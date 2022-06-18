use bevy::prelude::*;
use bevy_prototype_debug_lines::*;
use std::alloc;
use std::borrow::BorrowMut;
use std::cell::RefCell;
use std::mem;
use bevy_prototype_debug_lines::*;
use std::sync::Mutex;
use std::sync::Arc;
use std::u32::MIN;
use std::cmp;
use std::time::{Duration, Instant};

#[derive(Default, Clone, Copy)]
struct AABB {
    min : Vec2,
    max : Vec2
}

const SIMD_VECTOR_SIZE : usize = 4;
const UPDATE_INTERVAL : f32 = 1.0 / 120.0;
const COLLISION_ITERATIONS : usize = 2;
const COLLISION_RESPONSE_MULTIPLIER : f32 = 0.95;
const MAX_PENETRATION : f32 = 0.15;
const MAX_VELOCITY : f32 = 10.0;
const DUMPING : f32 = 0.99;
const CELL_SIZE : f32 = 5.0;
const GRAVITY : f32 = 450.0;

struct Particle {
    pub position_x : Vec<f32>,
    pub position_y : Vec<f32>,
    pub prev_position_x : Vec<f32>,
    pub prev_position_y : Vec<f32>,
    pub radius : Vec<f32>,
    pub alive_particles_count : u32,
}

impl Particle {
    fn new(capacity : usize) -> Self {
        const ALIGNMENT : usize = SIMD_VECTOR_SIZE * mem::size_of::<f32>();
        let size = mem::size_of::<f32>() * capacity;
        
        let mut p = Particle {
            alive_particles_count : 0,
            position_x: Vec::with_capacity(capacity),
            position_y: Vec::with_capacity(capacity),
            prev_position_x: Vec::with_capacity(capacity),
            prev_position_y: Vec::with_capacity(capacity),
            radius: Vec::with_capacity(capacity),
        };

        p.position_x.resize(capacity, 0.0);
        p.position_y.resize(capacity, 0.0);
        p.prev_position_x.resize(capacity, 0.0);
        p.prev_position_y.resize(capacity, 0.0);
        p.radius.resize(capacity, 0.0);

        return p;
    }

    fn add_particle(&mut self, position : Vec2, velocity : Vec2, radius : f32) {
        if self.alive_particles_count as usize == self.position_x.len() {
            panic!("Particle capacity exceeded");
        }

        let index = self.alive_particles_count as usize;
        let prev_pos = position - velocity;

        self.position_x[index] = position.x;
        self.position_y[index] = position.y;
        self.prev_position_x[index] = prev_pos.x;
        self.prev_position_y[index] = prev_pos.y;
        self.radius[index] = radius;
        self.alive_particles_count += 1;
    }

}

impl Default for Particle {
    fn default() -> Self { Particle::new(100) }
}

struct GridCell {
    pub particles : Vec<u32>
}

#[derive(Default)]
struct Grid {
    cell_dimensions : UVec2,
    cell_size : Vec2,
    cells : Vec<GridCell>,
    bounds : AABB
}

impl Grid {
    pub fn new(width : f32, height: f32, cell_size : f32, bounds : AABB) -> Self{
        let cell_dimensions = UVec2::new((width / cell_size).ceil() as u32, (height / cell_size).ceil() as u32);

        let mut grid = Grid {
            cell_size : Vec2::new(cell_size, cell_size),
            cell_dimensions,
            bounds,
            ..Default::default()
        };

        grid.cells.reserve_exact((cell_dimensions.x * cell_dimensions.y) as usize);
        for i in 0..grid.cells.capacity() {
            let mut cell = GridCell{ particles: Vec::new() };
            cell.particles.reserve(20);
            grid.cells.push(cell);
        }

        return grid;
    }

    pub fn position_to_cell(&self, position : Vec2) -> UVec2 {

        //return UVec2::new(self.cell_dimensions.x - 1, self.cell_dimensions.y - 1);

        let position = position - self.bounds.min;
        UVec2::new(
            (position.x / (self.cell_size.x as f32)) as u32,
            (position.y / (self.cell_size.y as f32)) as u32
        ).min(UVec2::new(self.cell_dimensions.x - 1, self.cell_dimensions.y - 1))
    }

    pub fn get_cell_index(&self, cell_position : UVec2) ->usize {
        return (cell_position.x + cell_position.y * self.cell_dimensions.x) as usize;
    }

    pub fn add_particle(&mut self, position : Vec2, index : u32) {
        let cell_index = self.get_cell_index(self.position_to_cell(position));
        //if self.cells[cell_index].particles.len() < 6 {
            self.cells[cell_index].particles.push(index);
        //}

        //print!("adding to grid {}, index = {}\n", self.position_to_cell(position), index);
    }

    pub fn clear(&mut self) {
        for cell in self.cells.iter_mut() {
            cell.particles.clear();
        }
    }

    pub fn get_adjacent_range(&self, cell_position : UVec2) -> (UVec2, UVec2) {
        let mut min = IVec2::new(cell_position.x as i32, cell_position.y as i32) - IVec2::new(1, 1);
        let mut max = IVec2::new(cell_position.x as i32, cell_position.y as i32) + IVec2::new(1, 1);

        min = min.max(IVec2::new(0, 0));
        max = max.min(IVec2::new(self.cell_dimensions.x as i32 - 1, self.cell_dimensions.y as i32 - 1));

        return (UVec2::new(min.x as u32, min.y as u32), UVec2::new(max.x as u32, max.y as u32));
    }
}

#[derive(Default)]
pub struct SimWorld {
    bounds : AABB,
    particles : Particle,
    grid : Grid,
    accumulated_dt : f32,
    tick_duration : f32
}

impl SimWorld {
    pub fn new(width : f32, height: f32, capacity: usize) -> Self {
        let bounds = AABB {
            min : Vec2::new(-width / 2.0, -height/2.0),
            max : Vec2::new(width / 2.0, height/2.0)
        };

        SimWorld {
            bounds,
            particles: Particle::new(capacity),
            grid : Grid::new(width, height, CELL_SIZE, bounds),
            ..Default::default()
        }
    }

    pub fn add_particle(&mut self, position : Vec2, velocity : Vec2, radius : f32) {
        self.particles.add_particle(position, velocity, radius);
    }

    fn sim_apply_gravity(&mut self, g : Vec2, dt : f32) {
        for i in 0..self.particles.alive_particles_count as usize {
            let g = g * dt * dt;
            self.particles.prev_position_x[i] -= g.x;
            self.particles.prev_position_y[i] -= g.y;
        }
    }

    fn sim_move(&mut self) {
        for i in 0..self.particles.alive_particles_count as usize {
            let mut velocity = Vec2::new(self.particles.position_x[i] - self.particles.prev_position_x[i], self.particles.position_y[i] - self.particles.prev_position_y[i]);
            let sq_len = velocity.length_squared();
            if (sq_len > MAX_VELOCITY * MAX_VELOCITY) {
                let corrected_velocity = velocity / sq_len.sqrt() * MAX_VELOCITY;
                self.particles.prev_position_x[i] = self.particles.position_x[i] - corrected_velocity.x;
                self.particles.prev_position_y[i] = self.particles.position_y[i] - corrected_velocity.y;
                velocity = corrected_velocity;
            }

            let prev_pos_x = self.particles.position_x[i];
            let prev_pos_y = self.particles.position_y[i];
            self.particles.position_x[i] += velocity.x * DUMPING;
            self.particles.position_y[i] += velocity.y * DUMPING;
            self.particles.prev_position_x[i] = prev_pos_x;
            self.particles.prev_position_y[i] = prev_pos_y;
        }
    }
    
    fn sim_restrict_bounds(&mut self) {
        for i in 0..self.particles.alive_particles_count as usize {
            self.particles.position_x[i] = self.bounds.max.x.min(self.particles.position_x[i]);
            self.particles.position_x[i] = self.bounds.min.x.max(self.particles.position_x[i]);
            self.particles.position_y[i] = self.bounds.max.y.min(self.particles.position_y[i]);
            self.particles.position_y[i] = self.bounds.min.y.max(self.particles.position_y[i]);
        }
    }

    fn sim_collide_slow(&mut self) {
        for i in 0..self.particles.alive_particles_count as usize {
            for j in i+1..self.particles.alive_particles_count as usize {
                let p1 = Vec2::new(self.particles.position_x[i], self.particles.position_y[i]);
                let p2 = Vec2::new(self.particles.position_x[j], self.particles.position_y[j]);
                let r1r2 = self.particles.radius[i] + self.particles.radius[j];
                let delta = p2 - p1;
                let sq_distance = delta.length_squared();
                if (sq_distance < r1r2 * r1r2) {
                    let distance = sq_distance.sqrt();
                    let mut delta_norm = Vec2::new(1.0, 0.0);
                    let shift = (r1r2 - distance).min(MAX_PENETRATION * r1r2) / 2.0 * COLLISION_RESPONSE_MULTIPLIER;
                    if distance > 1e-5 {
                        delta_norm = delta / distance;
                    }

                    let p1_new = p1 - delta_norm * shift;
                    let p2_new = p2 + delta_norm * shift;
                    self.particles.position_x[i] = p1_new.x;
                    self.particles.position_y[i] = p1_new.y;
                    self.particles.position_x[j] = p2_new.x;
                    self.particles.position_y[j] = p2_new.y;
                }
            }
        }
    }

    fn sim_collide(&mut self) {
        for i in 0..self.particles.alive_particles_count as usize {
            let p1 = Vec2::new(self.particles.position_x[i], self.particles.position_y[i]);
            let cell_position = self.grid.position_to_cell(p1);
            let range = self.grid.get_adjacent_range(cell_position);

            for x in range.0.x..=range.1.x {
                for y in range.0.y..=range.1.y {
                    //print!("Checking vs cell {}\n", UVec2::new(x, y));
                    
                    //let ux = x - range.0.x;
                    //let uy = y - range.0.y;
                    //if ux % 2 == 0 && uy % 2 == 0 { continue; }

                    for j in &self.grid.cells[self.grid.get_cell_index(UVec2::new(x, y))].particles {
                        let p1 = Vec2::new(self.particles.position_x[i], self.particles.position_y[i]);

                        let j = *j as usize;
                        if i != j {
                            let p2 = Vec2::new(self.particles.position_x[j], self.particles.position_y[j]);
                            let r1r2 = self.particles.radius[i] + self.particles.radius[j];
                            let delta = p2 - p1;
                            let sq_distance = delta.length_squared();
                            if (sq_distance < r1r2 * r1r2) {
                                let distance = sq_distance.sqrt();
                                let mut delta_norm = Vec2::new(1.0, 0.0);
                                let shift = ((r1r2 - distance) * COLLISION_RESPONSE_MULTIPLIER).min(MAX_PENETRATION * r1r2) / 2.0;
                                if distance > 1e-5 {
                                    delta_norm = delta / distance;
                                }
            
                                let p1_new = p1 - delta_norm * shift;
                                let p2_new = p2 + delta_norm * shift;
                                self.particles.position_x[i] = p1_new.x;
                                self.particles.position_y[i] = p1_new.y;
                                self.particles.position_x[j] = p2_new.x;
                                self.particles.position_y[j] = p2_new.y;
                            }
                       }

                    }
                }
            }
        }
    }

    pub fn get_particle_pos(&mut self, index : usize) -> Vec2 {
        if index >= self.particles.alive_particles_count as usize {
            panic!("particle index out of range");
        }
        let pos = Vec2::new(self.particles.position_x[index], self.particles.position_y[index]);
        let prev_pos = Vec2::new(self.particles.prev_position_x[index], self.particles.prev_position_y[index]);

        return prev_pos + (pos - prev_pos) * (self.accumulated_dt / UPDATE_INTERVAL);
    }

    pub fn get_collision_iterations(&self) -> u32 { COLLISION_ITERATIONS as u32 }
    pub fn get_particle_count(&self) -> usize { return self.particles.alive_particles_count as usize; }

    pub fn update_grid(&mut self) {
        self.grid.clear();
        for i in 0..self.particles.alive_particles_count as usize {
            self.grid.add_particle(Vec2::new(self.particles.position_x[i], self.particles.position_y[i]), i as u32);
        }
    }

    pub fn tick(&mut self, dt : f32) {
        self.accumulated_dt += dt;

        let now = Instant::now();
        let mut elapsed = 0.0;
        loop {
            if self.accumulated_dt < UPDATE_INTERVAL {
                break;
            }

            self.accumulated_dt -= UPDATE_INTERVAL;
            self.sim_apply_gravity(Vec2::new(0.0, -GRAVITY), UPDATE_INTERVAL);
            self.sim_move();

            self.update_grid();
            for i in 0..COLLISION_ITERATIONS {
                self.sim_collide();
                self.sim_restrict_bounds();
            }

            elapsed = (now.elapsed().as_secs() as f64 + now.elapsed().subsec_nanos() as f64 * 1e-9) as f32;
            if elapsed >= 1.0 / 30.0 {
                break;
            }
        }

        self.tick_duration = elapsed;
    }

    pub fn get_tick_duration(&self) -> f32 { return self.tick_duration; }

    pub fn apply_force_to_point(&mut self, position : Vec2, radius : f32, strength: f32, dt : f32) {
        let sq_radius = radius * radius;
        
        for i in 0..self.particles.alive_particles_count as usize {
            let delta = position - Vec2::new(self.particles.position_x[i], self.particles.position_y[i]);
            let sq_len = delta.length_squared();
            if (sq_len < sq_radius) {
                let mut force = Vec2::new(1.0, 0.0);
                if sq_len > 1e-5 {
                    force = delta * dt * strength;
                }

                self.particles.prev_position_x[i] -= force.x;
                self.particles.prev_position_y[i] -= force.y;
            }
        }
    }

    pub fn get_bounds(&self) -> (Vec2, Vec2) {
        return (self.bounds.min, self.bounds.max);
    }

    pub fn draw_grid(&self, debug_lines : &mut DebugLines) {

        let size = self.bounds.max - self.bounds.min;

        for i in 0..self.grid.cell_dimensions.x as usize {
            let p1 = self.bounds.min + Vec2::new(self.grid.cell_size.x * i as f32, 0.0);
            let p2 = p1 + Vec2::new(0.0, size.y);
            debug_lines.line(Vec3::new(p1.x, p1.y, 0.0), Vec3::new(p2.x, p2.y, 0.0), 0.0);
        }

        for i in 0..self.grid.cell_dimensions.y as usize {
            let p1 = self.bounds.min + Vec2::new(0.0, self.grid.cell_size.y * i as f32);
            let p2 = p1 + Vec2::new(size.x, 0.0);
            debug_lines.line(Vec3::new(p1.x, p1.y, 0.0), Vec3::new(p2.x, p2.y, 0.0), 0.0);
        }
    }
}