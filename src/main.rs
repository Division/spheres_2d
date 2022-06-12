use bevy::{
    prelude::*, render::camera::Camera2d,
    diagnostic::{Diagnostics, FrameTimeDiagnosticsPlugin},
};
use bevy_prototype_debug_lines::*;
use bevy::window::CursorMoved;

const SPHERE_COUNT : i32 = 8000;
const ARENA_SIZE : (f32, f32) = (1000.0, 700.0);
const SPRITE_RADIUS : f32 = 5.0;

mod simulation;

#[derive(Component)]
struct Sphere;

struct MouseLoc(Vec2);

#[derive(Component)]
struct TextChanges;

fn main() {
    let mut sim_world = simulation::SimWorld::new(ARENA_SIZE.0, ARENA_SIZE.1, 10000);

    App::new()
        .add_plugins(DefaultPlugins)
        .insert_resource(sim_world)
        .insert_resource(MouseLoc(Vec2::new(0.0, 0.0)))
        .add_plugin(DebugLinesPlugin::with_depth_test(false))
        .add_plugin(FrameTimeDiagnosticsPlugin)
        .add_startup_system(setup)
        .add_system(mouse_events)
        .add_system(simulation)
        .add_system(debug_info_system)
        .add_system(some_system)
        .run();
}

fn mouse_events(
    mut cursor_moved_events: EventReader<CursorMoved>,
    mut mouse_loc: ResMut<MouseLoc>,
    mut windows: ResMut<Windows>,
    mut camera: Query<&Transform, With<Camera2d>>
) {
    let window = windows.primary_mut();
    let size = Vec2::new(window.width() as f32, window.height() as f32);

    let camera = camera.get_single().unwrap();

    for event in cursor_moved_events.iter() {
        let mut n = event.position - size / 2.0;
        let mut v = *camera * Vec3::new(n.x, n.y, 0.0);
        mouse_loc.as_mut().0 = Vec2::new(v.x, v.y);
    }
}

fn setup(mut commands: Commands, mut sim_world: ResMut<simulation::SimWorld>, asset_server: Res<AssetServer>) {

    commands.spawn_bundle(OrthographicCameraBundle::new_2d());

    let sim_world = sim_world.as_mut();
    let bounds = sim_world.get_bounds();
    let half_size = (bounds.1 - bounds.0) / 2.0;
    let offset = Vec2::new(SPRITE_RADIUS * 5.0, SPRITE_RADIUS * 5.0);

    for i in 0..SPHERE_COUNT {
        let mut t = Transform::identity();
        let pos = Vec2::new((i % 50) as f32 * SPRITE_RADIUS * 2.0, (i / 50) as f32 * SPRITE_RADIUS * 2.0) - half_size + offset;
        t.translation.x = pos.x;
        t.translation.y = pos.y;

        sim_world.add_particle(pos, Vec2::new(0.5, -1.0), SPRITE_RADIUS);

        commands.spawn_bundle(SpriteBundle {
            sprite: Sprite { custom_size: Some(Vec2::new(SPRITE_RADIUS * 2.0, SPRITE_RADIUS * 2.0)), ..default()},
            texture: asset_server.load("icon.png"),
            ..default()
        }).insert(Sphere);
    }

    let font = asset_server.load("FiraSans-Bold.ttf");
    let font_size = 15.0;
    commands.spawn_bundle(UiCameraBundle::default());
    commands.spawn_bundle(TextBundle {
        style: Style {
            align_self: AlignSelf::FlexEnd,
            position_type: PositionType::Absolute,
            position: Rect {
                left: Val::Px(5.0),
                top: Val::Px(5.0),
                ..default()
            },
            ..default()
        },
        text: Text {
            sections: vec![
                TextSection {
                    value: "This text changes in the bottom right".to_string(),
                    style: TextStyle {
                        font: font.clone(),
                        font_size: font_size,
                        color: Color::RED,
                    },
                }
            ],
            alignment: Default::default(),
        },
        ..default()
    })
    .insert(TextChanges);
}

fn simulation(
    mouse_button_input: Res<Input<MouseButton>>,
    time: Res<Time>,
    mouse_loc : Res<MouseLoc>,
    mut sim_world: ResMut<simulation::SimWorld>,
    mut sprite_position: Query<&mut Transform, With<Sphere>>
) {
    let sim_world = sim_world.as_mut();

    let mouse_loc = mouse_loc.0;

    if mouse_button_input.pressed(MouseButton::Left) {
        sim_world.apply_force_to_point(mouse_loc, 50.0, 2.0, time.delta_seconds());
    }

    if mouse_button_input.pressed(MouseButton::Right) {
        sim_world.apply_force_to_point(mouse_loc, 50.0, -2.0, time.delta_seconds());
    }

    sim_world.tick(time.delta_seconds());

    let mut index = 0;
    for (mut transform) in sprite_position.iter_mut() {
        if index >= sim_world.get_particle_count() {
            break;
        }

        let pos = sim_world.get_particle_pos(index);
        transform.translation.x = pos.x;
        transform.translation.y = pos.y;
        index += 1;
    }
}

fn some_system(mut lines: ResMut<DebugLines>, mut sim_world: ResMut<simulation::SimWorld>) {
    let start = Vec3::new(-ARENA_SIZE.0 / 2.0, -ARENA_SIZE.1 / 2.0, 0.0);
    let end = Vec3::new(ARENA_SIZE.0 / 2.0, ARENA_SIZE.1 / 2.0, 0.0);
    lines.line(start, Vec3::new(start.x, end.y, 0.0), 0.0);
    lines.line(start, Vec3::new(end.x, start.y, 0.0), 0.0);
    lines.line(end, Vec3::new(start.x, end.y, 0.0), 0.0);
    lines.line(end, Vec3::new(end.x, start.y, 0.0), 0.0);

    let sim_world = sim_world.as_mut();
    sim_world.draw_grid(lines.as_mut());
}

fn debug_info_system(
    time: Res<Time>,
    diagnostics: Res<Diagnostics>,
    sim_world: Res<simulation::SimWorld>,
    mut query: Query<&mut Text, With<TextChanges>>,
) {
    for mut text in query.iter_mut() {
        let mut fps = 0.0;
        if let Some(fps_diagnostic) = diagnostics.get(FrameTimeDiagnosticsPlugin::FPS) {
            if let Some(fps_avg) = fps_diagnostic.average() {
                fps = fps_avg;
            }
        }

        let mut frame_time = time.delta_seconds_f64();
        if let Some(frame_time_diagnostic) = diagnostics.get(FrameTimeDiagnosticsPlugin::FRAME_TIME)
        {
            if let Some(frame_time_avg) = frame_time_diagnostic.average() {
                frame_time = frame_time_avg;
            }
        }

        text.sections[0].value = format!(
            "Tick {:.2}ms, {:.1} fps, {:.3} ms/frame\nParticle count: {}\nCollision iterations: {}",
            sim_world.get_tick_duration() * 1000.0,
            fps,
            frame_time * 1000.0,
            sim_world.get_particle_count(),
            sim_world.get_collision_iterations()
        );
    }
}