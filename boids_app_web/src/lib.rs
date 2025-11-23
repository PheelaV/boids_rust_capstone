/// Boids simulation using Nannou framework compiled to WebAssembly
///
/// This is a browser-based version of the desktop boids_app, using nannou's
/// experimental WASM support to render the simulation with WebGPU/WebGL.
///
/// Features:
/// - Full nannou rendering (same graphics as desktop)
/// - Interactive controls (mouse, keyboard)
/// - Dynamic boid count management
/// - Behavior toggles (separation, cohesion, alignment, wander)
///
/// Note: This uses nannou's experimental WASM support. Some features from the
/// desktop app may not be available (file I/O, audio, etc.).

use boids_lib::{
    boid::Boid,
    flock::Flock,
    options::{RunOptions, TrackerType},
};
use nannou::prelude::*;
use wasm_bindgen::prelude::*;

// When the `wee_alloc` feature is enabled, use `wee_alloc` as the global
// allocator.
#[cfg(feature = "wee_alloc")]
#[global_allocator]
static ALLOC: wee_alloc::WeeAlloc = wee_alloc::WeeAlloc::INIT;

/// Model struct holds the simulation state
pub struct Model<'a> {
    flock: Flock<'a>,
    run_options: RunOptions,
    color: Hsv,
    paused: bool,
    frame_count: u64,
}

/// Initialize the simulation
fn model(_app: &App) -> Model<'static> {
    let mut run_options = RunOptions::default();
    run_options.init_boids = 200;
    run_options.tracker_type = TrackerType::SpatHash;

    // Create the flock
    let flock = Flock::new(&run_options);

    Model {
        flock,
        run_options,
        color: hsv(0.6, 0.8, 1.0),
        paused: false,
        frame_count: 0,
    }
}

/// Update the simulation state
fn update(_app: &App, model: &mut Model, _update: Update) {
    if !model.paused {
        model.flock.update(&mut model.run_options);
        model.frame_count += 1;
    }
}

/// Render the simulation
fn view(app: &App, model: &Model, frame: Frame) {
    let draw = app.draw();

    // Clear the background
    draw.background().color(BLACK);

    // Draw all boids
    for (boid, metadata) in model.flock.view2() {
        draw_boid(&draw, boid, &model.color, &model.run_options);

        // Highlight selected boid if any
        if boid.id == model.run_options.clicked_boid_id {
            // Draw selection circle
            draw.ellipse()
                .xy(boid.position)
                .radius(model.run_options.size * 2.0)
                .stroke(HOTPINK)
                .stroke_weight(2.0)
                .no_fill()
                .z(0.5);
        }

        // Color neighbors differently
        if metadata.clicked_neighbour_id != std::usize::MAX {
            draw.ellipse()
                .xy(boid.position)
                .radius(model.run_options.size * 1.5)
                .color(rgba(0.5, 1.0, 0.5, 0.3))
                .z(0.4);
        }
    }

    // Draw FPS and info
    let fps_text = format!("FPS: {:.0} | Boids: {} | Frame: {}",
        app.fps(),
        model.flock.get_no_entities(),
        model.frame_count
    );
    draw.text(&fps_text)
        .xy(app.window_rect().pad(20.0).top_left())
        .color(WHITE)
        .font_size(16)
        .left_justify();

    if model.paused {
        draw.text("PAUSED")
            .xy(app.window_rect().xy())
            .color(RED)
            .font_size(48);
    }

    // Write to the frame
    draw.to_frame(app, &frame).unwrap();
}

/// Draw a single boid as a triangle pointing in its velocity direction
fn draw_boid(draw: &Draw, boid: &Boid, color: &Hsv, run_options: &RunOptions) {
    let position = boid.position;
    let mut theta = boid.velocity.angle();

    // Handle NaN angles
    if theta.is_nan() {
        theta = 0.0;
    }

    // Create triangle vertices (arrow shape)
    let size = run_options.size;
    let vertices = vec![
        pt2(-0.8 * size, 0.6 * size),
        pt2(size, 0.0),
        pt2(-0.8 * size, -0.6 * size),
        pt2(-0.5 * size, 0.0),
    ];

    draw.polygon()
        .stroke(AZURE)
        .points(vertices)
        .xy(position)
        .z(1.0)
        .rotate(2.0 * PI + theta)
        .hsv(
            color.hue.to_positive_degrees() / 360.0,
            color.saturation,
            color.value,
        );
}

/// Handle keyboard input
fn key_pressed(_app: &App, model: &mut Model, key: Key) {
    match key {
        // Space: pause/resume
        Key::Space => {
            model.paused = !model.paused;
        },
        // I: double boid count
        Key::I => {
            let current_count = model.flock.get_no_entities();
            for _ in 0..current_count {
                model.flock.insert_single(&model.run_options);
            }
        },
        // D: halve boid count
        Key::D => {
            let current_count = model.flock.get_no_entities();
            let to_remove = current_count / 2;
            for _ in 0..to_remove {
                model.flock.delete_last(&model.run_options);
            }
        },
        // R: reset simulation
        Key::R => {
            model.flock.restart(&model.run_options);
            model.frame_count = 0;
        },
        // 1: toggle separation
        Key::Key1 => {
            model.run_options.separation_on = !model.run_options.separation_on;
        },
        // 2: toggle cohesion
        Key::Key2 => {
            model.run_options.cohesion_on = !model.run_options.cohesion_on;
        },
        // 3: toggle alignment
        Key::Key3 => {
            model.run_options.alignment_on = !model.run_options.alignment_on;
        },
        // 4: toggle wander
        Key::Key4 => {
            model.run_options.wander_on = !model.run_options.wander_on;
        },
        _ => {}
    }
}

/// Handle mouse clicks to select boids
fn mouse_pressed(app: &App, model: &mut Model, button: MouseButton) {
    if button == MouseButton::Left {
        let mouse_pos = app.mouse.position();

        // Find closest boid to click position
        let mut closest_dist = f32::MAX;
        let mut closest_id = std::usize::MAX;

        for (boid, _) in model.flock.view2() {
            let dist = boid.position.distance(mouse_pos);
            if dist < closest_dist && dist < 50.0 {  // Within 50 pixels
                closest_dist = dist;
                closest_id = boid.id;
            }
        }

        model.run_options.clicked_boid_id = closest_id;
    }
}

/// Main entry point for WASM
#[wasm_bindgen(start)]
pub fn run() -> Result<(), JsValue> {
    // Set panic hook for better error messages in browser console
    #[cfg(feature = "console_error_panic_hook")]
    console_error_panic_hook::set_once();

    // Start the nannou app
    nannou::app(model)
        .update(update)
        .simple_window(view)
        .key_pressed(key_pressed)
        .mouse_pressed(mouse_pressed)
        .run();

    Ok(())
}
