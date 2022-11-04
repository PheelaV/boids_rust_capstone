// mod fv_capstone;
/*
    Boid algorithms: Crag Reynolds 1986, see http://www.red3d.com/cwr/
    Implementation inspired by:
    - Daniel Shiffman, Nature of Code, see https://natureofcode.com/book/
    - Nannou, creative coding library for Rust inspired by Processing,
      contains implemented examples of Nature of Code, see https://github.com/nannou-org/nannou/tree/master/nature_of_code
    - eisendaniel, see https://github.com/eisendaniel/boids/blob/master/src/main.rs
    - Craig Bester, see https://github.com/cycraig/rust-boids-wasm
*/
extern crate nannou;

use clap::Parser;
use cliargs::Args;
use nannou::{color::*, prelude::*};

use nannou_egui::{egui, Egui};

mod boid;
mod flock;
mod math_helpers;
mod options;
mod birdwatcher;
mod cliargs;

use flock::Flock;
use options::{Options, get_run_options, InitiationStrategy, Distance};
use birdwatcher::Birdwatcher;

use std::{sync::atomic::{AtomicU64, Ordering}, io, fs::OpenOptions};

// const NO_BOIDS: u32 = 10;
const WIDTH: u32 = 10000;
const HEIGHT: u32 = 1000;

fn main() {
    nannou::app(model)
    // .loop_mode(LoopMode::NTimes { number_of_updates: 60 })
    .update(update)
    .run();
}

struct ControlsState {
    execution_paused: bool,
    controls_open: bool,
}

struct Model {
    // boids: Vec<Boid>,
    egui: Egui,
    color: Hsv,
    flock: Flock,
    options: Options,
    last_update_millis: u128,
    fps: u8,
    fps_limit_on: bool,
    control_state: ControlsState,
    bird_watcher: Birdwatcher
}

fn model(app: &App) -> Model {
    // app.set_loop_mode(LoopMode::rate_fps(3.0));
    let args = Args::parse();
    let main_window = app
        .new_window()
        .key_pressed(key_pressed)
        .mouse_pressed(mouse_pressed)
        .closed(window_closed)
        .size(WIDTH, HEIGHT)
        .title("boids flocking")
        .raw_event(raw_window_event)
        .view(view) // This is where we forward all raw events for egui to process them
        .build()
        .unwrap();

    let window = app.window(main_window).unwrap();

    let run_options = get_run_options(&app.window_rect());

    let options = Options {
        init_boids: 100,
        init_width: 600,
        init_height: 600,
        run_options,
        initiation_strat: InitiationStrategy::RandomRandom,
        record_locations: true,
        no_boids: args.no_boids,
        // initiation_strat: InitiationStrategy::CircleCircumferenceIn
    };

    Model {
        egui: Egui::from_window(&window),
        color: Hsv::new(RgbHue::from_degrees(190.), 74., 32.),
        flock: Flock::new(args.no_boids, &options),
        options,
        last_update_millis: 0,
        control_state: ControlsState { 
            execution_paused: false,
            controls_open: false,
        },
        fps: 60,
        fps_limit_on: false,
        bird_watcher: Birdwatcher::new(),
    }
}

fn update(app: &App, model: &mut Model, update: Update) {
    if !should_render_update(app, model) {
        return
    }

    if app.keys.down.contains(&Key::Space) {
        return
    }

    let Model {
        ref mut egui,
        ref mut color,
        ref mut flock,
        ref mut options,
        ref mut bird_watcher,
        ..
    } = *model;

    let run_options = &mut options.run_options;

    // update window size as it could be resized
    let win = &app.window_rect();
    run_options.window.win_left =  win.left();
    run_options.window.win_right = win.right();
    run_options.window.win_top = win.top();
    run_options.window.win_bottom = win.bottom();
    run_options.window.win_h = win.h();
    run_options.window.win_w = win.w();

    // update controls UI
    egui.set_elapsed_time(update.since_start);
    let ctx = egui.begin_frame();
    egui::Window::new("controls")
        .default_size(egui::vec2(0.0, 200.0))
        .open(&mut model.control_state.controls_open)
        .show(&ctx, |ui| {
            ui.separator();
        
            ui.horizontal(|ui| {
                // ui.label("distancce");
                egui::ComboBox::from_label( "distance")
                    .selected_text(format!("{:?}", run_options.distance))
                    .show_ui(ui, |ui| {
                            ui.selectable_value(&mut run_options.distance, Distance::EucEnclosed, "EucEnclosed");
                            ui.selectable_value(&mut run_options.distance, Distance::EucThoroidal, "EucThoroidal");


                        }
                    );            
                });
            ui.horizontal(|ui| {
                ui.label("color");
                edit_hsv(ui, color);
            });
            
            ui.horizontal(|ui| {
                ui.label("baseline speed");
                ui.add(egui::Slider::new(&mut run_options.baseline_speed, 0.1..=20.0))
            });

            ui.horizontal(|ui| {
                ui.label("min speed");
                ui.add(egui::Slider::new(&mut run_options.min_speed, 0.0..=run_options.max_speed))
            });

            ui.horizontal(|ui| {
                ui.label("max speed");
                ui.add(egui::Slider::new(&mut run_options.max_speed, run_options.min_speed..=20.0))
            });

            ui.horizontal(|ui| {
                ui.label("max steering");
                ui.add(egui::Slider::new(&mut run_options.max_steering, 0.0..=10.0))
            });

            ui.separator();

            ui.horizontal(|ui| {
                ui.add(egui::Checkbox::new(&mut run_options.alignment_on, "allignment"));
                ui.add(egui::Checkbox::new(&mut run_options.cohesion_on, "cohesion"));
                ui.add(egui::Checkbox::new(&mut run_options.separation_on, "separation"));
            });

            ui.separator();

            ui.horizontal(|ui| {
                ui.label("allignment coef");
                ui.add(egui::Slider::new(&mut run_options.allignment_coefficient, 0.0..=10.0))
            });

            ui.horizontal(|ui| {
                ui.label("cohesion coef");
                ui.add(egui::Slider::new(&mut run_options.cohesion_coefficient, 0.0..=10.0))
            });

            ui.horizontal(|ui| {
                ui.label("separation coef");
                ui.add(egui::Slider::new(&mut run_options.separation_coefficient, 0.0..=10.0))
            });

            ui.separator();

            ui.horizontal(|ui| {
                ui.label("sensory distance");
                ui.add(egui::Slider::new(&mut run_options.sensory_distance, 50.0..=400.0))
            });

            ui.horizontal(|ui| {
                ui.label("allignment tresh coef");
                ui.add(egui::Slider::new(&mut run_options.allignment_treshold_coefficient, 0.0..=2.0))
            });

            ui.horizontal(|ui| {
                ui.label("cohesion tresh coef");
                ui.add(egui::Slider::new(&mut run_options.cohesion_treshold_coefficient, 0.0..=2.0))
            });

            ui.horizontal(|ui| {
                ui.label("separation tresh coef");
                ui.add(egui::Slider::new(&mut run_options.separation_treshold_coefficient, 0.0..=2.0))
            });

            ui.separator();

            ui.horizontal(|ui| {
                ui.label("size");
                ui.add(egui::Slider::new(&mut run_options.size, 5.0..=60.))
            });

            ui.horizontal(|ui| {
                ui.add(egui::Checkbox::new(&mut run_options.field_of_vision_on, "FOV (rad)"));

                ui.add_enabled(run_options.field_of_vision_on, egui::Slider::new(&mut run_options.field_of_vision_deg, 0.0..=360.))
            });

            ui.horizontal(|ui| {
                ui.add(egui::Checkbox::new(&mut model.fps_limit_on, "FPS limit"));

                ui.add_enabled(model.fps_limit_on, egui::Slider::new(&mut model.fps, 5..=120))
            });
        });

    // update sensory distances
    run_options.allignment_treshold_distance = run_options.sensory_distance * run_options.allignment_treshold_coefficient;
    run_options.cohesion_treshold_distance = run_options.sensory_distance * run_options.cohesion_treshold_coefficient;
    run_options.separation_treshold_distance = run_options.sensory_distance * run_options.separation_treshold_coefficient;
    run_options.field_of_vision_half_rad = run_options.field_of_vision_deg * PI/360.;

    // update model
    if model.control_state.execution_paused {
        return
    }
    flock.update(options);
    bird_watcher.watch(&flock);
    
    // println!("locations: {:#?}", flock.sample_locations());

    // let data = flock.sample_locations();

    // // let payload = serde_json::to_string(&data);
    // // let payload = csv::Writer::serialize(&mut self, record);
    // // let mut s = String::new();
    // let mut wtr = csv::Writer::from_writer(io::stdout());
    // wtr.serialize(data);
    // wtr.flush();
    // csv::Writer::from();

    // println!("json: {:#?}", payload);

}

fn raw_window_event(_app: &App, model: &mut Model, event: &nannou::winit::event::WindowEvent) {
    model.egui.handle_raw_event(event);
}

static RENDER_TICKS: AtomicU64 = AtomicU64::new(0);
static UPDATE_TICKS: AtomicU64 = AtomicU64::new(0);

fn should_render_view(model: & Model) -> bool {
    if !model.fps_limit_on {
        return true;
    }

    let render_ticks = RENDER_TICKS.load(Ordering::SeqCst);
    let update_ticks = UPDATE_TICKS.load(Ordering::SeqCst);

    let should_render = render_ticks != update_ticks;
    if should_render {
        RENDER_TICKS.store(update_ticks, Ordering::SeqCst)
    }

    should_render
}

fn should_render_update(app: &App, model: &mut Model) -> bool {
    if !model.fps_limit_on {
        return true;
    }

    let millis_now = app.duration.since_start.as_millis();
    
    let should_render = millis_now - model.last_update_millis > 1000 / model.fps as u128;

    if should_render {
        UPDATE_TICKS.fetch_add(1, Ordering::SeqCst);
        model.last_update_millis = millis_now;
    }

    should_render
}

fn key_pressed(_app: &App, model: &mut Model, key: Key) -> () {
    if key == Key::Space { 
        model.control_state.execution_paused = !model.control_state.execution_paused
    } else if key == Key::C {
        model.control_state.controls_open = !model.control_state.controls_open
    } else if key == Key::Key1 {
        model.options.run_options.alignment_on = !model.options.run_options.alignment_on;
    } else if key == Key::Key2 {
        model.options.run_options.cohesion_on = !model.options.run_options.cohesion_on;
    } else if key == Key::Key3 {
        model.options.run_options.separation_on = !model.options.run_options.separation_on;
    } else if key == Key::D {
        let id_delete = model.flock.delete_last();
        match id_delete {
            Some(x) => {
                if model.options.run_options.clicked_boid_id == x {
                    model.options.run_options.clicked_boid_id = std::u32::MAX;
                }
            },
            None => {},
        }

    } else if key == Key::I {
        model.flock.insert(&model.options);
    } else if key == Key::X {
        model.options.run_options.impl_mode = !model.options.run_options.impl_mode;
    }
}

fn mouse_pressed(app: &App, model: &mut Model, _button: MouseButton) -> () {
    if model.control_state.controls_open {
        return
    }
    let mouse_click = vec2(app.mouse.x, app.mouse.y);

    let clicked = model.flock.boids.iter()
    .filter(|b|{
        // get the distance from boid to click, if dist < size then boid is clicked
        let dist = b.position.distance(mouse_click);

        dist < model.options.run_options.size
    })
    .next();

    match clicked {
        Some(_) => {
            model.options.run_options.clicked_boid_id =  clicked.unwrap().id;
        },
        None => {},
    }

    // println!("mouse click: {:#?}", mouse_click);
}

fn window_closed(_app: &App, model: &mut Model) {
    let mut data = model.bird_watcher.pop_data();

    let file = OpenOptions::new()
        .write(true)
        .create(true)
        // .append(true)
        .open("./headings_data.csv")
        .unwrap();
    let mut wtr = csv::Writer::from_writer(file);
    // let mut wtr = csv::Writer::from_writer(io::stdout());

    data.sort_by(|a, b| a.id.cmp(&b.id));

    for data_point in data {
        let _ = wtr.serialize(data_point);
    }
    let _ = wtr.flush();
}


fn view(app: &App, model: &Model, frame: Frame) {
    if !should_render_view(model) {
        return
    }
    let draw = app.draw();


    draw.background().color(PLUM);

    // for b in &model.flock.boids{
    //     draw.ellipse()
    //     .xy(b.position)
    //     // .x_y(b.pos_x, b.pos_y)
    //     .radius(12.)
    //     .stroke(RED)
    //     .color(GREEN);
    // }
    // draw.ellipse()
    // .x(0.)
    // .y(0.)
    // .radius(60.)
    // .stroke(RED)
    // .color(RED);


    // let (a, b) = draw_animate(app, &draw, 1.);
    // let b = Boid::new(a, b, vec2(50., 50.), 10_000);
    // b.display(&draw, &c, &model.options.run_options);

    model.flock.draw(&draw, &model.color, &model.options.run_options);

    let c = hsv(0.5, 1., 1.);



    // draw.rect()
    // .xy(app.mouse.position())
    
    // .w_h(100.,100.)
    // .color(PINK);

    draw.to_frame(app, &frame).unwrap();
    model.egui.draw_to_frame(&frame).unwrap();
}

fn edit_hsv(ui: &mut egui::Ui, color: &mut Hsv) {
    let mut egui_hsv = egui::color::Hsva::new(
        color.hue.to_positive_radians() as f32 / (std::f32::consts::PI * 2.0),
        color.saturation,
        color.value,
        1.0,
    );

    if egui::color_picker::color_edit_button_hsva(
        ui,
        &mut egui_hsv,
        egui::color_picker::Alpha::Opaque,
    )
    .changed()
    {
        *color = nannou::color::hsv(egui_hsv.h, egui_hsv.s, egui_hsv.v);
    }
}

fn draw_animate(app: &App, draw: &Draw, offset: f32) -> (f32, f32) {
    // Generate sine wave data based on the time of the app
    let my_time = app.time + offset;

    let sine = (my_time * 1.0).sin();
    let slowersine = (my_time * 2.0).sin();

    // Get boundary of the window (to constrain the movements of our circle)
    let boundary = app.window_rect();
    let radius = 25.0;

    // Map the sine wave functions to ranges between the boundaries of the window
    let mut x = map_range(
        sine,
        -1.0,
        1.0,
        boundary.left() + radius,
        boundary.right() - radius,
    );
    let mut y = map_range(
        slowersine,
        -1.0,
        1.0,
        boundary.bottom() + radius,
        boundary.top() - radius,
    );

    if x > 0. {

        draw.ellipse().radius(radius).color(GREEN).x_y(0.0, 0.0);
    } else {
        draw.ellipse().radius(radius).color(GREEN).x_y(x, y);
    }

    // Clear the background to purple.
    // draw.background().color(PLUM);

    // Draw a blue ellipse at the x/y coordinates 0.0, 0.0
    // draw.ellipse().radius(radius).color(GREEN).x_y(x, y);

    (x, y)
}