extern crate nannou;
use boids_lib::birdwatcher::Birdwatcher;
use boids_lib::boid::{Boid, BoidMetadata};
use boids_lib::flock::Flock;
use boids_lib::flock::spathash_tracker::SpatHash1D;
use boids_lib::math_helpers::{distance_dyn, tor_vec_p};
use boids_lib::options::{self, Distance, RunOptions, SaveOptions, WindowSize, Boundary};
use clap_serde_derive::{clap::Parser, ClapSerde};
use nannou::draw::properties::ColorScalar;
use nannou::geom::{Ellipse, Tri};
use nannou::{color::*, prelude::*};
use nannou_egui::{egui, Egui};
use std::fs;
use std::sync::atomic::{AtomicU64, Ordering};

mod cliargs;
use cliargs::{Args, Config};

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

pub struct Model {
    egui: Egui,
    color: Hsv,
    flock: Flock,
    run_options: RunOptions,
    last_update_micros: u128,
    last_render_micros: u128,
    update_ticks: u128,
    since_last_update_micros: u128,
    fps: u8,
    fps_limit_on: bool,
    update_timing_on: bool,
    control_state: ControlsState,
    bird_watcher: Birdwatcher,
    window_id: WindowId,
    debug_grid: bool,
    debug_distance: bool,
    debug_labels: bool,
    repulsive_distance: f32,
    repulsive_force: f32
}

fn model(app: &App) -> Model {
    // let args = Args::parse();
    // Parse whole args with clap
    let mut args = Args::parse();

    let config = if let Ok(s) = fs::read_to_string(&args.config_path) {
        match toml::from_str::<<Config as ClapSerde>::Opt>(&s) {
            Ok(c) => Config::from(c).merge(&mut args.config),
            Err(err) => panic!("Error in configuration file:\n{}", err),
        }
    } else {
        Config::from(&mut args.config)
    };

    let mut run_options: RunOptions = Default::default();

    run_options.init_boids = config.no_boids;
    run_options.sample_rate = config.sample_rate;

    let save_options = SaveOptions {
        save_locations: config.save,
        save_locations_timestamp: config.save_timestamp,

        // default
        save_locations_path: run_options.save_options.save_locations_path,
    };

    run_options.window = options::get_window_size(config.init_width, config.init_height);

    run_options.sensory_distance = config.sensory_distance;
    run_options.allignment_coefficient = config.allignment_coefficient;
    run_options.cohesion_coefficient = config.cohesion_coefficient;
    run_options.separation_coefficient = config.separation_coefficient;
    run_options.allignment_treshold_coefficient = config.allignment_treshold_coefficient;
    run_options.cohesion_treshold_coefficient = config.cohesion_treshold_coefficient;
    run_options.separation_treshold_coefficient = config.separation_treshold_coefficient;
    run_options.min_speed = config.min_speed;
    run_options.max_speed = config.max_speed;
    run_options.max_steering = config.max_steering;
    run_options.field_of_vision_deg = config.field_of_vision;

    run_options.save_options = save_options;
    run_options.dbscan_flock_clustering_on = config.dbscan_flock_clustering_on;

    run_options.wander_on = config.wander_on;
    run_options.wander_coefficient = config.wander_coefficient;
    run_options.wander_distance = config.wander_distance;
    run_options.wander_radius = config.wander_radius;
    run_options.wander_rate = config.wander_rate;
    run_options.size = config.size;

    let main_window = app
        .new_window()
        .key_pressed(key_pressed)
        .key_released(key_released)
        .mouse_pressed(mouse_pressed)
        .closed(window_closed)
        .resizable(false)
        .size(
            run_options.window.win_w as u32,
            run_options.window.win_h as u32,
        )
        .title("boids flocking")
        .raw_event(raw_window_event)
        .view(view)
        .build()
        .unwrap();

    let window = app.window(main_window).unwrap();

    let birdwatcher = Birdwatcher::new(run_options.sample_rate);

    Model {
        egui: Egui::from_window(&window),
        color: Hsv::new(RgbHue::from_degrees(316.), 1., 0.94),
        flock: Flock::new(&run_options),
        run_options,
        last_update_micros: 0,
        last_render_micros: 0,
        update_ticks: 0,
        since_last_update_micros: 0,
        control_state: ControlsState {
            execution_paused: false,
            controls_open: false,
        },
        fps: 120,
        fps_limit_on: false,
        update_timing_on: true,
        bird_watcher: birdwatcher,
        window_id: main_window,
        debug_grid: false,
        debug_distance: false,
        repulsive_distance: 100.,
        repulsive_force: 0.05,
        debug_labels: false,
    }
}

fn update(app: &App, model: &mut Model, update: Update) {

    
    // Record update timing if we are looking at upate times or limiting frame times
    if model.fps_limit_on || model.update_timing_on {
        record_update_timing(app, model);

        // If we are limitting FPS, decide whether the update should run
        if model.fps_limit_on && !should_render_update(model) {
            return;
        }
    }

    //TODO: might not be needed
    if app.keys.down.contains(&Key::Space) {
        return;
    }

    let Model {
        ref mut egui,
        ref mut color,
        ref mut flock,
        ref mut run_options,
        ref mut bird_watcher,
        ref mut repulsive_distance,
        ref mut repulsive_force,
        ..
    } = *model;

    // update window size as it could be resized
    let win = &app.window_rect();

    let new_win = {
        let win_left  = win.left() as i32;
        let win_right = win.right() as i32;
        let wind_top = win.top() as i32;
        let win_bottom = win.bottom() as i32;
        let win_h = win.h() as i32;
        let win_w = win.w() as i32;
        WindowSize { win_left: win_left, win_right: win_right, win_top: wind_top, win_bottom: win_bottom, win_h: win_h, win_w: win_w }
    };

    // if the window size changed
    if new_win.win_h != run_options.window.win_h || new_win.win_w != run_options.window.win_w {
        // we just want to update the setting and wait till the window resizing stops
        run_options.window = new_win;
        return;
    }

    // update controls UI
    egui.set_elapsed_time(update.since_start);
    let ctx = egui.begin_frame();
    egui::Window::new("controls")
        .default_size(egui::vec2(0.0, 200.0))
        .open(&mut model.control_state.controls_open)
        .show(&ctx, |ui| {
            ui.separator();

            ui.horizontal(|ui| {
                egui::ComboBox::from_label("distance")
                    .selected_text(format!("{:?}", run_options.distance))
                    .show_ui(ui, |ui| {
                        ui.selectable_value(
                            &mut run_options.distance,
                            Distance::EucEnclosed,
                            "EucEnclosed",
                        );
                        ui.selectable_value(
                            &mut run_options.distance,
                            Distance::EucToroidal,
                            "EucThoroidal",
                        );
                    });
            });
            ui.horizontal(|ui| {
                egui::ComboBox::from_label("boundary")
                    .selected_text(format!("{:?}", run_options.boundary))
                    .show_ui(ui, |ui| {
                        ui.selectable_value(
                            &mut run_options.boundary,
                            Boundary::Reflective,
                            "Reflective",
                        );
                        ui.selectable_value(
                            &mut run_options.boundary,
                            Boundary::Toroidal,
                            "Toroidal",
                        );
                        ui.selectable_value(
                            &mut run_options.boundary,
                            Boundary::Absorbing,
                            "Absorbing",
                        );
                        ui.selectable_value(
                            &mut run_options.boundary,
                            Boundary::Repulsive { distance: *repulsive_distance, force: *repulsive_force },
                            "Repulsive",
                        );
                    });
                });

            if let Boundary::Repulsive { distance: _, force: _ } = run_options.boundary {
                ui.horizontal(|ui| {
                    ui.label("rep distance");
                    ui.add(egui::Slider::new(
                        repulsive_distance,
                        0.0..=300_f32,
                    ));
                });
                ui.horizontal(|ui| {
                    ui.label("rep force");
                    ui.add(egui::Slider::new(
                        repulsive_force,
                        0.001..=0.2,
                    ));
                    if ui.button("apply").clicked() {
                        run_options.boundary = Boundary::Repulsive { distance: *repulsive_distance, force: *repulsive_force }
                    }
                });
            }

            ui.horizontal(|ui| {
                ui.label("color");
                edit_hsv(ui, color);
            });

            ui.horizontal(|ui| {
                ui.label("baseline speed");
                ui.add(egui::Slider::new(
                    &mut run_options.baseline_speed,
                    0.1..=3.0,
                ))
            });

            ui.horizontal(|ui| {
                ui.label("min speed");
                ui.add(egui::Slider::new(
                    &mut run_options.min_speed,
                    0.0..=run_options.max_speed,
                ))
            });

            ui.horizontal(|ui| {
                ui.label("max speed");
                ui.add(egui::Slider::new(
                    &mut run_options.max_speed,
                    run_options.min_speed..=5.0,
                ))
            });

            ui.horizontal(|ui| {
                ui.label("max steering");
                ui.add(egui::Slider::new(&mut run_options.max_steering, 0.0..=10.0))
            });

            ui.separator();

            ui.horizontal(|ui| {
                ui.add(egui::Checkbox::new(
                    &mut run_options.alignment_on,
                    "allignment",
                ));
                ui.add(egui::Checkbox::new(
                    &mut run_options.cohesion_on,
                    "cohesion",
                ));
                ui.add(egui::Checkbox::new(
                    &mut run_options.separation_on,
                    "separation",
                ));
                ui.add(egui::Checkbox::new(&mut run_options.wander_on, "wander"));
            });

            ui.separator();

            ui.horizontal(|ui| {
                ui.label("wander coef");
                ui.add(egui::Slider::new(
                    &mut run_options.wander_coefficient,
                    0.0..=1.,
                ))
            });

            ui.horizontal(|ui| {
                ui.label("allignment coef");
                ui.add(egui::Slider::new(
                    &mut run_options.allignment_coefficient,
                    0.0..=10.0,
                ))
            });

            ui.horizontal(|ui| {
                ui.label("cohesion coef");
                ui.add(egui::Slider::new(
                    &mut run_options.cohesion_coefficient,
                    0.0..=10.0,
                ))
            });

            ui.horizontal(|ui| {
                ui.label("separation coef");
                ui.add(egui::Slider::new(
                    &mut run_options.separation_coefficient,
                    0.0..=10.0,
                ))
            });

            ui.separator();

            ui.horizontal(|ui| {
                ui.label("wander distance");
                ui.add(egui::Slider::new(
                    &mut run_options.wander_distance,
                    0.5_f32..=30.,
                ))
            });
            ui.horizontal(|ui| {
                ui.label("wander radius");
                ui.add(egui::Slider::new(
                    &mut run_options.wander_radius,
                    0.5_f32..=30.,
                ))
            });

            ui.horizontal(|ui| {
                ui.label("sensory distance");
                ui.add(egui::Slider::new(
                    &mut run_options.sensory_distance,
                    5.0..=60.0,
                ))
            });

            ui.horizontal(|ui| {
                ui.label("allignment tresh coef");
                ui.add(egui::Slider::new(
                    &mut run_options.allignment_treshold_coefficient,
                    0.0..=2.0,
                ))
            });

            ui.horizontal(|ui| {
                ui.label("cohesion tresh coef");
                ui.add(egui::Slider::new(
                    &mut run_options.cohesion_treshold_coefficient,
                    0.0..=2.0,
                ))
            });

            ui.horizontal(|ui| {
                ui.label("separation tresh coef");
                ui.add(egui::Slider::new(
                    &mut run_options.separation_treshold_coefficient,
                    0.0..=2.0,
                ))
            });

            ui.horizontal(|ui| {
                ui.label("wander rate");
                ui.add(egui::Slider::new(
                    &mut run_options.wander_rate,
                    0.000001_f32..=1_f32,
                ))
            });

            ui.separator();

            ui.horizontal(|ui| {
                ui.label("size");
                ui.add(egui::Slider::new(&mut run_options.size, 1.0..=60.))
            });

            ui.horizontal(|ui| {
                ui.label("no neighbours");
                ui.add(egui::Slider::new(
                    &mut run_options.neighbours_cosidered,
                    0..=60,
                ))
            });

            ui.horizontal(|ui| {
                ui.add(egui::Checkbox::new(
                    &mut run_options.field_of_vision_on,
                    "FOV (deg)",
                ));

                ui.add_enabled(
                    run_options.field_of_vision_on,
                    egui::Slider::new(&mut run_options.field_of_vision_deg, 0.0..=360.),
                )
            });

            ui.horizontal(|ui| {
                ui.add(egui::Checkbox::new(&mut model.fps_limit_on, "FPS limit"));

                ui.add_enabled(
                    model.fps_limit_on,
                    egui::Slider::new(&mut model.fps, 5..=120),
                )
            });

            ui.horizontal(|ui| {
                ui.add(egui::Checkbox::new(&mut model.debug_grid, "Debug grid"));
                ui.add(egui::Checkbox::new(
                    &mut model.debug_distance,
                    "Debug distance",
                ));
            });

            ui.horizontal(|ui| {
                ui.add(egui::Checkbox::new(
                    &mut run_options.dbscan_flock_clustering_on,
                    "Cluster flocks",
                ));
            });

            ui.horizontal(|ui| {
                ui.label(format!("No. boids: {n:3.}", n = run_options.init_boids));
            });

            // let ticks = UPDATE_TICKS.load(Ordering::SeqCst);
            ui.horizontal(|ui| {
                ui.label(format!("No. ticks: {t:3.}", t = model.update_ticks));
            });

            ui.horizontal(|ui| {
                ui.label(format!(
                    "Time since last update: {t:2.2}",
                    t = model.since_last_update_micros as f64 / 1000.
                ));
            });
        });

    // update sensory distances
    run_options.update_sensory_distances();
    run_options.field_of_vision_half_rad = run_options.field_of_vision_deg * PI / 360.;

    // update model
    if model.control_state.execution_paused {
        return;
    }
    flock.update(&run_options);
    bird_watcher.watch(&flock);

    model.update_ticks += 1;
}

fn raw_window_event(_app: &App, model: &mut Model, event: &nannou::winit::event::WindowEvent) {
    model.egui.handle_raw_event(event);
}

fn record_update_timing(app: &App, model: &mut Model) {
    let micros_now = app.duration.since_start.as_micros();
    let since_last_update = micros_now - model.last_update_micros;
    model.since_last_update_micros = since_last_update;
    model.last_update_micros = micros_now;
}

static RENDER_TICKS: AtomicU64 = AtomicU64::new(0);
static UPDATE_TICKS: AtomicU64 = AtomicU64::new(0);

fn should_render_view(model: &Model) -> bool {
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

fn should_render_update(model: &mut Model) -> bool {
    let should_render =
        model.last_update_micros - model.last_render_micros > 1_000_000 / model.fps as u128;

    if should_render {
        _ = UPDATE_TICKS.fetch_add(1, Ordering::SeqCst);
        model.last_render_micros = model.last_update_micros;
    }

    should_render
}

// impl From<glam::f32::Vec2> for nannou::geom::Vec2 {

// }

fn key_released(_: &App, model: &mut Model, key: Key) -> () {
    if key == Key::S {
        model.run_options.seek_target_on = false;
        model.run_options.seek_location = None;
    }
}
fn key_pressed(app: &App, model: &mut Model, key: Key) -> () {
    let Model {
        ref mut run_options,
        ..
    } = model;

    if model.control_state.controls_open && // if controls are open
    // allow only these actions
    !(key == Key::C || key == Key::R || key == Key::Space)
    {
        return;
    }

    if key == Key::Space {
        // pause the whole simulation
        model.control_state.execution_paused = !model.control_state.execution_paused;
        app.window(model.window_id)
            .unwrap()
            .set_resizable(model.control_state.execution_paused);
    } else if key == Key::C {
        // show/hide controls
        model.control_state.controls_open = !model.control_state.controls_open
    } else if key == Key::F {
        // show/hide controls
        model.run_options.dbscan_flock_clustering_on = !model.run_options.dbscan_flock_clustering_on
    } else if key == Key::M {
        // show/hide controls
        model.run_options.stop_movement = !model.run_options.stop_movement
    } else if key == Key::Key1 {
        // turn alignment on/off
        run_options.alignment_on = !run_options.alignment_on;
    } else if key == Key::Key2 {
        // turn cohesion on/off
        run_options.cohesion_on = !run_options.cohesion_on;
    } else if key == Key::Key3 {
        // turn separation on/off
        run_options.separation_on = !run_options.separation_on;
    } else if key == Key::Key4 {
        // turn separation on/off
        run_options.wander_on = !run_options.wander_on;
    } else if key == Key::D {
        // delete a boid from the flock
        if run_options.init_boids == 1 {
            return;
        }

        let ids_delete: Vec<usize> =
            ((run_options.init_boids / 2)..run_options.init_boids).collect();
        run_options.init_boids /= 2;
        model.flock.delete_multiple(&ids_delete, run_options);
        // let to_delete = run_options.init_boids / 2;
        // for _ in 0..to_delete {
        //     let deleted_boid = model.flock.delete_last();
        //     match deleted_boid {
        //         Some(boid) => {
        //             if run_options.clicked_boid_id == boid.id {
        //                 run_options.clicked_boid_id = std::usize::MAX;
        //             }
        //             // run_options.init_boids -= 1;
        //         }
        //         None => {}
        //     }
        // }
        // run_options.init_boids /= 2;
    } else if key == Key::I {
        // insert a boid into the flock
        // for _ in 0..run_options.init_boids {
        //     model.flock.insert_single(&run_options);
        // }
        if run_options.init_boids > 2_usize.pow(15) {
            return;
        }
        run_options.init_boids *= 2;

        model.flock.insert_multiple(run_options);

        // update the information about number of boids in the system
    } else if key == Key::F1 {
        // alternative implementation of allignment
        run_options.allignment_impl_mode = !run_options.allignment_impl_mode;
    } else if key == Key::F2 {
        // alternative implementation of cohesion
        run_options.cohesion_impl_mode = !run_options.cohesion_impl_mode;
    } else if key == Key::F3 {
        // alternative implementation of separationk
        run_options.separation_impl_mode = !run_options.separation_impl_mode;
    } else if key == Key::F10 {
        // alternative implementation of separationk
        run_options.clustering_impl = !run_options.clustering_impl;
    } else if key == Key::F9 {
        // alternative implementation of separationk
        model.debug_labels = !model.debug_labels;
    } else if key == Key::F12 {
        // switch on or of coloring of no_boids in vision
        run_options.col_by_neighbour = !run_options.col_by_neighbour;
    } else if key == Key::R {
        // restart the flock, not the simulation
        model.flock.restart(&model.run_options);
        model.update_ticks = 0;
    } else if key == Key::V {
        // switch vision on or off
        run_options.field_of_vision_on = !run_options.field_of_vision_on
    } else if key == Key::S {
        // switch vision on or off
        run_options.seek_target_on = !run_options.seek_target_on;
        run_options.seek_location = Some(app.mouse.position());
    } else if key == Key::X {
        run_options.clicked_boid_id = usize::MAX;
    }
}

fn mouse_pressed(app: &App, model: &mut Model, _button: MouseButton) -> () {
    if model.control_state.controls_open {
        return;
    }
    let Model {
        ref mut run_options,
        ..
    } = model;

    let mouse_click = vec2(app.mouse.x, app.mouse.y);

    // let distance_boid = model
    // .flock
    // .view()
    // .0
    // .iter()
    // .filter(|b| {
    //     b.id == run_options.clicked_boid_id
    // })
    // .next();
    // match distance_boid {
    //     Some(db) => println!("{:?}", distance_dyn(db.position.x, mouse_click.x, db.position.y, mouse_click.y, run_options)),
    //     None => ()
    //     ,
    // }

    // println!("mouse clicked: {:?}", mouse_click);
    let clicked = model
        .flock
        .view()
        .0
        .iter()
        .filter(|b| {
            // get the distance from boid to click, if dist < size then boid is clicked
            let dist = b.position.distance(mouse_click);

            dist < run_options.size
        })
        .next();

    match clicked {
        Some(_) => {
            run_options.clicked_boid_id = clicked.unwrap().id;
        }
        None => {}
    }
}

fn window_closed(_app: &App, model: &mut Model) {
    let _ = model
        .bird_watcher
        .pop_data_save(&model.run_options.save_options);
}

pub trait Drawable {
    fn draw(&self, draw: &Draw, color: &Hsv, run_options: &RunOptions, model: &Model);
}
impl Drawable for Flock {
    fn draw(&self, draw: &Draw, color: &Hsv, run_options: &RunOptions, model: &Model) {
        for e in self.view2() {
            e.0.draw(&draw, color, run_options, e.1, model);
        }

        let settings = SpatHash1D::get_tracker_settings(run_options);

        if model.debug_distance {
            let mut dist_max: f32 =
                ((run_options.window.win_w as f32).pow(2.) + (run_options.window.win_h as f32).pow(2.)).sqrt();
            
            if model.run_options.distance == Distance::EucToroidal {
                dist_max /= 2.;
            }
            let target = model.flock.view2().filter(|tuple| tuple.0.id == 0).next();

            match target {
                Some(b_tuple) => {
                    for x in ((run_options.window.win_left as i32)
                        ..(run_options.window.win_right as i32))
                        .step_by(10)
                    {
                        for y in ((run_options.window.win_bottom as i32)
                            ..(run_options.window.win_top as i32))
                            .step_by(10)
                        {
                            let dist = distance_dyn(b_tuple.0.position.x, x as f32, b_tuple.0.position.y, y as f32, run_options);
                            let dist_scaled = dist / dist_max;
                            // dbg!(dist_scaled);
                            draw.ellipse().x_y(x as f32, y as f32).radius(5.).color(hsv(
                                dist_scaled,
                                1.,
                                0.5,
                            ));
                        }
                    }
                }
                None => (),
            }
        }

        if model.debug_grid {
            let max_ined = (settings.x_cell_count * settings.y_cell_count) as f32;
            for x in ((run_options.window.win_left as i32)..(run_options.window.win_right as i32))
                .step_by(10)
            {
                for y in ((run_options.window.win_bottom as i32)
                    ..(run_options.window.win_top as i32))
                    .step_by(10)
                {
                    let index = SpatHash1D::get_table_index(
                        x as f32,
                        y as f32,
                        run_options.window.win_left as f32,
                        (run_options.window.win_right - 1) as f32,
                        (run_options.window.win_bottom as f32) as f32,
                        (run_options.window.win_top - 1) as f32,
                        settings.x_cell_res,
                        settings.y_cell_res,
                        settings.x_cell_count as f32,
                    );

                    let r = 1. * ((index + 1) as f32 / max_ined) * (1. - (index % 2) as f32);
                    let g = 1. - 1. * ((index + 1) as f32 / max_ined) * (0. + (index % 2) as f32);
                    draw.ellipse()
                        .x_y(x as f32, y as f32)
                        .rgb(r, g, 1.)
                        .z(-1.)
                        .radius(2.);
                }
            }
            let x_count = settings.x_cell_count as i32;
            let y_count = settings.y_cell_count as i32;

            for column in 0..=(x_count) {
                let line_x_point =
                    run_options.window.win_left as f32 + column as f32 * settings.x_cell_res;

                draw.line()
                    .start(Vec2::new(line_x_point, run_options.window.win_top as f32))
                    .end(Vec2::new(line_x_point, run_options.window.win_bottom as f32))
                    .weight(2.)
                    .z(-1.)
                    .color(DEEPPINK);
            }

            for row in 0..=(y_count) {
                let line_y_point = run_options.window.win_bottom as f32 + row as f32 * settings.y_cell_res;

                draw.line()
                    .start(Vec2::new(run_options.window.win_left as f32, line_y_point))
                    .end(Vec2::new(run_options.window.win_right as f32, line_y_point))
                    .weight(2.)
                    .z(-1.)
                    .color(DEEPPINK);
            }
        }
    }
}

trait MyRotate {
    fn myrotate(&self, rhs: Vec2) -> Self;
}

impl MyRotate for Vec2 {
    #[inline]
    fn myrotate(&self, rhs: Vec2) -> Self {
        Vec2::new(
            self.x * rhs.x - self.y * rhs.y,
            self.y * rhs.x + self.x * rhs.y,
        )
    }
}

pub trait DrawableBoid {
    fn draw(
        &self,
        draw: &Draw,
        color: &Hsv,
        run_options: &RunOptions,
        metadata: &BoidMetadata,
        model: &Model,
    );
}
impl DrawableBoid for Boid {
    fn draw(
        &self,
        draw: &Draw,
        color: &Hsv,
        run_options: &RunOptions,
        metadata: &BoidMetadata,
        model: &Model,
    ) {
        let Boid {
            position,
            // pos_x,
            // pos_y,
            velocity,
            ..
        } = self;
        // draw a triangle rotated in the direction of velocity
        let theta = velocity.angle();
        // // Draw Boid body as an arrow with a triangle cutout
        let vertices = vec![
            pt2(-0.8 * run_options.size, 0.6 * run_options.size),
            pt2(run_options.size, 0.),
            pt2(-0.8 * run_options.size, -0.6 * run_options.size),
            pt2(-0.5 * run_options.size, 0.),
        ];
        let drawing = draw
            .polygon()
            .stroke(AZURE)
            // .stroke_weight(1.0)
            .points(vertices)
            .xy(*position)
            .z(1.)
            // investigate this:
            .rotate(2. * PI + theta);
        // let drawing = draw.ellipse()
        // .xy(*position)
        // .z(1.)
        // .radius(run_options.size);

        // If this instance is selected, show diagnostics
        if self.id == run_options.clicked_boid_id {
            if run_options.wander_on {
                // the current velocity vector normalized
                let heading = self.velocity.normalize();
                // gives the center of the circle driving the locomotion
                // let loco_center = self.position +  heading * (run_options.wander_radius * (2_f32).sqrt());
                let loco_center =
                    self.position + heading * (run_options.wander_distance * (2_f32).sqrt());

                // vector pointing at the point on circumference
                let wander_point = heading.myrotate(Vec2::new(
                    run_options.wander_radius * metadata.wander_direction.cos(),
                    run_options.wander_radius * metadata.wander_direction.sin(),
                ));

                // places the point onto the locomotion circle with respect to agent's location
                // let wander_f = loco_center + wander_point;

                // dbg!(wander_f);

                draw.ellipse()
                    .radius(run_options.wander_radius)
                    .color(rgba(0.1, 0.1, 0.1, 0.5))
                    .xy(loco_center)
                    .z(30.);

                draw.ellipse()
                    .radius(run_options.wander_radius * 0.2)
                    .color(WHITE)
                    .xy(loco_center
                        + wander_point.clamp_length(
                            run_options.wander_radius * 0.9,
                            run_options.wander_radius * 0.9,
                        ))
                    // .xy(center + Vec2::new(dx, dy).clamp_length(run_options.size * 0.9, run_options.size * 0.9))
                    .z(31.);
            }
            //            draw.line()
            //            .start(self.position)
            //            .end(self.velocity * 100. + self.position)
            //            .color(LIME)
            //            .weight(4.0);
            //    } else if metadata.clicke
            // println!("vel {:#?}", velocity);
            // println!("rad {:#?}", theta);
            // println!("deg {:#?}", theta * 180. / PI);

            // Highlight the boid
            drawing.color(HOTPINK);

            if run_options.field_of_vision_on {
                // drawing a semi-sircle to demonstrate FOV
                // https://stackoverflow.com/a/72701981/9316685
                let radius = run_options.separation_treshold_distance * 2.;
                let section = Ellipse::new(Rect::from_x_y_w_h(0., 0., radius, radius), 120.)
                    .section(
                        theta - run_options.field_of_vision_half_rad,
                        run_options.field_of_vision_half_rad * 2.,
                    );
                let triangles = section.trangles();
                let mut tris = Vec::new();
                for t in triangles {
                    tris.push(Tri([
                        [t.0[0][0], t.0[0][1], 0.],
                        [t.0[1][0], t.0[1][1], 0.],
                        [t.0[2][0], t.0[2][1], 0.],
                    ]))
                }
                draw.mesh()
                    .tris(tris.into_iter())
                    .xy(self.position)
                    .color(RED)
                    .rgba(255., 0., 0., 40.)
                    .z(-1.);
            }

            // show alignment radius
            if run_options.alignment_on {
                show_radius(
                    &draw,
                    *position,
                    run_options.sensory_distance * run_options.allignment_treshold_coefficient,
                    RED,
                );
            }

            // show cohesion radius
            if run_options.cohesion_on {
                show_radius(
                    &draw,
                    *position,
                    run_options.sensory_distance * run_options.cohesion_treshold_coefficient,
                    GREEN,
                );
            }

            // show  separation radius
            if run_options.separation_on {
                show_radius(
                    &draw,
                    *position,
                    run_options.sensory_distance * run_options.separation_treshold_coefficient,
                    BLUE,
                );
            }

            show_radius(&draw, *position, run_options.sensory_distance, YELLOW);

            // velocity vector scaled 100x
            draw.line()
                .start(self.position)
                .end(self.velocity * 50. + self.position)
                .color(LIME)
                .weight(4.0)
                .z(20.);
        } else if metadata.clicked_neighbour_id != std::usize::MAX
            && metadata.clicked_neighbour_id == run_options.clicked_boid_id
        {
            drawing.color(CHARTREUSE);
        } else {
            let saturation = if run_options.col_by_neighbour {
                color.saturation / 100. * 3. * metadata.n_neighbours as f32
            } else {
                color.saturation
            };

            // color.hue.

            // Take the color from UI settings
            drawing.hsv(
                if metadata.cluster_id != 0 {
                    ((10. * metadata.cluster_id as f32) % 360.) / 360.
                } else {
                    300. / 360.
                },
                saturation,
                if metadata.cluster_id != 0 {
                    1f32
                } else {
                    color.value
                },
            );
            // drawing.hsv(
            //     color.hue.to_positive_degrees() / 360.0,
            //     saturation,
            //     color.value * ( if self.cluster_id != 0 {1f32} else {0f32}),
            // );
        }

        if model.debug_labels {
            draw.text(&self.id.to_string())
                .color(RED)
                .xy(self.position + vec2( if self.position.x < 0. { 20.} else { -10.}, if self.position.y < 0. { 20. } else { -10.}));
    
            draw.text(&metadata.cluster_id.to_string())
                .color(BROWN)
                .xy(self.position + vec2(if self.position.x < 0. { 10.} else { -20.}, if self.position.y < 0. { 10. } else { -20.}));
        }
    }
}

fn show_radius<C>(draw: &Draw, xy: Vec2, r: f32, color: C)
where
    C: IntoLinSrgba<ColorScalar>,
{
    draw.ellipse()
        .xy(xy)
        .radius(r)
        .no_fill()
        .stroke(color)
        .stroke_weight(2.);
}

fn view(app: &App, model: &Model, frame: Frame) {
    if !should_render_view(model) {
        return;
    }
    let draw = app.draw();

    let settings = SpatHash1D::get_tracker_settings(&model.run_options);

    let index = SpatHash1D::get_table_index(
        app.mouse.x as f32,
        app.mouse.y as f32,
        model.run_options.window.win_left as f32,
        (model.run_options.window.win_right - 1) as f32,
        model.run_options.window.win_bottom as f32,
        (model.run_options.window.win_top - 1) as f32,
        settings.x_cell_res,
        settings.y_cell_res,
        settings.x_cell_count as f32,
    );

    let mouse_label2 = format!("index: {}", index);

    let mouse_label = format!("{x:.2} | {y:.2}", x = app.mouse.x, y = app.mouse.y);

    draw.text(&mouse_label)
        .x_y(app.mouse.x + 50., app.mouse.y + 20.)
        .z(10.)
        .color(BLACK)
        .font_size(20);
    draw.text(&mouse_label2)
        .x_y(app.mouse.x + 50., app.mouse.y - 20.)
        .z(10.)
        .color(BLACK)
        .font_size(20);

    let distance_boid = model
        .flock
        .view2()
        .filter(|(b, _)| {
            b.id == model.run_options.clicked_boid_id
        })
        .next();

    match distance_boid {
            Some((db, dm)) => { 
                let distance = distance_dyn(db.position.x, app.mouse.x, db.position.y, app.mouse.y, &model.run_options);
                let vec_to = tor_vec_p(db.position.x, app.mouse.x, db.position.y, app.mouse.y, &model.run_options.window);
                draw.text(&format!("cd: {:.2}", distance))
                    .x_y(0., 0.)
                    .z(10.)
                    .color(WHITE)
                    .font_size(20);
                draw.text(&format!("vec to: {:.2}, {:.2}", vec_to.0, vec_to.1))
                    .x_y(0., -20.)
                    .z(10.)
                    .color(WHITE)
                    .font_size(20);
                draw.text(&format!("flock id: {:.2}", dm.cluster_id))
                    .x_y(0., -40.)
                    .z(10.)
                    .color(WHITE)
                    .font_size(20);
            },
            None => (),
        };

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

    model
        .flock
        .draw(&draw, &model.color, &model.run_options, model);

    // let c = hsv(0.5, 1., 1.);

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

// fn draw_animate(app: &App, draw: &Draw, offset: f32) -> (f32, f32) {
//     // Generate sine wave data based on the time of the app
//     let my_time = app.time + offset;

//     let sine = (my_time * 1.0).sin();
//     let slowersine = (my_time * 2.0).sin();

//     // Get boundary of the window (to constrain the movements of our circle)
//     let boundary = app.window_rect();
//     let radius = 25.0;

//     // Map the sine wave functions to ranges between the boundaries of the window
//     let x = map_range(
//         sine,
//         -1.0,
//         1.0,
//         boundary.left() + radius,
//         boundary.right() - radius,
//     );
//     let y = map_range(
//         slowersine,
//         -1.0,
//         1.0,
//         boundary.bottom() + radius,
//         boundary.top() - radius,
//     );

//     if x > 0. {
//         draw.ellipse().radius(radius).color(GREEN).x_y(0.0, 0.0);
//     } else {
//         draw.ellipse().radius(radius).color(GREEN).x_y(x, y);
//     }

//     // Clear the background to purple.
//     // draw.background().color(PLUM);

//     // Draw a blue ellipse at the x/y coordinates 0.0, 0.0
//     // draw.ellipse().radius(radius).color(GREEN).x_y(x, y);

//     (x, y)
// }
