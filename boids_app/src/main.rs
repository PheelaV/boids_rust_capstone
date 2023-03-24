// This is the one and only bit that drives the user interface, it has snowballed into a huge mass
// Todo: breakup into smaller pieces

extern crate nannou;
use boids_lib::{
    birdwatcher::Birdwatcher, 
    boid::{Boid, BoidMetadata}, 
    flock::{spathash_tracker::SpatHash1D, Flock, tracker::TrackerSignal, replay_tracker::ReplayTracker}, 
    math_helpers::{distance_dyn, tor_vec_p, tor_vec}, 
    options::{
        self, Boundary, Distance, NoiseModel, RunOptions, SaveOptions, WindowSize, TrackerType,
    }};
use circular_queue::CircularQueue;
use clap_serde_derive::{clap::Parser, ClapSerde};
use itertools::Itertools;
use nannou::{draw::{properties::ColorScalar}, geom::{Ellipse, Tri}, rand::{Rng, self}, color::*, prelude::*};
use nannou_egui::{egui, Egui};
use splines::{Interpolation, Spline, Key as sKey};
use std::{path::Path, fs, process::Command, sync::atomic::{AtomicU64, Ordering}};

mod cliargs;
use cliargs::{Args, Config};

fn main() {
    nannou::app(model)
        // https://github.com/nannou-org/nannou/issues/708 this won't work at the moment and a custom sollution had to be used instead
        // .loop_mode(LoopMode::NTimes { number_of_updates: 60 })
        .update(update)
        .run();
}

struct ControlsState {
    execution_paused: bool,
    controls_open: bool,
}

pub struct Model<'a> {
    egui: Egui,
    color: Hsv,
    flock: Flock<'a>,
    run_options: RunOptions,
    last_update_micros: u128,
    last_render_micros: u128,
    update_ticks: u64,
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
    debug_diagnostics: bool,
    repulsive_distance: f32,
    repulsive_force: f32,
    draw_video_frames: bool,
    ghost_buffer: Option<CircularQueue<(Boid, BoidMetadata)>>,
    ghost_buffer_ready: bool,
    ghost_number: usize,
    ghost_mode_on: bool,
}

fn model<'a>(app: &App) -> Model<'a> {
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
    run_options.alignment_coefficient = config.alignment_coefficient;
    run_options.cohesion_coefficient = config.cohesion_coefficient;
    run_options.separation_coefficient = config.separation_coefficient;
    run_options.alignment_treshold_coefficient = config.alignment_treshold_coefficient;
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
    run_options.rules_impl = config.rules_impl;

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

    let draw_video_frames = false;

    if draw_video_frames {
        empty_recording_folder("my_local_video_png".to_owned());
    }

    Model {
        egui: Egui::from_window(&window),
        color: Hsv::new(RgbHue::from_degrees(20.), 1., 0.0),
        flock: Flock::new(&run_options),
        last_update_micros: 0,
        last_render_micros: 0,
        update_ticks: 0,
        since_last_update_micros: 0,
        fps: 120,
        fps_limit_on: match run_options.tracker_type {
            TrackerType::Replay(_, _) => {
                true
            },
            _ => false
        },
        run_options,
        update_timing_on: true,
        control_state: ControlsState {
            execution_paused: false,
            controls_open: false,
        },
        bird_watcher: birdwatcher,
        window_id: main_window,
        debug_grid: false,
        debug_distance: false,
        debug_labels: false,
        repulsive_distance: 100.,
        repulsive_force: 0.05,
        debug_diagnostics: false,
        draw_video_frames,
        ghost_buffer: None,
        ghost_buffer_ready: false,
        ghost_mode_on: false,
        ghost_number: 16,
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
        ref mut update_ticks,
        ..
    } = *model;

    // update window size as it could be resized
    let win = &app.window_rect();

    let new_win = {
        let win_left = win.left() as i32;
        let win_right = win.right() as i32;
        let wind_top = win.top() as i32;
        let win_bottom = win.bottom() as i32;
        let win_h = win.h() as i32;
        let win_w = win.w() as i32;
        WindowSize {
            win_left: win_left,
            win_right: win_right,
            win_top: wind_top,
            win_bottom: win_bottom,
            win_h: win_h,
            win_w: win_w,
        }
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
            match run_options.tracker_type {
                TrackerType::Replay(_, ticks_max) => {
                    ui.horizontal(|ui| {
                        if ui.button("⏴").clicked() {
                            flock.tracker.signal(TrackerSignal::ReplayTickBackward)
                        } 
                        if ui.button("⏵").clicked() {
                            flock.tracker.signal(TrackerSignal::ReplayTickForward)
                        }
                        if ui.button("⏯").clicked() {
                            flock.tracker.signal(TrackerSignal::ReplayPlayPause)
                        } 
                        if ui.button("⟲⟳").clicked() {
                            flock.tracker.signal(TrackerSignal::ReplaySwitchDirection)
                        } 
                        ui.add(egui::Slider::new(update_ticks, 0..=ticks_max));
                    });

                    ui.horizontal(|ui| {
                        ui.add(egui::Checkbox::new(&mut model.fps_limit_on, "FPS limit"));
        
                        ui.add_enabled(
                            model.fps_limit_on,
                            egui::Slider::new(&mut model.fps, 5..=120),
                        )
                    });
                    ui.horizontal(|ui| {
                        ui.horizontal(|ui| {
                            ui.label("size");
                            ui.add(egui::Slider::new(&mut run_options.size, 1.0..=160.))
                        });
                        ui.horizontal(|ui| {
                            ui.label("color");
                            edit_hsv(ui, color);
                        });
                    });
                    return ();
                },
                _ => (),
                // TrackerType::SpatHash => todo!(),
                // TrackerType::Naive => todo!(),
            }
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
                            "EucToroidal",
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
                            Boundary::Repulsive {
                                distance: *repulsive_distance,
                                force: *repulsive_force,
                            },
                            "Repulsive",
                        );
                    });
            });

            if let Boundary::Repulsive {
                distance: _,
                force: _,
            } = run_options.boundary
            {
                ui.horizontal(|ui| {
                    ui.label("rep distance");
                    ui.add(egui::Slider::new(repulsive_distance, 0.0..=300_f32));
                });
                ui.horizontal(|ui| {
                    ui.label("rep force");
                    ui.add(egui::Slider::new(repulsive_force, 0.001..=0.2));
                    if ui.button("apply").clicked() {
                        run_options.boundary = Boundary::Repulsive {
                            distance: *repulsive_distance,
                            force: *repulsive_force,
                        }
                    }
                });
            }

            ui.horizontal(|ui| {
                egui::ComboBox::from_label("noise")
                    .selected_text(format!("{:?}", run_options.noise_model))
                    .show_ui(ui, |ui| {
                        ui.selectable_value(
                            &mut run_options.noise_model,
                            NoiseModel::Reynolds,
                            "Reynolds",
                        );
                        ui.selectable_value(
                            &mut run_options.noise_model,
                            NoiseModel::Vicsek,
                            "Vicsek",
                        );
                    });
            });

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
                    "alignment",
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
                ui.label("wander rate");
                ui.add(egui::Slider::new(
                    &mut run_options.wander_rate,
                    0.000001_f32..=1_f32,
                ))
            });

            ui.horizontal(|ui| {
                ui.label("sensory distance");
                ui.add(egui::Slider::new(
                    &mut run_options.sensory_distance,
                    5.0..=100.0,
                ))
            });

            ui.horizontal(|ui| {
                ui.label("alignment c");
                ui.add(egui::Slider::new(
                    &mut run_options.alignment_coefficient,
                    0.0..=10.0,
                ));
                ui.label("alignment tc");
                ui.add(egui::Slider::new(
                    &mut run_options.alignment_treshold_coefficient,
                    0.0..=2.0,
                ));
            });

            ui.horizontal(|ui| {
                ui.label("cohesion c");
                ui.add(egui::Slider::new(
                    &mut run_options.cohesion_coefficient,
                    0.0..=10.0,
                ));
                ui.label("cohesion tc");
                ui.add(egui::Slider::new(
                    &mut run_options.cohesion_treshold_coefficient,
                    0.0..=2.0,
                ));
            });

            ui.horizontal(|ui| {
                ui.label("separation c");
                ui.add(egui::Slider::new(
                    &mut run_options.separation_coefficient,
                    0.0..=10.0,
                ));
                ui.label("separation tc");
                ui.add(egui::Slider::new(
                    &mut run_options.separation_treshold_coefficient,
                    0.0..=2.0,
                ))
            });

            ui.separator();

            ui.horizontal(|ui| {
                ui.label("size");
                ui.add(egui::Slider::new(&mut run_options.size, 1.0..=160.))
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
                ui.add(egui::Checkbox::new(
                    &mut run_options.rules_impl,
                    "rules impl",
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
                ui.label(format!("No. ticks: {t:3.}", t = update_ticks));
            });

            ui.horizontal(|ui| {
                ui.label(format!(
                    "Time since last update: {t:2.2}",
                    t = model.since_last_update_micros as f64 / 1000.
                ));
            });
        });

    // update model
    if model.control_state.execution_paused {
        return;
    }
    flock.update(run_options);

    if model.ghost_mode_on {
        if let Some(gb) = &model.ghost_buffer {
            model.ghost_buffer_ready = gb.len() == (model.ghost_number + 1) // + 1 for current agents
            * run_options.init_boids 
        } else {
            let new_buff = CircularQueue::with_capacity(
                (model.ghost_number + 1) // + 1 for current agents
                * run_options.init_boids
            );
            model.ghost_buffer = Some(new_buff);
        }

        flock.view2()
            .sorted_by(|a, b| a.0.id.cmp(&b.0.id))
            .for_each(|em| 
                _ = model.ghost_buffer.as_mut()
                    .expect("ghost buffer uninitialized")
                    .push((*em.0, *em.1))
            );
    } else {
        if let Some(_) = model.ghost_buffer {
            model.ghost_buffer = None;
        }
        model.ghost_buffer_ready = false;
    }
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
        ref mut flock,
        ..
    } = model;

    // limit keys the app react's to when running in replay mode
    if let TrackerType::Replay(_, _) = run_options.tracker_type {
        if !(key == Key::Space 
            || key == Key::R // restart replay
            || key == Key::M // stop movement
            || key == Key::C // controls toggle
            || key == Key::X // deselect clicked agent
            || key == Key::Escape // escape app
            || key == Key::Return // switch direction
            || key == Key::Right // step forward
            || key == Key::Left // step backward
            || key == Key::F // clustering toggle
            || key == Key::F8 // diag
            || key == Key::F9 // labels
            || key == Key::F10 // debug distance
            || key == Key::F11 // debug grid
            || key == Key::F12 // density based colours
            ) {
            return;
        }
        // commands specific to replay tracker
        if key == Key::Right {
            flock.tracker.signal(TrackerSignal::ReplayTickForward);
        } else if key == Key::Left {
            flock.tracker.signal(TrackerSignal::ReplayTickBackward);
        } else if key == Key::Return {
            flock.tracker.signal(TrackerSignal::ReplaySwitchDirection)
        }
    }

    // if controls are open limit the keys only to the nescessary ones
    if model.control_state.controls_open &&
    // allow only these actions
    !(key == Key::C 
        || key == Key::R 
        || key == Key::F8
        || key == Key::F9 
        || key == Key::F10 
        || key == Key::F11 
        || key == Key::F12 
        || key == Key::Space)
    {
        return;
    }

    if key == Key::Space {
        // pause the whole simulation (movement and updates)
        // this is a requirement for resizing (because spathash will dislike you otherwise)
        model.control_state.execution_paused = !model.control_state.execution_paused;
        app.window(model.window_id)
            .unwrap()
            .set_resizable(model.control_state.execution_paused);
    } else if key == Key::C {
        // toggle controls
        model.control_state.controls_open = !model.control_state.controls_open
    } else if key == Key::F {
        // toggle flock clustering
        model.run_options.dbscan_flock_clustering_on = !model.run_options.dbscan_flock_clustering_on
    } else if key == Key::G {
        // toggle ghost mode
        model.ghost_mode_on = !model.ghost_mode_on;
    } else if key == Key::M {
        // stops movement, does not stop updates
        model.run_options.stop_movement = !model.run_options.stop_movement
    } else if key == Key::Key1 {
        // toggle alignment
        run_options.alignment_on = !run_options.alignment_on;
    } else if key == Key::Key2 {
        // toggle cohesion 
        run_options.cohesion_on = !run_options.cohesion_on;
    } else if key == Key::Key3 {
        // toggle separation 
        run_options.separation_on = !run_options.separation_on;
    } else if key == Key::Key4 {
        // toggle wander 
        run_options.wander_on = !run_options.wander_on;
    } else if key == Key::D && !model.ghost_mode_on {
        // delete a boid from the flock
        if run_options.init_boids == 1 {
            return;
        }
        let ids_delete: Vec<usize> =
            ((run_options.init_boids / 2)..run_options.init_boids).collect();
        run_options.init_boids /= 2;
        model.flock.delete_multiple(&ids_delete, run_options);
    } else if key == Key::I && !model.ghost_mode_on {
        if run_options.init_boids > 2_usize.pow(15) { // M1 Pro gets up to ~6 fps at 2^15 units
            return;
        }
        run_options.init_boids *= 2;

        model.flock.insert_multiple(run_options);
    } else if key == Key::F1 {
        // alternative implementation of alignment
        run_options.alignment_impl_mode = !run_options.alignment_impl_mode;
    } else if key == Key::F2 {
        // alternative implementation of cohesion
        run_options.cohesion_impl_mode = !run_options.cohesion_impl_mode;
    } else if key == Key::F3 {
         // alternative implementation of separation
        run_options.separation_impl_mode = !run_options.separation_impl_mode;
    } else if key == Key::F5 {
        // alternative implementation of run_rules
        run_options.rules_impl = !run_options.rules_impl;
    } else if key == Key::F8 {
        // show additional diagnostics
        model.debug_diagnostics = !model.debug_diagnostics;
    } else if key == Key::F9 {
        // show boid labels (id and flock id)
        model.debug_labels = !model.debug_labels;
    } else if key == Key::F10 {
        // show distance debugging interface
        model.debug_distance = !model.debug_distance;
    } else if key == Key::F11 {
        // show grid debugging interface
        model.debug_grid = !model.debug_grid;
    } else if key == Key::F12 {
        // switch on or of coloring of no_boids in vision
        run_options.col_by_neighbour = !run_options.col_by_neighbour;
    } else if key == Key::R {
        // restart the flock, not the simulation
        model.flock.restart(&model.run_options);
        model.update_ticks = 0;
        model.bird_watcher.restart();

        if model.ghost_mode_on {
            model.ghost_buffer.as_mut().expect("Ghost buffer should be present").clear();
            model.ghost_buffer_ready = false;
        }
        // TODO: under test, creating of video
        if model.draw_video_frames {
            let video_path = Path::new("/my_local_video_png");

            if video_path.exists() {
                fs::remove_dir_all("/my_local_video_png").expect("video folder deletion error");
                fs::create_dir("/my_local_video_png").expect("video folder creation error");
            }
        }

    } else if key == Key::V {
        // switch vision on or off
        run_options.field_of_vision_on = !run_options.field_of_vision_on
    } else if key == Key::S {
        // seek a target
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
    let clicked = model
        .flock
        .view2()
        .filter(|(e, _)| {
            // get the distance from boid to click, if dist < size then boid is clicked
            let dist = e.position.distance(mouse_click);

            dist < run_options.size
        })
        .map(|(e, _)| e)
        .next();

    match clicked {
        Some(_) => {
            run_options.clicked_boid_id = clicked.unwrap().id;
        }
        None => {}
    }
}

fn window_closed(_app: &App, model: &mut Model) {
    if model.draw_video_frames {
        // execute ffmpeg script
        let script_path = Path::new("my_local_helpers/create_video.sh");

        if !script_path.exists() {
            println!("draw_video_frames is on but video-compilation script has not been found")
        } else {
            let output = Command::new("my_local_helpers/create_video.sh")
            .spawn()
            .expect("Failed to call script")
            .wait()
            .expect("video script execution failed");
            
            if output.success() {
                println!("video created");
                empty_recording_folder("my_local_video_png".to_owned())
                
            } else {
                println!("video creation failed code: {}", output.code().unwrap())
            }
        }
    }

    let _ = model
        .bird_watcher
        .pop_data_save(&model.run_options.save_options);
}

fn empty_recording_folder(path: String) {

    let video_path = Path::new(&path);
    
    if video_path.exists() {
        fs::remove_dir_all(video_path).expect("video folder deletion error");
        fs::create_dir(video_path).expect("video folder creation error");
    }
}

fn draw_spline_path(draw: &Draw, spline: Vec<Vec2>) {
// <'a, I>
// where
//     I: Iterator<Item = &'a Vec2>,
    let spline_points = spline.len() as f32;

    let keys_iter = spline
        .iter()
        .enumerate()
        .map(|(i, p)| {
            sKey::new(i as f32 / spline_points, *p, Interpolation::CatmullRom)
        })
        .collect_vec();

    // the values above have been mapped to [0, ..., 1], as corner points are consumed by the
    // spline algorithm we can assess the range [1/spline_points, 1 - 1/spline_points]
    let spline = Spline::from_vec(keys_iter);
    let spline_parts  = 100f32;

    let spline_points = (1..(spline_parts as i32 - 1))
        .into_iter()
        .map(|i| spline.clamped_sample( i as f32 / (spline_parts)))
        .filter(|sp| sp.is_some())
        .map(|sp| sp.unwrap())
        .collect_vec();

    draw.path()
        .stroke()
        .weight(2.)     
        .color(BLACK)
        .points(spline_points);
}

 
fn crosses_boundary(v1: Vec2, v2: Vec2, run_options: &RunOptions) -> bool {
    let tor_vec = tor_vec( v1, v2, &run_options.window);
    let euc_vec = v2 - v1;

    (tor_vec - euc_vec).length_squared() > 0.00025 // five thousands
}

/// Takes in a sequence of spline segments, inclusive numerical ranges. 
/// Makes sure all the segments are at least 4 elements in length 
/// (minimum for Catmull-Rom spline).
/// Uses time notation to refer to previous and upcoming segments,
/// past and future respectively. Maybe a little counter-intuitive since
/// when going forward in time, the past elements have a higher index.
fn make_spline_sub_segments(sub_segments: &Vec<(usize, usize)>) 
    -> Vec<(usize, usize, Option<Vec<usize>>, Option<Vec<usize>>)> {

    // each one of these is the (start of spline, end of spline, points added to the start, points added to the end)
    let mut result: Vec<(usize, usize, Option<Vec<usize>>, Option<Vec<usize>>)> = Vec::new();

    for i in 0..sub_segments.len() {
        // keep track of current spline segment
        let (mut c_start, mut c_end) = sub_segments[i];
        let mut p_i = i + 1;
        let mut c_len: usize; // intentionally set to zero to make it enter condition below

        // here we always have past
        if c_start == 0 && i == 0 {
            let append_end;

            // begin by just joining up the past
            c_len = sub_segments[p_i].0 + 1;

            if c_len < 4 {
                // we can't do this straight away as we need to 
                // join the past and we did not know whether cur has >4 elements
                c_len = 4;
                append_end = Some(
                    ((c_end + 1)..c_len)
                        .into_iter()
                        .collect_vec()
                )
                
            } else {
                append_end = Some(vec![sub_segments[p_i].0])
            }
            result.push((0, c_len - 1, None, append_end));
        } else {
            let mut f_i = i - 1;
            // check whether we need to join up the past
            // future is implicit since we know we would not be joining
            // spline segments if there was not at least two
            let past_exists = i < sub_segments.len() - 1;
            
            let mut append_start = Vec::with_capacity(4);
            let mut append_end = Vec::with_capacity(if past_exists { 4 } else { 0});

            // stretch for future
            c_start = sub_segments[f_i].1;
            append_start.push(c_start);
            // stretch for past
            if past_exists { c_end = sub_segments[p_i].0;}
            append_end.push(c_end);
            // now that we are connected to potentially both sides
            // check length
            c_len = c_end - c_start + 1;

            let mut future_considered = false;
            let mut past_considered = false;

            while c_len < 4 {
                // have we tried expanding into the future?
                if !future_considered {
                    // try expanding future
                    future_considered = true;
                    if c_start > sub_segments[f_i].0 { 
                        c_start -= 1;
                        append_start.push(c_start);
                    }
                } else if past_exists && !past_considered {
                    past_considered = true;
                    if c_end < sub_segments[p_i].1 {
                        c_end += 1;
                        append_end.push(c_end);
                    }
                }
                c_len = c_end - c_start + 1;

                // if we have tried both future and past and are still under length
                // just check whether we can expand further or panic!
                if c_len < 4 // are we under length
                && future_considered && (past_considered || !past_exists) // have we considered both past (if exists) and future?
                {
                    if c_start != sub_segments[f_i].0 // (all negated) -> if we've hit the future end
                        || !(!past_exists || c_end == sub_segments[p_i].1) { // and past does not exist or we've hit its end
                            // try again
                            future_considered = false;
                            past_considered = false;
                            continue;
                    }

                    // or consider expansion
                    let mut lets_panic = true;
                    if f_i > 0 { 
                        f_i -= 1;
                        lets_panic = false;
                    }
                    if p_i < sub_segments.len() - 1 { 
                        p_i += 1; 
                        lets_panic = false;
                    }

                    if lets_panic {
                        panic!("Arrr, hit a dead end here mate.");
                    }
                }
            }

            if append_start.len() < 2 {
                if c_start != 0 {
                    c_start -= 1;
                    append_start.push(append_start.last().unwrap() - 1)
                    // we skip this one as it has already been take care of 
                    // continue; // maybe?
                }
            }

            if past_exists && append_end.len() < 3 {
                if sub_segments[p_i].1 > c_end {
                    c_end += 1;
                    append_end.push(c_end)
                }
            }

            append_start.reverse();

            result.push((c_start, c_end, Some(append_start), 
            if past_exists { Some(append_end) } else { None }));
        }
    }

    result
} 

/// Takes a point and projects it to its nearest representation according to the toroidal distance
fn project_this_to_that(this: Vec2, that: Vec2, run_options: &RunOptions) -> Vec2{
    let WindowSize {
        win_w,
        win_h,
        win_right,
        win_top,
        ..
    } = run_options.window;

    let x = {
        let x_a = this.x + win_w as f32 - that.x;
        let x_b = this.x - win_w as f32 - that.x;
        if x_a.abs() > win_right as f32 && x_b.abs() > win_right as f32 {
            // there is nothing to project in this dimension
            this.x
        } else if x_a.abs() < (x_b).abs() {
            x_a + that.x
        } else {
            x_b + that.x
        }
    };

    let y = {
        let y_a = this.y + win_h as f32 - that.y;
        let y_b = this.y - win_h as f32 - that.y;
        if y_a.abs() > win_top as f32 && y_b.abs() > win_top as f32 {
            // there is nothing to project in this dimension
            this.y
        } else if (y_a).abs() < (y_b).abs() {
            y_a + that.y
        } else {
            y_b + that.y
        }
    };

    return vec2(x, y);

}

// ghosts are representations of a given entity in the past (positive ghosts, negative would be future but that is not suppoerted atm)
fn draw_ghosts(entity: &Boid, _: &BoidMetadata, ghosts: Vec<&(&Boid, &BoidMetadata)>,
 draw: &Draw, run_options: &RunOptions) {
        // CatmullRom needs at least 4 points to draw a single line
        // it uses the first and last points as guides on how to angle the last (or only) segment(s)

        // accumulate all the spline points
        let mut spline_points: Vec<Vec2> = Vec::new();

        // add projected destination given current velocity, that is to extend
        // the spline to reach to the current position
        spline_points.push(entity.position + entity.velocity.normalize() 
            * (tor_vec(entity.position, ghosts[0].0.position
                ,&run_options.window)).length());
        spline_points.push(entity.position);
        // add all the ghosts
        spline_points.extend(ghosts.iter().map(|be| be.0.position));
        // add a tail extending from the last ghost
        spline_points.push(ghosts[ghosts.len() - 1].0.position - ghosts[ghosts.len() - 1].0.velocity);

        // detect periodic boundary conditions
        let mut spline_pieces: Vec<(usize, usize)> = Vec::new();

        let mut temp_end: usize = 0;
        for i in 0..(spline_points.len() - 1) {
            // get vector pointing to the previous step
            // check if it crosses boundary
            if crosses_boundary(spline_points[i], spline_points[i + 1], run_options) {
                if i == temp_end {
                    // only one point made it over
                    spline_pieces.push((i,i));
                } else {
                    spline_pieces.push(( temp_end, i));
                    temp_end = i + 1;
                }
            }
        }
        if temp_end == spline_points.len() {
            // only one point made it over
            spline_pieces.push((temp_end, temp_end));
        } else {
            spline_pieces.push((temp_end, spline_points.len() - 1));
        }
        
        let spline_segments = if spline_pieces.len() > 1 { 
            make_spline_sub_segments(&spline_pieces)
        } else {
            spline_pieces
                .iter()
                .map(|x| 
                    (x.0, x.1, None::<Vec<usize>>, None::<Vec<usize>>)
                )
                .collect_vec()
        };

        // map each segment into an iterator
        spline_segments
            .iter()
            .map(|s|{
                let mut sp = (s.0..=s.1)
                .into_iter()
                .map(|i| spline_points[i])
                .collect_vec();
            
                // in care there are appended points, we need to project them nearby
                let append_start = &s.2;
                if let Some(ap_s) = append_start {
                    for i in (0..ap_s.len()).rev() {
                            sp[i] = project_this_to_that(sp[i], sp[i + 1]
                                ,run_options)
                        }
                    }
                let append_end = &s.3;
                if let Some(ap_e) = append_end {
                    let sp_last = sp.len() - 1;
                    for i in (0..ap_e.len()).rev() {
                        sp[sp_last - i] = 
                        project_this_to_that(sp[sp_last - i]
                            ,sp[sp_last - i - 1], run_options)
                    }
                }

                sp
            })
            .for_each(|spline| {
                draw_spline_path(draw, spline)
            });

}

pub trait Drawable {
    fn draw(&self, draw: &Draw, color: &Hsv, run_options: &RunOptions, model: &Model);
}

impl<'a> Drawable for Flock<'a> {
    fn draw(&self, draw: &Draw, color: &Hsv, run_options: &RunOptions, model: &Model) {
        match run_options.tracker_type {
            // if it is a replay tracker, enable ghosting
            TrackerType::Replay(_, _) => {
                if let Some(rt) = model.flock.tracker.as_any().downcast_ref::<ReplayTracker>() {
                    let view_mat = self.view2().collect_vec();

                    (0..rt.no_boids)
                        .into_iter()
                        .for_each(|i| {
                            let (e, m) = view_mat[i]; 
                                let ghosts = view_mat.iter()
                                .skip(rt.no_boids + e.id)
                                .step_by(rt.no_boids)
                                .take(rt.no_ghosts.abs() as usize)
                                .collect_vec();
                            
                                draw_ghosts(e, m, ghosts, draw, &rt.run_options); 
                        });
                };
            },
            // for all other trackers proceed as standard
            _ => {

                if model.ghost_mode_on && model.ghost_buffer_ready {
                    // draw agents and ghosts
                    if let Some(gb) = &model.ghost_buffer {
                        gb.iter()
                            .take(run_options.init_boids)
                            .enumerate()
                            .for_each(|(i, em)| {
                                let ghosts = gb.iter()
                                    .skip(i)
                                    .step_by(run_options.init_boids)
                                    .map(|em| (&em.0, &em.1))
                                    .collect_vec();

                                draw_ghosts(&em.0, &em.1, ghosts.iter().collect_vec(), draw, run_options);
                                em.0.draw(draw, color, run_options, &em.1, model);
                            })

                    } else {
                        panic!("buffer should have been created in update")
                    }
                } else {
                    // just draw each boid
                    self.view2().for_each(|(e, m)|  e.draw(&draw, color, run_options, m, model))
                }
            }
        }

        let settings = SpatHash1D::get_tracker_settings(run_options);

        if model.debug_distance {
            let mut dist_max: f32 = ((run_options.window.win_w as f32).pow(2.)
                + (run_options.window.win_h as f32).pow(2.))
            .sqrt();

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
                            let dist = distance_dyn(
                                b_tuple.0.position.x,
                                x as f32,
                                b_tuple.0.position.y,
                                y as f32,
                                run_options,
                            );
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
                    .end(Vec2::new(
                        line_x_point,
                        run_options.window.win_bottom as f32,
                    ))
                    .weight(2.)
                    .z(-1.)
                    .color(DEEPPINK);
            }

            for row in 0..=(y_count) {
                let line_y_point =
                    run_options.window.win_bottom as f32 + row as f32 * settings.y_cell_res;

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
            velocity,
            ..
        } = self;
        // Draw a triangle rotated in the direction of velocity
        // self.position += self.velocity * run_options.baseline_speed;

        // let mut theta = (self.position self.position - self.velocity * run_options.baseline_speed).angle();
        let mut theta = velocity.angle();

        if theta.is_nan() { 
            theta = rand::thread_rng().gen_range(-PI..PI);
        }
        // Draw Boid body as an arrow with a triangle cutout
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
        
        // If this instance is selected, show diagnostics
        if self.id == run_options.clicked_boid_id {
            if run_options.wander_on {
                // the current velocity vector normalized
                let heading = self.velocity.normalize();
                // gives the center of the circle driving the locomotion
                let loco_center =
                    self.position + heading * (run_options.wander_distance * (2_f32).sqrt());

                // vector pointing at the point on circumference
                let wander_point = heading.myrotate(Vec2::new(
                    run_options.wander_radius * metadata.wander_direction.cos(),
                    run_options.wander_radius * metadata.wander_direction.sin(),
                ));

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
                    .z(31.);
            }

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
                    run_options.sensory_distance * run_options.alignment_treshold_coefficient,
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

            // velocity vector scaled 25x
            draw.line()
                .start(self.position)
                .end(self.velocity * 25. + self.position)
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

            // Take the color from UI settings
            drawing.hsv(
                if run_options.dbscan_flock_clustering_on && metadata.cluster_id != 0 {
                    ((10. * metadata.cluster_id as f32) % 360.) / 360.
                } else {
                    color.hue.to_positive_degrees() / 360.
                },
                saturation,
                if metadata.cluster_id != 0 {
                    1f32
                } else {
                    color.value
                },
            );
        }

        if model.debug_labels {
            draw.text(&self.id.to_string()).color(RED).xy(self.position
                + vec2(
                    if self.position.x < 0. { 20. } else { -10. },
                    if self.position.y < 0. { 20. } else { -10. },
                ));

            draw.text(&metadata.cluster_id.to_string())
                .color(BROWN)
                .xy(self.position
                    + vec2(
                        if self.position.x < 0. { 10. } else { -20. },
                        if self.position.y < 0. { 10. } else { -20. },
                    ));
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

    if model.debug_diagnostics {
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
            .filter(|(b, _)| b.id == model.run_options.clicked_boid_id)
            .next();
    
        match distance_boid {
            Some((db, dm)) => {
                let distance = distance_dyn(
                    db.position.x,
                    app.mouse.x,
                    db.position.y,
                    app.mouse.y,
                    &model.run_options,
                );
                let vec_from = tor_vec_p(
                    app.mouse.x,
                    db.position.x,
                    app.mouse.y,
                    db.position.y,
                    &model.run_options.window,
                );
                let vec_to = tor_vec_p(
                    db.position.x,
                    app.mouse.x,
                    db.position.y,
                    app.mouse.y,
                    &model.run_options.window,
                );
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
                draw.text(&format!(
                        "cosine test: {:.2}",
                        app.mouse
                        .position()
                        .normalize()
                        .dot(Vec2::new(1., 1.).normalize())
                    ))
                    .x_y(0., -60.)
                    .z(10.)
                    .color(WHITE)
                    .font_size(20);
                draw.text(&format!("vec fr: {:.2}, {:.2}", vec_from.0, vec_from.1))
                    .x_y(0., -90.)
                    .z(10.)
                    .color(WHITE)
                    .font_size(20);
            }
            None => (),
        };
    }

    draw.background().color(hsv(0.5, 0., 0.4));

    model
        .flock
        .draw(&draw, &model.color, &model.run_options, model);

    draw.to_frame(app, &frame).unwrap();
    
    model.egui.draw_to_frame(&frame).unwrap();

    if model.draw_video_frames {
        let video_frame_path = format!("my_local_video_png/{}.png", model.update_ticks);

        app.window(model.window_id)
        .expect("main app window")
        .capture_frame(video_frame_path);
    }
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