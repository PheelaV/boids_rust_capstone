extern crate nannou;
use std::{
    sync::atomic::{
        AtomicU64, 
        Ordering
    },
    io::BufReader,
    fs::File
};
use clap_serde_derive::{
    clap::Parser,
    ClapSerde,
};
use nannou::draw::properties::ColorScalar;
use nannou::geom::{
    Tri, 
    Ellipse
};
use nannou::{
    color::*, 
    prelude::*
};
use nannou_egui::{
    egui, 
    Egui
};
use boids_lib::boid::Boid;
use boids_lib::flock::Flock;
use boids_lib::options::{
        Distance, 
        RunOptions, 
        SaveOptions, 
        self
};
use boids_lib::birdwatcher::Birdwatcher;

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

struct Model {
    egui: Egui,
    color: Hsv,
    flock: Flock,
    run_options: RunOptions,
    last_update_micros: u128,
    last_render_micros: u128,
    since_last_update_micros: u128,
    fps: u8,
    fps_limit_on: bool,
    update_timing_on: bool,
    control_state: ControlsState,
    bird_watcher: Birdwatcher
}

fn model(app: &App) -> Model {
        // let args = Args::parse();
        // Parse whole args with clap
    let mut args = Args::parse();

    // Get config file
    let config = if let Ok(f) = File::open(&args.config_path) {
        // Parse config with serde
        match serde_yaml::from_reader::<_, <Config as ClapSerde>::Opt>(BufReader::new(f)) {
            // merge config already parsed from clap
            Ok(config) => Config::from(config).merge(&mut args.config),
            Err(err) => panic!("Error in configuration file:\n{}", err),
        }
    } else {
        // If there is not config file return only config parsed from clap
        Config::from(&mut args.config)
    };

    // let s  = serde_yaml::to_string(&config);

    // println!("{}", s.unwrap());

    let mut run_options: RunOptions = Default::default();

    run_options.init_boids = config.no_boids;
    run_options.sample_rate = config.sample_rate;

    let save_options = SaveOptions {
        save_locations: config.save,
        save_locations_timestamp: config.save_timestamp,

        // default
        save_locations_path: run_options.save_options.save_locations_path
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

    let main_window = app
        .new_window()
        .key_pressed(key_pressed)
        .mouse_pressed(mouse_pressed)
        .closed(window_closed)
        .size(run_options.window.win_w as u32, run_options.window. win_h as u32)
        .title("boids flocking")
        .raw_event(raw_window_event)
        .view(view)
        .build()
        .unwrap();

    let window = app.window(main_window).unwrap();

    let birdwatcher = Birdwatcher::new(run_options.sample_rate);

    Model {
        egui: Egui::from_window(&window),
        color: Hsv::new(RgbHue::from_degrees(190.), 74., 32.),
        flock: Flock::new(&run_options),
        run_options,
        last_update_micros: 0,
        last_render_micros: 0,
        since_last_update_micros: 0,
        control_state: ControlsState { 
            execution_paused: false,
            controls_open: false,
        },
        fps: 120,
        fps_limit_on: false,
        update_timing_on: true,
        bird_watcher: birdwatcher,
    }
}

fn update(app: &App, model: &mut Model, update: Update) {

    // Record update timing if we are looking at upate times or limiting frame times 
    if model.fps_limit_on || model.update_timing_on {
        record_update_timing(app, model);

        // If we are limitting FPS, decide whether the update should run
        if model.fps_limit_on && !should_render_update(model) {
            return
        }
    }

    //TODO: might not be needed
    if app.keys.down.contains(&Key::Space) {
        return
    }

    let Model {
        ref mut egui,
        ref mut color,
        ref mut flock,
        ref mut run_options,
        ref mut bird_watcher,
        ..
    } = *model;

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
                ui.add(egui::Slider::new(&mut run_options.sensory_distance, 10.0..=400.0))
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
                ui.add(egui::Slider::new(&mut run_options.size, 1.0..=60.))
            });

            ui.horizontal(|ui| {
                ui.label("no neighbours");
                ui.add(egui::Slider::new(&mut run_options.neighbours_cosidered, 1..=60))
            });

            ui.horizontal(|ui| {
                ui.add(egui::Checkbox::new(&mut run_options.field_of_vision_on, "FOV (deg)"));

                ui.add_enabled(run_options.field_of_vision_on, egui::Slider::new(&mut run_options.field_of_vision_deg, 0.0..=360.))
            });

            ui.horizontal(|ui| {
                ui.add(egui::Checkbox::new(&mut model.fps_limit_on, "FPS limit"));

                ui.add_enabled(model.fps_limit_on, egui::Slider::new(&mut model.fps, 5..=120))
            });

            ui.horizontal(|ui| {
                ui.label(format!("No. boids: {n:3.}", n = run_options.init_boids));
            });

            let ticks = UPDATE_TICKS.load(Ordering::SeqCst);
            ui.horizontal(|ui| {
                ui.label(format!("No. ticks: {t:3.}", t = ticks));
            });

            ui.horizontal(|ui| {
                ui.label(format!("Time since last update: {t:2.2}", t = model.since_last_update_micros as f64 / 1000.));
            });
        });

    // update sensory distances
    run_options.update_sensory_distances();
    run_options.field_of_vision_half_rad = run_options.field_of_vision_deg * PI/360.;

    // update model
    if model.control_state.execution_paused {
        return
    }
    flock.update(&run_options);
     bird_watcher.watch(&flock);
}

fn raw_window_event(_app: &App, model: &mut Model, event: &nannou::winit::event::WindowEvent) {
    model.egui.handle_raw_event(event);
}

fn record_update_timing(app: &App, model: &mut Model){
    let micros_now = app.duration.since_start.as_micros();
    let since_last_update = micros_now - model.last_update_micros;
    model.since_last_update_micros = since_last_update;
    model.last_update_micros = micros_now;
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

fn should_render_update(model: &mut Model) -> bool {
    let should_render = model.last_update_micros - model.last_render_micros > 1_000_000 / model.fps as u128;

    if should_render {
        UPDATE_TICKS.fetch_add(1, Ordering::SeqCst);
        model.last_render_micros = model.last_update_micros;
    }

    should_render
}

fn key_pressed(_app: &App, model: &mut Model, key: Key) -> () {
    let Model {
        ref mut run_options,
        ..
    } = model;

    if key == Key::Space {  // pause the whole simulation
        model.control_state.execution_paused = !model.control_state.execution_paused
    } else if key == Key::C { // show/hide controls
        model.control_state.controls_open = !model.control_state.controls_open
    } else if key == Key::Key1 { // turn alignment on/off
        run_options.alignment_on = !run_options.alignment_on;
    } else if key == Key::Key2 { // turn cohesion on/off
        run_options.cohesion_on = !run_options.cohesion_on;
    } else if key == Key::Key3 { // turn separation on/off
        run_options.separation_on = !run_options.separation_on;
    } else if key == Key::D { // delete a boid from the flock
        if run_options.init_boids == 1 {return}

        let to_delete = run_options.init_boids / 2;
        for _ in 0..to_delete{

            let id_delete = model.flock.delete_last();
            match id_delete {
                Some(x) => {
                    if run_options.clicked_boid_id == x {
                        run_options.clicked_boid_id = std::u32::MAX;
                    }
                    // run_options.init_boids -= 1;
                },
                None => {},
            }
        }
        run_options.init_boids /= 2;
            
    } else if key == Key::I { // insert a boid into the flock
        // actually make the boid insertion
        for _ in 0..run_options.init_boids {
            model.flock.insert(&run_options);
        }
        // update the information about number of boids in the system
        // run_options.init_boids += 1;
        run_options.init_boids *= 2;
        
    } else if key == Key::F1 { // alternative implementation of allignment
        run_options.allignment_impl_mode = !run_options.allignment_impl_mode;
    } else if key == Key::F2 { // alternative implementation of cohesion
        run_options.cohesion_impl_mode = !run_options.cohesion_impl_mode;
    } else if key == Key::F3 { // alternative implementation of separationk
        run_options.separation_impl_mode = !run_options.separation_impl_mode;
    } else if key == Key::F12 { // switch on or of coloring of no_boids in vision 
        run_options.col_by_neighbour = !run_options.col_by_neighbour;
    } else if key == Key::R { // restart the flock, not the simulation
        model.flock.restart(&model.run_options);
    } else if key == Key::V { // switch vision on or off
        run_options.field_of_vision_on = !run_options.field_of_vision_on
    }
}

fn mouse_pressed(app: &App, model: &mut Model, _button: MouseButton) -> () {
    if model.control_state.controls_open {
        return
    }
    let Model {
        ref mut run_options,
        flock,
        ..
    } = model;

    let mouse_click = vec2(app.mouse.x, app.mouse.y);

    let clicked = model.flock.view().iter()
    .filter(|b|{
        // get the distance from boid to click, if dist < size then boid is clicked
        let dist = b.position.distance(mouse_click);

        dist < run_options.size
    })
    .next();

    match clicked {
        Some(_) => {
            run_options.clicked_boid_id =  clicked.unwrap().id;
        },
        None => {},
    }
}

fn window_closed(_app: &App, model: &mut Model) {
    let _ = model.bird_watcher.pop_data_save(&model.run_options.save_options);
}

pub trait Drawable {
    fn draw(&self, draw: &Draw, color: &Hsv, run_options: &RunOptions);
}
impl Drawable for Flock {
    fn draw(&self, draw: &Draw, color: &Hsv, run_options: &RunOptions) {
        for b in self.view().iter() {
            b.draw(&draw, color, run_options);
        }

        for b in self.tracker.view().iter() {
            b.draw(&draw, &nannou::color::hsv(0., 0., 0.), run_options);
        }
    }
}

impl Drawable for Boid {
    fn draw(&self, draw: &Draw, color: &Hsv, run_uptions: &RunOptions) {
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
            pt2(-0.8 * run_uptions.size, 0.6 * run_uptions.size),
            pt2(run_uptions.size, 0.),
            pt2(-0.8 * run_uptions.size, -0.6 * run_uptions.size),
            pt2(-0.5 * run_uptions.size, 0.),
        ];
        let drawing = draw.polygon()
            .stroke(AZURE)
            .stroke_weight(1.0)
            .points(vertices)
            .xy(*position)
            .z(1.)
            // investigate this:
            .rotate(2. * PI + theta);        
        // let drawing = draw.ellipse()
        // .xy(*position)
        // .z(1.)
        // .radius(run_uptions.size);
        
        // If this instance is selected, show diagnostics
        if self.id == run_uptions.clicked_boid_id {
            // println!("vel {:#?}", velocity);
            // println!("rad {:#?}", theta);
            // println!("deg {:#?}", theta * 180. / PI);
            
            // Highlight the boid
            drawing.color(HOTPINK);

            // drawing a semi-sircle to demonstrate FOV
            // https://stackoverflow.com/a/72701981/9316685
            let radius = run_uptions.separation_treshold_distance * 2.;
            let section = Ellipse::new(
                Rect::from_x_y_w_h(0., 0., radius, radius), 120.
            )
            .section(
                theta - run_uptions.field_of_vision_half_rad,
                run_uptions.field_of_vision_half_rad * 2.
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

             // show alignment radius
            if run_uptions.alignment_on {
                show_radius(&draw, *position, run_uptions.sensory_distance * run_uptions.allignment_treshold_coefficient, RED);
            }
            
            // show cohesion radius
            if run_uptions.cohesion_on {
                show_radius(&draw, *position, run_uptions.sensory_distance * run_uptions.cohesion_treshold_coefficient, GREEN);
            }
            
            // show  separation radius
            if run_uptions.separation_on {
                show_radius(&draw, *position, run_uptions.sensory_distance * run_uptions.separation_treshold_coefficient, BLUE);
            }

            show_radius(&draw, *position, run_uptions.sensory_distance, YELLOW);

            // velocity vector scaled 100x
            draw.line()
            .start(self.position)
            .end(self.velocity * 100. + self.position)
            .color(LIME)
            .weight(4.0);
        } else {

        let saturation = if run_uptions.col_by_neighbour {
            color.saturation / 100. * 3. * self.n_neighbours as f32
        } else {
            color.saturation
        };


        // color.hue.


        // Take the color from UI settings
            drawing.hsv(
                if self.cluster_id != 0 {
                    ((10. * self.cluster_id as f32) % 360.) / 360.
                } else {
                    color.hue.to_positive_degrees() / 360.
                },
                saturation,
                if self.cluster_id != 0 {1f32} else {color.value}
            );
            // drawing.hsv(
            //     color.hue.to_positive_degrees() / 360.0,
            //     saturation,
            //     color.value * ( if self.cluster_id != 0 {1f32} else {0f32}),
            // );
        }

    }
}

fn show_radius<C>(draw: &Draw, xy: Vec2, r:f32, color: C) 
    where C: IntoLinSrgba<ColorScalar>,
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
        return
    }
    let draw = app.draw();
    
    let mouse_label = format!("{x:.2} - {y:.2}", x = app.mouse.x, y = app.mouse.y);

    draw.text(&mouse_label)
    .x_y(app.mouse.x, app.mouse.y)
    .color(RED);


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

    model.flock.draw(&draw, &model.color, &model.run_options);

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