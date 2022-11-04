use nannou::prelude::*;

pub struct Options {
    pub init_boids: usize,
    pub init_width: u32,
    pub init_height: u32,
    pub run_options: RunOptions,
    pub initiation_strat: InitiationStrategy,
    pub record_locations: bool,
    pub no_boids: u32,
}

// #[derive(Debug, Clone, Copy)]


#[derive(Debug)]
pub struct RunOptions {
    pub baseline_speed: f32,
    pub min_speed: f32,
    pub max_speed: f32,
    pub max_steering: f32,

    pub allignment_coefficient: f32,
    pub cohesion_coefficient: f32,
    pub separation_coefficient: f32,

    pub sensory_distance: f32,

    pub allignment_treshold_distance: f32,
    pub cohesion_treshold_distance: f32,
    pub separation_treshold_distance: f32,

    pub allignment_treshold_coefficient: f32,
    pub cohesion_treshold_coefficient: f32,
    pub separation_treshold_coefficient: f32,

    pub alignment_on: bool,
    pub cohesion_on: bool,
    pub separation_on: bool,

    pub window: WindowSize,

    pub size: f32,
    pub boundary: Boundary,
    pub distance: Distance,

    pub clicked_boid_id: u32,
    pub no_boids: u32,

    pub impl_mode: bool,

    pub field_of_vision_on: bool,
    /// Fiels of vision in radians, [0, π]
    pub field_of_vision_half_rad: f32,
    pub field_of_vision_deg: f32,
}

pub fn get_run_options(win: &Rect) -> RunOptions {
        let baseline_speed = 1.0;

        let min_speed = 0.65;
        let max_speed = 4.1;
        let max_steering = 0.7;

        let allignment_coefficient = 0.02;
        let cohesion_coefficient = 0.002;
        let separation_coefficient = 4.1;

        let sensory_distance = 60.;

        let allignment_treshold_distance = 200.;
        let cohesion_treshold_distance = 200.;
        let separation_treshold_distance = 200.;

        let allignment_treshold_coefficient = 1.15;
        let cohesion_treshold_coefficient = 0.95;
        let separation_treshold_coefficient = 0.35;

        let alignment_on = true;
        let cohesion_on = true;
        let separation_on = true;

        let field_of_vision_on = true;
        let field_of_vision_half_rad = 3./4. * PI / 2.;
        let field_of_vision_deg = 1. * 180. + 1.;

    RunOptions {
        baseline_speed,
        min_speed,
        max_speed,
        max_steering,
        allignment_coefficient,
        cohesion_coefficient,
        separation_coefficient,
        sensory_distance,
        allignment_treshold_distance,
        cohesion_treshold_distance,
        separation_treshold_distance,
        allignment_treshold_coefficient,
        cohesion_treshold_coefficient,
        separation_treshold_coefficient,
        alignment_on,
        cohesion_on,
        separation_on,        
        window: WindowSize {
            win_left: win.left(),
            win_right: win.right(),
            win_top: win.top(),
            win_bottom: win.bottom(),
            win_h: win.h(),
            win_w: win.w(),
        },
        size: 8.,
        boundary: Boundary::Thoroidal,
        distance: Distance::EucEnclosed,
        clicked_boid_id: std::u32::MAX,
        no_boids: 2,
        impl_mode: false,
        field_of_vision_on,
        field_of_vision_half_rad,
        field_of_vision_deg
    }
}

#[derive(Debug)]
pub struct WindowSize {
        /// lowest x value
        pub win_left: f32,
        /// highest x value
        pub win_right: f32,
        /// highest y value
        pub win_top: f32,
        /// lowest y value
        pub win_bottom: f32,
        /// height
        pub win_h: f32,
        /// width
        pub win_w: f32,
}

pub enum InitiationStrategy {
    CircleCenterOut,
    CircleCircumferenceIn,
    TwoWalls,
    RectangleIn,
    RandomIn,
    RandomRandom
}

#[derive(Debug)]
pub enum Boundary {
    Thoroidal,
    Absorbing,
    Reflective,
}

#[derive(Debug, PartialEq)]
pub enum Distance {
    EucThoroidal,
    EucEnclosed,
}