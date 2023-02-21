use serde::{Deserialize, Serialize};

#[derive(Debug)]
pub struct RunOptions {
    pub init_boids: usize,
    pub initiation_strat: InitiationStrategy,

    pub baseline_speed: f32,
    pub min_speed: f32,
    pub max_speed: f32,
    pub max_steering: f32,

    pub allignment_coefficient: f32,
    pub cohesion_coefficient: f32,
    pub separation_coefficient: f32,
    
    pub sensory_distance: f32,
    pub max_sensory_distance: f32,
    
    pub allignment_treshold_distance: f32,
    pub cohesion_treshold_distance: f32,
    pub separation_treshold_distance: f32,
    
    pub allignment_treshold_coefficient: f32,
    pub cohesion_treshold_coefficient: f32,
    pub separation_treshold_coefficient: f32,
    
    
    pub alignment_on: bool,
    pub cohesion_on: bool,
    pub separation_on: bool,

    pub wander_coefficient: f32,
    pub wander_rate: f32,
    pub wander_radius: f32,
    pub wander_distance: f32,
    pub wander_on: bool,
    
    pub window: WindowSize,
    pub save_options: SaveOptions,

    pub size: f32,
    pub boundary: Boundary,
    pub distance: Distance,

    pub clicked_boid_id: usize,

    // for testing
    pub allignment_impl_mode: bool,
    pub cohesion_impl_mode: bool,
    pub separation_impl_mode: bool,
    pub col_by_neighbour: bool,

    pub field_of_vision_on: bool,
    /// Fiels of vision in radians, [0, π]
    pub field_of_vision_half_rad: f32,
    pub field_of_vision_deg: f32,

    pub sample_rate: u64,
    pub dbscan_flock_clustering_on: bool,
    pub neighbours_cosidered: usize,
    pub stop_movement: bool
}

impl RunOptions {
    pub fn update_sensory_distances(&mut self) {
        self.allignment_treshold_distance =
            self.sensory_distance * self.allignment_treshold_coefficient;
        self.cohesion_treshold_distance =
            self.sensory_distance * self.cohesion_treshold_coefficient;
        self.separation_treshold_distance =
            self.sensory_distance * self.separation_treshold_coefficient;

        self.max_sensory_distance = self.allignment_treshold_distance.max(
            self.cohesion_treshold_distance
                .max(self.separation_treshold_distance),
        );
    }
}

impl Default for RunOptions {
    fn default() -> Self {
        let init_boids = 256;
        let init_height = 600.0;
        let init_width = 600.0;

        let baseline_speed = 1.0;

        let min_speed = 0.65;
        let max_speed = 4.1;
        let max_steering = 0.7;

        let allignment_coefficient = 0.02;
        let cohesion_coefficient = 0.002;
        let separation_coefficient = 4.1;

        let sensory_distance = 60.;
        let max_sensory_distance = 60.;

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
        let field_of_vision_half_rad = 3. / 4. * std::f32::consts::PI / 2.;
        let field_of_vision_deg = 1. * 180. + 1.;

        let sample_rate = 1u64;

        let mut res = RunOptions {
            init_boids,
            // initiation_strat: InitiationStrategy::CircleCircumferenceIn,
            initiation_strat: InitiationStrategy::RandomRandom,
            baseline_speed,
            min_speed,
            max_speed,
            max_steering,
            allignment_coefficient,
            cohesion_coefficient,
            separation_coefficient,
            sensory_distance,
            max_sensory_distance,
            allignment_treshold_distance,
            cohesion_treshold_distance,
            separation_treshold_distance,
            allignment_treshold_coefficient,
            cohesion_treshold_coefficient,
            separation_treshold_coefficient,
            alignment_on,
            cohesion_on,
            separation_on,
            window: self::get_window_size(init_width, init_height),
            save_options: SaveOptions {
                save_locations: false,
                save_locations_path: None,
                save_locations_timestamp: true,
            },
            size: 8.,
            boundary: Boundary::Thoroidal,
            // boundary: Boundary::Repulsive{distance: 100., force: 0.05},
            // boundary: Boundary::Reflective,
            distance: Distance::EucEnclosed,
            clicked_boid_id: std::usize::MAX,
            allignment_impl_mode: false,
            cohesion_impl_mode: false,
            separation_impl_mode: false,
            col_by_neighbour: false,
            field_of_vision_on,
            field_of_vision_half_rad,
            field_of_vision_deg,
            sample_rate,
            dbscan_flock_clustering_on: false,
            neighbours_cosidered: 0,
            stop_movement: false,
            wander_rate: 0.03, 
            wander_on: false,
            wander_radius: 5.2,
            wander_coefficient: 0.4,
            wander_distance: 21.5,
        };

        res.update_sensory_distances();

        res
    }
}

pub fn get_window_size(init_width: f32, init_height: f32) -> WindowSize {
    WindowSize {
        win_left: init_width / -2.,
        win_right: init_width / 2.,
        win_top: init_height / 2.,
        win_bottom: init_height / -2.,
        win_h: init_height,
        win_w: init_width,
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

impl WindowSize {
    pub fn new(win_left: f32, win_right: f32, wind_top: f32, win_bottom: f32, win_h: f32, win_w: f32) -> WindowSize {
        WindowSize { win_left: win_left, win_right: win_right, win_top: wind_top, win_bottom: win_bottom, win_h: win_h, win_w: win_w }
    }
}

#[derive(Debug)]
pub enum InitiationStrategy {
    CircleCenterOut,
    CircleCircumferenceIn,
    TwoWalls,
    RectangleIn,
    RandomIn,
    RandomRandom,
}

#[derive(Serialize, Deserialize)]
#[serde(tag = "type")]
#[derive(Debug)]
pub enum Boundary {
    Thoroidal,
    Absorbing,
    Reflective,
    Repulsive{distance: f32, force: f32},
}

// {"type": "Repulsive", "distance": 100, "force": 0.05}

#[derive(Debug, PartialEq)]
pub enum Distance {
    EucThoroidal,
    EucEnclosed,
}

#[derive(Debug)]
pub struct SaveOptions {
    pub save_locations: bool,
    pub save_locations_path: Option<String>,
    pub save_locations_timestamp: bool,
}
