use std::f32::consts::PI;
use glam::Vec2;
use serde::{Deserialize, Serialize};

use crate::math_helpers::deg_to_half_rad;

#[derive(Debug, Clone)]
pub struct RunOptions {
    pub init_boids: usize,
    pub initiation_strat: InitiationStrategy,

    pub baseline_speed: f32,
    pub min_speed: f32,
    pub min_speed_sq: f32,
    pub max_speed: f32,
    pub max_speed_sq: f32,
    pub max_steering: f32,
    pub max_steering_sq: f32,

    pub alignment_coefficient: f32,
    pub cohesion_coefficient: f32,
    pub separation_coefficient: f32,

    pub sensory_distance: f32,
    pub max_sensory_distance: f32,

    pub alignment_treshold_distance: f32,
    pub cohesion_treshold_distance: f32,
    pub separation_treshold_distance: f32,

    pub alignment_treshold_coefficient: f32,
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
    pub noise_model: NoiseModel,
    pub tracker_type: TrackerType,

    pub clicked_boid_id: usize,

    // for testing
    pub alignment_impl_mode: bool,
    pub cohesion_impl_mode: bool,
    pub separation_impl_mode: bool,
    pub col_by_neighbour: bool,

    pub rules_impl: bool,

    pub field_of_vision_on: bool,
    pub field_of_vision_deg: f32,
    /// field of vision in radians, [0, π]
    pub field_of_vision_half_rad: f32,
    /// field of vision in cos [-1, 1]
    pub field_of_vision_cos: f32,

    pub alignment_fov_deg: f32,
    pub cohesion_fov_deg: f32,
    pub separation_fov_deg: f32,

    pub alignment_fov_half_cos: f32,
    pub cohesion_fov_half_cos: f32,
    pub separation_fov_half_cos: f32,

    pub sample_rate: u16,
    pub dbscan_flock_clustering_on: bool,
    pub neighbours_cosidered: usize,
    pub stop_movement: bool,
    pub seek_target_on: bool,
    pub seek_location: Option<Vec2>,
    pub agent_steering: bool,
}

impl RunOptions {
    /// updates all treshold distances given sensory_distance and all treshold coefficients
    pub fn update_sensory_distances(&mut self) {
        self.alignment_treshold_distance =
            self.sensory_distance * self.alignment_treshold_coefficient;
        self.cohesion_treshold_distance =
            self.sensory_distance * self.cohesion_treshold_coefficient;
        self.separation_treshold_distance =
            self.sensory_distance * self.separation_treshold_coefficient;

        self.max_sensory_distance = self.alignment_treshold_distance.max(
            self.cohesion_treshold_distance
                .max(self.separation_treshold_distance),
        );
    }

    /// updates fov parameters given fov_deg
    pub fn update_fov(&mut self) {
        self.field_of_vision_half_rad = self.field_of_vision_deg * PI / 360.;
        self.field_of_vision_half_rad = deg_to_half_rad(self.field_of_vision_deg);
        self.field_of_vision_cos = deg_to_half_rad(self.field_of_vision_deg).cos();
        
        self.alignment_fov_half_cos = deg_to_half_rad(self.alignment_fov_deg).cos();
        self.cohesion_fov_half_cos = deg_to_half_rad(self.cohesion_fov_deg).cos();
        self.separation_fov_half_cos = deg_to_half_rad(self.separation_fov_deg).cos();
     }
}

impl Default for RunOptions {
    fn default() -> Self {
        let init_boids = 256;
        let init_height = 600;
        let init_width = 600;

        let baseline_speed = 1.0;

        let min_speed = 0.65;
        let max_speed = 4.1;
        let max_steering = 0.7;

        let alignment_coefficient = 0.02;
        let cohesion_coefficient = 0.002;
        let separation_coefficient = 4.1;

        let sensory_distance = 60.;
        let max_sensory_distance = 60.;

        let alignment_treshold_distance = 200.;
        let cohesion_treshold_distance = 200.;
        let separation_treshold_distance = 200.;

        let alignment_treshold_coefficient = 1.15;
        let cohesion_treshold_coefficient = 0.95;
        let separation_treshold_coefficient = 0.35;

        let alignment_on = true;
        let cohesion_on = true;
        let separation_on = true;

        let field_of_vision_on = true;
        let field_of_vision_half_rad = 3. / 4. * std::f32::consts::PI / 2.;
        let field_of_vision_deg = 1. * 180. + 1.;

        let sample_rate = 1_u16;

        let mut res = RunOptions {
            init_boids,
            // initiation_strat: InitiationStrategy::CircleCircumferenceIn,
            initiation_strat: InitiationStrategy::RandomRandom,
            baseline_speed,
            min_speed,
            min_speed_sq: min_speed.powf(2.),
            max_speed,
            max_speed_sq: max_speed.powf(2.),
            max_steering,
            max_steering_sq: max_steering.powf(2.),
            alignment_coefficient,
            cohesion_coefficient,
            separation_coefficient,
            sensory_distance,
            max_sensory_distance,
            alignment_treshold_distance,
            cohesion_treshold_distance,
            separation_treshold_distance,
            alignment_treshold_coefficient,
            cohesion_treshold_coefficient,
            separation_treshold_coefficient,
            alignment_on,
            cohesion_on,
            separation_on,
            window: self::get_window_size(init_width, init_height),
            save_options: SaveOptions {
                save_locations: false,
                save_locations_path: Some("./".to_owned()),
                save_locations_timestamp: true,
            },
            size: 8.,
            // boundary: Boundary::Absorbing,
            boundary: Boundary::Toroidal,
            // boundary: Boundary::Repulsive{distance: 100., force: 0.05},
            // boundary: Boundary::Reflective,
            // distance: Distance::EucToroidal,
            distance: Distance::EucEnclosed,
            noise_model: NoiseModel::Reynolds,
            // tracker_type: TrackerType::Naive,
            tracker_type: TrackerType::SpatHash,
            // tracker_type: TrackerType::Replay("boids-data.csv".to_owned(), 0), // demo
            // tracker_type: TrackerType::Replay("/Users/filipvlcek/Source/Repos/boids_rust/boidranalysis/Data/0326_experiment3_w/prepro_boids-data_1679828486948.csv".to_owned(), 0),
            // tracker_type: TrackerType::Replay("boids-data_1679316314043.csv".to_owned(), 0), // demo
            // tracker_type: TrackerType::Replay("boids-data_1.csv".to_owned(), 0), // demo
            // tracker_type: TrackerType::Replay("boids-data_2.csv".to_owned(), 0), // demo
            // tracker_type: TrackerType::Replay("boids-data_3.csv".to_owned(), 0), // demok
            // tracker_type: TrackerType::Replay("/Users/filipvlcek/Source/Repos/boids_rust/boidranalysis/Data/0324_experiment2_k_11/prepro_boids-data_1679672635630.csv".to_owned(), 0),
            // tracker_type: TrackerType::Replay("/Users/filipvlcek/Source/Repos/boids_rust/boidranalysis/Data/0324_experiment2_k_01/prepro_boids-data_1679666885599.csv".to_owned(), 0),
            clicked_boid_id: std::usize::MAX,
            alignment_impl_mode: false,
            cohesion_impl_mode: false,
            separation_impl_mode: false,
            col_by_neighbour: false,
            field_of_vision_on,
            field_of_vision_deg,
            field_of_vision_half_rad,
            field_of_vision_cos: 0.,
            sample_rate,
            dbscan_flock_clustering_on: false,
            neighbours_cosidered: 0,
            stop_movement: false,
            wander_rate: 0.03,
            wander_on: false,
            wander_radius: 5.2,
            wander_coefficient: 0.4,
            wander_distance: 21.5,
            seek_target_on: false,
            seek_location: None,
            alignment_fov_deg: 45.,
            cohesion_fov_deg: 100.,
            separation_fov_deg: 135.,
            alignment_fov_half_cos: 0.,
            cohesion_fov_half_cos: 0.,
            separation_fov_half_cos: 0.,
            rules_impl: false,
            agent_steering: true,
        };

        res.update_sensory_distances();
        res.update_fov();

        res
    }
}

pub fn get_window_size(init_width: u32, init_height: u32) -> WindowSize {
    WindowSize {
        win_left: init_width as i32 / -2,
        win_right: init_width as i32 / 2,
        win_top: init_height as i32 / 2,
        win_bottom: init_height as i32 / -2,
        win_h: init_height as i32,
        win_w: init_width as i32,
    }
}

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub struct WindowSize {
    /// lowest x value
    pub win_left: i32,
    /// highest x value
    pub win_right: i32,
    /// highest y value
    pub win_top: i32,
    /// lowest y value
    pub win_bottom: i32,
    /// height
    pub win_h: i32,
    /// width
    pub win_w: i32,
}

impl WindowSize {
    pub fn new(
        win_left: i32,
        win_right: i32,
        wind_top: i32,
        win_bottom: i32,
        win_h: i32,
        win_w: i32,
    ) -> WindowSize {
        WindowSize {
            win_left: win_left,
            win_right: win_right,
            win_top: wind_top,
            win_bottom: win_bottom,
            win_h: win_h,
            win_w: win_w,
        }
    }
}

#[derive(Debug, Clone)]
pub enum InitiationStrategy {
    CircleCenterOut,
    CircleCircumferenceIn,
    TwoWalls,
    RectangleIn,
    RandomIn,
    RandomRandom,
}

#[derive(Debug, PartialEq, Copy, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
// {"type": "Repulsive", "distance": 100, "force": 0.05}
pub enum Boundary {
    Toroidal,
    Absorbing,
    Reflective,
    Repulsive { distance: f32, force: f32 },
    RepulsiveCircle { radius: f32}
}

#[derive(Debug, PartialEq, Copy, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
// {"type": "EucToroidal"}
pub enum Distance {
    EucToroidal,
    EucEnclosed,
}

#[derive(Debug, PartialEq, Copy, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
// {"type": "EucToroidal"}
pub enum NoiseModel {
    Vicsek,
    Reynolds,
}

#[derive(Debug, Clone)]
pub enum TrackerType {
    SpatHash,
    Naive,
    Replay(String, u64)
}

#[derive(Debug, Clone)]
pub struct SaveOptions {
    pub save_locations: bool,
    pub save_locations_path: Option<String>,
    pub save_locations_timestamp: bool,
}
