mod utils;

use boids_lib::{
    flock::Flock,
    options::{self, RunOptions, WindowSize, TrackerType, InitiationStrategy, Distance, Boundary, NoiseModel, SaveOptions, NeighbourSampling},
};
use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

// When the `wee_alloc` feature is enabled, use `wee_alloc` as the global allocator.
#[cfg(feature = "wee_alloc")]
#[global_allocator]
static ALLOC: wee_alloc::WeeAlloc = wee_alloc::WeeAlloc::INIT;

/// Configuration for the simulation, serializable from JavaScript
#[derive(Serialize, Deserialize)]
#[wasm_bindgen(getter_with_clone)]
pub struct SimulationConfig {
    pub init_boids: usize,
    pub window_width: usize,
    pub window_height: usize,
    pub separation_coefficient: f32,
    pub cohesion_coefficient: f32,
    pub alignment_coefficient: f32,
    pub max_speed: f32,
    pub min_speed: f32,
    pub sensory_distance: f32,
    pub wander_on: bool,
    pub wander_coefficient: f32,
}

impl Default for SimulationConfig {
    fn default() -> Self {
        Self {
            init_boids: 200,
            window_width: 800,
            window_height: 600,
            separation_coefficient: 1.5,
            cohesion_coefficient: 1.0,
            alignment_coefficient: 1.0,
            max_speed: 4.0,
            min_speed: 2.0,
            sensory_distance: 60.0,
            wander_on: true,
            wander_coefficient: 0.5,
        }
    }
}

/// Represents a single boid's state for rendering
#[derive(Serialize, Deserialize)]
pub struct BoidState {
    pub id: usize,
    pub x: f32,
    pub y: f32,
    pub vx: f32,
    pub vy: f32,
}

/// Main WebAssembly simulation wrapper
#[wasm_bindgen]
pub struct WasmSimulation {
    flock: Flock<'static>,
    options: RunOptions,
    frame_count: u64,
    /// Original boid count for reset - doesn't change when adding/removing boids
    initial_boid_count: usize,
}

#[wasm_bindgen]
impl WasmSimulation {
    /// Create a new simulation with default configuration
    #[wasm_bindgen(constructor)]
    pub fn new() -> Result<WasmSimulation, JsValue> {
        utils::set_panic_hook();

        let config = SimulationConfig::default();
        Self::from_config(serde_wasm_bindgen::to_value(&config)?)
    }

    /// Create a simulation from a JavaScript configuration object
    pub fn from_config(config: JsValue) -> Result<WasmSimulation, JsValue> {
        utils::set_panic_hook();

        let config: SimulationConfig = serde_wasm_bindgen::from_value(config)?;
        let initial_boid_count = config.init_boids;
        let options = Self::config_to_run_options(config);
        let flock = Flock::new(&options);

        Ok(WasmSimulation {
            flock,
            options,
            frame_count: 0,
            initial_boid_count,
        })
    }

    /// Update the simulation by one step
    pub fn update(&mut self) {
        self.flock.update(&mut self.options);
        self.frame_count += 1;
    }

    /// Get the current state of all boids as a JSON object
    pub fn get_boids(&self) -> Result<JsValue, JsValue> {
        let boids: Vec<BoidState> = self
            .flock
            .view2()
            .map(|(boid, _)| BoidState {
                id: boid.id,
                x: boid.position.x,
                y: boid.position.y,
                vx: boid.velocity.x,
                vy: boid.velocity.y,
            })
            .collect();

        serde_wasm_bindgen::to_value(&boids).map_err(|e| JsValue::from_str(&e.to_string()))
    }

    /// Get the number of boids in the simulation
    pub fn get_boid_count(&self) -> usize {
        self.flock.tracker.get_no_entities()
    }

    /// Get current frame count
    pub fn get_frame_count(&self) -> u64 {
        self.frame_count
    }

    /// Set separation coefficient
    pub fn set_separation_coefficient(&mut self, value: f32) {
        self.options.separation_coefficient = value;
    }

    /// Set cohesion coefficient
    pub fn set_cohesion_coefficient(&mut self, value: f32) {
        self.options.cohesion_coefficient = value;
    }

    /// Set alignment coefficient
    pub fn set_alignment_coefficient(&mut self, value: f32) {
        self.options.alignment_coefficient = value;
    }

    /// Set max speed
    pub fn set_max_speed(&mut self, value: f32) {
        self.options.max_speed = value;
        self.options.max_speed_sq = value * value;
    }

    /// Toggle separation behavior
    pub fn toggle_separation(&mut self) {
        self.options.separation_on = !self.options.separation_on;
    }

    /// Toggle cohesion behavior
    pub fn toggle_cohesion(&mut self) {
        self.options.cohesion_on = !self.options.cohesion_on;
    }

    /// Toggle alignment behavior
    pub fn toggle_alignment(&mut self) {
        self.options.alignment_on = !self.options.alignment_on;
    }

    /// Toggle wander behavior
    pub fn toggle_wander(&mut self) {
        self.options.wander_on = !self.options.wander_on;
    }

    /// Set the window/canvas size
    pub fn set_window_size(&mut self, width: i32, height: i32) {
        self.options.window = options::get_window_size(width as u32, height as u32);
    }

    /// Get the current window width
    pub fn get_window_width(&self) -> i32 {
        self.options.window.win_w
    }

    /// Get the current window height
    pub fn get_window_height(&self) -> i32 {
        self.options.window.win_h
    }

    /// Add a new boid at a specific position with velocity
    pub fn add_boid(&mut self, _x: f32, _y: f32, _vx: f32, _vy: f32) {
        // Note: This is a simplified version - boids_lib doesn't expose a way to add boids
        // with custom positions/velocities through the Flock API yet.
        // For now, we just add a new boid with default initialization.
        self.flock.insert_single(&self.options);
    }

    /// Add multiple boids (for doubling, etc.)
    pub fn add_boids(&mut self, count: usize) {
        for _ in 0..count {
            self.flock.insert_single(&self.options);
        }
        self.options.init_boids = self.get_boid_count();
    }

    /// Remove the last boid
    pub fn remove_boid(&mut self) -> bool {
        let removed = self.flock.delete_last(&self.options).is_some();
        if removed {
            self.options.init_boids = self.get_boid_count();
        }
        removed
    }

    /// Remove multiple boids
    pub fn remove_boids(&mut self, count: usize) -> usize {
        let mut removed = 0;
        for _ in 0..count {
            if self.flock.delete_last(&self.options).is_some() {
                removed += 1;
            } else {
                break;
            }
        }
        if removed > 0 {
            self.options.init_boids = self.get_boid_count();
        }
        removed
    }

    /// Double the number of boids
    pub fn double_boids(&mut self) {
        let current_count = self.get_boid_count();
        self.add_boids(current_count);
    }

    /// Halve the number of boids
    pub fn halve_boids(&mut self) {
        let current_count = self.get_boid_count();
        let to_remove = current_count / 2;
        self.remove_boids(to_remove);
    }

    /// Set the number of boids to a specific count
    pub fn set_boid_count(&mut self, target_count: usize) {
        let current_count = self.get_boid_count();

        if target_count > current_count {
            let to_add = target_count - current_count;
            self.add_boids(to_add);
        } else if target_count < current_count {
            let to_remove = current_count - target_count;
            self.remove_boids(to_remove);
        }
    }

    /// Reset the simulation to initial state (original boid count)
    pub fn reset(&mut self) {
        // Restore original boid count before restart
        self.options.init_boids = self.initial_boid_count;
        self.flock.restart(&self.options);
        self.frame_count = 0;
    }

    /// Get simulation statistics as a JSON object
    pub fn get_stats(&self) -> Result<JsValue, JsValue> {
        #[derive(Serialize)]
        struct Stats {
            boid_count: usize,
            frame_count: u64,
            separation_on: bool,
            cohesion_on: bool,
            alignment_on: bool,
            wander_on: bool,
        }

        let stats = Stats {
            boid_count: self.flock.tracker.get_no_entities(),
            frame_count: self.frame_count,
            separation_on: self.options.separation_on,
            cohesion_on: self.options.cohesion_on,
            alignment_on: self.options.alignment_on,
            wander_on: self.options.wander_on,
        };

        serde_wasm_bindgen::to_value(&stats).map_err(|e| JsValue::from_str(&e.to_string()))
    }

    // Helper function to convert SimulationConfig to RunOptions
    fn config_to_run_options(config: SimulationConfig) -> RunOptions {
        let win_w = config.window_width as i32;
        let win_h = config.window_height as i32;

        RunOptions {
            init_boids: config.init_boids,
            initiation_strat: InitiationStrategy::RandomRandom,
            baseline_speed: 2.0,
            min_speed: config.min_speed,
            min_speed_sq: config.min_speed * config.min_speed,
            max_speed: config.max_speed,
            max_speed_sq: config.max_speed * config.max_speed,
            max_steering: 0.7,
            max_steering_sq: 0.49,
            alignment_coefficient: config.alignment_coefficient,
            cohesion_coefficient: config.cohesion_coefficient,
            separation_coefficient: config.separation_coefficient,
            sensory_distance: config.sensory_distance,
            max_sensory_distance: config.sensory_distance * 1.15,
            alignment_treshold_distance: config.sensory_distance * 1.15,
            cohesion_treshold_distance: config.sensory_distance * 0.95,
            separation_treshold_distance: config.sensory_distance * 0.35,
            alignment_treshold_coefficient: 1.15,
            cohesion_treshold_coefficient: 0.95,
            separation_treshold_coefficient: 0.35,
            alignment_on: true,
            cohesion_on: true,
            separation_on: true,
            wander_coefficient: config.wander_coefficient,
            wander_rate: 0.05,
            wander_radius: 1.0,
            wander_distance: 2.0,
            wander_on: config.wander_on,
            window: WindowSize {
                win_w,
                win_h,
                win_left: -win_w / 2,
                win_right: win_w / 2,
                win_bottom: -win_h / 2,
                win_top: win_h / 2,
            },
            save_options: SaveOptions {
                save_locations: false,
                save_locations_path: None,
                save_locations_timestamp: false,
            },
            size: 5.0,
            boundary: Boundary::Toroidal,
            distance: Distance::EucToroidal,
            noise_model: NoiseModel::Reynolds,
            tracker_type: TrackerType::SpatHash,
            clicked_boid_id: usize::MAX,
            alignment_impl_mode: false,
            cohesion_impl_mode: false,
            separation_impl_mode: false,
            col_by_neighbour: false,
            rules_impl: true,
            field_of_vision_on: true,
            field_of_vision_deg: 270.0,
            field_of_vision_half_rad: 3.0 * std::f32::consts::PI / 4.0,
            field_of_vision_cos: (3.0 * std::f32::consts::PI / 4.0).cos(),
            alignment_fov_deg: 270.0,
            cohesion_fov_deg: 270.0,
            separation_fov_deg: 270.0,
            alignment_fov_half_cos: (3.0 * std::f32::consts::PI / 4.0).cos(),
            cohesion_fov_half_cos: (3.0 * std::f32::consts::PI / 4.0).cos(),
            separation_fov_half_cos: (3.0 * std::f32::consts::PI / 4.0).cos(),
            sample_rate: 1,
            dbscan_flock_clustering_on: false,
            neighbours_cosidered: 10,
            stop_movement: false,
            seek_target_on: false,
            seek_location: None,
            agent_steering: true,
            neighbour_sampling: NeighbourSampling::default(),
            rng_seed: None, // Use random seed for web simulations
        }
    }
}
