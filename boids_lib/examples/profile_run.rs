//! Standalone binary for profiling boids simulation.
//!
//! Build and run with a profiler:
//! ```bash
//! cargo build --release --example profile_run
//! samply record ./target/release/examples/profile_run
//! ```

use boids_lib::flock_base;
use boids_lib::options::{self, RunOptions, SaveOptions};

fn main() {
    let mut ro = RunOptions::default();

    // Match benchmark configuration
    ro.init_boids = 1 << 15; // 32K boids - same as 2^15 benchmark
    ro.window = options::get_window_size(1000, 1000);
    ro.sensory_distance = 20.0;

    // Behavior coefficients (NOT scaled - max_steering handles SI conversion)
    ro.alignment_coefficient = 0.2;
    ro.cohesion_coefficient = 0.1;
    ro.separation_coefficient = 2.5;
    ro.alignment_treshold_coefficient = 1.0;
    ro.cohesion_treshold_coefficient = 0.8;
    ro.separation_treshold_coefficient = 0.3;

    // Speed/steering (SI units: pixels/second, pixels/second²)
    ro.max_speed = 120.0;  // pixels/second
    ro.max_speed_sq = ro.max_speed * ro.max_speed;
    ro.max_steering = 18.0;  // pixels/second²
    ro.max_steering_sq = ro.max_steering * ro.max_steering;

    // FOV and wander
    ro.field_of_vision_deg = 280.0;
    ro.wander_on = false;

    ro.agent_steering = true;
    ro.rules_impl = false;

    // Disable saving
    ro.save_options = SaveOptions {
        save_locations: false,
        save_locations_path: None,
        save_locations_timestamp: false,
    };

    // Deterministic seed for reproducibility
    ro.rng_seed = Some(42);

    ro.update_sensory_distances();
    ro.update_fov();

    println!("Running {} boids for 100 iterations...", ro.init_boids);
    flock_base(100, ro);
    println!("Done.");
}
