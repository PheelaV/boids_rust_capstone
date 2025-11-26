//! Scalability benchmark for boids simulation
//!
//! Tests performance across different boid counts (2^14, 2^15, 2^16)
//! with fixed sensory_distance=20, based on 2normal_s.toml config.
//!
//! Run with: cargo bench --bench flock_scalability

use std::time::Duration;

use boids_lib::flock_base;
use boids_lib::options::{self, RunOptions, SaveOptions};
use criterion::{
    criterion_group, criterion_main, BenchmarkId, Criterion, SamplingMode, Throughput,
};

/// Creates RunOptions based on 2normal_s.toml config with specified boid count
fn run_options_2normal_s(no_boids: usize) -> RunOptions {
    let mut ro: RunOptions = Default::default();

    // Core simulation params
    ro.init_boids = no_boids;
    ro.window = options::get_window_size(1000, 1000);
    ro.sensory_distance = 20.0; // Fixed at 20 for this benchmark

    // Behavior coefficients from 2normal_s.toml
    ro.alignment_coefficient = 0.2;
    ro.cohesion_coefficient = 0.1;
    ro.separation_coefficient = 2.5;
    ro.alignment_treshold_coefficient = 1.0;
    ro.cohesion_treshold_coefficient = 0.8;
    ro.separation_treshold_coefficient = 0.3;

    // Speed/steering from config
    ro.baseline_speed = 1.0;
    ro.min_speed = 2.0;
    ro.min_speed_sq = ro.min_speed * ro.min_speed;
    ro.max_speed = 2.0;
    ro.max_speed_sq = ro.max_speed * ro.max_speed;
    ro.max_steering = 0.3;
    ro.max_steering_sq = ro.max_steering * ro.max_steering;

    // FOV and wander from config
    ro.field_of_vision_deg = 280.0;
    ro.wander_on = false;
    ro.size = 5.0;
    ro.agent_steering = true;
    ro.rules_impl = false;

    // Disable saving for benchmarks
    ro.save_options = SaveOptions {
        save_locations: false,
        save_locations_path: None,
        save_locations_timestamp: false,
    };

    // Use deterministic seed for reproducibility
    ro.rng_seed = Some(42);

    ro.update_sensory_distances();
    ro.update_fov();
    ro
}

fn scalability_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("flock_scalability");

    // Test 2^14, 2^15, 2^16 boids with (exp, boid_count, iterations)
    let boid_counts: [(u32, usize, u64); 3] = [
        (14, 1 << 14, 256), // 16,384 boids
        (15, 1 << 15, 256), // 32,768 boids
        (16, 1 << 16, 25),  // 65,536 boids - fewer iterations due to cost
    ];

    for (exp, no_boids, iterations) in boid_counts {
        group.sample_size(10); // Criterion minimum is 10
        group.measurement_time(Duration::from_secs(30));
        group.sampling_mode(SamplingMode::Flat);
        group.throughput(Throughput::Elements(no_boids as u64 * iterations));

        group.bench_with_input(
            BenchmarkId::new("boids", format!("2^{}", exp)),
            &no_boids,
            |b, &n| {
                b.iter(|| flock_base(iterations, run_options_2normal_s(n)));
            },
        );
    }
    group.finish();
}

criterion_group!(benches, scalability_benchmark);
criterion_main!(benches);
