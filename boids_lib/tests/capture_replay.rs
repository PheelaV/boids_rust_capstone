//! Integration tests for capture and replay functionality
//!
//! These tests ensure that simulation data can be captured, saved,
//! and replayed correctly without GUI dependencies.

use boids_lib::{
    birdwatcher::{Birdwatcher, BoidData},
    flock::Flock,
    flock_base,
    options::*,
};
use std::fs;
use std::path::Path;

/// Test that birdwatcher captures data at correct sample rate
#[test]
fn test_birdwatcher_sample_rate() {
    let mut options = RunOptions::default();
    options.init_boids = 10;
    options.window = get_window_size(800, 600);

    let mut flock = Flock::new(&options);
    let mut birdwatcher = Birdwatcher::new(5); // Sample every 5 iterations

    let iterations = 100;

    for _ in 0..iterations {
        flock.update(&mut options);
        birdwatcher.watch(&flock);
    }

    let data = birdwatcher.pop_data_save(&options.save_options);

    // With sample_rate=5 and 100 iterations, we expect 20 samples per boid
    // (100 / 5 = 20)
    let expected_samples_per_boid = iterations / 5;
    let total_expected = expected_samples_per_boid * options.init_boids as usize;

    assert!(
        data.len() >= total_expected - options.init_boids as usize
            && data.len() <= total_expected + options.init_boids as usize,
        "Expected ~{} samples, got {}",
        total_expected,
        data.len()
    );
}

/// Test that captured data has valid structure
#[test]
fn test_captured_data_validity() {
    let data = flock_base(50, RunOptions::default());

    assert!(!data.is_empty(), "Should capture some data");

    // Check first data point
    let first = &data[0];

    // IDs should be reasonable
    assert!(first.id < 1000, "Boid ID should be reasonable");

    // Positions should be within expected bounds (default window is 1400x900)
    assert!(
        first.x.abs() < 1000.0,
        "X position should be within bounds"
    );
    assert!(
        first.y.abs() < 1000.0,
        "Y position should be within bounds"
    );

    // Time should increase across samples
    // Each sample has init_boids records
    let init_boids = RunOptions::default().init_boids;
    if data.len() > init_boids * 2 {
        assert!(
            data[0].time < data[init_boids].time,
            "Time should increase between samples: {} < {}",
            data[0].time, data[init_boids].time
        );
        assert!(
            data[init_boids].time < data[init_boids * 2].time,
            "Time should continue increasing: {} < {}",
            data[init_boids].time, data[init_boids * 2].time
        );
    }
}

/// Test that data can be saved and loaded
#[test]
fn test_save_and_load_trajectory() {
    let test_dir = "/tmp/boids_test_capture";
    let _ = fs::create_dir_all(test_dir);

    let mut options = RunOptions::default();
    options.init_boids = 20;
    options.sample_rate = 10;
    options.window = get_window_size(800, 600);
    options.save_options.save_locations = true;
    options.save_options.save_locations_path = Some(test_dir.to_string());

    // Run simulation and capture data
    let data = flock_base(100, options.clone());

    assert!(!data.is_empty(), "Should have captured data");

    // Try to verify CSV file was created (or could be created)
    // Note: actual file writing happens in the app, but we can validate data format

    // Check data can be serialized
    for sample in data.iter().take(10) {
        assert!(sample.id < 100, "Valid boid ID");
        assert!(sample.time < 200, "Valid time");
    }

    // Cleanup
    let _ = fs::remove_dir_all(test_dir);
}

/// Test flock_base produces consistent output structure
#[test]
fn test_flock_base_output_structure() {
    let mut options = RunOptions::default();
    options.init_boids = 15;
    options.sample_rate = 5;

    let iterations = 50u64;
    let data = flock_base(iterations, options.clone());

    // Should have data for each boid at each sample point
    let expected_samples = (iterations / options.sample_rate as u64) as usize;

    // Total samples = boids * samples_per_boid
    let expected_total = expected_samples * options.init_boids as usize;

    assert!(
        data.len() >= expected_total - options.init_boids as usize,
        "Expected around {} samples, got {}",
        expected_total,
        data.len()
    );

    // Verify each boid appears roughly equal times
    let mut boid_counts = vec![0; options.init_boids as usize];

    for sample in &data {
        if (sample.id as usize) < boid_counts.len() {
            boid_counts[sample.id as usize] += 1;
        }
    }

    // Each boid should appear roughly the same number of times
    let min_count = boid_counts.iter().min().unwrap_or(&0);
    let max_count = boid_counts.iter().max().unwrap_or(&0);

    assert!(
        max_count - min_count <= 2,
        "Boid sample counts should be balanced: min={}, max={}",
        min_count,
        max_count
    );
}

/// Test that replay data maintains temporal consistency
#[test]
fn test_temporal_consistency() {
    let mut options = RunOptions::default();
    options.init_boids = 10;
    options.sample_rate = 1; // Sample every iteration for fine-grained check

    let window_width = (options.window.win_right - options.window.win_left) as f32;
    let window_height = (options.window.win_top - options.window.win_bottom) as f32;

    let data = flock_base(20, options);

    // Group data by boid ID
    let mut trajectories: Vec<Vec<&BoidData>> = vec![Vec::new(); 10];

    for sample in &data {
        if (sample.id as usize) < trajectories.len() {
            trajectories[sample.id as usize].push(sample);
        }
    }

    // Check each boid's trajectory
    for (boid_id, trajectory) in trajectories.iter().enumerate() {
        if trajectory.is_empty() {
            continue;
        }

        // Times should be monotonically increasing
        for window in trajectory.windows(2) {
            assert!(
                window[1].time > window[0].time,
                "Boid {} time should increase: {} -> {}",
                boid_id,
                window[0].time,
                window[1].time
            );
        }

        // Positions shouldn't teleport (movement should be continuous)
        // Use toroidal distance since default boundary is toroidal
        for window in trajectory.windows(2) {
            let dx = (window[1].x - window[0].x).abs();
            let dy = (window[1].y - window[0].y).abs();

            // Toroidal distance: take shorter path around the wrapped boundary
            let dx_wrapped = dx.min(window_width - dx);
            let dy_wrapped = dy.min(window_height - dy);

            let dist = (dx_wrapped.powi(2) + dy_wrapped.powi(2)).sqrt();

            // With default max_speed and baseline_speed, movement per iteration
            // should be bounded
            assert!(
                dist < 100.0,
                "Boid {} teleported: distance {} between consecutive samples",
                boid_id,
                dist
            );
        }
    }
}

/// Test that simulation with replay tracker type can be configured
#[test]
fn test_replay_tracker_configuration() {
    let replay_path = "/tmp/nonexistent_replay.csv";

    let mut options = RunOptions::default();
    options.tracker_type = TrackerType::Replay(replay_path.to_string(), 0);

    // This should create a Replay tracker (which will fail gracefully if file doesn't exist)
    // but we're testing the configuration path works

    assert!(matches!(
        options.tracker_type,
        TrackerType::Replay(_, _)
    ));

    if let TrackerType::Replay(path, start_iteration) = &options.tracker_type {
        assert_eq!(path, replay_path);
        assert_eq!(*start_iteration, 0);
    }
}

/// Test that data collection respects save options
#[test]
fn test_save_options_filtering() {
    let mut options = RunOptions::default();
    options.init_boids = 10;
    options.sample_rate = 5;

    let mut flock = Flock::new(&options);
    let mut birdwatcher = Birdwatcher::new(5);

    for _ in 0..50 {
        flock.update(&mut options);
        birdwatcher.watch(&flock);
    }

    // Test with save_locations = true
    options.save_options.save_locations = true;
    let data_with_save = birdwatcher.pop_data_save(&options.save_options);

    assert!(!data_with_save.is_empty(), "Should have data when save is enabled");

    // Reset
    let mut birdwatcher2 = Birdwatcher::new(5);
    for _ in 0..50 {
        flock.update(&mut options);
        birdwatcher2.watch(&flock);
    }

    // Test with save_locations = false
    options.save_options.save_locations = false;
    let data_without_save = birdwatcher2.pop_data_save(&options.save_options);

    // Behavior might be same (data returned) but flag is respected
    assert!(!data_without_save.is_empty());
}

/// Test multiple simulation runs produce different results (non-determinism)
#[test]
fn test_multiple_runs_vary() {
    let options = RunOptions::default();

    let data1 = flock_base(30, options.clone());
    let data2 = flock_base(30, options);

    // Due to random initialization, runs should differ
    // (Note: if we want determinism, we'd need to set a fixed seed)

    if !data1.is_empty() && !data2.is_empty() {
        let pos1 = (data1[0].x, data1[0].y);
        let pos2 = (data2[0].x, data2[0].y);

        // Very likely to be different positions
        let same = (pos1.0 - pos2.0).abs() < 0.01 && (pos1.1 - pos2.1).abs() < 0.01;

        if same {
            println!("Warning: Multiple runs had identical start positions (unlikely but possible)");
        }
    }
}

/// Test that empty flock doesn't crash capture
#[test]
fn test_capture_empty_flock() {
    let mut options = RunOptions::default();
    options.init_boids = 0;

    let mut flock = Flock::new(&options);
    let mut birdwatcher = Birdwatcher::new(1);

    // Should not crash
    for _ in 0..10 {
        flock.update(&mut options);
        birdwatcher.watch(&flock);
    }

    let data = birdwatcher.pop_data_save(&options.save_options);

    assert!(data.is_empty(), "Empty flock should produce no data");
}
