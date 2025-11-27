//! Integration tests for simulation behavior and reproducibility
//!
//! These tests validate core behavioral properties of the simulation.
//! The determinism test is in a separate file (determinism_isolated.rs)
//! to avoid RNG interference from parallel tests.

use boids_lib::{flock_base, options::*, flock::Flock};
use glam::Vec2;

/// Test that simulation progresses (boids move)
#[test]
fn test_boids_move_over_time() {
    let mut options = RunOptions::default();
    options.init_boids = 10;
    options.window = get_window_size(800, 600);
    options.stop_movement = false;

    let mut flock = Flock::new(&options);

    // Capture initial positions
    let initial_positions: Vec<Vec2> = flock
        .view2()
        .map(|(boid, _)| boid.position)
        .collect();

    // Run simulation for a few steps
    for _ in 0..10 {
        flock.update(&mut options);
    }

    // Capture final positions
    let final_positions: Vec<Vec2> = flock
        .view2()
        .map(|(boid, _)| boid.position)
        .collect();

    // At least some boids should have moved
    let mut moved_count = 0;
    for (initial, final_pos) in initial_positions.iter().zip(final_positions.iter()) {
        let distance = (*final_pos - *initial).length();
        if distance > 0.1 {
            moved_count += 1;
        }
    }

    assert!(
        moved_count > 0,
        "At least some boids should move during simulation"
    );
}

/// Test that boid count remains stable (no creation/deletion)
#[test]
fn test_boid_count_stability() {
    let mut options = RunOptions::default();
    options.init_boids = 100;
    options.window = get_window_size(1000, 1000);

    let mut flock = Flock::new(&options);

    let initial_count = flock.view2().count();

    // Run for many iterations
    for _ in 0..1000 {
        flock.update(&mut options);
    }

    let final_count = flock.view2().count();

    assert_eq!(
        initial_count, final_count,
        "Boid count should remain stable during simulation"
    );
    assert_eq!(initial_count, 100, "Should have exactly 100 boids");
}

/// Test that velocities stay within configured bounds
#[test]
fn test_velocity_bounds() {
    let mut options = RunOptions::default();
    options.init_boids = 50;
    options.max_speed = 5.0;
    options.min_speed = 1.0;
    options.window = get_window_size(800, 600);

    let mut flock = Flock::new(&options);

    // Run simulation
    for _ in 0..100 {
        flock.update(&mut options);

        // Check all velocities are within bounds
        for (boid, _) in flock.view2() {
            let speed = boid.velocity.length();
            assert!(
                speed >= options.min_speed * 0.99, // Small tolerance for numerical errors
                "Speed {} should be >= min_speed {}",
                speed,
                options.min_speed
            );
            assert!(
                speed <= options.max_speed * 1.01,
                "Speed {} should be <= max_speed {}",
                speed,
                options.max_speed
            );
        }
    }
}

/// Test that sample rate correctly controls data collection
#[test]
fn test_sample_rate_collection() {
    let mut options = RunOptions::default();
    options.init_boids = 20;
    options.sample_rate = 5; // Collect every 5 iterations
    options.window = get_window_size(800, 600);

    let iterations = 100u64;
    let result = flock_base(iterations, options);

    // With sample_rate=5 and 100 iterations, we expect ~20 samples per boid
    // (iterations / sample_rate)
    let expected_samples_per_boid = (iterations / 5) as usize;
    let total_boids = 20;

    // Total samples should be approximately boids * samples_per_boid
    let expected_total = expected_samples_per_boid * total_boids;

    assert!(
        result.len() >= expected_total - total_boids
            && result.len() <= expected_total + total_boids,
        "Expected around {} samples, got {}",
        expected_total,
        result.len()
    );
}

/// Test that boids respect window boundaries in different modes
#[test]
fn test_boundary_containment_toroidal() {
    let mut options = RunOptions::default();
    options.init_boids = 30;
    options.window = get_window_size(400, 400);
    options.boundary = Boundary::Toroidal;
    options.max_speed = 10.0; // Fast movement to test boundaries

    let mut flock = Flock::new(&options);

    for _ in 0..200 {
        flock.update(&mut options);

        // In toroidal mode, boids should wrap around and stay within bounds
        for (boid, _) in flock.view2() {
            assert!(
                boid.position.x >= options.window.win_left as f32
                    && boid.position.x <= options.window.win_right as f32,
                "Boid x position {} outside bounds [{}, {}]",
                boid.position.x,
                options.window.win_left,
                options.window.win_right
            );
            assert!(
                boid.position.y >= options.window.win_bottom as f32
                    && boid.position.y <= options.window.win_top as f32,
                "Boid y position {} outside bounds [{}, {}]",
                boid.position.y,
                options.window.win_bottom,
                options.window.win_top
            );
        }
    }
}

#[test]
fn test_boundary_containment_reflective() {
    let mut options = RunOptions::default();
    options.init_boids = 30;
    options.window = get_window_size(400, 400);
    options.boundary = Boundary::Reflective;
    options.max_speed = 10.0;

    let mut flock = Flock::new(&options);

    for _ in 0..200 {
        flock.update(&mut options);

        // Boids should stay within bounds
        for (boid, _) in flock.view2() {
            // Allow small tolerance for numerical precision
            assert!(
                boid.position.x >= options.window.win_left as f32 - 0.1
                    && boid.position.x <= options.window.win_right as f32 + 0.1,
                "Boid x position {} outside bounds",
                boid.position.x
            );
            assert!(
                boid.position.y >= options.window.win_bottom as f32 - 0.1
                    && boid.position.y <= options.window.win_top as f32 + 0.1,
                "Boid y position {} outside bounds",
                boid.position.y
            );
        }
    }
}

/// Test that simulation runs without panicking for edge cases
#[test]
fn test_simulation_with_minimal_boids() {
    let mut options = RunOptions::default();
    options.init_boids = 1; // Single boid
    options.window = get_window_size(100, 100);

    let mut flock = Flock::new(&options);

    // Should not panic with single boid
    for _ in 0..100 {
        flock.update(&mut options);
    }

    assert_eq!(flock.view2().count(), 1);
}

#[test]
fn test_simulation_with_many_boids() {
    let mut options = RunOptions::default();
    options.init_boids = 1000;
    options.window = get_window_size(2000, 2000);

    let mut flock = Flock::new(&options);

    // Should handle large numbers without issues
    for _ in 0..10 {
        flock.update(&mut options);
    }

    assert_eq!(flock.view2().count(), 1000);
}

/// Test that restart properly changes boid count
#[test]
fn test_restart_changes_boid_count() {
    let mut options = RunOptions::default();
    options.init_boids = 200;
    options.window = get_window_size(800, 600);

    let mut flock = Flock::new(&options);
    assert_eq!(flock.view2().count(), 200, "Should start with 200 boids");

    // Add more boids (simulate doubling)
    for _ in 0..200 {
        flock.insert_single(&options);
    }
    assert_eq!(flock.view2().count(), 400, "Should have 400 after doubling");

    // Now restart with init_boids = 200
    options.init_boids = 200;
    flock.restart(&options);
    assert_eq!(flock.view2().count(), 200, "Should have 200 after restart");
}
