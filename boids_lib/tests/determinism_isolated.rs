//! Isolated determinism tests - runs as separate test binary
//!
//! These tests are isolated because they require exclusive access to the global RNG.
//! Running them in their own binary prevents interference from other tests.

use approx::assert_relative_eq;
use boids_lib::{flock_base, flock::Flock, options::*};
use glam::Vec2;

/// Test that simulation with same seed produces identical results.
#[test]
fn test_deterministic_simulation() {
    let mut options1 = RunOptions::default();
    options1.init_boids = 50;
    options1.window = get_window_size(800, 600);
    options1.rng_seed = Some(42);

    let options2 = options1.clone();

    // Run simulation twice with same settings and seed
    let result1 = flock_base(100, options1);
    let result2 = flock_base(100, options2);

    // Results should be identical for deterministic behavior
    assert_eq!(result1.len(), result2.len(), "Same number of data points");

    // Compare a few sample points
    if !result1.is_empty() {
        let sample1 = &result1[0];
        let sample2 = &result2[0];

        assert_eq!(sample1.id, sample2.id, "Boid IDs match");
        assert_relative_eq!(sample1.x, sample2.x, epsilon = 0.0001);
        assert_relative_eq!(sample1.y, sample2.y, epsilon = 0.0001);
    }
}

/// Helper to calculate average distance between all pairs of boids
fn calculate_average_neighbor_distance(flock: &Flock) -> f32 {
    let positions: Vec<Vec2> = flock.view2().map(|(b, _)| b.position).collect();
    if positions.len() < 2 {
        return 0.0;
    }

    let mut total_distance = 0.0;
    let mut count = 0;

    for i in 0..positions.len() {
        for j in (i + 1)..positions.len() {
            total_distance += (positions[i] - positions[j]).length();
            count += 1;
        }
    }

    total_distance / count as f32
}

/// Test that separation coefficient affects dispersion.
/// This test requires determinism to ensure both flocks start identically.
#[test]
fn test_coefficient_effects() {
    let mut options = RunOptions::default();
    options.init_boids = 20;
    options.window = get_window_size(800, 600);
    options.rng_seed = Some(12345);

    // Weak separation
    options.separation_coefficient = 0.1;
    options.cohesion_coefficient = 0.0;
    options.alignment_coefficient = 0.0;

    let mut flock_weak = Flock::new(&options);
    for _ in 0..50 {
        flock_weak.update(&mut options);
    }
    let dispersion_weak = calculate_average_neighbor_distance(&flock_weak);

    // Strong separation - reseed to get identical starting positions
    options.separation_coefficient = 5.0;
    let mut flock_strong = Flock::new(&options);
    for _ in 0..50 {
        flock_strong.update(&mut options);
    }
    let dispersion_strong = calculate_average_neighbor_distance(&flock_strong);

    // Stronger separation should lead to greater average distance
    assert!(
        dispersion_strong > dispersion_weak,
        "Strong separation ({}) should create more dispersion than weak ({})",
        dispersion_strong,
        dispersion_weak
    );
}
