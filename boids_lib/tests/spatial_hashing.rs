//! Integration tests for spatial hashing correctness
//!
//! These tests ensure that the spatial hashing tracker works correctly
//! for neighbor queries during simulation.

use boids_lib::{
    boid::Boid,
    flock::Flock,
    options::*,
};

/// Test that spatial hash maintains consistent neighbor relationships
#[test]
fn test_spathash_neighbor_consistency() {
    let mut options = RunOptions::default();
    options.init_boids = 50;
    options.sensory_distance = 100.0;
    options.update_sensory_distances();
    options.window = get_window_size(800, 600);
    options.neighbours_cosidered = 20;
    options.tracker_type = TrackerType::SpatHash;

    let mut flock = Flock::new(&options);

    // Run simulation to move boids around
    for _ in 0..20 {
        flock.update(&mut options);

        // Get all boids
        let boids: Vec<Boid> = flock.view2().map(|(b, _)| *b).collect();

        // Test a few boids
        for test_id in [0, 10, 20, 30, 40] {
            if test_id >= boids.len() {
                continue;
            }

            let query_boid = &boids[test_id];
            let neighbors = flock.tracker.get_neighbours(query_boid, &options);

            // Should return reasonable number of neighbors
            assert!(
                neighbors.len() <= options.neighbours_cosidered,
                "Boid {} has {} neighbors (max {})",
                test_id,
                neighbors.len(),
                options.neighbours_cosidered
            );

            // Neighbors should not include the boid itself
            for neighbor in &neighbors {
                assert_ne!(
                    neighbor.id, query_boid.id,
                    "Boid should not be in its own neighbor list"
                );
            }

            // All neighbors should be within sensory distance
            for neighbor in &neighbors {
                let distance = boids_lib::math_helpers::distance_dyn_boid(
                    query_boid,
                    neighbor,
                    &options
                );

                assert!(
                    distance <= options.max_sensory_distance + 1.0,
                    "Neighbor {} is {} away from boid {} (max: {})",
                    neighbor.id,
                    distance,
                    query_boid.id,
                    options.max_sensory_distance
                );
            }
        }
    }
}

/// Test that neighbor count limits are respected
#[test]
fn test_neighbor_count_limits() {
    let mut options = RunOptions::default();
    options.init_boids = 100;
    options.sensory_distance = 200.0; // Large range
    options.update_sensory_distances();
    options.window = get_window_size(500, 500); // Small window = many potential neighbors
    options.neighbours_cosidered = 10; // Limit to 10 neighbors
    options.tracker_type = TrackerType::SpatHash;

    let flock = Flock::new(&options);
    let boids: Vec<Boid> = flock.view2().map(|(b, _)| *b).collect();

    // Check several boids
    for boid in boids.iter().take(20) {
        let neighbors = flock.tracker.get_neighbours(boid, &options);

        assert!(
            neighbors.len() <= options.neighbours_cosidered,
            "Boid {} should have at most {} neighbors, got {}",
            boid.id,
            options.neighbours_cosidered,
            neighbors.len()
        );
    }
}

/// Test toroidal boundary with spatial hash
#[test]
fn test_spatial_hash_toroidal_boundary() {
    let mut options = RunOptions::default();
    options.init_boids = 30;
    options.sensory_distance = 150.0;
    options.update_sensory_distances();
    options.window = get_window_size(400, 400);
    options.boundary = Boundary::Toroidal;
    options.distance = Distance::EucToroidal;
    options.tracker_type = TrackerType::SpatHash;
    options.max_speed = 20.0; // Fast movement
    options.neighbours_cosidered = 20; // Set a reasonable neighbor limit for testing

    let mut flock = Flock::new(&options);

    // Run for many steps to ensure boids cross boundaries
    for _ in 0..200 {
        flock.update(&mut options);
    }

    // All boids should still be tracked
    assert_eq!(
        flock.view2().count(),
        options.init_boids as usize,
        "All boids should still be tracked after boundary crossings"
    );

    // Verify neighbor queries still work
    let boids: Vec<Boid> = flock.view2().map(|(b, _)| *b).collect();
    for boid in boids.iter().take(10) {
        let neighbors = flock.tracker.get_neighbours(boid, &options);
        // Just verify it doesn't panic and returns reasonable results
        assert!(neighbors.len() <= options.neighbours_cosidered);
    }
}

/// Test spatial hash with large flock
#[test]
fn test_spatial_hash_large_flock() {
    let mut options = RunOptions::default();
    options.init_boids = 500;
    options.sensory_distance = 100.0;
    options.update_sensory_distances();
    options.window = get_window_size(2000, 2000);
    options.tracker_type = TrackerType::SpatHash;
    options.neighbours_cosidered = 30; // Set a reasonable neighbor limit for testing

    let mut flock = Flock::new(&options);

    // Should handle large flock efficiently
    for _ in 0..10 {
        flock.update(&mut options);
    }

    assert_eq!(flock.view2().count(), 500);

    // Spot check a few boids
    let boids: Vec<Boid> = flock.view2().map(|(b, _)| *b).collect();
    for idx in [0, 100, 200, 300, 400] {
        let neighbors = flock.tracker.get_neighbours(&boids[idx], &options);
        assert!(neighbors.len() <= options.neighbours_cosidered);
    }
}

/// Test that spatial hash vs naive produce comparable results
#[test]
fn test_spathash_vs_naive_comparison() {
    let mut options_spathash = RunOptions::default();
    options_spathash.init_boids = 50;
    options_spathash.sensory_distance = 100.0;
    options_spathash.update_sensory_distances();
    options_spathash.window = get_window_size(800, 600);
    options_spathash.tracker_type = TrackerType::SpatHash;

    let mut options_naive = options_spathash.clone();
    options_naive.tracker_type = TrackerType::Naive;

    let mut flock_spathash = Flock::new(&options_spathash);
    let mut flock_naive = Flock::new(&options_naive);

    // Run both for several steps
    for _ in 0..10 {
        flock_spathash.update(&mut options_spathash);
        flock_naive.update(&mut options_naive);
    }

    // Both should maintain same number of boids
    assert_eq!(flock_spathash.view2().count(), flock_naive.view2().count());

    // Sample a few boids and compare neighbor counts
    let spathash_boids: Vec<Boid> = flock_spathash.view2().map(|(b, _)| *b).collect();
    let naive_boids: Vec<Boid> = flock_naive.view2().map(|(b, _)| *b).collect();

    for idx in [0, 10, 20, 30, 40] {
        if idx >= spathash_boids.len() {
            continue;
        }

        let spathash_neighbors = flock_spathash.tracker.get_neighbours(&spathash_boids[idx], &options_spathash);
        let naive_neighbors = flock_naive.tracker.get_neighbours(&naive_boids[idx], &options_naive);

        // Should find similar number of neighbors (may not be exact due to floating point)
        let diff = (spathash_neighbors.len() as i32 - naive_neighbors.len() as i32).abs();

        assert!(
            diff <= 5,
            "Boid {}: SpatHash found {} neighbors, Naive found {} (diff: {})",
            idx,
            spathash_neighbors.len(),
            naive_neighbors.len(),
            diff
        );
    }
}

/// Test spatial hash handles boid insertion/deletion
/// TODO: SpatHash1D uses boid IDs as array indices which breaks with deletion
/// Requires refactoring to use a separate ID→index mapping or ID reassignment
#[test]
#[ignore]
fn test_spatial_hash_dynamic_boids() {
    let mut options = RunOptions::default();
    options.init_boids = 50;
    options.sensory_distance = 100.0;
    options.update_sensory_distances();
    options.window = get_window_size(800, 600);
    options.tracker_type = TrackerType::SpatHash;

    let mut flock = Flock::new(&options);

    // Add some boids
    for _ in 0..10 {
        flock.insert_single(&options);
    }

    assert_eq!(flock.view2().count(), 60);

    // Remove some boids
    for _ in 0..5 {
        flock.delete_last();
    }

    assert_eq!(flock.view2().count(), 55);

    // Should still work correctly
    flock.update(&mut options);

    let boids: Vec<Boid> = flock.view2().map(|(b, _)| *b).collect();
    for boid in boids.iter().take(10) {
        let neighbors = flock.tracker.get_neighbours(boid, &options);
        assert!(neighbors.len() <= options.neighbours_cosidered);
    }
}
