//! Integration tests for individual boid behavior rules
//!
//! These tests validate that separation, cohesion, alignment, and wander
//! behaviors produce expected effects on boid movement.

use approx::assert_relative_eq;
use boids_lib::{boid::*, flock::Flock, options::*};
use glam::Vec2;

/// Test that separation pushes boids apart when too close
#[test]
fn test_separation_behavior() {
    let mut options = RunOptions::default();
    options.init_boids = 2;
    options.window = get_window_size(1000, 1000);

    // Enable only separation
    options.separation_on = true;
    options.cohesion_on = false;
    options.alignment_on = false;
    options.wander_on = false;

    options.separation_coefficient = 2.0;
    options.separation_treshold_distance = 100.0;

    // Create two boids very close together
    let mut boid1 = Boid::new(0.0, 0.0, Vec2::new(0.0, 0.0), 0);
    let boid2 = Boid::new(10.0, 0.0, Vec2::new(0.0, 0.0), 1);

    let boids = vec![&boid1, &boid2];
    let metadata = vec![BoidMetadata::new(&boid1), BoidMetadata::new(&boid2)];

    // Calculate separation force
    let force = boid1.run_rules(&boids, &metadata, &options);

    // Force should push boid1 away from boid2 (negative x direction)
    assert!(
        force.x < 0.0,
        "Separation should push boid1 left (away from boid2), got force.x = {}",
        force.x
    );
}

/// Test that cohesion pulls boids together when within range
#[test]
fn test_cohesion_behavior() {
    let mut options = RunOptions::default();
    options.init_boids = 3;
    options.window = get_window_size(1000, 1000);

    // Enable only cohesion
    options.separation_on = false;
    options.cohesion_on = true;
    options.alignment_on = false;
    options.wander_on = false;

    options.cohesion_coefficient = 1.0;
    options.cohesion_treshold_distance = 200.0;

    // Create a boid offset from center of mass
    let boid1 = Boid::new(-50.0, 0.0, Vec2::new(0.0, 0.0), 0);
    let boid2 = Boid::new(50.0, 0.0, Vec2::new(0.0, 0.0), 1);
    let boid3 = Boid::new(50.0, 50.0, Vec2::new(0.0, 0.0), 2);

    let boids = vec![&boid1, &boid2, &boid3];
    let metadata = vec![
        BoidMetadata::new(&boid1),
        BoidMetadata::new(&boid2),
        BoidMetadata::new(&boid3),
    ];

    // Calculate cohesion force for boid1
    let force = boid1.run_rules(&boids, &metadata, &options);

    // Force should pull boid1 towards the center of mass (positive x direction)
    assert!(
        force.x > 0.0,
        "Cohesion should pull boid1 right (towards center), got force.x = {}",
        force.x
    );
}

/// Test that alignment makes boids match velocity with neighbors
#[test]
fn test_alignment_behavior() {
    let mut options = RunOptions::default();
    options.init_boids = 3;
    options.window = get_window_size(1000, 1000);

    // Enable only alignment
    options.separation_on = false;
    options.cohesion_on = false;
    options.alignment_on = true;
    options.wander_on = false;

    options.alignment_coefficient = 1.0;
    options.alignment_treshold_distance = 150.0;
    options.agent_steering = true;

    // Create boids with different velocities
    let boid1 = Boid::new(0.0, 0.0, Vec2::new(1.0, 0.0), 0);
    let boid2 = Boid::new(50.0, 0.0, Vec2::new(0.0, 3.0), 1); // Moving up
    let boid3 = Boid::new(100.0, 0.0, Vec2::new(0.0, 3.0), 2); // Moving up

    let boids = vec![&boid1, &boid2, &boid3];
    let metadata = vec![
        BoidMetadata::new(&boid1),
        BoidMetadata::new(&boid2),
        BoidMetadata::new(&boid3),
    ];

    // Calculate alignment force for boid1
    let force = boid1.run_rules(&boids, &metadata, &options);

    // Force should have positive y component to align with neighbors moving up
    assert!(
        force.y > 0.0,
        "Alignment should steer boid1 upward to match neighbors, got force.y = {}",
        force.y
    );
}

/// Test that individual rules can be toggled on/off
#[test]
fn test_rules_toggle() {
    let mut options = RunOptions::default();
    options.init_boids = 10;
    options.window = get_window_size(500, 500);

    let mut flock = Flock::new(&options);

    // Capture state with all rules on
    options.separation_on = true;
    options.cohesion_on = true;
    options.alignment_on = true;

    for _ in 0..5 {
        flock.update(&mut options);
    }

    let positions_all_rules: Vec<Vec2> = flock.view2().map(|(b, _)| b.position).collect();

    // Reset and run with all rules off
    flock.restart(&options);

    options.separation_on = false;
    options.cohesion_on = false;
    options.alignment_on = false;

    for _ in 0..5 {
        flock.update(&mut options);
    }

    let positions_no_rules: Vec<Vec2> = flock.view2().map(|(b, _)| b.position).collect();

    // Results should be different (boids behave differently)
    let mut differences = 0;
    for (pos1, pos2) in positions_all_rules.iter().zip(positions_no_rules.iter()) {
        if (*pos1 - *pos2).length() > 1.0 {
            differences += 1;
        }
    }

    assert!(
        differences > 0,
        "Toggling rules should produce different behavior"
    );
}

/// Test that coefficients affect behavior strength
#[test]
fn test_coefficient_effects() {
    let mut options = RunOptions::default();
    options.init_boids = 20;
    options.window = get_window_size(800, 600);

    // Weak separation
    options.separation_coefficient = 0.1;
    options.cohesion_coefficient = 0.0;
    options.alignment_coefficient = 0.0;

    let mut flock_weak = Flock::new(&options);
    for _ in 0..50 {
        flock_weak.update(&mut options);
    }

    let dispersion_weak = calculate_average_neighbor_distance(&flock_weak);

    // Strong separation
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

/// Test field of vision constrains neighbor detection
#[test]
fn test_field_of_vision() {
    let mut options = RunOptions::default();
    options.window = get_window_size(1000, 1000);
    options.field_of_vision_on = true;
    options.field_of_vision_deg = 90.0; // 90 degree cone
    options.update_fov();

    // Create a boid facing right with neighbors in different positions
    let boid = Boid::new(0.0, 0.0, Vec2::new(1.0, 0.0), 0); // Facing right (+x)
    let neighbor_front = Boid::new(100.0, 0.0, Vec2::ZERO, 1); // In front (visible)
    let neighbor_behind = Boid::new(-100.0, 0.0, Vec2::ZERO, 2); // Behind (not visible)

    let all_neighbors = vec![&boid, &neighbor_front, &neighbor_behind];

    // Filter by field of vision
    let visible = boid.filter_sight2(&all_neighbors, &options);

    // Should only see the one in front
    assert_eq!(
        visible.len(),
        1,
        "With 90° FOV facing right, should only see front neighbor"
    );
    assert_eq!(
        visible[0].id, 1,
        "Should see the neighbor in front (id=1)"
    );
}

/// Test wander behavior adds variation
#[test]
fn test_wander_adds_variation() {
    let mut options = RunOptions::default();
    options.init_boids = 10;
    options.window = get_window_size(1000, 1000);

    // Disable other behaviors, enable wander
    options.separation_on = false;
    options.cohesion_on = false;
    options.alignment_on = false;
    options.wander_on = true;
    options.wander_coefficient = 1.0;
    options.wander_rate = 0.5;

    let mut flock = Flock::new(&options);

    // Record initial headings
    let initial_headings: Vec<Vec2> = flock.view2().map(|(b, _)| b.velocity.normalize()).collect();

    // Run simulation
    for _ in 0..100 {
        flock.update(&mut options);
    }

    // Record final headings
    let final_headings: Vec<Vec2> = flock.view2().map(|(b, _)| b.velocity.normalize()).collect();

    // Headings should have changed due to wander
    let mut heading_changes = 0;
    for (initial, final_h) in initial_headings.iter().zip(final_headings.iter()) {
        let angle_change = initial.dot(*final_h).acos();
        if angle_change > 0.1 {
            heading_changes += 1;
        }
    }

    assert!(
        heading_changes >= 5,
        "Wander should cause heading changes in multiple boids, got {} changes",
        heading_changes
    );
}

// Helper function
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
