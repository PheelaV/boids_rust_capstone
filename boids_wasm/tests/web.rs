//! Test suite for the Web and headless browsers.

#![cfg(target_arch = "wasm32")]

extern crate wasm_bindgen_test;
use wasm_bindgen_test::*;

use boids_wasm::*;
use wasm_bindgen::JsValue;

wasm_bindgen_test_configure!(run_in_browser);

#[wasm_bindgen_test]
fn test_simulation_creation() {
    let sim = WasmSimulation::new();
    assert!(sim.is_ok(), "Should be able to create a simulation");

    let sim = sim.unwrap();
    assert_eq!(sim.get_frame_count(), 0, "Frame count should start at 0");
    assert!(sim.get_boid_count() > 0, "Should have boids");
}

#[wasm_bindgen_test]
fn test_simulation_update() {
    let mut sim = WasmSimulation::new().unwrap();
    let initial_frame = sim.get_frame_count();

    sim.update();

    assert_eq!(
        sim.get_frame_count(),
        initial_frame + 1,
        "Frame count should increment after update"
    );
}

#[wasm_bindgen_test]
fn test_get_boids() {
    let sim = WasmSimulation::new().unwrap();
    let boids_result = sim.get_boids();

    assert!(boids_result.is_ok(), "Should be able to get boids");

    // Test that we can serialize and deserialize boids
    let boids_js = boids_result.unwrap();
    assert!(!boids_js.is_null(), "Boids data should not be null");
}

#[wasm_bindgen_test]
fn test_toggle_behaviors() {
    let mut sim = WasmSimulation::new().unwrap();

    // Test toggling behaviors doesn't crash
    sim.toggle_separation();
    sim.toggle_cohesion();
    sim.toggle_alignment();
    sim.toggle_wander();

    // Should still be able to update after toggling
    sim.update();
}

#[wasm_bindgen_test]
fn test_set_coefficients() {
    let mut sim = WasmSimulation::new().unwrap();

    sim.set_separation_coefficient(2.0);
    sim.set_cohesion_coefficient(1.5);
    sim.set_alignment_coefficient(1.2);
    sim.set_max_speed(5.0);

    // Should be able to update after setting coefficients
    sim.update();
}

#[wasm_bindgen_test]
fn test_reset() {
    let mut sim = WasmSimulation::new().unwrap();

    // Run a few frames
    for _ in 0..10 {
        sim.update();
    }

    assert_eq!(sim.get_frame_count(), 10, "Should have 10 frames");

    // Reset
    sim.reset();

    assert_eq!(sim.get_frame_count(), 0, "Frame count should reset to 0");
}

#[wasm_bindgen_test]
fn test_get_stats() {
    let sim = WasmSimulation::new().unwrap();
    let stats_result = sim.get_stats();

    assert!(stats_result.is_ok(), "Should be able to get stats");

    let stats = stats_result.unwrap();
    assert!(!stats.is_null(), "Stats should not be null");
}

#[wasm_bindgen_test]
fn test_custom_config() {
    use serde_json::json;

    let config = json!({
        "init_boids": 100,
        "window_width": 1024,
        "window_height": 768,
        "separation_coefficient": 2.0,
        "cohesion_coefficient": 1.5,
        "alignment_coefficient": 1.0,
        "max_speed": 5.0,
        "min_speed": 1.0,
        "sensory_distance": 80.0,
        "wander_on": true,
        "wander_coefficient": 0.8
    });

    let config_js = serde_wasm_bindgen::to_value(&config).unwrap();
    let sim = WasmSimulation::from_config(config_js);

    assert!(sim.is_ok(), "Should create simulation from custom config");

    let sim = sim.unwrap();
    assert_eq!(sim.get_boid_count(), 100, "Should have 100 boids");
}
