//! Python type wrappers for boids_lib structs.

use boids_lib::birdwatcher::BoidData;
use boids_lib::options::{get_window_size, RunOptions, WindowSize};
use pyo3::prelude::*;

use crate::enums::{
    PyBoundary, PyBoundaryRepulsive, PyBoundaryRepulsiveCircle, PyDistance, PyInitiationStrategy,
    PyNeighbourSampling, PyNoiseModel, PyTrackerType,
};

/// Window/simulation space dimensions.
#[pyclass]
#[derive(Clone, Debug)]
pub struct PyWindowSize {
    #[pyo3(get, set)]
    pub win_left: i32,
    #[pyo3(get, set)]
    pub win_right: i32,
    #[pyo3(get, set)]
    pub win_top: i32,
    #[pyo3(get, set)]
    pub win_bottom: i32,
    #[pyo3(get, set)]
    pub win_h: i32,
    #[pyo3(get, set)]
    pub win_w: i32,
}

#[pymethods]
impl PyWindowSize {
    /// Create a new WindowSize from width and height.
    #[new]
    #[pyo3(signature = (width, height))]
    fn new(width: u32, height: u32) -> Self {
        Self::from_rust(&get_window_size(width, height))
    }

    fn __repr__(&self) -> String {
        format!(
            "WindowSize(w={}, h={}, left={}, right={}, bottom={}, top={})",
            self.win_w, self.win_h, self.win_left, self.win_right, self.win_bottom, self.win_top
        )
    }
}

impl PyWindowSize {
    pub fn from_rust(ws: &WindowSize) -> Self {
        Self {
            win_left: ws.win_left,
            win_right: ws.win_right,
            win_top: ws.win_top,
            win_bottom: ws.win_bottom,
            win_h: ws.win_h,
            win_w: ws.win_w,
        }
    }

    pub fn to_rust(&self) -> WindowSize {
        WindowSize {
            win_left: self.win_left,
            win_right: self.win_right,
            win_top: self.win_top,
            win_bottom: self.win_bottom,
            win_h: self.win_h,
            win_w: self.win_w,
        }
    }
}

/// Recorded boid state data at a specific time step.
#[pyclass]
#[derive(Clone, Debug)]
pub struct PyBoidData {
    #[pyo3(get)]
    pub id: usize,
    #[pyo3(get)]
    pub x: f32,
    #[pyo3(get)]
    pub y: f32,
    #[pyo3(get)]
    pub cluster_id: usize,
    #[pyo3(get)]
    pub n_neighbours: usize,
    #[pyo3(get)]
    pub time: u64,
}

#[pymethods]
impl PyBoidData {
    fn __repr__(&self) -> String {
        format!(
            "BoidData(id={}, x={:.2}, y={:.2}, time={}, cluster_id={}, n_neighbours={})",
            self.id, self.x, self.y, self.time, self.cluster_id, self.n_neighbours
        )
    }
}

impl PyBoidData {
    pub fn from_rust(bd: &BoidData) -> Self {
        Self {
            id: bd.id,
            x: bd.x,
            y: bd.y,
            cluster_id: bd.cluster_id,
            n_neighbours: bd.n_neighbours,
            time: bd.time,
        }
    }
}

/// Main configuration object for boids simulation.
///
/// Contains all configurable parameters controlling the simulation behavior.
/// Use keyword arguments in the constructor for convenience, or set properties
/// individually for full control.
///
/// Example:
///     >>> opts = RunOptions(init_boids=128, width=800, height=600)
///     >>> opts.separation_coefficient = 5.0
///     >>> opts.rng_seed = 42
#[pyclass]
#[derive(Clone)]
pub struct PyRunOptions {
    pub(crate) inner: RunOptions,
    // Store repulsive boundary separately for Python access
    boundary_repulsive: Option<PyBoundaryRepulsive>,
    boundary_repulsive_circle: Option<PyBoundaryRepulsiveCircle>,
}

#[pymethods]
impl PyRunOptions {
    /// Create a new RunOptions with sensible defaults.
    ///
    /// All parameters are optional and will use the library defaults if not specified.
    ///
    /// Args:
    ///     init_boids: Number of boids to simulate (default: 256)
    ///     width: Simulation window width (default: 600)
    ///     height: Simulation window height (default: 600)
    ///     sample_rate: How often to sample boid positions (default: 1)
    ///     rng_seed: Random seed for reproducibility (default: None = random)
    ///     separation_coefficient: Strength of separation behavior (default: 4.1)
    ///     cohesion_coefficient: Strength of cohesion behavior (default: 0.002)
    ///     alignment_coefficient: Strength of alignment behavior (default: 0.02)
    ///     min_speed: Minimum boid speed (default: 0.65)
    ///     max_speed: Maximum boid speed (default: 4.1)
    ///     max_steering: Maximum steering force (default: 0.7)
    ///     sensory_distance: Base sensing radius (default: 60.0)
    ///     field_of_vision_deg: Field of vision in degrees (default: 181.0)
    ///     separation_on: Enable separation behavior (default: true)
    ///     cohesion_on: Enable cohesion behavior (default: true)
    ///     alignment_on: Enable alignment behavior (default: true)
    ///     wander_on: Enable wander behavior (default: false)
    ///     boundary: Boundary type (default: Toroidal)
    ///     distance: Distance calculation method (default: EucToroidal)
    ///     noise_model: Noise/wander model (default: Reynolds)
    ///     initiation_strategy: Initial boid placement (default: RandomRandom)
    ///     tracker_type: Spatial tracking algorithm (default: SpatHash)
    ///     neighbour_sampling: Neighbor sampling strategy (default: Uniform)
    #[new]
    #[pyo3(signature = (
        init_boids = 256,
        width = 600,
        height = 600,
        sample_rate = 1,
        rng_seed = None,
        separation_coefficient = None,
        cohesion_coefficient = None,
        alignment_coefficient = None,
        min_speed = None,
        max_speed = None,
        max_steering = None,
        sensory_distance = None,
        field_of_vision_deg = None,
        separation_on = None,
        cohesion_on = None,
        alignment_on = None,
        wander_on = None,
        boundary = None,
        distance = None,
        noise_model = None,
        initiation_strategy = None,
        tracker_type = None,
        neighbour_sampling = None,
    ))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        init_boids: usize,
        width: u32,
        height: u32,
        sample_rate: u16,
        rng_seed: Option<u64>,
        separation_coefficient: Option<f32>,
        cohesion_coefficient: Option<f32>,
        alignment_coefficient: Option<f32>,
        min_speed: Option<f32>,
        max_speed: Option<f32>,
        max_steering: Option<f32>,
        sensory_distance: Option<f32>,
        field_of_vision_deg: Option<f32>,
        separation_on: Option<bool>,
        cohesion_on: Option<bool>,
        alignment_on: Option<bool>,
        wander_on: Option<bool>,
        boundary: Option<PyBoundary>,
        distance: Option<PyDistance>,
        noise_model: Option<PyNoiseModel>,
        initiation_strategy: Option<PyInitiationStrategy>,
        tracker_type: Option<PyTrackerType>,
        neighbour_sampling: Option<PyNeighbourSampling>,
    ) -> Self {
        let mut inner = RunOptions::default();

        // Core parameters
        inner.init_boids = init_boids;
        inner.window = get_window_size(width, height);
        inner.sample_rate = sample_rate;
        inner.rng_seed = rng_seed;

        // Behavior coefficients
        if let Some(v) = separation_coefficient {
            inner.separation_coefficient = v;
        }
        if let Some(v) = cohesion_coefficient {
            inner.cohesion_coefficient = v;
        }
        if let Some(v) = alignment_coefficient {
            inner.alignment_coefficient = v;
        }

        // Speed parameters
        if let Some(v) = min_speed {
            inner.min_speed = v;
            inner.min_speed_sq = v.powi(2);
        }
        if let Some(v) = max_speed {
            inner.max_speed = v;
            inner.max_speed_sq = v.powi(2);
        }
        if let Some(v) = max_steering {
            inner.max_steering = v;
            inner.max_steering_sq = v.powi(2);
        }

        // Sensory parameters
        if let Some(v) = sensory_distance {
            inner.sensory_distance = v;
        }
        if let Some(v) = field_of_vision_deg {
            inner.field_of_vision_deg = v;
        }

        // Behavior toggles
        if let Some(v) = separation_on {
            inner.separation_on = v;
        }
        if let Some(v) = cohesion_on {
            inner.cohesion_on = v;
        }
        if let Some(v) = alignment_on {
            inner.alignment_on = v;
        }
        if let Some(v) = wander_on {
            inner.wander_on = v;
        }

        // Enums
        if let Some(v) = boundary {
            inner.boundary = v.to_rust();
        }
        if let Some(v) = distance {
            inner.distance = v.to_rust();
        }
        if let Some(v) = noise_model {
            inner.noise_model = v.to_rust();
        }
        if let Some(v) = initiation_strategy {
            inner.initiation_strat = v.to_rust();
        }
        if let Some(v) = tracker_type {
            inner.tracker_type = v.to_rust();
        }
        if let Some(v) = neighbour_sampling {
            inner.neighbour_sampling = v.to_rust();
        }

        // Update derived values
        inner.update_sensory_distances();
        inner.update_fov();

        Self {
            inner,
            boundary_repulsive: None,
            boundary_repulsive_circle: None,
        }
    }

    // ========== Core Properties ==========

    #[getter]
    fn init_boids(&self) -> usize {
        self.inner.init_boids
    }
    #[setter]
    fn set_init_boids(&mut self, v: usize) {
        self.inner.init_boids = v;
    }

    #[getter]
    fn sample_rate(&self) -> u16 {
        self.inner.sample_rate
    }
    #[setter]
    fn set_sample_rate(&mut self, v: u16) {
        self.inner.sample_rate = v;
    }

    #[getter]
    fn rng_seed(&self) -> Option<u64> {
        self.inner.rng_seed
    }
    #[setter]
    fn set_rng_seed(&mut self, v: Option<u64>) {
        self.inner.rng_seed = v;
    }

    // ========== Speed Properties ==========

    #[getter]
    fn baseline_speed(&self) -> f32 {
        self.inner.baseline_speed
    }
    #[setter]
    fn set_baseline_speed(&mut self, v: f32) {
        self.inner.baseline_speed = v;
    }

    #[getter]
    fn min_speed(&self) -> f32 {
        self.inner.min_speed
    }
    #[setter]
    fn set_min_speed(&mut self, v: f32) {
        self.inner.min_speed = v;
        self.inner.min_speed_sq = v.powi(2);
    }

    #[getter]
    fn max_speed(&self) -> f32 {
        self.inner.max_speed
    }
    #[setter]
    fn set_max_speed(&mut self, v: f32) {
        self.inner.max_speed = v;
        self.inner.max_speed_sq = v.powi(2);
    }

    #[getter]
    fn max_steering(&self) -> f32 {
        self.inner.max_steering
    }
    #[setter]
    fn set_max_steering(&mut self, v: f32) {
        self.inner.max_steering = v;
        self.inner.max_steering_sq = v.powi(2);
    }

    // ========== Behavior Coefficients ==========

    #[getter]
    fn separation_coefficient(&self) -> f32 {
        self.inner.separation_coefficient
    }
    #[setter]
    fn set_separation_coefficient(&mut self, v: f32) {
        self.inner.separation_coefficient = v;
    }

    #[getter]
    fn cohesion_coefficient(&self) -> f32 {
        self.inner.cohesion_coefficient
    }
    #[setter]
    fn set_cohesion_coefficient(&mut self, v: f32) {
        self.inner.cohesion_coefficient = v;
    }

    #[getter]
    fn alignment_coefficient(&self) -> f32 {
        self.inner.alignment_coefficient
    }
    #[setter]
    fn set_alignment_coefficient(&mut self, v: f32) {
        self.inner.alignment_coefficient = v;
    }

    // ========== Sensory Properties ==========

    #[getter]
    fn sensory_distance(&self) -> f32 {
        self.inner.sensory_distance
    }
    #[setter]
    fn set_sensory_distance(&mut self, v: f32) {
        self.inner.sensory_distance = v;
        self.inner.update_sensory_distances();
    }

    #[getter]
    fn alignment_treshold_coefficient(&self) -> f32 {
        self.inner.alignment_treshold_coefficient
    }
    #[setter]
    fn set_alignment_treshold_coefficient(&mut self, v: f32) {
        self.inner.alignment_treshold_coefficient = v;
        self.inner.update_sensory_distances();
    }

    #[getter]
    fn cohesion_treshold_coefficient(&self) -> f32 {
        self.inner.cohesion_treshold_coefficient
    }
    #[setter]
    fn set_cohesion_treshold_coefficient(&mut self, v: f32) {
        self.inner.cohesion_treshold_coefficient = v;
        self.inner.update_sensory_distances();
    }

    #[getter]
    fn separation_treshold_coefficient(&self) -> f32 {
        self.inner.separation_treshold_coefficient
    }
    #[setter]
    fn set_separation_treshold_coefficient(&mut self, v: f32) {
        self.inner.separation_treshold_coefficient = v;
        self.inner.update_sensory_distances();
    }

    // ========== Behavior Toggles ==========

    #[getter]
    fn separation_on(&self) -> bool {
        self.inner.separation_on
    }
    #[setter]
    fn set_separation_on(&mut self, v: bool) {
        self.inner.separation_on = v;
    }

    #[getter]
    fn cohesion_on(&self) -> bool {
        self.inner.cohesion_on
    }
    #[setter]
    fn set_cohesion_on(&mut self, v: bool) {
        self.inner.cohesion_on = v;
    }

    #[getter]
    fn alignment_on(&self) -> bool {
        self.inner.alignment_on
    }
    #[setter]
    fn set_alignment_on(&mut self, v: bool) {
        self.inner.alignment_on = v;
    }

    // ========== Wander Properties ==========

    #[getter]
    fn wander_on(&self) -> bool {
        self.inner.wander_on
    }
    #[setter]
    fn set_wander_on(&mut self, v: bool) {
        self.inner.wander_on = v;
    }

    #[getter]
    fn wander_coefficient(&self) -> f32 {
        self.inner.wander_coefficient
    }
    #[setter]
    fn set_wander_coefficient(&mut self, v: f32) {
        self.inner.wander_coefficient = v;
    }

    #[getter]
    fn wander_rate(&self) -> f32 {
        self.inner.wander_rate
    }
    #[setter]
    fn set_wander_rate(&mut self, v: f32) {
        self.inner.wander_rate = v;
    }

    #[getter]
    fn wander_radius(&self) -> f32 {
        self.inner.wander_radius
    }
    #[setter]
    fn set_wander_radius(&mut self, v: f32) {
        self.inner.wander_radius = v;
    }

    #[getter]
    fn wander_distance(&self) -> f32 {
        self.inner.wander_distance
    }
    #[setter]
    fn set_wander_distance(&mut self, v: f32) {
        self.inner.wander_distance = v;
    }

    // ========== Field of Vision Properties ==========

    #[getter]
    fn field_of_vision_on(&self) -> bool {
        self.inner.field_of_vision_on
    }
    #[setter]
    fn set_field_of_vision_on(&mut self, v: bool) {
        self.inner.field_of_vision_on = v;
    }

    #[getter]
    fn field_of_vision_deg(&self) -> f32 {
        self.inner.field_of_vision_deg
    }
    #[setter]
    fn set_field_of_vision_deg(&mut self, v: f32) {
        self.inner.field_of_vision_deg = v;
        self.inner.update_fov();
    }

    #[getter]
    fn alignment_fov_deg(&self) -> f32 {
        self.inner.alignment_fov_deg
    }
    #[setter]
    fn set_alignment_fov_deg(&mut self, v: f32) {
        self.inner.alignment_fov_deg = v;
        self.inner.update_fov();
    }

    #[getter]
    fn cohesion_fov_deg(&self) -> f32 {
        self.inner.cohesion_fov_deg
    }
    #[setter]
    fn set_cohesion_fov_deg(&mut self, v: f32) {
        self.inner.cohesion_fov_deg = v;
        self.inner.update_fov();
    }

    #[getter]
    fn separation_fov_deg(&self) -> f32 {
        self.inner.separation_fov_deg
    }
    #[setter]
    fn set_separation_fov_deg(&mut self, v: f32) {
        self.inner.separation_fov_deg = v;
        self.inner.update_fov();
    }

    // ========== Other Properties ==========

    #[getter]
    fn size(&self) -> f32 {
        self.inner.size
    }
    #[setter]
    fn set_size(&mut self, v: f32) {
        self.inner.size = v;
    }

    #[getter]
    fn dbscan_flock_clustering_on(&self) -> bool {
        self.inner.dbscan_flock_clustering_on
    }
    #[setter]
    fn set_dbscan_flock_clustering_on(&mut self, v: bool) {
        self.inner.dbscan_flock_clustering_on = v;
    }

    #[getter]
    fn neighbours_cosidered(&self) -> usize {
        self.inner.neighbours_cosidered
    }
    #[setter]
    fn set_neighbours_cosidered(&mut self, v: usize) {
        self.inner.neighbours_cosidered = v;
    }

    #[getter]
    fn agent_steering(&self) -> bool {
        self.inner.agent_steering
    }
    #[setter]
    fn set_agent_steering(&mut self, v: bool) {
        self.inner.agent_steering = v;
    }

    #[getter]
    fn stop_movement(&self) -> bool {
        self.inner.stop_movement
    }
    #[setter]
    fn set_stop_movement(&mut self, v: bool) {
        self.inner.stop_movement = v;
    }

    // ========== Window Property ==========

    #[getter]
    fn window(&self) -> PyWindowSize {
        PyWindowSize::from_rust(&self.inner.window)
    }
    #[setter]
    fn set_window(&mut self, v: PyWindowSize) {
        self.inner.window = v.to_rust();
    }

    /// Set window size using width and height.
    fn set_window_size(&mut self, width: u32, height: u32) {
        self.inner.window = get_window_size(width, height);
    }

    // ========== Enum Properties ==========

    #[getter]
    fn boundary(&self) -> PyBoundary {
        PyBoundary::from_rust(&self.inner.boundary)
    }
    #[setter]
    fn set_boundary(&mut self, v: PyBoundary) {
        self.inner.boundary = v.to_rust();
        self.boundary_repulsive = None;
        self.boundary_repulsive_circle = None;
    }

    /// Set a repulsive boundary that pushes boids away from edges.
    fn set_boundary_repulsive(&mut self, distance: f32, force: f32) {
        let repulsive = PyBoundaryRepulsive::new(distance, force);
        self.inner.boundary = repulsive.to_rust();
        self.boundary_repulsive = Some(repulsive);
        self.boundary_repulsive_circle = None;
    }

    /// Set a circular repulsive boundary.
    fn set_boundary_repulsive_circle(&mut self, radius: f32) {
        let repulsive = PyBoundaryRepulsiveCircle::new(radius);
        self.inner.boundary = repulsive.to_rust();
        self.boundary_repulsive_circle = Some(repulsive);
        self.boundary_repulsive = None;
    }

    #[getter]
    fn distance(&self) -> PyDistance {
        PyDistance::from_rust(&self.inner.distance)
    }
    #[setter]
    fn set_distance(&mut self, v: PyDistance) {
        self.inner.distance = v.to_rust();
    }

    #[getter]
    fn noise_model(&self) -> PyNoiseModel {
        PyNoiseModel::from_rust(&self.inner.noise_model)
    }
    #[setter]
    fn set_noise_model(&mut self, v: PyNoiseModel) {
        self.inner.noise_model = v.to_rust();
    }

    #[getter]
    fn initiation_strategy(&self) -> PyInitiationStrategy {
        PyInitiationStrategy::from_rust(&self.inner.initiation_strat)
    }
    #[setter]
    fn set_initiation_strategy(&mut self, v: PyInitiationStrategy) {
        self.inner.initiation_strat = v.to_rust();
    }

    #[getter]
    fn tracker_type(&self) -> PyTrackerType {
        PyTrackerType::from_rust(&self.inner.tracker_type)
    }
    #[setter]
    fn set_tracker_type(&mut self, v: PyTrackerType) {
        self.inner.tracker_type = v.to_rust();
    }

    #[getter]
    fn neighbour_sampling(&self) -> PyNeighbourSampling {
        PyNeighbourSampling::from_rust(&self.inner.neighbour_sampling)
    }
    #[setter]
    fn set_neighbour_sampling(&mut self, v: PyNeighbourSampling) {
        self.inner.neighbour_sampling = v.to_rust();
    }

    // ========== Repr ==========

    fn __repr__(&self) -> String {
        format!(
            "RunOptions(init_boids={}, window={}x{}, sensory_distance={:.1}, boundary={:?})",
            self.inner.init_boids,
            self.inner.window.win_w,
            self.inner.window.win_h,
            self.inner.sensory_distance,
            self.boundary()
        )
    }
}

impl PyRunOptions {
    /// Convert to the underlying Rust RunOptions.
    pub fn to_rust(&self) -> RunOptions {
        self.inner.clone()
    }

    /// Create from Rust RunOptions.
    pub fn from_rust(ro: RunOptions) -> Self {
        Self {
            inner: ro,
            boundary_repulsive: None,
            boundary_repulsive_circle: None,
        }
    }
}
