//! Python enum wrappers for boids_lib types.

use boids_lib::options::{
    Boundary, Distance, InitiationStrategy, NeighbourSampling, NoiseModel, TrackerType,
};
use pyo3::prelude::*;

/// Boundary condition type for the simulation space.
#[pyclass(eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum PyBoundary {
    /// Wrap-around boundaries (agents exiting one side appear on the opposite).
    Toroidal = 0,
    /// Agents stop at boundaries.
    Absorbing = 1,
    /// Agents bounce off boundaries.
    Reflective = 2,
}

#[pymethods]
impl PyBoundary {
    fn __repr__(&self) -> &'static str {
        match self {
            PyBoundary::Toroidal => "Boundary.Toroidal",
            PyBoundary::Absorbing => "Boundary.Absorbing",
            PyBoundary::Reflective => "Boundary.Reflective",
        }
    }
}

impl PyBoundary {
    pub fn to_rust(&self) -> Boundary {
        match self {
            PyBoundary::Toroidal => Boundary::Toroidal,
            PyBoundary::Absorbing => Boundary::Absorbing,
            PyBoundary::Reflective => Boundary::Reflective,
        }
    }

    pub fn from_rust(b: &Boundary) -> Self {
        match b {
            Boundary::Toroidal => PyBoundary::Toroidal,
            Boundary::Absorbing => PyBoundary::Absorbing,
            Boundary::Reflective => PyBoundary::Reflective,
            // Repulsive variants map to Reflective for enum purposes
            Boundary::Repulsive { .. } => PyBoundary::Reflective,
            Boundary::RepulsiveCircle { .. } => PyBoundary::Reflective,
        }
    }
}

/// Repulsive boundary configuration (agents are pushed away from edges).
#[pyclass]
#[derive(Clone, Debug)]
pub struct PyBoundaryRepulsive {
    #[pyo3(get, set)]
    pub distance: f32,
    #[pyo3(get, set)]
    pub force: f32,
}

#[pymethods]
impl PyBoundaryRepulsive {
    #[new]
    pub fn new(distance: f32, force: f32) -> Self {
        Self { distance, force }
    }

    fn __repr__(&self) -> String {
        format!(
            "BoundaryRepulsive(distance={}, force={})",
            self.distance, self.force
        )
    }
}

impl PyBoundaryRepulsive {
    pub fn to_rust(&self) -> Boundary {
        Boundary::Repulsive {
            distance: self.distance,
            force: self.force,
        }
    }
}

/// Repulsive circle boundary configuration.
#[pyclass]
#[derive(Clone, Debug)]
pub struct PyBoundaryRepulsiveCircle {
    #[pyo3(get, set)]
    pub radius: f32,
}

#[pymethods]
impl PyBoundaryRepulsiveCircle {
    #[new]
    pub fn new(radius: f32) -> Self {
        Self { radius }
    }

    fn __repr__(&self) -> String {
        format!("BoundaryRepulsiveCircle(radius={})", self.radius)
    }
}

impl PyBoundaryRepulsiveCircle {
    pub fn to_rust(&self) -> Boundary {
        Boundary::RepulsiveCircle {
            radius: self.radius,
        }
    }
}

/// Distance calculation method.
#[pyclass(eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum PyDistance {
    /// Euclidean distance with toroidal wrap-around.
    EucToroidal = 0,
    /// Standard Euclidean distance in enclosed space.
    EucEnclosed = 1,
}

#[pymethods]
impl PyDistance {
    fn __repr__(&self) -> &'static str {
        match self {
            PyDistance::EucToroidal => "Distance.EucToroidal",
            PyDistance::EucEnclosed => "Distance.EucEnclosed",
        }
    }
}

impl PyDistance {
    pub fn to_rust(&self) -> Distance {
        match self {
            PyDistance::EucToroidal => Distance::EucToroidal,
            PyDistance::EucEnclosed => Distance::EucEnclosed,
        }
    }

    pub fn from_rust(d: &Distance) -> Self {
        match d {
            Distance::EucToroidal => PyDistance::EucToroidal,
            Distance::EucEnclosed => PyDistance::EucEnclosed,
        }
    }
}

/// Noise/wander model for agent movement.
#[pyclass(eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum PyNoiseModel {
    /// Vicsek-style random noise.
    Vicsek = 0,
    /// Reynolds-style wander behavior.
    Reynolds = 1,
}

#[pymethods]
impl PyNoiseModel {
    fn __repr__(&self) -> &'static str {
        match self {
            PyNoiseModel::Vicsek => "NoiseModel.Vicsek",
            PyNoiseModel::Reynolds => "NoiseModel.Reynolds",
        }
    }
}

impl PyNoiseModel {
    pub fn to_rust(&self) -> NoiseModel {
        match self {
            PyNoiseModel::Vicsek => NoiseModel::Vicsek,
            PyNoiseModel::Reynolds => NoiseModel::Reynolds,
        }
    }

    pub fn from_rust(n: &NoiseModel) -> Self {
        match n {
            NoiseModel::Vicsek => PyNoiseModel::Vicsek,
            NoiseModel::Reynolds => PyNoiseModel::Reynolds,
        }
    }
}

/// Initial positioning strategy for boids.
#[pyclass(eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum PyInitiationStrategy {
    /// Boids start at center and move outward.
    CircleCenterOut = 0,
    /// Boids start at circumference and move inward.
    CircleCircumferenceIn = 1,
    /// Boids start at two opposite walls.
    TwoWalls = 2,
    /// Boids start within a rectangle.
    RectangleIn = 3,
    /// Boids start at random positions with aligned velocities.
    RandomIn = 4,
    /// Boids start at random positions with random velocities.
    RandomRandom = 5,
}

#[pymethods]
impl PyInitiationStrategy {
    fn __repr__(&self) -> &'static str {
        match self {
            PyInitiationStrategy::CircleCenterOut => "InitiationStrategy.CircleCenterOut",
            PyInitiationStrategy::CircleCircumferenceIn => "InitiationStrategy.CircleCircumferenceIn",
            PyInitiationStrategy::TwoWalls => "InitiationStrategy.TwoWalls",
            PyInitiationStrategy::RectangleIn => "InitiationStrategy.RectangleIn",
            PyInitiationStrategy::RandomIn => "InitiationStrategy.RandomIn",
            PyInitiationStrategy::RandomRandom => "InitiationStrategy.RandomRandom",
        }
    }
}

impl PyInitiationStrategy {
    pub fn to_rust(&self) -> InitiationStrategy {
        match self {
            PyInitiationStrategy::CircleCenterOut => InitiationStrategy::CircleCenterOut,
            PyInitiationStrategy::CircleCircumferenceIn => InitiationStrategy::CircleCircumferenceIn,
            PyInitiationStrategy::TwoWalls => InitiationStrategy::TwoWalls,
            PyInitiationStrategy::RectangleIn => InitiationStrategy::RectangleIn,
            PyInitiationStrategy::RandomIn => InitiationStrategy::RandomIn,
            PyInitiationStrategy::RandomRandom => InitiationStrategy::RandomRandom,
        }
    }

    pub fn from_rust(s: &InitiationStrategy) -> Self {
        match s {
            InitiationStrategy::CircleCenterOut => PyInitiationStrategy::CircleCenterOut,
            InitiationStrategy::CircleCircumferenceIn => PyInitiationStrategy::CircleCircumferenceIn,
            InitiationStrategy::TwoWalls => PyInitiationStrategy::TwoWalls,
            InitiationStrategy::RectangleIn => PyInitiationStrategy::RectangleIn,
            InitiationStrategy::RandomIn => PyInitiationStrategy::RandomIn,
            InitiationStrategy::RandomRandom => PyInitiationStrategy::RandomRandom,
        }
    }
}

/// Spatial tracking algorithm type.
#[pyclass(eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum PyTrackerType {
    /// Spatial hashing (O(n) average for neighbor queries).
    SpatHash = 0,
    /// Naive O(n^2) neighbor search.
    Naive = 1,
}

#[pymethods]
impl PyTrackerType {
    fn __repr__(&self) -> &'static str {
        match self {
            PyTrackerType::SpatHash => "TrackerType.SpatHash",
            PyTrackerType::Naive => "TrackerType.Naive",
        }
    }
}

impl PyTrackerType {
    pub fn to_rust(&self) -> TrackerType {
        match self {
            PyTrackerType::SpatHash => TrackerType::SpatHash,
            PyTrackerType::Naive => TrackerType::Naive,
        }
    }

    pub fn from_rust(t: &TrackerType) -> Self {
        match t {
            TrackerType::SpatHash => PyTrackerType::SpatHash,
            TrackerType::Naive => PyTrackerType::Naive,
            TrackerType::Replay(_, _) => PyTrackerType::SpatHash, // Default for replay
        }
    }
}

/// Neighbor sampling strategy.
#[pyclass(eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum PyNeighbourSampling {
    /// Biased sampling (early exit, faster but biased toward first cells).
    Biased = 0,
    /// Uniform strided sampling across all neighboring cells.
    Uniform = 1,
}

#[pymethods]
impl PyNeighbourSampling {
    fn __repr__(&self) -> &'static str {
        match self {
            PyNeighbourSampling::Biased => "NeighbourSampling.Biased",
            PyNeighbourSampling::Uniform => "NeighbourSampling.Uniform",
        }
    }
}

impl PyNeighbourSampling {
    pub fn to_rust(&self) -> NeighbourSampling {
        match self {
            PyNeighbourSampling::Biased => NeighbourSampling::Biased,
            PyNeighbourSampling::Uniform => NeighbourSampling::Uniform,
        }
    }

    pub fn from_rust(n: &NeighbourSampling) -> Self {
        match n {
            NeighbourSampling::Biased => PyNeighbourSampling::Biased,
            NeighbourSampling::Uniform => PyNeighbourSampling::Uniform,
        }
    }
}
