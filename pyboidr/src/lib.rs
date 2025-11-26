//! pyboidr - Python bindings for the boids simulation library.
//!
//! This crate provides Python bindings for boids_lib using PyO3, enabling
//! high-performance boid flocking simulations from Python.
//!
//! # Quick Start
//!
//! ```python
//! from pyboidr import run_simulation, RunOptions, Boundary
//!
//! # Create configuration
//! opts = RunOptions(
//!     init_boids=128,
//!     width=800,
//!     height=600,
//!     rng_seed=42,
//! )
//!
//! # Run simulation
//! df = run_simulation(1000, opts)
//! print(df.head())
//! ```

use pyo3::prelude::*;

mod enums;
mod simulation;
mod types;

use enums::{
    PyBoundary, PyBoundaryRepulsive, PyBoundaryRepulsiveCircle, PyDistance, PyInitiationStrategy,
    PyNeighbourSampling, PyNoiseModel, PyTrackerType,
};
use simulation::{flock_base, run_simulation, run_simulation_streaming};
use types::{PyBoidData, PyRunOptions, PyWindowSize};

/// Python module for boids simulation.
#[pymodule]
fn _pyboidr(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Register types
    m.add_class::<PyRunOptions>()?;
    m.add_class::<PyWindowSize>()?;
    m.add_class::<PyBoidData>()?;

    // Register enums
    m.add_class::<PyBoundary>()?;
    m.add_class::<PyBoundaryRepulsive>()?;
    m.add_class::<PyBoundaryRepulsiveCircle>()?;
    m.add_class::<PyDistance>()?;
    m.add_class::<PyNoiseModel>()?;
    m.add_class::<PyInitiationStrategy>()?;
    m.add_class::<PyTrackerType>()?;
    m.add_class::<PyNeighbourSampling>()?;

    // Register functions
    m.add_function(wrap_pyfunction!(run_simulation, m)?)?;
    m.add_function(wrap_pyfunction!(flock_base, m)?)?;
    m.add_function(wrap_pyfunction!(run_simulation_streaming, m)?)?;

    Ok(())
}
