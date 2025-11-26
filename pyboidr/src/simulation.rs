//! High-level simulation API for running boid simulations.

use pyo3::prelude::*;
use pyo3::types::PyDict;

use boids_lib::birdwatcher::{BoidData, Birdwatcher};
use boids_lib::flock::Flock;

use crate::types::{PyBoidData, PyRunOptions};

/// Run a boids simulation and return results as a Polars DataFrame.
///
/// This is the main entry point for running simulations. It runs the simulation
/// for the specified number of iterations and returns the recorded boid positions
/// and states as a Polars DataFrame.
///
/// Args:
///     no_iter: Number of simulation iterations to run.
///     options: Simulation configuration options.
///
/// Returns:
///     A Polars DataFrame with columns: id, x, y, time, cluster_id, n_neighbours.
///
/// Example:
///     >>> from pyboidr import run_simulation, RunOptions
///     >>> opts = RunOptions(init_boids=64, rng_seed=42)
///     >>> df = run_simulation(100, opts)
///     >>> print(df.head())
#[pyfunction]
#[pyo3(signature = (no_iter, options))]
pub fn run_simulation(py: Python<'_>, no_iter: u64, options: &PyRunOptions) -> PyResult<Py<PyAny>> {
    // Run the simulation
    let data = boids_lib::flock_base(no_iter, options.to_rust());

    // Convert to Polars DataFrame
    data_to_polars_dataframe(py, &data)
}

/// Run a boids simulation and return raw BoidData list.
///
/// This is a lower-level function that returns the raw simulation data
/// without converting to a DataFrame. Useful if you need custom processing.
///
/// Args:
///     no_iter: Number of simulation iterations to run.
///     options: Simulation configuration options.
///
/// Returns:
///     List of BoidData objects containing position/state data.
#[pyfunction]
#[pyo3(signature = (no_iter, options))]
pub fn flock_base(no_iter: u64, options: &PyRunOptions) -> PyResult<Vec<PyBoidData>> {
    let data = boids_lib::flock_base(no_iter, options.to_rust());
    Ok(data.iter().map(PyBoidData::from_rust).collect())
}

/// Run a simulation with a callback for streaming results.
///
/// This function allows you to process results as they are generated,
/// which can be useful for very long simulations or real-time visualization.
///
/// Args:
///     no_iter: Number of simulation iterations to run.
///     options: Simulation configuration options.
///     batch_size: Number of iterations between callbacks.
///     callback: Python callable that receives a Polars DataFrame for each batch.
///
/// Example:
///     >>> def process_batch(df):
///     ...     print(f"Got {len(df)} records")
///     >>> run_simulation_streaming(1000, opts, 100, process_batch)
#[pyfunction]
#[pyo3(signature = (no_iter, options, batch_size, callback))]
pub fn run_simulation_streaming(
    py: Python<'_>,
    no_iter: u64,
    options: &PyRunOptions,
    batch_size: u64,
    callback: Py<PyAny>,
) -> PyResult<()> {
    let mut ro = options.to_rust();
    let mut flock = Flock::new(&ro);
    let mut bird_watcher = Birdwatcher::new(ro.sample_rate);

    for i in 0..no_iter {
        flock.update(&mut ro);
        bird_watcher.watch(&flock);

        // Call callback with batch
        if (i + 1) % batch_size == 0 || i == no_iter - 1 {
            let data = bird_watcher.pop_data();
            if !data.is_empty() {
                let df = data_to_polars_dataframe(py, &data)?;
                callback.call1(py, (df,))?;
            }
        }
    }

    Ok(())
}

/// Convert BoidData to a Polars DataFrame.
fn data_to_polars_dataframe(py: Python<'_>, data: &[BoidData]) -> PyResult<Py<PyAny>> {
    // Build column vectors
    let ids: Vec<u64> = data.iter().map(|d| d.id as u64).collect();
    let xs: Vec<f32> = data.iter().map(|d| d.x).collect();
    let ys: Vec<f32> = data.iter().map(|d| d.y).collect();
    let times: Vec<u64> = data.iter().map(|d| d.time).collect();
    let cluster_ids: Vec<u64> = data.iter().map(|d| d.cluster_id as u64).collect();
    let n_neighbours: Vec<u64> = data.iter().map(|d| d.n_neighbours as u64).collect();

    // Create dict for polars.DataFrame constructor
    let dict = PyDict::new(py);
    dict.set_item("id", ids)?;
    dict.set_item("x", xs)?;
    dict.set_item("y", ys)?;
    dict.set_item("time", times)?;
    dict.set_item("cluster_id", cluster_ids)?;
    dict.set_item("n_neighbours", n_neighbours)?;

    // Import polars and create DataFrame
    let polars = py.import("polars")?;
    let df = polars.call_method1("DataFrame", (dict,))?;

    Ok(df.into())
}
