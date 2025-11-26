"""pyboidr - Python bindings for the boids simulation library.

This package provides high-performance boid flocking simulations powered by
the boids_lib Rust library.

Quick Start:
    >>> from pyboidr import run_simulation, RunOptions
    >>> opts = RunOptions(init_boids=128, rng_seed=42)
    >>> df = run_simulation(1000, opts)
    >>> print(df.head())

The main entry points are:
    - `run_simulation`: Run a simulation and get results as a Polars DataFrame
    - `flock_base`: Run a simulation and get raw BoidData list
    - `RunOptions`: Configure simulation parameters

Enums for configuration:
    - `Boundary`: Boundary condition type (Toroidal, Absorbing, Reflective)
    - `Distance`: Distance calculation method (EucToroidal, EucEnclosed)
    - `NoiseModel`: Noise/wander model (Vicsek, Reynolds)
    - `InitiationStrategy`: Initial boid placement strategy
    - `TrackerType`: Spatial tracking algorithm (SpatHash, Naive)
    - `NeighbourSampling`: Neighbor sampling strategy (Biased, Uniform)
"""

from pyboidr._pyboidr import (
    PyBoidData as BoidData,
)
from pyboidr._pyboidr import (
    PyBoundary as Boundary,
)
from pyboidr._pyboidr import (
    PyBoundaryRepulsive as BoundaryRepulsive,
)
from pyboidr._pyboidr import (
    PyBoundaryRepulsiveCircle as BoundaryRepulsiveCircle,
)
from pyboidr._pyboidr import (
    PyDistance as Distance,
)
from pyboidr._pyboidr import (
    PyInitiationStrategy as InitiationStrategy,
)
from pyboidr._pyboidr import (
    PyNeighbourSampling as NeighbourSampling,
)
from pyboidr._pyboidr import (
    PyNoiseModel as NoiseModel,
)
from pyboidr._pyboidr import (
    PyRunOptions as RunOptions,
)
from pyboidr._pyboidr import (
    PyTrackerType as TrackerType,
)
from pyboidr._pyboidr import (
    PyWindowSize as WindowSize,
)
from pyboidr._pyboidr import (
    flock_base,
    run_simulation,
    run_simulation_streaming,
)

__version__ = "0.1.0"

__all__ = [
    "BoidData",
    "Boundary",
    "BoundaryRepulsive",
    "BoundaryRepulsiveCircle",
    "Distance",
    "InitiationStrategy",
    "NeighbourSampling",
    "NoiseModel",
    "RunOptions",
    "TrackerType",
    "WindowSize",
    "flock_base",
    "run_simulation",
    "run_simulation_streaming",
]
