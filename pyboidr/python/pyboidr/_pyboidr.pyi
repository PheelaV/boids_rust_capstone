"""Type stubs for pyboidr._pyboidr Rust extension module."""

from collections.abc import Callable
from enum import IntEnum

import polars as pl

# ============ Enums ============

class PyBoundary(IntEnum):
    """Boundary condition type for the simulation space."""

    Toroidal = 0
    """Wrap-around boundaries (agents exiting one side appear on the opposite)."""
    Absorbing = 1
    """Agents stop at boundaries."""
    Reflective = 2
    """Agents bounce off boundaries."""

class PyBoundaryRepulsive:
    """Repulsive boundary configuration."""

    distance: float
    force: float

    def __init__(self, distance: float, force: float) -> None: ...

class PyBoundaryRepulsiveCircle:
    """Circular repulsive boundary configuration."""

    radius: float

    def __init__(self, radius: float) -> None: ...

class PyDistance(IntEnum):
    """Distance calculation method."""

    EucToroidal = 0
    """Euclidean distance with toroidal wrap-around."""
    EucEnclosed = 1
    """Standard Euclidean distance in enclosed space."""

class PyNoiseModel(IntEnum):
    """Noise/wander model for agent movement."""

    Vicsek = 0
    """Vicsek-style random noise."""
    Reynolds = 1
    """Reynolds-style wander behavior."""

class PyInitiationStrategy(IntEnum):
    """Initial positioning strategy for boids."""

    CircleCenterOut = 0
    CircleCircumferenceIn = 1
    TwoWalls = 2
    RectangleIn = 3
    RandomIn = 4
    RandomRandom = 5

class PyTrackerType(IntEnum):
    """Spatial tracking algorithm type."""

    SpatHash = 0
    """Spatial hashing (O(n) average for neighbor queries)."""
    Naive = 1
    """Naive O(n^2) neighbor search."""

class PyNeighbourSampling(IntEnum):
    """Neighbor sampling strategy."""

    Biased = 0
    """Biased sampling (early exit, faster but biased)."""
    Uniform = 1
    """Uniform strided sampling across all neighboring cells."""

# ============ Types ============

class PyWindowSize:
    """Window/simulation space dimensions."""

    win_left: int
    win_right: int
    win_top: int
    win_bottom: int
    win_h: int
    win_w: int

    def __init__(self, width: int, height: int) -> None: ...

class PyBoidData:
    """Recorded boid state data at a specific time step."""

    @property
    def id(self) -> int: ...
    @property
    def x(self) -> float: ...
    @property
    def y(self) -> float: ...
    @property
    def cluster_id(self) -> int: ...
    @property
    def n_neighbours(self) -> int: ...
    @property
    def time(self) -> int: ...

class PyRunOptions:
    """Main configuration object for boids simulation."""

    # Core properties
    init_boids: int
    sample_rate: int
    rng_seed: int | None

    # Speed properties
    baseline_speed: float
    min_speed: float
    max_speed: float
    max_steering: float

    # Behavior coefficients
    separation_coefficient: float
    cohesion_coefficient: float
    alignment_coefficient: float

    # Sensory properties
    sensory_distance: float
    alignment_treshold_coefficient: float
    cohesion_treshold_coefficient: float
    separation_treshold_coefficient: float

    # Behavior toggles
    separation_on: bool
    cohesion_on: bool
    alignment_on: bool

    # Wander properties
    wander_on: bool
    wander_coefficient: float
    wander_rate: float
    wander_radius: float
    wander_distance: float

    # Field of vision
    field_of_vision_on: bool
    field_of_vision_deg: float
    alignment_fov_deg: float
    cohesion_fov_deg: float
    separation_fov_deg: float

    # Other properties
    size: float
    dbscan_flock_clustering_on: bool
    neighbours_cosidered: int
    agent_steering: bool
    stop_movement: bool

    # Complex properties
    window: PyWindowSize
    boundary: PyBoundary
    distance: PyDistance
    noise_model: PyNoiseModel
    initiation_strategy: PyInitiationStrategy
    tracker_type: PyTrackerType
    neighbour_sampling: PyNeighbourSampling

    def __init__(
        self,
        init_boids: int = 256,
        width: int = 600,
        height: int = 600,
        sample_rate: int = 1,
        rng_seed: int | None = None,
        separation_coefficient: float | None = None,
        cohesion_coefficient: float | None = None,
        alignment_coefficient: float | None = None,
        min_speed: float | None = None,
        max_speed: float | None = None,
        max_steering: float | None = None,
        sensory_distance: float | None = None,
        field_of_vision_deg: float | None = None,
        separation_on: bool | None = None,
        cohesion_on: bool | None = None,
        alignment_on: bool | None = None,
        wander_on: bool | None = None,
        boundary: PyBoundary | None = None,
        distance: PyDistance | None = None,
        noise_model: PyNoiseModel | None = None,
        initiation_strategy: PyInitiationStrategy | None = None,
        tracker_type: PyTrackerType | None = None,
        neighbour_sampling: PyNeighbourSampling | None = None,
    ) -> None:
        """Create a new RunOptions with sensible defaults.

        Args:
            init_boids: Number of boids to simulate (default: 256)
            width: Simulation window width (default: 600)
            height: Simulation window height (default: 600)
            sample_rate: How often to sample boid positions (default: 1)
            rng_seed: Random seed for reproducibility (default: None = random)
            separation_coefficient: Strength of separation behavior
            cohesion_coefficient: Strength of cohesion behavior
            alignment_coefficient: Strength of alignment behavior
            min_speed: Minimum boid speed
            max_speed: Maximum boid speed
            max_steering: Maximum steering force
            sensory_distance: Base sensing radius
            field_of_vision_deg: Field of vision in degrees
            separation_on: Enable separation behavior
            cohesion_on: Enable cohesion behavior
            alignment_on: Enable alignment behavior
            wander_on: Enable wander behavior
            boundary: Boundary type
            distance: Distance calculation method
            noise_model: Noise/wander model
            initiation_strategy: Initial boid placement
            tracker_type: Spatial tracking algorithm
            neighbour_sampling: Neighbor sampling strategy
        """
        ...

    def set_window_size(self, width: int, height: int) -> None:
        """Set window size using width and height."""
        ...

    def set_boundary_repulsive(self, distance: float, force: float) -> None:
        """Set a repulsive boundary that pushes boids away from edges."""
        ...

    def set_boundary_repulsive_circle(self, radius: float) -> None:
        """Set a circular repulsive boundary."""
        ...

# ============ Functions ============

def run_simulation(no_iter: int, options: PyRunOptions) -> pl.DataFrame:
    """Run a boids simulation and return results as a Polars DataFrame.

    Args:
        no_iter: Number of simulation iterations to run.
        options: Simulation configuration options.

    Returns:
        A Polars DataFrame with columns: id, x, y, time, cluster_id, n_neighbours.

    Example:
        >>> from pyboidr import run_simulation, RunOptions
        >>> opts = RunOptions(init_boids=64, rng_seed=42)
        >>> df = run_simulation(100, opts)
        >>> print(df.head())
    """
    ...

def flock_base(no_iter: int, options: PyRunOptions) -> list[PyBoidData]:
    """Run a boids simulation and return raw BoidData list.

    Args:
        no_iter: Number of simulation iterations to run.
        options: Simulation configuration options.

    Returns:
        List of BoidData objects containing position/state data.
    """
    ...

def run_simulation_streaming(
    no_iter: int,
    options: PyRunOptions,
    batch_size: int,
    callback: Callable[[pl.DataFrame], None],
) -> None:
    """Run a simulation with a callback for streaming results.

    Args:
        no_iter: Number of simulation iterations to run.
        options: Simulation configuration options.
        batch_size: Number of iterations between callbacks.
        callback: Python callable that receives a Polars DataFrame for each batch.

    Example:
        >>> def process_batch(df):
        ...     print(f"Got {len(df)} records")
        >>> run_simulation_streaming(1000, opts, 100, process_batch)
    """
    ...
