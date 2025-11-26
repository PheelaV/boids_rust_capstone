# pyboidr

Python bindings for the boids simulation library.

## Installation

```bash
cd pyboidr
uv sync
```

## Quick Start

```python
from pyboidr import run_simulation, RunOptions

# Create configuration
opts = RunOptions(
    init_boids=128,
    width=800,
    height=600,
    rng_seed=42,
)

# Run simulation - returns Polars DataFrame
df = run_simulation(1000, opts)
print(df.head())
```

## Configuration

The `RunOptions` class provides full control over simulation parameters:

```python
opts = RunOptions(
    init_boids=256,           # Number of boids
    width=600,                # Simulation width
    height=600,               # Simulation height
    sample_rate=1,            # Data sampling frequency
    rng_seed=None,            # Random seed (None = random)

    # Behavior coefficients
    separation_coefficient=4.1,
    cohesion_coefficient=0.002,
    alignment_coefficient=0.02,

    # Speed parameters
    min_speed=0.65,
    max_speed=4.1,
    max_steering=0.7,

    # Sensory parameters
    sensory_distance=60.0,
    field_of_vision_deg=181.0,
)
```

## Enums

```python
from pyboidr import Boundary, Distance, NoiseModel

# Boundary types
opts.boundary = Boundary.Toroidal    # Wrap-around
opts.boundary = Boundary.Absorbing   # Stop at edges
opts.boundary = Boundary.Reflective  # Bounce off edges

# Distance calculation
opts.distance = Distance.EucToroidal  # With wrap-around
opts.distance = Distance.EucEnclosed  # Standard Euclidean

# Noise model
opts.noise_model = NoiseModel.Reynolds  # Reynolds wander
opts.noise_model = NoiseModel.Vicsek    # Vicsek noise
```

## Development

```bash
# Install dev dependencies
uv sync --dev

# Run tests
pytest -v

# Lint
ruff check .
ruff format .
```
