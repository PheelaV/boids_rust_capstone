"""Pytest fixtures for pyboidr tests."""

import pytest

from pyboidr import RunOptions


@pytest.fixture
def default_options() -> RunOptions:
    """Create RunOptions with default values."""
    return RunOptions()


@pytest.fixture
def small_simulation_options() -> RunOptions:
    """Create RunOptions for small, fast simulations."""
    return RunOptions(
        init_boids=16,
        width=200,
        height=200,
        sample_rate=10,
        rng_seed=42,
    )


@pytest.fixture
def deterministic_options() -> RunOptions:
    """Create deterministic RunOptions for reproducible tests."""
    return RunOptions(
        init_boids=32,
        rng_seed=12345,
        sample_rate=1,
    )
