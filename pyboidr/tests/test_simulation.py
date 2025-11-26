"""Tests for simulation execution and DataFrame output."""

import polars as pl

from pyboidr import RunOptions, flock_base, run_simulation


class TestRunSimulation:
    """Tests for the run_simulation function."""

    def test_returns_polars_dataframe(self, small_simulation_options: RunOptions) -> None:
        """Verify run_simulation returns a Polars DataFrame."""
        df = run_simulation(10, small_simulation_options)
        assert isinstance(df, pl.DataFrame)

    def test_has_expected_columns(self, small_simulation_options: RunOptions) -> None:
        """Verify DataFrame has all expected columns."""
        df = run_simulation(10, small_simulation_options)
        expected_columns = {"id", "x", "y", "time", "cluster_id", "n_neighbours"}
        assert set(df.columns) == expected_columns

    def test_non_empty_result(self, small_simulation_options: RunOptions) -> None:
        """Verify simulation produces non-empty results."""
        df = run_simulation(10, small_simulation_options)
        assert len(df) > 0

    def test_deterministic_with_seed(self, deterministic_options: RunOptions) -> None:
        """Verify same seed produces identical results."""
        df1 = run_simulation(50, deterministic_options)
        df2 = run_simulation(50, deterministic_options)

        assert len(df1) == len(df2)
        assert df1.equals(df2)

    def test_different_seeds_produce_different_results(self) -> None:
        """Verify different seeds produce different results."""
        opts1 = RunOptions(init_boids=16, rng_seed=1)
        opts2 = RunOptions(init_boids=16, rng_seed=2)

        df1 = run_simulation(50, opts1)
        df2 = run_simulation(50, opts2)

        # Results should differ (very unlikely to be identical)
        assert not df1.equals(df2)


class TestFlockBase:
    """Tests for the flock_base function."""

    def test_returns_list_of_boid_data(self, small_simulation_options: RunOptions) -> None:
        """Verify flock_base returns a list of BoidData."""
        data = flock_base(10, small_simulation_options)
        assert isinstance(data, list)
        assert len(data) > 0

    def test_boid_data_has_attributes(self, small_simulation_options: RunOptions) -> None:
        """Verify BoidData objects have expected attributes."""
        data = flock_base(10, small_simulation_options)
        boid = data[0]

        assert hasattr(boid, "id")
        assert hasattr(boid, "x")
        assert hasattr(boid, "y")
        assert hasattr(boid, "time")
        assert hasattr(boid, "cluster_id")
        assert hasattr(boid, "n_neighbours")


class TestSampleRate:
    """Tests for sample_rate configuration."""

    def test_sample_rate_affects_data_size(self) -> None:
        """Verify sample_rate controls how often data is collected."""
        opts_frequent = RunOptions(init_boids=16, sample_rate=1, rng_seed=42)
        opts_sparse = RunOptions(init_boids=16, sample_rate=10, rng_seed=42)

        df_frequent = run_simulation(100, opts_frequent)
        df_sparse = run_simulation(100, opts_sparse)

        # Frequent sampling should produce ~10x more records
        assert len(df_frequent) > len(df_sparse)


class TestStreamingSimulation:
    """Tests for streaming simulation."""

    def test_streaming_collects_all_data(self, small_simulation_options: RunOptions) -> None:
        """Verify streaming simulation collects all data via callbacks."""
        from pyboidr import run_simulation_streaming

        collected: list[pl.DataFrame] = []

        def collector(df: pl.DataFrame) -> None:
            collected.append(df)

        run_simulation_streaming(100, small_simulation_options, 25, collector)

        # Should have 4 batches (100 / 25)
        assert len(collected) >= 1

        # Total records should match non-streaming
        total_streaming = sum(len(df) for df in collected)
        df_regular = run_simulation(100, small_simulation_options)
        assert total_streaming == len(df_regular)
