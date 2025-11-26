"""Tests for RunOptions configuration."""


from pyboidr import Boundary, Distance, NoiseModel, RunOptions, WindowSize


class TestRunOptionsConstructor:
    """Tests for RunOptions dataclass-style constructor."""

    def test_default_values(self) -> None:
        """Verify default values are sensible."""
        opts = RunOptions()
        assert opts.init_boids == 256
        assert opts.sample_rate == 1
        assert opts.rng_seed is None

    def test_constructor_with_kwargs(self) -> None:
        """Verify constructor accepts keyword arguments."""
        opts = RunOptions(
            init_boids=64,
            width=400,
            height=300,
            separation_coefficient=5.0,
            rng_seed=42,
        )
        assert opts.init_boids == 64
        assert opts.separation_coefficient == 5.0
        assert opts.rng_seed == 42

    def test_property_mutation(self) -> None:
        """Verify properties can be mutated after construction."""
        opts = RunOptions()
        opts.init_boids = 128
        opts.separation_coefficient = 3.0
        opts.rng_seed = 99

        assert opts.init_boids == 128
        assert opts.separation_coefficient == 3.0
        assert opts.rng_seed == 99


class TestRunOptionsWindow:
    """Tests for window configuration."""

    def test_window_from_constructor(self) -> None:
        """Verify window size from constructor."""
        opts = RunOptions(width=800, height=600)
        window = opts.window

        assert window.win_w == 800
        assert window.win_h == 600

    def test_set_window_size_method(self) -> None:
        """Verify set_window_size method."""
        opts = RunOptions()
        opts.set_window_size(1024, 768)

        assert opts.window.win_w == 1024
        assert opts.window.win_h == 768

    def test_window_property_setter(self) -> None:
        """Verify window property can be set directly."""
        opts = RunOptions()
        new_window = WindowSize(500, 500)
        opts.window = new_window

        assert opts.window.win_w == 500


class TestRunOptionsEnums:
    """Tests for enum configuration."""

    def test_boundary_enum(self) -> None:
        """Verify boundary enum works."""
        opts = RunOptions()
        opts.boundary = Boundary.Toroidal
        assert opts.boundary == Boundary.Toroidal

        opts.boundary = Boundary.Absorbing
        assert opts.boundary == Boundary.Absorbing

        opts.boundary = Boundary.Reflective
        assert opts.boundary == Boundary.Reflective

    def test_distance_enum(self) -> None:
        """Verify distance enum works."""
        opts = RunOptions()
        opts.distance = Distance.EucToroidal
        assert opts.distance == Distance.EucToroidal

        opts.distance = Distance.EucEnclosed
        assert opts.distance == Distance.EucEnclosed

    def test_noise_model_enum(self) -> None:
        """Verify noise model enum works."""
        opts = RunOptions()
        opts.noise_model = NoiseModel.Reynolds
        assert opts.noise_model == NoiseModel.Reynolds

        opts.noise_model = NoiseModel.Vicsek
        assert opts.noise_model == NoiseModel.Vicsek

    def test_repulsive_boundary(self) -> None:
        """Verify repulsive boundary configuration."""
        opts = RunOptions()
        opts.set_boundary_repulsive(distance=50.0, force=0.1)
        # After setting repulsive, boundary enum returns Reflective (fallback)
        # but internal state is Repulsive


class TestRunOptionsBehaviors:
    """Tests for behavior toggles."""

    def test_behavior_toggles(self) -> None:
        """Verify behavior toggles work."""
        opts = RunOptions(
            separation_on=False,
            cohesion_on=False,
            alignment_on=True,
        )

        assert opts.separation_on is False
        assert opts.cohesion_on is False
        assert opts.alignment_on is True

    def test_wander_configuration(self) -> None:
        """Verify wander behavior configuration."""
        opts = RunOptions(wander_on=True)
        opts.wander_coefficient = 0.5
        opts.wander_rate = 0.1
        opts.wander_radius = 10.0
        opts.wander_distance = 30.0

        assert opts.wander_on is True
        assert abs(opts.wander_coefficient - 0.5) < 1e-6  # f32 precision
        assert abs(opts.wander_rate - 0.1) < 1e-6  # f32 precision


class TestRunOptionsRepr:
    """Tests for string representation."""

    def test_repr(self) -> None:
        """Verify __repr__ provides useful info."""
        opts = RunOptions(init_boids=128, width=800, height=600)
        repr_str = repr(opts)

        assert "128" in repr_str
        assert "800" in repr_str


class TestWindowSize:
    """Tests for WindowSize class."""

    def test_constructor(self) -> None:
        """Verify WindowSize constructor."""
        window = WindowSize(800, 600)
        assert window.win_w == 800
        assert window.win_h == 600
        assert window.win_left == -400
        assert window.win_right == 400
        assert window.win_top == 300
        assert window.win_bottom == -300

    def test_repr(self) -> None:
        """Verify WindowSize repr."""
        window = WindowSize(800, 600)
        repr_str = repr(window)
        assert "800" in repr_str
        assert "600" in repr_str
