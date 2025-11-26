"""Tests for enum types."""

from pyboidr import (
    Boundary,
    BoundaryRepulsive,
    BoundaryRepulsiveCircle,
    Distance,
    InitiationStrategy,
    NeighbourSampling,
    NoiseModel,
    TrackerType,
)


class TestBoundary:
    """Tests for Boundary enum."""

    def test_values_exist(self) -> None:
        """Verify all boundary variants exist."""
        assert Boundary.Toroidal is not None
        assert Boundary.Absorbing is not None
        assert Boundary.Reflective is not None

    def test_repr(self) -> None:
        """Verify repr is informative."""
        assert "Toroidal" in repr(Boundary.Toroidal)

    def test_equality(self) -> None:
        """Verify enum equality."""
        assert Boundary.Toroidal == Boundary.Toroidal
        assert Boundary.Toroidal != Boundary.Absorbing


class TestBoundaryRepulsive:
    """Tests for BoundaryRepulsive type."""

    def test_constructor(self) -> None:
        """Verify BoundaryRepulsive constructor."""
        repulsive = BoundaryRepulsive(distance=50.0, force=0.1)
        assert repulsive.distance == 50.0
        assert abs(repulsive.force - 0.1) < 1e-6  # f32 precision

    def test_repr(self) -> None:
        """Verify repr is informative."""
        repulsive = BoundaryRepulsive(distance=50.0, force=0.1)
        repr_str = repr(repulsive)
        assert "50" in repr_str
        assert "0.1" in repr_str


class TestBoundaryRepulsiveCircle:
    """Tests for BoundaryRepulsiveCircle type."""

    def test_constructor(self) -> None:
        """Verify BoundaryRepulsiveCircle constructor."""
        circle = BoundaryRepulsiveCircle(radius=100.0)
        assert circle.radius == 100.0


class TestDistance:
    """Tests for Distance enum."""

    def test_values_exist(self) -> None:
        """Verify all distance variants exist."""
        assert Distance.EucToroidal is not None
        assert Distance.EucEnclosed is not None


class TestNoiseModel:
    """Tests for NoiseModel enum."""

    def test_values_exist(self) -> None:
        """Verify all noise model variants exist."""
        assert NoiseModel.Vicsek is not None
        assert NoiseModel.Reynolds is not None


class TestInitiationStrategy:
    """Tests for InitiationStrategy enum."""

    def test_values_exist(self) -> None:
        """Verify all initiation strategy variants exist."""
        assert InitiationStrategy.CircleCenterOut is not None
        assert InitiationStrategy.CircleCircumferenceIn is not None
        assert InitiationStrategy.TwoWalls is not None
        assert InitiationStrategy.RectangleIn is not None
        assert InitiationStrategy.RandomIn is not None
        assert InitiationStrategy.RandomRandom is not None


class TestTrackerType:
    """Tests for TrackerType enum."""

    def test_values_exist(self) -> None:
        """Verify all tracker type variants exist."""
        assert TrackerType.SpatHash is not None
        assert TrackerType.Naive is not None


class TestNeighbourSampling:
    """Tests for NeighbourSampling enum."""

    def test_values_exist(self) -> None:
        """Verify all neighbour sampling variants exist."""
        assert NeighbourSampling.Biased is not None
        assert NeighbourSampling.Uniform is not None
