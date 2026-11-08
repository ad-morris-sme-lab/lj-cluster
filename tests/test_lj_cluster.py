"""Tests for the Lennard-Jones cluster model."""
import numpy as np
import pytest
from lj_cluster import LennardJonesParameters, lennard_jones_energy, lennard_jones_gradient
from lj_cluster.geometry import validate_positions

EPS = 2.0e-21
SIGMA = 3.0e-10
PARAMS = LennardJonesParameters(EPS, SIGMA)

def test_parameter_validation() -> None:
    """Parameters accept positive SI values and reject zero."""
    assert LennardJonesParameters(EPS, SIGMA).sigma_meter == SIGMA
    with pytest.raises(ValueError):
        LennardJonesParameters(0.0, SIGMA)

def test_geometry_validation() -> None:
    """Geometry requires two or more finite xyz rows."""
    assert validate_positions([[0, 0, 0], [1, 0, 0]]).shape == (2, 3)
    with pytest.raises(ValueError):
        validate_positions([[0, 0]])
    with pytest.raises(ValueError):
        validate_positions([[0, 0, 0], [float("nan"), 0, 0]])

def test_pair_energy_at_minimum() -> None:
    """A pair at 2^(1/6) sigma has energy -epsilon."""
    r = 2 ** (1 / 6) * SIGMA
    value = lennard_jones_energy([[0, 0, 0], [r, 0, 0]], PARAMS)
    assert value == pytest.approx(-EPS, rel=1e-12)

def test_pair_energy_at_sigma() -> None:
    """Potential crosses zero at sigma."""
    assert lennard_jones_energy([[0, 0, 0], [SIGMA, 0, 0]], PARAMS) == pytest.approx(0.0, abs=1e-34)

def test_overlapping_atoms_rejected() -> None:
    """Coincident particles raise rather than divide by zero."""
    with pytest.raises(ValueError, match="identical"):
        lennard_jones_energy([[0, 0, 0], [0, 0, 0]], PARAMS)

def test_gradient_finite_difference() -> None:
    """Analytic gradient agrees with a central finite difference."""
    coords = np.array([[0., 0., 0.], [1.3 * SIGMA, 0.2 * SIGMA, 0.]])
    analytic = lennard_jones_gradient(coords, PARAMS)[1, 0]
    step = 1e-15
    plus, minus = coords.copy(), coords.copy()
    plus[1, 0] += step
    minus[1, 0] -= step
    numeric = (lennard_jones_energy(plus, PARAMS) - lennard_jones_energy(minus, PARAMS)) / (2 * step)
    assert analytic == pytest.approx(numeric, rel=1e-7)

def test_gradient_conserves_total_force() -> None:
    """Translation invariance makes the summed gradient zero."""
    grad = lennard_jones_gradient([[0, 0, 0], [1.2 * SIGMA, 0, 0]], PARAMS)
    assert np.allclose(np.sum(grad, axis=0), 0.0, atol=1e-30)


def test_extreme_repulsive_energy_reports_overflow() -> None:
    """Unrepresentable short-range repulsion raises a clear error."""
    with pytest.raises(ValueError, match="float64 range"):
        lennard_jones_energy([[0, 0, 0], [1e-100, 0, 0]], PARAMS)

def test_extreme_gradient_reports_overflow() -> None:
    """Unrepresentable short-range gradient raises a clear error."""
    with pytest.raises(ValueError, match="float64 range"):
        lennard_jones_gradient([[0, 0, 0], [1e-100, 0, 0]], PARAMS)
