"""Cartesian gradient of the Lennard-Jones potential."""
import numpy as np
from numpy.typing import NDArray
from .geometry import validate_positions
from .parameters import LennardJonesParameters

def lennard_jones_gradient(positions_meter: NDArray[np.float64], parameters: LennardJonesParameters) -> NDArray[np.float64]:
    """Compute the Cartesian energy gradient in joules per metre.

    Parameters
    ----------
    positions_meter : array_like
        N by 3 Cartesian coordinates in metres.
    parameters : LennardJonesParameters
        Homogeneous pair parameters.

    Returns
    -------
    numpy.ndarray
        Gradient array with shape (N, 3), in J/m.

    Raises
    ------
    ValueError
        If coordinates overlap or the requested gradient exceeds float64 range.
    """
    coords = validate_positions(positions_meter)
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise", under="ignore"):
            delta = coords[:, None, :] - coords[None, :, :]
            distances = np.hypot(np.hypot(delta[..., 0], delta[..., 1]), delta[..., 2])
            i, j = np.triu_indices(len(coords), k=1)
            rij = distances[i, j]
            if np.any(rij == 0.0):
                raise ValueError("distinct atoms must not occupy identical positions (overlap at float64 precision)")
            ratio6 = (parameters.sigma_meter / rij) ** 6
            factor = 24.0 * parameters.epsilon_joule * (ratio6 - 2.0 * ratio6**2) / rij**2
            pair_grad = factor[:, None] * delta[i, j]
            grad = np.zeros_like(coords)
            np.add.at(grad, i, pair_grad)
            np.add.at(grad, j, -pair_grad)
    except FloatingPointError as exc:
        raise ValueError("gradient is outside the representable float64 range") from exc
    if not np.all(np.isfinite(grad)):
        raise ValueError("gradient is outside the representable float64 range")
    return grad
