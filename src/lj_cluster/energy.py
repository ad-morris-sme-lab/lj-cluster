"""Lennard-Jones potential energy calculation."""
import numpy as np
from numpy.typing import NDArray
from .geometry import validate_positions
from .parameters import LennardJonesParameters

def lennard_jones_energy(positions_meter: NDArray[np.float64], parameters: LennardJonesParameters) -> float:
    """Compute total pairwise Lennard-Jones energy in joules.

    Parameters
    ----------
    positions_meter : array_like
        N by 3 Cartesian coordinates in metres.
    parameters : LennardJonesParameters
        Homogeneous pair parameters.

    Returns
    -------
    float
        Total potential energy in joules.

    Raises
    ------
    ValueError
        If coordinates overlap or the requested energy exceeds float64 range.
    """
    coords = validate_positions(positions_meter)
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise", under="ignore"):
            delta = coords[:, None, :] - coords[None, :, :]
            distances = np.hypot(np.hypot(delta[..., 0], delta[..., 1]), delta[..., 2])
            i, j = np.triu_indices(len(coords), k=1)
            pair_distances = distances[i, j]
            if np.any(pair_distances == 0.0):
                raise ValueError("distinct atoms must not occupy identical positions (overlap at float64 precision)")
            ratio6 = (parameters.sigma_meter / pair_distances) ** 6
            energy = 4.0 * parameters.epsilon_joule * np.sum(ratio6 * (ratio6 - 1.0))
    except FloatingPointError as exc:
        raise ValueError("energy is outside the representable float64 range") from exc
    if not np.isfinite(energy):
        raise ValueError("energy is outside the representable float64 range")
    return float(energy)
