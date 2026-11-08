"""Geometry validation helpers."""
import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]

def validate_positions(positions_meter: FloatArray) -> FloatArray:
    """Return validated Cartesian coordinates.

    Parameters
    ----------
    positions_meter : array_like
        N by 3 finite Cartesian coordinates in metres.

    Returns
    -------
    numpy.ndarray
        A float64 view or copy of shape (N, 3).

    Raises
    ------
    ValueError
        If coordinates are malformed, non-finite, or contain fewer than two atoms.
    """
    coords = np.asarray(positions_meter, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 3 or coords.shape[0] < 2:
        raise ValueError("positions_meter must have shape (N, 3), N >= 2")
    if not np.all(np.isfinite(coords)):
        raise ValueError("positions_meter must contain only finite values")
    return coords
