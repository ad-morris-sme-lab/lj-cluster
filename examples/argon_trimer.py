"""Evaluate a simple argon trimer geometry."""
import numpy as np
from lj_cluster import LennardJonesParameters, lennard_jones_energy, lennard_jones_gradient


def main() -> None:
    """Print energy and gradient norm for an equilateral trimer."""
    epsilon = 1.654e-21  # J, representative Ar pair well depth
    sigma = 3.405e-10   # m, representative Ar zero-crossing distance
    params = LennardJonesParameters(epsilon, sigma)
    side = 1.10 * 2 ** (1 / 6) * sigma
    positions = np.array([[0.0, 0.0, 0.0], [side, 0.0, 0.0], [side / 2, side * np.sqrt(3) / 2, 0.0]])
    energy = lennard_jones_energy(positions, params)
    grad_norm = float(np.linalg.norm(lennard_jones_gradient(positions, params)))
    print(f"Argon trimer energy: {energy:.6e} J")
    print(f"Gradient norm: {grad_norm:.6e} J/m")

if __name__ == "__main__":
    main()
