# Lennard-Jones Cluster Energy

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue)](https://www.python.org/) [![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

A compact Python package for evaluating pairwise Lennard-Jones energies and Cartesian gradients for finite clusters. Inputs and outputs use SI units.

## Why this problem

The Lennard-Jones cluster calculator provides a concrete molecular energy problem with transparent mathematics and straightforward numerical checks. It stays small enough to audit while supporting an analytic gradient and explicit edge-case handling.

## How to run this locally

From a fresh Conda environment, run these commands at the repository root:

```bash
conda create -n lj-cluster -c conda-forge python=3.12 numpy=2.5.4 pytest=9.2.0 pip
conda activate lj-cluster
python -m pip install -e .
python examples/argon_trimer.py
python -m pytest -q
```

Example output:

```text
Argon trimer energy: -4.020793e-21 J
Gradient norm: 3.481886e-11 J/m
```

## Usage

Coordinates are an N by 3 array in metres. Pair parameters are in joules and metres.

```python
import numpy as np
from lj_cluster import LennardJonesParameters, lennard_jones_energy

parameters = LennardJonesParameters(epsilon_joule=1.654e-21, sigma_meter=3.405e-10)
positions = np.array([[0., 0., 0.], [4.0e-10, 0., 0.]])
print(lennard_jones_energy(positions, parameters), "J")
```

## Cite this work

### Abstract

This project provides a small, reproducible Python implementation of the pairwise Lennard-Jones potential for finite molecular clusters. Given Cartesian coordinates and homogeneous pair parameters, it computes the total potential energy and analytic Cartesian gradient in SI units. Input validation rejects malformed geometries and coincident particles, while floating-point overflow is reported explicitly. The code is organized as an installable package with a runnable argon trimer example and unit tests covering known potential values, gradient consistency, translation invariance, and numerical edge cases. The implementation is intended for teaching and small-cluster calculations. It assumes classical, isotropic pair interactions and does not model chemical bonding, electrostatics, many-body interactions, periodic boundaries, or thermal sampling.

### BibTeX

```bibtex
@software{lj_cluster_2026,
  author = {{Lennard-Jones Cluster Energy contributors}},
  title = {Lennard-Jones Cluster Energy},
  version = {1.0.0},
  year = {2026},
  url = {https://github.com/ad-morris-sme-lab/lj-cluster}
}
```

The potential is due to J. E. Lennard-Jones, "On the determination of molecular fields. II. From the equation of state of a gas," *Proceedings of the Royal Society A* **106**, 463-477 (1924), https://doi.org/10.1098/rspa.1924.0082.

## License

MIT, see [LICENSE](LICENSE).
