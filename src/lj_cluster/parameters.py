"""Validated Lennard-Jones pair parameters."""
from dataclasses import dataclass
import math

@dataclass(frozen=True)
class LennardJonesParameters:
    """Pair parameters in SI units.

    Parameters
    ----------
    epsilon_joule : float
        Well depth in joules; must be finite and positive.
    sigma_meter : float
        Zero-crossing distance in metres; must be finite and positive.
    """
    epsilon_joule: float
    sigma_meter: float

    def __post_init__(self) -> None:
        """Reject nonphysical or non-finite pair parameters."""
        if not math.isfinite(self.epsilon_joule) or self.epsilon_joule <= 0:
            raise ValueError("epsilon_joule must be finite and positive")
        if not math.isfinite(self.sigma_meter) or self.sigma_meter <= 0:
            raise ValueError("sigma_meter must be finite and positive")
