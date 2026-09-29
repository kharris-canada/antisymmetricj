"""Tools for antisymmetric J-coupling powder-distribution simulations."""

from antisymmetricj.distribution import Distribution, weighted_distribution
from antisymmetricj.geometry import f_zy
from antisymmetricj.orientations import (
    OrientationSet,
    generate_zcw,
    load_orientations,
    write_orientations,
)
from antisymmetricj.spectrum import Spectrum, Transitions, calculate_spectrum

__all__ = [
    "Distribution",
    "OrientationSet",
    "Spectrum",
    "Transitions",
    "calculate_spectrum",
    "f_zy",
    "generate_zcw",
    "load_orientations",
    "weighted_distribution",
    "write_orientations",
]

__version__ = "0.1.0"
