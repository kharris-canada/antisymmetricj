"""Weighted distributions for :math:`f_{zy}` and :math:`f_{zy}^2`."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray


@dataclass(frozen=True)
class Distribution:
    """Histogram bins and fraction-per-bin populations."""

    left: NDArray[np.float64]
    center: NDArray[np.float64]
    right: NDArray[np.float64]
    fraction: NDArray[np.float64]


def weighted_distribution(
    values: ArrayLike,
    weights: ArrayLike,
    *,
    squared: bool = False,
    bins: int = 200,
) -> Distribution:
    """Return weighted fraction-per-bin populations.

    The domain is fixed to ``[-1, 1]`` for ``f_zy`` and ``[0, 1]`` for
    ``f_zy**2`` so multiple parameter cases can be compared directly.
    """

    if bins < 1:
        raise ValueError("bins must be positive")

    values_array = np.asarray(values, dtype=float)
    weights_array = np.asarray(weights, dtype=float)
    if squared:
        values_array = values_array**2
        domain = (0.0, 1.0)
    else:
        domain = (-1.0, 1.0)

    if values_array.shape != weights_array.shape:
        raise ValueError("values and weights must have the same shape")
    if not np.isfinite(values_array).all():
        raise ValueError("values must be finite")
    if not np.isfinite(weights_array).all():
        raise ValueError("weights must be finite")
    if np.any(weights_array < 0.0):
        raise ValueError("weights must be non-negative")

    total_weight = float(np.sum(weights_array))
    if total_weight <= 0.0:
        raise ValueError("total weight must be positive")

    counts, edges = np.histogram(
        values_array,
        bins=bins,
        range=domain,
        weights=weights_array,
    )
    fraction = counts / total_weight
    left = edges[:-1]
    right = edges[1:]
    center = 0.5 * (left + right)
    return Distribution(left=left, center=center, right=right, fraction=fraction)
