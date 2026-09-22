from __future__ import annotations

import numpy as np
import pytest

from antisymmetricj.distribution import weighted_distribution
from antisymmetricj.geometry import f_zy
from antisymmetricj.orientations import generate_zcw


def test_weighted_distribution_fractions_sum_to_one() -> None:
    values = np.array([-0.5, 0.0, 0.5])
    weights = np.array([1.0, 2.0, 1.0])
    distribution = weighted_distribution(values, weights, bins=10)
    assert float(np.sum(distribution.fraction)) == pytest.approx(1.0)


def test_squared_distribution_uses_zero_to_one_domain() -> None:
    distribution = weighted_distribution(
        np.array([-1.0, 0.0, 1.0]),
        np.ones(3),
        squared=True,
        bins=5,
    )
    assert distribution.left[0] == pytest.approx(0.0)
    assert distribution.right[-1] == pytest.approx(1.0)
    assert float(np.sum(distribution.fraction)) == pytest.approx(1.0)


def test_rejects_bad_distribution_inputs() -> None:
    with pytest.raises(ValueError, match="same shape"):
        weighted_distribution([0.0], [1.0, 2.0])
    with pytest.raises(ValueError, match="non-negative"):
        weighted_distribution([0.0], [-1.0])


def test_zcw_distribution_is_close_to_flat_for_full_support_case() -> None:
    orientations = generate_zcw(28_656)
    values = f_zy(orientations.alpha, orientations.beta, c_yx=1.0, c_zx=1.0)
    distribution = weighted_distribution(values, orientations.weights, bins=50)

    occupied = distribution.fraction[distribution.fraction > 0.0]
    assert len(occupied) == 50
    assert float(np.max(occupied) / np.min(occupied)) < 1.25
