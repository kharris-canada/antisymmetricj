from __future__ import annotations

import numpy as np
import pytest

from antisymmetricj.geometry import SQRT_3, f_zy, simplified_f_zy


def test_mathematica_regression_value() -> None:
    assert f_zy(np.pi / 4.0, np.pi / 4.0, c_yx=1.0, c_zx=1.0) == pytest.approx(
        0.1691019787257628,
        abs=1e-15,
    )


def test_scalar_and_array_broadcasting() -> None:
    beta = np.array([0.0, np.pi / 2.0, np.pi])
    values = f_zy(0.0, beta, c_yx=1.0, c_zx=1.0)
    expected = simplified_f_zy(0.0, beta, c_yx=1.0, c_zx=1.0)
    assert isinstance(values, np.ndarray)
    np.testing.assert_allclose(values, expected)


def test_geometric_calculation_matches_simplified_oracle() -> None:
    rng = np.random.default_rng(20260801)
    alpha = rng.uniform(0.0, 2.0 * np.pi, size=1000)
    beta = rng.uniform(0.0, np.pi, size=1000)
    for c_yx, c_zx in [(1.0, 1.0), (0.0, 0.0), (0.0, 0.8), (-0.4, 1.7)]:
        np.testing.assert_allclose(
            f_zy(alpha, beta, c_yx=c_yx, c_zx=c_zx),
            simplified_f_zy(alpha, beta, c_yx=c_yx, c_zx=c_zx),
            atol=2e-15,
        )


def test_pole_orientations_are_finite() -> None:
    values = f_zy(
        np.array([0.0, np.pi / 3.0, np.pi]),
        np.array([0.0, np.pi, np.pi]),
    )
    assert np.isfinite(values).all()


def test_support_bound_for_figure2a_cases() -> None:
    alpha = np.linspace(0.0, 2.0 * np.pi, 301)
    beta = np.linspace(0.0, np.pi, 201)
    alpha_grid, beta_grid = np.meshgrid(alpha, beta)

    for c_yx, c_zx in [(0.0, 0.0), (1.0, 1.0), (0.0, 0.8)]:
        values = f_zy(alpha_grid, beta_grid, c_yx=c_yx, c_zx=c_zx)
        bound = np.sqrt(c_yx**2 + c_zx**2 + 1.0) / SQRT_3
        assert float(np.max(np.abs(values))) <= bound + 1e-12


def test_rejects_non_finite_angles() -> None:
    with pytest.raises(ValueError, match="finite"):
        f_zy(np.nan, 0.0)
