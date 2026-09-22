from __future__ import annotations

from pathlib import Path

import numpy as np
from PIL import Image

from antisymmetricj.distribution import Distribution
from antisymmetricj.plotting import (
    FZY_PLOT_XLIM,
    _padded_fzy_plot_points,
    plot_fzy_overlay,
)


def test_fzy_plot_padding_does_not_mutate_distribution(tmp_path: Path) -> None:
    distribution = Distribution(
        left=np.array([-1.0, 0.0]),
        center=np.array([-0.5, 0.5]),
        right=np.array([0.0, 1.0]),
        fraction=np.array([0.25, 0.75]),
    )
    original_center = distribution.center.copy()
    original_fraction = distribution.fraction.copy()

    plot_x, plot_y = _padded_fzy_plot_points(distribution, FZY_PLOT_XLIM)

    np.testing.assert_allclose(plot_x[[0, -1]], FZY_PLOT_XLIM)
    np.testing.assert_allclose(plot_y[[0, -1]], [0.0, 0.0])
    np.testing.assert_allclose(distribution.center, original_center)
    np.testing.assert_allclose(distribution.fraction, original_fraction)

    output = tmp_path / "fzy_plot.jpg"
    plot_fzy_overlay(
        output,
        [distribution],
        labels=["test"],
        xlabel=r"$f_{zy}$",
    )

    with Image.open(output) as image:
        assert image.size[0] > 0
        assert image.size[1] > 0
