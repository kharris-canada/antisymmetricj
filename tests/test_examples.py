from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np


def test_reproduce_figure2_example(tmp_path: Path) -> None:
    example_path = Path("examples/reproduce_figure2.py")
    spec = importlib.util.spec_from_file_location("reproduce_figure2", example_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    output_dir = tmp_path / "output"
    module.main(output_dir=output_dir, min_orientations=28_656, bins=50)

    expected = [
        "figure2a_intensities.txt",
        "figure2a.jpg",
        "figure2b_intensities.txt",
        "figure2b.jpg",
    ]
    for name in expected:
        assert (output_dir / name).stat().st_size > 0

    figure2a = np.loadtxt(output_dir / "figure2a_intensities.txt")
    assert figure2a[0, 0] == -1.0
    assert figure2a[-1, 2] == 1.0
    np.testing.assert_allclose(np.sum(figure2a[:, 3:], axis=0), np.ones(3))

    figure2b = np.loadtxt(output_dir / "figure2b_intensities.txt")
    assert figure2b[0, 0] == 0.0
    assert figure2b[-1, 2] == 1.0
    assert figure2b[0, 3] > figure2b[-1, 3]
