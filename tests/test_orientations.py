from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from antisymmetricj.orientations import (
    OrientationSet,
    generate_zcw,
    load_orientations,
    write_orientations,
)


def test_orientation_set_is_immutable() -> None:
    orientations = OrientationSet(
        alpha=np.array([0.0]),
        beta=np.array([0.0]),
        weights=np.array([1.0]),
    )
    with pytest.raises(ValueError):
        orientations.alpha[0] = 1.0


def test_zcw_legal_counts_and_ranges() -> None:
    expected_counts = [20, 33, 54, 28_656]
    for count in expected_counts:
        orientations = generate_zcw(count)
        assert len(orientations) == count
        assert np.all(orientations.alpha >= 0.0)
        assert np.all(orientations.alpha < 2.0 * np.pi)
        assert np.all(orientations.beta > 0.0)
        assert np.all(orientations.beta < np.pi)
        assert np.all(orientations.weights > 0.0)


def test_zcw_default_paper_count() -> None:
    assert len(generate_zcw()) == 832_039


def test_zcw_is_deterministic() -> None:
    first = generate_zcw(54)
    second = generate_zcw(54)
    np.testing.assert_array_equal(first.alpha, second.alpha)
    np.testing.assert_array_equal(first.beta, second.beta)
    np.testing.assert_array_equal(first.weights, second.weights)


def test_plain_two_column_file_defaults_to_equal_radian_weights(tmp_path: Path) -> None:
    path = tmp_path / "plain.txt"
    path.write_text("# comment\n0 0\n1.5707963267948966 1.5707963267948966\n")
    orientations = load_orientations(path)
    np.testing.assert_allclose(orientations.weights, [1.0, 1.0])
    np.testing.assert_allclose(orientations.beta, [0.0, np.pi / 2.0])


def test_plain_three_column_file_uses_weights(tmp_path: Path) -> None:
    path = tmp_path / "weighted.txt"
    path.write_text("0 0 2\n1 1 3\n")
    orientations = load_orientations(path)
    np.testing.assert_allclose(orientations.weights, [2.0, 3.0])


def test_explicit_degree_conversion(tmp_path: Path) -> None:
    path = tmp_path / "degrees.txt"
    path.write_text("0 0 1\n90 90 1\n")
    orientations = load_orientations(path, angle_unit="degrees")
    np.testing.assert_allclose(orientations.alpha, [0.0, np.pi / 2.0])
    np.testing.assert_allclose(orientations.beta, [0.0, np.pi / 2.0])


def test_simpson_style_file_defaults_to_degrees(tmp_path: Path) -> None:
    path = tmp_path / "simpson.txt"
    path.write_text("2\n0 0 0 1\n90 90 0 2\n")
    orientations = load_orientations(path, format="simpson")
    np.testing.assert_allclose(orientations.alpha, [0.0, np.pi / 2.0])
    np.testing.assert_allclose(orientations.beta, [0.0, np.pi / 2.0])
    np.testing.assert_allclose(orientations.weights, [1.0, 2.0])


def test_rejects_invalid_orientations(tmp_path: Path) -> None:
    path = tmp_path / "bad_beta.txt"
    path.write_text("0 4 1\n")
    with pytest.raises(ValueError, match=r"\[0, pi\]"):
        load_orientations(path)

    negative_weight = tmp_path / "negative_weight.txt"
    negative_weight.write_text("0 0 -1\n")
    with pytest.raises(ValueError, match="non-negative"):
        load_orientations(negative_weight)


def test_write_then_load_roundtrip(tmp_path: Path) -> None:
    original = generate_zcw(20)
    path = tmp_path / "orientations.txt"
    write_orientations(original, path)
    loaded = load_orientations(path)
    np.testing.assert_allclose(loaded.alpha, original.alpha)
    np.testing.assert_allclose(loaded.beta, original.beta)
    np.testing.assert_allclose(loaded.weights, original.weights)
