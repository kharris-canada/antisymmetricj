from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from PIL import Image

from antisymmetricj.cli import build_parser, main


def test_cli_zcw_and_distribution(tmp_path: Path) -> None:
    orientations = tmp_path / "orientations_outfile.txt"
    intensities = tmp_path / "distribution_outfile.txt"
    plot = tmp_path / "plot.jpg"

    assert main(["make_zcw_set", str(orientations), "--min-orientations", "54"]) == 0
    assert orientations.exists()

    assert (
        main(
            [
                "calc_fzy",
                str(orientations),
                str(intensities),
                "--plot_fzy",
                str(plot),
                "--bins",
                "20",
                "--c-yx",
                "1",
                "--c-zx",
                "1",
            ]
        )
        == 0
    )
    assert intensities.stat().st_size > 0
    assert plot.stat().st_size > 0
    with Image.open(plot) as image:
        assert image.size[0] > 0
        assert image.size[1] > 0


def test_cli_distribution_can_generate_default_orientations(tmp_path: Path) -> None:
    intensities = tmp_path / "distribution_outfile.txt"
    plot = tmp_path / "default_orientations.jpg"

    assert (
        main(
            [
                "calc_fzy",
                str(intensities),
                "--plot_fzy",
                str(plot),
                "--bins",
                "20",
                "--c-yx",
                "1",
                "--c-zx",
                "1",
            ]
        )
        == 0
    )

    assert intensities.stat().st_size > 0
    assert plot.stat().st_size > 0
    data = np.loadtxt(intensities)
    assert data.shape == (20, 4)
    np.testing.assert_allclose(np.sum(data[:, 3]), 1.0)
    assert "orientations 28656" in intensities.read_text()


def test_cli_distribution_default_orientation_count_can_be_overridden(
    tmp_path: Path,
) -> None:
    intensities = tmp_path / "distribution_outfile.txt"

    assert (
        main(
            [
                "calc_fzy",
                str(intensities),
                "--min-orientations",
                "54",
            ]
        )
        == 0
    )

    assert "orientations 54" in intensities.read_text()


def test_cli_distribution_help_names_distribution_outfile(
    capsys: pytest.CaptureFixture[str],
) -> None:
    parser = build_parser()

    with pytest.raises(SystemExit):
        parser.parse_args(["calc_fzy", "--help"])

    captured = capsys.readouterr()
    assert "DISTRIBUTION_OUTFILE" in captured.out
    assert "ORIENTATIONS_INFILE" in captured.out

    with pytest.raises(SystemExit):
        parser.parse_args(["make_zcw_set", "--help"])

    captured = capsys.readouterr()
    assert "ORIENTATIONS_OUTFILE" in captured.out


def test_cli_rejects_plot_abbreviation(tmp_path: Path) -> None:
    parser = build_parser()

    with pytest.raises(SystemExit):
        parser.parse_args(
            [
                "calc_fzy",
                str(tmp_path / "distribution_outfile.txt"),
                "--plot",
                str(tmp_path / "plot.jpg"),
            ]
        )
