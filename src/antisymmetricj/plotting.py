"""Plot and text-output helpers for f_zy distribution examples."""

from __future__ import annotations

import os
import tempfile
from collections.abc import Sequence
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(tempfile.gettempdir()) / "antisymmetricj-matplotlib"),
)
os.environ.setdefault(
    "XDG_CACHE_HOME",
    str(Path(tempfile.gettempdir()) / "antisymmetricj-cache"),
)

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np

from antisymmetricj.distribution import Distribution
from antisymmetricj.spectrum import Spectrum

FZY_PLOT_XLIM = (-1.5, 1.5)


def write_fzy_distribution_table(
    path: str | Path,
    distribution: Distribution,
    *,
    metadata: Sequence[str],
    fraction_label: str = "fraction",
) -> None:
    """Write one f_zy distribution table with metadata comments."""

    file_path = Path(path)
    file_path.parent.mkdir(parents=True, exist_ok=True)
    header = "\n".join([*metadata, f"bin_left bin_center bin_right {fraction_label}"])
    data = np.column_stack(
        (
            distribution.left,
            distribution.center,
            distribution.right,
            distribution.fraction,
        )
    )
    np.savetxt(file_path, data, fmt="%.17g", header=header)


def write_fzy_overlay_distribution_table(
    path: str | Path,
    distributions: Sequence[Distribution],
    *,
    labels: Sequence[str],
    metadata: Sequence[str],
) -> None:
    """Write a shared-bin f_zy table with one fraction column per distribution."""

    if len(distributions) == 0:
        raise ValueError("at least one distribution is required")
    reference = distributions[0]
    for distribution in distributions[1:]:
        if not np.allclose(distribution.left, reference.left):
            raise ValueError("distributions must share bin edges")
        if not np.allclose(distribution.right, reference.right):
            raise ValueError("distributions must share bin edges")

    column_labels = " ".join(_column_safe(label) for label in labels)
    header = "\n".join([*metadata, f"bin_left bin_center bin_right {column_labels}"])
    data = np.column_stack(
        (
            reference.left,
            reference.center,
            reference.right,
            *(distribution.fraction for distribution in distributions),
        )
    )
    file_path = Path(path)
    file_path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(file_path, data, fmt="%.17g", header=header)


def write_spectrum_table(
    path: str | Path,
    spectrum: Spectrum,
    *,
    metadata: Sequence[str],
) -> None:
    """Write one NMR spectrum table with metadata comments."""

    file_path = Path(path)
    file_path.parent.mkdir(parents=True, exist_ok=True)
    header = "\n".join([*metadata, "bin_left bin_center bin_right intensity"])
    data = np.column_stack(
        (
            spectrum.left,
            spectrum.center,
            spectrum.right,
            spectrum.intensity,
        )
    )
    np.savetxt(file_path, data, fmt="%.17g", header=header)


def plot_fzy_overlay(
    path: str | Path,
    distributions: Sequence[Distribution],
    *,
    labels: Sequence[str],
    xlabel: str,
    ylabel: str = "fraction per bin",
    title: str | None = None,
    colors: Sequence[str] = ("#0072B2", "#D55E00", "#009E73"),
    linestyles: Sequence[str] = ("-", "--", "-."),
    xlim: tuple[float, float] | None = FZY_PLOT_XLIM,
) -> None:
    """Write an f_zy overlay line plot as a 300-dpi JPEG."""

    file_path = Path(path)
    file_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    for index, distribution in enumerate(distributions):
        plot_x, plot_y = _padded_fzy_plot_points(distribution, xlim)
        ax.plot(
            plot_x,
            plot_y,
            color=colors[index % len(colors)],
            linestyle=linestyles[index % len(linestyles)],
            linewidth=1.6,
            label=labels[index],
        )

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if xlim is not None:
        ax.set_xlim(xlim)
    if title is not None:
        ax.set_title(title)
    ax.legend(frameon=False)
    ax.grid(True, alpha=0.25, linewidth=0.6)
    fig.tight_layout()
    fig.savefig(file_path, dpi=300, facecolor="white")
    plt.close(fig)


def plot_spectrum(
    path: str | Path,
    spectrum: Spectrum,
    *,
    xlabel: str = "frequency / Hz",
    ylabel: str = "raw relative intensity",
    title: str | None = None,
    color: str = "#0072B2",
    linewidth: float = 0.8,
) -> None:
    """Write a Gaussian-broadened NMR spectrum plot as a 300-dpi JPEG."""

    file_path = Path(path)
    file_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    ax.plot(
        spectrum.center,
        spectrum.intensity,
        color=color,
        linestyle="-",
        linewidth=linewidth,
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim((float(spectrum.left[0]), float(spectrum.right[-1])))
    if title is not None:
        ax.set_title(title)
    ax.grid(True, alpha=0.25, linewidth=0.6)
    fig.tight_layout()
    fig.savefig(file_path, dpi=300, facecolor="white")
    plt.close(fig)


def _padded_fzy_plot_points(
    distribution: Distribution,
    xlim: tuple[float, float] | None,
) -> tuple[np.ndarray, np.ndarray]:
    if xlim is None:
        return distribution.center, distribution.fraction

    left_limit, right_limit = xlim
    x_values = [left_limit]
    y_values = [0.0]

    if distribution.left[0] > left_limit:
        x_values.append(float(distribution.left[0]))
        y_values.append(0.0)

    x_values.extend(distribution.center)
    y_values.extend(distribution.fraction)

    if distribution.right[-1] < right_limit:
        x_values.append(float(distribution.right[-1]))
        y_values.append(0.0)

    x_values.append(right_limit)
    y_values.append(0.0)

    return np.asarray(x_values, dtype=float), np.asarray(y_values, dtype=float)


def _column_safe(label: str) -> str:
    safe = []
    previous_was_separator = False
    for character in label.lower():
        if character.isalnum():
            safe.append(character)
            previous_was_separator = False
        elif not previous_was_separator:
            safe.append("_")
            previous_was_separator = True
    return "".join(safe).strip("_")
