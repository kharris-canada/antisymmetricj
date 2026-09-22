"""Crystallite orientation generation and text-file I/O."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class OrientationSet:
    """Immutable alpha/beta/weight arrays for powder averaging."""

    alpha: NDArray[np.float64]
    beta: NDArray[np.float64]
    weights: NDArray[np.float64]

    def __post_init__(self) -> None:
        alpha = np.asarray(self.alpha, dtype=float)
        beta = np.asarray(self.beta, dtype=float)
        weights = np.asarray(self.weights, dtype=float)

        if alpha.ndim != 1 or beta.ndim != 1 or weights.ndim != 1:
            raise ValueError("alpha, beta, and weights must be one-dimensional")
        if not (len(alpha) == len(beta) == len(weights)):
            raise ValueError("alpha, beta, and weights must have the same length")
        if len(alpha) == 0:
            raise ValueError("orientation sets cannot be empty")
        if not np.isfinite(alpha).all():
            raise ValueError("alpha values must be finite")
        if not np.isfinite(beta).all():
            raise ValueError("beta values must be finite")
        if not np.isfinite(weights).all():
            raise ValueError("weights must be finite")
        if np.any(beta < 0.0) or np.any(beta > np.pi):
            raise ValueError("beta values must lie in [0, pi] radians")
        if np.any(weights < 0.0):
            raise ValueError("weights must be non-negative")
        if float(np.sum(weights)) <= 0.0:
            raise ValueError("at least one weight must be positive")

        alpha.setflags(write=False)
        beta.setflags(write=False)
        weights.setflags(write=False)
        object.__setattr__(self, "alpha", alpha)
        object.__setattr__(self, "beta", beta)
        object.__setattr__(self, "weights", weights)

    def __len__(self) -> int:
        return len(self.alpha)


def generate_zcw(min_orientations: int = 832_039) -> OrientationSet:
    """Generate the legacy full-sphere ZCW orientation set.

    The denominator is the smallest Fibonacci number ``F_m`` for which
    ``F_m - 1`` meets the requested size. Orientations use indices
    ``1..F_m-1``. This reproduces the paper's default count because
    ``832040 - 1 == 832039``.
    """

    if min_orientations < 1:
        raise ValueError("min_orientations must be positive")

    previous, denominator = _zcw_denominator(min_orientations)
    indices = np.arange(1, denominator, dtype=np.float64)
    alpha = np.mod(indices * previous, denominator) * (2.0 * np.pi / denominator)
    beta = indices * (np.pi / denominator)
    weights = np.sin(beta)
    return OrientationSet(alpha=alpha, beta=beta, weights=weights)


def load_orientations(
    path: str | Path,
    *,
    format: str = "auto",
    angle_unit: str = "radians",
) -> OrientationSet:
    """Load orientations from plain or SIMPSON-style whitespace text."""

    file_path = Path(path)
    if format not in {"auto", "plain", "simpson"}:
        raise ValueError("format must be 'auto', 'plain', or 'simpson'")
    unit = _normalize_angle_unit(angle_unit)

    rows: list[list[float]] = []
    metadata_unit: str | None = None
    for raw_line in file_path.read_text(encoding="utf-8").splitlines():
        line = raw_line.split("#", maxsplit=1)[0].strip()
        if not line:
            lower = raw_line.lower()
            if "degree" in lower or "degrees" in lower:
                metadata_unit = "degrees"
            elif "radian" in lower or "radians" in lower:
                metadata_unit = "radians"
            continue
        rows.append([float(part) for part in line.split()])

    if not rows:
        raise ValueError(f"no orientation rows found in {file_path}")

    detected = _detect_format(rows, format)
    if detected == "simpson":
        expected = int(rows[0][0])
        data_rows = rows[1:]
        if expected != len(data_rows):
            raise ValueError(
                f"SIMPSON row count says {expected}, but {len(data_rows)} rows follow"
            )
        if angle_unit == "radians":
            unit = metadata_unit or "degrees"
    else:
        data_rows = rows
        if metadata_unit is not None and angle_unit == "radians":
            unit = metadata_unit

    data = np.asarray(data_rows, dtype=float)
    if data.ndim != 2 or data.shape[1] not in {2, 3, 4}:
        raise ValueError("orientation rows must be alpha beta [weight]")

    alpha = data[:, 0]
    beta = data[:, 1]
    if data.shape[1] == 2:
        weights = np.ones(len(data), dtype=float)
    elif data.shape[1] == 3:
        weights = data[:, 2]
    else:
        weights = data[:, 3]

    if unit == "degrees":
        alpha = np.deg2rad(alpha)
        beta = np.deg2rad(beta)

    return OrientationSet(alpha=alpha, beta=beta, weights=weights)


def write_orientations(orientation_set: OrientationSet, path: str | Path) -> None:
    """Write a simple radian-valued ``alpha beta weight`` orientation file."""

    file_path = Path(path)
    file_path.parent.mkdir(parents=True, exist_ok=True)
    data = np.column_stack(
        (orientation_set.alpha, orientation_set.beta, orientation_set.weights)
    )
    header = "alpha beta weight\nangle_unit radians"
    np.savetxt(file_path, data, fmt="%.17g", header=header)


def _zcw_denominator(min_orientations: int) -> tuple[int, int]:
    previous, current = 1, 1
    while current - 1 < min_orientations:
        previous, current = current, previous + current
    return previous, current


def _detect_format(rows: list[list[float]], requested_format: str) -> str:
    if requested_format != "auto":
        return requested_format
    if len(rows[0]) == 1 and float(rows[0][0]).is_integer():
        row_count = int(rows[0][0])
        if row_count == len(rows) - 1:
            return "simpson"
    return "plain"


def _normalize_angle_unit(angle_unit: str) -> str:
    normalized = angle_unit.lower()
    if normalized in {"radian", "radians", "rad"}:
        return "radians"
    if normalized in {"degree", "degrees", "deg"}:
        return "degrees"
    raise ValueError("angle_unit must be radians or degrees")
