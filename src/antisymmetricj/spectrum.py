"""AB spectrum calculations including antisymmetric J-coupling effects."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray


@dataclass(frozen=True)
class Spectrum:
    """Binned NMR spectrum with raw relative intensities."""

    left: NDArray[np.float64]
    center: NDArray[np.float64]
    right: NDArray[np.float64]
    intensity: NDArray[np.float64]


@dataclass(frozen=True)
class Transitions:
    """Four Table 2 transition frequencies and relative intensities."""

    frequency: NDArray[np.float64]
    intensity: NDArray[np.float64]


def table2_transitions(
    fzy_values: ArrayLike,
    *,
    sigma_i_hz: float,
    sigma_s_hz: float,
    j_iso_hz: float,
    j_zy_anti_hz: float,
) -> Transitions:
    """Return Table 2 transition frequencies and relative intensities.

    The common center frequency term, ``nu_L [1 - 1/2 (sigma_I + sigma_S)]``,
    is intentionally omitted.  The returned arrays have shape ``(4, ...)`` in
    the Table 2 order ``v_34``, ``v_12``, ``v_24``, ``v_13``.
    """

    fzy_array = np.asarray(fzy_values, dtype=float)
    if not np.isfinite(fzy_array).all():
        raise ValueError("f_zy values must be finite")

    delta_s_hz = sigma_i_hz - sigma_s_hz
    a_hz = j_zy_anti_hz * fzy_array
    c_hz = np.sqrt(j_iso_hz**2 + a_hz**2 + delta_s_hz**2)
    if np.any(c_hz <= 0.0):
        raise ValueError("C must be positive for Table 2 intensities")

    half_c = 0.5 * c_hz
    half_j = 0.5 * j_iso_hz
    frequency = np.stack(
        (
            half_c + half_j,
            half_c - half_j,
            -half_c + half_j,
            -half_c - half_j,
        )
    )

    j_over_c = j_iso_hz / c_hz
    intensity = np.stack(
        (
            1.0 - j_over_c,
            1.0 + j_over_c,
            1.0 + j_over_c,
            1.0 - j_over_c,
        )
    )

    return Transitions(frequency=frequency, intensity=intensity)


def calculate_spectrum(
    fzy_values: ArrayLike,
    weights: ArrayLike,
    *,
    sigma_i_hz: float,
    sigma_s_hz: float,
    j_iso_hz: float,
    j_zy_anti_hz: float,
    bins: int = 4096,
    fwhm_hz: float = 1.0,
) -> Spectrum:
    """Return a Gaussian-broadened raw relative-intensity spectrum."""

    if bins < 1:
        raise ValueError("bins must be positive")
    if fwhm_hz <= 0.0:
        raise ValueError("fwhm_hz must be positive")

    fzy_array = np.asarray(fzy_values, dtype=float)
    weights_array = np.asarray(weights, dtype=float)
    if fzy_array.shape != weights_array.shape:
        raise ValueError("f_zy values and weights must have the same shape")
    if not np.isfinite(weights_array).all():
        raise ValueError("weights must be finite")
    if np.any(weights_array < 0.0):
        raise ValueError("weights must be non-negative")

    spectrum_limit = abs(j_iso_hz) + abs(j_zy_anti_hz)
    if spectrum_limit <= 0.0:
        raise ValueError("|J_iso_hz| + |J_zy_anti_hz| must be positive")

    transitions = table2_transitions(
        fzy_array,
        sigma_i_hz=sigma_i_hz,
        sigma_s_hz=sigma_s_hz,
        j_iso_hz=j_iso_hz,
        j_zy_anti_hz=j_zy_anti_hz,
    )
    weighted_intensity = transitions.intensity * weights_array

    counts, edges = np.histogram(
        transitions.frequency.ravel(),
        bins=bins,
        range=(-spectrum_limit, spectrum_limit),
        weights=weighted_intensity.ravel(),
    )
    broadened = _gaussian_broaden_counts(
        counts.astype(float),
        bin_width=edges[1] - edges[0],
        fwhm_hz=fwhm_hz,
    )

    left = edges[:-1]
    right = edges[1:]
    center = 0.5 * (left + right)
    return Spectrum(left=left, center=center, right=right, intensity=broadened)


def _gaussian_broaden_counts(
    counts: NDArray[np.float64],
    *,
    bin_width: float,
    fwhm_hz: float,
) -> NDArray[np.float64]:
    sigma_hz = fwhm_hz / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    sigma_bins = sigma_hz / bin_width
    radius = max(1, int(np.ceil(4.0 * sigma_bins)))
    offsets = np.arange(-radius, radius + 1, dtype=float)
    kernel = np.exp(-0.5 * (offsets / sigma_bins) ** 2)
    kernel /= np.sum(kernel)
    return np.convolve(counts, kernel, mode="same")
