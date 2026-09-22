"""Geometric evaluation of the equation (11) :math:`f_{zy}` factor."""

from __future__ import annotations

from typing import cast, overload

import numpy as np
from numpy.typing import ArrayLike, NDArray

SQRT_3 = float(np.sqrt(3.0))


@overload
def f_zy(
    alpha: float,
    beta: float,
    *,
    c_yx: float = 1.0,
    c_zx: float = 1.0,
) -> float: ...


@overload
def f_zy(
    alpha: ArrayLike,
    beta: ArrayLike,
    *,
    c_yx: float = 1.0,
    c_zx: float = 1.0,
) -> NDArray[np.float64]: ...


def f_zy(
    alpha: ArrayLike,
    beta: ArrayLike,
    *,
    c_yx: float = 1.0,
    c_zx: float = 1.0,
) -> float | NDArray[np.float64]:
    """Return the equation (11) ``f_zy`` value for crystallite orientations.

    We need to get the powder distribution of this function, and we have several choices.
    The natural choice generally used would be to select the lab frame as the origin
    and then use 3 rotation angles to define each crystallite angle (probably J^PAS).
    We would then use linear algebra to get the projections.

    However, if we use the rotor frame as the origin, we can define the crystallite orientation with only 2 angles.
    What we are doing is generating each crystallite axis system *starting* from the rotor frame.
    The equation is invariant to rotations around the rotor Z axis, so the third angle does nothing
    and we don't need to include it. Note: I'm using the ZYZ convention of Arfken, Spiess, Mehring, SIMPSON, WSOLIDS, etc.

    ``alpha`` and ``beta`` are the two crystallite orientation angles in the
    rotor frame, in radians. Scalars and broadcast-compatible arrays are both
    accepted.

    The calculation deliberately follows an instructive flow for the
    geometric construction instead of jumping straight to a single simplified
    expression: rotate the crystallite axes, project them into the plane
    perpendicular to the rotor-frame ``z`` axis, prime the needed projections by
    a ``-pi/2`` rotor-axis rotation, and evaluate the projected dot products.
    """

    alpha_array = np.asarray(alpha, dtype=float)
    beta_array = np.asarray(beta, dtype=float)
    alpha_b, beta_b = np.broadcast_arrays(alpha_array, beta_array)

    if not np.isfinite(alpha_b).all() or not np.isfinite(beta_b).all():
        raise ValueError("alpha and beta must be finite")

    sin_alpha = np.sin(alpha_b)
    cos_alpha = np.cos(alpha_b)
    sin_beta = np.sin(beta_b)
    cos_beta = np.cos(beta_b)

    # Columns of R_y(beta) @ R_z(alpha). These are the crystallite x, y, and z
    # axes expressed in the rotor frame.
    x_axis = np.stack(
        (cos_alpha * cos_beta, sin_alpha, -cos_alpha * sin_beta),
        axis=-1,
    )
    y_axis = np.stack(
        (-sin_alpha * cos_beta, cos_alpha, sin_alpha * sin_beta),
        axis=-1,
    )
    z_axis = np.stack((sin_beta, np.zeros_like(sin_beta), cos_beta), axis=-1)

    # Project into the rotor-frame xy plane. Leaving the vectors unnormalized
    # preserves the original sin(xi_a) sin(xi_b) prefactor and avoids pole
    # singularities when a projected vector has zero length.
    p_x = _project_perpendicular_to_rotor_z(x_axis)
    p_y = _project_perpendicular_to_rotor_z(y_axis)
    p_z = _project_perpendicular_to_rotor_z(z_axis)

    # Priming means rotating the projected vector by -90 degrees around the
    # rotor-frame z axis: (x, y, 0) -> (y, -x, 0).
    p_x_prime = _prime_projection(p_x)
    p_y_prime = _prime_projection(p_y)

    result = (
        c_yx * _dot(p_y, p_x_prime) + c_zx * _dot(p_z, p_x_prime) + _dot(p_z, p_y_prime)
    ) / SQRT_3

    if result.shape == ():
        return float(result)
    return result


def _project_perpendicular_to_rotor_z(
    vectors: NDArray[np.float64],
) -> NDArray[np.float64]:
    projected = vectors.copy()
    projected[..., 2] = 0.0
    return projected


def _prime_projection(vectors: NDArray[np.float64]) -> NDArray[np.float64]:
    return np.stack((vectors[..., 1], -vectors[..., 0], vectors[..., 2]), axis=-1)


def _dot(left: NDArray[np.float64], right: NDArray[np.float64]) -> NDArray[np.float64]:
    return np.sum(left * right, axis=-1)


def simplified_f_zy(
    alpha: ArrayLike,
    beta: ArrayLike,
    *,
    c_yx: float = 1.0,
    c_zx: float = 1.0,
) -> NDArray[np.float64]:
    """Independent closed-form oracle used by tests."""

    alpha_array = np.asarray(alpha, dtype=float)
    beta_array = np.asarray(beta, dtype=float)
    alpha_b, beta_b = np.broadcast_arrays(alpha_array, beta_array)
    result = (
        -c_yx * np.cos(beta_b)
        + c_zx * np.sin(beta_b) * np.sin(alpha_b)
        + np.sin(beta_b) * np.cos(alpha_b)
    ) / SQRT_3
    return cast("NDArray[np.float64]", result)
