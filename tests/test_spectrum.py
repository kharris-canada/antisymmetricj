from __future__ import annotations

import numpy as np
import pytest

from antisymmetricj.spectrum import calculate_spectrum, table2_transitions


def test_table2_transitions_match_hand_calculated_values() -> None:
    transitions = table2_transitions(
        np.array([0.5]),
        sigma_i_hz=60.0,
        sigma_s_hz=10.0,
        j_iso_hz=50.0,
        j_zy_anti_hz=100.0,
    )

    expected_c = np.sqrt(50.0**2 + 50.0**2 + 50.0**2)
    expected_frequencies = np.array(
        [
            0.5 * expected_c + 25.0,
            0.5 * expected_c - 25.0,
            -0.5 * expected_c + 25.0,
            -0.5 * expected_c - 25.0,
        ]
    )
    expected_intensities = np.array(
        [
            1.0 - 50.0 / expected_c,
            1.0 + 50.0 / expected_c,
            1.0 + 50.0 / expected_c,
            1.0 - 50.0 / expected_c,
        ]
    )

    np.testing.assert_allclose(transitions.frequency[:, 0], expected_frequencies)
    np.testing.assert_allclose(transitions.intensity[:, 0], expected_intensities)


def test_table2_intensities_sum_to_four_per_crystallite() -> None:
    transitions = table2_transitions(
        np.array([-1.0, -0.2, 0.4, 1.0]),
        sigma_i_hz=60.0,
        sigma_s_hz=10.0,
        j_iso_hz=50.0,
        j_zy_anti_hz=100.0,
    )

    np.testing.assert_allclose(np.sum(transitions.intensity, axis=0), 4.0)
    assert np.all(transitions.intensity >= 0.0)


def test_calculate_spectrum_uses_requested_domain_and_preserves_area() -> None:
    fzy_values = np.array([-0.5, 0.5])
    weights = np.array([2.0, 3.0])

    spectrum = calculate_spectrum(
        fzy_values,
        weights,
        sigma_i_hz=0.0,
        sigma_s_hz=0.0,
        j_iso_hz=50.0,
        j_zy_anti_hz=100.0,
        bins=512,
        fwhm_hz=1.0,
    )

    assert spectrum.left[0] == pytest.approx(-150.0)
    assert spectrum.right[-1] == pytest.approx(150.0)
    assert len(spectrum.center) == 512
    assert float(np.sum(spectrum.intensity)) == pytest.approx(4.0 * np.sum(weights))


def test_calculate_spectrum_rejects_invalid_inputs() -> None:
    with pytest.raises(ValueError, match="same shape"):
        calculate_spectrum(
            [0.0],
            [1.0, 2.0],
            sigma_i_hz=0.0,
            sigma_s_hz=0.0,
            j_iso_hz=1.0,
            j_zy_anti_hz=1.0,
        )

    with pytest.raises(ValueError, match="must be positive"):
        calculate_spectrum(
            [0.0],
            [1.0],
            sigma_i_hz=0.0,
            sigma_s_hz=0.0,
            j_iso_hz=0.0,
            j_zy_anti_hz=0.0,
        )
