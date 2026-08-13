"""Materials, samples and measurement handling.

The physics tests cover the matrix. These cover everything around it,
where the failures are quieter: a table read wrong, a wavelength silently
extrapolated, a filename parsed into the wrong sample.
"""

from __future__ import annotations

import numpy as np
import pytest

from plasmonics import SAMPLES, Dispersion, default_registry
from plasmonics.experiment import (
    align_to_model,
    angle_axis,
    find_critical_angle,
    load_measurement,
    parse_filename,
)
from plasmonics.samples import ANGLE_SWEEP_DEG, THETA_CRITICAL_DEG


# --------------------------------------------------------------------
# optical constants
# --------------------------------------------------------------------

def test_interpolation_is_exact_at_tabulated_points():
    d = Dispersion(np.array([0.4, 0.6, 0.8]), np.array([1.5, 1.6, 1.7]), name="toy")
    for wl, n in [(0.4, 1.5), (0.6, 1.6), (0.8, 1.7)]:
        assert d.n(wl) == pytest.approx(n)


def test_interpolation_is_linear_between_points():
    d = Dispersion(np.array([0.4, 0.6]), np.array([1.5, 1.7]), name="toy")
    assert d.n(0.5) == pytest.approx(1.6)


def test_lossless_material_returns_a_real_index():
    """No spurious imaginary part: it would drag the whole stack into
    complex arithmetic and hide sign errors."""
    d = Dispersion(np.array([0.4, 0.8]), np.array([1.5, 1.5]), name="toy")
    assert isinstance(d.n(0.6), float)


def test_absorbing_material_returns_a_complex_index():
    d = Dispersion(
        np.array([0.4, 0.8]), np.array([0.2, 0.3]),
        np.array([0.4, 0.8]), np.array([3.0, 3.5]), name="toy-metal",
    )
    n = d.n(0.6)
    assert isinstance(n, complex) and n.imag > 0


def test_zero_extinction_table_is_treated_as_lossless():
    d = Dispersion(
        np.array([0.4, 0.8]), np.array([1.5, 1.5]),
        np.array([0.4, 0.8]), np.array([0.0, 0.0]), name="toy",
    )
    assert isinstance(d.n(0.6), float)


def test_out_of_range_wavelength_raises_instead_of_extrapolating():
    """Silent extrapolation is the dangerous failure: numpy's interp
    clamps at the ends, which would return a plausible wrong index."""
    d = Dispersion(np.array([0.5, 0.8]), np.array([2.1, 2.0]), name="toy")
    with pytest.raises(ValueError):
        d.n(0.4)
    with pytest.raises(ValueError):
        d.n(0.9)


def test_constant_material_is_flat():
    air = Dispersion.constant(1.0, name="air")
    assert air.n(0.4) == pytest.approx(1.0)
    assert air.n(9.0) == pytest.approx(1.0)


# --------------------------------------------------------------------
# the shipped tables
# --------------------------------------------------------------------

def test_registry_loads_every_material_the_samples_need():
    registry = default_registry()
    needed = {m for s in SAMPLES.values() for m in s.media}
    assert needed <= set(registry)


def test_metals_absorb_and_dielectrics_do_not():
    registry = default_registry()
    for metal in ("gold", "copper"):
        assert isinstance(registry[metal].n(0.633), complex)
    assert registry["bk7"].n(0.633) == pytest.approx(1.5151, abs=2e-3)


def test_ta2o5_is_high_index_in_the_visible():
    """The mechanism this work points at needs a high-index dielectric;
    if this table were wrong the whole interpretation would be too."""
    assert default_registry()["ta2o5"].n(0.550).real > 2.0


def test_ta2o5_table_starts_at_500nm():
    """Documented limit: it bounds every simulated map in the repository."""
    registry = default_registry()
    assert registry["ta2o5"].lam_min == pytest.approx(0.500, abs=1e-6)
    with pytest.raises(ValueError):
        registry["ta2o5"].n(0.450)


def test_quartz_and_sio2_are_the_same_table():
    registry = default_registry()
    assert registry["quartz"].n(0.633) == registry["sio2"].n(0.633)


# --------------------------------------------------------------------
# sample definitions
# --------------------------------------------------------------------

@pytest.mark.parametrize("key", list(SAMPLES))
def test_sample_thicknesses_line_up_with_media(key):
    sample = SAMPLES[key]
    assert len(sample.thicknesses) == len(sample.media)
    assert sample.thicknesses[0] == 0.0
    assert sample.thicknesses[-1] == 0.0
    assert all(t > 0 for t in sample.layer_thicknesses)


@pytest.mark.parametrize("key", list(SAMPLES))
def test_every_sample_is_glass_in_and_air_out(key):
    sample = SAMPLES[key]
    assert sample.media[0] == "bk7"
    assert sample.media[-1] == "air"
    assert sample.media[1] == "quartz"
    assert sample.layer_thicknesses[0] == pytest.approx(0.500)


def test_mismatched_sample_definition_is_rejected():
    from plasmonics.samples import Sample

    with pytest.raises(ValueError):
        Sample(key="bad", media=("bk7", "gold", "air"), layer_thicknesses=(0.02, 0.03), label="x")


def test_the_thick_dielectric_samples_are_the_ones_the_result_is_about():
    """M6, M7 and M8 carry 70 nm or more of Ta2O5; M1, M4 and M5 do not.
    Guards the claim the report and the site are built on."""
    thick = {"M6", "M7", "M8"}
    for key, sample in SAMPLES.items():
        final_dielectric_nm = sample.layer_thicknesses[-1] * 1000
        assert (final_dielectric_nm >= 70) == (key in thick), key


# --------------------------------------------------------------------
# measurement handling
# --------------------------------------------------------------------

def test_angle_axis_matches_the_goniometer_sweep():
    theta = angle_axis()
    assert theta.size == 601
    assert theta[0] == pytest.approx(ANGLE_SWEEP_DEG[0])
    assert theta[-1] == pytest.approx(ANGLE_SWEEP_DEG[1])
    assert np.allclose(np.diff(theta), ANGLE_SWEEP_DEG[2])


@pytest.mark.parametrize(
    "name,sample,pol,lam",
    [
        ("m7S_700_Norm_Air.txt", "M7", "s", 700),
        ("m1P_550_Norm_Air.txt", "M1", "p", 550),
        ("m8S_633_Norm_Air.txt", "M8", "s", 633),
    ],
)
def test_filename_parsing(name, sample, pol, lam):
    got = parse_filename(name)
    assert got["sample"] == sample
    assert got["pol"] == pol
    assert got["wavelength_nm"] == lam


def test_filename_without_a_polarisation_is_rejected():
    with pytest.raises(ValueError):
        parse_filename("m7X_700_Norm_Air.txt")


def _synthetic_curve(theta, theta_c=39.4):
    """A measurement-shaped curve: flat, a steep rise at the critical
    angle, then a slow climb, plus a little noise."""
    rng = np.random.default_rng(0)
    R = 0.45 + 0.5 / (1 + np.exp(-(theta - theta_c) * 6))
    return R + rng.normal(0, 0.004, theta.size)


def test_critical_angle_is_found_at_the_steepest_point():
    theta = angle_axis()
    R = _synthetic_curve(theta, theta_c=39.4)
    assert find_critical_angle(theta, R) == pytest.approx(39.4, abs=0.3)


def test_alignment_puts_the_critical_angle_at_41_degrees():
    theta = angle_axis()
    R = _synthetic_curve(theta, theta_c=39.4)
    shifted, offset = align_to_model(theta, R)
    assert offset == pytest.approx(THETA_CRITICAL_DEG - 39.4, abs=0.3)
    assert find_critical_angle(shifted, R) == pytest.approx(THETA_CRITICAL_DEG, abs=1e-9)


def test_alignment_shifts_the_axis_and_leaves_the_data_alone():
    theta = angle_axis()
    R = _synthetic_curve(theta)
    shifted, offset = align_to_model(theta, R)
    assert np.allclose(shifted - theta, offset)


def test_a_deep_minimum_elsewhere_does_not_steal_the_critical_angle():
    """A plasmon dip is steeper than the critical edge; the search window
    is what stops it from being mistaken for one."""
    theta = angle_axis()
    R = _synthetic_curve(theta, theta_c=39.4)
    R -= 0.8 * np.exp(-((theta - 55.0) ** 2) / 0.5)
    assert find_critical_angle(theta, R) == pytest.approx(39.4, abs=0.3)


def test_measurement_with_the_wrong_point_count_is_rejected(tmp_path):
    """A file from a different sweep must fail loudly, not be silently
    paired with the wrong angles."""
    path = tmp_path / "short_Norm_Air.txt"
    path.write_text(" ".join("0.5" for _ in range(500)))
    with pytest.raises(ValueError):
        load_measurement(path)


def test_measurement_round_trip(tmp_path):
    theta_expected = angle_axis()
    values = np.linspace(0.4, 0.9, theta_expected.size)
    path = tmp_path / "m9S_633_Norm_Air.txt"
    path.write_text("\t".join(f"{v:.6f}" for v in values))

    theta, R = load_measurement(path)
    assert np.allclose(theta, theta_expected)
    assert np.allclose(R, values, atol=1e-6)


def test_empty_measurement_is_rejected(tmp_path):
    path = tmp_path / "empty_Norm_Air.txt"
    path.write_text("\n")
    with pytest.raises(ValueError):
        load_measurement(path)
