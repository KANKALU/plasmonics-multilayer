"""The transfer matrix, checked against cases with closed-form answers.

The point of these is that none of the expected values come from running
this code. They come from Fresnel's equations, from the Airy formula for
a single film, and from symmetries the method has to respect no matter
how it is implemented. A regression test that compares today's output to
yesterday's output would pass just as happily with the physics wrong.

`test_reference.py` carries the other half: the same problem solved by a
second, independent algorithm, which is what covers the multilayer cases
where no closed form exists.

Mutation-tested. Every change that alters the physics — the branch of
kz, the TM admittance, the round-trip phase factor, which matrix element
the amplitude is read from, sin/cos in the tangential wavevector, the
sign of the extinction coefficient — makes tests fail. Two deliberate
changes survive, and should:

  * flipping the sign of r at every interface multiplies the total
    amplitude by -1, and
  * dropping the factor 2 from t scales the whole stack matrix by a
    constant,

and reflectance reads |r|^2 = |M21/M11|^2, which is blind to both. They
would matter for transmittance, which this module does not compute.

One caution learned the hard way: when swapping implementations to check
coverage, run with PYTHONDONTWRITEBYTECODE=1. Edits that keep the file
the same size can leave a stale .pyc in place, and the suite then passes
against bytecode that no longer matches the source.
"""

from __future__ import annotations

import numpy as np
import pytest

from plasmonics import reflectance, reflectance_vs_angle

N_GLASS = 1.5168  # BK7 near 633 nm
N_AIR = 1.0
LAM = 0.633


# --------------------------------------------------------------------
# helpers: the analytic answers
# --------------------------------------------------------------------

def fresnel(n1: float, n2: float, theta_deg: float, pol: str) -> float:
    """Reflectance of a single interface, straight from Fresnel."""
    th_i = np.deg2rad(theta_deg)
    sin_t = (n1 / n2) * np.sin(th_i)
    cos_t = np.sqrt(1 - sin_t**2 + 0j)
    cos_i = np.cos(th_i)

    if pol == "s":
        r = (n1 * cos_i - n2 * cos_t) / (n1 * cos_i + n2 * cos_t)
    else:
        r = (n2 * cos_i - n1 * cos_t) / (n2 * cos_i + n1 * cos_t)
    return float(abs(r) ** 2)


def airy_single_film(n0: float, n1: float, n2: float, d_um: float, lam_um: float) -> float:
    """Reflectance of one film between two media, at normal incidence.

    The standard two-interface interference result: the film's two
    reflections beat against each other with a round-trip phase.
    """
    r01 = (n0 - n1) / (n0 + n1)
    r12 = (n1 - n2) / (n1 + n2)
    beta = 2 * np.pi * n1 * d_um / lam_um
    r = (r01 + r12 * np.exp(2j * beta)) / (1 + r01 * r12 * np.exp(2j * beta))
    return float(abs(r) ** 2)


def _fresnel_r(n1, n2, cos1, cos2, pol):
    if pol == "s":
        return (n1 * cos1 - n2 * cos2) / (n1 * cos1 + n2 * cos2)
    return (n2 * cos1 - n1 * cos2) / (n2 * cos1 + n1 * cos2)


def airy_at_angle(n0, n1, n2, d_um, lam_um, theta_deg, pol) -> float:
    """Airy formula for one film at arbitrary incidence.

    Same interference of two interface reflections, but the coefficients
    are now angle- and polarisation-dependent and the round-trip phase
    carries cos(theta) inside the film. Unlike the normal-incidence case,
    this distinguishes TE from TM and is sensitive to which matrix
    element the reflection amplitude is read from.
    """
    cos0 = np.cos(np.deg2rad(theta_deg))
    sin0 = np.sin(np.deg2rad(theta_deg))
    cos1 = np.sqrt(1 - (n0 / n1 * sin0) ** 2 + 0j)
    cos2 = np.sqrt(1 - (n0 / n2 * sin0) ** 2 + 0j)

    r01 = _fresnel_r(n0, n1, cos0, cos1, pol)
    r12 = _fresnel_r(n1, n2, cos1, cos2, pol)
    beta = 2 * np.pi * n1 * cos1 * d_um / lam_um
    r = (r01 + r12 * np.exp(2j * beta)) / (1 + r01 * r12 * np.exp(2j * beta))
    return float(abs(r) ** 2)


BREWSTER = np.rad2deg(np.arctan(N_AIR / N_GLASS))
CRITICAL = np.rad2deg(np.arcsin(N_AIR / N_GLASS))


# --------------------------------------------------------------------
# single interface
# --------------------------------------------------------------------

@pytest.mark.parametrize("pol", ["s", "p"])
@pytest.mark.parametrize("theta", [0.0, 10.0, 25.0, 38.0, 41.5, 55.0, 75.0])
def test_bare_interface_matches_fresnel(pol, theta):
    """With no layers at all, the method must reduce to Fresnel exactly."""
    got = reflectance([N_GLASS, N_AIR], [0.0, 0.0], LAM, theta, pol)
    assert got == pytest.approx(fresnel(N_GLASS, N_AIR, theta, pol), abs=1e-9)


def test_brewster_angle_kills_tm():
    """TM reflectance vanishes at the Brewster angle; TE does not."""
    r_tm = reflectance([N_GLASS, N_AIR], [0.0, 0.0], LAM, BREWSTER, "p")
    r_te = reflectance([N_GLASS, N_AIR], [0.0, 0.0], LAM, BREWSTER, "s")
    assert r_tm == pytest.approx(0.0, abs=1e-12)
    assert r_te > 0.1


def test_no_brewster_analogue_in_te():
    """TE has no zero anywhere: the minimum sits at normal incidence."""
    theta = np.linspace(0, 89.9, 900)
    R = reflectance_vs_angle([N_GLASS, N_AIR], [0.0, 0.0], LAM, theta, "s")
    assert R.min() == pytest.approx(R[0], abs=1e-12)
    assert R.min() > 0.0


@pytest.mark.parametrize("pol", ["s", "p"])
def test_total_internal_reflection(pol):
    """Past the critical angle a lossless interface reflects everything."""
    for theta in (CRITICAL + 0.5, 55.0, 70.0, 85.0):
        R = reflectance([N_GLASS, N_AIR], [0.0, 0.0], LAM, theta, pol)
        assert R == pytest.approx(1.0, abs=1e-9)


def test_below_critical_angle_transmits():
    """Below it, some light must get through."""
    for theta in (0.0, 20.0, CRITICAL - 0.5):
        for pol in ("s", "p"):
            assert reflectance([N_GLASS, N_AIR], [0.0, 0.0], LAM, theta, pol) < 1.0


def test_polarisations_agree_at_normal_incidence():
    """At theta = 0 the two polarisations are indistinguishable."""
    r_te = reflectance([N_GLASS, N_AIR], [0.0, 0.0], LAM, 0.0, "s")
    r_tm = reflectance([N_GLASS, N_AIR], [0.0, 0.0], LAM, 0.0, "p")
    expected = ((N_GLASS - N_AIR) / (N_GLASS + N_AIR)) ** 2
    assert r_te == pytest.approx(expected, abs=1e-12)
    assert r_tm == pytest.approx(expected, abs=1e-12)


# --------------------------------------------------------------------
# one film
# --------------------------------------------------------------------

@pytest.mark.parametrize("d_um", [0.02, 0.07, 0.1055, 0.25])
def test_single_film_matches_airy(d_um):
    """One dielectric film at normal incidence, against the Airy formula."""
    n_film = 2.1  # Ta2O5 in the visible, near enough
    got = reflectance([N_GLASS, n_film, N_AIR], [0.0, d_um, 0.0], LAM, 0.0, "s")
    assert got == pytest.approx(airy_single_film(N_GLASS, n_film, N_AIR, d_um, LAM), abs=1e-9)


def test_quarter_wave_layer_is_antireflective_when_index_matched():
    """A quarter-wave film of index sqrt(n0 n2) reflects nothing.

    This is the textbook single-layer AR coating, and it exercises the
    phase term rather than just the interface coefficients.
    """
    n0, n2 = N_GLASS, N_AIR
    n_film = np.sqrt(n0 * n2)
    d = LAM / (4 * n_film)
    for pol in ("s", "p"):
        R = reflectance([n0, n_film, n2], [0.0, d, 0.0], LAM, 0.0, pol)
        assert R == pytest.approx(0.0, abs=1e-12)


def test_half_wave_layer_is_invisible():
    """A half-wave film restores the bare-interface result."""
    n_film = 2.1
    d = LAM / (2 * n_film)
    bare = reflectance([N_GLASS, N_AIR], [0.0, 0.0], LAM, 0.0, "s")
    with_film = reflectance([N_GLASS, n_film, N_AIR], [0.0, d, 0.0], LAM, 0.0, "s")
    assert with_film == pytest.approx(bare, abs=1e-9)


@pytest.mark.parametrize("pol", ["s", "p"])
@pytest.mark.parametrize("theta", [5.0, 18.0, 33.0, 40.0, 47.0, 62.0, 78.0])
@pytest.mark.parametrize("d_um", [0.020, 0.070, 0.190])
def test_single_film_matches_airy_at_angle(pol, theta, d_um):
    """One film at oblique incidence, against the angle-resolved Airy formula.

    The strongest analytic check here: it separates the two
    polarisations, spans both sides of the critical angle, and fixes the
    reflection amplitude down to which matrix element it comes from.
    """
    n_film = 2.1
    got = reflectance([N_GLASS, n_film, N_AIR], [0.0, d_um, 0.0], LAM, theta, pol)
    assert got == pytest.approx(airy_at_angle(N_GLASS, n_film, N_AIR, d_um, LAM, theta, pol), abs=1e-9)


# --------------------------------------------------------------------
# structural invariances
# --------------------------------------------------------------------

@pytest.mark.parametrize("theta", [0.0, 30.0, 55.0])
@pytest.mark.parametrize("pol", ["s", "p"])
def test_zero_thickness_layer_is_a_no_op(theta, pol):
    """Inserting a layer of zero thickness must change nothing."""
    without = reflectance([N_GLASS, N_AIR], [0.0, 0.0], LAM, theta, pol)
    with_ghost = reflectance([N_GLASS, 3.4, N_AIR], [0.0, 0.0, 0.0], LAM, theta, pol)
    assert with_ghost == pytest.approx(without, abs=1e-9)


@pytest.mark.parametrize("theta", [0.0, 20.0, 35.0])
@pytest.mark.parametrize("pol", ["s", "p"])
def test_index_matched_layer_is_invisible(theta, pol):
    """A layer with its neighbour's index is not an interface at all."""
    plain = reflectance([N_GLASS, N_AIR], [0.0, 0.0], LAM, theta, pol)
    padded = reflectance([N_GLASS, N_GLASS, N_AIR], [0.0, 0.3, 0.0], LAM, theta, pol)
    assert padded == pytest.approx(plain, abs=1e-9)


@pytest.mark.parametrize("pol", ["s", "p"])
def test_splitting_a_layer_in_two_changes_nothing(pol):
    """One 70 nm layer equals two 35 nm layers of the same material."""
    one = reflectance([N_GLASS, 2.1, N_AIR], [0.0, 0.070, 0.0], LAM, 45.0, pol)
    two = reflectance([N_GLASS, 2.1, 2.1, N_AIR], [0.0, 0.035, 0.035, 0.0], LAM, 45.0, pol)
    assert two == pytest.approx(one, abs=1e-9)


# --------------------------------------------------------------------
# absorbing layers: the regime the experiment actually lives in
# --------------------------------------------------------------------

CU_633 = complex(0.21, 3.67)  # copper, Johnson and Christy, near 633 nm


@pytest.mark.parametrize("theta", np.arange(20.0, 89.5, 2.5))
@pytest.mark.parametrize("pol", ["s", "p"])
def test_metal_stack_stays_physical_past_critical_angle(theta, pol):
    """A passive stack can never reflect more than it receives.

    Necessary but far from sufficient: several wrong implementations
    still land inside [0, 1]. The two tests below are what actually pin
    the branch choice down.
    """
    R = reflectance([N_GLASS, CU_633, 2.1, N_AIR], [0.0, 0.020, 0.070, 0.0], LAM, theta, pol)
    assert np.isfinite(R)
    assert 0.0 <= R <= 1.0


@pytest.mark.parametrize("theta", [0.0, 30.0, 55.0, 75.0])
@pytest.mark.parametrize("pol", ["s", "p"])
def test_thick_absorbing_film_converges_to_bulk_interface(theta, pol):
    """A metal film thick enough to swallow the light behaves like bulk metal.

    Once the film is many skin depths deep, nothing reaches the far
    interface, so the answer must collapse onto the bare glass/metal
    Fresnel value. This is the test that pins the kz branch: take the
    other root and the evanescent field grows through the film instead
    of decaying, so a thicker film moves the result *away* from bulk.
    """
    bulk = fresnel(N_GLASS, CU_633, theta, pol)
    thick = reflectance([N_GLASS, CU_633, N_AIR], [0.0, 2.0, 0.0], LAM, theta, pol)
    assert thick == pytest.approx(bulk, abs=1e-6)


@pytest.mark.parametrize("pol", ["s", "p"])
def test_frustrated_total_reflection_decays_with_gap(pol):
    """Widening an air gap past the critical angle drives R monotonically to 1.

    The evanescent tail is what couples across the gap, so its decay has
    to show up as a monotone approach to total reflection. With the wrong
    branch the tail grows and this sequence runs the other way.
    """
    theta = 55.0  # comfortably past the 41 deg critical angle
    gaps = [0.05, 0.10, 0.20, 0.40, 0.80]
    R = [
        reflectance([N_GLASS, N_AIR, N_GLASS], [0.0, g, 0.0], LAM, theta, pol)
        for g in gaps
    ]
    assert all(b > a for a, b in zip(R, R[1:])), R
    assert R[-1] == pytest.approx(1.0, abs=1e-4)
    assert R[0] < 0.99


def test_thick_metal_reflects_almost_everything():
    """Bulk copper is a mirror: 200 nm leaves nothing to transmit."""
    R = reflectance([N_GLASS, CU_633, N_AIR], [0.0, 0.200, 0.0], LAM, 0.0, "s")
    assert R > 0.9


def test_absorbing_layer_breaks_total_internal_reflection():
    """With a lossy film present, R < 1 even above the critical angle."""
    R = reflectance([N_GLASS, CU_633, N_AIR], [0.0, 0.020, 0.0], LAM, 55.0, "p")
    assert R < 1.0


# --------------------------------------------------------------------
# input validation
# --------------------------------------------------------------------

def test_outer_media_must_be_semi_infinite():
    with pytest.raises(ValueError):
        reflectance([N_GLASS, N_AIR], [0.5, 0.0], LAM, 45.0, "s")


def test_thickness_count_must_match_media():
    with pytest.raises(ValueError):
        reflectance([N_GLASS, 2.1, N_AIR], [0.0, 0.0], LAM, 45.0, "s")


def test_unknown_polarisation_is_rejected():
    with pytest.raises(ValueError):
        reflectance([N_GLASS, N_AIR], [0.0, 0.0], LAM, 45.0, "circular")


def test_sweep_returns_one_value_per_angle():
    theta = np.linspace(20, 80, 37)
    R = reflectance_vs_angle([N_GLASS, N_AIR], [0.0, 0.0], LAM, theta, "s")
    assert R.shape == theta.shape
