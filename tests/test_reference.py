"""The matrix method against a second, independent algorithm.

Everything in `test_tmm.py` checks the code against closed-form answers,
which only exist for one interface and one film. The samples in this work
have up to six media, and there the analytic route runs out.

So this file solves the same problem a different way. The recursive
Fresnel method walks the stack from the far side inwards, folding each
layer into a running reflection amplitude, and never forms a matrix at
all. Two algorithms this different agreeing to twelve digits on hundreds
of randomly generated stacks is much stronger evidence than either one
agreeing with itself.

The reference implementation lives here rather than in the package
because its only job is to be independent. If a bug ever lands in
`tmm.py`, importing shared code here would hide it.
"""

from __future__ import annotations

import numpy as np
import pytest

from plasmonics import SAMPLES, default_registry, reflectance


def reflectance_recursive(n_list, thicknesses, lam_um, theta_deg, pol) -> float:
    """Recursive Fresnel reflectance, written from scratch.

    Starting at the last interface and working back toward the incident
    medium, each layer folds in through the Airy relation

        r ← (r_{j-1,j} + r e^{2iβ_j}) / (1 + r_{j-1,j} r e^{2iβ_j})

    with β_j = k₀ n_j cos θ_j d_j. No transfer matrices anywhere.
    """
    n = [complex(x) for x in n_list]
    d = [float(x) for x in thicknesses]
    k0 = 2 * np.pi / lam_um

    kx = n[0] * np.sin(np.deg2rad(theta_deg))
    cos = []
    for nj in n:
        c = np.sqrt(1 - (kx / nj) ** 2 + 0j)
        cos.append(c if c.imag >= 0 else -c)

    def r_interface(i, j):
        if pol == "s":
            return (n[i] * cos[i] - n[j] * cos[j]) / (n[i] * cos[i] + n[j] * cos[j])
        return (n[j] * cos[i] - n[i] * cos[j]) / (n[j] * cos[i] + n[i] * cos[j])

    N = len(n)
    r = r_interface(N - 2, N - 1)
    for j in range(N - 2, 0, -1):
        e = np.exp(2j * k0 * n[j] * cos[j] * d[j])
        r_prev = r_interface(j - 1, j)
        r = (r_prev + r * e) / (1 + r_prev * r * e)

    return float(abs(r) ** 2)


CU = complex(0.21, 3.67)   # copper near 633 nm
AU = complex(0.19, 3.09)   # gold near 633 nm
GLASS = 1.5168
TA2O5 = 2.1
LAM = 0.633


HAND_PICKED = [
    ([GLASS, 1.0], [0, 0]),
    ([GLASS, TA2O5, 1.0], [0, 0.070, 0]),
    ([GLASS, CU, 1.0], [0, 0.020, 0]),
    ([GLASS, CU, TA2O5, 1.0], [0, 0.020, 0.070, 0]),
    ([GLASS, AU, TA2O5, 1.0], [0, 0.020, 0.070, 0]),
    ([GLASS, 1.4585, CU, TA2O5, 1.0], [0, 0.500, 0.020, 0.070, 0]),
    ([GLASS, 1.4585, CU, TA2O5, CU, TA2O5, 1.0], [0, 0.500, 0.020, 0.070, 0.020, 0.091, 0]),
]


@pytest.mark.parametrize("n_list,thicknesses", HAND_PICKED)
@pytest.mark.parametrize("theta", [0.0, 15.0, 33.0, 41.0, 48.0, 60.0, 72.0, 85.0])
@pytest.mark.parametrize("pol", ["s", "p"])
def test_matches_recursive_method(n_list, thicknesses, theta, pol):
    """Every stack in this work, both polarisations, across the sweep."""
    matrix = reflectance(n_list, thicknesses, LAM, theta, pol)
    recursive = reflectance_recursive(n_list, thicknesses, LAM, theta, pol)
    assert matrix == pytest.approx(recursive, abs=1e-11)


@pytest.mark.parametrize("seed", range(40))
def test_matches_recursive_method_on_random_stacks(seed):
    """Randomly generated stacks, absorbing layers included.

    Random inputs reach corners hand-picked cases miss: near-degenerate
    indices, layers thin enough to be nearly transparent, angles right on
    a critical angle.
    """
    rng = np.random.default_rng(seed)
    n_layers = int(rng.integers(1, 6))

    n_list = [float(rng.uniform(1.3, 1.9))]
    thicknesses = [0.0]
    for _ in range(n_layers):
        if rng.random() < 0.4:                       # an absorbing layer
            n_list.append(complex(rng.uniform(0.1, 1.5), rng.uniform(0.5, 5.0)))
        else:                                        # a transparent one
            n_list.append(float(rng.uniform(1.2, 2.6)))
        thicknesses.append(float(rng.uniform(0.005, 0.4)))
    n_list.append(1.0)
    thicknesses.append(0.0)

    for theta in rng.uniform(0, 89, size=6):
        for pol in ("s", "p"):
            matrix = reflectance(n_list, thicknesses, LAM, float(theta), pol)
            recursive = reflectance_recursive(n_list, thicknesses, LAM, float(theta), pol)
            assert matrix == pytest.approx(recursive, abs=1e-10), (
                f"seed={seed} theta={theta:.2f} pol={pol} n={n_list} d={thicknesses}"
            )


@pytest.mark.parametrize("key", list(SAMPLES))
@pytest.mark.parametrize("lam_nm", [550, 633, 750])
def test_real_samples_match_recursive_method(key, lam_nm):
    """The six measured stacks, with their tabulated optical constants."""
    registry = default_registry()
    sample = SAMPLES[key]
    lam_um = lam_nm / 1000.0
    n_list = registry.indices(sample.media, lam_um)

    for theta in (25.0, 41.0, 45.0, 55.0, 65.0, 78.0):
        for pol in ("s", "p"):
            matrix = reflectance(n_list, sample.thicknesses, lam_um, theta, pol)
            recursive = reflectance_recursive(n_list, sample.thicknesses, lam_um, theta, pol)
            assert matrix == pytest.approx(recursive, abs=1e-10), f"{key} {lam_nm}nm {theta}deg {pol}"
