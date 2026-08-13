"""Transfer-matrix method for the reflectance of a planar multilayer.

Geometry: ``N`` media stacked along ``z``.  The first and last are
semi-infinite; only the interior layers carry a thickness.  The
tangential wavevector is conserved across every interface, so with
``k0 = 2*pi/lambda`` and reduced quantities ``kx = n0 sin(theta0)`` the
normal component in medium ``l`` is ``kz_l = sqrt(n_l^2 - kx^2)``.

Each layer contributes

    M_l = exp(-i k0 kz_l d_l)/t_l * [[1,               r_l          ],
                                     [r_l e^{2i k0 kz_l d_l}, e^{2i k0 kz_l d_l}]]

and the stack matrix is the ordered product ``M = M_1 M_2 ... M_N``,
from which ``r = M21/M11`` and ``R = |r|^2``.  The method is exact; no
thin-film or weak-absorption approximation is made.

Polarisation is written ``'s'`` for TE and ``'p'`` for TM throughout,
matching the labels on the measured data files.
"""

from __future__ import annotations

import numpy as np

__all__ = ["reflectance", "reflectance_vs_angle", "reflectance_map", "kz_reduced", "admittances"]

_EPS = 1e-15


def _safe_div(a, b, eps: float = _EPS):
    """Divide, nudging a vanishing denominator off zero.

    Denominators here vanish only at grazing incidence and at exact
    branch points, where the reflectance is already saturated; the nudge
    keeps a sweep from producing NaNs at isolated angles.
    """
    if np.abs(b) < eps:
        b = eps + 0j if (np.iscomplexobj(a) or np.iscomplexobj(b)) else eps
    return a / b


def kz_reduced(n_list, theta0_deg: float, eps: float = 1e-12) -> list[complex]:
    """Normal wavevector components divided by ``k0``.

    The square root branch is chosen so that evanescent waves decay
    rather than grow (``Im kz >= 0``), which is what makes the method
    stable for metal layers well past the critical angle.
    """
    kx = n_list[0] * np.sin(np.deg2rad(theta0_deg))
    out = []
    for n in n_list:
        kz = np.sqrt(n**2 - kx**2 + 0j)
        if np.imag(kz) < 0 or (abs(np.imag(kz)) < 1e-15 and np.real(kz) < 0):
            kz = -kz
        if abs(kz) < eps:
            kz = eps + 0j
        out.append(kz)
    return out


def admittances(n_list, kz, pol: str) -> list[complex]:
    """Reduced optical admittances (mu = 1).

    TE sees ``Y = kz``; TM sees ``Y = eps/kz = n^2/kz``.  Casting Fresnel
    coefficients as admittance ratios keeps one interface formula for
    both polarisations.
    """
    pol = pol.lower()
    if pol == "s":
        return list(kz)
    if pol == "p":
        return [_safe_div(n**2, kzj) for n, kzj in zip(n_list, kz)]
    raise ValueError("pol must be 's' (TE) or 'p' (TM)")


def _interface(Y_l, Y_next) -> tuple[complex, complex]:
    denom = Y_l + Y_next
    return _safe_div(Y_l - Y_next, denom), _safe_div(2.0 * Y_l, denom)


def reflectance(n_list, thicknesses, lam_um: float, theta0_deg: float, pol: str) -> float:
    """Reflectance of one stack at one angle and wavelength.

    Parameters
    ----------
    n_list : sequence of complex
        Refractive index of every medium, incident side first.
    thicknesses : sequence of float
        Thickness in micrometres, same length as ``n_list``, with zeros
        at both ends for the semi-infinite media.
    lam_um : float
        Vacuum wavelength in micrometres.
    theta0_deg : float
        Angle of incidence inside the first medium, in degrees.
    pol : {'s', 'p'}
        TE or TM.
    """
    n_list = [np.complex128(n) for n in n_list]
    d = [float(x) for x in thicknesses]

    if len(d) != len(n_list):
        raise ValueError("thicknesses must have the same length as n_list")
    if d[0] != 0.0 or d[-1] != 0.0:
        raise ValueError("outer media are semi-infinite: thicknesses[0] and [-1] must be 0")
    for j, n in enumerate(n_list):
        if not np.isfinite(n.real) or not np.isfinite(n.imag) or abs(n) == 0:
            raise ValueError(f"invalid index n[{j}]={n} at lambda={lam_um} um")

    kz = kz_reduced(n_list, theta0_deg)
    Y = admittances(n_list, kz, pol)
    k0 = 2.0 * np.pi / lam_um

    M11, M12, M21, M22 = 1 + 0j, 0 + 0j, 0 + 0j, 1 + 0j
    for l in range(len(n_list) - 1):
        r_l, t_l = _interface(Y[l], Y[l + 1])
        phase = np.exp(-1j * k0 * kz[l] * d[l])
        phase2 = np.exp(2j * k0 * kz[l] * d[l])
        A = _safe_div(phase, t_l)

        m11, m12 = A, A * r_l
        m21, m22 = A * r_l * phase2, A * phase2

        M11, M12, M21, M22 = (
            M11 * m11 + M12 * m21,
            M11 * m12 + M12 * m22,
            M21 * m11 + M22 * m21,
            M21 * m12 + M22 * m22,
        )

    r = _safe_div(M21, M11)
    R = float(np.abs(r) ** 2)
    return min(max(R, 0.0), 1.0)


def reflectance_vs_angle(n_list, thicknesses, lam_um: float, thetas_deg, pol: str) -> np.ndarray:
    """Angular sweep at fixed wavelength."""
    thetas_deg = np.atleast_1d(thetas_deg)
    return np.array([reflectance(n_list, thicknesses, lam_um, th, pol) for th in thetas_deg])


def reflectance_map(registry, media, thicknesses, thetas_deg, lambdas_um, pol: str) -> np.ndarray:
    """Reflectance over the angle-wavelength plane.

    Returns an array of shape ``(len(lambdas_um), len(thetas_deg))``.  The
    indices are re-evaluated once per wavelength, which is where nearly
    all of the cost sits.
    """
    R = np.zeros((len(lambdas_um), len(thetas_deg)), dtype=float)
    for i, lam in enumerate(lambdas_um):
        n_list = registry.indices(media, lam)
        R[i, :] = reflectance_vs_angle(n_list, thicknesses, lam, thetas_deg, pol)
    return R
