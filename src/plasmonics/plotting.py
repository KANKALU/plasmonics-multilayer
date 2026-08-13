"""Figures: angular curves, angle-wavelength maps, and stack diagrams."""

from __future__ import annotations

import numpy as np

from .samples import OIL_ARTEFACT_ONSET_DEG, THETA_CRITICAL_DEG

__all__ = ["plot_curves", "plot_map", "plot_stack"]

POL_LABEL = {"s": "TE (s)", "p": "TM (p)"}
POL_COLOR = {"s": "#2C6E8F", "p": "#B4643C"}


def _mark_reference_angles(ax, oil: bool = True) -> None:
    ax.axvline(THETA_CRITICAL_DEG, color="#555555", ls="--", lw=1.2, zorder=1)
    ax.annotate(
        rf"$\theta_c \approx {THETA_CRITICAL_DEG:.0f}\degree$",
        xy=(THETA_CRITICAL_DEG, 0.02),
        xytext=(3, 0),
        textcoords="offset points",
        fontsize=8,
        color="#555555",
    )
    if oil:
        ax.axvspan(OIL_ARTEFACT_ONSET_DEG, 90, color="#999999", alpha=0.12, lw=0, zorder=0)


def plot_curves(ax, theta_deg, curves: dict, title: str = "", mark_oil: bool = True):
    """Plot one or more reflectance curves against angle.

    ``curves`` maps a legend label to either an array or a
    ``(array, style_dict)`` pair.
    """
    for label, entry in curves.items():
        R, style = entry if isinstance(entry, tuple) else (entry, {})
        ax.plot(theta_deg, R, lw=1.8, label=label, **style)

    _mark_reference_angles(ax, oil=mark_oil)
    ax.set_xlabel(r"Angle of incidence $\theta$ (deg)")
    ax.set_ylabel("Reflectance")
    ax.set_ylim(0, 1.05)
    ax.grid(alpha=0.25, lw=0.6)
    ax.legend(frameon=False, fontsize=9)
    if title:
        ax.set_title(title, fontsize=11)
    return ax


def plot_map(ax, theta_deg, lambdas_nm, R, title: str = "", cmap: str = "magma"):
    """Colour map of reflectance over the angle-wavelength plane.

    Resonances read as dark bands that move to smaller angle as the
    wavelength grows, which is the signature the report reads off M7
    and M8 in TE.
    """
    TH, LAM = np.meshgrid(theta_deg, lambdas_nm)
    im = ax.pcolormesh(TH, LAM, R, shading="auto", cmap=cmap, vmin=0, vmax=1)
    ax.set_xlabel(r"Angle of incidence $\theta$ (deg)")
    ax.set_ylabel(r"Wavelength $\lambda$ (nm)")
    if title:
        ax.set_title(title, fontsize=11)
    return im


LAYER_COLORS = {
    "gold": "#C9962B",
    "copper": "#B4643C",
    "ta2o5": "#7A8CA3",
    "sio2": "#C4CCD6",
    "quartz": "#DEE3E9",
    "bk7": "#EDF0F3",
    "air": "#FFFFFF",
}


def plot_stack(ax, sample, substrate_scale: float = 0.15):
    """Draw a sample as a to-scale stack of layers.

    The 500 nm substrate is compressed by ``substrate_scale`` so the
    20-90 nm deposited layers stay legible; those layers keep their true
    relative thicknesses.
    """
    media = sample.media[1:-1]
    thicknesses = list(sample.layer_thicknesses)
    drawn = [t * substrate_scale if m == "quartz" else t for m, t in zip(media, thicknesses)]

    y = 0.0
    for medium, t_real, t_draw in zip(media, thicknesses, drawn):
        ax.add_patch(
            __import__("matplotlib.patches", fromlist=["Rectangle"]).Rectangle(
                (0, y), 1, t_draw, facecolor=LAYER_COLORS.get(medium, "#CCCCCC"),
                edgecolor="#33393F", lw=0.7,
            )
        )
        ax.text(
            1.05, y + t_draw / 2,
            f"{medium} ({t_real * 1000:.0f} nm)",
            va="center", fontsize=8,
        )
        y += t_draw

    ax.set_xlim(0, 2.2)
    ax.set_ylim(-0.02, y * 1.08)
    ax.text(0.5, y + y * 0.02, "air", ha="center", fontsize=8, color="#555555")
    ax.set_title(f"{sample.key}: {sample.label}", fontsize=10)
    ax.axis("off")
    return ax
