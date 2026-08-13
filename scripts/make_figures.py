"""Regenerate the simulated figures from the transfer-matrix code.

    python scripts/make_figures.py            # all samples
    python scripts/make_figures.py M7 M8      # a subset

Writes to figures/generated/.  The experimental figures under
figures/maps/ and figures/comparisons/ came from the measurement run and
are not regenerated here: they need the raw goniometer files, which are
not distributed with this repository.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from plasmonics import (  # noqa: E402
    SAMPLES,
    default_registry,
    reflectance_map,
    reflectance_vs_angle,
)
from plasmonics.plotting import POL_COLOR, POL_LABEL, plot_curves, plot_map, plot_stack  # noqa: E402

OUT = ROOT / "figures" / "generated"

THETA = np.arange(20.0, 80.01, 0.2)
LAMBDAS_NM = np.arange(400.0, 1000.1, 4.0)


def valid_wavelengths(registry, sample, lambdas_nm=LAMBDAS_NM) -> np.ndarray:
    """Clip a wavelength grid to where every material in the stack is tabulated.

    Ta2O5 is the binding constraint: its table starts at 500 nm, so maps
    of any Ta2O5 sample cannot be drawn below that, however far the
    measurement itself extends.
    """
    lo = max(registry[m].lam_min for m in sample.media) * 1000.0
    hi = min(registry[m].lam_max for m in sample.media) * 1000.0
    return lambdas_nm[(lambdas_nm >= lo) & (lambdas_nm <= hi)]


def curve_panel(registry, sample, lam_nm: float) -> None:
    """TE and TM reflectance at one wavelength, on one axis."""
    lam_um = lam_nm / 1000.0
    n_list = registry.indices(sample.media, lam_um)
    curves = {
        POL_LABEL[pol]: (
            reflectance_vs_angle(n_list, sample.thicknesses, lam_um, THETA, pol),
            {"color": POL_COLOR[pol]},
        )
        for pol in ("p", "s")
    }

    fig, ax = plt.subplots(figsize=(7, 4.2), dpi=150)
    plot_curves(
        ax, THETA, curves,
        title=f"{sample.key} - {sample.label} - {lam_nm:.0f} nm (simulation)",
        mark_oil=False,
    )
    fig.tight_layout()
    fig.savefig(OUT / f"{sample.key.lower()}_curves_{lam_nm:.0f}nm.png")
    plt.close(fig)


def map_panel(registry, sample) -> None:
    """TE and TM reflectance maps side by side."""
    lambdas_nm = valid_wavelengths(registry, sample)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), dpi=150, sharey=True)
    for ax, pol in zip(axes, ("p", "s")):
        R = reflectance_map(registry, sample.media, sample.thicknesses, THETA, lambdas_nm / 1000.0, pol)
        im = plot_map(ax, THETA, lambdas_nm, R, title=f"{POL_LABEL[pol]}")
    fig.colorbar(im, ax=axes, label="Reflectance", pad=0.02)
    fig.suptitle(f"{sample.key} - {sample.label} - simulated reflectance", fontsize=12)
    fig.savefig(OUT / f"{sample.key.lower()}_maps.png", bbox_inches="tight")
    plt.close(fig)


def stack_sheet() -> None:
    """All six stacks drawn to scale on one sheet."""
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), dpi=150)
    for ax, sample in zip(axes.ravel(), SAMPLES.values()):
        plot_stack(ax, sample)
    fig.suptitle("Measured samples, deposited layers to scale", fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT / "sample_stacks.png", bbox_inches="tight")
    plt.close(fig)


def main(keys=None) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    registry = default_registry()
    keys = keys or list(SAMPLES)

    for key in keys:
        sample = SAMPLES[key.upper()]
        print(f"{sample.key} ...", flush=True)
        for lam_nm in (550, 633, 700):
            curve_panel(registry, sample, lam_nm)
        map_panel(registry, sample)

    stack_sheet()
    print(f"Wrote {len(list(OUT.glob('*.png')))} figures to {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main(sys.argv[1:] or None)
