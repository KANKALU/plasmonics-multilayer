"""Precompute the reflectance curves the GitHub Pages figure reads.

The page plots simulated TE and TM reflectance for every sample at every
measured wavelength.  Running the transfer matrix in the browser would
mean reimplementing it in JavaScript, so instead the curves are computed
once here and written to a JSON file the page loads.

    python scripts/export_web_data.py

Output: docs/data/curves.json
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from plasmonics import (  # noqa: E402
    MEASURED_WAVELENGTHS_NM,
    SAMPLES,
    default_registry,
    reflectance_vs_angle,
)

# Coarser than the 0.1 deg measurement grid: the page needs a smooth
# line, not the full resolution, and the payload stays small.
THETA = np.round(np.arange(20.0, 80.001, 0.25), 4)
DECIMALS = 4


def main() -> None:
    registry = default_registry()
    out = {
        "theta_deg": THETA.tolist(),
        "wavelengths_nm": list(MEASURED_WAVELENGTHS_NM),
        "samples": {},
    }

    for key, sample in SAMPLES.items():
        entry = {
            "label": sample.label,
            "note": sample.note,
            "media": list(sample.media),
            "layers_nm": [round(t * 1000) for t in sample.layer_thicknesses],
            "curves": {},
        }
        for lam_nm in MEASURED_WAVELENGTHS_NM:
            lam_um = lam_nm / 1000.0
            n_list = registry.indices(sample.media, lam_um)
            entry["curves"][str(lam_nm)] = {
                pol: np.round(
                    reflectance_vs_angle(n_list, sample.thicknesses, lam_um, THETA, pol),
                    DECIMALS,
                ).tolist()
                for pol in ("s", "p")
            }
            print(f"  {key} {lam_nm} nm")
        out["samples"][key] = entry

    dest = ROOT / "docs" / "data" / "curves.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, separators=(",", ":")), encoding="utf-8")
    print(f"\nWrote {dest.relative_to(ROOT)} ({dest.stat().st_size / 1024:.0f} kB)")


if __name__ == "__main__":
    main()
