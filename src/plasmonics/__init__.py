"""Reflectance of periodic metal-oxide multilayers.

Transfer-matrix simulation and measurement analysis for surface plasmon
resonance in the Kretschmann geometry, and for the TE-polarised
reflectance minimum that appears once the final dielectric layer is
thick enough.

Typical use::

    from plasmonics import default_registry, SAMPLES, reflectance_vs_angle
    import numpy as np

    reg = default_registry()
    m7 = SAMPLES["M7"]
    theta = np.arange(20, 80.1, 0.1)
    n = reg.indices(m7.media, 0.700)
    R_te = reflectance_vs_angle(n, m7.thicknesses, 0.700, theta, "s")
"""

from .materials import Dispersion, MaterialRegistry, default_registry, load_nk_file
from .samples import (
    ANGLE_SWEEP_DEG,
    MEASURED_WAVELENGTHS_NM,
    OIL_ARTEFACT_ONSET_DEG,
    SAMPLES,
    THETA_CRITICAL_DEG,
    Sample,
)
from .tmm import admittances, kz_reduced, reflectance, reflectance_map, reflectance_vs_angle

__version__ = "1.0.0"

__all__ = [
    "Dispersion",
    "MaterialRegistry",
    "default_registry",
    "load_nk_file",
    "Sample",
    "SAMPLES",
    "MEASURED_WAVELENGTHS_NM",
    "ANGLE_SWEEP_DEG",
    "THETA_CRITICAL_DEG",
    "OIL_ARTEFACT_ONSET_DEG",
    "reflectance",
    "reflectance_vs_angle",
    "reflectance_map",
    "kz_reduced",
    "admittances",
    "__version__",
]
