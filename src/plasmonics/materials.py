"""Optical constants: loading tabulated n,k data and interpolating n(lambda).

Tabulated data files follow the plain-text export format of
RefractiveIndex.INFO: a ``wl  n`` header followed by wavelength/index
pairs, optionally followed by a ``wl  k`` header and the extinction
coefficient pairs.  Wavelengths are in micrometres.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

__all__ = ["Dispersion", "MaterialRegistry", "load_nk_file", "DATA_DIR"]

DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "refractive_index"


def load_nk_file(path: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray | None, np.ndarray | None]:
    """Parse a RefractiveIndex.INFO text export.

    Returns ``(wl_n, n, wl_k, k)``.  The ``k`` arrays are ``None`` for
    non-absorbing materials, whose files carry no ``wl k`` block.
    """
    lines = [ln.strip() for ln in Path(path).read_text(encoding="utf-8").splitlines() if ln.strip()]

    idx_n = idx_k = None
    for i, line in enumerate(lines):
        tokens = line.lower().replace("\t", " ").split()
        if len(tokens) >= 2 and tokens[0] == "wl":
            if tokens[1] == "n":
                idx_n = i
            elif tokens[1] == "k":
                idx_k = i

    if idx_n is None:
        raise ValueError(f"{path}: no 'wl n' header found")

    def _block(start: int, stop: int | None) -> tuple[np.ndarray, np.ndarray]:
        xs, ys = [], []
        for line in lines[start + 1 : stop]:
            parts = line.split()
            if len(parts) >= 2:
                xs.append(float(parts[0]))
                ys.append(float(parts[1]))
        return np.asarray(xs, float), np.asarray(ys, float)

    wl_n, n = _block(idx_n, idx_k)
    wl_k = k = None
    if idx_k is not None:
        wl_k, k = _block(idx_k, None)
        if wl_k.size == 0:
            wl_k = k = None

    return wl_n, n, wl_k, k


@dataclass
class Dispersion:
    """Wavelength-dependent complex refractive index of one material.

    ``n(lam)`` linearly interpolates the tabulated data.  A material whose
    extinction coefficient is absent or identically zero returns a real
    index, which lets the transfer-matrix code stay in real arithmetic
    when no absorbing layer is present.
    """

    wl_n: np.ndarray
    n_vals: np.ndarray
    wl_k: np.ndarray | None = None
    k_vals: np.ndarray | None = None
    name: str = "unnamed"

    _lossless: bool = field(init=False, default=True)

    def __post_init__(self) -> None:
        self.wl_n = np.asarray(self.wl_n, float)
        self.n_vals = np.asarray(self.n_vals, float)
        if self.wl_k is None or self.k_vals is None:
            self._lossless = True
        else:
            self.wl_k = np.asarray(self.wl_k, float)
            self.k_vals = np.asarray(self.k_vals, float)
            self._lossless = bool(np.allclose(self.k_vals, 0.0, atol=1e-14))

    @classmethod
    def from_file(cls, path: str | Path, name: str | None = None) -> "Dispersion":
        wl_n, n, wl_k, k = load_nk_file(path)
        return cls(wl_n, n, wl_k, k, name=name or Path(path).stem)

    @classmethod
    def constant(cls, value: float, name: str = "constant") -> "Dispersion":
        """A non-dispersive medium, e.g. air."""
        wl = np.linspace(0.1, 100.0, 2)
        return cls(wl, np.full(2, float(value)), name=name)

    @property
    def lam_min(self) -> float:
        return self.wl_n.min() if self.wl_k is None else max(self.wl_n.min(), self.wl_k.min())

    @property
    def lam_max(self) -> float:
        return self.wl_n.max() if self.wl_k is None else min(self.wl_n.max(), self.wl_k.max())

    def n(self, lam_um: float) -> complex | float:
        """Refractive index at ``lam_um`` micrometres."""
        if not (self.lam_min <= lam_um <= self.lam_max):
            raise ValueError(
                f"{self.name}: lambda={lam_um} um outside tabulated range "
                f"[{self.lam_min:.3f}, {self.lam_max:.3f}] um"
            )
        n = float(np.interp(lam_um, self.wl_n, self.n_vals))
        if self._lossless:
            return n
        k = float(np.interp(lam_um, self.wl_k, self.k_vals))
        return complex(n, k)

    def __repr__(self) -> str:
        kind = "lossless" if self._lossless else "absorbing"
        return f"Dispersion({self.name!r}, {kind}, {self.lam_min:.3f}-{self.lam_max:.3f} um)"


class MaterialRegistry(dict):
    """Named lookup of :class:`Dispersion` models, keyed by material name."""

    def register_file(self, name: str, path: str | Path) -> Dispersion:
        model = Dispersion.from_file(path, name=name)
        self[name] = model
        return model

    def indices(self, names, lam_um: float) -> list[complex | float]:
        """Refractive indices of a list of materials at one wavelength."""
        return [self[name].n(lam_um) for name in names]


def default_registry(data_dir: str | Path = DATA_DIR) -> MaterialRegistry:
    """Registry preloaded with the materials used in this work.

    The substrate quartz and the M1 dielectric SiO2 share one data file:
    both are fused silica, and the original notebooks loaded the same
    table under two different names.
    """
    data_dir = Path(data_dir)
    reg = MaterialRegistry()
    reg["air"] = Dispersion.constant(1.0, name="air")
    for name, filename in {
        "bk7": "bk7.txt",
        "gold": "gold.txt",
        "copper": "copper.txt",
        "ta2o5": "ta2o5.txt",
        "sio2": "sio2.txt",
    }.items():
        reg.register_file(name, data_dir / filename)
    reg["quartz"] = reg["sio2"]
    return reg
