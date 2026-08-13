"""The six measured samples, as stacks the transfer-matrix code can read.

Every sample is a multilayer deposited on a 500 nm quartz substrate,
index-matched to the BK7 prism with immersion oil.  The oil is not part
of the modelled stack: it has its own critical angle near 70-75 deg, and
the steep drop it produces at high angle is a mounting artefact, not
sample physics.  Comparisons are therefore read between the BK7/air
critical angle (~41 deg) and roughly 70 deg.

Thicknesses are in micrometres, ordered from the prism outwards.  The
two semi-infinite media (BK7 and air) carry no thickness.
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = ["Sample", "SAMPLES", "stack_thicknesses"]


@dataclass(frozen=True)
class Sample:
    """One measured multilayer.

    Attributes
    ----------
    media : tuple of str
        Material names, incident medium first, matching the registry keys.
    layer_thicknesses : tuple of float
        Thicknesses in micrometres of the interior layers only, in the
        same order as ``media[1:-1]``.
    label : str
        How the deposited stack is written in the report.
    """

    key: str
    media: tuple[str, ...]
    layer_thicknesses: tuple[float, ...]
    label: str
    note: str = ""

    @property
    def thicknesses(self) -> list[float]:
        """Thicknesses padded with zeros for the semi-infinite media."""
        return [0.0, *self.layer_thicknesses, 0.0]

    def __post_init__(self) -> None:
        if len(self.layer_thicknesses) != len(self.media) - 2:
            raise ValueError(
                f"{self.key}: {len(self.media)} media need "
                f"{len(self.media) - 2} interior thicknesses"
            )


def stack_thicknesses(layer_thicknesses) -> list[float]:
    """Pad interior thicknesses with the two semi-infinite zeros."""
    return [0.0, *layer_thicknesses, 0.0]


SAMPLES: dict[str, Sample] = {
    "M1": Sample(
        key="M1",
        media=("bk7", "quartz", "gold", "sio2", "air"),
        layer_thicknesses=(0.500, 0.022, 0.010),
        label="Au(22 nm) / SiO₂(10 nm)",
        note="Clearest TM plasmon minimum of the set; TE stays smooth.",
    ),
    "M4": Sample(
        key="M4",
        media=("bk7", "quartz", "copper", "ta2o5", "air"),
        layer_thicknesses=(0.500, 0.020, 0.020),
        label="Cu(20 nm) / Ta₂O₅(20 nm)",
        note="Best overall agreement between measurement and model.",
    ),
    "M5": Sample(
        key="M5",
        media=("bk7", "quartz", "copper", "ta2o5", "air"),
        layer_thicknesses=(0.500, 0.020, 0.030),
        label="Cu(20 nm) / Ta₂O₅(30 nm)",
        note="M4 with a thicker dielectric; TM valley deeper than modelled.",
    ),
    "M6": Sample(
        key="M6",
        media=("bk7", "quartz", "copper", "ta2o5", "copper", "ta2o5", "air"),
        layer_thicknesses=(0.500, 0.020, 0.070, 0.020, 0.091),
        label="Cu(20 nm) / Ta₂O₅(70 nm) / Cu(20 nm) / Ta₂O₅(91 nm)",
        note="Two metal/dielectric periods; TE minimum predicted and measured.",
    ),
    "M7": Sample(
        key="M7",
        media=("bk7", "quartz", "copper", "ta2o5", "air"),
        layer_thicknesses=(0.500, 0.020, 0.070),
        label="Cu(20 nm) / Ta₂O₅(70 nm)",
        note="Narrow TE minimum in both curves; no TM plasmon resonance.",
    ),
    "M8": Sample(
        key="M8",
        media=("bk7", "quartz", "gold", "ta2o5", "air"),
        layer_thicknesses=(0.500, 0.020, 0.070),
        label="Au(20 nm) / Ta₂O₅(70 nm)",
        note="M7 with gold instead of copper; the TE minimum survives.",
    ),
}

#: Wavelengths measured, in nanometres.
MEASURED_WAVELENGTHS_NM = (550, 633, 650, 700, 750)

#: Angular sweep of the goniometer: 20 to 80 degrees in 0.1 degree steps.
ANGLE_SWEEP_DEG = (20.0, 80.0, 0.1)

#: Critical angle of the BK7/air interface, used to align measurement to model.
THETA_CRITICAL_DEG = 41.0

#: Above this angle the immersion oil dominates; comparisons stop here.
OIL_ARTEFACT_ONSET_DEG = 70.0
