"""Crystal-structure metadata used by the CRIPO Python API.

This module introduces a small abstraction layer above the legacy integer
``crystal_type`` codes. The current physics backend still uses the historical
Bragg-cut generators, but the public API can now carry richer crystal metadata
and a more general unit-cell description.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


ANGSTROM_TO_CM = 1.0e-8


class CrystalSystem(str, Enum):
    """Crystal-system family for a supported CRIPO structure."""

    HEXAGONAL = "hexagonal"
    CUBIC = "cubic"
    TETRAGONAL = "tetragonal"


class CrystalStructure(Enum):
    """Supported legacy CRIPO crystal structures.

    The enum values match the historical Fortran ``crystal_type`` codes so the
    current backend can still use the original branching logic.
    """

    GRAPHITE = 1
    HEXAGONAL = 2
    BCC = 3
    FCC = 4
    SIMPLE_CUBIC = 5
    LEGACY_HEXAGONAL_SPECIAL = 6
    TETRAGONAL = 7

    @property
    def legacy_code(self) -> int:
        """Return the original CRIPO integer code for the structure."""
        return int(self.value)

    @property
    def crystal_system(self) -> CrystalSystem:
        """Return the crystal-system family for the structure."""
        if self in {CrystalStructure.GRAPHITE, CrystalStructure.HEXAGONAL, CrystalStructure.LEGACY_HEXAGONAL_SPECIAL}:
            return CrystalSystem.HEXAGONAL
        if self in {CrystalStructure.BCC, CrystalStructure.FCC, CrystalStructure.SIMPLE_CUBIC}:
            return CrystalSystem.CUBIC
        return CrystalSystem.TETRAGONAL

    @property
    def display_name(self) -> str:
        """Return a human-readable structure label."""
        labels = {
            CrystalStructure.GRAPHITE: "graphite",
            CrystalStructure.HEXAGONAL: "hexagonal",
            CrystalStructure.BCC: "body-centered cubic",
            CrystalStructure.FCC: "face-centered cubic",
            CrystalStructure.SIMPLE_CUBIC: "simple cubic",
            CrystalStructure.LEGACY_HEXAGONAL_SPECIAL: "legacy hexagonal type 6",
            CrystalStructure.TETRAGONAL: "tetragonal",
        }
        return labels[self]


@dataclass(frozen=True)
class UnitCell:
    """Geometric description of a crystal unit cell.

    The current CRIPO backend only uses ``a`` and ``c`` because that is what
    the historical algorithms support. ``b`` and the cell angles are stored now
    so the public API can evolve toward more general crystal structures.
    """

    a_angstrom: float
    b_angstrom: float | None = None
    c_angstrom: float
    alpha_deg: float = 90.0
    beta_deg: float = 90.0
    gamma_deg: float = 90.0

    def __post_init__(self) -> None:
        """Fill in defaults for omitted lattice parameters."""
        if self.b_angstrom is None:
            object.__setattr__(self, "b_angstrom", self.a_angstrom)

    @property
    def a_cm(self) -> float:
        """Return the ``a`` lattice parameter in centimeters."""
        return self.a_angstrom * ANGSTROM_TO_CM

    @property
    def b_cm(self) -> float:
        """Return the ``b`` lattice parameter in centimeters."""
        return self.b_angstrom * ANGSTROM_TO_CM

    @property
    def c_cm(self) -> float:
        """Return the ``c`` lattice parameter in centimeters."""
        return self.c_angstrom * ANGSTROM_TO_CM

    @classmethod
    def from_legacy(cls, a_angstrom: float, c_angstrom: float) -> UnitCell:
        """Build a unit cell from the legacy CRIPO ``a`` and ``c`` inputs."""
        return cls(a_angstrom=a_angstrom, c_angstrom=c_angstrom)


def normalize_crystal_structure(value: int | str | CrystalStructure) -> CrystalStructure:
    """Normalize a crystal structure identifier to :class:`CrystalStructure`."""
    if isinstance(value, CrystalStructure):
        return value
    if isinstance(value, int):
        return CrystalStructure(value)

    normalized = value.strip().lower().replace("-", "_").replace(" ", "_")
    aliases = {
        "graphite": CrystalStructure.GRAPHITE,
        "hexagonal": CrystalStructure.HEXAGONAL,
        "bcc": CrystalStructure.BCC,
        "fcc": CrystalStructure.FCC,
        "cubic": CrystalStructure.SIMPLE_CUBIC,
        "simple_cubic": CrystalStructure.SIMPLE_CUBIC,
        "tetragonal": CrystalStructure.TETRAGONAL,
        "legacy_hexagonal_special": CrystalStructure.LEGACY_HEXAGONAL_SPECIAL,
        "hexagonal_special": CrystalStructure.LEGACY_HEXAGONAL_SPECIAL,
        "type_6": CrystalStructure.LEGACY_HEXAGONAL_SPECIAL,
    }
    try:
        return aliases[normalized]
    except KeyError as exc:
        raise ValueError(f"Unsupported crystal structure: {value!r}") from exc
