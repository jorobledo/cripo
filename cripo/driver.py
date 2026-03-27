"""Top-level CRIPO driver translated from Cripog.for."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

from .crystal import CrystalStructure, UnitCell, normalize_crystal_structure
from .efes import EfesState, efes
from .facfor import facfor1


@dataclass
class CripoInput:
    """Input parameters required to run one CRIPO cross-section calculation."""

    name: str
    sigma_incoherent: float
    b_coherent_fm: float
    sigma_thermal: float
    debye_temperature: float
    temperature: float
    atomic_mass_amu: float
    crystal_structure: CrystalStructure | int | str
    unit_cell: UnitCell
    lowest_lethargy_exponent: int
    highest_lethargy_exponent: int
    points_per_lethargy_decade: int
    electrons_z: int
    neutron_electron_length_fm: float = 0.0013

    def resolved_crystal_structure(self) -> CrystalStructure:
        """Return the normalized crystal structure for the current input."""
        return normalize_crystal_structure(self.crystal_structure)

    def resolved_unit_cell(self) -> UnitCell:
        """Return the unit cell that should be used by the current backend."""
        return self.unit_cell


def _mantissa_from_energy(energy: float) -> float:
    """Return the lethargy mantissa used by the original output tables."""
    al = 10.0 - math.log10(energy)
    return (al - math.floor(al)) * 0.1


def run_cripo(inputs: CripoInput, output_dir: str | Path = ".") -> dict[str, object]:
    """Run the CRIPO calculation and write the legacy `.CRI` and `.DAT` files.

    Parameters match the original interactive inputs through ``CripoInput``.
    The returned dictionary contains the in-memory cross-section rows, Debye-
    Waller factor, Bragg information, and the generated output file paths.
    """
    output_dir = Path(output_dir)
    structure = inputs.resolved_crystal_structure()
    unit_cell = inputs.resolved_unit_cell()
    cri_path = output_dir / "CRIPOOUT.CRI"
    dat_path = output_dir / "CRIPOOUT.DAT"

    bcoh = inputs.b_coherent_fm / 10.0
    sigc = 4.0 * math.pi * bcoh**2
    if sigc == 0.0:
        sigc = 1.0
    isc = 0
    bne = inputs.neutron_electron_length_fm

    a = unit_cell.a_cm
    c = unit_cell.c_cm
    elim = 1.0e5
    idim = 1000
    il0 = 8
    appd = float(inputs.points_per_lethargy_decade)
    bl = il0 - inputs.lowest_lethargy_exponent + 1.0 / appd
    ndec = inputs.highest_lethargy_exponent - inputs.lowest_lethargy_exponent
    ipun = ndec * inputs.points_per_lethargy_decade + 1
    ak = math.sqrt(sigc / 4.0 / math.pi) + bne * inputs.electrons_z / 10.0

    state = EfesState(idim=idim, elim=elim)
    rows: list[dict[str, float]] = []

    with cri_path.open("w", encoding="utf-8") as cri, dat_path.open("w", encoding="utf-8") as dat:
        cri.write("\n" + "*" * 130 + "\n\n")
        cri.write(f"{'':54}ELEMENTO - {inputs.name[:24]}\n")
        dat.write("    ENERGY        EL.COH.       EL.INC.      INEL.INC.     INEL.COH.     ABSORPTION      TOTAL\n")
        dat.write("     (EV)         (BARN)        (BARN)        (BARN)        (BARN)        (BARN)        (BARN)\n")
        cri.write(
            f" SIGMA I = {inputs.sigma_incoherent:8.4f} - SIGMA C = {sigc:8.4f} - SIGMA TH = {inputs.sigma_thermal:8.4f}"
            f" - DEBYE TEMPERATURE  = {inputs.debye_temperature:6.1f} - SAMPLE TEMPERATURE = {inputs.temperature:6.1f}\n"
        )
        cri.write(
            f" PARAMETROS DEL CRISTAL A = {unit_cell.a_angstrom:7.5f} , C = {unit_cell.c_angstrom:7.5f}"
            f" - MASA DEL DISPERSOR EN AMU = {inputs.atomic_mass_amu:8.4f}"
            f" - LONG.SCATT. N-E = {bne:8.4f} - CRYSTAL TYPE = {structure.legacy_code:4d} ({structure.display_name})\n\n"
        )
        cri.write(
            f"{'':10}ENERGY{'':8}LETARGY{'':7}EL.COH.{'':7}EL.INC.{'':6}INEL.INC.{'':5}INEL.COH."
            f"{'':5}ABSORPTION{'':6}TOTAL\n"
        )

        ind = 1
        for i in range(1, ipun + 1):
            b = bl - i / appd
            g = il0 - b
            energ = 10.0**g
            (efe1, efe2, efe3, efe4, efe5, efe6), ier, fdw, idim, elim = efes(
                energ,
                inputs.atomic_mass_amu,
                inputs.debye_temperature,
                inputs.temperature,
                ind,
                idim,
                elim,
                structure.legacy_code,
                a,
                c,
                state,
            )
            if ier >= 100:
                raise RuntimeError(f"CRIPO calculation failed with IER={ier}")

            mantissa = _mantissa_from_energy(energ)
            fbar = facfor1(energ, inputs.electrons_z)
            sigce = 4.0 * math.pi * (ak - bne * fbar * inputs.electrons_z / 10.0) ** 2
            f1 = efe2 * sigce
            f2 = efe3 * inputs.sigma_incoherent
            f3 = efe5 * inputs.sigma_incoherent
            f4 = efe5 * sigce
            f5 = efe4 * inputs.sigma_thermal
            f6 = efe1 * inputs.sigma_incoherent + efe6 * sigce + efe4 * inputs.sigma_thermal

            cri.write(
                f"{i:6d}{energ:14.6E}{mantissa:11.4f}   {f1:14.6E}{f2:14.6E}{f3:14.6E}{f4:14.6E}{f5:14.6E}{f6:14.6E}\n"
            )
            dat.write(f"{energ:14.6E}{f1:14.6E}{f2:14.6E}{f3:14.6E}{f4:14.6E}{f5:14.6E}{f6:14.6E}\n")
            rows.append(
                {
                    "energy_ev": energ,
                    "lethargy_mantissa": mantissa,
                    "elastic_coherent": f1,
                    "elastic_incoherent": f2,
                    "inelastic_incoherent": f3,
                    "inelastic_coherent": f4,
                    "absorption": f5,
                    "total": f6,
                }
            )
            ind = 0

        cri.write("\n" + "*" * 130 + "\n")
        cri.write(f" FACTOR DE DEBYE-WALLER CALCULADO = {state.fdw:13.6E}\n")
        cri.write(f" SE CALCULARON {state.idim:4d} CORTES DE BRAGG,CORRESPONDIENDO EL ULTIMO A UNA ENERGIA DE {state.elim:11.4E} EV\n")

        bragg_energies = list(state.elcoh_state.raw_a1[:51])
        bragg_jumps = [value * sigc for value in state.elcoh_state.raw_b1[:51]]

        cri.write("\n" + "*" * 130 + "\n")
        cri.write(f"{'':52}PRIMEROS CORTES DE BRAGG\n\n")
        cri.write(" NRO. ENERGIA LETARGIA       SALTO" * 3 + "\n")
        fulti = 0.0
        eulti = 0.0
        for idx in range(0, len(bragg_energies), 3):
            chunk = []
            for j in range(idx, min(idx + 3, len(bragg_energies))):
                finf = 0.0 if j == 0 else fulti * eulti / bragg_energies[j]
                leth = _mantissa_from_energy(bragg_energies[j])
                chunk.append(f"{j+1:4d}{bragg_energies[j]:10.3E}{leth:7.4f}{finf:10.3E}{bragg_jumps[j]:10.3E} *")
                fulti = bragg_jumps[j]
                eulti = bragg_energies[j]
            cri.write("".join(chunk) + "\n")

    return {
        "rows": rows,
        "fdw": state.fdw,
        "idim": state.idim,
        "elim": state.elim,
        "output_files": {"cri": str(cri_path), "dat": str(dat_path)},
        "bragg_energies": state.elcoh_state.a1,
        "bragg_jumps": [value * sigc for value in state.elcoh_state.b1],
    }
