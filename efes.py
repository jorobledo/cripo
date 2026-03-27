"""Combined scattering functions translated from EFES."""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .elcoh import ElcohState, elcoh
from .fadewa import fadewa
from .totinc import TotincState, totinc


@dataclass
class EfesState:
    """Mutable state shared across repeated EFES evaluations in one run."""

    fdw: float = 0.0
    idim: int = 0
    elim: float = 0.0
    totinc_state: TotincState = field(default_factory=TotincState)
    elcoh_state: ElcohState = field(default_factory=ElcohState)
    fdw_initialized: bool = False


def efes(
    energ: float,
    amu: float,
    debye: float,
    temp: float,
    ind: int,
    idim: int,
    elim: float,
    indic: int,
    a: float,
    c: float,
    state: EfesState,
) -> tuple[tuple[float, float, float, float, float, float], int, float, int, float]:
    """Evaluate the six normalized cross-section components at one energy.

    The returned tuple contains the EFES values used by the driver to build the
    physical coherent, incoherent, absorption, and total cross sections.
    """
    ier = 0
    if ind == 1 or not state.fdw_initialized:
        state.fdw, ier = fadewa(amu, debye, temp)
        state.fdw_initialized = True
        if ier == 102:
            return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0), ier, state.fdw, idim, elim

    efe1, ier = totinc(energ, amu, debye, temp, ind, state.totinc_state)
    if ier >= 100:
        return (efe1, 0.0, 0.0, 0.0, 0.0, 0.0), ier, state.fdw, idim, elim

    efe2, ier, idim_out, elim_out = elcoh(energ, state.fdw, idim, elim, indic, ind, a, c, state.elcoh_state)
    state.idim = idim_out
    state.elim = elim_out
    if ier == 103:
        return (efe1, efe2, 0.0, 0.0, 0.0, 0.0), ier, state.fdw, idim_out, elim_out

    p = 19.303706e18 * energ * state.fdw
    efe3 = (1.0 - math.exp(-p)) / p
    efe4 = 15.905676e-2 / math.sqrt(energ)
    efe5 = efe1 - efe3
    efe6 = efe1 + efe2 - efe3
    return (efe1, efe2, efe3, efe4, efe5, efe6), ier, state.fdw, idim_out, elim_out
