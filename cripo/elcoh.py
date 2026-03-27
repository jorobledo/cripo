"""Elastic coherent scattering translated from ELCOH."""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .multi2 import Multi2Result, multi2


@dataclass
class ElcohState:
    """Cached Bragg-cut data and lookup state for elastic coherent scattering."""

    a1: list[float] = field(default_factory=list)
    b1: list[float] = field(default_factory=list)
    raw_a1: list[float] = field(default_factory=list)
    raw_b1: list[float] = field(default_factory=list)
    fdw: float = 0.0
    idim: int = 0
    elim: float = 0.0
    ier: int = 0
    v0: float = 0.0
    iult: int = 0
    initialized: bool = False


def elcoh(
    energ: float,
    fdw: float,
    idim: int,
    elim: float,
    indic: int,
    ind: int,
    a: float,
    c: float,
    state: ElcohState,
) -> tuple[float, int, int, float]:
    """Compute the elastic coherent contribution for a single neutron energy.

    On the first call the Bragg-cut table is generated and converted into the
    cumulative step representation used by later lookups.
    """
    if ind == 1 or not state.initialized:
        c1 = 4.3935985e9
        pi = math.pi
        const2 = -(c1 * c1) * fdw
        result = multi2(a, c, elim, indic, idim)
        if result.ier == 103:
            state.ier = result.ier
            state.a1 = result.a1
            state.b1 = result.b1
            state.idim = result.idim
            state.elim = result.elim
            state.v0 = result.v0
            state.fdw = fdw
            state.initialized = True
            return 0.0, result.ier, result.idim, result.elim

        const1 = 4.0 * pi * pi / result.v0 / c1 / c1 / c1
        eulti = 0.0
        fulti = 0.0
        b1_vals = list(result.b1)
        for i in range(result.idim):
            finf = eulti * fulti / result.a1[i]
            salto = const1 * b1_vals[i] * math.exp(const2 * result.a1[i]) / math.sqrt(result.a1[i] ** 3)
            b1_vals[i] = finf + salto
            eulti = result.a1[i]
            fulti = b1_vals[i]

        state.a1 = list(result.a1)
        state.b1 = b1_vals
        state.raw_a1 = list(result.a1)
        state.raw_b1 = list(result.b1)
        state.fdw = fdw
        state.idim = result.idim
        state.elim = result.elim
        state.ier = result.ier
        state.v0 = result.v0
        state.iult = 0
        state.initialized = True

    if energ > state.elim:
        const2 = -(4.3935985e9 ** 2) * state.fdw
        const3 = const2 * energ
        efe2 = -(1.0 - math.exp(const3)) / const3
        return efe2, state.ier, state.idim, state.elim

    if energ < state.a1[0]:
        return 0.0, state.ier, state.idim, state.elim

    if state.iult >= len(state.a1):
        state.iult = len(state.a1) - 1

    while state.iult > 0 and energ < state.a1[state.iult]:
        state.iult -= 1
    while state.iult + 1 < len(state.a1) and energ > state.a1[state.iult + 1]:
        state.iult += 1

    efe2 = state.b1[state.iult] * state.a1[state.iult] / energ
    return efe2, state.ier, state.idim, state.elim
