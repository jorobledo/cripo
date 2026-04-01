"""Incoherent total scattering translated from TOTINC."""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .scipy_interpolation import interpolate_local_polynomial
from .totinc_data import TOTINC_A, TOTINC_B, TOTINC_C


@dataclass
class TotincState:
    """Cached interpolation state for the incoherent total-scattering model."""

    f: list[float] = field(default_factory=list)
    flim: float = 0.0
    a1: float = 0.0
    b1: float = 0.0
    initialized: bool = False


def totinc(
    e: float, amu: float, debye: float, temp: float, ind: int, state: TotincState
) -> tuple[float, int]:
    """Compute the normalized incoherent total scattering term.

    The first call initializes the temperature-dependent interpolation table.
    Later calls reuse the cached state for fast evaluation over an energy grid.
    """
    d = [0.0] * 15
    eps = 1.0e-4
    ier = 0
    boltz = 0.8617346e-4
    x = math.sqrt(e / boltz / debye)

    if ind == 1 or not state.initialized:
        coef = 1.008665 / amu
        theta = temp / debye
        f = [1.0] * 36
        for i in range(1, 4):
            for j in range(1, 37):
                for k in range(1, 16):
                    l = 540 * (i - 1) + j + 36 * (k - 1)
                    if 1080 < l < 1297:
                        d[k - 1] = 0.0
                    else:
                        if l > 1080:
                            l -= 216
                        d[k - 1] = TOTINC_A[l - 1]
                p = min(d)
                dlog = [math.log10(value - p + 0.1) for value in d]
                y = interpolate_local_polynomial(theta, TOTINC_C, dlog, 10)
                y = p + 10.0**y - 0.1
                f[j - 1] += y * coef**i

        flim = (amu / (amu + 1.008665)) ** 2
        b1 = math.log((f[34] - flim) / (f[33] - flim)) / math.log10(
            (TOTINC_B[33] / TOTINC_B[34]) ** 2
        )
        a1 = (f[34] - flim) * math.exp(
            b1 * math.log10(TOTINC_B[34] ** 2 * boltz * debye)
        )
        state.f = [math.log10(value) for value in f]
        state.flim = flim
        state.a1 = a1
        state.b1 = b1
        state.initialized = True

    if x < TOTINC_B[0]:
        return 0.0, 101
    if x <= TOTINC_B[34]:
        efe1 = interpolate_local_polynomial(x, TOTINC_B, state.f, 10)
        return 10.0**efe1, ier

    # The original file uses an undefined UFA1 symbol here. Based on the
    # routine comments and preceding algebra, this is the intended expression.
    efe1 = state.flim + state.a1 * math.exp(-state.b1 * math.log10(e))
    return efe1, ier
