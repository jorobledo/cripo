"""Unused cubic multiplicity helper translated from MULCUB."""

from __future__ import annotations

import math


def mulcub(is_value: int, indic: int) -> tuple[float, int]:
    """Compute cubic reciprocal-lattice multiplicities for a given index.

    This is a direct translation of the unused legacy helper and is retained
    for completeness relative to the original Fortran code base.
    """
    n = 0
    if indic != 5 and indic != 3:
        iresto = is_value - (is_value // 8) * 8
        if iresto not in (0, 3, 4):
            return 0.0, 0

    i2is = is_value if indic != 3 else 2 * is_value
    a2s = float(i2is)
    iq1sup = int(math.sqrt(a2s))
    a2s = i2is / 3.0 if (i2is // 3) * 3 != i2is else i2is / 3
    i1 = int(math.sqrt(a2s))
    iq1inf = i1 if i2is == i1 * i1 * 3 else i1 + 1

    for iq1 in range(iq1inf, iq1sup + 1):
        k2 = i2is - iq1 * iq1
        a2 = float(k2)
        iq2sup = min(iq1, int(math.sqrt(a2)))
        a2 = k2 / 2.0 if (k2 // 2) * 2 != k2 else k2 / 2
        i2 = int(math.sqrt(a2))
        iq2inf = i2 if k2 == i2 * i2 * 2 else i2 + 1
        if iq2inf > iq2sup:
            continue
        for iq2 in range(iq2inf, iq2sup + 1):
            m = 0
            iq3cua = k2 - iq2 * iq2
            iq3 = int(math.sqrt(float(iq3cua)))
            if iq3cua != iq3 * iq3 or iq3 > iq2:
                continue
            if iq3 == iq2:
                m += 1
                if iq2 == iq1:
                    n += 8
                    continue
            elif iq2 == iq1:
                m += 1
            if iq3 == 0:
                m += 1
                if iq2 == 0:
                    m += 1
            n += 48 // (2**m)
    return float(n), n
