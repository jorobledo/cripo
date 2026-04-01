"""Form-factor helpers translated from FACFOR/FACFOR1."""

from __future__ import annotations

import math

from .scipy_interpolation import interpolate_local_polynomial


def facfor(energ: float) -> float:
    """Interpolate the tabulated electron form-factor correction."""
    z = [0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 4.4, 4.8, 5.2, 5.6, 6.0]
    f = [
        0.00178781,
        0.00463235,
        0.0110622,
        0.0236917,
        0.0441111,
        0.0774218,
        0.13411,
        0.213965,
        0.315336,
        0.465072,
        0.651999,
        0.807181,
        0.901653,
        0.954104,
        0.983089,
        0.997777,
    ]
    al = 2.0 - math.log10(energ)
    if al < 0.0:
        return 0.0
    if al > 6.0:
        return 1.0
    return interpolate_local_polynomial(al, z, f, 5)


def facfor1(energ: float, iz: int) -> float:
    """Compute the analytic neutron-electron form-factor approximation."""
    aq0 = 1.9 * iz**0.333333
    akk = 21.968 * math.sqrt(energ)
    axx = 12.0 * (akk / aq0) ** 2
    return 2.0 * (math.sqrt(1.0 + axx) - 1.0) / axx
