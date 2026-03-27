"""Interpolation helpers translated from ATSM/ALI."""

from __future__ import annotations

from typing import Sequence


def atsm(
    x: float,
    z: Sequence[float],
    f: Sequence[float],
    irow: int,
    icol: int,
    ndim: int,
) -> tuple[list[float], list[float]]:
    """Select interpolation nodes around ``x`` using the ATSM logic.

    This mirrors the node-selection strategy of the original IBM routine used
    by the Fortran source before applying polynomial interpolation.
    """
    arg = [0.0] * ndim
    val = [0.0] * (ndim if icol <= 1 else 2 * ndim)
    if irow <= 1:
        arg[0] = z[0]
        val[0] = f[0]
        if icol == 2 and len(f) > 1:
            val[1] = f[1]
        return arg, val

    n = min(ndim, irow)
    if z[irow - 1] >= z[0]:
        lo = 0
        hi = irow - 1
    else:
        lo = irow - 1
        hi = 0

    while abs(hi - lo) > 1:
        k = (hi + lo) // 2
        if x <= z[k]:
            hi = k
        else:
            lo = k
    j = hi if abs(z[hi] - x) <= abs(z[lo] - x) else lo

    jl = 0
    jr = 0
    k = j
    for i in range(n):
        arg[i] = z[k]
        if icol > 1:
            val[2 * i] = f[k]
            kk = k + irow
            val[2 * i + 1] = f[kk]
        else:
            val[i] = f[k]

        jjr = j + jr
        if jjr >= irow - 1:
            jl += 1
            k = j - jl
            continue
        jjl = j - jl
        if jjl <= 0 or abs(z[jjr + 1] - x) <= abs(z[jjl - 1] - x):
            jl += 1
            k = j - jl
        else:
            jr += 1
            k = j + jr
    return arg, val


def ali(
    x: float,
    arg: Sequence[float],
    val: Sequence[float],
    ndim: int,
    eps: float,
) -> tuple[float, int]:
    """Evaluate a local polynomial interpolation with the ALI algorithm.

    Returns the interpolated value and the legacy status code describing early
    convergence or interpolation issues.
    """
    values = list(val[:ndim])
    ier = 2
    delt2 = 0.0
    if ndim < 1:
        return 0.0, ier
    if ndim == 1:
        return values[0], ier

    j_result = ndim
    for j in range(1, ndim):
        delt1 = delt2
        iend = j
        for i in range(iend):
            h = arg[i] - arg[j]
            if h == 0.0:
                return values[iend - 1], 3
            values[j] = (values[i] * (x - arg[j]) - values[j] * (x - arg[i])) / h
        delt2 = abs(values[j] - values[iend - 1])
        if j < 2:
            j_result = j + 1
            continue
        if delt2 <= eps:
            return values[j], 0
        if j >= 4 and delt2 >= delt1:
            return values[iend - 1], 1
        j_result = j + 1

    return values[j_result - 1], ier
