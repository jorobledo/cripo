"""Debye-Waller factor translated from FADEWA."""

from __future__ import annotations

from .interpolation import ali, atsm


def fadewa(amu: float, debye: float, temp: float) -> tuple[float, int]:
    """Compute the Debye-Waller factor for the requested material state."""
    t = [
        0.0,
        0.05,
        0.10,
        0.15,
        0.2,
        0.3,
        0.4,
        0.5,
        0.6,
        0.7,
        0.8,
        0.9,
        1.0,
        1.1,
        1.2,
        1.3,
        1.4,
        1.5,
        1.6,
        1.8,
        2.0,
        2.2,
        2.5,
        3.0,
        3.5,
        4.0,
        5.0,
        6.0,
        7.0,
    ]
    f = [
        0.5,
        0.508225,
        0.532889,
        0.573583,
        0.62835,
        0.767813,
        0.931032,
        1.106947,
        1.290136,
        1.477799,
        1.668387,
        1.860982,
        2.055009,
        2.250093,
        2.445978,
        2.642484,
        2.839481,
        3.036873,
        3.234587,
        3.630769,
        4.027708,
        4.4252,
        5.022186,
        6.018497,
        7.015859,
        8.01388,
        10.011106,
        12.009257,
        14.007935,
    ]
    theta = temp / debye
    ier = 0
    if theta > 7.0:
        fi1 = 2.0 * theta
    else:
        arg, val = atsm(theta, t, f, 26, 1, 10)
        fi1, ier = ali(theta, arg, val, 10, 1.0e-4)
    fdw = 0.7276358e-14 * fi1 / amu / debye
    return fdw, ier
