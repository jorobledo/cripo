"""Packed-digit helper translated from CS1."""

from __future__ import annotations

from typing import Sequence


def cs1(data: list[int], n: int, ind: int, j1: int = 0) -> tuple[int, list[int]]:
    """Decode or update a two-digit entry in a packed integer table.

    The original Fortran stores several small integers inside a single base-10
    integer. This helper extracts the ``n``-th packed value, and when ``ind``
    is positive it also writes back an updated value using ``j1`` as the
    previous content.
    """
    m1 = (n - 1) // 4
    m2 = n - m1 * 4
    j2 = data[m1]
    m3 = 10 ** (2 * (m2 - 1))
    j3 = j2 // m3
    j4 = j3 - (j3 // 100) * 100
    if ind <= 0:
        return j4, data

    j4 = j4 - j1
    updated = list(data)
    updated[m1] = j2 - j4 * m3
    return j4, updated
