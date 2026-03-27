"""SciPy-backed interpolation utilities used by the CRIPO physics routines."""

from __future__ import annotations

from bisect import bisect_left
from typing import Sequence

from scipy.interpolate import BarycentricInterpolator


def select_nearest_window(x: float, grid: Sequence[float], values: Sequence[float], npoints: int) -> tuple[list[float], list[float]]:
    """Select a locally centered interpolation window around ``x``.

    The returned points follow the same spirit as the legacy ATSM routine:
    choose nearby nodes around the target location and keep them ordered by the
    original grid so they can be passed directly into a SciPy interpolator.
    """
    if len(grid) != len(values):
        raise ValueError("grid and values must have the same length")
    if not grid:
        raise ValueError("grid must not be empty")

    n = min(npoints, len(grid))
    if n == len(grid):
        return list(grid), list(values)

    idx = bisect_left(grid, x)
    left = idx - 1
    right = idx
    chosen: list[int] = []

    while len(chosen) < n:
        if left < 0:
            chosen.append(right)
            right += 1
            continue
        if right >= len(grid):
            chosen.append(left)
            left -= 1
            continue
        if abs(grid[left] - x) <= abs(grid[right] - x):
            chosen.append(left)
            left -= 1
        else:
            chosen.append(right)
            right += 1

    chosen.sort()
    return [grid[i] for i in chosen], [values[i] for i in chosen]


def interpolate_local_polynomial(x: float, grid: Sequence[float], values: Sequence[float], npoints: int) -> float:
    """Interpolate ``values(grid)`` at ``x`` using SciPy barycentric interpolation."""
    x_nodes, y_nodes = select_nearest_window(x, grid, values, npoints)
    interpolator = BarycentricInterpolator(x_nodes, y_nodes)
    return float(interpolator(x))
