from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable

import mpmath as mp


DEFAULT_DPS = 120


def configure_precision(dps: int = DEFAULT_DPS) -> None:
    """Set the global mpmath working precision."""
    mp.mp.dps = int(dps)


def to_mpf(value: str | float | int | mp.mpf) -> mp.mpf:
    """Convert a scalar to mpmath's arbitrary-precision floating type."""
    return mp.mpf(value)


@dataclass(frozen=True)
class NumberFormatOptions:
    """User-facing numerical formatting options for terminal tables."""

    digits: int = 24
    min_fixed_exponent: int = -6
    max_fixed_exponent: int = 24


DEFAULT_NUMBER_FORMAT = NumberFormatOptions()


def format_mpf(
    value: mp.mpf | float | int,
    digits: int = DEFAULT_NUMBER_FORMAT.digits,
    *,
    min_fixed_exponent: int = DEFAULT_NUMBER_FORMAT.min_fixed_exponent,
    max_fixed_exponent: int | None = None,
) -> str:
    """
    Format a scalar without losing the leading zero for sub-unit numbers.

    The previous implementation delegated too aggressively to scientific
    notation, which becomes visually dangerous when a narrow terminal truncates
    the exponent field.  The present formatter keeps ordinary fixed-point
    output for physically typical magnitudes while still falling back to a
    compact scientific string for extremely small or extremely large values.
    """
    x = mp.mpf(value)
    if x == 0:
        return "0"

    if max_fixed_exponent is None:
        max_fixed_exponent = int(digits)

    try:
        exponent = int(mp.floor(mp.log10(abs(x))))
    except ValueError:
        exponent = 0

    if min_fixed_exponent <= exponent <= max_fixed_exponent:
        return mp.nstr(
            x,
            n=int(digits),
            strip_zeros=False,
            min_fixed=int(min_fixed_exponent),
            max_fixed=int(max_fixed_exponent),
        )
    return mp.nstr(x, n=int(digits), strip_zeros=False)


def format_percent(value: mp.mpf | float | int, digits: int = 6) -> str:
    """Format a dimensionless fraction as a percentage string."""
    return f"{float(mp.mpf(value) * 100):.{digits}f}%"


def quadratic_form(matrix: mp.matrix, vector: mp.matrix) -> mp.mpf:
    """Return v^T M v for real column vectors."""
    return (vector.T * matrix * vector)[0]


def bracket_root(
    function: Callable[[mp.mpf], mp.mpf],
    left: mp.mpf,
    right: mp.mpf,
    growth: mp.mpf = mp.mpf("2.0"),
    max_expansions: int = 128,
) -> tuple[mp.mpf, mp.mpf]:
    """
    Expand a positive bracket until a sign change is found.

    The routine assumes the physical root is positive and starts from the
    user-supplied interval [left, right].
    """
    f_left = function(left)
    f_right = function(right)
    expansions = 0
    while f_left * f_right > 0:
        right *= growth
        f_right = function(right)
        expansions += 1
        if expansions > max_expansions:
            raise RuntimeError(
                "Unable to bracket a positive root inside the allowed search window."
            )
    return left, right


def bisect_root(
    function: Callable[[mp.mpf], mp.mpf],
    left: mp.mpf,
    right: mp.mpf,
    *,
    absolute_tolerance: mp.mpf = mp.mpf("1e-70"),
    max_iterations: int = 1024,
) -> mp.mpf:
    """High-precision bisection for a sign-changing interval."""
    f_left = function(left)
    f_right = function(right)
    if f_left == 0:
        return left
    if f_right == 0:
        return right
    if f_left * f_right > 0:
        raise ValueError("Bisection requires a sign-changing bracket.")
    for _ in range(max_iterations):
        midpoint = (left + right) / 2
        f_mid = function(midpoint)
        if abs(f_mid) < absolute_tolerance or abs(right - left) < absolute_tolerance:
            return midpoint
        if f_left * f_mid <= 0:
            right = midpoint
            f_right = f_mid
        else:
            left = midpoint
            f_left = f_mid
    return (left + right) / 2


def sum_mpf(values: Iterable[mp.mpf]) -> mp.mpf:
    """Stable mpmath-based sum."""
    total = mp.mpf("0")
    for value in values:
        total += value
    return total


@dataclass(frozen=True)
class ResidualCheck:
    """Small container used in printed diagnostics."""

    value: mp.mpf
    description: str
