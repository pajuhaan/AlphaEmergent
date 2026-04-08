from __future__ import annotations

import mpmath as mp


Series = list[mp.mpf]


def zeros(order: int) -> Series:
    return [mp.mpf("0") for _ in range(order + 1)]


def one(order: int) -> Series:
    out = zeros(order)
    out[0] = mp.mpf("1")
    return out


def add(a: Series, b: Series, order: int) -> Series:
    return [
        (a[i] if i < len(a) else mp.mpf("0"))
        + (b[i] if i < len(b) else mp.mpf("0"))
        for i in range(order + 1)
    ]


def sub(a: Series, b: Series, order: int) -> Series:
    return [
        (a[i] if i < len(a) else mp.mpf("0"))
        - (b[i] if i < len(b) else mp.mpf("0"))
        for i in range(order + 1)
    ]


def scale(a: Series, scalar: mp.mpf, order: int) -> Series:
    return [
        (a[i] if i < len(a) else mp.mpf("0")) * scalar
        for i in range(order + 1)
    ]


def multiply(a: Series, b: Series, order: int) -> Series:
    out = zeros(order)
    for i in range(min(order, len(a) - 1) + 1):
        ai = a[i]
        if ai == 0:
            continue
        for j in range(min(order - i, len(b) - 1) + 1):
            out[i + j] += ai * b[j]
    return out


def shift(a: Series, shift_by: int, order: int) -> Series:
    out = zeros(order)
    for i in range(min(len(a), order + 1)):
        j = i + shift_by
        if 0 <= j <= order:
            out[j] = a[i]
    return out


def reciprocal(a: Series, order: int) -> Series:
    if a[0] == 0:
        raise ZeroDivisionError("A series reciprocal requires a nonzero constant term.")
    out = zeros(order)
    out[0] = 1 / a[0]
    for m in range(1, order + 1):
        accumulator = mp.mpf("0")
        for k in range(1, m + 1):
            accumulator += (a[k] if k < len(a) else mp.mpf("0")) * out[m - k]
        out[m] = -accumulator / a[0]
    return out


def sqrt(a: Series, order: int) -> Series:
    if a[0] == 0:
        raise ZeroDivisionError("A square-root series requires a nonzero constant term.")
    out = zeros(order)
    out[0] = mp.sqrt(a[0])
    for m in range(1, order + 1):
        accumulator = a[m] if m < len(a) else mp.mpf("0")
        for k in range(1, m):
            accumulator -= out[k] * out[m - k]
        out[m] = accumulator / (2 * out[0])
    return out


def inverse_sqrt(a: Series, order: int) -> Series:
    return reciprocal(sqrt(a, order), order)
