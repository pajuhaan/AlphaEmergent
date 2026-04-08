from __future__ import annotations

from dataclasses import dataclass
from typing import Dict

import mpmath as mp

from .article_constants import PI, QED_A1_COEFFICIENTS


@dataclass(frozen=True)
class QEDInducedCoefficients:
    """
    QED-induced scalar coefficients reconstructed from the universal A_1^(2n).
    """
    a1_coefficients: Dict[int, mp.mpf]
    d_coefficients: Dict[int, mp.mpf]
    c_coefficients: Dict[int, mp.mpf]
    order: int


def d_coefficients_from_qed_a1(a_coefficients: Dict[int, mp.mpf]) -> Dict[int, mp.mpf]:
    """
    Eqs. (d12-qed)–(d5-qed).
    """
    a1 = a_coefficients[1]
    a2 = a_coefficients[2]
    a3 = a_coefficients[3]
    a4 = a_coefficients[4]
    a5 = a_coefficients[5]
    return {
        1: 2 * a1,
        2: 2 * a2 - 3 * a1**2,
        3: 2 * a3 - 6 * a1 * a2 + 4 * a1**3,
        4: 2 * a4 - 6 * a1 * a3 - 3 * a2**2 + 12 * a1**2 * a2 - 5 * a1**4,
        5: (
            2 * a5
            - 6 * a1 * a4
            - 6 * a2 * a3
            + 12 * a1**2 * a3
            + 12 * a1 * a2**2
            - 20 * a1**3 * a2
            + 6 * a1**5
        ),
    }


def c_coefficients_from_qed_d(d_coefficients: Dict[int, mp.mpf]) -> Dict[int, mp.mpf]:
    """
    Eqs. (c12-recursion)–(c5-recursion).
    """
    d1 = d_coefficients[1]
    d2 = d_coefficients[2]
    d3 = d_coefficients[3]
    d4 = d_coefficients[4]
    d5 = d_coefficients[5]
    return {
        1: d1,
        2: d2 + 1 / (3 * PI**2),
        3: d3 + d2 / (3 * PI**2) - 1 / (36 * PI**4),
        4: d4 + d3 / (3 * PI**2) - 5 * d2 / (36 * PI**4),
        5: (
            d5
            + d4 / (3 * PI**2)
            - 5 * d3 / (36 * PI**4)
            - d2**2 / (9 * PI**4)
            + d2 / (18 * PI**6)
            + 5 / (1296 * PI**8)
        ),
    }


def build_qed_induced_coefficients(order: int = 5) -> QEDInducedCoefficients:
    """
    Reconstruct the QED-induced scalar branch retained in the manuscript.

    The present implementation keeps exactly the first five universal
    pure-photonic coefficients used in the article.
    """
    if order < 1 or order > 5:
        raise ValueError("The current implementation supports 1 <= order <= 5.")
    d_coefficients = d_coefficients_from_qed_a1(QED_A1_COEFFICIENTS)
    c_coefficients = c_coefficients_from_qed_d(d_coefficients)
    return QEDInducedCoefficients(
        a1_coefficients=QED_A1_COEFFICIENTS,
        d_coefficients=d_coefficients,
        c_coefficients=c_coefficients,
        order=order,
    )


def dc_qed_truncated(x_value: mp.mpf, coefficients: QEDInducedCoefficients) -> mp.mpf:
    """
    Truncated QED-induced scalar branch
        D_C^QED,[M](x) = sum_{n=1}^M c_n^QED x^n.
    """
    x_value = mp.mpf(x_value)
    total = mp.mpf("0")
    for n in range(1, coefficients.order + 1):
        total += coefficients.c_coefficients[n] * x_value**n
    return total


def dc_qed_from_alpha(alpha: mp.mpf, coefficients: QEDInducedCoefficients) -> mp.mpf:
    """
    Evaluate the QED-induced scalar branch at a direct α input.
    """
    return dc_qed_truncated(mp.mpf(alpha) / PI, coefficients)
