from __future__ import annotations

from dataclasses import dataclass

import mpmath as mp

from .article_constants import K_ALP, PI, ScalarModelPreset
from .scalar_article import scalar_series_coefficients, solve_dc_article
from .series import add, inverse_sqrt, multiply, one, reciprocal, scale, shift, sqrt, sub


@dataclass(frozen=True)
class G2PhysicalPointResult:
    """Direct physical-point evaluation of the Relator g-factor proxy."""

    alpha: mp.mpf
    x: mp.mpf
    d_c: mp.mpf
    zeta: mp.mpf
    lambda_pi: mp.mpf
    d_total: mp.mpf
    g_e2: mp.mpf
    a_e: mp.mpf



def zeta_from_dc(d_c: mp.mpf) -> mp.mpf:
    """ALP inversion ζ(D) = [sqrt(1 + 4D/(3π²)) - 1] / 2."""
    d_c = mp.mpf(d_c)
    return (mp.sqrt(1 + 4 * d_c / (3 * PI**2)) - 1) / 2



def lambda_pi_from_dc(d_c: mp.mpf) -> mp.mpf:
    """Auxiliary shell branch Λ_π(D) = 2π² ζ(D) / K."""
    return 2 * PI**2 * zeta_from_dc(d_c) / K_ALP



def d_total_from_x_and_dc(x_value: mp.mpf, d_c: mp.mpf) -> mp.mpf:
    """Operational orbital slowdown map used in the article's g-2 cross-check."""
    x_value = mp.mpf(x_value)
    d_c = mp.mpf(d_c)
    zeta_value = zeta_from_dc(d_c)
    return d_c - x_value * zeta_value + x_value**2 * zeta_value**2 / (4 * d_c)



def g_e2_from_d_total(d_total: mp.mpf) -> mp.mpf:
    """g_e / 2 = (1 - D_total)^(-1/2)."""
    d_total = mp.mpf(d_total)
    return 1 / mp.sqrt(1 - d_total)



def a_e_from_d_total(d_total: mp.mpf) -> mp.mpf:
    """Anomalous part a_e = g_e/2 - 1."""
    return g_e2_from_d_total(d_total) - 1



def evaluate_g2_at_alpha(alpha: mp.mpf, model: ScalarModelPreset) -> G2PhysicalPointResult:
    """Evaluate the full Relator g-2 cross-check chain at a direct α input."""
    alpha = mp.mpf(alpha)
    x_value = alpha / PI
    scalar = solve_dc_article(alpha, model)
    d_total = d_total_from_x_and_dc(x_value, scalar.d_c)
    g_e2 = g_e2_from_d_total(d_total)
    return G2PhysicalPointResult(
        alpha=alpha,
        x=x_value,
        d_c=scalar.d_c,
        zeta=zeta_from_dc(scalar.d_c),
        lambda_pi=lambda_pi_from_dc(scalar.d_c),
        d_total=d_total,
        g_e2=g_e2,
        a_e=g_e2 - 1,
    )



def relator_g2_coefficients(model: ScalarModelPreset, order: int = 10) -> dict[int, mp.mpf]:
    """
    Taylor coefficients of a_e^Rel(x) = Σ Â_1^(2n) x^n.

    The nontrivial coefficients of g_e/2 coincide with those of a_e.
    """
    if order < 1:
        raise ValueError("The requested order must be positive.")

    d_coefficients = scalar_series_coefficients(model, order)
    d_series = [mp.mpf("0")] * (order + 1)
    for n in range(1, order + 1):
        d_series[n] = d_coefficients[n]

    zeta_series = scale(
        sub(
            sqrt(add(one(order), scale(d_series, 4 / (3 * PI**2), order), order), order),
            one(order),
            order,
        ),
        mp.mpf("0.5"),
        order,
    )

    y_series = [mp.mpf("0")] * (order + 1)
    y_series[0] = mp.mpf("1")
    for power in range(1, order + 1):
        if power + 1 <= order:
            y_series[power] = d_series[power + 1]
    inverse_y = reciprocal(y_series, order)

    d_total_series = sub(d_series, shift(zeta_series, 1, order), order)
    third_term = shift(
        scale(multiply(multiply(zeta_series, zeta_series, order), inverse_y, order), mp.mpf("0.25"), order),
        1,
        order,
    )
    d_total_series = add(d_total_series, third_term, order)

    g_series = inverse_sqrt(sub(one(order), d_total_series, order), order)
    a_series = sub(g_series, one(order), order)
    return {n: a_series[n] for n in range(1, order + 1)}
