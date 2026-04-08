from __future__ import annotations

from dataclasses import dataclass

import mpmath as mp

from .article_constants import PI, S_IR, S_UV, ScalarModelPreset
from .common import bisect_root, bracket_root
from .series import add, multiply, one, reciprocal, scale, sqrt, sub


@dataclass(frozen=True)
class ArticleScalarResult:
    """Evaluated article scalar branch D_C(α)."""

    alpha: mp.mpf
    x: mp.mpf
    model: ScalarModelPreset
    d_c: mp.mpf
    theta1_effective: mp.mpf
    phi_dyn: mp.mpf
    radicand: mp.mpf
    residual: mp.mpf



def dynamic_numerator(d_c: mp.mpf, model: ScalarModelPreset) -> mp.mpf:
    """Eq. (app-ND-def): N(D)."""
    d_c = mp.mpf(d_c)
    return (
        model.a_uv * (1 + S_IR * d_c)
        + model.a_ir * (1 + S_UV * d_c)
        - 2 * model.chi * d_c * mp.sqrt(model.a_uv * model.a_ir)
    )



def dynamic_denominator(d_c: mp.mpf, model: ScalarModelPreset) -> mp.mpf:
    """Eq. (app-QD-def): Q(D)."""
    d_c = mp.mpf(d_c)
    return (1 + S_UV * d_c) * (1 + S_IR * d_c) - model.chi**2 * d_c**2



def dynamic_numerator_prime(model: ScalarModelPreset) -> mp.mpf:
    """Eq. (app-Nprime-def): N'(D)."""
    return (
        model.a_uv * S_IR
        + model.a_ir * S_UV
        - 2 * model.chi * mp.sqrt(model.a_uv * model.a_ir)
    )



def dynamic_denominator_prime(d_c: mp.mpf, model: ScalarModelPreset) -> mp.mpf:
    """Eq. (app-Qprime-def): Q'(D)."""
    d_c = mp.mpf(d_c)
    return S_UV + S_IR + 2 * d_c * (S_UV * S_IR - model.chi**2)



def phi_dyn(d_c: mp.mpf, model: ScalarModelPreset) -> mp.mpf:
    """Visible dynamic response Φ_dyn(D) in the realized scalar truncation."""
    return (
        model.source_renorm_factor
        * dynamic_numerator(d_c, model)
        / dynamic_denominator(d_c, model)
    )



def phi_dyn_prime(d_c: mp.mpf, model: ScalarModelPreset) -> mp.mpf:
    """Derivative Φ'_dyn(D) entering the appendix sensitivity identities."""
    d_c = mp.mpf(d_c)
    numerator = dynamic_numerator(d_c, model)
    denominator = dynamic_denominator(d_c, model)
    return model.source_renorm_factor * (
        dynamic_numerator_prime(model) * denominator
        - numerator * dynamic_denominator_prime(d_c, model)
    ) / denominator**2



def radicand(d_c: mp.mpf, model: ScalarModelPreset) -> mp.mpf:
    """Scalar mother radicand R_moth(D) = 1 - Θ_eff D + D² Φ_dyn(D)."""
    d_c = mp.mpf(d_c)
    return 1 - model.theta1_effective * d_c + d_c**2 * phi_dyn(d_c, model)



def radicand_prime(d_c: mp.mpf, model: ScalarModelPreset) -> mp.mpf:
    """Derivative of the scalar mother radicand with respect to D."""
    d_c = mp.mpf(d_c)
    return -model.theta1_effective + 2 * d_c * phi_dyn(d_c, model) + d_c**2 * phi_dyn_prime(d_c, model)



def dc_equation(d_c: mp.mpf, alpha: mp.mpf, model: ScalarModelPreset) -> mp.mpf:
    """Scalar closure written as D - x sqrt(R_moth(D)) = 0, x = α / π."""
    x_value = mp.mpf(alpha) / PI
    return d_c - x_value * mp.sqrt(radicand(d_c, model))



def solve_dc_article(alpha: mp.mpf, model: ScalarModelPreset) -> ArticleScalarResult:
    """Solve the article scalar relation D_C(α) on the physical positive branch."""
    alpha = mp.mpf(alpha)
    x_value = alpha / PI
    equation = lambda d_value: dc_equation(d_value, alpha, model)
    left = mp.mpf("0")
    right = max(mp.mpf("0.1"), 4 * x_value)
    left, right = bracket_root(equation, left, right)
    d_c = bisect_root(equation, left, right)
    rad = radicand(d_c, model)
    residual = dc_equation(d_c, alpha, model)
    return ArticleScalarResult(
        alpha=alpha,
        x=x_value,
        model=model,
        d_c=d_c,
        theta1_effective=model.theta1_effective,
        phi_dyn=phi_dyn(d_c, model),
        radicand=rad,
        residual=residual,
    )



def alpha_from_dc_article(d_c: mp.mpf, model: ScalarModelPreset) -> mp.mpf:
    """Algebraic inversion α(D) = π D / sqrt(R_moth(D))."""
    d_c = mp.mpf(d_c)
    return PI * d_c / mp.sqrt(radicand(d_c, model))



def alpha_prime_from_d(d_c: mp.mpf, model: ScalarModelPreset) -> mp.mpf:
    """Closed derivative dα/dD used in the appendix sensitivity table."""
    d_c = mp.mpf(d_c)
    r_value = radicand(d_c, model)
    return PI * (1 - (d_c / 2) * radicand_prime(d_c, model) / r_value) / mp.sqrt(r_value)



def scalar_series_coefficients(model: ScalarModelPreset, order: int = 5) -> dict[int, mp.mpf]:
    """
    Formal coefficients of D_C(x) = sum_{n>=1} c_n x^n.

    The recursion follows directly from the implicit scalar mother law
        D(x) = x sqrt(R_moth(D(x))).
    Because R_moth(0) = 1, the coefficient c_n depends only on the already
    determined coefficients c_1, ..., c_{n-1}.
    """
    if order < 1:
        raise ValueError("Series order must be positive.")

    coefficients = [mp.mpf("0") for _ in range(order + 1)]
    coefficients[1] = mp.mpf("1")

    q1 = S_UV + S_IR
    q2 = S_UV * S_IR - model.chi**2
    n_const = model.a_uv + model.a_ir
    n_linear = dynamic_numerator_prime(model)

    for n in range(2, order + 1):
        truncation_order = n - 1
        d_series = coefficients[: truncation_order + 1]
        d_square = multiply(d_series, d_series, truncation_order)

        n_series = add(
            [n_const] + [mp.mpf("0")] * truncation_order,
            scale(d_series, n_linear, truncation_order),
            truncation_order,
        )
        q_series = add(
            one(truncation_order),
            scale(d_series, q1, truncation_order),
            truncation_order,
        )
        q_series = add(
            q_series,
            scale(d_square, q2, truncation_order),
            truncation_order,
        )
        inverse_q = reciprocal(q_series, truncation_order)

        r_series = add(
            one(truncation_order),
            scale(d_series, -model.theta1_effective, truncation_order),
            truncation_order,
        )
        dynamic_series = multiply(d_square, multiply(n_series, inverse_q, truncation_order), truncation_order)
        r_series = add(
            r_series,
            scale(dynamic_series, model.source_renorm_factor, truncation_order),
            truncation_order,
        )
        sqrt_r = sqrt(r_series, truncation_order)
        coefficients[n] = sqrt_r[n - 1]

    return {index: coefficients[index] for index in range(1, order + 1)}
