from __future__ import annotations

from dataclasses import dataclass

import mpmath as mp

from .article_constants import K_ALP, PI, ScalarModelPreset
from .common import bisect_root, bracket_root
from .scalar_article import alpha_from_dc_article, solve_dc_article
from .scalar_qed import QEDInducedCoefficients, dc_qed_truncated


@dataclass(frozen=True)
class AlphaLockResult:
    """Generic ALP solve result."""
    branch_label: str
    lambda_geom: mp.mpf
    zeta: mp.mpf
    d_lock: mp.mpf
    alpha: mp.mpf
    inverse_alpha: mp.mpf
    residual: mp.mpf
    note: str


def dc_lock_from_lambda(lambda_geom: mp.mpf) -> mp.mpf:
    """
    Eq. (Dlock-of-Lambda):
        D_C^lock(Λ) = (3/2) K Λ [1 + K Λ / (2π²)].
    """
    lambda_geom = mp.mpf(lambda_geom)
    return mp.mpf("1.5") * K_ALP * lambda_geom * (1 + K_ALP * lambda_geom / (2 * PI**2))


def zeta_from_lambda(lambda_geom: mp.mpf) -> mp.mpf:
    """
    ζ = K Λ / (2π²) in the baseline shell-lock closure.
    """
    lambda_geom = mp.mpf(lambda_geom)
    return K_ALP * lambda_geom / (2 * PI**2)


def zeta_from_dc(d_c: mp.mpf) -> mp.mpf:
    """
    Auxiliary inverse-ALP relation used in the pure-photonic cross-check:
        ζ(D) = [sqrt(1 + 4D/(3π²)) - 1] / 2.
    """
    d_c = mp.mpf(d_c)
    return (mp.sqrt(1 + 4 * d_c / (3 * PI**2)) - 1) / 2


def d_total_from_d_and_alpha(d_c: mp.mpf, alpha: mp.mpf) -> mp.mpf:
    """
    Relator bridge:
        D_total = D - x ζ(D) + x² ζ(D)² / (4D), x = α/π.
    """
    d_c = mp.mpf(d_c)
    alpha = mp.mpf(alpha)
    x_value = alpha / PI
    zeta_value = zeta_from_dc(d_c)
    return d_c - x_value * zeta_value + x_value**2 * zeta_value**2 / (4 * d_c)


def solve_alpha_from_article_branch(
    lambda_geom: mp.mpf,
    model: ScalarModelPreset,
) -> AlphaLockResult:
    """
    Solve α from the article scalar branch by direct algebraic inversion.

    Once Λ_geom is fixed, no coupled fixed-point solve remains:
        α = π D_lock / sqrt(R_moth(D_lock)).
    """
    lambda_geom = mp.mpf(lambda_geom)
    d_lock = dc_lock_from_lambda(lambda_geom)
    alpha = alpha_from_dc_article(d_lock, model)
    direct_check = solve_dc_article(alpha, model)
    residual = direct_check.d_c - d_lock
    return AlphaLockResult(
        branch_label=f"article:{model.key}",
        lambda_geom=lambda_geom,
        zeta=zeta_from_lambda(lambda_geom),
        d_lock=d_lock,
        alpha=alpha,
        inverse_alpha=1 / alpha,
        residual=residual,
        note="direct algebraic inversion of the article scalar mother law",
    )


def solve_alpha_from_qed_branch(
    lambda_geom: mp.mpf,
    coefficients: QEDInducedCoefficients,
) -> AlphaLockResult:
    """
    Solve α from the truncated universal-QED-induced scalar branch.

    This is the one-way cross-check defined by
        D_C^QED,[M](α/π) = D_C^lock(Λ_geom).
    """
    lambda_geom = mp.mpf(lambda_geom)
    d_lock = dc_lock_from_lambda(lambda_geom)

    function = lambda x_value: dc_qed_truncated(x_value, coefficients) - d_lock
    left = mp.mpf("0")
    right = mp.mpf("0.01")
    left, right = bracket_root(function, left, right)
    x_root = bisect_root(function, left, right)
    alpha = PI * x_root
    residual = dc_qed_truncated(x_root, coefficients) - d_lock
    return AlphaLockResult(
        branch_label=f"qed:M={coefficients.order}",
        lambda_geom=lambda_geom,
        zeta=zeta_from_lambda(lambda_geom),
        d_lock=d_lock,
        alpha=alpha,
        inverse_alpha=1 / alpha,
        residual=residual,
        note="one-way universal-QED cross-check over the same ALP target",
    )
