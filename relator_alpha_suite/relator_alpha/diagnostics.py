from __future__ import annotations

from dataclasses import dataclass

import mpmath as mp

from .alp import dc_lock_from_lambda
from .article_constants import K_ALP, LAMBDA_GEOM_MEAN, PI, REALIZED_MODELS, S_IR, S_UV, ScalarModelPreset
from .common import configure_precision
from .scalar_article import alpha_from_dc_article, alpha_prime_from_d, radicand


@dataclass(frozen=True)
class ExactSensitivityReport:
    d_lock: mp.mpf
    radicand_at_lock: mp.mpf
    d_dlock_d_lambda: mp.mpf
    d_alpha_d_d: mp.mpf
    d_alpha_d_lambda: mp.mpf
    d_alpha_d_k: mp.mpf
    dlog_alpha_dlog_lambda: mp.mpf
    dlog_alpha_dlog_k: mp.mpf
    d_alpha_d_theta_eff: mp.mpf
    dlog_alpha_dlog_theta_eff: mp.mpf
    dlog_alpha_dlog_b11: mp.mpf


@dataclass(frozen=True)
class NumericalSensitivityReport:
    theta_eff: mp.mpf
    a_uv: mp.mpf
    a_ir: mp.mpf
    chi: mp.mpf
    s_uv: mp.mpf
    s_ir: mp.mpf



def realized_alpha_prediction(model: ScalarModelPreset, lambda_geom: mp.mpf = LAMBDA_GEOM_MEAN) -> mp.mpf:
    """Direct geometric-lock prediction α_pred = α(D_lock(Λ_geom))."""
    d_lock = dc_lock_from_lambda(lambda_geom)
    return alpha_from_dc_article(d_lock, model)



def exact_sensitivities(
    model: ScalarModelPreset,
    lambda_geom: mp.mpf = LAMBDA_GEOM_MEAN,
) -> ExactSensitivityReport:
    """Closed appendix sensitivities evaluated at the realized lock point."""
    lambda_geom = mp.mpf(lambda_geom)
    d_lock = dc_lock_from_lambda(lambda_geom)
    alpha_pred = alpha_from_dc_article(d_lock, model)
    r_star = radicand(d_lock, model)

    d_dlock_d_lambda = mp.mpf("1.5") * K_ALP * (1 + K_ALP * lambda_geom / PI**2)
    d_dlock_d_k = mp.mpf("1.5") * lambda_geom * (1 + K_ALP * lambda_geom / PI**2)
    d_alpha_d_d = alpha_prime_from_d(d_lock, model)
    d_alpha_d_lambda = d_alpha_d_d * d_dlock_d_lambda
    d_alpha_d_k = d_alpha_d_d * d_dlock_d_k
    d_alpha_d_theta_eff = alpha_pred * d_lock / (2 * r_star)
    return ExactSensitivityReport(
        d_lock=d_lock,
        radicand_at_lock=r_star,
        d_dlock_d_lambda=d_dlock_d_lambda,
        d_alpha_d_d=d_alpha_d_d,
        d_alpha_d_lambda=d_alpha_d_lambda,
        d_alpha_d_k=d_alpha_d_k,
        dlog_alpha_dlog_lambda=lambda_geom * d_alpha_d_lambda / alpha_pred,
        dlog_alpha_dlog_k=K_ALP * d_alpha_d_k / alpha_pred,
        d_alpha_d_theta_eff=d_alpha_d_theta_eff,
        dlog_alpha_dlog_theta_eff=model.theta1_effective * d_alpha_d_theta_eff / alpha_pred,
        dlog_alpha_dlog_b11=model.b11_effective * d_alpha_d_theta_eff / alpha_pred,
    )



def appendix_source_prefactor(
    chi: mp.mpf,
    s_uv: mp.mpf,
    s_ir: mp.mpf,
) -> mp.mpf:
    """
    Source prefactor used by the appendix realized-lock convention.

    The appendix locally writes λ_u = (1 - ρ_dyn²)^(-1/2), so λ_u² becomes
    1 / (1 - ρ_dyn²).  This is numerically the same as the main-text λ_u.
    """
    chi = mp.mpf(chi)
    s_uv = mp.mpf(s_uv)
    s_ir = mp.mpf(s_ir)
    rho = chi / mp.sqrt(s_uv * s_ir)
    return 1 / (1 - rho**2)



def realized_alpha_from_parameters(
    *,
    d_lock: mp.mpf,
    theta_eff: mp.mpf,
    a_uv: mp.mpf,
    a_ir: mp.mpf,
    chi: mp.mpf,
    s_uv: mp.mpf,
    s_ir: mp.mpf,
) -> mp.mpf:
    """Direct α(D_lock) evaluation for numerical sensitivity studies."""
    d_lock = mp.mpf(d_lock)
    theta_eff = mp.mpf(theta_eff)
    a_uv = mp.mpf(a_uv)
    a_ir = mp.mpf(a_ir)
    chi = mp.mpf(chi)
    s_uv = mp.mpf(s_uv)
    s_ir = mp.mpf(s_ir)
    prefactor = appendix_source_prefactor(chi, s_uv, s_ir)
    numerator = a_uv * (1 + s_ir * d_lock) + a_ir * (1 + s_uv * d_lock) - 2 * chi * d_lock * mp.sqrt(a_uv * a_ir)
    denominator = (1 + s_uv * d_lock) * (1 + s_ir * d_lock) - chi**2 * d_lock**2
    radicand_value = 1 - theta_eff * d_lock + prefactor * d_lock**2 * numerator / denominator
    return PI * d_lock / mp.sqrt(radicand_value)



def symmetric_log_sensitivity(
    model: ScalarModelPreset,
    parameter_name: str,
    *,
    lambda_geom: mp.mpf = LAMBDA_GEOM_MEAN,
    epsilon: mp.mpf = mp.mpf("1e-8"),
) -> mp.mpf:
    """Symmetric logarithmic sensitivity estimator used in the appendix table."""
    if model.family != "realized_lock":
        raise ValueError("Numerical sensitivities are defined for the realized appendix convention.")

    d_lock = dc_lock_from_lambda(lambda_geom)
    base = {
        "theta_eff": model.theta1_effective,
        "a_uv": model.a_uv,
        "a_ir": model.a_ir,
        "chi": model.chi,
        "s_uv": S_UV,
        "s_ir": S_IR,
    }
    if parameter_name not in base:
        raise KeyError(f"Unknown sensitivity parameter: {parameter_name}")

    plus = dict(base)
    minus = dict(base)
    plus[parameter_name] = base[parameter_name] * (1 + epsilon)
    minus[parameter_name] = base[parameter_name] * (1 - epsilon)

    alpha_plus = realized_alpha_from_parameters(d_lock=d_lock, **plus)
    alpha_minus = realized_alpha_from_parameters(d_lock=d_lock, **minus)
    return (mp.log(alpha_plus) - mp.log(alpha_minus)) / (mp.log(1 + epsilon) - mp.log(1 - epsilon))



def numerical_sensitivities(
    model: ScalarModelPreset | None = None,
    *,
    lambda_geom: mp.mpf = LAMBDA_GEOM_MEAN,
    epsilon: mp.mpf = mp.mpf("1e-8"),
) -> NumericalSensitivityReport:
    """Compute the numerical appendix sensitivity audit for the realized rank-100 model."""
    if model is None:
        model = REALIZED_MODELS[100]
    return NumericalSensitivityReport(
        theta_eff=symmetric_log_sensitivity(model, "theta_eff", lambda_geom=lambda_geom, epsilon=epsilon),
        a_uv=symmetric_log_sensitivity(model, "a_uv", lambda_geom=lambda_geom, epsilon=epsilon),
        a_ir=symmetric_log_sensitivity(model, "a_ir", lambda_geom=lambda_geom, epsilon=epsilon),
        chi=symmetric_log_sensitivity(model, "chi", lambda_geom=lambda_geom, epsilon=epsilon),
        s_uv=symmetric_log_sensitivity(model, "s_uv", lambda_geom=lambda_geom, epsilon=epsilon),
        s_ir=symmetric_log_sensitivity(model, "s_ir", lambda_geom=lambda_geom, epsilon=epsilon),
    )



def rank_stability(
    ranks: tuple[int, ...] = (1, 3, 5, 100),
    *,
    lambda_geom: mp.mpf = LAMBDA_GEOM_MEAN,
) -> dict[int, tuple[mp.mpf, mp.mpf]]:
    """Return α_pred^(N) and drift from the N=100 realized-lock value."""
    alphas = {rank: realized_alpha_prediction(REALIZED_MODELS[rank], lambda_geom) for rank in ranks}
    reference = alphas[100]
    return {rank: (alphas[rank], alphas[rank] - reference) for rank in ranks}



def precision_audit(
    precisions: tuple[int, ...] = (50, 80, 120, 200, 300),
    *,
    lambda_geom: mp.mpf = LAMBDA_GEOM_MEAN,
) -> dict[int, tuple[mp.mpf, mp.mpf]]:
    """Return α_pred at several working precisions and drift from 300 dps."""
    original_dps = mp.mp.dps
    results: dict[int, mp.mpf] = {}
    for dps in precisions:
        configure_precision(dps)
        results[dps] = realized_alpha_prediction(REALIZED_MODELS[100], lambda_geom)
    reference = results[300]
    out = {dps: (alpha_value, alpha_value - reference) for dps, alpha_value in results.items()}
    configure_precision(original_dps)
    return out
