from __future__ import annotations

from dataclasses import dataclass

import mpmath as mp

from .article_constants import (
    C_UV_GAUSS,
    DELTA_LAMBDA_OUT_EXACT,
    DELTA_LAMBDA_UV_TO_IR_MEAN,
    EXTERIOR_MEMORY_CONTRIBUTION_MEAN,
    ETA_0,
    ELL_0,
    LAMBDA_GEOM_MEAN,
    LAMBDA_IND,
    P_IR_CHI_MEAN,
    R_CHI_MEAN,
)


@dataclass(frozen=True)
class LambdaGeomResult:
    """Container for the vector-channel geometric representative."""
    eta_0: mp.mpf
    ell_0: mp.mpf
    lambda_ind: mp.mpf
    c_uv_gauss: mp.mpf
    delta_lambda_out: mp.mpf
    p_ir_chi: mp.mpf
    delta_lambda_uv_to_ir: mp.mpf
    r_chi: mp.mpf
    exterior_memory_contribution: mp.mpf
    lambda_geom: mp.mpf
    status: str


def build_lambda_geom_representative(
    *,
    p_ir_chi: mp.mpf | None = None,
    r_chi: mp.mpf | None = None,
    delta_lambda_out: mp.mpf | None = None,
) -> LambdaGeomResult:
    """
    Build Λ_geom from the article's load-bearing decomposition.

    Important honesty note:
    the full operator-exact Λ_geom is *not* reduced to a one-line closed form
    in the manuscript. The default returned here is therefore the article's
    diagnostic representative Λ_geom^(mean), assembled from:
        Λ_ind + ΔΛ_UV→IR^(mean)(ℓ0) + ΔΛ_out(η0) R_χ^(mean).
    The function also accepts custom P_IR^(χ) and R_χ values.
    """
    delta_lambda_out = (
        DELTA_LAMBDA_OUT_EXACT if delta_lambda_out is None else mp.mpf(delta_lambda_out)
    )
    if p_ir_chi is None:
        p_ir_value = P_IR_CHI_MEAN
        delta_lambda_uv_to_ir = DELTA_LAMBDA_UV_TO_IR_MEAN
    else:
        p_ir_value = mp.mpf(p_ir_chi)
        delta_lambda_uv_to_ir = C_UV_GAUSS * p_ir_value
    r_chi_value = R_CHI_MEAN if r_chi is None else mp.mpf(r_chi)
    exterior_memory_contribution = delta_lambda_out * r_chi_value
    lambda_geom = LAMBDA_IND + delta_lambda_uv_to_ir + exterior_memory_contribution
    status = (
        "diagnostic representative Λ_geom^(mean)"
        if p_ir_chi is None and r_chi is None and delta_lambda_out == DELTA_LAMBDA_OUT_EXACT
        else "custom representative assembled from supplied shell inputs"
    )
    return LambdaGeomResult(
        eta_0=ETA_0,
        ell_0=ELL_0,
        lambda_ind=LAMBDA_IND,
        c_uv_gauss=C_UV_GAUSS,
        delta_lambda_out=delta_lambda_out,
        p_ir_chi=p_ir_value,
        delta_lambda_uv_to_ir=delta_lambda_uv_to_ir,
        r_chi=r_chi_value,
        exterior_memory_contribution=exterior_memory_contribution,
        lambda_geom=lambda_geom,
        status=status,
    )
