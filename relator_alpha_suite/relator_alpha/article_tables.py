
from __future__ import annotations

from functools import lru_cache

import mpmath as mp

from .alp import (
    d_total_from_d_and_alpha,
    dc_lock_from_lambda,
    solve_alpha_from_article_branch,
    solve_alpha_from_qed_branch,
    zeta_from_lambda,
)
from .article_constants import (
    ALPHA_BENCHMARK_DIRECT,
    ALPHA_BENCHMARK_DIRECT_U,
    ALPHA_BENCHMARK_INVERSE,
    ALPHA_BENCHMARK_INVERSE_U,
    ALPHA_CODATA22_DIRECT,
    B11_3RHO,
    B11_GLUE,
    C_UNI,
    C_UV_GAUSS,
    CURRENT_MODELS,
    DELTA_LAMBDA_OUT_EXACT,
    DELTA_LAMBDA_UV_TO_IR_MEAN,
    ELL_0,
    ETA_0,
    EXPERIMENTAL_ALPHA_REFERENCES,
    EXTERIOR_MEMORY_CONTRIBUTION_MEAN,
    GAMMA_CL_0,
    K_ALP,
    LAMBDA_GEOM_MEAN,
    LAMBDA_H,
    LAMBDA_IND,
    N_C,
    PI,
    P_IR_CHI_MEAN,
    QED_A1_COEFFICIENTS,
    QED_A_E_BENCHMARK_XSTAR,
    QED_G_E2_BENCHMARK_XSTAR,
    RHO_STAT,
    R_CHI_MEAN,
    REALIZED_MODELS,
    REFINED_MODELS,
    S_IR,
    S_UV,
    SCALAR_RANK_DATA,
    THETA_1_3RHO,
    THETA_1_COLL,
    THETA_1_CORE,
    THETA_1_LOG,
    VECTOR_SHELL_SOURCE_WEIGHTS,
    compute_b11_3rho_from_formula,
    compute_b11_glue_from_formula,
)
from .common import format_mpf, format_percent
from .diagnostics import exact_sensitivities, numerical_sensitivities, precision_audit, rank_stability
from .g2 import evaluate_g2_at_alpha, relator_g2_coefficients
from .presentation import ArticleTable, TableColumn
from .scalar_article import (
    dynamic_denominator,
    dynamic_numerator,
    solve_dc_article,
    scalar_series_coefficients,
)
from .scalar_qed import build_qed_induced_coefficients, dc_qed_from_alpha


RANKS = (1, 3, 5, 100)



def _num(value: mp.mpf | float | int, digits: int) -> str:
    return format_mpf(value, digits)



def _signed_num(value: mp.mpf | float | int, digits: int) -> str:
    x = mp.mpf(value)
    return ("+" if x > 0 else "") + _num(x, digits)



def _ppt(value: mp.mpf, reference: mp.mpf, digits: int = 8) -> str:
    return f"{float(mp.mpf('1e12') * value / reference):.{digits}f} ppt"


@lru_cache(maxsize=None)
def _current_lock_result():
    return solve_alpha_from_article_branch(LAMBDA_GEOM_MEAN, CURRENT_MODELS[5])


@lru_cache(maxsize=None)
def _qed_lock_result():
    coefficients = build_qed_induced_coefficients(order=5)
    return solve_alpha_from_qed_branch(LAMBDA_GEOM_MEAN, coefficients)


@lru_cache(maxsize=None)
def _qed_coefficients(order: int = 5):
    return build_qed_induced_coefficients(order=order)


@lru_cache(maxsize=None)
def _scalar_coefficients(model_key: str, order: int) -> dict[int, mp.mpf]:
    model_lookup = {
        **{m.key: m for m in CURRENT_MODELS.values()},
        **{m.key: m for m in REFINED_MODELS.values()},
        **{m.key: m for m in REALIZED_MODELS.values()},
    }
    return scalar_series_coefficients(model_lookup[model_key], order)


@lru_cache(maxsize=None)
def _g2_coefficients(model_key: str, order: int) -> dict[int, mp.mpf]:
    model_lookup = {
        **{m.key: m for m in CURRENT_MODELS.values()},
        **{m.key: m for m in REFINED_MODELS.values()},
        **{m.key: m for m in REALIZED_MODELS.values()},
    }
    return relator_g2_coefficients(model_lookup[model_key], order)



def scalar_static_fixed_data_table(digits: int) -> ArticleTable:
    reconstructed_glue = compute_b11_glue_from_formula()
    reconstructed_3rho = compute_b11_3rho_from_formula()
    rows = [
        ["η₀", "1/π", _num(ETA_0, digits)],
        ["ℓ₀", "1/(π√π)", _num(ELL_0, digits)],
        ["C_uni", "(1/π)(4/3 + 1/(4π²))", _num(C_UNI, digits)],
        ["S_UV", "log 2", _num(S_UV, digits)],
        ["S_IR", "1/(8π²)", _num(S_IR, digits)],
        ["K", "(150π² - 8π⁴ - 315)/(180π⁶)", _num(K_ALP, digits)],
        ["γ_cℓ(0)", "static collar–log overlap", _num(GAMMA_CL_0, digits)],
        ["N_c", "√(π/8) ℓ₀", _num(N_C, digits)],
        ["ρ_stat", "γ_cℓ(0)/√(N_c N_ℓ)", _num(RHO_STAT, digits)],
        ["Λ_h", "(1-ρ_stat²)⁻¹", _num(LAMBDA_H, digits)],
        ["Θ₁,core", "2π C_uni", _num(THETA_1_CORE, digits)],
        ["Θ₁^(coll)", "(γ_E/π²) C₁₁(ℓ₀)", _num(THETA_1_COLL, digits)],
        ["Θ₁^(log)", "η₀²/[8(1-η₀²/4)]", _num(THETA_1_LOG, digits)],
        ["B₁₁^glue", "visible 2-channel static load", _num(B11_GLUE, digits)],
        ["B₁₁^glue[formula]", "t_visᵀ (K_vis^(2))⁻¹ t_vis", _num(reconstructed_glue, digits)],
        ["ΔB₁₁^glue", "reconstructed - stored", _num(reconstructed_glue - B11_GLUE, digits)],
        ["B₁₁^(3,ρ)", "3-channel static completion", _num(B11_3RHO, digits)],
        ["B₁₁^(3,ρ)[formula]", "(t⊥)ᵀ (K_stat^(3,ρ))⁻¹ t⊥", _num(reconstructed_3rho, digits)],
        ["ΔB₁₁^(3,ρ)", "reconstructed - stored", _num(reconstructed_3rho - B11_3RHO, digits)],
        ["Θ₁^(3,ρ)", "Θ₁,core + B₁₁^(3,ρ)", _num(THETA_1_3RHO, digits)],
    ]
    return ArticleTable(
        label="tab:dc-static-fixed-data",
        title="Fixed shell geometry, ALP constant K, and static hidden-mode data",
        columns=[TableColumn("Symbol"), TableColumn("Def."), TableColumn("Value", "number")],
        rows=rows,
    )




def scalar_source_prefactor_convention_table(digits: int) -> ArticleTable:
    rows = []
    for rank in RANKS:
        rows.append(
            [
                str(rank),
                _num(CURRENT_MODELS[rank].lambda_u_display, digits),
                _num(CURRENT_MODELS[rank].source_renorm_factor, digits),
                _num(REFINED_MODELS[rank].source_renorm_factor, digits),
                _num(REALIZED_MODELS[rank].source_renorm_factor, digits),
            ]
        )
    return ArticleTable(
        label="tab:dc-source-prefactor-conventions",
        title="Rank-dependent source-prefactor conventions used by the current, refined, and realized scalar evaluators",
        columns=[
            TableColumn("N", "number"),
            TableColumn("λ_u^(main text)", "number"),
            TableColumn("current factor", "number"),
            TableColumn("refined factor = (λ_u^(main text))²", "number"),
            TableColumn("realized factor = λ_u^(main text)", "number"),
        ],
        rows=rows,
    )


def scalar_visible_dynamic_amplitudes_table(digits: int) -> ArticleTable:
    rows = []
    for rank in RANKS:
        data = SCALAR_RANK_DATA[rank]
        rows.append(
            [
                str(rank),
                _num(data.delta_a_mix, digits),
                _num(data.delta_a0, digits),
                _num(data.chi_cross, digits),
            ]
        )
    return ArticleTable(
        label="tab:dc-visible-dynamic-amplitudes",
        title="Rank-dependent visible dynamic diagnostics",
        columns=[
            TableColumn("N", "number"),
            TableColumn("ΔA_mix^(N)", "number"),
            TableColumn("δA₀^(N)", "number"),
            TableColumn("χ_cross^(N)", "number"),
        ],
        rows=rows,
    )



def scalar_visible_dynamic_block_table(digits: int) -> ArticleTable:
    rows = []
    for rank in RANKS:
        data = SCALAR_RANK_DATA[rank]
        rows.append(
            [
                str(rank),
                _num(data.a_uv, digits),
                _num(data.a_ir, digits),
                _num(data.rho_dyn, digits),
                _num(data.lambda_u_display, digits),
            ]
        )
    return ArticleTable(
        label="tab:dc-visible-dynamic-block",
        title="Visible dynamic two-channel parameters",
        columns=[
            TableColumn("N", "number"),
            TableColumn("A_uv^(N)", "number"),
            TableColumn("A_ir^(N)", "number"),
            TableColumn("ρ_dyn^(N)", "number"),
            TableColumn("λ_u^(N)", "number"),
        ],
        rows=rows,
    )





def scalar_current_physical_intermediates_table(digits: int) -> ArticleTable:
    rows = []
    for rank in RANKS:
        model = CURRENT_MODELS[rank]
        scalar = solve_dc_article(ALPHA_CODATA22_DIRECT, model)
        numerator = dynamic_numerator(scalar.d_c, model)
        denominator = dynamic_denominator(scalar.d_c, model)
        rows.append(
            [
                str(rank),
                _num(scalar.d_c, digits),
                _num(numerator, digits),
                _num(denominator, digits),
                _num(model.source_renorm_factor, digits),
                _num(scalar.phi_dyn, digits),
                _num(scalar.radicand, digits),
                _num(mp.sqrt(scalar.radicand), digits),
            ]
        )
    return ArticleTable(
        label="tab:dc-current-physical-intermediates",
        title="Current visible scalar intermediates at the physical reference point",
        columns=[
            TableColumn("N", "number"),
            TableColumn("D_C(x_ref)", "number"),
            TableColumn("N(D_C)", "number"),
            TableColumn("Q(D_C)", "number"),
            TableColumn("source factor", "number"),
            TableColumn("Φ_dyn(D_C)", "number"),
            TableColumn("R_moth(D_C)", "number"),
            TableColumn("√R_moth(D_C)", "number"),
        ],
        rows=rows,
        note="All quantities are evaluated on the current visible scalar branch with α = α_ref (CODATA-like reference point).",
    )


def scalar_hidden_dynamic_jacobi_table(digits: int) -> ArticleTable:
    data = SCALAR_RANK_DATA[100]
    rows = [
        ["∥v^(100)∥²", "hidden dynamic source norm", _num(data.hidden_source_norm, digits)],
        ["α₀^(100)", "first Jacobi diagonal coefficient", _num(data.jacobi_alpha0, digits)],
        ["β₁^(100)", "first Jacobi off-diagonal coefficient", _num(data.jacobi_beta1, digits)],
        ["α₁^(100)", "second Jacobi diagonal coefficient", _num(data.jacobi_alpha1, digits)],
        ["β₂^(100)", "second Jacobi off-diagonal coefficient", _num(data.jacobi_beta2, digits)],
        ["α₂^(100)", "third Jacobi diagonal coefficient", _num(data.jacobi_alpha2, digits)],
    ]
    return ArticleTable(
        label="tab:dc-hidden-dynamic-jacobi",
        title="First coefficients of the finite hidden dynamic Jacobi chain (N = 100)",
        columns=[TableColumn("Symbol"), TableColumn("Def."), TableColumn("Value", "number")],
        rows=rows,
    )



def scalar_physical_roots_table(digits: int) -> ArticleTable:
    rows = []
    for rank in RANKS:
        visible = solve_dc_article(ALPHA_CODATA22_DIRECT, CURRENT_MODELS[rank]).d_c
        refined = solve_dc_article(ALPHA_CODATA22_DIRECT, REFINED_MODELS[rank]).d_c
        rows.append(
            [
                str(rank),
                _num(visible, digits),
                _num(refined, digits),
                _num(refined - visible, digits),
            ]
        )
    return ArticleTable(
        label="tab:dc-physical-roots",
        title="Scalar roots at the physical reference point",
        columns=[
            TableColumn("N", "number"),
            TableColumn("D_C^(vis,N)(x_ref)", "number"),
            TableColumn("D_C^(ref,N)(x_ref)", "number"),
            TableColumn("ΔD_C^(N)", "number"),
        ],
        rows=rows,
    )



def scalar_coefficients_table(digits: int) -> ArticleTable:
    current_rank1 = _scalar_coefficients(CURRENT_MODELS[1].key, 5)
    current_rank5 = _scalar_coefficients(CURRENT_MODELS[5].key, 5)
    current_rank100 = _scalar_coefficients(CURRENT_MODELS[100].key, 5)
    # The published coefficient table follows the appendix / series convention,
    # numerically represented by the realized-lock model at N = 100.
    refined_series_rank100 = _scalar_coefficients(REALIZED_MODELS[100].key, 5)

    rows = []
    for n in range(1, 6):
        rows.append(
            [
                f"c_{n}",
                _num(current_rank1[n], digits),
                _num(current_rank5[n], digits),
                _num(current_rank100[n], digits),
                _num(refined_series_rank100[n], digits),
            ]
        )
    return ArticleTable(
        label="tab:dc-scalar-coefficients",
        title="First coefficients of the scalar branch D_C(x)",
        columns=[
            TableColumn("Coefficient"),
            TableColumn("N = 1 (vis)", "number"),
            TableColumn("N = 5 (vis)", "number"),
            TableColumn("N = 100 (vis)", "number"),
            TableColumn("N = 100 (ref)", "number"),
        ],
        rows=rows,
    )




def vector_exact_status_table(digits: int) -> ArticleTable:
    rows = [
        ["Minimal branch ratio", "η₀ = 1/π", _num(ETA_0, digits)],
        ["Minimal shell width", "ℓ₀ = 1/(π√π)", _num(ELL_0, digits)],
        ["Filamentary reference core", "Λ_ind", _num(LAMBDA_IND, digits)],
        ["Gaussian UV constant", "C_UV^(Gauss)", _num(C_UV_GAUSS, digits)],
        ["Mean shell-admission quotient", "P_IR^(χ)(ℓ₀)", _num(P_IR_CHI_MEAN, digits)],
        ["Mean admitted UV load", "ΔΛ_UV→IR^(mean)(ℓ₀)", _num(DELTA_LAMBDA_UV_TO_IR_MEAN, digits)],
        ["Exterior shell subtraction", "ΔΛ_out(η₀)", _num(DELTA_LAMBDA_OUT_EXACT, digits)],
        ["Mean shell-memory factor", "R_χ^(mean)", _num(R_CHI_MEAN, digits)],
        ["Exterior-memory contribution", "ΔΛ_out(η₀) R_χ^(mean)", _num(EXTERIOR_MEMORY_CONTRIBUTION_MEAN, digits)],
        ["Representative vector constant", "Λ_geom^(mean)", _num(LAMBDA_GEOM_MEAN, digits)],
    ]
    return ArticleTable(
        label="tab:vector-exact-status",
        title="Minimal-branch quantities fixed in closed form or direct shell evaluation",
        columns=[TableColumn("Item"), TableColumn("Symbol / factor"), TableColumn("Value", "number")],
        rows=rows,
    )


def vector_exact_source_table(digits: int) -> ArticleTable:
    cumulative = sum(VECTOR_SHELL_SOURCE_WEIGHTS.values())
    tail = 1 - cumulative
    rows = [
        ["Mode-1 source weight", "w₁", _num(VECTOR_SHELL_SOURCE_WEIGHTS[1], digits)],
        ["Mode-3 source weight", "w₃", _num(VECTOR_SHELL_SOURCE_WEIGHTS[3], digits)],
        ["Mode-5 source weight", "w₅", _num(VECTOR_SHELL_SOURCE_WEIGHTS[5], digits)],
        ["Mode-7 source weight", "w₇", _num(VECTOR_SHELL_SOURCE_WEIGHTS[7], digits)],
        ["Cumulative weight through j = 7", "w₁ + w₃ + w₅ + w₇", _num(cumulative, digits)],
        ["Tail above j = 7", "1 - (w₁ + w₃ + w₅ + w₇)", _num(tail, digits)],
    ]
    return ArticleTable(
        label="tab:vector-exact-source",
        title="Normalized odd-toroidal shell-source weights on the minimal branch",
        columns=[TableColumn("Item"), TableColumn("Symbol / factor"), TableColumn("Value", "number")],
        rows=rows,
    )





def vector_mean_input_intermediates_table(digits: int) -> ArticleTable:
    bare = DELTA_LAMBDA_OUT_EXACT
    delta_uv = C_UV_GAUSS * P_IR_CHI_MEAN
    memory_only = bare * (R_CHI_MEAN - 1)
    closure_residual = LAMBDA_IND + delta_uv + bare * R_CHI_MEAN - LAMBDA_GEOM_MEAN
    rows = [
        ["UV shell-admission quotient", "P_IR^(χ)(ℓ₀)", _num(P_IR_CHI_MEAN, digits)],
        ["Gaussian UV constant", "C_UV^(Gauss)", _num(C_UV_GAUSS, digits)],
        ["Admitted UV load", "C_UV^(Gauss) P_IR^(χ)(ℓ₀)", _num(delta_uv, digits)],
        ["Bare exterior subtraction", "ΔΛ_out(η₀)", _num(bare, digits)],
        ["Shell-memory factor", "R_χ^(mean)", _num(R_CHI_MEAN, digits)],
        ["Pure memory correction", "ΔΛ_out(η₀)[R_χ^(mean)-1]", _num(memory_only, digits)],
        ["Exterior-memory contribution", "ΔΛ_out(η₀) R_χ^(mean)", _num(bare * R_CHI_MEAN, digits)],
        ["Closure residual", "Λ_ind + ΔΛ_UV→IR + ΔΛ_out R_χ - Λ_geom", _num(closure_residual, digits)],
    ]
    return ArticleTable(
        label="tab:vector-mean-input-intermediates",
        title="Intermediate representative inputs entering Λ_geom^(mean)",
        columns=[TableColumn("Item"), TableColumn("Formula / factor"), TableColumn("Value", "number")],
        rows=rows,
    )


def vector_lambda_ledger_table(digits: int) -> ArticleTable:
    total = LAMBDA_GEOM_MEAN
    rows = [
        ["Intrinsic core", "Λ_ind", _num(LAMBDA_IND, digits), format_percent(LAMBDA_IND / total, 6)],
        [
            "UV-to-IR shell transfer",
            "ΔΛ_UV→IR^(mean)(ℓ₀)",
            _num(DELTA_LAMBDA_UV_TO_IR_MEAN, digits),
            format_percent(DELTA_LAMBDA_UV_TO_IR_MEAN / total, 6),
        ],
        [
            "Exterior-memory contribution",
            "ΔΛ_out(η₀) R_χ^(mean)",
            _num(EXTERIOR_MEMORY_CONTRIBUTION_MEAN, digits),
            format_percent(EXTERIOR_MEMORY_CONTRIBUTION_MEAN / total, 6),
        ],
        [
            "Total representative geometric constant",
            "Λ_geom^(mean)",
            _num(LAMBDA_GEOM_MEAN, digits),
            format_percent(mp.mpf("1"), 6),
        ],
    ]
    return ArticleTable(
        label="tab:vector-lambda-ledger",
        title="Top-level contribution ledger for Λ_geom^(mean)",
        columns=[
            TableColumn("Contribution"),
            TableColumn("Formula / factor"),
            TableColumn("Value", "number"),
            TableColumn("Share of Λ_geom^(mean)", "number"),
        ],
        rows=rows,
    )



def vector_lambda_out_ledger_table(digits: int) -> ArticleTable:
    bare = DELTA_LAMBDA_OUT_EXACT
    total = EXTERIOR_MEMORY_CONTRIBUTION_MEAN
    memory = total - bare
    rows = [
        [
            "Bare exterior subtraction",
            "ΔΛ_out(η₀)",
            _num(bare, digits),
            format_percent(bare / total, 6),
        ],
        [
            "Memory correction above bare exterior",
            "ΔΛ_out(η₀) [R_χ^(mean) - 1]",
            _num(memory, digits),
            format_percent(memory / total, 6),
        ],
        [
            "Total exterior-memory contribution",
            "ΔΛ_out(η₀) R_χ^(mean)",
            _num(total, digits),
            format_percent(mp.mpf("1"), 6),
        ],
    ]
    return ArticleTable(
        label="tab:vector-lambda-out-ledger",
        title="Internal decomposition of the representative exterior-memory contribution",
        columns=[
            TableColumn("Contribution"),
            TableColumn("Formula / factor"),
            TableColumn("Value", "number"),
            TableColumn("Share of ΔΛ_out(η₀) R_χ^(mean)", "number"),
        ],
        rows=rows,
    )





def alpha_lock_construction_table(digits: int) -> ArticleTable:
    lambda_geom = LAMBDA_GEOM_MEAN
    k_lambda = K_ALP * lambda_geom
    zeta = zeta_from_lambda(lambda_geom)
    linear = mp.mpf('1.5') * k_lambda
    quadratic = linear * zeta
    d_lock = dc_lock_from_lambda(lambda_geom)
    rows = [
        ["Pinned ALP constant", "K", _num(K_ALP, digits)],
        ["Representative vector input", "Λ_geom^(mean)", _num(lambda_geom, digits)],
        ["Linear shell product", "K Λ_geom", _num(k_lambda, digits)],
        ["Quadratic lock ratio", "K Λ_geom / (2π²)", _num(zeta, digits)],
        ["ALP bracket", "1 + K Λ_geom / (2π²)", _num(1 + zeta, digits)],
        ["Linear D-lock contribution", "(3/2) K Λ_geom", _num(linear, digits)],
        ["Quadratic D-lock contribution", "(3/2) K Λ_geom · K Λ_geom/(2π²)", _num(quadratic, digits)],
        ["Locked scalar target", "D_C^lock(Λ_geom^(mean))", _num(d_lock, digits)],
    ]
    return ArticleTable(
        label="tab:alpha-lock-construction",
        title="Closed ALP construction of the locked scalar target D_C^lock",
        columns=[TableColumn("Item"), TableColumn("Formula / factor"), TableColumn("Value", "number")],
        rows=rows,
    )



def alpha_internal_root_table(digits: int) -> ArticleTable:
    model = CURRENT_MODELS[5]
    result = _current_lock_result()
    scalar = solve_dc_article(result.alpha, model)
    d_c = scalar.d_c
    x_value = result.alpha / PI
    numerator = dynamic_numerator(d_c, model)
    denominator = dynamic_denominator(d_c, model)
    sqrt_r = mp.sqrt(scalar.radicand)
    rows = [
        ["Current article lock output", "α_lock^(art)", _num(result.alpha, digits)],
        ["Inverse lock output", "(α_lock^(art))⁻¹", _num(result.inverse_alpha, digits)],
        ["Sommerfeld variable", "x_lock = α_lock^(art)/π", _num(x_value, digits)],
        ["Pinned ALP constant", "K", _num(K_ALP, digits)],
        ["Geometric vector input", "Λ_geom^(mean)", _num(LAMBDA_GEOM_MEAN, digits)],
        ["Vector overlap strength", "ζ_* = K Λ_geom/(2π²)", _num(result.zeta, digits)],
        ["Locked scalar target", "D_C^lock", _num(result.d_lock, digits)],
        ["Source prefactor", "λ_src", _num(model.source_renorm_factor, digits)],
        ["Effective static load", "Θ₁,eff", _num(model.theta1_effective, digits)],
        ["Dynamic numerator", "N(D_C)", _num(numerator, digits)],
        ["Dynamic denominator", "Q(D_C)", _num(denominator, digits)],
        ["Visible dynamic response", "Φ_dyn(D_C)", _num(scalar.phi_dyn, digits)],
        ["Scalar mother radicand", "R_moth(D_C)", _num(scalar.radicand, digits)],
        ["Square root of radicand", "√R_moth(D_C)", _num(sqrt_r, digits)],
        ["Scalar block at the lock output", "D_C(α_lock^(art))", _num(d_c, digits)],
        ["Closure RHS", "x_lock √R_moth(D_C)", _num(x_value * sqrt_r, digits)],
        ["Bridge value", "D_total(α_lock^(art))", _num(d_total_from_d_and_alpha(d_c, result.alpha), digits)],
        ["Scalar lock residual", "D_C(α_lock^(art)) - D_C^lock", _num(d_c - result.d_lock, digits)],
        ["Closure residual", "D_C - x√R_moth", _num(result.residual, digits)],
    ]
    return ArticleTable(
        label="tab:alpha-internal-root",
        title="Internal current-rank-5 values entering the article geometric-lock output",
        columns=[TableColumn("Item"), TableColumn("Symbol"), TableColumn("Value", "number")],
        rows=rows,
    )


def alpha_reference_comparison_table(digits: int) -> ArticleTable:
    result = _current_lock_result()
    delta_direct = result.alpha - ALPHA_BENCHMARK_DIRECT
    delta_inverse = result.inverse_alpha - ALPHA_BENCHMARK_INVERSE
    rows = [
        ["Adopted benchmark in inverse form", "α_exp⁻¹", _num(ALPHA_BENCHMARK_INVERSE, digits)],
        ["Quoted one-sigma uncertainty", "u(α⁻¹)", _num(ALPHA_BENCHMARK_INVERSE_U, digits)],
        ["Derived benchmark in direct form", "α_exp", _num(ALPHA_BENCHMARK_DIRECT, digits)],
        ["Difference in inverse form", "Δ(α⁻¹)", _num(delta_inverse, digits)],
        ["Difference in direct form", "Δα", _num(delta_direct, digits)],
        ["Relative difference", "10¹² Δα / α_lock^(art)", _ppt(delta_direct, result.alpha)],
        ["Deviation in sigma units", "z_σ(α)", _num(delta_direct / ALPHA_BENCHMARK_DIRECT_U, digits)],
    ]
    return ArticleTable(
        label="tab:alpha-reference-comparison",
        title="Diagnostic external comparison of the current article lock output",
        columns=[TableColumn("Item"), TableColumn("Symbol"), TableColumn("Value", "number")],
        rows=rows,
    )



def alpha_experimental_benchmark_table(digits: int) -> ArticleTable:
    result = _current_lock_result()
    rows = []
    for reference in EXPERIMENTAL_ALPHA_REFERENCES:
        delta = result.alpha - reference.alpha_ref
        rows.append(
            [
                reference.year,
                reference.method,
                reference.reference_label,
                _num(reference.alpha_ref, digits),
                _num(reference.u_ref, digits),
                _signed_num(delta, digits),
                f"{float(delta / reference.u_ref):+.2f}",
            ]
        )
    return ArticleTable(
        label="tab:alpha-expt-benchmark",
        title="Experimental and adjusted reference values of α against the current article lock output",
        columns=[
            TableColumn("Year"),
            TableColumn("Method"),
            TableColumn("Ref."),
            TableColumn("α_ref", "number"),
            TableColumn("u_ref", "number"),
            TableColumn("Δα", "number"),
            TableColumn("Δα / u_ref", "number"),
        ],
        rows=rows,
    )





def qed_universal_input_coefficients_table(digits: int) -> ArticleTable:
    coefficients = _qed_coefficients(5)
    rows = [[str(n), _num(coefficients.a1_coefficients[n], digits)] for n in range(1, coefficients.order + 1)]
    return ArticleTable(
        label="tab:qed-universal-input-coefficients",
        title="Universal pure-photonic QED coefficients used in the one-way scalar cross-check",
        columns=[TableColumn("n", "number"), TableColumn("A₁^(2n)", "number")],
        rows=rows,
    )



def qed_induced_scalar_coefficients_table(digits: int) -> ArticleTable:
    coefficients = _qed_coefficients(5)
    rows = [
        [str(n), _num(coefficients.d_coefficients[n], digits), _num(coefficients.c_coefficients[n], digits)]
        for n in range(1, coefficients.order + 1)
    ]
    return ArticleTable(
        label="tab:qed-induced-scalar-coefficients",
        title="Intermediate QED-induced scalar coefficients entering D_C^QED,[5](x)",
        columns=[
            TableColumn("n", "number"),
            TableColumn("d_n^QED", "number"),
            TableColumn("c_n^QED", "number"),
        ],
        rows=rows,
    )



def qed_lock_evaluation_table(digits: int) -> ArticleTable:
    coefficients = _qed_coefficients(5)
    result = _qed_lock_result()
    d_qed = dc_qed_from_alpha(result.alpha, coefficients)
    rows = [
        ["Cross-check solution", "α_QED→sc^[5]", _num(result.alpha, digits)],
        ["Inverse cross-check solution", "(α_QED→sc^[5])⁻¹", _num(result.inverse_alpha, digits)],
        ["Sommerfeld variable", "x_QED = α_QED→sc^[5]/π", _num(result.alpha / PI, digits)],
        ["Locked scalar target", "D_C^lock", _num(result.d_lock, digits)],
        ["QED-induced scalar value", "D_C^QED,[5](x_QED)", _num(d_qed, digits)],
        ["Cross-check residual", "D_C^QED,[5](x_QED) - D_C^lock", _num(d_qed - result.d_lock, digits)],
    ]
    return ArticleTable(
        label="tab:qed-lock-evaluation",
        title="Internal values of the universal-QED-induced one-way lock cross-check",
        columns=[TableColumn("Item"), TableColumn("Symbol"), TableColumn("Value", "number")],
        rows=rows,
    )


def alpha_qed_crosscheck_table(digits: int) -> ArticleTable:
    article = _current_lock_result()
    qed = _qed_lock_result()
    rows = [
        [
            "Current article lock output",
            _num(article.alpha, digits),
            _num(article.alpha - ALPHA_BENCHMARK_DIRECT, digits),
            f"{float((article.alpha - ALPHA_BENCHMARK_DIRECT) / ALPHA_BENCHMARK_DIRECT_U):.3f}",
        ],
        [
            "Universal-QED cross-check",
            _num(qed.alpha, digits),
            _num(qed.alpha - ALPHA_BENCHMARK_DIRECT, digits),
            f"{float((qed.alpha - ALPHA_BENCHMARK_DIRECT) / ALPHA_BENCHMARK_DIRECT_U):.3f}",
        ],
    ]
    note = (
        "The current spread between the article scalar branch and the universal-QED-induced "
        f"branch at the same geometric lock is Δα = {_num(article.alpha - qed.alpha, digits)} "
        f"({_ppt(article.alpha - qed.alpha, article.alpha)})."
    )
    return ArticleTable(
        label="tab:alpha-qed-crosscheck",
        title="Comparison of the current article lock output and the universal-QED inverse-α cross-check",
        columns=[
            TableColumn("Realization"),
            TableColumn("α", "number"),
            TableColumn("Δα_exp", "number"),
            TableColumn("n_σ", "number"),
        ],
        rows=rows,
        note=note,
    )



def logical_status_summary_table() -> ArticleTable:
    rows = [
        ["Standing input", "ℂ ⊕ ℝ³, Rω = c, shell admissibility, analyticity, gauge regularity"],
        ["Derived branch data", "σ_ℂ = R/√π, r_* = πR, η = 1/π, ℓ = 1/(π√π)"],
        ["Derived infrared data", "Q = n q₀, Dirac limit, g = 2, dim-5 Pauli slot"],
        ["Derived scalar block", "D_C(α) from the cavity–shell pencil"],
        ["Derived vector geometry", "R_χ, Λ_geom"],
        ["Derived shell lock", "C_log = 1/3, ALP"],
        ["Primary output", "D_C^lock(Λ_geom), α_pred"],
        ["Secondary consequence", "pure-photonic Â₁^(2n) from the realized evaluator"],
    ]
    return ArticleTable(
        label="tab:logical-status-summary",
        title="Logical status of the construction",
        columns=[TableColumn("Logical status"), TableColumn("Content")],
        rows=rows,
    )



def g2_rank_summary_table(digits: int) -> ArticleTable:
    rows = []
    models = [
        ("Current rank-1", CURRENT_MODELS[1]),
        ("Refined rank-1", REFINED_MODELS[1]),
        ("Current rank-3", CURRENT_MODELS[3]),
        ("Refined rank-3", REFINED_MODELS[3]),
        ("Current rank-5", CURRENT_MODELS[5]),
        ("Refined rank-5", REFINED_MODELS[5]),
        ("Current rank-100", CURRENT_MODELS[100]),
        ("Refined rank-100", REFINED_MODELS[100]),
    ]
    for label, model in models:
        physical = evaluate_g2_at_alpha(ALPHA_BENCHMARK_DIRECT, model)
        coefficients = _g2_coefficients(model.key, 5)
        max_delta = max(abs(coefficients[n] - QED_A1_COEFFICIENTS[n]) for n in range(1, 6))
        rows.append(
            [
                label,
                _num(physical.d_c, digits),
                _num(physical.g_e2, digits),
                _num(physical.g_e2 - QED_G_E2_BENCHMARK_XSTAR, digits),
                _num(max_delta, digits),
            ]
        )
    return ArticleTable(
        label="tab:g2-rank-summary",
        title="Pure-photonic rank summary for the current and refined scalar realizations",
        columns=[
            TableColumn("Model"),
            TableColumn("D_C(x_*)", "number"),
            TableColumn("g_e,2(x_*)", "number"),
            TableColumn("Δg_e,2(x_*)", "number"),
            TableColumn("max₁≤n≤5 |ΔA₁^(2n)|", "number"),
        ],
        rows=rows,
    )




def g2_physical_point_table(digits: int) -> ArticleTable:
    physical = evaluate_g2_at_alpha(ALPHA_BENCHMARK_DIRECT, REFINED_MODELS[100])
    linear_drag = physical.x * physical.zeta
    quadratic_backreaction = physical.x**2 * physical.zeta**2 / (4 * physical.d_c)
    one_minus_dtotal = 1 - physical.d_total
    rows = [
        ["α_ref", _num(physical.alpha, digits)],
        ["x_* = α_ref/π", _num(physical.x, digits)],
        ["D_C(x_*)", _num(physical.d_c, digits)],
        ["ζ(D_C(x_*))", _num(physical.zeta, digits)],
        ["x_* ζ(D_C(x_*))", _num(linear_drag, digits)],
        ["x_*² ζ(D_C(x_*))² /(4 D_C)", _num(quadratic_backreaction, digits)],
        ["Λ_π(x_*)", _num(physical.lambda_pi, digits)],
        ["D_total(x_*)", _num(physical.d_total, digits)],
        ["1 - D_total(x_*)", _num(one_minus_dtotal, digits)],
        ["√(1 - D_total(x_*))", _num(mp.sqrt(one_minus_dtotal), digits)],
        ["g_e,2^Rel(x_*)", _num(physical.g_e2, digits)],
        ["a_e^Rel(x_*)", _num(physical.a_e, digits)],
    ]
    return ArticleTable(
        label="tab:g2-physical-point",
        title="Physical-point bridge intermediates entering the pure-photonic Relator cross-check (refined rank-100)",
        columns=[TableColumn("Quantity"), TableColumn("Value", "number")],
        rows=rows,
    )




def g2_bridge_intermediates_table(digits: int) -> ArticleTable:
    physical = evaluate_g2_at_alpha(ALPHA_BENCHMARK_DIRECT, REFINED_MODELS[100])
    rows = [
        ["Benchmark g_e / 2", _num(QED_G_E2_BENCHMARK_XSTAR, digits), "pure-QED reference"],
        ["Relator g_e / 2", _num(physical.g_e2, digits), "refined rank-100"],
        ["Δ(g_e / 2)", _num(physical.g_e2 - QED_G_E2_BENCHMARK_XSTAR, digits), "Relator - benchmark"],
        ["Benchmark a_e", _num(QED_A_E_BENCHMARK_XSTAR, digits), "pure-QED reference"],
        ["Relator a_e", _num(physical.a_e, digits), "refined rank-100"],
        ["Δa_e", _num(physical.a_e - QED_A_E_BENCHMARK_XSTAR, digits), "Relator - benchmark"],
    ]
    return ArticleTable(
        label="tab:g2-bridge-intermediates",
        title="Final physical-point comparison of the Relator bridge against the pure-photonic benchmark",
        columns=[TableColumn("Quantity"), TableColumn("Value", "number"), TableColumn("Comment")],
        rows=rows,
    )


def g2_known_coefficients_table(digits: int) -> ArticleTable:
    current = _g2_coefficients(CURRENT_MODELS[5].key, 5)
    refined = _g2_coefficients(REFINED_MODELS[100].key, 5)
    rows = []
    for n in range(1, 6):
        pure = QED_A1_COEFFICIENTS[n]
        rows.append(
            [
                str(n),
                _num(pure, digits),
                _num(current[n], digits),
                _num(refined[n], digits),
                _num(current[n] - pure, digits),
                _num(refined[n] - pure, digits),
            ]
        )
    return ArticleTable(
        label="tab:g2-known-coefficients",
        title="Comparison of the first five nontrivial pure-photonic electron-anomaly coefficients",
        columns=[
            TableColumn("n", "number"),
            TableColumn("A₁^(2n) (pure QED)", "number"),
            TableColumn("Â₁^(2n) (current rank-5)", "number"),
            TableColumn("Â₁^(2n) (refined rank-100)", "number"),
            TableColumn("ΔA₁^(2n) (current rank-5)", "number"),
            TableColumn("ΔA₁^(2n) (refined rank-100)", "number"),
        ],
        rows=rows,
    )



def g2_first_ten_predictions_table(digits: int) -> ArticleTable:
    coefficients = _g2_coefficients(REFINED_MODELS[100].key, 10)
    rows = []
    for n in range(1, 11):
        if n in QED_A1_COEFFICIENTS:
            pure = _num(QED_A1_COEFFICIENTS[n], digits)
            status = "benchmarked"
        else:
            pure = "---"
            status = "Relator prediction"
        rows.append([str(n), _num(coefficients[n], digits), pure, status])
    return ArticleTable(
        label="tab:g2-first-ten-predictions",
        title="First ten nontrivial coefficients of the refined rank-100 Relator prediction",
        columns=[
            TableColumn("n", "number"),
            TableColumn("Â₁^(2n) (Relator)", "number"),
            TableColumn("A₁^(2n) (pure QED)"),
            TableColumn("Status"),
        ],
        rows=rows,
    )



def app_closed_vs_truncated_table() -> ArticleTable:
    rows = [
        ["η₀, ℓ₀, K", "exact closed", "no truncation"],
        ["C₁₁(ℓ₀), γ_cℓ(D)", "exact closed", "error-function formulas"],
        ["B₁₁^glue", "exact closed", "visible 2-channel static block"],
        ["B₁₁^(3,ρ)", "closed in the chosen completion", "canonical 3-channel static lift"],
        ["D_C^lock(Λ)", "exact closed", "no root search"],
        ["α(D)", "exact closed", "once R_moth is fixed"],
        ["ζ(D), D_total(D)", "exact closed", "algebraic Relator bridge"],
        ["Λ_geom^(mean)", "representative", "target-free vector-shell input"],
        ["(A_uv^(N), A_ir^(N), χ_N)", "finite-rank", "visible dynamic block"],
        ["α_QED→sc^[5]", "external truncated", "QED-induced cross-check"],
    ]
    return ArticleTable(
        label="tab:app-closed-vs-truncated",
        title="Status of the ingredients entering the realized geometric lock",
        columns=[TableColumn("Ingredient"), TableColumn("Status"), TableColumn("Comment")],
        rows=rows,
    )





def app_sensitivity_basepoint_table(digits: int) -> ArticleTable:
    model = REALIZED_MODELS[100]
    rows = [
        ["K", _num(K_ALP, digits)],
        ["Λ_geom^(mean)", _num(LAMBDA_GEOM_MEAN, digits)],
        ["λ_src", _num(model.source_renorm_factor, digits)],
        ["Θ₁,eff", _num(model.theta1_effective, digits)],
        ["B₁₁^(3,ρ)", _num(model.b11_effective, digits)],
        ["A_uv^(100)", _num(model.a_uv, digits)],
        ["A_ir^(100)", _num(model.a_ir, digits)],
        ["χ_100", _num(model.chi, digits)],
        ["S_UV", _num(S_UV, digits)],
        ["S_IR", _num(S_IR, digits)],
    ]
    return ArticleTable(
        label="tab:app-sensitivity-basepoint",
        title="Base-point parameters used in the appendix sensitivity audit",
        columns=[TableColumn("Quantity"), TableColumn("Value", "number")],
        rows=rows,
    )



def app_exact_sensitivities_table(digits: int) -> ArticleTable:
    model = REALIZED_MODELS[100]
    report = exact_sensitivities(model)
    rows = [
        ["K", _num(K_ALP, digits)],
        ["Λ_geom^(mean)", _num(LAMBDA_GEOM_MEAN, digits)],
        ["source prefactor", _num(model.source_renorm_factor, digits)],
        ["Θ₁,eff", _num(model.theta1_effective, digits)],
        ["B₁₁^(3,ρ)", _num(model.b11_effective, digits)],
        ["D_* = D_C^lock(Λ_geom^(mean))", _num(report.d_lock, digits)],
        ["R_* = R_moth(D_*)", _num(report.radicand_at_lock, digits)],
        ["d D_C^lock / dΛ", _num(report.d_dlock_d_lambda, digits)],
        ["d D_C^lock / dK", _num(mp.mpf('1.5') * LAMBDA_GEOM_MEAN * (1 + K_ALP * LAMBDA_GEOM_MEAN / PI**2), digits)],
        ["dα / dD at D_*", _num(report.d_alpha_d_d, digits)],
        ["dα_pred / dΛ_geom", _num(report.d_alpha_d_lambda, digits)],
        ["dα_pred / dK", _num(report.d_alpha_d_k, digits)],
        ["∂ ln α_pred / ∂ ln Λ_geom", _num(report.dlog_alpha_dlog_lambda, digits)],
        ["∂ ln α_pred / ∂ ln K", _num(report.dlog_alpha_dlog_k, digits)],
        ["∂ α_pred / ∂ Θ₁,eff", _num(report.d_alpha_d_theta_eff, digits)],
        ["∂ ln α_pred / ∂ ln Θ₁,eff", _num(report.dlog_alpha_dlog_theta_eff, digits)],
        ["∂ ln α_pred / ∂ ln B₁₁^(3,ρ)", _num(report.dlog_alpha_dlog_b11, digits)],
    ]
    return ArticleTable(
        label="tab:app-exact-sensitivities",
        title="Exact closed sensitivities at the realized geometric lock",
        columns=[TableColumn("Quantity"), TableColumn("Value", "number")],
        rows=rows,
    )



def app_numerical_sensitivities_table(digits: int) -> ArticleTable:
    model = REALIZED_MODELS[100]
    report = numerical_sensitivities(model)
    rows = [
        ["Θ₁,eff", _num(model.theta1_effective, digits), _num(report.theta_eff, digits)],
        ["A_uv^(100)", _num(model.a_uv, digits), _num(report.a_uv, digits)],
        ["A_ir^(100)", _num(model.a_ir, digits), _num(report.a_ir, digits)],
        ["χ_100", _num(model.chi, digits), _num(report.chi, digits)],
        ["S_UV", _num(S_UV, digits), _num(report.s_uv, digits)],
        ["S_IR", _num(S_IR, digits), _num(report.s_ir, digits)],
    ]
    return ArticleTable(
        label="tab:app-numerical-sensitivities",
        title="Numerical logarithmic sensitivities of the realized scalar truncation",
        columns=[TableColumn("Parameter"), TableColumn("Base value", "number"), TableColumn("S_q^(10⁻⁸)", "number")],
        rows=rows,
    )


def app_rank_stability_table(digits: int) -> ArticleTable:
    stability = rank_stability()
    rows = []
    for rank in RANKS:
        alpha_value, drift = stability[rank]
        rows.append([str(rank), _num(alpha_value, digits), _num(drift, digits)])
    return ArticleTable(
        label="tab:app-rank-stability",
        title="Rank-stability of the realized geometric-lock prediction",
        columns=[
            TableColumn("Rank N", "number"),
            TableColumn("α_pred^(N)", "number"),
            TableColumn("α_pred^(N) - α_pred^(100)", "number"),
        ],
        rows=rows,
    )



def app_precision_stability_table(digits: int) -> ArticleTable:
    audit = precision_audit()
    rows = []
    for dps in (50, 80, 120, 200, 300):
        alpha_value, drift = audit[dps]
        rows.append([f"{dps} dps", _num(alpha_value, digits), _num(drift, digits)])
    return ArticleTable(
        label="tab:app-precision-stability",
        title="Arithmetic-precision audit of the realized geometric-lock prediction",
        columns=[TableColumn("Working precision"), TableColumn("α_pred", "number"), TableColumn("drift from 300 dps", "number")],
        rows=rows,
    )





def all_article_tables(digits: int = 24) -> list[ArticleTable]:
    return [
        scalar_static_fixed_data_table(digits),
        scalar_source_prefactor_convention_table(digits),
        scalar_visible_dynamic_amplitudes_table(digits),
        scalar_visible_dynamic_block_table(digits),
        scalar_current_physical_intermediates_table(digits),
        scalar_hidden_dynamic_jacobi_table(digits),
        scalar_physical_roots_table(digits),
        scalar_coefficients_table(digits),
        vector_exact_status_table(digits),
        vector_exact_source_table(digits),
        vector_mean_input_intermediates_table(digits),
        vector_lambda_ledger_table(digits),
        vector_lambda_out_ledger_table(digits),
        alpha_lock_construction_table(digits),
        alpha_internal_root_table(digits),
        alpha_reference_comparison_table(digits),
        alpha_experimental_benchmark_table(digits),
        qed_universal_input_coefficients_table(digits),
        qed_induced_scalar_coefficients_table(digits),
        qed_lock_evaluation_table(digits),
        alpha_qed_crosscheck_table(digits),
        logical_status_summary_table(),
        g2_rank_summary_table(digits),
        g2_physical_point_table(digits),
        g2_bridge_intermediates_table(digits),
        g2_known_coefficients_table(digits),
        g2_first_ten_predictions_table(digits),
        app_closed_vs_truncated_table(),
        app_sensitivity_basepoint_table(digits),
        app_exact_sensitivities_table(digits),
        app_numerical_sensitivities_table(digits),
        app_rank_stability_table(digits),
        app_precision_stability_table(digits),
    ]
