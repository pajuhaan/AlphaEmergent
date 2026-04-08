from __future__ import annotations

from functools import lru_cache

import mpmath as mp

from .alp import d_total_from_d_and_alpha, solve_alpha_from_article_branch, solve_alpha_from_qed_branch
from .article_constants import (
    ALPHA_BENCHMARK_DIRECT,
    ALPHA_BENCHMARK_DIRECT_U,
    ALPHA_BENCHMARK_INVERSE,
    ALPHA_BENCHMARK_INVERSE_U,
    ALPHA_CODATA22_DIRECT,
    B11_3RHO,
    B11_GLUE,
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
)
from .common import format_mpf, format_percent
from .diagnostics import exact_sensitivities, numerical_sensitivities, precision_audit, rank_stability
from .g2 import evaluate_g2_at_alpha, relator_g2_coefficients
from .presentation import ArticleTable, TableColumn
from .scalar_article import scalar_series_coefficients, solve_dc_article
from .scalar_qed import build_qed_induced_coefficients


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
def _scalar_coefficients(model_key: str, order: int) -> dict[int, mp.mpf]:
    model_lookup = {**{m.key: m for m in CURRENT_MODELS.values()}, **{m.key: m for m in REFINED_MODELS.values()}, **{m.key: m for m in REALIZED_MODELS.values()}}
    return scalar_series_coefficients(model_lookup[model_key], order)


@lru_cache(maxsize=None)
def _g2_coefficients(model_key: str, order: int) -> dict[int, mp.mpf]:
    model_lookup = {**{m.key: m for m in CURRENT_MODELS.values()}, **{m.key: m for m in REFINED_MODELS.values()}, **{m.key: m for m in REALIZED_MODELS.values()}}
    return relator_g2_coefficients(model_lookup[model_key], order)



def scalar_static_fixed_data_table(digits: int) -> ArticleTable:
    rows = [
        ["η₀", "1/π", _num(ETA_0, digits)],
        ["ℓ₀", "1/(π√π)", _num(ELL_0, digits)],
        ["Θ₁,core", "2π C_uni", _num(THETA_1_CORE, digits)],
        ["Θ₁^(coll)", "(γ_E/π²) C₁₁(ℓ₀)", _num(THETA_1_COLL, digits)],
        ["Θ₁^(log)", "η₀²/[8(1-η₀²/4)]", _num(THETA_1_LOG, digits)],
        ["γ_cℓ(0)", "static collar–log overlap", _num(GAMMA_CL_0, digits)],
        ["ρ_stat", "γ_cℓ(0)/√(N_c N_ℓ)", _num(RHO_STAT, digits)],
        ["Λ_h", "(1-ρ_stat²)⁻¹", _num(LAMBDA_H, digits)],
        ["B₁₁^glue", "visible 2-channel static load", _num(B11_GLUE, digits)],
        ["B₁₁^(3,ρ)", "3-channel static completion", _num(B11_3RHO, digits)],
        ["Θ₁^(3,ρ)", "Θ₁,core + B₁₁^(3,ρ)", _num(THETA_1_3RHO, digits)],
    ]
    return ArticleTable(
        label="tab:dc-static-fixed-data",
        title="Fixed shell geometry and static hidden-mode data",
        columns=[TableColumn("Symbol"), TableColumn("Def."), TableColumn("Value", "number")],
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
        ["Exterior shell subtraction", "ΔΛ_out(η₀)", _num(DELTA_LAMBDA_OUT_EXACT, digits)],
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



def alpha_internal_root_table(digits: int) -> ArticleTable:
    result = _current_lock_result()
    d_c = solve_dc_article(result.alpha, CURRENT_MODELS[5]).d_c
    rows = [
        ["Current article lock output", "α_lock^(art)", _num(result.alpha, digits)],
        ["Inverse lock output", "(α_lock^(art))⁻¹", _num(result.inverse_alpha, digits)],
        ["Geometric vector input", "Λ_geom^(mean)", _num(LAMBDA_GEOM_MEAN, digits)],
        ["Vector overlap strength", "ζ_*", _num(result.zeta, digits)],
        ["Locked scalar value", "D_C^lock", _num(result.d_lock, digits)],
        ["Scalar block at the lock output", "D_C(α_lock^(art))", _num(d_c, digits)],
        ["Bridge value", "D_total(α_lock^(art))", _num(d_total_from_d_and_alpha(d_c, result.alpha), digits)],
        ["Scalar lock residual", "D_C(α_lock^(art)) - D_C^lock", _num(result.residual, digits)],
    ]
    return ArticleTable(
        label="tab:alpha-internal-root",
        title="Internal values of the current article geometric-lock output",
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
    rows = [
        ["D_C(x_*)", _num(physical.d_c, digits)],
        ["ζ(D_C(x_*))", _num(physical.zeta, digits)],
        ["Λ_π(x_*)", _num(physical.lambda_pi, digits)],
        ["D_total(x_*)", _num(physical.d_total, digits)],
        ["g_e,2^Rel(x_*)", _num(physical.g_e2, digits)],
        ["a_e^Rel(x_*)", _num(physical.a_e, digits)],
    ]
    return ArticleTable(
        label="tab:g2-physical-point",
        title="Physical-point values entering the pure-photonic Relator cross-check (refined rank-100)",
        columns=[TableColumn("Quantity"), TableColumn("Value", "number")],
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



def app_exact_sensitivities_table(digits: int) -> ArticleTable:
    report = exact_sensitivities(REALIZED_MODELS[100])
    rows = [
        ["D_* = D_C^lock(Λ_geom^(mean))", _num(report.d_lock, digits)],
        ["R_* = R_moth(D_*)", _num(report.radicand_at_lock, digits)],
        ["d D_C^lock / dΛ", _num(report.d_dlock_d_lambda, digits)],
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
        title="Exact closed sensitivities at the geometric lock",
        columns=[TableColumn("Quantity"), TableColumn("Value", "number")],
        rows=rows,
    )



def app_numerical_sensitivities_table(digits: int) -> ArticleTable:
    report = numerical_sensitivities(REALIZED_MODELS[100])
    rows = [
        ["Θ₁,eff", _num(report.theta_eff, digits)],
        ["A_uv^(100)", _num(report.a_uv, digits)],
        ["A_ir^(100)", _num(report.a_ir, digits)],
        ["χ_100", _num(report.chi, digits)],
        ["s_uv", _num(report.s_uv, digits)],
        ["s_ir", _num(report.s_ir, digits)],
    ]
    return ArticleTable(
        label="tab:app-numerical-sensitivities",
        title="Numerical logarithmic sensitivities of the realized scalar truncation",
        columns=[TableColumn("Parameter"), TableColumn("S_q^(10⁻⁸)", "number")],
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
        scalar_visible_dynamic_amplitudes_table(digits),
        scalar_visible_dynamic_block_table(digits),
        scalar_hidden_dynamic_jacobi_table(digits),
        scalar_physical_roots_table(digits),
        scalar_coefficients_table(digits),
        vector_exact_status_table(digits),
        vector_exact_source_table(digits),
        vector_lambda_ledger_table(digits),
        vector_lambda_out_ledger_table(digits),
        alpha_internal_root_table(digits),
        alpha_reference_comparison_table(digits),
        alpha_experimental_benchmark_table(digits),
        alpha_qed_crosscheck_table(digits),
        logical_status_summary_table(),
        g2_rank_summary_table(digits),
        g2_physical_point_table(digits),
        g2_known_coefficients_table(digits),
        g2_first_ten_predictions_table(digits),
        app_closed_vs_truncated_table(),
        app_exact_sensitivities_table(digits),
        app_numerical_sensitivities_table(digits),
        app_rank_stability_table(digits),
        app_precision_stability_table(digits),
    ]
