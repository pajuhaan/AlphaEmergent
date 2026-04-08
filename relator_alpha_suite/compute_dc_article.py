#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import mpmath as mp
from rich.table import Table

from relator_alpha.article_constants import LAMBDA_GEOM_MEAN, SCALAR_MODEL_PRESETS
from relator_alpha.alp import solve_alpha_from_article_branch
from relator_alpha.common import configure_precision, format_mpf
from relator_alpha.presentation import build_console
from relator_alpha.scalar_article import solve_dc_article



def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compute D_C(α) from one of the article scalar evaluators."
    )
    parser.add_argument(
        "--model",
        choices=sorted(SCALAR_MODEL_PRESETS.keys()),
        default="current_rank5",
        help="Scalar evaluator preset.",
    )
    parser.add_argument(
        "--alpha",
        type=str,
        default=None,
        help="Direct α input. If omitted, α is derived internally from the ALP target.",
    )
    parser.add_argument(
        "--lambda-geom",
        type=str,
        default=str(LAMBDA_GEOM_MEAN),
        help="Λ_geom used when α is derived internally.",
    )
    parser.add_argument("--digits", type=int, default=24)
    parser.add_argument("--dps", type=int, default=160)
    parser.add_argument("--width", type=int, default=200)
    return parser



def main() -> None:
    args = build_parser().parse_args()
    configure_precision(args.dps)
    console = build_console(args.width)
    model = SCALAR_MODEL_PRESETS[args.model]

    if args.alpha is None:
        alpha_result = solve_alpha_from_article_branch(mp.mpf(args.lambda_geom), model)
        alpha_value = alpha_result.alpha
        alpha_source = "internally derived from the ALP target"
    else:
        alpha_value = mp.mpf(args.alpha)
        alpha_source = "user-supplied α"

    result = solve_dc_article(alpha_value, model)

    parameter_table = Table(title=f"Article scalar preset: {model.key}")
    parameter_table.add_column("Parameter")
    parameter_table.add_column("Value", justify="right", no_wrap=True)
    parameter_table.add_column("Comment")
    parameter_table.add_row("family", model.family, model.description)
    parameter_table.add_row("rank N", str(model.rank_n), "visible excited-shell truncation rank")
    parameter_table.add_row("B₁₁,eff", format_mpf(model.b11_effective, args.digits), "effective static shell load")
    parameter_table.add_row("Θ₁,eff", format_mpf(model.theta1_effective, args.digits), "Θ₁,core + B₁₁,eff")
    parameter_table.add_row("A_uv^(N)", format_mpf(model.a_uv, args.digits), "visible UV pole weight")
    parameter_table.add_row("A_ir^(N)", format_mpf(model.a_ir, args.digits), "visible IR pole weight")
    parameter_table.add_row("χ_N", format_mpf(model.chi, args.digits), "mixing parameter")
    parameter_table.add_row("ρ_dyn^(N)", format_mpf(model.rho_dyn, args.digits), "dynamic correlation ratio")
    parameter_table.add_row("λ_u^(main text)", format_mpf(model.lambda_u_display, args.digits), "published visible-source factor")
    parameter_table.add_row("source prefactor", format_mpf(model.source_renorm_factor, args.digits), model.source_renorm_label)

    result_table = Table(title="Article scalar evaluation")
    result_table.add_column("Quantity")
    result_table.add_column("Value", justify="right", no_wrap=True)
    result_table.add_column("Comment")
    result_table.add_row("α", format_mpf(result.alpha, args.digits), alpha_source)
    result_table.add_row("x = α/π", format_mpf(result.x, args.digits), "Sommerfeld variable")
    result_table.add_row("D_C(α)", format_mpf(result.d_c, args.digits), "physical positive branch")
    result_table.add_row("Φ_dyn(D_C)", format_mpf(result.phi_dyn, args.digits), "visible dynamic response")
    result_table.add_row("R_moth(D_C)", format_mpf(result.radicand, args.digits), "scalar mother radicand")
    result_table.add_row("closure residual", format_mpf(result.residual, args.digits), "D_C - x sqrt(R_moth)")

    console.print(parameter_table)
    console.print()
    console.print(result_table)


if __name__ == "__main__":
    main()
