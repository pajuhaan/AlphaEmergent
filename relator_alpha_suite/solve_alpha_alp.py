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

from relator_alpha.common import configure_precision, format_mpf
from relator_alpha.presentation import build_console
from relator_alpha.article_constants import K_ALP, LAMBDA_GEOM_MEAN, SCALAR_MODEL_PRESETS
from relator_alpha.alp import (
    d_total_from_d_and_alpha,
    dc_lock_from_lambda,
    solve_alpha_from_article_branch,
    solve_alpha_from_qed_branch,
    zeta_from_lambda,
)
from relator_alpha.scalar_article import dynamic_denominator, dynamic_numerator, solve_dc_article
from relator_alpha.scalar_qed import build_qed_induced_coefficients, dc_qed_from_alpha



def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Solve the ALP closure for α using either the article scalar branch or the one-way QED cross-check."
    )
    parser.add_argument(
        "--branch",
        choices=["article", "qed"],
        default="article",
        help="Scalar branch used in the ALP inversion.",
    )
    parser.add_argument(
        "--model",
        choices=sorted(SCALAR_MODEL_PRESETS.keys()),
        default="current_rank5",
        help="Article scalar preset used when --branch=article.",
    )
    parser.add_argument(
        "--order",
        type=int,
        default=5,
        help="QED truncation order used when --branch=qed.",
    )
    parser.add_argument(
        "--lambda-geom",
        type=str,
        default=str(LAMBDA_GEOM_MEAN),
        help="Vector-channel geometric input Λ_geom.",
    )
    parser.add_argument("--digits", type=int, default=24, help="Displayed significant digits.")
    parser.add_argument("--dps", type=int, default=160, help="Working precision in decimal digits.")
    parser.add_argument("--width", type=int, default=280, help="Console width used for Rich tables.")
    return parser



def main() -> None:
    args = build_parser().parse_args()
    configure_precision(args.dps)
    console = build_console(args.width)
    lambda_geom = mp.mpf(args.lambda_geom)
    d_lock = dc_lock_from_lambda(lambda_geom)
    k_lambda = K_ALP * lambda_geom
    zeta_lambda = zeta_from_lambda(lambda_geom)

    if args.branch == "article":
        model = SCALAR_MODEL_PRESETS[args.model]
        result = solve_alpha_from_article_branch(lambda_geom, model)
        scalar = solve_dc_article(result.alpha, model)
        dc_value = scalar.d_c
        numerator = dynamic_numerator(dc_value, model)
        denominator = dynamic_denominator(dc_value, model)
        branch_comment = model.description
        branch_descriptor = f"article scalar branch / {model.key}"
        source_prefactor = model.source_renorm_factor
        source_label = model.source_renorm_label
    else:
        coefficients = build_qed_induced_coefficients(order=args.order)
        result = solve_alpha_from_qed_branch(lambda_geom, coefficients)
        dc_value = dc_qed_from_alpha(result.alpha, coefficients)
        numerator = None
        denominator = None
        branch_comment = f"truncated universal-QED branch, order M = {coefficients.order}"
        branch_descriptor = f"QED-induced cross-check / M = {coefficients.order}"
        source_prefactor = None
        source_label = None

    meta_table = Table(title="ALP solve configuration")
    meta_table.add_column("Field")
    meta_table.add_column("Value", justify="right", no_wrap=True)
    meta_table.add_column("Comment")
    meta_table.add_row("branch", branch_descriptor, branch_comment)
    meta_table.add_row("K", format_mpf(K_ALP, args.digits), "closed ALP shell constant")
    meta_table.add_row("Λ_geom", format_mpf(lambda_geom, args.digits), "vector-channel geometric input")
    meta_table.add_row("K Λ_geom", format_mpf(k_lambda, args.digits), "linear ALP shell product")
    meta_table.add_row("K Λ_geom /(2π²)", format_mpf(zeta_lambda, args.digits), "closed overlap ratio")
    if source_prefactor is not None:
        meta_table.add_row("source prefactor", format_mpf(source_prefactor, args.digits), source_label)

    alp_table = Table(title="ALP target and α solution")
    alp_table.add_column("Quantity")
    alp_table.add_column("Value", justify="right", no_wrap=True)
    alp_table.add_column("Comment")
    alp_table.add_row("ζ", format_mpf(result.zeta, args.digits), "baseline overlap strength")
    alp_table.add_row("D_C^lock", format_mpf(result.d_lock, args.digits), "ALP target")
    alp_table.add_row("α", format_mpf(result.alpha, args.digits), branch_comment)
    alp_table.add_row("α⁻¹", format_mpf(result.inverse_alpha, args.digits), "inverse fine-structure constant")
    alp_table.add_row("x = α/π", format_mpf(result.alpha / mp.pi, args.digits), "Sommerfeld variable")
    alp_table.add_row("D_C(branch)", format_mpf(dc_value, args.digits), "scalar-branch value at the solved α")
    if numerator is not None and denominator is not None:
        alp_table.add_row("N(D_C)", format_mpf(numerator, args.digits), "dynamic numerator")
        alp_table.add_row("Q(D_C)", format_mpf(denominator, args.digits), "dynamic denominator")
        alp_table.add_row("R_moth(D_C)", format_mpf(scalar.radicand, args.digits), "scalar mother radicand")
    alp_table.add_row(
        "D_total",
        format_mpf(d_total_from_d_and_alpha(dc_value, result.alpha), args.digits),
        "Relator bridge value",
    )
    alp_table.add_row("lock residual", format_mpf(result.residual, args.digits), result.note)

    console.print(meta_table)
    console.print()
    console.print(alp_table)


if __name__ == "__main__":
    main()
