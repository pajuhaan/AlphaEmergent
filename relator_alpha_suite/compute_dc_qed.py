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

from relator_alpha.article_constants import LAMBDA_GEOM_MEAN
from relator_alpha.alp import solve_alpha_from_qed_branch
from relator_alpha.common import configure_precision, format_mpf
from relator_alpha.presentation import build_console
from relator_alpha.scalar_qed import build_qed_induced_coefficients, dc_qed_from_alpha



def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compute the QED-induced scalar branch D_C^QED from universal A_1^(2n)."
    )
    parser.add_argument("--order", type=int, default=5, help="Truncation order M for D_C^QED,[M](x).")
    parser.add_argument(
        "--alpha",
        type=str,
        default=None,
        help="Direct α input. If omitted, α is derived from the ALP cross-check.",
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
    coefficients = build_qed_induced_coefficients(order=args.order)

    if args.alpha is None:
        alpha_result = solve_alpha_from_qed_branch(mp.mpf(args.lambda_geom), coefficients)
        alpha_value = alpha_result.alpha
        alpha_source = "internally derived from the ALP QED cross-check"
    else:
        alpha_value = mp.mpf(args.alpha)
        alpha_source = "user-supplied α"

    d_c_value = dc_qed_from_alpha(alpha_value, coefficients)

    q_table = Table(title="Universal pure-photonic QED input coefficients")
    q_table.add_column("n", justify="right")
    q_table.add_column("A₁^(2n)", justify="right", no_wrap=True)
    for n, value in coefficients.a1_coefficients.items():
        q_table.add_row(str(n), format_mpf(value, args.digits))

    coeff_table = Table(title=f"QED-induced scalar coefficients through order M = {coefficients.order}")
    coeff_table.add_column("n", justify="right")
    coeff_table.add_column("d_n^QED", justify="right", no_wrap=True)
    coeff_table.add_column("c_n^QED", justify="right", no_wrap=True)
    for n in range(1, coefficients.order + 1):
        coeff_table.add_row(
            str(n),
            format_mpf(coefficients.d_coefficients[n], args.digits),
            format_mpf(coefficients.c_coefficients[n], args.digits),
        )

    eval_table = Table(title="QED-induced scalar evaluation")
    eval_table.add_column("Quantity")
    eval_table.add_column("Value", justify="right", no_wrap=True)
    eval_table.add_column("Comment")
    eval_table.add_row("α", format_mpf(alpha_value, args.digits), alpha_source)
    eval_table.add_row("x = α/π", format_mpf(alpha_value / mp.pi, args.digits), "Sommerfeld variable")
    eval_table.add_row(f"D_C^QED,[{coefficients.order}](α/π)", format_mpf(d_c_value, args.digits), "truncated induced scalar branch")

    console.print(q_table)
    console.print()
    console.print(coeff_table)
    console.print()
    console.print(eval_table)


if __name__ == "__main__":
    main()
