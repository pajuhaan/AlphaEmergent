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
from relator_alpha.geometry import build_lambda_geom_representative
from relator_alpha.presentation import build_console



def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compute the vector-channel representative Λ_geom."
    )
    parser.add_argument("--p-ir-chi", type=str, default=None, help="Optional custom P_IR^(χ)(ℓ₀).")
    parser.add_argument("--r-chi", type=str, default=None, help="Optional custom shell-memory factor R_χ.")
    parser.add_argument("--delta-lambda-out", type=str, default=None, help="Optional custom ΔΛ_out(η₀).")
    parser.add_argument("--digits", type=int, default=24, help="Displayed significant digits.")
    parser.add_argument("--dps", type=int, default=160, help="Working precision in decimal digits.")
    parser.add_argument("--width", type=int, default=180, help="Console width used for Rich tables.")
    return parser



def main() -> None:
    args = build_parser().parse_args()
    configure_precision(args.dps)
    console = build_console(args.width)

    result = build_lambda_geom_representative(
        p_ir_chi=mp.mpf(args.p_ir_chi) if args.p_ir_chi is not None else None,
        r_chi=mp.mpf(args.r_chi) if args.r_chi is not None else None,
        delta_lambda_out=mp.mpf(args.delta_lambda_out) if args.delta_lambda_out is not None else None,
    )

    table = Table(title="Vector-channel geometric representative")
    table.add_column("Symbol")
    table.add_column("Value", justify="right", no_wrap=True)
    table.add_column("Comment")
    table.add_row("η₀", format_mpf(result.eta_0, args.digits), "minimal branch ratio")
    table.add_row("ℓ₀", format_mpf(result.ell_0, args.digits), "minimal shell width")
    table.add_row("Λ_ind", format_mpf(result.lambda_ind, args.digits), "filamentary reference core")
    table.add_row("C_UV^(Gauss)", format_mpf(result.c_uv_gauss, args.digits), "Gaussian UV constant")
    table.add_row("P_IR^(χ)(ℓ₀)", format_mpf(result.p_ir_chi, args.digits), "shell-admission quotient")
    table.add_row("ΔΛ_UV→IR(ℓ₀)", format_mpf(result.delta_lambda_uv_to_ir, args.digits), "admitted UV load")
    table.add_row("ΔΛ_out(η₀)", format_mpf(result.delta_lambda_out, args.digits), "exterior subtraction")
    table.add_row("R_χ", format_mpf(result.r_chi, args.digits), "shell-memory factor")
    table.add_row("ΔΛ_out(η₀) R_χ", format_mpf(result.exterior_memory_contribution, args.digits), "exterior-memory contribution")
    table.add_row("Λ_geom", format_mpf(result.lambda_geom, args.digits), result.status)
    console.print(table)
    console.print(
        "Note: the default output is the article's diagnostic representative "
        "Λ_geom^(mean), not the unresolved operator-exact Λ_geom."
    )


if __name__ == "__main__":
    main()
