#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from relator_alpha.article_tables import all_article_tables
from relator_alpha.common import configure_precision
from relator_alpha.presentation import build_console, print_article_tables



def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run the full Relator-alpha pipeline and print the numerical tables "
            "of the article in publication order, including intermediate values."
        )
    )
    parser.add_argument("--digits", type=int, default=24, help="Displayed significant digits.")
    parser.add_argument("--dps", type=int, default=160, help="Working precision in decimal digits.")
    parser.add_argument("--width", type=int, default=300, help="Console width used for Rich tables.")
    return parser



def main() -> None:
    args = build_parser().parse_args()
    configure_precision(args.dps)
    console = build_console(args.width)

    console.print("[bold]Relator Alpha Suite[/bold]")
    console.print(
        "Default pipeline mode: compute and print the article-style numerical tables, "
        "including intermediate values such as K, ALP-lock pieces, scalar N(D)/Q(D), "
        "QED-induced coefficients, bridge terms, and appendix sensitivities."
    )
    console.print()
    print_article_tables(console, all_article_tables(digits=args.digits))
    console.print()
    console.print(
        "[bold green]Done.[/bold green] "
        "All default article tables and their intermediate diagnostics were evaluated."
    )


if __name__ == "__main__":
    main()
