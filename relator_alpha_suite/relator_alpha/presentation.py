from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

from rich import box
from rich.console import Console
from rich.table import Table
from rich.text import Text


@dataclass(frozen=True)
class TableColumn:
    header: str
    kind: str = "text"  # 'text' | 'number'


@dataclass(frozen=True)
class ArticleTable:
    label: str
    title: str
    columns: Sequence[TableColumn]
    rows: Sequence[Sequence[str]]
    note: str | None = None


def build_console(width: int = 220) -> Console:
    """Build a wide console to prevent terminal-side ellipsis of scientific notation."""
    return Console(width=width, highlight=False)


def build_rich_table(spec: ArticleTable) -> Table:
    table = Table(
        title=f"{spec.label} — {spec.title}",
        box=box.SIMPLE_HEAVY,
        show_lines=False,
        expand=False,
        pad_edge=True,
    )
    for column in spec.columns:
        if column.kind == "number":
            table.add_column(column.header, justify="right", no_wrap=True)
        else:
            table.add_column(column.header, justify="left", overflow="fold")
    for row in spec.rows:
        table.add_row(*row)
    return table


def print_article_tables(console: Console, tables: Iterable[ArticleTable]) -> None:
    first = True
    for spec in tables:
        if not first:
            console.print()
        first = False
        console.print(build_rich_table(spec))
        if spec.note:
            console.print(Text(spec.note))
