#!/usr/bin/env python3
"""Task 8: GSEA results for every MSigDB collection in one table.

GSEA writes a directory per collection per phenotype and buries the statistics in files whose
names embed a timestamp. One row per gene set, with ``collection`` distinguishing them, is the
only shape that can answer "is this pathway up in any of our experiments".
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from outputs import RunContext, write_table
from schemas import GENE_SET_ENRICHMENT, NA

REPORT_COLUMNS = {
    "gene_set": ("NAME",),
    "size": ("SIZE",),
    "es": ("ES",),
    "nes": ("NES",),
    "pvalue": ("NOM p-val", "NOM.p.val", "NOM p-value"),
    "fdr_qvalue": ("FDR q-val", "FDR.q.val"),
    "fwer_pvalue": ("FWER p-val", "FWER.p.val"),
    "rank_at_max": ("RANK AT MAX",),
}
SYMBOL_COLUMNS = ("SYMBOL", "GENE SYMBOL", "PROBE")
_VERSION_SUFFIX = re.compile(r"\.v\d.*$")


def collection_name(report_dir_name: str) -> str:
    """``preranked_h.all.v2023.2.Hs.symbols.GseaPreranked.170`` -> ``h.all``."""
    name = re.sub(r"\.(Gsea|GseaPreranked)\.\d+$", "", report_dir_name)
    name = re.sub(r"^preranked_", "", name)
    return _VERSION_SUFFIX.sub("", name) or name


def _pick(frame: pd.DataFrame, candidates: tuple[str, ...]) -> pd.Series | None:
    for candidate in candidates:
        if candidate in frame.columns:
            return frame[candidate]
    return None


def _leading_edge(report_dir: Path, gene_set: str) -> str:
    """Core-enrichment symbols from the per-set detail file GSEA writes for plotted sets."""
    for suffix in (".tsv", ".xls"):
        detail = report_dir / f"{gene_set}{suffix}"
        if not detail.is_file():
            continue
        try:
            frame = pd.read_csv(detail, sep="\t", dtype=str)
        except (OSError, ValueError):
            return NA
        core = frame
        for column in frame.columns:
            if column.strip().upper().startswith("CORE ENRICHMENT"):
                core = frame[frame[column].astype(str).str.strip().str.lower() == "yes"]
                break
        symbols = _pick(core, SYMBOL_COLUMNS)
        if symbols is None:
            return NA
        values = [s for s in symbols.astype(str) if s and s.lower() != "nan"]
        return ",".join(values) if values else NA
    return NA


def write_gene_set_enrichment(ctx: RunContext, results_dir: Path) -> tuple[int, str]:
    """Returns ``(row_count, status)``; status is ``ok``, ``partial`` or ``skipped``."""
    gsea_root = Path(results_dir) / "gsea"
    report_dirs = sorted(
        [d for d in gsea_root.glob("*.Gsea.*") if d.is_dir()] +
        [d for d in gsea_root.glob("*.GseaPreranked.*") if d.is_dir()]
    ) if gsea_root.is_dir() else []

    rows: list[pd.DataFrame] = []
    empty_dirs: list[str] = []
    for report_dir in report_dirs:
        collection = collection_name(report_dir.name)
        reports = sorted(
            p for p in report_dir.glob("gsea_report_for_*")
            if p.suffix in (".tsv", ".xls")
        )
        if not reports:
            empty_dirs.append(report_dir.name)
            continue
        for report in reports:
            frame = pd.read_csv(report, sep="\t")
            names = _pick(frame, REPORT_COLUMNS["gene_set"])
            if names is None:
                empty_dirs.append(report.name)
                continue
            frame = frame[names.notna()]
            names = names[names.notna()].astype(str)
            out = pd.DataFrame({
                "run_id": ctx.run_id,
                "comparison_id": ctx.comparison_id,
                "collection": collection,
                "gene_set": names,
                "leading_edge_genes": [_leading_edge(report_dir, name) for name in names],
            })
            for field in ("size", "es", "nes", "pvalue", "fdr_qvalue", "fwer_pvalue",
                          "rank_at_max"):
                values = _pick(frame, REPORT_COLUMNS[field])
                out[field] = values.values if values is not None else None
            rows.append(out)

    if not rows:
        print("WARNING: no GSEA reports found under "
              f"{gsea_root}. Writing gene_set_enrichment.tsv with headers only.")
        write_table(results_dir, "gene_set_enrichment",
                    pd.DataFrame(columns=GENE_SET_ENRICHMENT.column_names))
        return 0, "skipped"

    combined = pd.concat(rows, ignore_index=True)
    combined = combined.sort_values(["collection", "gene_set"], kind="stable")
    write_table(results_dir, "gene_set_enrichment",
                combined[GENE_SET_ENRICHMENT.column_names])
    if empty_dirs:
        print(f"WARNING: GSEA produced no report table in: {', '.join(empty_dirs)}")
    return len(combined), "partial" if empty_dirs else "ok"
