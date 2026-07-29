#!/usr/bin/env python3
"""Tasks 5 and 7: differential expression as statistics-only tables.

``de_gene.tsv`` deliberately contains no per-sample values. Expression lives in
``expression_gene.tsv``, keyed by ``sample_id``, which is what makes a query portable across
experiments.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from annotation import Annotation
from outputs import RunContext, write_table
from schemas import DE_GENE, NA, TRANSCRIPT_DE

DESEQ2_RESULTS = ("deseq2", "deseq2_results.csv")
SLEUTH_FILES = {
    "wald": ("sleuth", "sleuth_DETranscripts_wt.csv"),
    "lrt": ("sleuth", "sleuth_DETranscripts_lrt.csv"),
}


def write_de_gene(ctx: RunContext, results_dir: Path, annotation: Annotation) -> int:
    """Task 5. Every gene DESeq2 tested, unfiltered -- filtering is the query's job."""
    path = Path(results_dir).joinpath(*DESEQ2_RESULTS)
    if not path.is_file():
        raise FileNotFoundError(f"DESeq2 results missing: {path}")

    # write.csv() emits the row names (the gene labels) as the first, unnamed column.
    frame = pd.read_csv(path)
    label_column = frame.columns[0]
    labels = frame[label_column].astype(str)
    resolved = [annotation.resolve_gene(label) for label in labels]

    out = pd.DataFrame({
        "run_id": ctx.run_id,
        "comparison_id": ctx.comparison_id,
        "gene_id": [gene_id for gene_id, _ in resolved],
        "gene_symbol": [symbol for _, symbol in resolved],
        "base_mean": frame.get("baseMean"),
        "log2_fold_change": frame.get("log2FoldChange"),
        "lfc_se": frame.get("lfcSE"),
        "stat": frame.get("stat"),
        "pvalue": frame.get("pvalue"),
        "padj": frame.get("padj"),
    })
    missing = [c for c in ("base_mean", "log2_fold_change", "pvalue") if out[c].isna().all()]
    if missing:
        raise ValueError(f"{path} did not contain usable {missing}; columns were "
                         f"{list(frame.columns)}")
    write_table(results_dir, "de_gene", out[DE_GENE.column_names])
    return len(out)


def write_de_gene_xlsx(results_dir: Path) -> Path | None:
    """Human-facing copy of ``de_gene.tsv``. Generated *from* the TSV, never instead of it."""
    source = Path(results_dir) / DE_GENE.relpath
    if not source.is_file():
        return None
    target = Path(results_dir) / "reports" / "de_gene.xlsx"
    target.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(source, sep="\t")
    frame = frame.sort_values("padj", ascending=True, na_position="last")
    try:
        frame.to_excel(target, index=False, sheet_name="DE genes")
    except (ImportError, ValueError) as error:
        print(f"Skipping de_gene.xlsx ({error})")
        return None
    return target


def _load_sleuth(path: Path, test_type: str, ctx: RunContext,
                 annotation: Annotation) -> pd.DataFrame:
    frame = pd.read_csv(path)
    # write.csv() prefixes an unnamed row-number column for these tables.
    if frame.columns[0].startswith("Unnamed"):
        frame = frame.drop(columns=frame.columns[0])
    if "target_id" not in frame.columns:
        raise ValueError(f"{path} has no target_id column; columns were {list(frame.columns)}")

    targets = frame["target_id"].astype(str)
    mapped = annotation.annotate_targets(targets)
    symbols = frame["gene"].astype(str) if "gene" in frame.columns else mapped["gene_symbol"]

    if "test_stat" in frame.columns:
        test_stat = frame["test_stat"]
    elif {"b", "se_b"}.issubset(frame.columns):
        # Sleuth's Wald table reports the effect size and its standard error but not the
        # statistic itself; b / se_b is that statistic.
        test_stat = frame["b"] / frame["se_b"]
    else:
        test_stat = _absent(frame)

    return pd.DataFrame({
        "run_id": ctx.run_id,
        "comparison_id": ctx.comparison_id,
        "transcript_id": targets,
        "gene_id": mapped["gene_id"],
        "gene_symbol": symbols,
        "test_stat": test_stat,
        "pval": _numeric(frame, "pval"),
        "qval": _numeric(frame, "qval"),
        "b": _numeric(frame, "b"),
        "se_b": _numeric(frame, "se_b"),
        "mean_obs": _numeric(frame, "mean_obs"),
        "test_type": test_type,
    })


def _absent(frame: pd.DataFrame) -> pd.Series:
    """A float column of NaN, so that concatenating LRT and Wald keeps a numeric dtype."""
    return pd.Series(float("nan"), index=frame.index, dtype="float64")


def _numeric(frame: pd.DataFrame, name: str) -> pd.Series:
    return frame[name] if name in frame.columns else _absent(frame)


def write_transcript_de(ctx: RunContext, results_dir: Path,
                        annotation: Annotation) -> tuple[int, str]:
    """Task 7. Returns ``(row_count, status)`` where status is ``ok`` or ``skipped``."""
    parts = []
    missing = []
    for test_type, relpath in SLEUTH_FILES.items():
        path = Path(results_dir).joinpath(*relpath)
        if path.is_file():
            parts.append(_load_sleuth(path, test_type, ctx, annotation))
        else:
            missing.append(str(path))

    if not parts:
        print("WARNING: no Sleuth results found (" + "; ".join(missing) + "). Writing "
              "transcript_de.tsv with headers only.")
        write_table(results_dir, "transcript_de",
                    pd.DataFrame(columns=TRANSCRIPT_DE.column_names))
        return 0, "skipped"

    if missing:
        print(f"WARNING: Sleuth produced only some result tables; missing {missing}.")
    combined = pd.concat(parts, ignore_index=True)
    combined["gene_symbol"] = combined["gene_symbol"].fillna(NA)
    write_table(results_dir, "transcript_de", combined[TRANSCRIPT_DE.column_names])
    return len(combined), "ok" if not missing else "partial"
