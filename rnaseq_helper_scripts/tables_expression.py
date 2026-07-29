#!/usr/bin/env python3
"""Task 6: kallisto abundances in long format.

Transcript rows come straight from ``abundance.tsv``; gene rows are summarized with the same
tx2gene mapping DESeq2 consumes, so ``expression_gene`` and ``de_gene`` cannot disagree about
what a gene is.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from annotation import Annotation
from outputs import RunContext, write_table
from schemas import EXPRESSION_GENE, EXPRESSION_TRANSCRIPT, NA


def abundance_path(results_dir: Path, sample_name: str, genome_build: str) -> Path:
    return (Path(results_dir) / "kallisto" / sample_name /
            f"{sample_name}_{genome_build}" / "abundance.tsv")


def _load_sample(ctx: RunContext, results_dir: Path, sample_name: str, genome_build: str,
                 annotation: Annotation) -> pd.DataFrame:
    path = abundance_path(results_dir, sample_name, genome_build)
    if not path.is_file():
        raise FileNotFoundError(
            f"kallisto abundance.tsv missing for {sample_name}: {path}"
        )
    frame = pd.read_csv(path, sep="\t")
    mapped = annotation.annotate_targets(frame["target_id"].astype(str))
    return pd.DataFrame({
        "run_id": ctx.run_id,
        "sample_id": ctx.sample_id(sample_name),
        "transcript_id": frame["target_id"].astype(str),
        "gene_id": mapped["gene_id"],
        "gene_symbol": mapped["gene_symbol"],
        "length": frame["length"],
        "eff_length": frame["eff_length"],
        "est_counts": frame["est_counts"],
        "tpm": frame["tpm"],
    })


def write_expression_tables(ctx: RunContext, results_dir: Path, sample_names: list[str],
                            genome_build: str, annotation: Annotation,
                            write_transcripts: bool = True) -> dict[str, int]:
    """Write ``expression_transcript.tsv`` and ``expression_gene.tsv``; return row counts."""
    per_sample = [
        _load_sample(ctx, results_dir, sample_name, genome_build, annotation)
        for sample_name in sample_names
    ]
    combined = pd.concat(per_sample, ignore_index=True)

    if write_transcripts:
        transcripts = combined[EXPRESSION_TRANSCRIPT.column_names]
    else:
        # Kept behind a flag because this table is by far the largest; the header-only file
        # keeps the contract surface complete for the loader.
        transcripts = pd.DataFrame(columns=EXPRESSION_TRANSCRIPT.column_names)
    write_table(results_dir, "expression_transcript", transcripts)

    mapped = combined[combined["gene_id"] != NA]
    genes = (
        mapped.groupby(["run_id", "sample_id", "gene_id", "gene_symbol"], as_index=False)
        .agg(est_counts=("est_counts", "sum"), tpm=("tpm", "sum"))
        .sort_values(["sample_id", "gene_id"], kind="stable")
    )
    write_table(results_dir, "expression_gene", genes[EXPRESSION_GENE.column_names])

    unmapped = int((combined["gene_id"] == NA).sum())
    if unmapped:
        print(f"expression_gene: {unmapped} transcript rows had no gene mapping and were "
              f"excluded from gene-level summarization (they remain in expression_transcript).")

    return {
        "expression_transcript": int(len(transcripts)),
        "expression_gene": int(len(genes)),
        "unmapped_transcript_rows": unmapped,
    }
