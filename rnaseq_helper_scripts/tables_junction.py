#!/usr/bin/env python3
"""Task 11: splice junctions as a table instead of 476 strand-split BED files.

Coordinate convention: STAR reports the first and last base of the intron, both 1-based
inclusive. ``start`` is therefore ``col2 - 1`` and ``end`` is ``col3`` unchanged, which makes
every coordinate in this table 0-based half-open, the same as BED and the same as rMATS.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from outputs import RunContext, write_table
from schemas import JUNCTION

SJ_COLUMNS = ["chr", "intron_start_1based", "intron_end_1based", "star_strand",
              "intron_motif", "is_annotated", "n_uniquely_mapped", "n_multi_mapped",
              "max_overhang"]
STAR_STRAND_SYMBOL = {0: ".", 1: "+", 2: "-"}


def sj_path(results_dir: Path, sample_name: str) -> Path:
    return Path(results_dir) / "star" / sample_name / f"{sample_name}_SJ.out.tab"


def read_junctions(ctx: RunContext, results_dir: Path, sample_name: str) -> pd.DataFrame:
    path = sj_path(results_dir, sample_name)
    if not path.is_file():
        raise FileNotFoundError(f"STAR SJ.out.tab missing for {sample_name}: {path}")
    frame = pd.read_csv(path, sep="\t", header=None, names=SJ_COLUMNS)
    return pd.DataFrame({
        "run_id": ctx.run_id,
        "sample_id": ctx.sample_id(sample_name),
        "chr": frame["chr"].astype(str),
        "start": frame["intron_start_1based"].astype("int64") - 1,
        "end": frame["intron_end_1based"].astype("int64"),
        "strand": frame["star_strand"].map(STAR_STRAND_SYMBOL).fillna("."),
        "intron_motif": frame["intron_motif"],
        "is_annotated": frame["is_annotated"].astype(bool),
        "n_uniquely_mapped": frame["n_uniquely_mapped"],
        "n_multi_mapped": frame["n_multi_mapped"],
        "max_overhang": frame["max_overhang"],
    })


def write_junctions(ctx: RunContext, results_dir: Path,
                    sample_names: list[str]) -> dict[str, int]:
    """Write ``junction.tsv``; returns the per-sample row counts for the acceptance check."""
    per_sample = {}
    frames = []
    for sample_name in sample_names:
        frame = read_junctions(ctx, results_dir, sample_name)
        per_sample[sample_name] = len(frame)
        frames.append(frame)

    combined = (pd.concat(frames, ignore_index=True) if frames
                else pd.DataFrame(columns=JUNCTION.column_names))
    write_table(results_dir, "junction", combined[JUNCTION.column_names])
    return per_sample
