#!/usr/bin/env python3
"""Task 14: every QC number in one long table.

Long format on purpose: adding a metric is a new row, never a new column, so the warehouse
schema survives the next tool upgrade.
"""

from __future__ import annotations

import json
import re
import zipfile
from pathlib import Path

import pandas as pd

import strandedness
from outputs import RunContext, write_table
from schemas import NA, QC_METRIC


def _row(ctx: RunContext, sample_id: str, tool: str, name: str,
         value=None, value_str=None) -> dict:
    return {
        "run_id": ctx.run_id,
        "sample_id": sample_id,
        "tool": tool,
        "metric_name": name,
        "metric_value": value if value is not None else NA,
        "metric_value_str": value_str if value_str is not None else NA,
    }


def _fastp_metrics(ctx: RunContext, sample_name: str, sample_id: str,
                   results_dir: Path) -> list[dict]:
    path = Path(results_dir) / "fastp" / sample_name / f"{sample_name}.json"
    if not path.is_file():
        return []
    data = json.loads(path.read_text())
    summary = data.get("summary", {})
    before = summary.get("before_filtering", {})
    after = summary.get("after_filtering", {})
    wanted = {
        "reads_before_filtering": before.get("total_reads"),
        "reads_after_filtering": after.get("total_reads"),
        "bases_before_filtering": before.get("total_bases"),
        "q30_rate": before.get("q30_rate"),
        "gc_content": before.get("gc_content"),
        "read1_mean_length": before.get("read1_mean_length"),
        "duplication_rate": data.get("duplication", {}).get("rate"),
        "adapter_trimmed_reads": data.get("adapter_cutting", {}).get("adapter_trimmed_reads"),
    }
    return [_row(ctx, sample_id, "fastp", name, value=value)
            for name, value in wanted.items() if value is not None]


_STAR_FIELDS = {
    "input_reads": "Number of input reads",
    "uniquely_mapped_reads": "Uniquely mapped reads number",
    "uniquely_mapped_pct": "Uniquely mapped reads %",
    "multi_mapped_pct": "% of reads mapped to multiple loci",
    "too_short_pct": "% of reads unmapped: too short",
    "mismatch_rate_per_base": "Mismatch rate per base, %",
    "splices_total": "Number of splices: Total",
    "splices_annotated": "Number of splices: Annotated (sjdb)",
    "average_input_read_length": "Average input read length",
    "average_mapped_length": "Average mapped length",
}


def _star_metrics(ctx: RunContext, sample_name: str, sample_id: str,
                  results_dir: Path) -> list[dict]:
    path = Path(results_dir) / "star" / sample_name / f"{sample_name}_Log.final.out"
    if not path.is_file():
        return []
    text = path.read_text()
    values: dict[str, float] = {}
    for metric, label in _STAR_FIELDS.items():
        match = re.search(re.escape(label) + r"\s*\|\s*([\d.]+)", text)
        if match:
            values[metric] = float(match.group(1))

    rows = [_row(ctx, sample_id, "star", name, value=value) for name, value in values.items()]
    total, annotated = values.get("splices_total"), values.get("splices_annotated")
    if total and annotated is not None:
        rows.append(_row(ctx, sample_id, "star", "splices_annotated_pct",
                         value=100.0 * annotated / total))
    return rows


def _rseqc_metrics(ctx: RunContext, sample_name: str, sample_id: str,
                   results_dir: Path) -> list[dict]:
    path = strandedness.rseqc_path(results_dir, sample_name)
    if not path.is_file():
        return []
    parsed = strandedness.parse_rseqc(path.read_text())
    rows = [_row(ctx, sample_id, "rseqc", "strandedness",
                 value_str=parsed["library_type"])]
    for metric in ("fraction_failed", "fraction_forward", "fraction_reverse"):
        if parsed[metric] is not None:
            rows.append(_row(ctx, sample_id, "rseqc", metric, value=parsed[metric]))
    return rows


def _kallisto_metrics(ctx: RunContext, sample_name: str, sample_id: str,
                      results_dir: Path, genome_build: str) -> list[dict]:
    path = (Path(results_dir) / "kallisto" / sample_name /
            f"{sample_name}_{genome_build}" / "run_info.json")
    if not path.is_file():
        return []
    data = json.loads(path.read_text())
    return [_row(ctx, sample_id, "kallisto", metric, value=data[metric])
            for metric in ("n_processed", "n_pseudoaligned", "p_pseudoaligned")
            if metric in data]


def _fastqc_metrics(ctx: RunContext, sample_name: str, sample_id: str,
                    results_dir: Path) -> list[dict]:
    """Per-module PASS/WARN/FAIL, read from the summary.txt inside each FastQC zip."""
    fastqc_dir = Path(results_dir) / "fastqc" / sample_name
    if not fastqc_dir.is_dir():
        return []
    rows = []
    for archive_path in sorted(fastqc_dir.glob("*_fastqc.zip")):
        read = "r2" if re.search(r"(_2|_R2)(_001)?_fastqc\.zip$", archive_path.name) else "r1"
        try:
            with zipfile.ZipFile(archive_path) as archive:
                names = [n for n in archive.namelist() if n.endswith("summary.txt")]
                if not names:
                    continue
                content = archive.read(names[0]).decode("utf-8", errors="replace")
        except (zipfile.BadZipFile, OSError):
            continue
        for line in content.splitlines():
            fields = line.split("\t")
            if len(fields) < 2:
                continue
            status, module = fields[0].strip(), fields[1].strip()
            metric = re.sub(r"[^a-z0-9]+", "_", module.lower()).strip("_")
            rows.append(_row(ctx, sample_id, "fastqc", f"{metric}.{read}", value_str=status))
    return rows


def collect_qc_metrics(ctx: RunContext, results_dir: Path, sample_names: list[str],
                       genome_build: str) -> pd.DataFrame:
    rows: list[dict] = []
    for sample_name in sample_names:
        sample_id = ctx.sample_id(sample_name)
        rows.extend(_fastp_metrics(ctx, sample_name, sample_id, results_dir))
        rows.extend(_fastqc_metrics(ctx, sample_name, sample_id, results_dir))
        rows.extend(_star_metrics(ctx, sample_name, sample_id, results_dir))
        rows.extend(_rseqc_metrics(ctx, sample_name, sample_id, results_dir))
        rows.extend(_kallisto_metrics(ctx, sample_name, sample_id, results_dir, genome_build))
    return pd.DataFrame(rows, columns=QC_METRIC.column_names)


def write_qc_metrics(ctx: RunContext, results_dir: Path, sample_names: list[str],
                     genome_build: str) -> pd.DataFrame:
    frame = collect_qc_metrics(ctx, results_dir, sample_names, genome_build)
    write_table(results_dir, "qc_metric", frame)
    return frame


def update_metadata_strandedness(results_dir: Path, library_types: dict[str, str]) -> Path:
    """Write the inferred strandedness back into ``metadata.tsv``.

    Task 12 resolves STAR's str1/str2 from this column, so it has to land in the metadata
    rather than staying inside a per-sample RSeQC text file.
    """
    path = Path(results_dir) / "metadata.tsv"
    if not path.is_file():
        raise FileNotFoundError(f"metadata.tsv missing from {results_dir}")
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if "Sample name" not in frame.columns:
        raise ValueError(f"{path} has no 'Sample name' column")
    frame["rseqc_measured_strandedness"] = frame["Sample name"].map(library_types).fillna(NA)
    frame.to_csv(path, sep="\t", index=False)
    return path
