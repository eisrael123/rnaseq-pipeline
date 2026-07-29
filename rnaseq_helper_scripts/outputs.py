#!/usr/bin/env python3
"""Canonical output layout and the only sanctioned way to write a results table.

Every table goes through :func:`write_table`, which enforces the column list and order from
``schemas.py``, renders missing values as ``NA`` (never ``0``, never empty) and refuses to write
a column the schema does not declare.
"""

from __future__ import annotations

import csv
import hashlib
import math
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from schemas import NA, Table, table as get_table

# Layout from section 2 of PIPELINE_CHANGES.md. `tables/` is the contract surface: anything the
# warehouse reads lives there and nowhere else.
SUBDIRS = (
    "tables",
    "artifacts",
    "artifacts/bigwig",
    "artifacts/bam",
    "artifacts/star_logs",
    "artifacts/junction_bed",
    "artifacts/kallisto",
    "artifacts/rmats_raw",
    "artifacts/deseq2",
    "artifacts/sleuth",
    "artifacts/gsea",
    "artifacts/qc",
    "reports",
    "logs",
)


def create_layout(results_dir: Path) -> Path:
    """Create the canonical directory tree. Idempotent, safe to call on a resumed run."""
    results_dir = Path(results_dir)
    for sub in SUBDIRS:
        (results_dir / sub).mkdir(parents=True, exist_ok=True)
    return results_dir


def tables_dir(results_dir: Path) -> Path:
    return Path(results_dir) / "tables"


def table_path(results_dir: Path, name: str) -> Path:
    return Path(results_dir) / get_table(name).relpath


def stage_log_path(results_dir: Path, stage: str) -> Path:
    log_dir = Path(results_dir) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    return log_dir / f"{stage}.log"


@dataclass
class RunContext:
    """Everything a table writer needs to stamp provenance onto its rows."""

    run_id: str
    experiment_id: str
    results_dir: Path
    genome_build: str
    # internal STAR/kallisto sample name -> warehouse sample_id
    sample_ids: dict[str, str] = field(default_factory=dict)
    test_group: str = "test"
    cntl_group: str = "cntl"

    @property
    def comparison_id(self) -> str:
        return f"{self.experiment_id}__{self.test_group}_vs_{self.cntl_group}"

    def sample_id(self, sample_name: str) -> str:
        """Map an internal sample name to its warehouse ``sample_id``.

        Raises rather than inventing an id: an unmapped sample would silently break the
        metadata.tsv foreign key the validator checks.
        """
        try:
            return self.sample_ids[sample_name]
        except KeyError:
            raise KeyError(
                f"sample {sample_name!r} is not in metadata.tsv; known samples: "
                f"{sorted(self.sample_ids)}"
            ) from None


def _clean_str(value) -> str:
    text = str(value)
    if text == "":
        return NA
    # TSV is the contract; a stray tab or newline in a free-text field would shift every
    # downstream column.
    return text.replace("\t", " ").replace("\r", " ").replace("\n", " ")


def _is_missing(value) -> bool:
    if value is None:
        return True
    if isinstance(value, float) and math.isnan(value):
        return True
    if isinstance(value, str) and value.strip() in ("", "NA", "nan", "NaN", "None", "."):
        return True
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False


def _format_scalar(value, dtype: str) -> str:
    if _is_missing(value):
        return NA
    if dtype == "int":
        try:
            return str(int(round(float(value))))
        except (TypeError, ValueError):
            return NA
    if dtype == "float":
        try:
            number = float(value)
        except (TypeError, ValueError):
            return NA
        if math.isnan(number):
            return NA
        # repr() round-trips, which matters for p-values near the float floor.
        return repr(number)
    if dtype == "bool":
        if isinstance(value, str):
            return "true" if value.strip().lower() in ("1", "true", "yes", "t") else "false"
        return "true" if bool(value) else "false"
    return _clean_str(value)


def format_frame(name: str, frame: pd.DataFrame) -> pd.DataFrame:
    """Coerce ``frame`` to exactly the schema's columns, in order, as formatted strings."""
    spec: Table = get_table(name)
    missing = [c for c in spec.column_names if c not in frame.columns]
    if missing:
        raise ValueError(f"{name}: writer did not produce required column(s) {missing}")
    extra = [c for c in frame.columns if c not in spec.column_names]
    if extra:
        raise ValueError(f"{name}: writer produced undeclared column(s) {extra}")

    out = pd.DataFrame(index=frame.index)
    for column in spec.columns:
        series = frame[column.name]
        out[column.name] = [_format_scalar(v, column.dtype) for v in series]
    return out[spec.column_names]


def write_table(results_dir: Path, name: str, frame: pd.DataFrame | None = None,
                rows: list[dict] | None = None) -> Path:
    """Write ``tables/<name>.tsv``. A zero-row table is still written, with its header."""
    spec = get_table(name)
    if frame is None:
        frame = pd.DataFrame(rows or [], columns=spec.column_names)
    if frame.empty:
        frame = pd.DataFrame(columns=spec.column_names)
    formatted = format_frame(name, frame)

    path = Path(results_dir) / spec.relpath
    path.parent.mkdir(parents=True, exist_ok=True)
    formatted.to_csv(path, sep="\t", index=False, quoting=csv.QUOTE_NONE,
                     escapechar="\\", lineterminator="\n")
    return path


def sha256_file(path: Path, chunk_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def md5_file(path: Path, chunk_size: int = 1 << 20) -> str:
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()
