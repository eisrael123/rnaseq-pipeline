#!/usr/bin/env python3
"""Render ``docs/SCHEMA.md`` from ``schemas.py``.

The document is generated rather than hand-written because a stale schema reference is worse
than none: someone will trust it. ``--check`` fails if the committed file is out of date, and
the test suite runs that check.

    python rnaseq_helper_scripts/generate_schema_doc.py [--check]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import schemas
from outputs import SUBDIRS
from provenance import MANIFEST_FIELDS, OPTIONAL_MANIFEST_FIELDS

REPO_ROOT = Path(__file__).resolve().parent.parent
DOC_PATH = REPO_ROOT / "docs" / "SCHEMA.md"

PREAMBLE = f"""# Output schema

Generated from `rnaseq_helper_scripts/schemas.py` by
`rnaseq_helper_scripts/generate_schema_doc.py`. Do not edit by hand; edit the schema and
regenerate.

Pipeline version `{schemas.PIPELINE_VERSION}`.

## Conventions

- Every table is a tab-separated file with a single header line and Unix line endings.
- **Column order is part of the contract.** `validate_outputs.py` compares names and order.
- Missing values are the literal string `{schemas.NA}`. Never an empty cell, never `0`.
- Every row carries `run_id`, which joins to `run_manifest.json`.
- Coordinates are **{schemas.COORDINATE_CONVENTION}**, matching BED. rMATS already emits
  0-based starts; STAR `SJ.out.tab` is 1-based and is converted when it is read.
- Gene IDs are unversioned Ensembl accessions (`{schemas.GENE_ID_PATTERN}`). The exceptions
  are ERCC spike-ins and viral loci, which are listed in `docs/gene_id_exceptions.txt`.
- `comparison_id` is `{{experiment_id}}__{{test_group}}_vs_{{cntl_group}}`.
"""


def _layout_section() -> str:
    lines = ["## Directory layout", "",
             "```", "<results_dir>/",
             "├── run_manifest.json",
             "├── metadata.tsv",
             "├── checksums.sha256"]
    tops = [s for s in SUBDIRS if "/" not in s]
    for i, sub in enumerate(tops):
        last_top = i == len(tops) - 1
        lines.append(f"{'└──' if last_top else '├──'} {sub}/")
        children = [s for s in SUBDIRS if s.startswith(f"{sub}/")]
        for j, child in enumerate(children):
            spine = "    " if last_top else "│   "
            branch = "└──" if j == len(children) - 1 else "├──"
            lines.append(f"{spine}{branch} {child.split('/', 1)[1]}/")
    lines += ["```", "",
              "`tables/` is the contract surface: everything the warehouse reads lives there and "
              "nowhere else. `artifacts/` holds raw tool output for re-analysis, `reports/` the "
              "human-readable XLSX, figures and HTML, `logs/` the per-stage logs.", ""]
    return "\n".join(lines)


def _manifest_section() -> str:
    lines = ["## `run_manifest.json`", "",
             "Written at the start of the run with `exit_status: running` and rewritten at the "
             "end. A run that crashes still leaves a manifest, with `exit_status: failed` and "
             "the stage that failed.", "",
             "| Field | Description |", "|---|---|"]
    lines += [f"| `{name}` | {description} |" for name, description in MANIFEST_FIELDS]
    lines += ["",
              "Added as the run progresses, so absent from a manifest written at startup:", "",
              "| Field | Description |", "|---|---|"]
    lines += [f"| `{name}` | {description} |" for name, description in OPTIONAL_MANIFEST_FIELDS]
    lines.append("")
    return "\n".join(lines)


def _table_section(spec: schemas.Table) -> str:
    lines = [f"### `{spec.relpath}`", ""]
    if spec.description:
        lines += [spec.description, ""]
    if spec.grain:
        lines += [f"**Grain:** {spec.grain}", ""]
    if spec.key_columns:
        lines += ["**Key:** " + ", ".join(f"`{c}`" for c in spec.key_columns), ""]
    lines += ["| Column | Type | Notes |", "|---|---|---|"]
    for column in spec.columns:
        lines.append(f"| `{column.name}` | {_cell(column.type_label)} | {_cell(column.notes)} |")
    lines.append("")
    return "\n".join(lines)


def _cell(text: str) -> str:
    """Enum labels separate values with `|`, which would end the table cell."""
    return text.replace("|", "\\|")


def _event_key_section() -> str:
    lines = [
        "## `event_key`", "",
        "rMATS assigns its `ID` column per invocation, so it cannot identify an event across "
        "comparisons. `event_key` is built only from values intrinsic to the locus:", "",
        "```", "<event_type>:<chr>:<strand>:<coord_1>-<coord_2>-...-<coord_n>", "```", "",
        "The coordinate slots per event type are fixed, because the five rMATS files use "
        "different column names for the same positional roles:", "",
        "| Event type | Slots, in order |", "|---|---|",
    ]
    for event_type, slots in schemas.RMATS_COORD_SLOTS.items():
        lines.append(f"| `{event_type}` | " + ", ".join(f"`{s}`" for s in slots) + " |")
    lines += ["",
              f"Unused slots up to `coord_{schemas.N_COORD_SLOTS}` are `{schemas.NA}`. Only MXE "
              f"uses all {schemas.N_COORD_SLOTS}.", ""]
    return "\n".join(lines)


def _bigwig_section() -> str:
    return "\n".join([
        "## bigWig naming", "",
        "```", "<sample_id>.<content>.<genome_build>.<strand>.bw", "```", "",
        "`content` is one of " + ", ".join(f"`{c}`" for c in schemas.BIGWIG_CONTENTS) +
        "; `strand` is one of " + ", ".join(f"`{s}`" for s in schemas.STRANDS) + ".", "",
        "STAR's `str1`/`str2` never reaches a filename: which one is the plus strand depends on "
        "the library chemistry, and the mapping is resolved once from the RSeQC call. Values are "
        "**positive on both strands**, unlike the pre-1.4 files, which stored minus-strand "
        "coverage as negative numbers. Tracks are CPM-normalized and the factor actually applied "
        "is recorded in `bigwig_manifest.tsv`; the manifest, not the filename, is authoritative.",
        "",
    ])


def render() -> str:
    parts = [PREAMBLE, "", _layout_section(), _manifest_section(),
             "## Tables", "",
             "Tables every successful run must contain:", "",
             "\n".join(f"- [`{schemas.table(n).relpath}`](#{schemas.table(n).relpath.replace('/', '').replace('.', '')})"
                       for n in schemas.REQUIRED_TABLES),
             ""]
    parts += [_table_section(schemas.table(name)) for name in schemas.REQUIRED_TABLES]
    parts += ["## Inputs", "",
              "`metadata.tsv` is produced by `metadata.py`, copied into `<results_dir>/` and "
              "updated in place once strandedness is inferred. It additionally carries the "
              "legacy columns "
              + ", ".join(f"`{c}`" for c in schemas.METADATA_LEGACY_COLUMNS) +
              ", which `rnaseq.py` still reads by name, so the validator requires these columns "
              "to be present rather than requiring an exact match.", "",
              _table_section(schemas.METADATA),
              _event_key_section(), _bigwig_section()]
    return "\n".join(parts).rstrip() + "\n"


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true",
                        help="exit non-zero if docs/SCHEMA.md is out of date")
    args = parser.parse_args(argv)

    content = render()
    if args.check:
        current = DOC_PATH.read_text() if DOC_PATH.is_file() else ""
        if current != content:
            print(f"ERROR: {DOC_PATH.relative_to(REPO_ROOT)} is out of date. Regenerate with:\n"
                  f"  python rnaseq_helper_scripts/generate_schema_doc.py", file=sys.stderr)
            return 1
        return 0

    DOC_PATH.parent.mkdir(parents=True, exist_ok=True)
    DOC_PATH.write_text(content)
    print(f"wrote {DOC_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
