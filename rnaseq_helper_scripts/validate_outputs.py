#!/usr/bin/env python3
"""Task 16: validate a results directory against the output contract.

    python validate_outputs.py <results_dir>

Exits non-zero with a numbered list of failures. Called automatically at the end of
``rnaseq.py``, which records the outcome in the manifest as ``"validation": "pass" | "fail"``.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import pandas as pd

import checksums
from schemas import (GENE_ID_PATTERN, METADATA, REQUIRED_TABLES, VIRAL_GENE_ID_PREFIXES,
                     table as get_table)

GENE_ID_RE = re.compile(GENE_ID_PATTERN)
EXCEPTIONS_FILE = Path(__file__).resolve().parent.parent / "docs" / "gene_id_exceptions.txt"
GENE_ID_TABLES = ("de_gene", "expression_gene", "expression_transcript", "transcript_de",
                  "splicing_event", "signal_over_gene")


class Report:
    def __init__(self) -> None:
        self.failures: list[str] = []
        self.warnings: list[str] = []

    def fail(self, message: str) -> None:
        self.failures.append(message)

    def warn(self, message: str) -> None:
        self.warnings.append(message)

    @property
    def ok(self) -> bool:
        return not self.failures


def _find_nulls(value, path: str = "") -> list[str]:
    if value is None:
        return [path or "(root)"]
    found = []
    if isinstance(value, dict):
        for key, child in value.items():
            found.extend(_find_nulls(child, f"{path}.{key}" if path else key))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            found.extend(_find_nulls(child, f"{path}[{index}]"))
    return found


def load_gene_id_exceptions() -> tuple[set[str], tuple[str, ...]]:
    """Documented non-Ensembl gene ids from ``docs/gene_id_exceptions.txt``.

    One entry per line, ``#`` starts a comment, a trailing ``*`` makes the entry a prefix.
    Returns ``(exact, prefixes)``.
    """
    exact: set[str] = set()
    prefixes: list[str] = list(VIRAL_GENE_ID_PREFIXES)
    if EXCEPTIONS_FILE.is_file():
        for line in EXCEPTIONS_FILE.read_text().splitlines():
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            if line.endswith("*"):
                prefixes.append(line[:-1])
            else:
                exact.add(line)
    return exact, tuple(prefixes)


def _is_acceptable_gene_id(gene_id: str, exact: set[str], prefixes: tuple[str, ...]) -> bool:
    return bool(GENE_ID_RE.match(gene_id)) or gene_id in exact or gene_id.startswith(prefixes)


def _read_columns(path: Path, columns: list[str]) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", usecols=columns, dtype=str, keep_default_na=False)


def validate(results_dir: Path) -> Report:
    results_dir = Path(results_dir)
    report = Report()

    # 1. Manifest.
    manifest_path = results_dir / "run_manifest.json"
    manifest = {}
    if not manifest_path.is_file():
        report.fail("run_manifest.json is missing")
    else:
        try:
            manifest = json.loads(manifest_path.read_text())
        except json.JSONDecodeError as error:
            report.fail(f"run_manifest.json is not valid JSON: {error}")
        nulls = _find_nulls(manifest)
        if nulls:
            report.fail(f"run_manifest.json has null value(s) at: {', '.join(sorted(nulls))}")
        if manifest.get("exit_status") != "success":
            report.fail(f"run_manifest.json exit_status is {manifest.get('exit_status')!r}, "
                        "expected 'success'")

    manifest_run_id = manifest.get("run_id")

    # 2 + 3. Every table exists with exactly the declared columns, in order.
    present: dict[str, Path] = {}
    for name in REQUIRED_TABLES:
        spec = get_table(name)
        path = results_dir / spec.relpath
        if not path.is_file():
            report.fail(f"{spec.relpath} is missing")
            continue
        header = path.read_text().split("\n", 1)[0].rstrip("\r")
        if not header:
            report.fail(f"{spec.relpath} has no header row")
            continue
        found = header.split("\t")
        if found != spec.column_names:
            report.fail(
                f"{spec.relpath} columns do not match the schema\n"
                f"      expected: {spec.column_names}\n"
                f"      found:    {found}"
            )
            continue
        present[name] = path

    # 5. metadata.tsv is the sample_id authority.
    metadata_path = results_dir / "metadata.tsv"
    known_sample_ids: set[str] = set()
    if not metadata_path.is_file():
        report.fail("metadata.tsv is missing")
    else:
        metadata = pd.read_csv(metadata_path, sep="\t", dtype=str, keep_default_na=False)
        missing_columns = [c for c in METADATA.column_names if c not in metadata.columns]
        if missing_columns:
            report.fail(f"metadata.tsv is missing column(s) {missing_columns}")
        if "sample_id" in metadata.columns:
            known_sample_ids = set(metadata["sample_id"])

    exact_exceptions, exception_prefixes = load_gene_id_exceptions()
    gene_id_source = manifest.get("gene_id_source", "ensembl")

    for name, path in present.items():
        spec = get_table(name)

        # 4. run_id present, single-valued, matching the manifest.
        run_ids = set(_read_columns(path, ["run_id"])["run_id"])
        if len(run_ids) > 1:
            report.fail(f"{spec.relpath} contains {len(run_ids)} distinct run_id values")
        elif run_ids and manifest_run_id and run_ids != {manifest_run_id}:
            report.fail(f"{spec.relpath} run_id {run_ids.pop()!r} does not match the manifest "
                        f"run_id {manifest_run_id!r}")

        # 5. sample_id foreign key.
        if "sample_id" in spec.column_names and known_sample_ids:
            found = set(_read_columns(path, ["sample_id"])["sample_id"])
            unknown = sorted(found - known_sample_ids)
            if unknown:
                report.fail(f"{spec.relpath} references sample_id(s) not in metadata.tsv: "
                            f"{unknown[:5]}{'...' if len(unknown) > 5 else ''}")

        # 6. gene_id shape.
        if name in GENE_ID_TABLES:
            gene_ids = set(_read_columns(path, ["gene_id"])["gene_id"]) - {"NA", ""}
            bad = sorted(g for g in gene_ids
                         if not _is_acceptable_gene_id(g, exact_exceptions, exception_prefixes))
            if bad:
                message = (f"{spec.relpath} has {len(bad)} gene_id value(s) that are neither "
                           f"Ensembl nor documented exceptions, e.g. {bad[:5]}")
                if gene_id_source == "ensembl":
                    report.fail(message)
                else:
                    # The reference biomart export has no Ensembl gene column, so gene_id holds
                    # symbols. Open question 3 in PIPELINE_CHANGES.md has to be resolved before
                    # this can be an error.
                    report.warn(message + f"\n      (manifest gene_id_source="
                                          f"{gene_id_source!r}; see docs/SCHEMA.md)")

    # 7. Replicate event keys are a subset of the event keys.
    if "splicing_event" in present and "splicing_event_replicate" in present:
        events = set(_read_columns(present["splicing_event"], ["event_key"])["event_key"])
        replicates = set(
            _read_columns(present["splicing_event_replicate"], ["event_key"])["event_key"])
        orphans = sorted(replicates - events)
        if orphans:
            report.fail(f"splicing_event_replicate.tsv has {len(orphans)} event_key value(s) "
                        f"absent from splicing_event.tsv, e.g. {orphans[:3]}")

    # 8. Checksums.
    for problem in checksums.verify_checksums(results_dir):
        report.fail(f"checksums: {problem}")

    # 9. Housekeeping.
    for path in sorted(results_dir.rglob("*")):
        relative = path.relative_to(results_dir)
        if path.name == ".DS_Store":
            report.fail(f"{relative}: .DS_Store must not be shipped")
        elif path.is_file() and path.stat().st_size == 0:
            # An empty log means the stage had nothing to say, which is normal. An empty
            # table or manifest means a stage produced nothing, which is not.
            if relative.parts[0] == "logs" or path.suffix in (".log", ".err"):
                report.warn(f"{relative}: zero-byte file")
            else:
                report.fail(f"{relative}: zero-byte file")

    return report


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        description="Check a results directory against the output contract in docs/SCHEMA.md.")
    parser.add_argument("results_dir", type=Path)
    args = parser.parse_args(argv)

    results_dir = args.results_dir
    if not results_dir.is_dir():
        parser.error(f"no such directory: {results_dir}")

    report = validate(results_dir)
    for index, warning in enumerate(report.warnings, start=1):
        print(f"warning {index}: {warning}")
    if report.ok:
        print(f"PASS: {results_dir} satisfies the output contract "
              f"({len(report.warnings)} warning(s)).")
        return 0
    print(f"FAIL: {len(report.failures)} problem(s) in {results_dir}")
    for index, failure in enumerate(report.failures, start=1):
        print(f"  {index}. {failure}")
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
