#!/usr/bin/env python3
"""Recompute signal_over_gene.tsv (+ strand check) for an existing results directory.

Use after fixing the annotation BED / bigwig.py dedupe logic, without re-running STAR.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "rnaseq_helper_scripts"))

import bigwig  # noqa: E402
import checksums  # noqa: E402
import validate_outputs  # noqa: E402
from annotation import Annotation  # noqa: E402
from outputs import RunContext  # noqa: E402


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print(f"Usage: {argv[0]} <results_dir> <reference_dir>", file=sys.stderr)
        return 2

    results = Path(argv[1]).resolve()
    ref = Path(argv[2]).resolve()
    manifest_path = results / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    meta = pd.read_csv(results / "metadata.tsv", sep="\t", keep_default_na=False)
    bw_man = pd.read_csv(results / "tables" / "bigwig_manifest.tsv", sep="\t",
                         keep_default_na=False)

    genome_build = manifest["genome_build"]
    ctx = RunContext(
        run_id=manifest["run_id"],
        experiment_id=manifest["experiment_id"],
        results_dir=results,
        genome_build=genome_build,
        sample_ids={row["Sample name"]: row["sample_id"] for _, row in meta.iterrows()},
    )

    mart = ref / genome_build / "biomart" / f"{genome_build}.mart_export.txt"
    if not mart.is_file():
        raise FileNotFoundError(mart)
    annotation = Annotation(mart)
    print(f"gene_id_source={annotation.source}")

    rows, status = bigwig.write_signal_over_gene(
        ctx, results, ref, annotation, bw_man.to_dict(orient="records"))
    print(f"signal_over_gene rows={rows} status={status}")

    library_type = meta["strandedness"].iloc[0]
    organism = meta["organism"].iloc[0]
    strand_status = bigwig.assert_strand_assignment(
        results, organism, library_type, annotation)
    print(f"strand_check={strand_status}")

    for junk in results.rglob(".DS_Store"):
        junk.unlink(missing_ok=True)

    known = {row["file_path"]: row["sha256"] for _, row in bw_man.iterrows()}
    checksums.write_checksums(results, known=known)

    manifest.setdefault("stage_status", {})
    manifest["stage_status"]["signal_over_gene_rows"] = rows
    manifest["stage_status"]["signal_over_gene"] = status
    manifest["stage_status"]["strand_check"] = strand_status

    report = validate_outputs.validate(results)
    manifest["validation"] = "pass" if report.ok else "fail"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    for failure in report.failures:
        print(f"FAIL: {failure}")
    for warning in report.warnings[:12]:
        print(f"WARN: {warning}")
    print(f"validation={manifest['validation']} "
          f"({len(report.failures)} failure(s), {len(report.warnings)} warning(s))")
    return 0 if report.ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
