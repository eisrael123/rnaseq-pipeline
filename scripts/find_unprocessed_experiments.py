#!/usr/bin/env python3
"""Survey the lab mount: which experiments hold FASTQ data, and which still need an
``input_args.json``.

    python3 scripts/find_unprocessed_experiments.py

Prints paths relative to 1_RNA_seq, grouped by what has to happen to each one.

An *experiment* is the directory that directly contains the ``cntl``/``test``
subdirectories -- that is the directory the pipeline is pointed at, and the one
``input_args.json`` has to sit in.

"Already processed" is decided by FASTQ *filename*, not by directory name: the archive's
folder names do not match the local ``Model_Experiment`` names, but the FASTQ filenames
recorded across the processed runs' metadata.tsv identify them unambiguously.
"""
from __future__ import annotations

import csv
import glob
import json
import os
import sys
from collections import defaultdict
from pathlib import Path

LAB_ROOT = Path("/Volumes/FlemingtonLabMain1/2b_Flemington_Lab_Experiments/1_RNA_seq")
PROCESSED = "/Volumes/TUNGSACore3/rnaseq_runs/Flemington_RNAseq_analysis_output/*/metadata.tsv"
FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
CONDITION_DIRS = {"cntl", "test", "control"}

# The fields metadata_form.html writes today. A file missing any of these predates a schema
# change and needs re-downloading from the form rather than hand-editing.
EXPECTED_FIELDS = {
    "reference_dir", "species_name", "investigator_name", "experiment_type", "cell_line",
    "perturbation_type", "perturbation_agent", "perturbation_target", "perturbation_dose",
    "co_treatment", "co_treatment_target", "co_treatment_arm", "facs_purified",
    "facs_gfp_promoter", "timepoint_hours", "library_selection", "library_strandedness",
    "sequencing_run_date", "notes", "series_label", "series_variance",
}


def processed_fastq_names() -> set[str]:
    names: set[str] = set()
    files = glob.glob(PROCESSED)
    if not files:
        sys.exit(f"No processed metadata.tsv found under {PROCESSED} — is the drive mounted?")
    for path in files:
        with open(path, newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                for column in ("fastq_r1", "fastq_r2"):
                    value = row.get(column, "NA")
                    if value and value != "NA":
                        names.add(Path(value).name)
    return names


def describe_json(path: Path) -> str:
    """One line on whether this input_args.json matches the current form output."""
    try:
        data = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        return f"UNREADABLE ({type(error).__name__})"
    missing = sorted(EXPECTED_FIELDS - set(data))
    extra = sorted(set(data) - EXPECTED_FIELDS)
    if not missing and not extra:
        return "current format"
    parts = []
    if missing:
        parts.append(f"missing {', '.join(missing)}")
    if extra:
        parts.append(f"unrecognized {', '.join(extra)}")
    return "NEEDS UPDATE: " + "; ".join(parts)


def main() -> int:
    if not LAB_ROOT.is_dir():
        sys.exit(f"Not reachable: {LAB_ROOT}")
    known = processed_fastq_names()
    print(f"{len(known)} FASTQ filenames already accounted for.\n", file=sys.stderr)

    # experiment dir -> {condition dir name: [fastq names]}
    experiments: dict[Path, dict[str, list[str]]] = defaultdict(dict)
    unreadable: list[str] = []
    for dirpath, dirnames, filenames in os.walk(
            LAB_ROOT, onerror=lambda e: unreadable.append(str(e))):
        dirnames[:] = [d for d in dirnames if not d.startswith((".", "_"))]
        fastqs = [f for f in filenames
                  if f.endswith(FASTQ_SUFFIXES) and not f.startswith("._")]
        if not fastqs:
            continue
        here = Path(dirpath)
        # FASTQs normally sit in a cntl/ or test/ dir; the experiment is its parent.
        if here.name.lower() in CONDITION_DIRS:
            experiments[here.parent][here.name] = fastqs
        else:
            experiments[here]["(loose)"] = fastqs

    todo_no_json, todo_with_json, done = [], [], []
    for directory in sorted(experiments):
        conditions = experiments[directory]
        all_fastqs = [f for group in conditions.values() for f in group]
        unprocessed = [f for f in all_fastqs if f not in known]
        json_path = directory / "input_args.json"
        has_json = json_path.is_file()
        record = (directory, conditions, len(all_fastqs), has_json, json_path)
        if not unprocessed:
            done.append(record)
        elif has_json:
            todo_with_json.append(record)
        else:
            todo_no_json.append(record)

    def show(record, note_json):
        directory, conditions, total, has_json, json_path = record
        rel = directory.relative_to(LAB_ROOT)
        shape = ", ".join(f"{name}={len(files)}" for name, files in sorted(conditions.items()))
        print(f"{rel}")
        print(f"    {total} FASTQ file(s)  [{shape}]")
        missing_side = CONDITION_DIRS.intersection({"cntl", "test"}) - {
            n.lower() for n in conditions}
        if missing_side:
            print(f"    *** INCOMPLETE: no {'/'.join(sorted(missing_side))} directory")
        if note_json and has_json:
            print(f"    input_args.json: {describe_json(json_path)}")

    print(f"=== NEED an input_args.json written ({len(todo_no_json)}) ===")
    for record in todo_no_json:
        show(record, note_json=False)

    print(f"\n=== HAVE an input_args.json, not yet processed ({len(todo_with_json)}) ===")
    for record in todo_with_json:
        show(record, note_json=True)

    print(f"\n=== already processed ({len(done)}) ===")
    for directory, conditions, total, has_json, json_path in done:
        marker = describe_json(json_path) if has_json else "no input_args.json on the mount"
        print(f"{directory.relative_to(LAB_ROOT)}  ({total} FASTQ)  [{marker}]")

    if unreadable:
        print(f"\n=== {len(unreadable)} directories could not be read ===")
        for message in unreadable[:20]:
            print(f"  {message}")
        print("  (the lists above may be incomplete)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
