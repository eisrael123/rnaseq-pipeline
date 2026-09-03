#!/usr/bin/env python3
"""Check every input_args.json under a root for correctness before anything is run.

    python3 scripts/validate_input_args.py [root ...]

Default roots are the lab mount and the local copy. Checks each file on its own, then
checks series members against each other, and finally cross-references what the file
claims against what the directory path and the FASTQ layout say.

Exits non-zero if anything is wrong.
"""
from __future__ import annotations

import json
import re
import sys
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "rnaseq_helper_scripts"))
sys.path.insert(0, str(REPO))
import vocab  # noqa: E402

# metadata.py guards against being *run* outside the pipeline env; we only import it for its
# FASTQ-naming patterns, so satisfy the guard rather than duplicating the patterns here. A
# local copy is exactly how this check drifted out of agreement with the pipeline once already.
import os  # noqa: E402
os.environ.setdefault("CONDA_DEFAULT_ENV", "rnaseqpipeline")
import metadata  # noqa: E402

DEFAULT_ROOTS = [
    Path("/Volumes/FlemingtonLabMain1/2b_Flemington_Lab_Experiments/1_RNA_seq"),
    Path("/Volumes/TUNGSACore3/rnaseq_runs"),
]
REQUIRED = {
    "reference_dir", "species_name", "investigator_name", "experiment_type", "cell_line",
    "perturbation_type", "perturbation_agent", "perturbation_target", "perturbation_dose",
    "co_treatment", "co_treatment_target", "co_treatment_arm", "facs_purified",
    "facs_gfp_promoter", "timepoint_hours", "library_selection", "library_strandedness",
    "sequencing_run_date", "notes", "series_label", "series_variance",
}
# field -> vocabulary, for the ones that are closed sets. NA is allowed where the pipeline
# writes NA itself.
CLOSED = {
    "cell_line": vocab.CELL_LINES,
    "perturbation_type": vocab.PERTURBATION_TYPES,
    "co_treatment": vocab.CO_TREATMENTS,
    "co_treatment_target": vocab.CO_TREATMENT_TARGETS | {"NA"},
    "co_treatment_arm": vocab.CO_TREATMENT_ARMS,
    "facs_purified": vocab.FACS_PURIFIED_VALUES,
    "facs_gfp_promoter": vocab.FACS_GFP_PROMOTERS,
    "library_selection": vocab.LIBRARY_SELECTIONS,
    "library_strandedness": vocab.LIBRARY_STRANDEDNESS_VALUES,
    "series_variance": vocab.SERIES_VARIANCES,
    "experiment_type": vocab.LIBRARY_LAYOUTS,
    "species_name": set(vocab.GENOME_BUILD_ORGANISM),
}
# Hints the directory path gives, which the file should agree with.
PATH_HINTS = [
    (r"(?<![A-Za-z])(\d+)\s*hr(?![A-Za-z])", "timepoint_hours", lambda m: str(float(m.group(1)))),
    (r"(?<![A-Za-z0-9])(\d+)\s*nM(?![A-Za-z])", "perturbation_dose", lambda m: m.group(0)),
    (r"(?<![A-Za-z0-9])(\d+)\s*uM(?![A-Za-z])", "perturbation_dose", lambda m: m.group(0)),
    (r"ribodeplet", "library_selection", lambda m: "ribodepleted"),
    (r"polyA", "library_selection", lambda m: "polyA"),
    (r"stranded", "library_strandedness", lambda m: "stranded"),
]
FASTQ = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


def check_one(path: Path, root: Path) -> tuple[dict | None, list[str]]:
    rel = path.parent.relative_to(root)
    problems: list[str] = []
    try:
        data = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        return None, [f"unreadable/invalid JSON: {error}"]
    if not isinstance(data, dict):
        return None, ["top level is not a JSON object"]

    missing = sorted(REQUIRED - set(data))
    extra = sorted(set(data) - REQUIRED)
    if missing:
        problems.append(f"missing field(s): {', '.join(missing)}")
    if extra:
        problems.append(f"unrecognized field(s): {', '.join(extra)} "
                        f"(metadata.py will refuse the file)")

    for field, allowed in CLOSED.items():
        value = data.get(field)
        if value is not None and value not in allowed:
            problems.append(f"{field}={value!r} is not in the allowed set "
                            f"({', '.join(sorted(allowed))})")

    # --- the path is evidence; the file should not contradict it ---
    text = str(rel)
    # First hint wins per field. Some archive folders name both alternatives -- e.g.
    # "2_Mutu_polyA_Zta_24hrs_ribodepletion" contains both "polyA" and "ribodepletion" -- and
    # the more specific pattern is listed first, so later ones must not re-flag the same field.
    hinted: set[str] = set()
    for pattern, field, extract in PATH_HINTS:
        match = re.search(pattern, text, re.I)
        if not match or field not in data or field in hinted:
            continue
        hinted.add(field)
        expected, actual = extract(match), str(data[field])
        norm = lambda s: s.replace(" ", "").lower()
        if field == "timepoint_hours":
            try:
                same = float(actual) == float(expected)
            except ValueError:
                same = False
        else:
            same = norm(actual) == norm(expected)
        if not same:
            problems.append(f"path says {field}~{expected!r} but file says {actual!r}")

    # cell_line must appear in the path (the archive groups by cell line at the top level)
    cell = data.get("cell_line")
    if cell and cell.lower() not in text.lower():
        problems.append(f"cell_line={cell!r} does not appear anywhere in the path {text!r}")

    # genome build vs. cell line: KSHV lines need a KSHV build, EBV lines an EBV build
    build = str(data.get("species_name", ""))
    if cell == "BCBL1" and "KSHV" not in build:
        problems.append(f"BCBL1 is a KSHV line but species_name={build!r}")
    if cell in {"Akata", "Mutu", "P3HR1", "Raji", "SNU719", "YCCEL1"} and "KSHV" in build:
        problems.append(f"{cell} is an EBV line but species_name={build!r}")

    # --- series rules, mirroring metadata.py ---
    label, variance = data.get("series_label", "NA"), data.get("series_variance", "NA")
    co = data.get("co_treatment", "none")
    if (label == "NA") != (variance == "NA"):
        problems.append(f"series_label={label!r} and series_variance={variance!r} are mutually "
                        f"inclusive: both real, or both NA")
    if co != "none":
        if variance != "co_treatment":
            problems.append(f"co_treatment={co!r} forces series_variance='co_treatment', "
                            f"got {variance!r}")
        if label == "NA":
            problems.append(f"co_treatment={co!r} requires a series_label")
        if data.get("co_treatment_arm") not in {"test", "cntl"}:
            problems.append(f"co_treatment={co!r} requires co_treatment_arm test|cntl, "
                            f"got {data.get('co_treatment_arm')!r}")
        if data.get("co_treatment_target", "NA") == "NA":
            problems.append(f"co_treatment={co!r} but co_treatment_target is NA")
    else:
        if variance == "co_treatment":
            problems.append("series_variance='co_treatment' but co_treatment='none'")
        for field in ("co_treatment_target", "co_treatment_arm"):
            if data.get(field, "NA") != "NA":
                problems.append(f"co_treatment='none' so {field} must be NA, "
                                f"got {data.get(field)!r}")

    # --- the FASTQ layout the pipeline requires ---
    subdirs = {d.name.lower() for d in path.parent.iterdir() if d.is_dir()}
    for needed in ("cntl", "test"):
        if needed not in subdirs:
            problems.append(f"no {needed}/ subdirectory beside this file")
    layout = data.get("experiment_type")
    fastqs = [f for d in path.parent.iterdir() if d.is_dir()
              for f in d.iterdir() if f.name.endswith(FASTQ)]
    if fastqs and layout in {"PE", "SE"}:
        looks_paired = any(pattern.match(f.name)
                           for f in fastqs for pattern in metadata.READ2_PATTERNS)
        if layout == "PE" and not looks_paired:
            problems.append("experiment_type=PE but no R2/_2 files found")
        if layout == "SE" and looks_paired:
            problems.append("experiment_type=SE but R2/_2 files are present")

    # the directory name becomes experiment_id, so it must parse as Model_Experiment
    name = path.parent.name
    model, _, descriptor = name.partition("_")
    if not model or not descriptor:
        problems.append(f"directory name {name!r} is not 'Model_Experiment' -- it becomes "
                        f"experiment_id and every sample_id is built from it")
    elif model not in vocab.CELL_LINES:
        problems.append(f"directory name {name!r} starts with {model!r}, which is not a known "
                        f"cell line; experiment_id would disagree with cell_line={cell!r}")
    return data, problems


def main(argv: list[str]) -> int:
    roots = [Path(a) for a in argv] or DEFAULT_ROOTS
    found: list[tuple[Path, Path, dict]] = []
    total_problems = 0

    for root in roots:
        if not root.is_dir():
            print(f"(skipping unreachable root {root})\n")
            continue
        try:
            paths = sorted(root.rglob("input_args.json"))
        except PermissionError as error:
            print(f"(cannot read {root}: {error})\n")
            continue
        print(f"### {root}  —  {len(paths)} file(s)\n")
        for path in paths:
            data, problems = check_one(path, root)
            rel = path.parent.relative_to(root)
            if problems:
                total_problems += len(problems)
                print(f"  {rel}")
                for problem in problems:
                    print(f"      ✗ {problem}")
            else:
                print(f"  {rel}   ok")
            if data:
                found.append((root, path, data))
        print()

    # --- cross-file: series members must agree with each other ---
    series: dict[str, list[tuple[Path, dict]]] = defaultdict(list)
    for _, path, data in found:
        label = data.get("series_label", "NA")
        if label != "NA":
            series[label].append((path, data))

    print("### series membership\n")
    if not series:
        print("  (no series declared)")
    for label, members in sorted(series.items()):
        variances = {d.get("series_variance") for _, d in members}
        print(f"  {label}  —  {len(members)} member(s), variance={', '.join(sorted(map(str, variances)))}")
        for path, _ in members:
            print(f"      {path.parent}")
        if len(members) < 2:
            print("      ✗ only one member: a series needs at least two "
                  "(typo in the label, or the partner has no input_args.json yet)")
            total_problems += 1
        if len(variances) > 1:
            print(f"      ✗ members disagree about series_variance: {sorted(map(str, variances))}")
            total_problems += 1
        axis = next(iter(variances), None)
        field = {"timepoint": "timepoint_hours", "dose": "perturbation_dose",
                 "cell_line": "cell_line", "library_prep": "library_selection",
                 "co_treatment": "co_treatment_arm"}.get(str(axis))
        if field and len(members) > 1:
            values = [str(d.get(field)) for _, d in members]
            if len(set(values)) != len(values):
                print(f"      ✗ members do not differ on {field}: {values} — the axis says they "
                      f"should be distinct")
                total_problems += 1
            else:
                print(f"      varies on {field}: {', '.join(values)}")

    print(f"\n{'=' * 60}\n{total_problems} problem(s) found across {len(found)} file(s).")
    return 1 if total_problems else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
