#!/usr/bin/env python3
# metadata.py
#
# Usage: metadata.py --root-fastq-dir DIR --output-dir DIR [--quick-input | run/experiment flags]
#
# Every argument is a named flag. --quick-input reads the remaining ones from the input_args.json
# that metadata_form.html downloads, so the two directories are all that has to be typed.
#
# Emits <results_dir>/<investigator>_metadata_<timestamp>.tsv (the file to hand to rnaseq.py)
# and <results_dir>/metadata.tsv (the warehouse copy, identical content).

import argparse
import hashlib
import json
import logging
import os
import re
import sys
from datetime import date, datetime
from itertools import chain
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / "rnaseq_helper_scripts"))

from outputs import create_layout  # noqa: E402
from schemas import METADATA, METADATA_LEGACY_COLUMNS, NA  # noqa: E402
import vocab  # noqa: E402

REQUIRED_ENV = "rnaseqpipeline"

conda_env = os.environ.get("CONDA_DEFAULT_ENV")
if conda_env != REQUIRED_ENV:
    sys.stderr.write(
        f"\nERROR: This script must be run inside the '{REQUIRED_ENV}' conda environment.\n"
        f"Currently active environment: {conda_env or 'None'}\n"
        f"Please run:\n\n    conda activate {REQUIRED_ENV}\n\n"
    )
    sys.exit(1)

READ1_PATTERNS = [
    re.compile(r'(.*)_1\.fastq\.gz'),
    re.compile(r'(.*)_R1_001\.fastq\.gz'),
    re.compile(r'(.*)_1\.fq\.gz'),
    re.compile(r'(.*)_R1_001\.fq\.gz'),
    re.compile(r'(.*)_R1\.fastq\.gz'),
    re.compile(r'(.*)_R1\.fq\.gz'),
]
READ2_PATTERNS = [
    re.compile(r'(.*)_2\.fastq\.gz'),
    re.compile(r'(.*)_R2_001\.fastq\.gz'),
    re.compile(r'(.*)_2\.fq\.gz'),
    re.compile(r'(.*)_R2_001\.fq\.gz'),
    re.compile(r'(.*)_R2\.fastq\.gz'),
    re.compile(r'(.*)_R2\.fq\.gz'),
]
SINGLE_PATTERNS = [
    re.compile(r'(.*)\.fastq\.gz'),
    re.compile(r'(.*)\.fq\.gz'),
]

# Names that mean somebody pasted an example instead of substituting their own path. Cheap to
# check, and the failure it prevents is silent and permanent.
PLACEHOLDER_NAMES = frozenset({"Model_Experiment", "model_experiment", "CellLine_Experiment"})

# Written by metadata_form.html and dropped at the top level of the FASTQ directory, next to the
# test/ and cntl/ subdirectories. The name is fixed because that is the whole interface: the form
# names the download, --quick-input looks for it, and nobody has to type or agree on a path.
QUICK_INPUT_FILENAME = "input_args.json"
QUICK_INPUT_FIELDS = frozenset({
    "reference_dir", "species_name", "investigator_name", "experiment_type",
    "cell_line", "perturbation_type", "induced_program", "library_selection",
    "perturbation_target", "perturbation_dose", "timepoint_hours",
    "sequencing_run_date", "notes",
})


def fail(message: str) -> None:
    """Exit non-zero with a readable message. No placeholders, no partial metadata."""
    sys.stderr.write(f"\nERROR: {message}\n\n")
    logging.error(message)
    sys.exit(1)


def md5_file(path: Path, chunk_size: int = 1 << 22) -> str:
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def condition_from_subdir(subdir_name: str) -> str:
    """`test` or `cntl` from the condition subdirectory name. A hard requirement."""
    lowered = subdir_name.lower()
    if lowered.startswith("test"):
        return "test"
    if lowered.startswith("cntl"):
        return "cntl"
    fail(
        f"condition subdirectory {subdir_name!r} is neither 'test' nor 'cntl'. The FASTQ parent "
        f"directory must contain exactly two subdirectories named 'test' and 'cntl'.\n"
        f"If {subdir_name!r} looks like an experiment name, you pointed one level too high: pass "
        f"the Model_Experiment directory itself, not its parent."
    )


def available_references(reference_dir: Path) -> list[str]:
    if not reference_dir.is_dir():
        fail(f"reference directory not found: {reference_dir}")
    return sorted(d.name for d in reference_dir.iterdir() if d.is_dir())


def parse_fastq_files(fastq_root_dir: Path, layout: str) -> list[dict]:
    """One record per sample: condition, replicate ordering and read paths."""
    records: list[dict] = []
    for subdir in sorted(fastq_root_dir.iterdir()):
        if not subdir.is_dir():
            continue
        condition = condition_from_subdir(subdir.name)
        fastqs = sorted(chain(subdir.rglob('*.fastq.gz'), subdir.rglob('*.fq.gz')))

        if layout == "PE":
            paired: dict[str, dict] = {}
            for fastq in fastqs:
                for pattern in READ1_PATTERNS:
                    match = pattern.match(fastq.name)
                    if match:
                        paired.setdefault(match.group(1), {})["r1"] = fastq
                        break
                for pattern in READ2_PATTERNS:
                    match = pattern.match(fastq.name)
                    if match:
                        paired.setdefault(match.group(1), {})["r2"] = fastq
                        break
            for sample_name in sorted(paired):
                reads = paired[sample_name]
                if "r1" not in reads or "r2" not in reads:
                    fail(
                        f"sample {sample_name!r} in {subdir} is missing its "
                        f"{'R2' if 'r1' in reads else 'R1'} file. Paired-end runs need both."
                    )
                records.append({
                    "sample_name": sample_name,
                    "condition": condition,
                    "fastq_r1": reads["r1"],
                    "fastq_r2": reads["r2"],
                })
        else:
            for fastq in fastqs:
                for pattern in SINGLE_PATTERNS:
                    match = pattern.match(fastq.name)
                    if match:
                        records.append({
                            "sample_name": match.group(1),
                            "condition": condition,
                            "fastq_r1": fastq,
                            "fastq_r2": None,
                        })
                        break
                else:
                    fail(f"FASTQ {fastq.name} does not match any expected single-end pattern.")
    return records


def load_quick_input(fastq_root_dir: Path) -> dict[str, str]:
    """Read the form's ``input_args.json`` from the top level of the FASTQ directory."""
    path = fastq_root_dir / QUICK_INPUT_FILENAME
    if not path.is_file():
        fail(
            f"--quick-input was passed but no {QUICK_INPUT_FILENAME} was found at the top level "
            f"of {fastq_root_dir}.\n"
            f"Fill out metadata_form.html, then drag the file it downloads into that directory, "
            f"alongside the test/ and cntl/ subdirectories.\n"
            f"To supply the values on the command line instead, drop --quick-input and pass the "
            f"run and experiment flags (see --help)."
        )
    try:
        data = json.loads(path.read_text())
    except json.JSONDecodeError as error:
        fail(f"{path} is not valid JSON ({error}). Download a fresh copy from "
             f"metadata_form.html rather than editing it by hand.")
    if not isinstance(data, dict):
        fail(f"{path} must hold a JSON object, got {type(data).__name__}.")
    # A misspelled key would otherwise read as "not supplied" and send the run to a prompt that
    # nothing is watching, or to a missing-value error naming a field the file appears to set.
    unknown = sorted(set(data) - QUICK_INPUT_FIELDS)
    if unknown:
        fail(f"{path} has unrecognized field(s): {', '.join(unknown)}.\n"
             f"recognized fields: {', '.join(sorted(QUICK_INPUT_FIELDS))}")
    return {key: str(value).strip() for key, value in data.items() if value is not None}


def resolve(field: str, value: str | None, *, prompt: str, allowed=None,
            interactive: bool, flag: str, default: str | None = None,
            from_form: str | None = None) -> str:
    """Return a validated value for ``field``, prompting only when a terminal is attached.

    Nothing here depends on being interactive: every field can be supplied by flag, which is
    how run_pipeline_one_shot.sh drives it.

    ``from_form`` is the value read from ``input_args.json``. It is used only when no flag was
    passed, and it is deliberately not re-checked against the vocabulary: the form's <select>
    is what constrains it, and a second copy of the allowed values here could only drift from
    the first. An explicit flag still wins, so a one-off run can override the file without
    editing it.
    """
    if value is None and from_form is not None:
        return from_form
    if value is None and default is not None:
        value = default
    while value is None:
        if not interactive:
            fail(
                f"{field} was not supplied and there is no terminal to prompt on. "
                f"Pass {flag}, or --quick-input to read it from {QUICK_INPUT_FILENAME}."
                + (f"\nallowed values: {', '.join(sorted(allowed))}" if allowed else "")
            )
        entered = input(f"{prompt}: ").strip()
        value = entered or None
    if allowed is not None:
        try:
            vocab.validate(field, value)
        except vocab.VocabError as error:
            fail(str(error))
    return value


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="metadata.py",
        description="Generate the metadata TSV that drives rnaseq.py.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    # The two directories nobody else can know: they name the local copy of the data and where
    # its output goes, so they stay on the command line and out of the form.
    parser.add_argument("--root-fastq-dir", dest="fastq_root_dir", required=True,
                        help="directory containing the test/ and cntl/ subdirectories")
    parser.add_argument("--output-dir", dest="results_dir", required=True,
                        help="output directory")
    parser.add_argument("--quick-input", action="store_true",
                        help=f"read the remaining values from {QUICK_INPUT_FILENAME} at the top "
                             f"level of --root-fastq-dir, as written by metadata_form.html. "
                             f"Without this flag the file is never consulted; with it, an "
                             f"explicit flag still overrides what the file says.")

    run = parser.add_argument_group("run (from --quick-input, a flag, or a prompt)")
    run.add_argument("--reference-dir", dest="reference_dir", help="path to referenceFiles")
    run.add_argument("--species-name", dest="species_name",
                     help="genome build; must match a directory under --reference-dir")
    run.add_argument("--investigator-name", dest="investigator_name",
                     help="investigator name; whitespace is stripped, since it becomes part of "
                          "the output filename")
    run.add_argument("--experiment-type", dest="experiment_type",
                     choices=["PE", "SE", "pe", "se"], help="PE or SE")

    experiment = parser.add_argument_group("experiment structure (prompted if omitted)")
    experiment.add_argument("--cell-line", help=f"one of: {', '.join(sorted(vocab.CELL_LINES))}")
    experiment.add_argument("--organism", help=f"one of: {', '.join(sorted(vocab.ORGANISMS))}")
    experiment.add_argument("--perturbation-type",
                            help=f"one of: {', '.join(sorted(vocab.PERTURBATION_TYPES))}")
    experiment.add_argument("--induced-program",
                            help=f"biological program induced; one of: "
                                 f"{', '.join(sorted(vocab.INDUCED_PROGRAMS))}")
    experiment.add_argument("--library-selection",
                            help=f"RNA selection before library prep; one of: "
                                 f"{', '.join(sorted(vocab.LIBRARY_SELECTIONS))}")
    experiment.add_argument("--perturbation-target", help="e.g. BMRF1, SRSF1, CC115")
    experiment.add_argument("--perturbation-dose", help="e.g. 100nM")
    experiment.add_argument("--timepoint-hours", help="hours post perturbation, or NA")
    experiment.add_argument("--sequencing-run-date", help="ISO 8601 date, or NA")
    experiment.add_argument("--notes", help="free text")
    experiment.add_argument("--non-interactive", action="store_true",
                            help="never prompt; missing values are an error")
    return parser.parse_args(argv)


def main(argv: list[str]) -> None:
    args = parse_args(argv)

    fastq_root_dir = Path(args.fastq_root_dir).resolve()
    results_dir = Path(args.results_dir).resolve()

    if not fastq_root_dir.is_dir():
        fail(f"FASTQ directory not found: {fastq_root_dir}")

    form = load_quick_input(fastq_root_dir) if args.quick_input else {}

    def run_arg(field: str, flag: str) -> str:
        value = getattr(args, field) or form.get(field)
        if not value:
            fail(f"{field} was not supplied. Pass {flag}, or --quick-input to read it from "
                 f"{QUICK_INPUT_FILENAME} at the top level of {fastq_root_dir}.")
        return value

    reference_dir = Path(run_arg("reference_dir", "--reference-dir")).resolve()
    genome_build = run_arg("species_name", "--species-name")
    # Whitespace would produce an output filename that the pipeline's own
    # <investigator>_metadata_*.tsv glob no longer matches.
    investigator = re.sub(r"\s+", "", run_arg("investigator_name", "--investigator-name"))
    if not investigator:
        fail("investigator_name is blank once whitespace is stripped.")
    library_layout = run_arg("experiment_type", "--experiment-type").upper()
    # Not vocabulary policing: this picks which FASTQ pairing branch runs, so an unrecognized
    # value would quietly parse a paired-end run as single-end.
    if library_layout not in vocab.LIBRARY_LAYOUTS:
        fail(f"experiment_type must be PE or SE, got {library_layout!r}.")
    # Only the first underscore is structural: it separates the cell model from the experiment
    # descriptor. Everything after it is free, so a uniquifying suffix can use either separator.
    model, _, descriptor = fastq_root_dir.name.partition('_')
    if not model or not descriptor:
        fail(
            f"FASTQ directory name {fastq_root_dir.name!r} must be in 'Model_Experiment' format, "
            f"e.g. 'SNU719_Rta-Zta-2025-04'. The first underscore separates the cell model from "
            f"the experiment descriptor; both parts must be non-empty."
        )
    experiment_id = fastq_root_dir.name
    if experiment_id in PLACEHOLDER_NAMES:
        fail(
            f"the FASTQ directory is named {experiment_id!r}, which is the documentation "
            f"placeholder rather than a real experiment. This name becomes experiment_id and is "
            f"baked into every sample_id, so it must be the actual experiment directory name.\n"
            f"If you are running under Docker, mount to the real name:\n"
            f"    -v \"/host/path/SNU719_Rta-Zta-2025-04-10:/data/SNU719_Rta-Zta-2025-04-10:ro\""
        )

    references = available_references(reference_dir)
    if genome_build not in references:
        fail(f"reference {genome_build!r} not available. Available: {', '.join(references)}")

    create_layout(results_dir)
    timestamp = datetime.now().strftime("%m%d%Y_%H%M%S")
    logging.basicConfig(
        filename=str(results_dir / "logs" / f"metadata_{timestamp}.log"),
        level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s',
    )
    logging.info("Starting metadata generation for %s", experiment_id)

    interactive = sys.stdin.isatty() and not args.non_interactive

    cell_line = resolve(
        "cell_line", args.cell_line, prompt=f"Cell line for {experiment_id}",
        allowed=vocab.CELL_LINES, interactive=interactive, flag="--cell-line",
        # The Model half of Model_Experiment is the cell line by convention.
        default=model if model in vocab.CELL_LINES else None,
        from_form=form.get("cell_line"),
    )
    organism = resolve(
        "organism", args.organism, prompt=f"Organism for {genome_build}",
        allowed=vocab.ORGANISMS, interactive=interactive, flag="--organism",
        default=vocab.organism_for_build(genome_build),
    )
    perturbation_type = resolve(
        "perturbation_type", args.perturbation_type,
        prompt=f"Perturbation type ({'/'.join(sorted(vocab.PERTURBATION_TYPES))})",
        allowed=vocab.PERTURBATION_TYPES, interactive=interactive, flag="--perturbation-type",
        from_form=form.get("perturbation_type"),
    )
    induced_program = resolve(
        "induced_program", args.induced_program,
        prompt=f"Induced program ({'/'.join(sorted(vocab.INDUCED_PROGRAMS))})",
        allowed=vocab.INDUCED_PROGRAMS, interactive=interactive, flag="--induced-program",
        from_form=form.get("induced_program"),
    )
    if perturbation_type == "none":
        perturbation_target = NA
        perturbation_dose = NA
    else:
        perturbation_target = resolve(
            "perturbation_target", args.perturbation_target,
            prompt="Perturbation target (e.g. BMRF1)", interactive=interactive,
            flag="--perturbation-target", from_form=form.get("perturbation_target"),
        )
        perturbation_dose = resolve(
            "perturbation_dose", args.perturbation_dose,
            prompt="Perturbation dose (e.g. 100nM, or NA)", interactive=interactive,
            flag="--perturbation-dose", from_form=form.get("perturbation_dose"),
        )
    timepoint_hours = resolve(
        "timepoint_hours", args.timepoint_hours, prompt="Timepoint in hours (or NA)",
        interactive=interactive, flag="--timepoint-hours",
        from_form=form.get("timepoint_hours"),
    )
    # Canonicalisation rather than a second verifier: the column has to read the same whether a
    # value arrived as "24" from the form or "24.0" from a flag, or a query grouping on it splits
    # one timepoint in two.
    if timepoint_hours != NA:
        try:
            timepoint_hours = str(float(timepoint_hours))
        except ValueError:
            fail(f"timepoint_hours must be a number or NA, got {timepoint_hours!r}")
    # No default, deliberately. Most of the archive is polyA, so defaulting would be right most
    # of the time and silently wrong for the ribodepleted runs -- and wrong in the direction that
    # makes a query look answerable when it is not.
    library_selection = resolve(
        "library_selection", args.library_selection,
        prompt=f"Library selection ({'/'.join(sorted(vocab.LIBRARY_SELECTIONS))})",
        allowed=vocab.LIBRARY_SELECTIONS, interactive=interactive, flag="--library-selection",
        from_form=form.get("library_selection"),
    )
    sequencing_run_date = resolve(
        "sequencing_run_date", args.sequencing_run_date,
        prompt="Sequencing run date (YYYY-MM-DD, or NA)", interactive=interactive,
        flag="--sequencing-run-date", from_form=form.get("sequencing_run_date"),
    )
    if sequencing_run_date != NA:
        try:
            sequencing_run_date = date.fromisoformat(sequencing_run_date).isoformat()
        except ValueError:
            fail(f"sequencing_run_date must be ISO 8601 (YYYY-MM-DD) or NA, "
                 f"got {sequencing_run_date!r}")

    # The one field with nothing to say when it is empty, so it never prompts.
    notes = args.notes if args.notes is not None else form.get("notes", "")

    records = parse_fastq_files(fastq_root_dir, library_layout)
    if not records:
        fail(f"no {library_layout} FASTQ files found under {fastq_root_dir}")

    conditions_found = {record["condition"] for record in records}
    missing_conditions = sorted(vocab.CONDITIONS - conditions_found)
    if missing_conditions:
        fail(f"no samples found for condition(s): {', '.join(missing_conditions)}")

    # Controls first, then tests, each ordered by filename. Replicate indices follow that order
    # so they are stable across reruns of the same input directory.
    records.sort(key=lambda r: (r["condition"] != "cntl", r["fastq_r1"].name))
    replicate_counter: dict[str, int] = {}

    rows = []
    for record in records:
        condition = record["condition"]
        replicate_counter[condition] = replicate_counter.get(condition, 0) + 1
        replicate_index = replicate_counter[condition]
        r1, r2 = record["fastq_r1"], record["fastq_r2"]

        print(f"Checksumming {r1.name}" + (f" and {r2.name}" if r2 else ""), flush=True)
        rows.append({
            # Spec columns.
            "sample_id": f"{experiment_id}_{condition}{replicate_index}",
            "experiment_id": experiment_id,
            "condition": condition,
            "replicate_index": replicate_index,
            "cell_line": cell_line,
            "organism": organism,
            "genome_build": genome_build,
            "perturbation_type": perturbation_type if condition == "test" else "none",
            # Experiment-level, unlike the perturbation_* fields: it describes the design, so it
            # is on control rows too. "Every sample from a reactivation experiment" is then one
            # filter rather than a subquery back through experiment_id.
            "induced_program": induced_program,
            "perturbation_target": perturbation_target if condition == "test" else NA,
            "perturbation_dose": perturbation_dose if condition == "test" else NA,
            "timepoint_hours": timepoint_hours,
            "library_layout": library_layout,
            "library_selection": library_selection,
            # Filled in by rnaseq.py once RSeQC has run (Task 14).
            "strandedness": NA,
            "fastq_r1": str(r1),
            "fastq_r2": str(r2) if r2 else NA,
            "fastq_r1_md5": md5_file(r1),
            "fastq_r2_md5": md5_file(r2) if r2 else NA,
            "investigator": investigator,
            "sequencing_run_date": sequencing_run_date,
            "notes": notes or NA,
            # Legacy columns rnaseq.py still reads by name.
            "Path Read 1": str(r1),
            "Path Read 2": str(r2) if r2 else NA,
            "Species": genome_build,
            "Sample name": record["sample_name"],
            "Condition": "test" if condition == "test" else "control",
            "Control?": "no" if condition == "test" else "yes",
            "ConditionReplicate": f"{condition}{replicate_index}",
            "Experiment name": experiment_id,
        })

    duplicates = pd.Series([r["Sample name"] for r in rows])
    if duplicates.duplicated().any():
        fail(
            "duplicate sample names across conditions: "
            f"{sorted(duplicates[duplicates.duplicated()].unique())}. Sample names must be "
            "unique because they become output directory names."
        )

    column_order = list(METADATA.column_names) + list(METADATA_LEGACY_COLUMNS)
    metadata_df = pd.DataFrame(rows)[column_order]

    metadata_filename = results_dir / f"{investigator}_metadata_{timestamp}.tsv"
    metadata_df.to_csv(metadata_filename, sep='\t', index=False)
    metadata_df.to_csv(results_dir / "metadata.tsv", sep='\t', index=False)

    logging.info("Metadata file generated: %s", metadata_filename)
    print(f"Metadata file successfully generated: {metadata_filename}")
    print(f"Warehouse copy: {results_dir / 'metadata.tsv'}")
    print(f"{len(rows)} samples: " + ", ".join(r["sample_id"] for r in rows))


if __name__ == "__main__":
    main(sys.argv[1:])
