"""--quick-input reads the form's input_args.json instead of a dozen flags.

The point of the file is that Erik never types a path or a flag: he fills out metadata_form.html,
drags the download into the experiment folder, and the two directory arguments are all that is
left for whoever runs the pipeline.
"""

import json
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd

from test_metadata_vocab import GENOME_BUILD, build_experiment

REPO_ROOT = Path(__file__).resolve().parent.parent

FORM_OUTPUT = {
    "reference_dir": None,  # filled in per experiment; the form hardcodes the container path
    "species_name": GENOME_BUILD,
    "investigator_name": "ethan",
    "experiment_type": "PE",
    "cell_line": "SNU719",
    "perturbation_type": "transfection",
    "perturbation_agent": "Zta",
    "perturbation_target": "NA",
    "perturbation_dose": "NA",
    "co_treatment": "none",
    "co_treatment_target": "NA",
    "facs_purified": "yes",
    "facs_gfp_promoter": "pCMV",
    "timepoint_hours": "24",
    "library_selection": "polyA",
    "library_strandedness": "stranded",
    "sequencing_run_date": "2026-01-15",
    "notes": "",
}


def write_form_output(experiment, **overrides):
    """Drop an input_args.json at the top level of the FASTQ directory, as the form's download."""
    payload = dict(FORM_OUTPUT, reference_dir=str(experiment["reference_dir"]))
    payload.update(overrides)
    payload = {key: value for key, value in payload.items() if value is not None}
    path = experiment["fastq_root"] / "input_args.json"
    path.write_text(json.dumps(payload, indent=2) + "\n")
    return path


def run(experiment, *extra):
    experiment["results_dir"].mkdir(exist_ok=True)
    environment = dict(os.environ, CONDA_DEFAULT_ENV="rnaseqpipeline")
    return subprocess.run(
        [sys.executable, str(REPO_ROOT / "metadata.py"),
         "--root-fastq-dir", str(experiment["fastq_root"]),
         "--output-dir", str(experiment["results_dir"]),
         "--non-interactive", *extra],
        capture_output=True, text=True, env=environment,
    )


def test_quick_input_supplies_everything_but_the_two_directories(tmp_path):
    experiment = build_experiment(tmp_path)
    write_form_output(experiment)

    result = run(experiment, "--quick-input")
    assert result.returncode == 0, result.stderr

    frame = pd.read_csv(experiment["results_dir"] / "metadata.tsv", sep="\t", dtype=str,
                        keep_default_na=False)
    assert set(frame["cell_line"]) == {"SNU719"}
    assert set(frame["library_selection"]) == {"polyA"}
    assert set(frame["co_treatment"]) == {"none"}
    assert set(frame["facs_gfp_promoter"]) == {"pCMV"}
    assert set(frame["library_strandedness"]) == {"stranded"}
    assert set(frame["investigator"]) == {"ethan"}
    assert set(frame["library_layout"]) == {"PE"}
    # Derived from the genome build rather than asked for, which is why the form has no field.
    assert set(frame["organism"]) == {"human"}
    # Normalised the same way a flag-supplied value is, so grouping on the column cannot split
    # one timepoint into "24" and "24.0".
    assert set(frame["timepoint_hours"]) == {"24.0"}


def test_the_file_is_ignored_without_the_flag(tmp_path):
    """--quick-input is opt-in: a file lying in the directory must not take effect by itself."""
    experiment = build_experiment(tmp_path)
    write_form_output(experiment)

    result = run(experiment)
    assert result.returncode != 0
    assert "--quick-input" in result.stderr


def test_a_missing_file_says_where_it_belongs(tmp_path):
    experiment = build_experiment(tmp_path)

    result = run(experiment, "--quick-input")
    assert result.returncode != 0
    assert "input_args.json" in result.stderr
    assert "metadata_form.html" in result.stderr
    assert not (experiment["results_dir"] / "metadata.tsv").exists()


def test_a_misspelled_field_is_rejected_rather_than_ignored(tmp_path):
    """Silently ignoring it would send the run to a missing-value error naming a field the file
    appears to set, which is the least debuggable outcome available."""
    experiment = build_experiment(tmp_path)
    path = write_form_output(experiment)
    payload = json.loads(path.read_text())
    payload["cell_lines"] = payload.pop("cell_line")
    path.write_text(json.dumps(payload))

    result = run(experiment, "--quick-input")
    assert result.returncode != 0
    assert "cell_lines" in result.stderr
    assert "unrecognized" in result.stderr


def test_malformed_json_points_back_at_the_form(tmp_path):
    experiment = build_experiment(tmp_path)
    (experiment["fastq_root"] / "input_args.json").write_text("{not json")

    result = run(experiment, "--quick-input")
    assert result.returncode != 0
    assert "not valid JSON" in result.stderr


def test_an_explicit_flag_overrides_the_file(tmp_path):
    """A one-off rerun should not need the file edited and re-downloaded."""
    experiment = build_experiment(tmp_path)
    write_form_output(experiment)

    result = run(experiment, "--quick-input", "--library-selection", "ribodepleted")
    assert result.returncode == 0, result.stderr

    frame = pd.read_csv(experiment["results_dir"] / "metadata.tsv", sep="\t", dtype=str,
                        keep_default_na=False)
    assert set(frame["library_selection"]) == {"ribodepleted"}


def test_the_file_does_not_disturb_sample_discovery(tmp_path):
    """It sits alongside test/ and cntl/, and parse_fastq_files walks directories only."""
    experiment = build_experiment(tmp_path)
    write_form_output(experiment)

    result = run(experiment, "--quick-input")
    assert result.returncode == 0, result.stderr

    frame = pd.read_csv(experiment["results_dir"] / "metadata.tsv", sep="\t", dtype=str,
                        keep_default_na=False)
    assert list(frame["sample_id"]) == ["SNU719_Zta-plus-Rta_cntl1",
                                        "SNU719_Zta-plus-Rta_test1"]


def test_investigator_whitespace_is_stripped(tmp_path):
    """Kept whitespace-free for consistency as a queryable metadata.tsv column value."""
    experiment = build_experiment(tmp_path)
    write_form_output(experiment, investigator_name="ethan israel")

    result = run(experiment, "--quick-input")
    assert result.returncode == 0, result.stderr

    frame = pd.read_csv(experiment["results_dir"] / "metadata.tsv", sep="\t", dtype=str,
                        keep_default_na=False)
    assert set(frame["investigator"]) == {"ethanisrael"}


def test_only_metadata_tsv_is_written(tmp_path):
    """There is exactly one metadata file; no separate timestamped copy."""
    experiment = build_experiment(tmp_path)
    write_form_output(experiment)

    result = run(experiment, "--quick-input")
    assert result.returncode == 0, result.stderr

    tsvs = sorted(p.name for p in experiment["results_dir"].glob("*.tsv"))
    assert tsvs == ["metadata.tsv"]
