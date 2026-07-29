"""metadata.py rejects values outside the controlled vocabulary, and says what is allowed."""

import gzip
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

import vocab

REPO_ROOT = Path(__file__).resolve().parent.parent
GENOME_BUILD = "hg38plusAkataInverted"


def test_vocab_rejects_unknown_value():
    with pytest.raises(vocab.VocabError) as error:
        vocab.validate("cell_line", "NotARealCellLine")
    message = str(error.value)
    assert "NotARealCellLine" in message
    # The error has to be actionable: it lists the allowed values and where to add one.
    assert "Mutu" in message
    assert "vocab.py" in message


def test_vocab_accepts_known_values():
    assert vocab.validate("cell_line", "SNU719") == "SNU719"
    assert vocab.validate("perturbation_type", "siRNA") == "siRNA"
    assert vocab.organism_for_build(GENOME_BUILD) == "human"


def build_experiment(tmp_path, name="SNU719_Zta-plus-Rta"):
    """A minimal Model_Experiment input directory with two conditions."""
    root = tmp_path / name
    for condition in ("cntl", "test"):
        directory = root / condition
        directory.mkdir(parents=True)
        for read in ("1", "2"):
            with gzip.open(directory / f"{condition}_a_{read}.fq.gz", "wt") as handle:
                handle.write("@read1\nACGT\n+\nIIII\n")
    reference = tmp_path / "referenceFiles" / GENOME_BUILD
    (reference / "annotations").mkdir(parents=True, exist_ok=True)
    return {"fastq_root": root, "reference_dir": tmp_path / "referenceFiles",
            "results_dir": tmp_path / "out"}


@pytest.fixture
def experiment(tmp_path):
    return build_experiment(tmp_path)


def run_metadata(experiment, cell_line):
    (experiment["results_dir"]).mkdir(exist_ok=True)
    environment = dict(os.environ, CONDA_DEFAULT_ENV="rnaseqpipeline")
    return subprocess.run(
        [sys.executable, str(REPO_ROOT / "metadata.py"),
         str(experiment["fastq_root"]), str(experiment["reference_dir"]), GENOME_BUILD,
         "ethan", "PE", str(experiment["results_dir"]),
         "--non-interactive", "--cell-line", cell_line,
         "--perturbation-type", "transfection", "--perturbation-target", "Zta",
         "--perturbation-dose", "NA", "--timepoint-hours", "24",
         "--sequencing-run-date", "2026-01-15"],
        capture_output=True, text=True, env=environment,
    )


def test_metadata_vocab_rejects_unknown(experiment):
    result = run_metadata(experiment, "HeLa")
    assert result.returncode != 0
    assert "invalid cell_line" in result.stderr
    assert "SNU719" in result.stderr
    assert not (experiment["results_dir"] / "metadata.tsv").exists(), \
        "a rejected run must not leave metadata behind"


def test_metadata_accepts_known_vocab_and_derives_sample_ids(experiment):
    result = run_metadata(experiment, "SNU719")
    assert result.returncode == 0, result.stderr

    frame = pd.read_csv(experiment["results_dir"] / "metadata.tsv", sep="\t", dtype=str,
                        keep_default_na=False)
    assert list(frame["sample_id"]) == ["SNU719_Zta-plus-Rta_cntl1", "SNU719_Zta-plus-Rta_test1"]
    assert list(frame["condition"]) == ["cntl", "test"]
    assert list(frame["replicate_index"]) == ["1", "1"]
    assert set(frame["organism"]) == {"human"}
    # Controls carry no perturbation.
    assert frame.loc[frame["condition"] == "cntl", "perturbation_target"].iloc[0] == "NA"
    # Legacy columns survive for rnaseq.py.
    assert {"Sample name", "Control?", "Experiment name"} <= set(frame.columns)
    assert frame["fastq_r1_md5"].str.len().eq(32).all()


def test_experiment_name_may_carry_a_uniquifying_suffix(tmp_path):
    """Only the first underscore is structural, so a date suffix can use either separator.

    Repeating an experiment must be able to produce a distinct experiment_id, because
    experiment_id is what sample_id and comparison_id are built from.
    """
    experiment = build_experiment(tmp_path, "SNU719_Rta-Zta_2025-04")
    result = run_metadata(experiment, "SNU719")
    assert result.returncode == 0, result.stderr

    frame = pd.read_csv(experiment["results_dir"] / "metadata.tsv", sep="\t", dtype=str,
                        keep_default_na=False)
    assert set(frame["experiment_id"]) == {"SNU719_Rta-Zta_2025-04"}
    assert list(frame["sample_id"]) == ["SNU719_Rta-Zta_2025-04_cntl1",
                                        "SNU719_Rta-Zta_2025-04_test1"]


def test_cell_line_still_defaults_from_the_first_segment(tmp_path):
    """The cell line is derived from the model half even when a suffix follows."""
    experiment = build_experiment(tmp_path, "SNU719_Rta-Zta_2025-04")
    environment = dict(os.environ, CONDA_DEFAULT_ENV="rnaseqpipeline")
    experiment["results_dir"].mkdir(exist_ok=True)
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "metadata.py"),
         str(experiment["fastq_root"]), str(experiment["reference_dir"]), GENOME_BUILD,
         "ethan", "PE", str(experiment["results_dir"]),
         "--non-interactive", "--perturbation-type", "none",
         "--timepoint-hours", "NA", "--sequencing-run-date", "NA"],
        capture_output=True, text=True, env=environment,
    )
    assert result.returncode == 0, result.stderr
    frame = pd.read_csv(experiment["results_dir"] / "metadata.tsv", sep="\t", dtype=str)
    assert set(frame["cell_line"]) == {"SNU719"}


def test_experiment_name_needs_both_halves(tmp_path):
    experiment = build_experiment(tmp_path, "SNU719")
    result = run_metadata(experiment, "SNU719")
    assert result.returncode != 0
    assert "Model_Experiment" in result.stderr


def test_pointing_one_level_too_high_says_so(tmp_path):
    """The old rule rejected '1_fastq' by accident; this is the replacement guardrail."""
    build_experiment(tmp_path / "1_fastq", "SNU719_Rta-Zta-2025-04")
    experiment = {"fastq_root": tmp_path / "1_fastq",
                  "reference_dir": tmp_path / "1_fastq" / "referenceFiles",
                  "results_dir": tmp_path / "out"}
    result = run_metadata(experiment, "SNU719")
    assert result.returncode != 0
    assert "one level too high" in result.stderr


def test_metadata_requires_a_terminal_or_a_flag(experiment):
    """Nothing may depend on an interactive prompt: --non-interactive must fail loudly.

    cell_line is not named here because it is derived from the input directory name; the
    first field with no default is perturbation_type.
    """
    environment = dict(os.environ, CONDA_DEFAULT_ENV="rnaseqpipeline")
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "metadata.py"),
         str(experiment["fastq_root"]), str(experiment["reference_dir"]), GENOME_BUILD,
         "ethan", "PE", str(experiment["results_dir"]), "--non-interactive"],
        capture_output=True, text=True, env=environment,
    )
    assert result.returncode != 0
    assert "no terminal to prompt on" in result.stderr
    assert "--perturbation-type" in result.stderr
