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
    assert vocab.validate("perturbation_type", "chemical") == "chemical"
    assert vocab.validate("library_selection", "ribodepleted") == "ribodepleted"
    assert vocab.validate("library_strandedness", "stranded") == "stranded"
    assert vocab.validate("facs_gfp_promoter", "BMRF1p") == "BMRF1p"
    assert vocab.organism_for_build(GENOME_BUILD) == "human"


def test_sirna_is_a_co_treatment_not_a_perturbation_type():
    """It describes a second treatment alongside the first, not how the first was delivered."""
    assert vocab.validate("co_treatment", "siRNA") == "siRNA"
    with pytest.raises(vocab.VocabError):
        vocab.validate("perturbation_type", "siRNA")


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


def run_args(experiment):
    """The four run flags that used to be positional arguments."""
    return ["--root-fastq-dir", str(experiment["fastq_root"]),
            "--output-dir", str(experiment["results_dir"]),
            "--reference-dir", str(experiment["reference_dir"]),
            "--species-name", GENOME_BUILD,
            "--investigator-name", "ethan",
            "--experiment-type", "PE"]


def run_metadata(experiment, cell_line, library_selection="polyA"):
    (experiment["results_dir"]).mkdir(exist_ok=True)
    environment = dict(os.environ, CONDA_DEFAULT_ENV="rnaseqpipeline")
    return subprocess.run(
        [sys.executable, str(REPO_ROOT / "metadata.py"), *run_args(experiment),
         "--non-interactive", "--cell-line", cell_line,
         "--perturbation-type", "transfection", "--perturbation-agent", "Zta",
         "--perturbation-target", "NA",
         "--perturbation-dose", "NA", "--timepoint-hours", "24",
         "--co-treatment", "none", "--facs-purified", "yes",
         "--facs-gfp-promoter", "pCMV",
         "--library-selection", library_selection,
         "--library-strandedness", "stranded",
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
    # "Sample name" survives for rnaseq.py; the other historical mirrors were retired as
    # pure duplicates of condition/experiment_id/fastq_r1/fastq_r2/replicate_index.
    assert {"Sample name"} <= set(frame.columns)
    assert not {"Path Read 1", "Path Read 2", "Species", "Condition", "Control?",
                "ConditionReplicate", "Experiment name"} & set(frame.columns)
    assert frame["fastq_r1_md5"].str.len().eq(32).all()


def test_co_treatment_is_recorded_apart_from_the_perturbation(experiment):
    """Zta+PAA has to stay two facts, or "every PAA experiment" stops being a query."""
    environment = dict(os.environ, CONDA_DEFAULT_ENV="rnaseqpipeline")
    experiment["results_dir"].mkdir(exist_ok=True)
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "metadata.py"), *run_args(experiment),
         "--non-interactive", "--cell-line", "SNU719",
         "--perturbation-type", "transfection", "--perturbation-agent", "Zta",
         "--perturbation-target", "NA", "--perturbation-dose", "NA",
         "--timepoint-hours", "24", "--co-treatment", "PAA",
         "--co-treatment-target", "PAA_replication",
         "--co-treatment-arm", "test",
         "--series-label", "SNU719_Zta-plus-Rta_PAA_2025-04-10",
         "--facs-purified", "yes", "--facs-gfp-promoter", "pCMV",
         "--library-selection", "polyA", "--library-strandedness", "stranded",
         "--sequencing-run-date", "2025-04-10"],
        capture_output=True, text=True, env=environment,
    )
    assert result.returncode == 0, result.stderr

    frame = pd.read_csv(experiment["results_dir"] / "metadata.tsv", sep="\t", dtype=str,
                        keep_default_na=False)
    # A co-treatment (real or control/mock) is a fact about the whole experiment, applied
    # uniformly -- unlike the primary perturbation, it is not test-condition-only.
    assert set(frame["co_treatment"]) == {"PAA"}
    assert set(frame["co_treatment_target"]) == {"PAA_replication"}
    assert set(frame["co_treatment_arm"]) == {"test"}
    # A co-treatment pair is a series, and its variance is forced rather than asked for.
    assert set(frame["series_label"]) == {"SNU719_Zta-plus-Rta_PAA_2025-04-10"}
    assert set(frame["series_variance"]) == {"co_treatment"}
    test_row = frame.loc[frame["condition"] == "test"].iloc[0]
    assert test_row["perturbation_agent"] == "Zta"
    cntl_row = frame.loc[frame["condition"] == "cntl"].iloc[0]
    assert cntl_row["perturbation_agent"] == "NA"
    # FACS describes the material, so it is on both arms.
    assert set(frame["facs_purified"]) == {"yes"}


def _co_treatment_args(experiment):
    return [
        sys.executable, str(REPO_ROOT / "metadata.py"), *run_args(experiment),
        "--non-interactive", "--cell-line", "SNU719",
        "--perturbation-type", "transfection", "--perturbation-agent", "Zta",
        "--perturbation-target", "NA", "--perturbation-dose", "NA",
        "--timepoint-hours", "24", "--co-treatment", "CRISPR",
        "--co-treatment-target", "EXOSC3",
        "--facs-purified", "yes", "--facs-gfp-promoter", "pCMV",
        "--library-selection", "polyA", "--library-strandedness", "stranded",
        "--sequencing-run-date", "2025-08-31",
    ]


def _plain_args(experiment):
    return [
        sys.executable, str(REPO_ROOT / "metadata.py"), *run_args(experiment),
        "--non-interactive", "--cell-line", "SNU719",
        "--perturbation-type", "transfection", "--perturbation-agent", "Zta",
        "--perturbation-target", "NA", "--perturbation-dose", "NA",
        "--timepoint-hours", "24", "--co-treatment", "none",
        "--facs-purified", "yes", "--facs-gfp-promoter", "pCMV",
        "--library-selection", "polyA", "--library-strandedness", "stranded",
        "--sequencing-run-date", "2025-08-31",
    ]


def test_co_treatment_requires_arm_and_series_label(experiment):
    """A co-treatment always comes as a pair, so both halves of that fact are required."""
    environment = dict(os.environ, CONDA_DEFAULT_ENV="rnaseqpipeline")
    base = _co_treatment_args(experiment)

    # --co-treatment-arm omitted.
    result = subprocess.run([*base, "--series-label", "Mutu_Zta_EXOSC3_2025-08-31"],
                            capture_output=True, text=True, env=environment)
    assert result.returncode != 0

    # --series-label omitted.
    result = subprocess.run([*base, "--co-treatment-arm", "test"],
                            capture_output=True, text=True, env=environment)
    assert result.returncode != 0


def test_co_treatment_forces_series_variance(experiment):
    """The one right answer is filled in, and an attempt to override it is refused."""
    environment = dict(os.environ, CONDA_DEFAULT_ENV="rnaseqpipeline")
    base = [*_co_treatment_args(experiment), "--co-treatment-arm", "test",
            "--series-label", "Mutu_Zta_EXOSC3_2025-08-31"]

    # Not supplied at all: it is still set, without prompting.
    result = subprocess.run(base, capture_output=True, text=True, env=environment)
    assert result.returncode == 0, result.stderr
    frame = pd.read_csv(experiment["results_dir"] / "metadata.tsv", sep="\t", dtype=str,
                        keep_default_na=False)
    assert set(frame["series_variance"]) == {"co_treatment"}

    # Contradicting it is an error rather than a silent overwrite.
    result = subprocess.run([*base, "--series-variance", "timepoint"],
                            capture_output=True, text=True, env=environment)
    assert result.returncode != 0
    assert "forces it to 'co_treatment'" in result.stderr


def test_series_label_and_variance_are_mutually_inclusive(experiment):
    """Either both carry a real value or both are NA -- never one without the other."""
    environment = dict(os.environ, CONDA_DEFAULT_ENV="rnaseqpipeline")
    base = _plain_args(experiment)

    # Label without variance.
    result = subprocess.run([*base, "--series-label", "Akata_anti-IgG_2022-12-08"],
                            capture_output=True, text=True, env=environment)
    assert result.returncode != 0
    assert "mutually" in result.stderr

    # Variance without label.
    result = subprocess.run([*base, "--series-variance", "timepoint"],
                            capture_output=True, text=True, env=environment)
    assert result.returncode != 0
    assert "mutually" in result.stderr

    # Neither: both land as NA, and no series is invented.
    result = subprocess.run(base, capture_output=True, text=True, env=environment)
    assert result.returncode == 0, result.stderr
    frame = pd.read_csv(experiment["results_dir"] / "metadata.tsv", sep="\t", dtype=str,
                        keep_default_na=False)
    assert set(frame["series_label"]) == {"NA"}
    assert set(frame["series_variance"]) == {"NA"}

    # Both: recorded on every row, since series membership is a whole-experiment fact.
    result = subprocess.run(
        [*base, "--series-label", "Akata_anti-IgG_2022-12-08", "--series-variance", "timepoint"],
        capture_output=True, text=True, env=environment)
    assert result.returncode == 0, result.stderr
    frame = pd.read_csv(experiment["results_dir"] / "metadata.tsv", sep="\t", dtype=str,
                        keep_default_na=False)
    assert set(frame["series_label"]) == {"Akata_anti-IgG_2022-12-08"}
    assert set(frame["series_variance"]) == {"timepoint"}


def test_series_variance_co_treatment_rejected_without_a_co_treatment(experiment):
    """That value is set automatically; claiming it by hand would be a lie about the design."""
    result = subprocess.run(
        [*_plain_args(experiment), "--series-label", "X_2020-01-01",
         "--series-variance", "co_treatment"],
        capture_output=True, text=True, env=dict(os.environ, CONDA_DEFAULT_ENV="rnaseqpipeline"),
    )
    assert result.returncode != 0
    assert "co_treatment is 'none'" in result.stderr


def test_co_treatment_target_rejects_free_text(experiment):
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "metadata.py"), *run_args(experiment),
         "--non-interactive", "--cell-line", "SNU719",
         "--perturbation-type", "transfection", "--perturbation-agent", "Zta",
         "--perturbation-target", "NA", "--perturbation-dose", "NA",
         "--timepoint-hours", "24", "--co-treatment", "siRNA",
         "--co-treatment-target", "CNOT-1",
         "--facs-purified", "yes", "--facs-gfp-promoter", "pCMV",
         "--library-selection", "polyA", "--library-strandedness", "stranded",
         "--sequencing-run-date", "2025-04-10"],
        capture_output=True, text=True,
        env=dict(os.environ, CONDA_DEFAULT_ENV="rnaseqpipeline"),
    )
    assert result.returncode != 0
    assert "invalid co_treatment_target" in result.stderr
    assert "CNOT1" in result.stderr


def test_library_selection_is_recorded_on_every_sample(experiment, tmp_path):
    """polyA and ribodepleted libraries see different transcriptomes, so it is a stored field.

    It is per sample rather than encoded in experiment_id: the archive is mostly polyA with a
    few ribodepleted runs, and a query that pools them silently compares a transcript that was
    selected away against one that was retained.
    """
    result = run_metadata(experiment, "SNU719")
    assert result.returncode == 0, result.stderr
    frame = pd.read_csv(experiment["results_dir"] / "metadata.tsv", sep="\t", dtype=str,
                        keep_default_na=False)
    assert set(frame["library_selection"]) == {"polyA"}

    ribo = build_experiment(tmp_path / "ribo", "Mutu_total-RNA_2025-06")
    result = run_metadata(ribo, "Mutu", library_selection="ribodepleted")
    assert result.returncode == 0, result.stderr
    frame = pd.read_csv(ribo["results_dir"] / "metadata.tsv", sep="\t", dtype=str,
                        keep_default_na=False)
    assert set(frame["library_selection"]) == {"ribodepleted"}


def test_library_selection_rejects_free_text(experiment):
    result = run_metadata(experiment, "SNU719", library_selection="rRNA-depleted")
    assert result.returncode != 0
    assert "invalid library_selection" in result.stderr
    assert "ribodepleted" in result.stderr


def test_library_selection_has_no_default(experiment):
    """Most of the archive is polyA, which is exactly why guessing is not allowed.

    A default would be right most of the time and wrong for the ribodepleted runs, and wrong
    in the direction where nothing downstream can detect it.
    """
    environment = dict(os.environ, CONDA_DEFAULT_ENV="rnaseqpipeline")
    experiment["results_dir"].mkdir(exist_ok=True)
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "metadata.py"), *run_args(experiment),
         "--non-interactive", "--cell-line", "SNU719",
         "--perturbation-type", "transfection", "--perturbation-agent", "Zta",
         "--perturbation-target", "NA", "--perturbation-dose", "NA",
         "--timepoint-hours", "24", "--co-treatment", "none", "--facs-purified", "yes",
         "--facs-gfp-promoter", "pCMV", "--library-strandedness", "stranded",
         "--sequencing-run-date", "2025-04-10"],
        capture_output=True, text=True, env=environment,
    )
    assert result.returncode != 0
    assert "--library-selection" in result.stderr
    assert not (experiment["results_dir"] / "metadata.tsv").exists()


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
        [sys.executable, str(REPO_ROOT / "metadata.py"), *run_args(experiment),
         "--non-interactive", "--perturbation-type", "none",
         "--library-selection", "polyA", "--library-strandedness", "unknown",
         "--co-treatment", "none", "--facs-purified", "no",
         "--facs-gfp-promoter", "none",
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


def test_placeholder_experiment_name_is_rejected(tmp_path):
    """A container mounted at a fixed /data/Model_Experiment would stamp every run alike."""
    experiment = build_experiment(tmp_path, "Model_Experiment")
    result = run_metadata(experiment, "SNU719")
    assert result.returncode != 0
    assert "documentation placeholder" in result.stderr


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
        [sys.executable, str(REPO_ROOT / "metadata.py"), *run_args(experiment),
         "--non-interactive"],
        capture_output=True, text=True, env=environment,
    )
    assert result.returncode != 0
    assert "no terminal to prompt on" in result.stderr
    assert "--perturbation-type" in result.stderr
