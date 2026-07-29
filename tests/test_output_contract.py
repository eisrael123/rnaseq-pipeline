"""Contract tests: the shape of a results directory, not the biology in it."""

import json

import pandas as pd
import pytest

import validate_outputs
from schemas import REQUIRED_TABLES, table as get_table

from conftest import RUN_ID


def test_manifest_complete(results_dir):
    """Every provenance field is present and populated for a finished run."""
    manifest = json.loads((results_dir / "run_manifest.json").read_text())
    required = [
        "run_id", "experiment_id", "pipeline_version", "git_commit", "git_dirty",
        "docker_image_digest", "genome_build", "annotation_version", "reference_dir",
        "reference_dir_sha256", "investigator", "library_layout", "run_start_utc",
        "run_end_utc", "exit_status", "tool_versions", "parameters",
    ]
    assert [key for key in required if key not in manifest] == []
    assert validate_outputs._find_nulls(manifest) == []
    assert manifest["exit_status"] == "success"
    assert manifest["tool_versions"], "tool_versions must not be empty"


@pytest.mark.parametrize("name", REQUIRED_TABLES)
def test_schema_columns(results_dir, name):
    """Each table's header is exactly the schema's columns, in order."""
    spec = get_table(name)
    path = results_dir / spec.relpath
    assert path.is_file(), f"{spec.relpath} was not written"
    header = path.read_text().split("\n", 1)[0].rstrip("\r").split("\t")
    assert header == spec.column_names


@pytest.mark.parametrize("name", REQUIRED_TABLES)
def test_run_id_present(results_dir, name):
    """run_id is on every row of every table and matches the manifest."""
    spec = get_table(name)
    frame = pd.read_csv(results_dir / spec.relpath, sep="\t", dtype=str)
    assert "run_id" in frame.columns
    if frame.empty:
        return
    assert set(frame["run_id"]) == {RUN_ID}


def test_validator_passes_on_a_good_run(results_dir):
    report = validate_outputs.validate(results_dir)
    assert report.failures == []


def test_validator_catches_missing_column(results_dir):
    """Dropping a declared column is a failure, not a silently narrower table."""
    path = results_dir / get_table("de_gene").relpath
    frame = pd.read_csv(path, sep="\t", dtype=str).drop(columns=["lfc_se"])
    frame.to_csv(path, sep="\t", index=False)

    report = validate_outputs.validate(results_dir)
    assert any("de_gene.tsv columns do not match" in failure for failure in report.failures)


def test_validator_catches_unknown_sample_id(results_dir):
    path = results_dir / get_table("qc_metric").relpath
    frame = pd.read_csv(path, sep="\t", dtype=str)
    frame.loc[0, "sample_id"] = "not_a_sample"
    frame.to_csv(path, sep="\t", index=False)

    report = validate_outputs.validate(results_dir)
    assert any("not in metadata.tsv" in failure for failure in report.failures)


def test_validator_catches_checksum_mismatch(results_dir):
    target = results_dir / get_table("junction").relpath
    target.write_text(target.read_text() + "\n")

    report = validate_outputs.validate(results_dir)
    assert any("sha256 mismatch" in failure for failure in report.failures)


def test_de_gene_holds_no_per_sample_columns(results_dir):
    """The whole point of Task 5: statistics only, no sample names in the header."""
    header = (results_dir / get_table("de_gene").relpath).read_text().split("\n")[0]
    metadata = pd.read_csv(results_dir / "metadata.tsv", sep="\t", dtype=str)
    for sample_id in metadata["sample_id"]:
        assert sample_id not in header


def test_expression_tpm_sums_to_a_million(results_dir):
    frame = pd.read_csv(results_dir / get_table("expression_transcript").relpath, sep="\t")
    for sample_id, group in frame.groupby("sample_id"):
        assert abs(group["tpm"].sum() - 1e6) < 1.0, sample_id


def test_expression_gene_row_count_is_rectangular(results_dir):
    frame = pd.read_csv(results_dir / get_table("expression_gene").relpath, sep="\t")
    per_sample = frame.groupby("sample_id")["gene_id"].nunique()
    assert per_sample.nunique() == 1
    assert len(frame) == int(per_sample.iloc[0]) * len(per_sample)


def test_missing_values_are_na_not_blank(results_dir):
    """padj is NA for the independently-filtered gene, and it is spelled NA."""
    text = (results_dir / get_table("de_gene").relpath).read_text()
    assert "\tNA" in text
    for line in text.strip().split("\n"):
        assert "\t\t" not in line, "an empty cell escaped the NA convention"
