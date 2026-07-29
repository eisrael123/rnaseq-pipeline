"""rMATS normalization: coordinate slots, event_key stability, replicate explosion."""

import pandas as pd
import pytest

import tables_rmats
from outputs import RunContext
from schemas import EVENT_TYPES, RMATS_COORD_SLOTS, event_key, table as get_table

from conftest import EXPERIMENT_ID, GENOME_BUILD, write_rmats


def test_every_event_type_has_a_coordinate_mapping():
    assert set(RMATS_COORD_SLOTS) == set(EVENT_TYPES)
    assert len(RMATS_COORD_SLOTS["MXE"]) == 8
    assert all(len(slots) == 6 for name, slots in RMATS_COORD_SLOTS.items() if name != "MXE")


def test_event_key_stability(tmp_path):
    """The same event in two different comparisons gets byte-identical keys.

    rMATS assigns IDs per invocation, so a key built from the ID would make cross-experiment
    joins impossible. Here the second run has different IDs for the same loci.
    """
    keys = []
    for offset, run in ((0, "runA"), (500, "runB")):
        rmats_out = write_rmats(tmp_path / run / "rmats_out", id_offset=offset)
        results_dir = tmp_path / run / "results"
        (results_dir / "tables").mkdir(parents=True)
        ctx = RunContext(run_id=run, experiment_id=EXPERIMENT_ID, results_dir=results_dir,
                         genome_build=GENOME_BUILD)
        tables_rmats.write_splicing_tables(ctx, results_dir, rmats_out)
        frame = pd.read_csv(results_dir / get_table("splicing_event").relpath, sep="\t")
        keys.append(frame.sort_values(["counting_mode", "event_key"])["event_key"].tolist())

    assert keys[0] == keys[1]
    assert keys[0], "no events were written"
    # And the IDs really did differ, so the test is not vacuous.
    first = pd.read_csv(tmp_path / "runA" / "results" / get_table("splicing_event").relpath,
                        sep="\t")
    second = pd.read_csv(tmp_path / "runB" / "results" / get_table("splicing_event").relpath,
                         sep="\t")
    assert set(first["rmats_event_id"]) != set(second["rmats_event_id"])


def test_event_key_encodes_type_locus_and_coordinates():
    key = event_key("SE", "chr8", "+", [1000, 1200, 500, 700, 1500, 1700])
    assert key == "SE:chr8:+:1000-1200-500-700-1500-1700"
    # A5SS and A3SS share coordinate column names; the type prefix keeps them distinct.
    assert event_key("A5SS", "chr1", "+", [1, 2, 3, 4, 5, 6]) != \
           event_key("A3SS", "chr1", "+", [1, 2, 3, 4, 5, 6])


def test_event_key_refuses_missing_coordinates():
    with pytest.raises(ValueError):
        event_key("MXE", "chr1", "+", [1, 2, 3, 4, 5, 6])


def test_unused_coordinate_slots_are_na(results_dir):
    frame = pd.read_csv(results_dir / get_table("splicing_event").relpath, sep="\t",
                        dtype=str, keep_default_na=False)
    se_rows = frame[frame["event_type"] == "SE"]
    assert (se_rows["coord_7"] == "NA").all()
    assert (se_rows["coord_8"] == "NA").all()
    mxe_rows = frame[frame["event_type"] == "MXE"]
    assert (mxe_rows["coord_8"] != "NA").all()


def test_replicates_are_exploded_per_sample_group(results_dir):
    events = pd.read_csv(results_dir / get_table("splicing_event").relpath, sep="\t")
    replicates = pd.read_csv(
        results_dir / get_table("splicing_event_replicate").relpath, sep="\t")

    assert set(replicates["sample_group"]) == {"test", "cntl"}
    assert set(replicates["replicate_index"]) == {1, 2}
    # Two groups x two replicates for every event in every counting mode.
    assert len(replicates) == len(events) * 4
    assert set(replicates["event_key"]) <= set(events["event_key"])


def test_group_1_is_the_test_group(results_dir):
    """--b1 is the test BAM list, so IncLevel1 must land on sample_group=test."""
    events = pd.read_csv(results_dir / get_table("splicing_event").relpath, sep="\t")
    replicates = pd.read_csv(
        results_dir / get_table("splicing_event_replicate").relpath, sep="\t")
    row = events[events["gene_symbol"] == "MYC"].iloc[0]
    test_levels = replicates[(replicates["event_key"] == row["event_key"]) &
                             (replicates["counting_mode"] == row["counting_mode"]) &
                             (replicates["sample_group"] == "test")]["inc_level"]
    assert abs(test_levels.mean() - row["inc_level_1_mean"]) < 1e-9
    assert row["inc_level_1_mean"] > row["inc_level_2_mean"]


def test_summary_records_its_thresholds(results_dir):
    summary = pd.read_csv(results_dir / get_table("splicing_summary").relpath, sep="\t")
    assert set(summary["fdr_threshold"]) == {tables_rmats.RMATS_FDR_THRESHOLD}
    assert set(summary["inclusion_diff_threshold"]) == {tables_rmats.RMATS_INC_DIFF_THRESHOLD}

    se_jc = summary[(summary["event_type"] == "SE") & (summary["counting_mode"] == "JC")].iloc[0]
    assert se_jc["total_events"] == 2
    # Only the MYC event clears both FDR <= 0.05 and |dPSI| >= 0.1, and it is test-high.
    assert se_jc["significant_events"] == 1
    assert se_jc["sig_higher_inclusion_test"] == 1
    assert se_jc["sig_higher_inclusion_cntl"] == 0


def test_gene_ids_are_unversioned(results_dir):
    frame = pd.read_csv(results_dir / get_table("splicing_event").relpath, sep="\t")
    assert not frame["gene_id"].astype(str).str.contains(r"\.").any()
