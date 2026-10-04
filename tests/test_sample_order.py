"""Sample order used for the DESeq2 counts matrix, rMATS replicate lists and GSEA classes.

Regression: Mutu_DDX5-KO_2026-09-30 failed in DESeq2 because samples were sorted by the trailing
integer of the name -- the lane in ``DDX5-2X-KO-1_L05`` -- putting replicate 2 (L04) before
replicate 1 (L05), while deseq2_metadata.tsv followed metadata.tsv's replicate_index.
"""

import pandas as pd
import pytest

pytest.importorskip("matplotlib")  # rnaseq.py imports it at module level
rnaseq = pytest.importorskip("rnaseq")


def _metadata(rows):
    return pd.DataFrame(rows, columns=["sample_id", "condition", "replicate_index", "Sample name"])


def test_order_follows_replicate_index_not_lane_number():
    md = _metadata([
        ("X_cntl1", "cntl", 1, "CNTL-2X-KO-1_L02"),
        ("X_cntl2", "cntl", 2, "CNTL-2X-KO-2_L05"),
        ("X_test1", "test", 1, "DDX5-2X-KO-1_L05"),
        ("X_test2", "test", 2, "DDX5-2X-KO-2_L04"),
    ])
    assert rnaseq.ordered_sample_names(md) == [
        "CNTL-2X-KO-1_L02", "CNTL-2X-KO-2_L05", "DDX5-2X-KO-1_L05", "DDX5-2X-KO-2_L04",
    ]


def test_controls_first_regardless_of_row_order():
    md = _metadata([
        ("X_test2", "test", 2, "B2"),
        ("X_cntl1", "cntl", 1, "A1"),
        ("X_test1", "test", 1, "B1"),
        ("X_cntl2", "cntl", 2, "A2"),
    ])
    assert rnaseq.ordered_sample_names(md) == ["A1", "A2", "B1", "B2"]


def test_falls_back_to_natural_sort_without_replicate_index():
    md = pd.DataFrame({
        "Sample name": ["S10", "S2", "C1"],
        "condition": ["test", "test", "cntl"],
    })
    assert rnaseq.ordered_sample_names(md) == ["C1", "S2", "S10"]
