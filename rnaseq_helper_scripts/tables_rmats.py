#!/usr/bin/env python3
"""Tasks 9 and 10: rMATS output normalized into queryable tables.

rMATS emits five files per counting mode with five different coordinate column names for the
same positional roles, per-replicate values packed into comma-separated strings, and a summary
that puts JC and JCEC side by side in one row. This module turns that into three tidy tables
keyed by a stable ``event_key``.

Run standalone for reprocessing an existing rMATS directory:

    python tables_rmats.py <rmats_out_dir> <results_dir> <run_id> <comparison_id>
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

from annotation import strip_version
from outputs import RunContext, write_table
from schemas import (COUNTING_MODES, EVENT_TYPES, NA, N_COORD_SLOTS, RMATS_COORD_SLOTS,
                     RMATS_FDR_THRESHOLD, RMATS_INC_DIFF_THRESHOLD, SPLICING_EVENT,
                     SPLICING_EVENT_REPLICATE, SPLICING_SUMMARY, event_key)

# rMATS is invoked with --b1 = test BAMs and --b2 = control BAMs, so SAMPLE_1/IncLevel1 is the
# test group. Changing the invocation without changing this mapping would silently invert every
# inclusion difference in the warehouse.
GROUP_COLUMNS = {
    "test": {"inclusion": "IJC_SAMPLE_1", "skipping": "SJC_SAMPLE_1", "inc_level": "IncLevel1"},
    "cntl": {"inclusion": "IJC_SAMPLE_2", "skipping": "SJC_SAMPLE_2", "inc_level": "IncLevel2"},
}


def mats_path(rmats_out_dir: Path, event_type: str, counting_mode: str) -> Path:
    return Path(rmats_out_dir) / f"{event_type}.MATS.{counting_mode}.txt"


def _unquote(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.strip('"')


def _split_cell(cell) -> list[str]:
    if cell is None or (isinstance(cell, float) and pd.isna(cell)):
        return []
    return [part.strip() for part in str(cell).split(",")]


def _numeric_or_none(value: str):
    if value in ("", "NA", "nan", "NaN", "null", "None"):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _replicate_mean(cell):
    """Mean inclusion level across replicates, ignoring the replicates rMATS could not score."""
    values = [v for v in (_numeric_or_none(x) for x in _split_cell(cell)) if v is not None]
    return sum(values) / len(values) if values else None


def read_event_table(rmats_out_dir: Path, event_type: str, counting_mode: str) -> pd.DataFrame:
    """Load one rMATS file and map it onto the canonical schema."""
    path = mats_path(rmats_out_dir, event_type, counting_mode)
    frame = pd.read_csv(path, sep="\t", dtype=str)

    slots = RMATS_COORD_SLOTS[event_type]
    missing = [column for column in slots if column not in frame.columns]
    if missing:
        raise ValueError(
            f"{path} is missing coordinate column(s) {missing}; the canonical slot mapping in "
            f"schemas.RMATS_COORD_SLOTS does not match this rMATS version."
        )

    out = pd.DataFrame(index=frame.index)
    out["event_type"] = event_type
    out["counting_mode"] = counting_mode
    out["rmats_event_id"] = frame["ID"]
    out["gene_id"] = _unquote(frame["GeneID"]).map(strip_version)
    out["gene_symbol"] = _unquote(frame["geneSymbol"])
    out["chr"] = _unquote(frame["chr"])
    out["strand"] = _unquote(frame["strand"])

    for index in range(1, N_COORD_SLOTS + 1):
        column = f"coord_{index}"
        out[column] = frame[slots[index - 1]] if index <= len(slots) else None

    out["inc_form_len"] = frame.get("IncFormLen")
    out["skip_form_len"] = frame.get("SkipFormLen")
    out["pvalue"] = frame.get("PValue")
    out["fdr"] = frame.get("FDR")
    out["inc_level_1_mean"] = frame["IncLevel1"].map(_replicate_mean)
    out["inc_level_2_mean"] = frame["IncLevel2"].map(_replicate_mean)
    out["inc_level_difference"] = frame.get("IncLevelDifference")

    coord_columns = [out[f"coord_{i}"] for i in range(1, len(slots) + 1)]
    out["event_key"] = [
        event_key(event_type, chrom, strand, list(coords))
        for chrom, strand, *coords in zip(out["chr"], out["strand"], *coord_columns)
    ]

    # Keep the raw per-replicate strings alongside so the replicate table can be exploded from
    # the same parse rather than re-reading the file. write_table drops them again because it
    # only writes columns the schema declares.
    for group, columns in GROUP_COLUMNS.items():
        for role, column in columns.items():
            out[f"rep_{group}_{role}"] = frame.get(column)
    return out


def _explode_replicates(events: pd.DataFrame, ctx: RunContext) -> list[dict]:
    rows: list[dict] = []
    for event in events.itertuples(index=False):
        for group in GROUP_COLUMNS:
            inclusion = _split_cell(getattr(event, f"rep_{group}_inclusion"))
            skipping = _split_cell(getattr(event, f"rep_{group}_skipping"))
            inc_levels = _split_cell(getattr(event, f"rep_{group}_inc_level"))
            n_replicates = max(len(inclusion), len(skipping), len(inc_levels))
            for index in range(n_replicates):
                rows.append({
                    "run_id": ctx.run_id,
                    "comparison_id": ctx.comparison_id,
                    "counting_mode": event.counting_mode,
                    "event_key": event.event_key,
                    "sample_group": group,
                    "replicate_index": index + 1,
                    "inclusion_junction_count":
                        inclusion[index] if index < len(inclusion) else None,
                    "skipping_junction_count":
                        skipping[index] if index < len(skipping) else None,
                    "inc_level":
                        _numeric_or_none(inc_levels[index]) if index < len(inc_levels) else None,
                })
    return rows


def _summarize(events: pd.DataFrame, ctx: RunContext) -> list[dict]:
    """Task 10. Significance is applied here, with the thresholds recorded in the table.

    rMATS' own summary.txt leaves its cut-offs implicit, which is what makes the archived
    counts uninterpretable.
    """
    rows = []
    fdr = pd.to_numeric(events["fdr"], errors="coerce")
    difference = pd.to_numeric(events["inc_level_difference"], errors="coerce")
    significant = (fdr <= RMATS_FDR_THRESHOLD) & (difference.abs() >= RMATS_INC_DIFF_THRESHOLD)

    for (event_type, counting_mode), group in events.groupby(
            ["event_type", "counting_mode"], sort=False):
        mask = significant.loc[group.index]
        group_difference = difference.loc[group.index]
        rows.append({
            "run_id": ctx.run_id,
            "comparison_id": ctx.comparison_id,
            "event_type": event_type,
            "counting_mode": counting_mode,
            "total_events": len(group),
            "significant_events": int(mask.sum()),
            "sig_higher_inclusion_test": int((mask & (group_difference > 0)).sum()),
            "sig_higher_inclusion_cntl": int((mask & (group_difference < 0)).sum()),
            "fdr_threshold": RMATS_FDR_THRESHOLD,
            "inclusion_diff_threshold": RMATS_INC_DIFF_THRESHOLD,
        })
    return rows


def _cross_check_summary(rmats_out_dir: Path, summary_rows: list[dict]) -> None:
    """Compare our totals against rMATS' summary.txt and warn on disagreement."""
    path = Path(rmats_out_dir) / "summary.txt"
    if not path.is_file():
        return
    try:
        reference = pd.read_csv(path, sep="\t")
    except (OSError, ValueError):
        return
    if "EventType" not in reference.columns:
        return
    totals = {str(row.EventType): row for row in reference.itertuples(index=False)}
    for row in summary_rows:
        expected = totals.get(row["event_type"])
        if expected is None:
            continue
        column = f"TotalEvents{row['counting_mode']}"
        if hasattr(expected, column) and int(getattr(expected, column)) != row["total_events"]:
            print(f"WARNING: {row['event_type']} {row['counting_mode']} total is "
                  f"{row['total_events']} but summary.txt says {getattr(expected, column)}.")


def write_splicing_tables(ctx: RunContext, results_dir: Path,
                          rmats_out_dir: Path) -> dict:
    """Write splicing_event, splicing_event_replicate and splicing_summary."""
    rmats_out_dir = Path(rmats_out_dir)
    frames = []
    missing = []
    for counting_mode in COUNTING_MODES:
        for event_type in EVENT_TYPES:
            path = mats_path(rmats_out_dir, event_type, counting_mode)
            if not path.is_file():
                missing.append(path.name)
                continue
            frames.append(read_event_table(rmats_out_dir, event_type, counting_mode))

    if not frames:
        print(f"WARNING: no rMATS tables found in {rmats_out_dir}. Writing splicing tables "
              "with headers only.")
        for name, spec in (("splicing_event", SPLICING_EVENT),
                           ("splicing_event_replicate", SPLICING_EVENT_REPLICATE),
                           ("splicing_summary", SPLICING_SUMMARY)):
            write_table(results_dir, name, pd.DataFrame(columns=spec.column_names))
        return {"status": "skipped", "splicing_event": 0, "splicing_event_replicate": 0,
                "missing_files": missing}

    events = pd.concat(frames, ignore_index=True)
    events["run_id"] = ctx.run_id
    events["comparison_id"] = ctx.comparison_id

    replicates = pd.DataFrame(_explode_replicates(events, ctx),
                              columns=SPLICING_EVENT_REPLICATE.column_names)
    summary_rows = _summarize(events, ctx)
    _cross_check_summary(rmats_out_dir, summary_rows)

    write_table(results_dir, "splicing_event", events[SPLICING_EVENT.column_names])
    write_table(results_dir, "splicing_event_replicate", replicates)
    write_table(results_dir, "splicing_summary",
                pd.DataFrame(summary_rows, columns=SPLICING_SUMMARY.column_names))

    duplicated = events.duplicated(subset=["counting_mode", "event_key"]).sum()
    if duplicated:
        print(f"WARNING: {duplicated} duplicate (counting_mode, event_key) pairs in rMATS "
              "output; event_key is meant to be unique within a counting mode.")
    if missing:
        print(f"WARNING: rMATS did not produce {', '.join(missing)}.")

    return {
        "status": "partial" if missing else "ok",
        "splicing_event": int(len(events)),
        "splicing_event_replicate": int(len(replicates)),
        "distinct_event_keys": int(events["event_key"].nunique()),
        "missing_files": missing,
    }


def main(argv: list[str]) -> int:
    if len(argv) != 4:
        print(__doc__)
        return 1
    rmats_out_dir, results_dir, run_id, comparison_id = argv
    experiment_id, _, groups = comparison_id.partition("__")
    test_group, _, cntl_group = groups.partition("_vs_")
    ctx = RunContext(
        run_id=run_id,
        experiment_id=experiment_id,
        results_dir=Path(results_dir),
        genome_build=NA,
        test_group=test_group or "test",
        cntl_group=cntl_group or "cntl",
    )
    result = write_splicing_tables(ctx, Path(results_dir), Path(rmats_out_dir))
    print(result)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
