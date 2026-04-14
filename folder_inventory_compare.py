#!/usr/bin/env python3
"""
Inventory a directory tree to TSV (relative path + human-readable size) and
compare two such inventories or two folders end-to-end.

Compare mode (default) normalizes relative paths before diffing: run-specific
path segments (e.g. GSEA ``apr13``, ISO dates, UUIDs, ``tmp…`` temp dirs) are
replaced with stable placeholders so two outputs from different run dates can
still match. Raw paths in written TSVs are unchanged. Use ``--no-normalize`` for a strict path-for-path comparison.

Use ``--strip-all-numerics`` to remove every digit from the compare-key path
(after optional normalization), so e.g. GSEA run id folders align across runs.
"""

from __future__ import annotations

import argparse
import re
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path


TSV_HEADER = "relative_path\tsize"

# Placeholders for path normalization (compare step only; inventories stay raw).
_NORM_MON_DAY = "__MON_DAY__"
_NORM_ISO_DATE = "__ISO_DATE__"
_NORM_UUID = "__UUID__"
_NORM_TMPDIR = "__TMPDIR__"

# Whole path segment: GSEA-style run date folder (e.g. apr13, jan01).
_RE_SEG_MON_DAY = re.compile(r"(?i)^(jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec)\d{1,2}$")
# Whole segment: macOS/Linux style random temp dir name (e.g. tmpoagw1kfc).
_RE_SEG_TMP = re.compile(r"(?i)^tmp[a-z0-9]{6,}$")
# Substrings within a segment (filenames that embed dates / hashes).
_RE_UUID = re.compile(
    r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}",
    re.I,
)
_RE_ISO_DATE = re.compile(r"\d{4}[-_]\d{2}[-_]\d{2}")
# Month+day embedded in a longer name (not just whole segment).
_RE_EMBED_MON_DAY = re.compile(
    r"(?i)(?<![0-9a-z])(jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec)\d{1,2}(?![0-9a-z])"
)


def normalize_rel_path(rel: str) -> str:
    """
    Strip common run-specific tokens from a POSIX relative path for comparison.

    Applied only when comparing inventories; written TSVs are unchanged.
    """
    parts = rel.split("/")
    out: list[str] = []
    for seg in parts:
        s = seg
        if _RE_SEG_TMP.match(s):
            s = _NORM_TMPDIR
        elif _RE_SEG_MON_DAY.match(s):
            s = _NORM_MON_DAY
        else:
            s = _RE_UUID.sub(_NORM_UUID, s)
            s = _RE_ISO_DATE.sub(_NORM_ISO_DATE, s)
            s = _RE_EMBED_MON_DAY.sub(_NORM_MON_DAY, s)
        out.append(s)
    return "/".join(out)


def _strip_all_numerics_from_path(s: str) -> str:
    """Remove ASCII digits from *s*, then collapse repeated slashes."""
    out = "".join(c for c in s if not c.isdigit())
    out = re.sub(r"/+", "/", out)
    return out


def _path_key(
    path: str, *, normalize_paths: bool, strip_all_numerics: bool = False
) -> str:
    s = normalize_rel_path(path) if normalize_paths else path
    if strip_all_numerics:
        s = _strip_all_numerics_from_path(s)
    return s


def _path_key_counter(
    path_to_size: dict[str, str],
    *,
    normalize_paths: bool,
    strip_all_numerics: bool = False,
) -> Counter[str]:
    c: Counter[str] = Counter()
    for path in path_to_size:
        c[
            _path_key(
                path,
                normalize_paths=normalize_paths,
                strip_all_numerics=strip_all_numerics,
            )
        ] += 1
    return c


def _entries_by_path_key(
    path_to_size: dict[str, str],
    *,
    normalize_paths: bool,
    strip_all_numerics: bool = False,
) -> dict[str, list[tuple[str, str]]]:
    """Map compare-key -> list of (raw_relative_path, size_str), sorted by raw path."""
    d: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for raw, sz in path_to_size.items():
        d[
            _path_key(
                raw,
                normalize_paths=normalize_paths,
                strip_all_numerics=strip_all_numerics,
            )
        ].append((raw, sz))
    for k in d:
        d[k].sort(key=lambda t: t[0])
    return dict(d)


def _normalization_conflicts(
    path_to_size: dict[str, str],
    *,
    normalize_paths: bool,
    strip_all_numerics: bool = False,
) -> list[str]:
    """Same compare-key but different reported sizes within one inventory."""
    norm_to_sizes: dict[str, set[str]] = defaultdict(set)
    norm_to_raw: dict[str, list[str]] = defaultdict(list)
    for raw, sz in path_to_size.items():
        n = _path_key(
            raw,
            normalize_paths=normalize_paths,
            strip_all_numerics=strip_all_numerics,
        )
        norm_to_sizes[n].add(sz)
        norm_to_raw[n].append(raw)
    issues: list[str] = []
    for n, sizes in norm_to_sizes.items():
        if len(sizes) > 1:
            pairs = sorted((r, path_to_size[r]) for r in norm_to_raw[n])
            detail = "; ".join(f"{r!r} -> {sz!r}" for r, sz in pairs)
            issues.append(
                f"Inconsistent sizes for compare key {n!r}: {detail}"
            )
    return issues


def human_readable_size(size_bytes: int) -> str:
    """Format byte size as e.g. ``1.1 Gb``, ``512.0 Kb``."""
    if size_bytes < 0:
        raise ValueError("size_bytes must be non-negative")
    if size_bytes == 0:
        return "0 B"
    kb = 1024
    if size_bytes >= kb**4:
        return f"{size_bytes / kb**4:.1f} Tb"
    if size_bytes >= kb**3:
        return f"{size_bytes / kb**3:.1f} Gb"
    if size_bytes >= kb**2:
        return f"{size_bytes / kb**2:.1f} Mb"
    if size_bytes >= kb:
        return f"{size_bytes / kb:.1f} Kb"
    return f"{size_bytes} B"


def parse_human_readable_size(s: str) -> int | None:
    """
    Parse strings produced by :func:`human_readable_size` back to whole bytes.

    Returns None if the string cannot be parsed.
    """
    s2 = s.strip()
    if not s2:
        return None
    parts = s2.split()
    if len(parts) != 2:
        return None
    val_s, unit_s = parts[0], parts[1]
    try:
        val = float(val_s)
    except ValueError:
        return None
    u = unit_s.lower()
    mult = {"b": 1, "kb": 1024, "mb": 1024**2, "gb": 1024**3, "tb": 1024**4}
    if u not in mult:
        return None
    return int(round(val * mult[u]))


def format_signed_size_delta(delta_bytes: int) -> str:
    """Human-readable signed delta, e.g. ``+1.2 Kb``, ``-0.5 Mb``."""
    if delta_bytes == 0:
        return "0 B"
    sign = "+" if delta_bytes > 0 else "-"
    mag = human_readable_size(abs(delta_bytes))
    return f"{sign}{mag}"


def folder_to_tsv(folder_abs: Path | str, output_tsv: Path | str) -> Path:
    """
    Walk *folder_abs* (must exist and be a directory), write one row per file:
    column 1 = path relative to *folder_abs* (POSIX-style), column 2 = size string.

    Returns the path to the written TSV.
    """
    root = Path(folder_abs).expanduser().resolve(strict=True)
    if not root.is_dir():
        raise NotADirectoryError(f"Not a directory or does not exist: {root}")

    out_path = Path(output_tsv).expanduser().resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[tuple[str, str]] = []
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        try:
            rel = path.relative_to(root).as_posix()
        except ValueError:
            continue
        size_b = path.stat().st_size
        rows.append((rel, human_readable_size(size_b)))

    with open(out_path, "w", encoding="utf-8", newline="\n") as f:
        f.write(TSV_HEADER + "\n")
        for rel, size in rows:
            f.write(f"{rel}\t{size}\n")
    return out_path


def _parse_inventory_tsv(tsv_path: Path) -> dict[str, str]:
    """Load path -> size string; skip blank lines and the header row."""
    text = tsv_path.read_text(encoding="utf-8")
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    if not lines:
        return {}
    start = 0
    if lines[0].replace(" ", "").lower().startswith("relative_path"):
        start = 1
    out: dict[str, str] = {}
    for line in lines[start:]:
        if "\t" not in line:
            raise ValueError(f"Bad line (no tab): {line!r} in {tsv_path}")
        rel, size = line.split("\t", 1)
        out[rel] = size.strip()
    return out


def _write_log(log_path: Path | str, lines: list[str]) -> Path:
    """Write comparison output lines to a .txt log file."""
    out_path = Path(log_path).expanduser().resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8", newline="\n") as f:
        for line in lines:
            f.write(line + "\n")
    return out_path


def compare_inventory_tsvs(
    tsv_a: Path | str,
    tsv_b: Path | str,
    output_log: Path | str,
    *,
    normalize_paths: bool = True,
    strip_all_numerics: bool = False,
) -> bool:
    """
    Compare two TSV files produced by :func:`folder_to_tsv`.

    Two phases: (1) compare path keys only — report ``Only in first`` /
    ``Only in second`` when a key is missing on one side or has different
    multiplicity; (2) for keys present on both with the same count, compare
    sizes pairwise and report ``Size Mismatch:`` with signed human-readable delta
    (first minus second) and the path key.

    When *normalize_paths* is True, keys use :func:`normalize_rel_path`.
    When *strip_all_numerics* is True, all digits are removed from the compare key
    (after optional normalization), then slash runs are collapsed.
    """
    pa, pb = Path(tsv_a), Path(tsv_b)
    a = _parse_inventory_tsv(pa)
    b = _parse_inventory_tsv(pb)
    report_lines: list[str] = []

    report_lines.append(
        f"Path normalization for compare: {'enabled' if normalize_paths else 'disabled'}"
    )
    report_lines.append(
        f"Strip all numerics from compare keys: {'enabled' if strip_all_numerics else 'disabled'}"
    )

    conflict_a = _normalization_conflicts(
        a, normalize_paths=normalize_paths, strip_all_numerics=strip_all_numerics
    )
    conflict_b = _normalization_conflicts(
        b, normalize_paths=normalize_paths, strip_all_numerics=strip_all_numerics
    )
    for line in conflict_a:
        report_lines.append(f"[first inventory] {line}")
    for line in conflict_b:
        report_lines.append(f"[second inventory] {line}")

    n_a = len(a)
    n_b = len(b)
    if n_a != n_b:
        report_lines.append(
            f"Raw file count mismatch: {pa} has {n_a} files, {pb} has {n_b} files."
        )

    cnt_a = _path_key_counter(
        a, normalize_paths=normalize_paths, strip_all_numerics=strip_all_numerics
    )
    cnt_b = _path_key_counter(
        b, normalize_paths=normalize_paths, strip_all_numerics=strip_all_numerics
    )

    report_lines.append("--- Path / name comparison ---")
    only_path_a = cnt_a - cnt_b
    only_path_b = cnt_b - cnt_a
    for path_key, count in sorted(only_path_a.items()):
        report_lines.append(f"Only in first (x{count}): {path_key!r}")
    for path_key, count in sorted(only_path_b.items()):
        report_lines.append(f"Only in second (x{count}): {path_key!r}")
    if not only_path_a and not only_path_b:
        report_lines.append("(No path-only differences.)")

    entries_a = _entries_by_path_key(
        a, normalize_paths=normalize_paths, strip_all_numerics=strip_all_numerics
    )
    entries_b = _entries_by_path_key(
        b, normalize_paths=normalize_paths, strip_all_numerics=strip_all_numerics
    )

    report_lines.append("--- Size comparison (paths appearing on both sides) ---")
    size_mismatch_lines: list[str] = []
    common_keys = set(cnt_a) & set(cnt_b)
    for path_key in sorted(common_keys):
        if cnt_a[path_key] != cnt_b[path_key]:
            continue
        ea = entries_a.get(path_key, [])
        eb = entries_b.get(path_key, [])
        for (raw_a, sa), (raw_b, sb) in zip(ea, eb):
            if sa == sb:
                continue
            ba = parse_human_readable_size(sa)
            bb = parse_human_readable_size(sb)
            if ba is None or bb is None:
                size_mismatch_lines.append(
                    f"Size Mismatch (unparseable sizes; first={sa!r} second={sb!r}): {path_key}"
                )
            else:
                delta = ba - bb
                size_mismatch_lines.append(
                    f"Size Mismatch: {format_signed_size_delta(delta)}  {path_key}"
                )

    if not size_mismatch_lines:
        report_lines.append("(No size differences for shared paths.)")
    else:
        report_lines.extend(size_mismatch_lines)

    has_conflicts = bool(conflict_a or conflict_b)
    has_path_diff = bool(only_path_a or only_path_b)
    has_size_diff = bool(size_mismatch_lines)
    ok = not has_conflicts and not has_path_diff and not has_size_diff
    if ok:
        report_lines.append(
            "Inventories match: same paths and sizes (per current compare-key options)."
        )
    _write_log(output_log, report_lines)
    return ok


def compare_folders(
    folder_a: Path | str,
    folder_b: Path | str,
    output_log: Path | str,
    *,
    normalize_paths: bool = True,
    strip_all_numerics: bool = False,
) -> bool:
    """
    Build TSV inventories for two folders (temporary directory), compare them,
    then remove the temp tree. Returns the same boolean as
    :func:`compare_inventory_tsvs`.
    """
    fa, fb = Path(folder_a), Path(folder_b)
    with tempfile.TemporaryDirectory() as td:
        td_path = Path(td)
        path_a = td_path / "inventory_a.tsv"
        path_b = td_path / "inventory_b.tsv"
        folder_to_tsv(fa, path_a)
        folder_to_tsv(fb, path_b)
        return compare_inventory_tsvs(
            path_a,
            path_b,
            output_log,
            normalize_paths=normalize_paths,
            strip_all_numerics=strip_all_numerics,
        )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Inventory a folder to TSV, or compare two TSVs / two folders."
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_inv = sub.add_parser("inventory", help="Write TSV for one folder")
    p_inv.add_argument(
        "folder",
        type=Path,
        help="Absolute (or resolvable) path to the root folder",
    )
    p_inv.add_argument(
        "output_tsv",
        type=Path,
        help="Path for the output .tsv file",
    )

    p_cmp = sub.add_parser("compare-tsv", help="Compare two inventory TSV files")
    p_cmp.add_argument("tsv_a", type=Path)
    p_cmp.add_argument("tsv_b", type=Path)
    p_cmp.add_argument("output_log", type=Path, help="Path to output .txt comparison log")
    p_cmp.add_argument(
        "--no-normalize",
        action="store_true",
        help="Compare raw relative paths (no date/tmp token stripping)",
    )
    p_cmp.add_argument(
        "--strip-all-numerics",
        action="store_true",
        help="Remove all digits from compare-key paths (after optional normalization)",
    )

    p_dirs = sub.add_parser("compare-dirs", help="Compare two folders (full test)")
    p_dirs.add_argument("folder_a", type=Path)
    p_dirs.add_argument("folder_b", type=Path)
    p_dirs.add_argument("output_log", type=Path, help="Path to output .txt comparison log")
    p_dirs.add_argument(
        "--no-normalize",
        action="store_true",
        help="Compare raw relative paths (no date/tmp token stripping)",
    )
    p_dirs.add_argument(
        "--strip-all-numerics",
        action="store_true",
        help="Remove all digits from compare-key paths (after optional normalization)",
    )

    args = parser.parse_args(argv)

    if args.command == "inventory":
        out = folder_to_tsv(args.folder, args.output_tsv)
        print(f"Wrote {out}", file=sys.stdout)
        return 0

    if args.command == "compare-tsv":
        ok = compare_inventory_tsvs(
            args.tsv_a,
            args.tsv_b,
            args.output_log,
            normalize_paths=not args.no_normalize,
            strip_all_numerics=args.strip_all_numerics,
        )
        print(f"Wrote comparison log: {Path(args.output_log).expanduser().resolve()}")
        return 0 if ok else 1

    if args.command == "compare-dirs":
        ok = compare_folders(
            args.folder_a,
            args.folder_b,
            args.output_log,
            normalize_paths=not args.no_normalize,
            strip_all_numerics=args.strip_all_numerics,
        )
        print(f"Wrote comparison log: {Path(args.output_log).expanduser().resolve()}")
        return 0 if ok else 1

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
