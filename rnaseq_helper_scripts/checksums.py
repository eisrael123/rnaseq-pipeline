#!/usr/bin/env python3
"""Task 15: checksums over everything the warehouse ingests.

This is how a reprocessed experiment is verified before the legacy originals are deleted, so it
is load-bearing. Output is standard ``sha256sum`` format, relative to ``results_dir``:

    cd <results_dir> && sha256sum -c checksums.sha256
"""

from __future__ import annotations

from pathlib import Path

from outputs import sha256_file

CHECKSUM_FILE = "checksums.sha256"
COVERED_SUBDIRS = ("tables", "artifacts")
IGNORED_NAMES = {".DS_Store"}


def covered_files(results_dir: Path) -> list[Path]:
    results_dir = Path(results_dir)
    files: list[Path] = []
    for subdir in COVERED_SUBDIRS:
        root = results_dir / subdir
        if not root.is_dir():
            continue
        for path in sorted(root.rglob("*")):
            if path.is_file() and path.name not in IGNORED_NAMES:
                files.append(path)
    return files


def write_checksums(results_dir: Path, known: dict[str, str] | None = None) -> Path:
    """Hash every file under ``tables/`` and ``artifacts/``.

    ``known`` lets callers pass hashes they already computed (the bigWig manifest), so a 243 GB
    signal directory is not read twice.
    """
    results_dir = Path(results_dir)
    known = known or {}
    lines = []
    for path in covered_files(results_dir):
        rel = path.relative_to(results_dir).as_posix()
        digest = known.get(rel) or sha256_file(path)
        lines.append(f"{digest}  {rel}")

    target = results_dir / CHECKSUM_FILE
    target.write_text("\n".join(lines) + ("\n" if lines else ""))
    return target


def verify_checksums(results_dir: Path) -> list[str]:
    """Return a list of human-readable problems; empty means everything verified."""
    results_dir = Path(results_dir)
    target = results_dir / CHECKSUM_FILE
    if not target.is_file():
        return [f"{CHECKSUM_FILE} is missing"]

    problems = []
    listed = set()
    for line_number, line in enumerate(target.read_text().splitlines(), start=1):
        if not line.strip():
            continue
        digest, _, rel = line.partition("  ")
        if not rel:
            problems.append(f"{CHECKSUM_FILE}:{line_number}: malformed line")
            continue
        listed.add(rel)
        path = results_dir / rel
        if not path.is_file():
            problems.append(f"{rel}: listed in {CHECKSUM_FILE} but missing")
            continue
        actual = sha256_file(path)
        if actual != digest:
            problems.append(f"{rel}: sha256 mismatch")

    for path in covered_files(results_dir):
        rel = path.relative_to(results_dir).as_posix()
        if rel not in listed:
            problems.append(f"{rel}: not listed in {CHECKSUM_FILE}")
    return problems
