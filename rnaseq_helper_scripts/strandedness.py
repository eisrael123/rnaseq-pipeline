#!/usr/bin/env python3
"""Library strandedness: one parse of the RSeQC output, one vocabulary, four consumers.

Previously the RSeQC free-text file was re-parsed by string matching in four places, each with
its own spelling of the answer. Everything now goes through :func:`parse_rseqc` and the
canonical labels below.
"""

from __future__ import annotations

import re
from pathlib import Path

FR_FIRSTSTRAND = "fr-firststrand"   # dUTP: read 1 is antisense to the transcript
FR_SECONDSTRAND = "fr-secondstrand"  # read 1 is sense to the transcript
UNSTRANDED = "unstranded"
UNDETERMINED = "undetermined"

LIBRARY_TYPES = (FR_FIRSTSTRAND, FR_SECONDSTRAND, UNSTRANDED, UNDETERMINED)

# Below this difference between the forward and reverse fractions the library is not stranded.
STRANDED_MARGIN = 0.1
# Above this failure fraction RSeQC could not classify enough reads to trust the call.
MAX_FRACTION_FAILED = 0.5

_FAILED_RE = re.compile(r"failed to determine:\s*([\d.eE+-]+)")
_EXPLAINED_RE = re.compile(r'explained by\s*"([^"]+)":\s*([\d.eE+-]+)')

CANONICAL_LINE = "Library Type:"


def parse_rseqc(text: str) -> dict:
    """Interpret ``infer_experiment.py`` output.

    Returns ``library_type``, ``layout`` and the three fractions RSeQC reports. A pattern whose
    read 1 is on the same strand as the transcript (``1++,1--,2+-,2-+`` for paired data,
    ``++,--`` for single) means a second-strand library.
    """
    fraction_failed = None
    match = _FAILED_RE.search(text)
    if match:
        fraction_failed = float(match.group(1))

    fraction_forward = None
    fraction_reverse = None
    for pattern, value in _EXPLAINED_RE.findall(text):
        fraction = float(value)
        if pattern.startswith("1++") or pattern.startswith("++"):
            fraction_forward = fraction
        else:
            fraction_reverse = fraction

    if "PairEnd" in text:
        layout = "PE"
    elif "SingleEnd" in text:
        layout = "SE"
    else:
        layout = UNDETERMINED

    if fraction_forward is None or fraction_reverse is None:
        library_type = UNDETERMINED
    elif fraction_failed is not None and fraction_failed > MAX_FRACTION_FAILED:
        library_type = UNDETERMINED
    elif abs(fraction_forward - fraction_reverse) < STRANDED_MARGIN:
        library_type = UNSTRANDED
    elif fraction_forward > fraction_reverse:
        library_type = FR_SECONDSTRAND
    else:
        library_type = FR_FIRSTSTRAND

    return {
        "library_type": library_type,
        "layout": layout,
        "fraction_failed": fraction_failed,
        "fraction_forward": fraction_forward,
        "fraction_reverse": fraction_reverse,
    }


def rseqc_path(results_dir: Path, sample_name: str) -> Path:
    return (Path(results_dir) / "rseqc" /
            f"{sample_name}_Aligned.sortedByCoord.out_strandedness.txt")


def read_library_type(results_dir: Path, sample_name: str) -> str:
    path = rseqc_path(results_dir, sample_name)
    if not path.is_file():
        raise FileNotFoundError(f"RSeQC strandedness file missing for {sample_name}: {path}")
    return parse_rseqc(path.read_text())["library_type"]


def kallisto_option(library_type: str) -> str:
    return {
        FR_SECONDSTRAND: "--fr-stranded",
        FR_FIRSTSTRAND: "--rf-stranded",
    }.get(library_type, "")


def rmats_lib_type(library_type: str) -> str:
    return {
        FR_SECONDSTRAND: "fr-secondstrand",
        FR_FIRSTSTRAND: "fr-firststrand",
    }.get(library_type, "fr-unstranded")


# STAR emits two signal tracks per sample, str1 and str2, in read-strand order rather than
# genomic-strand order, so the meaning of str1 depends on the library chemistry. For a
# first-strand (dUTP) library read 1 is antisense to the transcript, so str1 tracks the minus
# strand; for a second-strand library read 1 is sense, so str1 tracks the plus strand. This is
# the same relationship the pipeline has always relied on implicitly -- the legacy code negated
# str1 values for first-strand libraries and str2 values for second-strand ones, i.e. it drew
# whichever track it considered the minus strand below the axis. Making it explicit is the point
# of this table; Task 12's assertion checks it against known-strand control genes at runtime.
STAR_STRAND_LABELS = {
    FR_FIRSTSTRAND: {"str1": "minus", "str2": "plus"},
    FR_SECONDSTRAND: {"str1": "plus", "str2": "minus"},
    UNSTRANDED: {"str1": "unstranded"},
}


def star_strand_label(library_type: str, star_strand: str) -> str:
    """Genomic strand for a STAR ``str1``/``str2`` signal track.

    Raises for ``undetermined`` libraries: guessing here is exactly how the legacy archive
    ended up with unrecoverable strand information.
    """
    if library_type == UNDETERMINED:
        raise ValueError(
            "strandedness is undetermined, so STAR str1/str2 cannot be resolved to a genomic "
            "strand. Refusing to guess -- inspect the RSeQC output for this sample."
        )
    try:
        return STAR_STRAND_LABELS[library_type][star_strand]
    except KeyError:
        raise ValueError(
            f"no strand label for library_type={library_type!r} track={star_strand!r}. "
            f"An unstranded library must not produce a str2 track."
        ) from None
