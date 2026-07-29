#!/usr/bin/env python3
"""Controlled vocabularies for metadata fields.

Adding a value here is a one-line change; that friction is intentional. Free-text metadata is
what made the legacy archive unqueryable, so every field except ``notes`` is constrained.
"""

CELL_LINES = frozenset({
    "Mutu", "Akata", "DG75", "HepG2", "Raji", "SNU719", "BCBL1", "HEK293",
})
PERTURBATION_TYPES = frozenset({
    "transfection", "siRNA", "drug", "BCR-crosslink", "none",
})
ORGANISMS = frozenset({"human", "mouse"})
CONDITIONS = frozenset({"test", "cntl"})
LIBRARY_LAYOUTS = frozenset({"PE", "SE"})

# Genome build -> organism, so metadata.py can derive organism instead of asking for it.
GENOME_BUILD_ORGANISM = {
    "hg38": "human",
    "hg38plusAkataInverted": "human",
    "hg38plusKSHV": "human",
    "hg38plusKSHVALT": "human",
    "mm39": "mouse",
    "mm39plusMHV68": "mouse",
}

# Known-strand control genes used to assert that bigWig strand labelling is not inverted.
# GAPDH/Gapdh are on the plus strand, ACTB/Actb on the minus strand in both builds.
STRAND_CONTROL_GENES = {
    "human": {"plus": ("GAPDH",), "minus": ("ACTB",)},
    "mouse": {"plus": ("Gapdh",), "minus": ("Actb",)},
}

VOCABULARIES = {
    "cell_line": CELL_LINES,
    "perturbation_type": PERTURBATION_TYPES,
    "organism": ORGANISMS,
    "condition": CONDITIONS,
    "library_layout": LIBRARY_LAYOUTS,
}


class VocabError(ValueError):
    """Raised when a metadata value is outside its controlled vocabulary."""


def validate(field: str, value: str) -> str:
    """Return ``value`` if it is allowed for ``field``, else raise with the allowed values."""
    try:
        allowed = VOCABULARIES[field]
    except KeyError:
        raise VocabError(f"no vocabulary defined for field {field!r}") from None
    if value not in allowed:
        raise VocabError(
            f"invalid {field}: {value!r}\n"
            f"allowed values: {', '.join(sorted(allowed))}\n"
            f"to add one, edit CELL_LINES/PERTURBATION_TYPES/... in "
            f"rnaseq_helper_scripts/vocab.py"
        )
    return value


def organism_for_build(genome_build: str) -> str | None:
    return GENOME_BUILD_ORGANISM.get(genome_build)
