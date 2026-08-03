#!/usr/bin/env python3
"""Controlled vocabularies for metadata fields.

Adding a value here is a one-line change; that friction is intentional. Free-text metadata is
what made the legacy archive unqueryable, so every field except ``notes`` is constrained.
"""

CELL_LINES = frozenset({
    "Mutu", "Akata", "DG75", "SNU719", "YCCEL1", "BCBL1", "Raji", "P3HR1", "HepG2", "HEK293",
})
# How the perturbation was delivered. What was delivered is `perturbation_agent`, and the gene it
# acts on -- when it has one -- is `perturbation_target`.
PERTURBATION_TYPES = frozenset({
    "transfection", "chemical", "BCR-crosslink", "none",
})
# A second perturbation applied alongside the first, e.g. Zta transfection under PAA. Recorded
# separately so "Zta+PAA" stays two queryable facts rather than one opaque string.
CO_TREATMENTS = frozenset({
    "none", "PAA", "siRNA", "CRISPR", "unknown",
})
# What the co-treatment acts on. `PAA_replication` is the odd one out: PAA blocks viral DNA
# replication rather than targeting a host gene, so the "target" is the process it blocks.
CO_TREATMENT_TARGETS = frozenset({
    "PAA_replication", "CNOT1", "CNOT9", "UPF1", "EXOSC3", "NAT10", "TET1", "BRRF1",
})
# Whether the sequenced population was sorted. An unsorted transfection is a mixture of
# transfected and untransfected cells, so this changes what a expression value means.
FACS_PURIFIED_VALUES = frozenset({"yes", "no", "unknown"})
# Which promoter drove the GFP that sorting selected on. pCMV reports transfection; BMRF1p
# reports that the lytic cycle actually started, which is a different population.
FACS_GFP_PROMOTERS = frozenset({"pCMV", "BMRF1p", "none", "unknown"})
ORGANISMS = frozenset({"human", "mouse"})
CONDITIONS = frozenset({"test", "cntl"})
LIBRARY_LAYOUTS = frozenset({"PE", "SE"})
# What was intended at library prep, which is not the same as the `strandedness` column: that one
# is measured per sample by RSeQC after alignment. The two normally agree, and a disagreement is
# a useful flag for a mislabelled sample or the wrong kit, which is why both are kept.
LIBRARY_STRANDEDNESS_VALUES = frozenset({"stranded", "unstranded", "unknown"})
# How RNA was selected before library construction. polyA selection and rRNA depletion see
# different transcriptomes -- non-polyadenylated and unprocessed RNA is present in one and
# absent by construction in the other -- so this is not a comparable axis: a gene that looks
# absent may simply have been selected away. Recording it is what keeps a cross-experiment
# query from reading a library-prep difference as biology. Values track the SRA/ENA
# library_selection vocabulary (`PolyA`, `Inverse rRNA`) without its inconsistent casing.
LIBRARY_SELECTIONS = frozenset({"polyA", "ribodepleted", "unknown"})

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
    "co_treatment": CO_TREATMENTS,
    "co_treatment_target": CO_TREATMENT_TARGETS,
    "facs_purified": FACS_PURIFIED_VALUES,
    "facs_gfp_promoter": FACS_GFP_PROMOTERS,
    "organism": ORGANISMS,
    "condition": CONDITIONS,
    "library_layout": LIBRARY_LAYOUTS,
    "library_selection": LIBRARY_SELECTIONS,
    "library_strandedness": LIBRARY_STRANDEDNESS_VALUES,
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
