#!/usr/bin/env python3
"""Single source of truth for every machine-readable table the pipeline emits.

The table writers, ``validate_outputs.py`` and ``docs/SCHEMA.md`` are all generated from or
checked against the definitions here, so code and documentation cannot drift apart.

Column order is part of the contract: the warehouse loader reads positionally as a fallback and
the validator compares names *and* order.
"""

from dataclasses import dataclass, field

# Written for every value that cannot be determined. Never an empty string, never 0.
NA = "NA"

PIPELINE_VERSION = "1.4.0"

# Coordinates in every table are 0-based half-open, matching BED. rMATS already emits
# 0-based starts; STAR SJ.out.tab is 1-based and is converted on read.
COORDINATE_CONVENTION = "0-based half-open"

CONDITIONS = ("test", "cntl")
COUNTING_MODES = ("JC", "JCEC")
EVENT_TYPES = ("SE", "A5SS", "A3SS", "MXE", "RI")
STRANDS = ("plus", "minus", "unstranded")
BIGWIG_CONTENTS = ("unique", "uniquemulti")
NORMALIZATIONS = ("raw", "cpm", "rpkm", "bpm", "spike_in")
QC_TOOLS = ("fastqc", "fastp", "rseqc", "star", "kallisto")
TEST_TYPES = ("wald", "lrt")

# rMATS significance thresholds. These are applied by us rather than inherited implicitly from
# rMATS' own summary.txt so that splicing_summary counts are interpretable.
RMATS_FDR_THRESHOLD = 0.05
RMATS_INC_DIFF_THRESHOLD = 0.1

# Canonical coordinate slots per rMATS event type. Explicit by design: the five rMATS files use
# five different coordinate column names for the same positional roles, and event_key stability
# depends entirely on this mapping never being inferred.
RMATS_COORD_SLOTS: dict[str, tuple[str, ...]] = {
    "SE": (
        "exonStart_0base", "exonEnd",
        "upstreamES", "upstreamEE",
        "downstreamES", "downstreamEE",
    ),
    "A5SS": (
        "longExonStart_0base", "longExonEnd",
        "shortES", "shortEE",
        "flankingES", "flankingEE",
    ),
    "A3SS": (
        "longExonStart_0base", "longExonEnd",
        "shortES", "shortEE",
        "flankingES", "flankingEE",
    ),
    "RI": (
        "riExonStart_0base", "riExonEnd",
        "upstreamES", "upstreamEE",
        "downstreamES", "downstreamEE",
    ),
    "MXE": (
        "1stExonStart_0base", "1stExonEnd",
        "2ndExonStart_0base", "2ndExonEnd",
        "upstreamES", "upstreamEE",
        "downstreamES", "downstreamEE",
    ),
}

N_COORD_SLOTS = 8

# gene_id values that are legitimately not Ensembl accessions. Viral contigs are annotated with
# their historical locus names, which is what every downstream lab query uses.
GENE_ID_PATTERN = r"^ENS[A-Z]*G\d+$"
VIRAL_GENE_ID_PREFIXES = ("BHLF", "BMRF", "BZLF", "BRLF", "BALF", "BLLF", "BNRF", "BCRF",
                          "BDLF", "BGLF", "BKRF", "BLRF", "BFRF", "BFLF", "BPLF", "BOLF",
                          "BSRF", "BTRF", "BVLF", "BVRF", "BXLF", "BXRF", "BILF", "BARF",
                          "BaRF", "BcRF", "EBNA", "LMP", "RPMS", "A73", "K1", "K2", "K3",
                          "K4", "K5", "K6", "K7", "K8", "K9", "K10", "K11", "K12", "K14",
                          "K15", "ORF", "vFLIP", "LANA", "PAN", "M1", "M2", "M3", "M4",
                          "M7", "M8", "M9", "M10", "M11", "M12")


@dataclass(frozen=True)
class Column:
    name: str
    dtype: str  # str | int | float | bool | enum
    enum: tuple[str, ...] | None = None
    notes: str = ""

    @property
    def type_label(self) -> str:
        if self.dtype == "enum" and self.enum:
            return "enum: " + " | ".join(f"`{v}`" for v in self.enum)
        return self.dtype


@dataclass(frozen=True)
class Table:
    name: str
    subdir: str
    columns: tuple[Column, ...]
    description: str = ""
    grain: str = ""
    key_columns: tuple[str, ...] = field(default_factory=tuple)

    @property
    def filename(self) -> str:
        return f"{self.name}.tsv"

    @property
    def relpath(self) -> str:
        return f"{self.subdir}/{self.filename}" if self.subdir else self.filename

    @property
    def column_names(self) -> list[str]:
        return [c.name for c in self.columns]

    def column(self, name: str) -> Column:
        for col in self.columns:
            if col.name == name:
                return col
        raise KeyError(f"{self.name} has no column {name!r}")


def _c(name: str, dtype: str, enum: tuple[str, ...] | None = None, notes: str = "") -> Column:
    return Column(name=name, dtype=dtype, enum=enum, notes=notes)


_RUN_ID = _c("run_id", "str", notes="Foreign key to run_manifest.json. Every row carries it.")
_SAMPLE_ID = _c("sample_id", "str", notes="Foreign key to metadata.tsv.")
_COMPARISON_ID = _c("comparison_id", "str", notes="`{experiment_id}__{test_group}_vs_{cntl_group}`.")


def _coord_columns() -> tuple[Column, ...]:
    return tuple(
        _c(f"coord_{i}", "int", notes=f"Positional slot {i}; `NA` when the event type has no slot {i}.")
        for i in range(1, N_COORD_SLOTS + 1)
    )


TABLES: dict[str, Table] = {}


def _register(table: Table) -> Table:
    TABLES[table.name] = table
    return table


DE_GENE = _register(Table(
    name="de_gene",
    subdir="tables",
    description="DESeq2 gene-level statistics. Statistics only — no per-sample expression.",
    grain="one row per gene per comparison",
    key_columns=("run_id", "comparison_id", "gene_id"),
    columns=(
        _RUN_ID,
        _COMPARISON_ID,
        _c("gene_id", "str", notes="Ensembl gene accession, version suffix stripped."),
        _c("gene_symbol", "str"),
        _c("base_mean", "float"),
        _c("log2_fold_change", "float"),
        _c("lfc_se", "float"),
        _c("stat", "float"),
        _c("pvalue", "float"),
        _c("padj", "float", notes="`NA` where DESeq2 independent filtering removed the gene."),
    ),
))

EXPRESSION_TRANSCRIPT = _register(Table(
    name="expression_transcript",
    subdir="tables",
    description="kallisto transcript abundances, long format.",
    grain="one row per transcript per sample",
    key_columns=("run_id", "sample_id", "transcript_id"),
    columns=(
        _RUN_ID,
        _SAMPLE_ID,
        _c("transcript_id", "str"),
        _c("gene_id", "str"),
        _c("length", "int"),
        _c("eff_length", "float"),
        _c("est_counts", "float"),
        _c("tpm", "float"),
    ),
))

EXPRESSION_GENE = _register(Table(
    name="expression_gene",
    subdir="tables",
    description="Gene-level summarization of kallisto abundances, long format. Uses the same "
                "tx2gene mapping as the DESeq2 input so the two agree.",
    grain="one row per gene per sample",
    key_columns=("run_id", "sample_id", "gene_id"),
    columns=(
        _RUN_ID,
        _SAMPLE_ID,
        _c("gene_id", "str"),
        _c("gene_symbol", "str"),
        _c("est_counts", "float"),
        _c("tpm", "float"),
    ),
))

TRANSCRIPT_DE = _register(Table(
    name="transcript_de",
    subdir="tables",
    description="Sleuth transcript-level differential expression, both test types stacked.",
    grain="one row per transcript per comparison per test_type",
    key_columns=("run_id", "comparison_id", "transcript_id", "test_type"),
    columns=(
        _RUN_ID,
        _COMPARISON_ID,
        _c("transcript_id", "str"),
        _c("gene_id", "str"),
        _c("gene_symbol", "str"),
        _c("test_stat", "float", notes="LRT statistic, or the Wald statistic `b / se_b`."),
        _c("pval", "float"),
        _c("qval", "float"),
        _c("b", "float", notes="Beta / effect size. `NA` for LRT."),
        _c("se_b", "float", notes="`NA` for LRT."),
        _c("mean_obs", "float"),
        _c("test_type", "enum", TEST_TYPES),
    ),
))

GENE_SET_ENRICHMENT = _register(Table(
    name="gene_set_enrichment",
    subdir="tables",
    description="GSEA results for every MSigDB collection concatenated into one table.",
    grain="one row per gene set per collection per comparison",
    key_columns=("run_id", "comparison_id", "collection", "gene_set"),
    columns=(
        _RUN_ID,
        _COMPARISON_ID,
        _c("collection", "str", notes="MSigDB collection, e.g. `h.all`, `c2.cp.reactome`."),
        _c("gene_set", "str"),
        _c("size", "int"),
        _c("es", "float"),
        _c("nes", "float"),
        _c("pvalue", "float"),
        _c("fdr_qvalue", "float"),
        _c("fwer_pvalue", "float"),
        _c("rank_at_max", "int"),
        _c("leading_edge_genes", "str",
           notes="Comma-separated gene symbols. `NA` when GSEA wrote no per-set detail file."),
    ),
))

SPLICING_EVENT = _register(Table(
    name="splicing_event",
    subdir="tables",
    description="All five rMATS event types normalized onto one schema, both counting modes.",
    grain="one row per event per counting_mode per comparison",
    key_columns=("run_id", "comparison_id", "counting_mode", "event_key"),
    columns=(
        _RUN_ID,
        _COMPARISON_ID,
        _c("counting_mode", "enum", COUNTING_MODES),
        _c("event_type", "enum", EVENT_TYPES),
        _c("event_key", "str", notes="Stable cross-comparison identifier; see docs/SCHEMA.md."),
        _c("rmats_event_id", "int", notes="rMATS `ID`, for traceability back to artifacts/rmats_raw/."),
        _c("gene_id", "str"),
        _c("gene_symbol", "str"),
        _c("chr", "str"),
        _c("strand", "enum", ("+", "-")),
        *_coord_columns(),
        _c("inc_form_len", "int"),
        _c("skip_form_len", "int"),
        _c("pvalue", "float"),
        _c("fdr", "float"),
        _c("inc_level_1_mean", "float", notes="Mean across test replicates."),
        _c("inc_level_2_mean", "float", notes="Mean across cntl replicates."),
        _c("inc_level_difference", "float"),
    ),
))

SPLICING_EVENT_REPLICATE = _register(Table(
    name="splicing_event_replicate",
    subdir="tables",
    description="rMATS per-replicate junction counts and inclusion levels, exploded out of the "
                "comma-separated strings rMATS packs into single cells.",
    grain="one row per event per counting_mode per replicate",
    key_columns=("run_id", "comparison_id", "counting_mode", "event_key", "sample_group",
                 "replicate_index"),
    columns=(
        _RUN_ID,
        _COMPARISON_ID,
        _c("counting_mode", "enum", COUNTING_MODES),
        _c("event_key", "str", notes="Subset of splicing_event.event_key."),
        _c("sample_group", "enum", CONDITIONS, notes="rMATS `SAMPLE_1` is `test` (`--b1`)."),
        _c("replicate_index", "int", notes="1-based, in BAM list order."),
        _c("inclusion_junction_count", "int"),
        _c("skipping_junction_count", "int"),
        _c("inc_level", "float"),
    ),
))

SPLICING_SUMMARY = _register(Table(
    name="splicing_summary",
    subdir="tables",
    description="Event counts per type, split by counting mode, with the significance "
                "thresholds recorded explicitly.",
    grain="one row per event_type per counting_mode per comparison",
    key_columns=("run_id", "comparison_id", "event_type", "counting_mode"),
    columns=(
        _RUN_ID,
        _COMPARISON_ID,
        _c("event_type", "enum", EVENT_TYPES),
        _c("counting_mode", "enum", COUNTING_MODES),
        _c("total_events", "int"),
        _c("significant_events", "int"),
        _c("sig_higher_inclusion_test", "int"),
        _c("sig_higher_inclusion_cntl", "int"),
        _c("fdr_threshold", "float"),
        _c("inclusion_diff_threshold", "float"),
    ),
))

JUNCTION = _register(Table(
    name="junction",
    subdir="tables",
    description="STAR SJ.out.tab parsed directly, replacing the strand-split BED files.",
    grain="one row per junction per sample",
    key_columns=("run_id", "sample_id", "chr", "start", "end"),
    columns=(
        _RUN_ID,
        _SAMPLE_ID,
        _c("chr", "str"),
        _c("start", "int", notes="Intron start, converted from STAR's 1-based value to 0-based."),
        _c("end", "int", notes="Intron end, exclusive."),
        _c("strand", "enum", ("+", "-", "."), notes="STAR 1/2/0."),
        _c("intron_motif", "int", notes="0 non-canonical, 1 GT/AG, 2 CT/AC, 3 GC/AG, 4 CT/GC, 5 AT/AC, 6 GT/AT."),
        _c("is_annotated", "bool"),
        _c("n_uniquely_mapped", "int"),
        _c("n_multi_mapped", "int"),
        _c("max_overhang", "int"),
    ),
))

BIGWIG_MANIFEST = _register(Table(
    name="bigwig_manifest",
    subdir="tables",
    description="Authoritative description of every bigWig. Filenames are cosmetic; this is not.",
    grain="one row per bigWig file",
    key_columns=("run_id", "file_path"),
    columns=(
        _RUN_ID,
        _SAMPLE_ID,
        _c("file_path", "str", notes="Relative to results_dir."),
        _c("content", "enum", BIGWIG_CONTENTS),
        _c("strand", "enum", STRANDS),
        _c("genome_build", "str"),
        _c("normalization", "enum", NORMALIZATIONS),
        _c("scale_factor", "float", notes="Factor actually applied to the raw STAR wiggle values."),
        _c("sha256", "str"),
        _c("bytes", "int"),
    ),
))

SIGNAL_OVER_GENE = _register(Table(
    name="signal_over_gene",
    subdir="tables",
    description="Precomputed per-feature bigWig aggregates, so signal questions do not require "
                "reading 243 GB of bigWig.",
    grain="one row per feature per strand per sample per feature_set",
    key_columns=("run_id", "sample_id", "strand", "gene_id", "feature_set"),
    columns=(
        _RUN_ID,
        _SAMPLE_ID,
        _c("strand", "enum", STRANDS),
        _c("gene_id", "str"),
        _c("feature_set", "str", notes="`gene_body` today; extensible (`promoter_tss_500`, `exon`, `dog_10kb`)."),
        _c("mean_coverage", "float", notes="sum_coverage / feature_length."),
        _c("max_coverage", "float"),
        _c("sum_coverage", "float"),
        _c("covered_bases", "int"),
        _c("feature_length", "int"),
    ),
))

QC_METRIC = _register(Table(
    name="qc_metric",
    subdir="tables",
    description="All QC output in long format, so a new metric never changes the schema.",
    grain="one row per metric per sample per tool",
    key_columns=("run_id", "sample_id", "tool", "metric_name"),
    columns=(
        _RUN_ID,
        _SAMPLE_ID,
        _c("tool", "enum", QC_TOOLS),
        _c("metric_name", "str"),
        _c("metric_value", "float", notes="`NA` for non-numeric metrics."),
        _c("metric_value_str", "str", notes="`NA` for numeric metrics."),
    ),
))

# metadata.tsv is not written by a table writer (metadata.py owns it) and it keeps legacy
# columns for backwards compatibility, so the validator requires these columns to be present
# rather than requiring an exact match.
METADATA = Table(
    name="metadata",
    subdir="",
    description="One row per FASTQ sample. Produced by metadata.py and copied into results_dir.",
    grain="one row per sample",
    key_columns=("sample_id",),
    columns=(
        _c("sample_id", "str", notes="`{experiment_id}_{condition}{replicate_index}`."),
        _c("experiment_id", "str",
           notes="Parent directory name, `Model_Experiment`. The first underscore is the only "
                 "structural one. Must be unique across the archive, since `sample_id` and "
                 "`comparison_id` are derived from it."),
        _c("condition", "enum", CONDITIONS),
        _c("replicate_index", "int", notes="1-based, stable ordering by filename."),
        _c("cell_line", "str", notes="Controlled vocabulary, see vocab.py."),
        _c("organism", "enum", ("human", "mouse")),
        _c("genome_build", "str", notes="Must match a directory under reference_dir."),
        _c("perturbation_type", "enum",
           ("transfection", "siRNA", "drug", "BCR-crosslink", "none"),
           notes="How the perturbation was delivered."),
        _c("induced_program", "enum",
           ("lytic_reactivation", "latency", "none", "unknown"),
           notes="The biological program the perturbation was meant to induce, independent of "
                 "how. `unknown` is a missing record; `none` is a recorded absence. Describes "
                 "the experiment, so unlike `perturbation_*` it is set on control rows too."),
        _c("perturbation_target", "str", notes="`NA` for controls."),
        _c("perturbation_dose", "str", notes="`NA` if not applicable."),
        _c("timepoint_hours", "float", notes="`NA` if not applicable."),
        _c("library_layout", "enum", ("PE", "SE")),
        _c("strandedness", "str", notes="Filled in by the RSeQC step; `NA` until then."),
        _c("fastq_r1", "str", notes="Absolute path, as mounted."),
        _c("fastq_r2", "str", notes="`NA` for SE."),
        _c("fastq_r1_md5", "str"),
        _c("fastq_r2_md5", "str", notes="`NA` for SE."),
        _c("investigator", "str"),
        _c("sequencing_run_date", "str", notes="ISO 8601 date, or `NA`."),
        _c("notes", "str", notes="Free text. The only free-text field."),
    ),
)

# Columns rnaseq.py still reads by their historical names. Kept in metadata.tsv so the
# pipeline internals do not have to change in the same commit as the schema.
METADATA_LEGACY_COLUMNS = (
    "Path Read 1", "Path Read 2", "Species", "Sample name", "Condition", "Control?",
    "ConditionReplicate", "Experiment name",
)

# Tables every successful run must contain, in the order docs list them.
REQUIRED_TABLES = (
    "de_gene",
    "expression_gene",
    "expression_transcript",
    "transcript_de",
    "gene_set_enrichment",
    "splicing_event",
    "splicing_event_replicate",
    "splicing_summary",
    "junction",
    "bigwig_manifest",
    "signal_over_gene",
    "qc_metric",
)


def table(name: str) -> Table:
    try:
        return TABLES[name]
    except KeyError:
        raise KeyError(f"unknown table {name!r}; known: {sorted(TABLES)}") from None


def event_key(event_type: str, chrom: str, strand: str, coords: list) -> str:
    """Stable identifier for "the same splicing event" across comparisons and runs.

    Two comparisons that detect the same event must produce byte-identical keys, so this is
    built only from the event type, locus and canonical coordinate slots -- never from rMATS'
    per-run ``ID``, which is assigned per invocation.
    """
    if event_type not in RMATS_COORD_SLOTS:
        raise ValueError(f"unknown rMATS event type {event_type!r}")
    n_slots = len(RMATS_COORD_SLOTS[event_type])
    used = coords[:n_slots]
    if len(used) != n_slots or any(c is None or c == NA for c in used):
        raise ValueError(
            f"{event_type} event_key needs {n_slots} coordinates, got {used!r}"
        )
    return f"{event_type}:{chrom}:{strand}:" + "-".join(str(int(c)) for c in used)
