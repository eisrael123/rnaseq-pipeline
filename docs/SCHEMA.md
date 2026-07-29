# Output schema

Generated from `rnaseq_helper_scripts/schemas.py` by
`rnaseq_helper_scripts/generate_schema_doc.py`. Do not edit by hand; edit the schema and
regenerate.

Pipeline version `1.4.0`.

## Conventions

- Every table is a tab-separated file with a single header line and Unix line endings.
- **Column order is part of the contract.** `validate_outputs.py` compares names and order.
- Missing values are the literal string `NA`. Never an empty cell, never `0`.
- Every row carries `run_id`, which joins to `run_manifest.json`.
- Coordinates are **0-based half-open**, matching BED. rMATS already emits
  0-based starts; STAR `SJ.out.tab` is 1-based and is converted when it is read.
- Gene IDs are unversioned Ensembl accessions (`^ENS[A-Z]*G\d+$`). The exceptions
  are ERCC spike-ins and viral loci, which are listed in `docs/gene_id_exceptions.txt`.
- `comparison_id` is `{experiment_id}__{test_group}_vs_{cntl_group}`.


## Directory layout

```
<results_dir>/
├── run_manifest.json
├── metadata.tsv
├── checksums.sha256
├── tables/
├── artifacts/
│   ├── bigwig/
│   ├── bam/
│   ├── star_logs/
│   ├── junction_bed/
│   ├── kallisto/
│   ├── rmats_raw/
│   ├── deseq2/
│   ├── sleuth/
│   ├── gsea/
│   └── qc/
├── reports/
└── logs/
```

`tables/` is the contract surface: everything the warehouse reads lives there and nowhere else. `artifacts/` holds raw tool output for re-analysis, `reports/` the human-readable XLSX, figures and HTML, `logs/` the per-stage logs.

## `run_manifest.json`

Written at the start of the run with `exit_status: running` and rewritten at the end. A run that crashes still leaves a manifest, with `exit_status: failed` and the stage that failed.

| Field | Description |
|---|---|
| `run_id` | `R` + 12 hex characters, derived from experiment, start time and commit. Joins to the `run_id` column of every table. |
| `experiment_id` | Input directory name, `Model_Experiment`. |
| `pipeline_version` | Version of this pipeline's output contract. |
| `git_commit` | Commit the code was run from. |
| `git_dirty` | `true` if the working tree had uncommitted changes. Results from a dirty tree are not reproducible from the commit alone. |
| `docker_image_digest` | Image the run executed in. |
| `docker_image_digest_source` | How the digest was obtained: the runner's environment, `/work/VERSION`, or the local git commit for a bare run. |
| `genome_build` | Reference build name, matching a directory under `reference_dir`. |
| `annotation_version` | GTF version, e.g. `gencode_v44`. Read from an `ANNOTATION_VERSION` marker, the GTF header, or `--annotation-version`. |
| `reference_dir` | Path to the reference directory as mounted. |
| `reference_dir_sha256` | Digest over the reference file inventory, so two runs against different references are distinguishable. |
| `investigator` | Who ran it. |
| `library_layout` | `PE` or `SE`. |
| `run_start_utc` | ISO 8601, UTC. |
| `run_end_utc` | ISO 8601, UTC. `null` while the run is in progress. |
| `exit_status` | `running`, `success` or `failed`. |
| `tool_versions` | Object mapping tool name to version string, captured at runtime rather than hardcoded. |
| `parameters` | Object recording the non-default arguments passed to each tool. |

Added as the run progresses, so absent from a manifest written at startup:

| Field | Description |
|---|---|
| `stage_status` | Object mapping stage name to `ok`, `skipped` or `partial`. |
| `failed_stage` | Stage that raised, present only when `exit_status` is `failed`. |
| `error` | Exception type and message, present only when `exit_status` is `failed`. |
| `validation` | Result of `validate_outputs.py`: `status`, `failures` and `warnings`. |

## Tables

Tables every successful run must contain:

- [`tables/de_gene.tsv`](#tablesde_genetsv)
- [`tables/expression_gene.tsv`](#tablesexpression_genetsv)
- [`tables/expression_transcript.tsv`](#tablesexpression_transcripttsv)
- [`tables/transcript_de.tsv`](#tablestranscript_detsv)
- [`tables/gene_set_enrichment.tsv`](#tablesgene_set_enrichmenttsv)
- [`tables/splicing_event.tsv`](#tablessplicing_eventtsv)
- [`tables/splicing_event_replicate.tsv`](#tablessplicing_event_replicatetsv)
- [`tables/splicing_summary.tsv`](#tablessplicing_summarytsv)
- [`tables/junction.tsv`](#tablesjunctiontsv)
- [`tables/bigwig_manifest.tsv`](#tablesbigwig_manifesttsv)
- [`tables/signal_over_gene.tsv`](#tablessignal_over_genetsv)
- [`tables/qc_metric.tsv`](#tablesqc_metrictsv)

### `tables/de_gene.tsv`

DESeq2 gene-level statistics. Statistics only — no per-sample expression.

**Grain:** one row per gene per comparison

**Key:** `run_id`, `comparison_id`, `gene_id`

| Column | Type | Notes |
|---|---|---|
| `run_id` | str | Foreign key to run_manifest.json. Every row carries it. |
| `comparison_id` | str | `{experiment_id}__{test_group}_vs_{cntl_group}`. |
| `gene_id` | str | Ensembl gene accession, version suffix stripped. |
| `gene_symbol` | str |  |
| `base_mean` | float |  |
| `log2_fold_change` | float |  |
| `lfc_se` | float |  |
| `stat` | float |  |
| `pvalue` | float |  |
| `padj` | float | `NA` where DESeq2 independent filtering removed the gene. |

### `tables/expression_gene.tsv`

Gene-level summarization of kallisto abundances, long format. Uses the same tx2gene mapping as the DESeq2 input so the two agree.

**Grain:** one row per gene per sample

**Key:** `run_id`, `sample_id`, `gene_id`

| Column | Type | Notes |
|---|---|---|
| `run_id` | str | Foreign key to run_manifest.json. Every row carries it. |
| `sample_id` | str | Foreign key to metadata.tsv. |
| `gene_id` | str |  |
| `gene_symbol` | str |  |
| `est_counts` | float |  |
| `tpm` | float |  |

### `tables/expression_transcript.tsv`

kallisto transcript abundances, long format.

**Grain:** one row per transcript per sample

**Key:** `run_id`, `sample_id`, `transcript_id`

| Column | Type | Notes |
|---|---|---|
| `run_id` | str | Foreign key to run_manifest.json. Every row carries it. |
| `sample_id` | str | Foreign key to metadata.tsv. |
| `transcript_id` | str |  |
| `gene_id` | str |  |
| `length` | int |  |
| `eff_length` | float |  |
| `est_counts` | float |  |
| `tpm` | float |  |

### `tables/transcript_de.tsv`

Sleuth transcript-level differential expression, both test types stacked.

**Grain:** one row per transcript per comparison per test_type

**Key:** `run_id`, `comparison_id`, `transcript_id`, `test_type`

| Column | Type | Notes |
|---|---|---|
| `run_id` | str | Foreign key to run_manifest.json. Every row carries it. |
| `comparison_id` | str | `{experiment_id}__{test_group}_vs_{cntl_group}`. |
| `transcript_id` | str |  |
| `gene_id` | str |  |
| `gene_symbol` | str |  |
| `test_stat` | float | LRT statistic, or the Wald statistic `b / se_b`. |
| `pval` | float |  |
| `qval` | float |  |
| `b` | float | Beta / effect size. `NA` for LRT. |
| `se_b` | float | `NA` for LRT. |
| `mean_obs` | float |  |
| `test_type` | enum: `wald` \| `lrt` |  |

### `tables/gene_set_enrichment.tsv`

GSEA results for every MSigDB collection concatenated into one table.

**Grain:** one row per gene set per collection per comparison

**Key:** `run_id`, `comparison_id`, `collection`, `gene_set`

| Column | Type | Notes |
|---|---|---|
| `run_id` | str | Foreign key to run_manifest.json. Every row carries it. |
| `comparison_id` | str | `{experiment_id}__{test_group}_vs_{cntl_group}`. |
| `collection` | str | MSigDB collection, e.g. `h.all`, `c2.cp.reactome`. |
| `gene_set` | str |  |
| `size` | int |  |
| `es` | float |  |
| `nes` | float |  |
| `pvalue` | float |  |
| `fdr_qvalue` | float |  |
| `fwer_pvalue` | float |  |
| `rank_at_max` | int |  |
| `leading_edge_genes` | str | Comma-separated gene symbols. `NA` when GSEA wrote no per-set detail file. |

### `tables/splicing_event.tsv`

All five rMATS event types normalized onto one schema, both counting modes.

**Grain:** one row per event per counting_mode per comparison

**Key:** `run_id`, `comparison_id`, `counting_mode`, `event_key`

| Column | Type | Notes |
|---|---|---|
| `run_id` | str | Foreign key to run_manifest.json. Every row carries it. |
| `comparison_id` | str | `{experiment_id}__{test_group}_vs_{cntl_group}`. |
| `counting_mode` | enum: `JC` \| `JCEC` |  |
| `event_type` | enum: `SE` \| `A5SS` \| `A3SS` \| `MXE` \| `RI` |  |
| `event_key` | str | Stable cross-comparison identifier; see docs/SCHEMA.md. |
| `rmats_event_id` | int | rMATS `ID`, for traceability back to artifacts/rmats_raw/. |
| `gene_id` | str |  |
| `gene_symbol` | str |  |
| `chr` | str |  |
| `strand` | enum: `+` \| `-` |  |
| `coord_1` | int | Positional slot 1; `NA` when the event type has no slot 1. |
| `coord_2` | int | Positional slot 2; `NA` when the event type has no slot 2. |
| `coord_3` | int | Positional slot 3; `NA` when the event type has no slot 3. |
| `coord_4` | int | Positional slot 4; `NA` when the event type has no slot 4. |
| `coord_5` | int | Positional slot 5; `NA` when the event type has no slot 5. |
| `coord_6` | int | Positional slot 6; `NA` when the event type has no slot 6. |
| `coord_7` | int | Positional slot 7; `NA` when the event type has no slot 7. |
| `coord_8` | int | Positional slot 8; `NA` when the event type has no slot 8. |
| `inc_form_len` | int |  |
| `skip_form_len` | int |  |
| `pvalue` | float |  |
| `fdr` | float |  |
| `inc_level_1_mean` | float | Mean across test replicates. |
| `inc_level_2_mean` | float | Mean across cntl replicates. |
| `inc_level_difference` | float |  |

### `tables/splicing_event_replicate.tsv`

rMATS per-replicate junction counts and inclusion levels, exploded out of the comma-separated strings rMATS packs into single cells.

**Grain:** one row per event per counting_mode per replicate

**Key:** `run_id`, `comparison_id`, `counting_mode`, `event_key`, `sample_group`, `replicate_index`

| Column | Type | Notes |
|---|---|---|
| `run_id` | str | Foreign key to run_manifest.json. Every row carries it. |
| `comparison_id` | str | `{experiment_id}__{test_group}_vs_{cntl_group}`. |
| `counting_mode` | enum: `JC` \| `JCEC` |  |
| `event_key` | str | Subset of splicing_event.event_key. |
| `sample_group` | enum: `test` \| `cntl` | rMATS `SAMPLE_1` is `test` (`--b1`). |
| `replicate_index` | int | 1-based, in BAM list order. |
| `inclusion_junction_count` | int |  |
| `skipping_junction_count` | int |  |
| `inc_level` | float |  |

### `tables/splicing_summary.tsv`

Event counts per type, split by counting mode, with the significance thresholds recorded explicitly.

**Grain:** one row per event_type per counting_mode per comparison

**Key:** `run_id`, `comparison_id`, `event_type`, `counting_mode`

| Column | Type | Notes |
|---|---|---|
| `run_id` | str | Foreign key to run_manifest.json. Every row carries it. |
| `comparison_id` | str | `{experiment_id}__{test_group}_vs_{cntl_group}`. |
| `event_type` | enum: `SE` \| `A5SS` \| `A3SS` \| `MXE` \| `RI` |  |
| `counting_mode` | enum: `JC` \| `JCEC` |  |
| `total_events` | int |  |
| `significant_events` | int |  |
| `sig_higher_inclusion_test` | int |  |
| `sig_higher_inclusion_cntl` | int |  |
| `fdr_threshold` | float |  |
| `inclusion_diff_threshold` | float |  |

### `tables/junction.tsv`

STAR SJ.out.tab parsed directly, replacing the strand-split BED files.

**Grain:** one row per junction per sample

**Key:** `run_id`, `sample_id`, `chr`, `start`, `end`

| Column | Type | Notes |
|---|---|---|
| `run_id` | str | Foreign key to run_manifest.json. Every row carries it. |
| `sample_id` | str | Foreign key to metadata.tsv. |
| `chr` | str |  |
| `start` | int | Intron start, converted from STAR's 1-based value to 0-based. |
| `end` | int | Intron end, exclusive. |
| `strand` | enum: `+` \| `-` \| `.` | STAR 1/2/0. |
| `intron_motif` | int | 0 non-canonical, 1 GT/AG, 2 CT/AC, 3 GC/AG, 4 CT/GC, 5 AT/AC, 6 GT/AT. |
| `is_annotated` | bool |  |
| `n_uniquely_mapped` | int |  |
| `n_multi_mapped` | int |  |
| `max_overhang` | int |  |

### `tables/bigwig_manifest.tsv`

Authoritative description of every bigWig. Filenames are cosmetic; this is not.

**Grain:** one row per bigWig file

**Key:** `run_id`, `file_path`

| Column | Type | Notes |
|---|---|---|
| `run_id` | str | Foreign key to run_manifest.json. Every row carries it. |
| `sample_id` | str | Foreign key to metadata.tsv. |
| `file_path` | str | Relative to results_dir. |
| `content` | enum: `unique` \| `uniquemulti` |  |
| `strand` | enum: `plus` \| `minus` \| `unstranded` |  |
| `genome_build` | str |  |
| `normalization` | enum: `raw` \| `cpm` \| `rpkm` \| `bpm` \| `spike_in` |  |
| `scale_factor` | float | Factor actually applied to the raw STAR wiggle values. |
| `sha256` | str |  |
| `bytes` | int |  |

### `tables/signal_over_gene.tsv`

Precomputed per-feature bigWig aggregates, so signal questions do not require reading 243 GB of bigWig.

**Grain:** one row per feature per strand per sample per feature_set

**Key:** `run_id`, `sample_id`, `strand`, `gene_id`, `feature_set`

| Column | Type | Notes |
|---|---|---|
| `run_id` | str | Foreign key to run_manifest.json. Every row carries it. |
| `sample_id` | str | Foreign key to metadata.tsv. |
| `strand` | enum: `plus` \| `minus` \| `unstranded` |  |
| `gene_id` | str |  |
| `feature_set` | str | `gene_body` today; extensible (`promoter_tss_500`, `exon`, `dog_10kb`). |
| `mean_coverage` | float | sum_coverage / feature_length. |
| `max_coverage` | float |  |
| `sum_coverage` | float |  |
| `covered_bases` | int |  |
| `feature_length` | int |  |

### `tables/qc_metric.tsv`

All QC output in long format, so a new metric never changes the schema.

**Grain:** one row per metric per sample per tool

**Key:** `run_id`, `sample_id`, `tool`, `metric_name`

| Column | Type | Notes |
|---|---|---|
| `run_id` | str | Foreign key to run_manifest.json. Every row carries it. |
| `sample_id` | str | Foreign key to metadata.tsv. |
| `tool` | enum: `fastqc` \| `fastp` \| `rseqc` \| `star` \| `kallisto` |  |
| `metric_name` | str |  |
| `metric_value` | float | `NA` for non-numeric metrics. |
| `metric_value_str` | str | `NA` for numeric metrics. |

## Inputs

`metadata.tsv` is produced by `metadata.py`, copied into `<results_dir>/` and updated in place once strandedness is inferred. It additionally carries the legacy columns `Path Read 1`, `Path Read 2`, `Species`, `Sample name`, `Condition`, `Control?`, `ConditionReplicate`, `Experiment name`, which `rnaseq.py` still reads by name, so the validator requires these columns to be present rather than requiring an exact match.

### `metadata.tsv`

One row per FASTQ sample. Produced by metadata.py and copied into results_dir.

**Grain:** one row per sample

**Key:** `sample_id`

| Column | Type | Notes |
|---|---|---|
| `sample_id` | str | `{experiment_id}_{condition}{replicate_index}`. |
| `experiment_id` | str | Parent directory name, `Model_Experiment`. The first underscore is the only structural one. Must be unique across the archive, since `sample_id` and `comparison_id` are derived from it. |
| `condition` | enum: `test` \| `cntl` |  |
| `replicate_index` | int | 1-based, stable ordering by filename. |
| `cell_line` | str | Controlled vocabulary, see vocab.py. |
| `organism` | enum: `human` \| `mouse` |  |
| `genome_build` | str | Must match a directory under reference_dir. |
| `perturbation_type` | enum: `transfection` \| `siRNA` \| `drug` \| `BCR-crosslink` \| `none` |  |
| `perturbation_target` | str | `NA` for controls. |
| `perturbation_dose` | str | `NA` if not applicable. |
| `timepoint_hours` | float | `NA` if not applicable. |
| `library_layout` | enum: `PE` \| `SE` |  |
| `strandedness` | str | Filled in by the RSeQC step; `NA` until then. |
| `fastq_r1` | str | Absolute path, as mounted. |
| `fastq_r2` | str | `NA` for SE. |
| `fastq_r1_md5` | str |  |
| `fastq_r2_md5` | str | `NA` for SE. |
| `investigator` | str |  |
| `sequencing_run_date` | str | ISO 8601 date, or `NA`. |
| `notes` | str | Free text. The only free-text field. |

## `event_key`

rMATS assigns its `ID` column per invocation, so it cannot identify an event across comparisons. `event_key` is built only from values intrinsic to the locus:

```
<event_type>:<chr>:<strand>:<coord_1>-<coord_2>-...-<coord_n>
```

The coordinate slots per event type are fixed, because the five rMATS files use different column names for the same positional roles:

| Event type | Slots, in order |
|---|---|
| `SE` | `exonStart_0base`, `exonEnd`, `upstreamES`, `upstreamEE`, `downstreamES`, `downstreamEE` |
| `A5SS` | `longExonStart_0base`, `longExonEnd`, `shortES`, `shortEE`, `flankingES`, `flankingEE` |
| `A3SS` | `longExonStart_0base`, `longExonEnd`, `shortES`, `shortEE`, `flankingES`, `flankingEE` |
| `RI` | `riExonStart_0base`, `riExonEnd`, `upstreamES`, `upstreamEE`, `downstreamES`, `downstreamEE` |
| `MXE` | `1stExonStart_0base`, `1stExonEnd`, `2ndExonStart_0base`, `2ndExonEnd`, `upstreamES`, `upstreamEE`, `downstreamES`, `downstreamEE` |

Unused slots up to `coord_8` are `NA`. Only MXE uses all 8.

## bigWig naming

```
<sample_id>.<content>.<genome_build>.<strand>.bw
```

`content` is one of `unique`, `uniquemulti`; `strand` is one of `plus`, `minus`, `unstranded`.

STAR's `str1`/`str2` never reaches a filename: which one is the plus strand depends on the library chemistry, and the mapping is resolved once from the RSeQC call. Values are **positive on both strands**, unlike the pre-1.4 files, which stored minus-strand coverage as negative numbers. Tracks are CPM-normalized and the factor actually applied is recorded in `bigwig_manifest.tsv`; the manifest, not the filename, is authoritative.
