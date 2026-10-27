# RNA-Seq Pipeline

## Contents
- [What it does](#what-it-does)
- [Input](#input)
- [Setup](#setup)
- [Running the pipeline](#running-the-pipeline)
- [Arguments](#arguments)
- [Outputs](#outputs)
- [Validating a run](#validating-a-run)
- [Development](#development)

## What it does

Takes the FASTQs of a control-vs-test experiment through to differential analysis and a report:

- Read QC (`FastQC`, `fastp`)
- Alignment (`STAR`) and strandedness inference (`RSeQC`)
- Transcript quantification (`kallisto`) and gene-level summaries
- Differential expression (`DESeq2`, with ERCC normalization when spike-ins are present) and
  transcript-level analysis (`Sleuth`)
- Gene set enrichment (`GSEA`, every MSigDB collection in the reference folder)
- Alternative splicing (`rMATS`)
- Raw-count bigWigs, splice-junction BEDs, and signal over each gene
- Figures, spreadsheets and `report.html`

Every step also writes a machine-readable table to `<results_dir>/tables/`, described in
[docs/SCHEMA.md](docs/SCHEMA.md). Those tables, not the spreadsheets or the report, are what the
lab's database loads.

## Input

```text
SNU719_Rta-Zta_2025-04-10/
├── input_args.json      from the form (see below)
├── cntl/
│   ├── SNU719_Rta-Zta_cntl1_1.fq.gz
│   ├── SNU719_Rta-Zta_cntl1_2.fq.gz
│   └── ...
└── test/
    ├── SNU719_Rta-Zta_test1_1.fq.gz
    ├── SNU719_Rta-Zta_test1_2.fq.gz
    └── ...
```

- **Folder name is `Model_Experiment`.** The first underscore separates the cell model from the
  descriptor; after that the descriptor is free-form (`SNU719_Rta-Zta_2025-04-10` is model
  `SNU719`).
- **The folder name is permanent.** It becomes `experiment_id`, and every `sample_id` and
  `comparison_id` is built from it, so it must be unique across the archive. If an experiment is
  repeated, add a distinguishing suffix such as a date, or the two runs get identical sample IDs.
- **Exactly two subfolders, `cntl` and `test`.** The condition comes from the folder, not the
  filename. Replicates are numbered by sorting filenames within each folder, so zero-pad past
  nine (`cntl10` sorts before `cntl2`).

### Filling in `input_args.json` with the form

`metadata_form.html` is a self-contained page (double-click to open; no server or install) that
asks for each experiment detail with a dropdown and downloads the answers as `input_args.json`.

1. Fill in every field. The download button stays disabled, listing what is missing, until all
   are done.
2. Put the downloaded file in the experiment folder next to `cntl/` and `test/`, **keeping the
   name `input_args.json`**. `metadata.py` ignores loose files when it looks for samples, so it
   does not interfere.
3. Pass `--quick-input` to `metadata.py` (the one-shot script does this for you).

- The file is read only with `--quick-input`.
- An explicit flag overrides the file, e.g. `--quick-input --library-selection ribodepleted`.
- Values from the file are not re-checked against `vocab.py`; the dropdowns are the check. Don't
  hand-edit the JSON; re-download it from the form.
- The dropdowns are a hand-maintained copy of `rnaseq_helper_scripts/vocab.py`. When you add a
  value to `vocab.py`, add the matching `<option>` to `metadata_form.html` in the same commit.

### Reference files

`referenceFiles` is not in this repository. On the Mac15 it is at
`/Applications/ngs/pipelines/docker-rnaseq/referenceFiles` and holds, per genome build, the STAR
index, annotations, kallisto index, Biomart export and GSEA gene sets. Supported builds:
`hg38`, `hg38plusAkataInverted`, `hg38plusKSHV`, `hg38plusKSHVALT`, `mm39`, `mm39plusMHV68`.

## Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/eisrael123/rnaseq-pipeline.git
   cd rnaseq-pipeline
   ```
2. Install and start [Docker Desktop](https://www.docker.com/products/docker-desktop/). `docker ps`
   should run without a "cannot connect to the Docker daemon" error.
3. In Docker Desktop, **Settings → Resources**:
   - **File Sharing:** add every parent folder you will mount, including `/Volumes` for external
     drives and the lab share.
   - **CPUs and Memory:** give Docker as much as the machine can spare. STAR may use up to 60 GB to
     sort BAMs, and GSEA runs with a 16 GB Java heap.
4. Build the image from the repository folder:
   ```bash
   docker build --build-arg GIT_COMMIT="$(git rev-parse HEAD)" -t rnaseqpipeline:latest .
   ```
   `GIT_COMMIT` is required. Every run records the commit it came from in `run_manifest.json`,
   and the image cannot read it itself (`.git` is not copied in), so a build without it fails.
   Rebuild after every `git pull`.

Run one pipeline at a time. It is CPU- and memory-heavy, and competing jobs slow it sharply and
can cause out-of-memory failures.

## Running the pipeline

### One-shot script (normal use)

Requires `input_args.json` in the experiment folder. The script checks for it and stops before
starting Docker if it is missing.

1. Edit the two variables at the top of `run_pipeline_one_shot.sh`:
   - `HOST_FASTQ_ROOT_DIR`: the experiment folder
   - `HOST_OUTPUT_DIR`: where the output goes; must be empty
2. Run it:
   ```bash
   ./run_pipeline_one_shot.sh
   ```

`HOST_REFERENCE_DIR` is already set for the Mac15; change it only on another machine. To override
one experiment detail for a single run, add the flag to the `metadata.py` call in the script.

### Step by step in an interactive container

Useful for debugging, or for running without the form.

```bash
docker run --rm -it \
  -v "/path/to/referenceFiles:/data/referenceFiles:ro" \
  -v "/path/to/SNU719_Rta-Zta_2025-04-10:/data/SNU719_Rta-Zta_2025-04-10:ro" \
  -v "/path/to/output_folder:/data/output" \
  -w /work \
  rnaseqpipeline:latest bash
```

**Mount the FASTQ folder under its real name.** `metadata.py` takes `experiment_id` from the
folder name it is given, so a generic mount point such as `/data/Model_Experiment` would give
every experiment the same ID. `metadata.py` rejects that name for this reason.

Inside the container, use the container paths (the right-hand side of each `-v`):

```bash
conda activate rnaseqpipeline

# With input_args.json in the experiment folder:
python metadata.py \
  --root-fastq-dir /data/SNU719_Rta-Zta_2025-04-10 \
  --output-dir /data/output \
  --quick-input

# Or without it: give the run arguments and answer the prompts (or pass them as flags):
python metadata.py \
  --root-fastq-dir /data/SNU719_Rta-Zta_2025-04-10 \
  --output-dir /data/output \
  --reference-dir /data/referenceFiles \
  --species-name hg38plusAkataInverted \
  --investigator-name ethan \
  --experiment-type PE

python rnaseq.py /data/referenceFiles /work/rnaseq_helper_scripts /data/output
```

## Arguments

### metadata.py

All arguments are named flags. Two are always required, because only you know where the data is:

- `--root-fastq-dir`: the experiment folder (contains `cntl/` and `test/`)
- `--output-dir`: where the output goes; must be empty

The rest come from `input_args.json` (`--quick-input`), from flags, or from interactive prompts.
Add `--non-interactive` to make a missing value an error instead of a prompt.

Run settings:

- `--reference-dir`: path to `referenceFiles`
- `--species-name`: genome build; must be a folder in `referenceFiles`
- `--investigator-name`: whitespace is stripped, since it becomes part of output filenames
- `--experiment-type`: `PE` (paired-end) or `SE` (single-end)

Experiment description, kept so results stay interpretable years later. Controlled-vocabulary
fields reject unknown values and list the allowed ones; to add a value, edit
`rnaseq_helper_scripts/vocab.py` (and the form).

| Flag | Vocabulary | Meaning |
|---|---|---|
| `--cell-line`, `--organism` | controlled | |
| `--perturbation-type` | controlled | *How* it was delivered: `transfection`, `chemical`, `BCR-crosslink`, `none` |
| `--perturbation-agent` | free text | *What* was delivered: `Zta`, `Zta+Rta`, `anti-IgG`, `CC115` |
| `--perturbation-target` | free text | The gene the agent acts on (`BMRF1` for an siRNA), or `NA` when it has none, as for Zta or Rta overexpression |
| `--perturbation-dose` | free text | e.g. `5ug+5ug`, `100nM` |
| `--co-treatment`, `--co-treatment-target` | controlled | A second treatment and what it acts on. Zta under PAA is agent `Zta` with co-treatment `PAA`, target `PAA_replication`, so "every PAA experiment" stays one filter |
| `--facs-purified`, `--facs-gfp-promoter` | controlled | Whether, and on what, cells were sorted. Recorded on control rows too. `pCMV`-GFP selects transfected cells; `BMRF1p`-GFP selects cells that entered the lytic cycle |
| `--library-selection` | controlled | `polyA`, `ribodepleted`, or `unknown` for backfilled runs. No default on purpose: the two methods are not directly comparable, so it must be recorded, not assumed |
| `--library-strandedness` | controlled | What was ordered at library prep. RSeQC's measured value is stored separately in `metadata.tsv`; a disagreement points to a mislabeled sample or the wrong kit |
| `--timepoint-hours`, `--sequencing-run-date`, `--notes` | free text | |

### rnaseq.py

```bash
python rnaseq.py <reference_dir> <scripts_dir> <results_dir>
```

- `<reference_dir>`: same as `--reference-dir`
- `<scripts_dir>`: `rnaseq_helper_scripts` (`/work/rnaseq_helper_scripts` in the container)
- `<results_dir>`: same as `--output-dir`. `metadata.py` must have written `metadata.tsv` there
  first; `rnaseq.py` updates it in place (e.g. with the RSeQC strandedness).

Optional:

- `--annotation-version`: annotation release to record in the manifest, needed only when the GTF
  header and `ANNOTATION_VERSION` marker file give none
- `--skip-transcript-expression`: write `expression_transcript.tsv` with headers only (it is by
  far the largest table)
- `--normalization`: normalization of `signal_over_gene`. Leave the default, so values stay
  comparable across runs

## Outputs

Column-by-column reference: [docs/SCHEMA.md](docs/SCHEMA.md). What the pipeline wrote before
this layout existed is recorded in [docs/CURRENT_OUTPUTS.md](docs/CURRENT_OUTPUTS.md).

```text
<results_dir>/
├── run_manifest.json   provenance: commit, image, tool versions, reference, status
├── metadata.tsv        one row per sample, with measured strandedness filled in
├── checksums.sha256
├── tables/             machine-readable results; what the database reads
├── artifacts/          raw tool output, BAMs and bigWigs, for re-analysis
│   ├── bigwig/unique/        raw-count bigWigs, uniquely mapped reads
│   ├── bigwig/uniquemulti/   raw-count bigWigs, unique + multi-mapped reads
│   ├── rmats_raw/JCEC/       the five *.MATS.JCEC.txt files and summary.txt
│   ├── deseq2/, sleuth/      XLSX at the top, intermediates in prep_files/
│   ├── kallisto/<sample>/
│   └── gsea/                 one folder per gene set collection
├── reports/            XLSX, figures and report.html, for reading
└── logs/               one log per stage; errors.log if any output failed
```

- **Tables** are long-format TSVs: one row per gene, transcript, event or metric, never one column
  per sample. Missing values are `NA`. Every row has a `run_id` that joins to `run_manifest.json`.
- **bigWigs** hold STAR's raw counts (whole reads for unique; 1/N fractions for multi-mappers),
  named `<sample_id>.<content>.<build>.<strand>.bw`. The strand is resolved from RSeQC, not
  STAR's `str1`/`str2`. Minus-strand tracks are **negative** so they mirror below the axis in a
  genome browser; use absolute values to compare coverage. Unstranded libraries get one
  `unstranded` track per content type. No normalized bigWigs are written.
- **`signal_over_gene`** is in CPM (raw signal × 1e6 / mapped reads), with positive values on
  every strand, so it is comparable between samples.
- **Figures** (volcano, PCA, splicing bar charts, ERCC) are publication-ready: Arial, large text,
  1 pt black axes and outlines, 300 dpi, colours from the lab's 2026 bioRxiv preprint (Nguyen et
  al.). See `rnaseq_helper_scripts/figure_style.py`.
- **File names** people open carry the experiment name: `<experiment>_test_vs_cntl_` on figures
  and JCEC files, `<experiment>_deseq2_` / `<experiment>_sleuth_` on DESeq2 and Sleuth output
  (except `deseq2_metadata.tsv`), and `<experiment>_de_gene.xlsx`. Summary tables and
  `report.html` are unprefixed.
- **Strandedness:** `report.html` shows `rf-stranded (fr-firststrand)`,
  `fr-stranded (fr-secondstrand)` or `unstranded`. If RSeQC cannot decide, the sample is analysed
  as unstranded and the run says so on the terminal, in the log and RSeQC file, as
  `strandedness_fallback` in the manifest, and in `report.html`; `metadata.tsv` keeps
  `undetermined`.
- **GSEA** runs every collection in the reference folder, including the full `gene_ontology` GMT
  (the only one with cellular-component sets), with a 16 GB Java heap set through
  `JVM_MEM_OPTS`.
- **kallisto** in Docker uses an index rebuilt in April 2026 with kallisto 0.52.0; on macOS it
  uses the original `*_ORIGINAL.kallisto.idx`. Abundances can differ slightly, so the manifest
  records the index path and SHA-256 under `parameters.kallisto_index`.
- **Failures:** if a spreadsheet, figure, junction BED, `report.html` or GSEA collection fails,
  the run finishes but reports it as an error on the terminal, in `logs/errors.log` and in the
  manifest's `output_errors`. Validation then fails and the run exits non-zero.

### Older runs in the archive

Runs made between the refactor and October 2026 (all 70 archived experiments) differ from the
above:

- CPM bigWigs in `artifacts/bigwig/`, with raw-count copies in a separate
  `raw_count_bigwigs/<experiment>/` folder on the share. The raw copies were reconstructed from
  the CPM tracks (CPM ÷ `scale_factor`, rounded), not taken from STAR.
- Unprefixed file names, and kallisto output one level deeper, in
  `artifacts/kallisto/kallisto/<sample>/`.
- No `gene_ontology` GSEA (skipped from commit 9ea83ea, July 2026). GSEA ran with the conda
  `gsea-cli`'s 4 GB heap, too little for that file. (The original pipeline's `gsea-cli.sh` had
  `-Xmx12g` written into it.)
- Older figure style.

Compared with the original pipeline's output, `deseq2_results_genes.tsv` is split into
`tables/de_gene.tsv` (statistics) and `tables/expression_gene.tsv` (per-sample TPMs), and
`results.html` is now `reports/report.html`.

## Validating a run

`rnaseq.py` validates its own output and records the result in the manifest. To re-check a run
later, or check an older one:

```bash
python rnaseq_helper_scripts/validate_outputs.py /data/output
```

It fails on:
- a missing or incomplete manifest
- a table whose columns do not match the schema
- a `sample_id` not in `metadata.tsv`
- a gene ID that is neither an unversioned Ensembl ID nor listed in
  `docs/gene_id_exceptions.txt` (where ERCC spike-ins and EBV/KSHV loci belong)
- a splicing replicate with no parent event
- a checksum mismatch
- an output the run recorded as failed (`output_errors`)
- a missing DESeq2 XLSX, volcano plot or `report.html`

## Development

```bash
pip install pandas openpyxl pytest
python -m pytest
```

The tests synthesize tool output and run the real writers over it, so no bioinformatics tools are
needed. `docs/SCHEMA.md` is generated from `rnaseq_helper_scripts/schemas.py`; after changing a
schema, regenerate it, or the tests fail:

```bash
python rnaseq_helper_scripts/generate_schema_doc.py
```
