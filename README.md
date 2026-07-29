# RNA-Seq Pipeline Usage Guide 

## Table of Contents
- [General Information](#general-information)
- [Getting Started](#getting-started)
- [Running the Pipeline on Docker](#running-the-pipeline-on-docker)
- [Outputs](#outputs)
- [Validating a run](#validating-a-run)
- [Development](#development)

## General Information

### What this pipeline does
This RNA-seq pipeline processes raw sequencing data through differential analysis and reporting. Key analytical outputs include:

- Read quality assessment (`FastQC`, `fastp`).
- Genome alignment (`STAR`).
- Strandedness inference (`RSeQC`).
- Transcript abundance quantification (`kallisto`) and gene-level summarization.
- Differential expression analysis (`DESeq2`).
- Optional ERCC-based normalization (if ERCC spike-ins are present).
- Transcript-level differential analysis (`Sleuth`).
- Gene set enrichment analysis (`GSEA`, MSigDB gene set collections).
- Alternative splicing analysis (`rMATS`).
- Final report generation (`report.html`) and organized output directories.

Every one of those steps also writes a machine-readable table under `<results_dir>/tables/`,
described in [docs/SCHEMA.md](docs/SCHEMA.md). Those tables, not the spreadsheets or the HTML
report, are what gets loaded into the lab's data warehouse.

### Input structure 
- Input parent directory name must be in `Model_Experiment` format. The **first** underscore
  separates the cell model from the experiment descriptor; both parts must be non-empty.
- The descriptor is free-form after that, so `SNU719_Rta-Zta-2025-04` and `SNU719_Rta-Zta_2025-04`
  are both valid and both describe the `SNU719` model.
- **This name is permanent.** It becomes `experiment_id`, and `sample_id` and `comparison_id` are
  built from it, so it has to be unique across the whole archive. If an experiment is ever
  repeated, give the second one a distinguishing suffix such as a date — otherwise the two
  produce identical sample IDs and any query grouped by sample silently pools them.
- Parent directory must contain two condition subdirectories: `cntl` and `test`.
- FASTQ files are placed inside those two condition folders. The condition comes from the
  directory, not from the filename, and replicate numbers are assigned by sorting filenames
  within each condition (zero-pad past nine, or `cntl10` will sort before `cntl2`).

Example:

```text
SNU719_Zta-plus-Rta/
├── cntl
│   ├── SNU_Cntl1_1.fq.gz
│   ├── SNU_Cntl1_2.fq.gz
│   ├── SNU_Cntl2_1.fq.gz
│   ├── SNU_Cntl2_2.fq.gz
│   ├── SNU_Cntl3_1.fq.gz
│   └── SNU_Cntl3_2.fq.gz
└── test
    ├── SNU_Zta1_1.fq.gz
    ├── SNU_Zta1_2.fq.gz
    ├── SNU_Zta2_1.fq.gz
    ├── SNU_Zta2_2.fq.gz
    ├── SNU_Zta3_1.fq.gz
    └── SNU_Zta3_2.fq.gz
```

### Reference files and supported references
This repository provides a folder named `/referenceFiles`, which contains assets used by the pipeline (for example STAR index, annotation files, kallisto index, Biomart export, and GSEA gene sets). 

When using the pipeline, you must provide a reference directory path you pass in commands. In this README, examples from inside the Docker Container use:

`/data/referenceFiles`

Supported species currently include:
- `hg38`
- `hg38plusAkataInverted`
- `hg38plusKSHV`
- `hg38plusKSHVALT`
- `mm39`
- `mm39plusMHV68`

## Getting Started

### 1) Clone this repository: 
  ```bash
  cd /path/you/want/the/folder/to/exist
  git clone <PASTE_GITHUB_REPO_URL_HERE> rnaseq-pipeline
  ```

### 2) Ensure [Docker Desktop](https://www.docker.com/products/docker-desktop/) is installed and running. 

After Docker Desktop opens, verify it is running from Terminal:
```bash
docker --version
docker ps
```
If `docker ps` shows an error about connecting to the Docker daemon, Docker Desktop is not running yet.


## Running the Pipeline on Docker

### Prerequisites

#### 1) Docker Desktop File Sharing includes host directories you will mount:
  - Open Docker Desktop, then go to `Settings` -> `Resources` -> `File Sharing`.
  - Add or confirm the parent folders you will mount (for example `/Users`, `~/Documents`, etc).
  - Click `Apply & Restart` if Docker prompts you.

#### 2) Move into the rnaseq-pipeline folder (do this every new terminal session):
  ```bash
  cd /path/to/rnaseq-pipeline
  ```
- Confirm you are in the correct folder:
  ```bash
  pwd
  ls
  ```
  You should see `Dockerfile`, `metadata.py`, and `rnaseq.py` in the `ls` output.

#### 3) Build the Docker image
```bash
docker build --build-arg GIT_COMMIT="$(git rev-parse HEAD)" -t rnaseqpipeline:latest .
```

`GIT_COMMIT` is required. The build context excludes `.git`, so the commit cannot be read from
inside the image, and every result records the commit it came from in `run_manifest.json`. A
build without it fails rather than producing an image whose output cannot be traced back to a
version of the code.

### CPU, memory, and concurrent runs

To give the container as much compute headroom as your host allows, open Docker Desktop and go to **Settings** → **Resources**. Increase **CPUs** and **Memory** toward what your machine can spare (STAR alignment, abundance estimation, and downstream steps are demanding), then click **Apply & Restart** if Docker asks you to.

This pipeline is **computationally and memory intensive**. Do **not** run multiple instances of the pipeline—or other heavy workloads—at the same time on the same machine unless you know you have spare capacity. Competing processes slow runs sharply and can trigger out-of-memory or other downstream failures.

### Pipeline Arguments 
#### For metadata.py: 
- `<fastq_root_dir>`: Filepath of Fastq Directory (the one that contains `cntl_*` and `test_*` subdirectories) 
- `<reference_dir>`: Filepath of `referenceFiles`. 
- `<species_name>`: e.g. hg38, mm39, etc. It must match an existing directory under `referenceFiles`.
- `<investigator_name>`: Name of investigator, should not contain spaces.
- `<PE | SE>`: use `PE` for paired-end data or `SE` for single-end data.
- `<results_dir>`: Filepath of where you want the output to exist. The folder must be empty.

`metadata.py` also records what the experiment actually was, so results stay interpretable
years later. It prompts for each of these; pass them as flags to skip the prompts, and add
`--non-interactive` to make a missing value an error instead:

- `--cell-line`, `--organism`, `--perturbation-type`, `--induced-program`: controlled
  vocabularies. An unrecognized value is rejected with the list of allowed values. To add one,
  edit `rnaseq_helper_scripts/vocab.py`.
- `--perturbation-target` (e.g. `BMRF1`), `--perturbation-dose` (e.g. `100nM`),
  `--timepoint-hours`, `--sequencing-run-date`, `--notes`.

`--perturbation-type` is *how* the perturbation was delivered; `--induced-program` is *what it
was meant to induce*. They are separate axes because lytic reactivation can be driven by Zta or
Rta transfection, by TPA/butyrate, or by BCR crosslinking — keeping them apart is what makes
"every reactivation experiment, regardless of method" a single query. Unlike the
`--perturbation-*` fields, `--induced-program` describes the experiment and so is recorded on
control samples too. Use `unknown` for an unrecoverable value and `none` for a recorded absence;
they are not the same thing.

#### For rnaseq.py: 
- `<metadata_file>`: The file path of generated metadata tsv file `output_dir/*.tsv`.
- `<reference_dir>`: Same argument as metadata.py.
- `<scripts_dir>`: The `./rnaseq_helper_scripts` folder in this directory.
- `<results_dir>`: Same argument as metadata.py. 

Optional:
- `--annotation-version`: override the annotation release recorded in the manifest. Needed only
  when the GTF has no version in its header and no `ANNOTATION_VERSION` marker file.
- `--skip-transcript-expression`: write `expression_transcript.tsv` with headers only. It is by
  far the largest table.
- `--normalization`: bigWig normalization. Leave at the default; it has to be uniform across
  runs for tracks to be comparable.


### 1) Run interactively (recommended first run)
#### Mount references, input FASTQs, and output directory:

```bash
docker run --rm -it \
  -v "/path/on/your/computer/to/referenceFiles/:/data/referenceFiles:ro" \
  -v "/path/on/your/computer/to/Model_Experiment:/data/Model_Experiment:ro" \
  -v "/path/on/your/computer/for/output_folder:/data/output" \
  -w /work \
  rnaseqpipeline:latest \
  bash
```

#### **IMPORTANT NOTE**: How paths change inside the container
When you use `-v` to mount folders, Docker maps folders from your computer to new paths inside the container.    

Use the **container paths** (right side of every colon in each `-v` line) when calling `metadata.py` and `rnaseq.py` inside the container.

- `fastq_root_dir`: `/data/Model_Experiment`
- `reference_dir`: `/data/referenceFiles`
- `results_dir`: `/data/output`
- `scripts_dir`: `/work/rnaseq_helper_scripts`

Quick mapping examples from the command above:
- `/path/on/your/computer/to/Model_Experiment` -> `/data/Model_Experiment`
- `/path/on/your/computer/to/referenceFiles` -> `/data/referenceFiles`
- `/path/on/your/computer/for/output_folder` -> `/data/output`

#### Inside the container:

#### 1. Activate the pipeline environment:
   ```bash
   conda activate rnaseqpipeline
   ```  

#### 2. Generate metadata:
   ```bash
   metadata.py <fastq_root_dir> <reference_dir> <species_name> <investigator_name> <PE|SE> <results_dir>
   ``` 
   
   Example:
   ```bash
   python metadata.py /data/SNU719_Zta-plus-Rta /data/referenceFiles hg38plusAkataInverted ethan PE /data/output
   ```

#### 3. Run the pipeline
   ```bash
   python rnaseq.py <metadata_file> <reference_dir> <scripts_dir> <results_dir>
   ```
   Example:
   ```bash
   python rnaseq.py /data/output/ethan_metadata_04232026_174629.tsv /data/referenceFiles /work/rnaseq_helper_scripts /data/output
   ```
 

### 2) One-shot Docker command (non-interactive)
```bash
./run_pipeline_one_shot.sh
```

Edit variables at the top of `run_pipeline_one_shot.sh` before running:
- `HOST_REFERENCE_DIR`
- `HOST_FASTQ_ROOT_DIR`
- `HOST_OUTPUT_DIR`
- `SPECIES_NAME`
- `INVESTIGATOR_NAME`
- `EXPERIMENT_TYPE` (`PE` or `SE`)
- `CELL_LINE`, `PERTURBATION_TYPE`, `PERTURBATION_TARGET`, `PERTURBATION_DOSE`,
  `TIMEPOINT_HOURS`, `SEQUENCING_RUN_DATE` — the script runs `metadata.py --non-interactive`,
  so all of these must be filled in.

### Checklist
- Parent folder naming follows `Model_Experiment` (single underscore).
- Subdirectories are exactly `cntl` and `test`.
- Reference name exists in `referenceFiles`.
- Environment is `rnaseqpipeline` (native or container).
- Metadata and rnaseq are run with matching paths in the selected environment.

## Outputs

Full column-by-column reference: **[docs/SCHEMA.md](docs/SCHEMA.md)**. What the pipeline wrote
before this layout existed is recorded in [docs/CURRENT_OUTPUTS.md](docs/CURRENT_OUTPUTS.md).

```text
<results_dir>/
├── run_manifest.json   provenance: commit, image, tool versions, reference, exit status
├── metadata.tsv        one row per sample, with the inferred strandedness filled in
├── checksums.sha256
├── tables/             the machine-readable results; this is what the warehouse reads
├── artifacts/          raw tool output, bigWigs and BAMs, kept for re-analysis
├── reports/            XLSX, figures and report.html, for reading rather than querying
└── logs/               one log per stage
```

Everything in `tables/` is a TSV in long format: one row per gene, transcript, event or metric,
never one column per sample. Missing values are the literal string `NA`. Every row carries a
`run_id` that joins back to `run_manifest.json`.

### Differences from pre-1.4 output

- `deseq2_results_genes.tsv` is split. Statistics live in `tables/de_gene.tsv`; per-sample TPMs
  live in `tables/expression_gene.tsv`. Neither has sample names in its header. The XLSX is
  still produced for the bench.
- bigWigs are named `<sample_id>.<content>.<build>.<strand>.bw`. STAR's `str1`/`str2` no longer
  appears anywhere: which one is the plus strand depends on the library chemistry, so it is
  resolved from the RSeQC call and written into the filename and the manifest.
- bigWig values are CPM-normalized and **positive on both strands**. Minus-strand tracks used to
  hold negative values. Genome browser sessions that relied on the old sign will need their
  track ranges adjusted.
- Output filenames are no longer rewritten after the fact. The four `rename_*` passes are gone,
  along with the naming inconsistencies they caused.

## Validating a run

`rnaseq.py` runs the validator itself and records the result in the manifest. To re-check a
directory later, or to check an older run:

```bash
python rnaseq_helper_scripts/validate_outputs.py /data/output
```

It fails on a missing or incomplete manifest, a table whose columns do not match the schema, a
`sample_id` that is not in `metadata.tsv`, a gene ID that is neither an unversioned Ensembl
accession nor a listed exception, a splicing replicate with no parent event, and any checksum
that no longer matches. Legitimate non-Ensembl IDs — ERCC spike-ins, EBV and KSHV loci — belong
in `docs/gene_id_exceptions.txt`.

## Development

```bash
pip install pandas openpyxl pytest
python -m pytest
```

The tests synthesize tool output and run the real writers over it, so no bioinformatics tools
are needed. `docs/SCHEMA.md` is generated from `rnaseq_helper_scripts/schemas.py`; after
changing a schema, regenerate it, or the test suite will fail:

```bash
python rnaseq_helper_scripts/generate_schema_doc.py
```

