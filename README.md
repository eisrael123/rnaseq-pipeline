# RNA-Seq Pipeline Usage Guide 

## Table of Contents
- [General Information](#general-information)
- [Getting Started](#getting-started)
- [Running the Pipeline on Docker](#running-the-pipeline-on-docker)
  - [Filling in the parameters with the form](#filling-in-the-parameters-with-the-form)
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

Every argument is a named flag; there are no positional arguments. There are two ways to supply
them, and they can be mixed:

| | What you type | Where the rest comes from |
|---|---|---|
| **With the form** (normal) | Two directories, plus `--quick-input` | `input_args.json` in the FASTQ folder |
| **By hand** | Every flag, or answer the prompts | You |

Two flags are always required, because they name where the data sits on *this* machine and
nothing else can know that:

- `--root-fastq-dir`: Filepath of the FASTQ directory (the one that contains the `cntl` and
  `test` subdirectories).
- `--output-dir`: Filepath of where you want the output to exist. The folder must be empty.

The rest can come from a flag, from an interactive prompt, or from
[the form](#filling-in-the-parameters-with-the-form) via `--quick-input`:

- `--reference-dir`: Filepath of `referenceFiles`.
- `--species-name`: e.g. hg38, mm39, etc. It must match an existing directory under
  `referenceFiles`.
- `--investigator-name`: Name of investigator. Whitespace is stripped, since the name becomes
  part of the output filename.
- `--experiment-type`: `PE` for paired-end data or `SE` for single-end data.

`metadata.py` also records what the experiment actually was, so results stay interpretable
years later. The form asks for all of these; without it, `metadata.py` prompts for each. Pass
them as flags to skip the prompts, and add `--non-interactive` to make a missing value an error
instead:

- `--cell-line`, `--organism`, `--perturbation-type`, `--co-treatment`,
  `--co-treatment-target`, `--facs-purified`, `--facs-gfp-promoter`, `--library-selection`,
  `--library-strandedness`: controlled vocabularies. An unrecognized value is rejected with the
  list of allowed values. To add one, edit `rnaseq_helper_scripts/vocab.py`.
- `--perturbation-agent` (e.g. `Zta`, `Zta+Rta`, `anti-IgG`), `--perturbation-target` (e.g.
  `BMRF1`, or `NA`), `--perturbation-dose` (e.g. `5ug+5ug`, `100nM`), `--timepoint-hours`,
  `--sequencing-run-date`, `--notes`: free text.

The perturbation is recorded as three separate axes, because collapsing them is what makes an
archive unqueryable. `--perturbation-type` is *how* it was delivered (`transfection`,
`chemical`, `BCR-crosslink`, `none`); `--perturbation-agent` is *what* was delivered (`Zta`,
`anti-IgG`, `CC115`); `--perturbation-target` is the gene the agent acts on, which is `NA` when
it has none — Zta or Rta overexpression has no separate molecular target, but an siRNA does.

`--co-treatment` is a second treatment applied alongside the first, with
`--co-treatment-target` naming what it acts on. Zta under PAA is a Zta perturbation with a `PAA`
co-treatment, not a single combined value, so "every PAA experiment" stays one filter. Use
`PAA_replication` as the target where the co-treatment blocks viral DNA replication rather than
acting on a host gene.

`--facs-purified` and `--facs-gfp-promoter` describe the sequenced material, so they are
recorded on control rows too. An unsorted transfection is a mixture of transfected and
untransfected cells, which changes what an expression value means; and sorting on `pCMV`-driven
GFP selects cells that were transfected, while `BMRF1p` selects cells where the lytic cycle
actually started. Those are different populations.

`--library-strandedness` is what was ordered at library prep. It is deliberately *not* the same
field as the `strandedness` column in `metadata.tsv`, which RSeQC measures per sample after
alignment. The two normally agree, and a disagreement is a useful flag for a mislabeled sample
or the wrong kit, which is why both are kept.

`--library-selection` is `polyA` or `ribodepleted` (or `unknown` when backfilling older runs),
and it has no default on purpose. It is not a comparable axis: non-polyadenylated and
unprocessed transcripts are absent from a polyA library by construction, so the same gene can
read as absent in one library and abundant in another for reasons that have nothing to do with
the biology. Recording it per sample is what lets a cross-experiment query either filter to one
selection method or state that it is mixing them; defaulting it would make the archive quietly
claim otherwise. It sits alongside `library_layout` (`PE`/`SE`) and the RSeQC-inferred
`strandedness` in `metadata.tsv`.

### Filling in the parameters with the form

Typing a dozen flags correctly is not the job of whoever ran the experiment. `metadata_form.html`
is a single self-contained page — open it by double-clicking, no server and no install — that
asks for each value with a dropdown, then downloads the answers as `input_args.json`.

1. Open `metadata_form.html` in a browser and fill in every field. The download button stays
   disabled until nothing is missing, and it lists what is still outstanding.
2. Drag the downloaded `input_args.json` into the experiment's FASTQ folder, alongside the
   `cntl` and `test` subdirectories. **Keep the filename exactly as downloaded** — that fixed
   name is how the pipeline finds it.
   ```text
    SNU719_Zta-plus-Rta/
    ├── input_args.json
    ├── cntl 
      └── ...
    └── test
      └── ...
    ```
3. Then, pass `--quick-input`, and the two directory arguments:

```bash
python metadata.py \
  --root-fastq-dir /data/SNU719_Zta-plus-Rta \
  --output-dir /data/output \
  --quick-input
```

The file sits next to the condition folders rather than inside them, so it does not disturb
sample discovery — `metadata.py` walks directories and ignores loose files.

Notes on how it behaves:

- **`--quick-input` is opt-in.** Without the flag the file is never read, even if it is sitting
  right there. Nothing about the existing flag-driven or interactive workflows changes.
- **An explicit flag still wins.** `--quick-input --library-selection ribodepleted` overrides
  what the file says, so a one-off rerun does not need the file edited and re-downloaded.
- **Values that came from the form are not re-checked against `vocab.py`.** The dropdown is
  what constrains them, and a second copy of the allowed values in `metadata.py` could only
  drift from the first. The consequence is worth knowing: a **hand-edited** `input_args.json`
  can put a value into the archive that the vocabulary would have rejected. Re-download from
  the form rather than editing the JSON.
- **The form's dropdowns are a hand-maintained copy of `rnaseq_helper_scripts/vocab.py`.** They
  agree today. If you add a cell line, genome build, or any other vocabulary value to
  `vocab.py`, add the matching `<option>` to `metadata_form.html` in the same commit, or the
  form will silently be unable to offer a value the pipeline supports.

#### For rnaseq.py: 
- `<reference_dir>`: Same argument as metadata.py.
- `<scripts_dir>`: The `./rnaseq_helper_scripts` folder in this directory.
- `<results_dir>`: Same argument as metadata.py. rnaseq.py reads `<results_dir>/metadata.tsv`
  directly -- metadata.py must have already generated it there -- and updates it in place (e.g.
  once RSeQC infers `rseqc_measured_strandedness`).

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
  -v "/path/on/your/computer/to/SNU719_Rta-Zta-2025-04-10:/data/SNU719_Rta-Zta-2025-04-10:ro" \
  -v "/path/on/your/computer/for/output_folder:/data/output" \
  -w /work \
  rnaseqpipeline:latest \
  bash
```

**The FASTQ mount point must carry the real experiment directory name**, not a generic one.
`metadata.py` reads `experiment_id` from the basename of the path you hand it, and that id is
baked into every `sample_id`. Mounting to a fixed `/data/Model_Experiment` would label every
experiment in the archive identically; `metadata.py` rejects that name outright to make the
mistake loud rather than silent.

#### **IMPORTANT NOTE**: How paths change inside the container
When you use `-v` to mount folders, Docker maps folders from your computer to new paths inside the container.    

Use the **container paths** (right side of every colon in each `-v` line) when calling `metadata.py` and `rnaseq.py` inside the container.

- `--root-fastq-dir`: `/data/<your experiment directory name>`
- `--reference-dir`: `/data/referenceFiles`
- `--output-dir`: `/data/output`
- `scripts_dir` (rnaseq.py): `/work/rnaseq_helper_scripts`

Quick mapping examples from the command above:
- `/path/on/your/computer/to/SNU719_Rta-Zta-2025-04-10` -> `/data/SNU719_Rta-Zta-2025-04-10`
- `/path/on/your/computer/to/referenceFiles` -> `/data/referenceFiles`
- `/path/on/your/computer/for/output_folder` -> `/data/output`

#### Inside the container:

#### 1. Activate the pipeline environment:
   ```bash
   conda activate rnaseqpipeline
   ```  

#### 2. Generate metadata:

   If the experiment folder has an `input_args.json` from
   [the form](#filling-in-the-parameters-with-the-form):
   ```bash
   python metadata.py \
     --root-fastq-dir /data/SNU719_Zta-plus-Rta \
     --output-dir /data/output \
     --quick-input
   ```

   Otherwise pass the run arguments yourself, and answer the experiment prompts (or pass those
   as flags too — see [Pipeline Arguments](#for-metadatapy)):
   ```bash
   python metadata.py \
     --root-fastq-dir /data/SNU719_Zta-plus-Rta \
     --output-dir /data/output \
     --reference-dir /data/referenceFiles \
     --species-name hg38plusAkataInverted \
     --investigator-name ethan \
     --experiment-type PE
   ```

#### 3. Run the pipeline
   ```bash
   python rnaseq.py <reference_dir> <scripts_dir> <results_dir>
   ```
   Example:
   ```bash
   python rnaseq.py /data/referenceFiles /work/rnaseq_helper_scripts /data/output
   ```
 

### 2) One-shot Docker command (non-interactive)

This is the normal way to run the pipeline. It builds the metadata and runs the analysis in a
single container, start to finish.

**Prerequisite:** the experiment's FASTQ folder must already contain an `input_args.json` from
[the form](#filling-in-the-parameters-with-the-form). The script reads every experiment detail
from it and exits immediately, before starting Docker, if it is not there.

Edit exactly two variables at the top of `run_pipeline_one_shot.sh`:

- `HOST_FASTQ_ROOT_DIR` — the experiment folder (holds `cntl/`, `test/`, and `input_args.json`)
- `HOST_OUTPUT_DIR` — where the output goes; must be empty

Then run it:

```bash
./run_pipeline_one_shot.sh
```

Nothing else in the script needs touching. `HOST_REFERENCE_DIR` is already set to the Mac15
location of `referenceFiles` and only changes on a machine where it lives somewhere else; the
container paths below it are wired to the mounts.

The species, investigator, layout, cell line, perturbation and library fields that used to be
variables in this script are gone — they live in `input_args.json` now, and the script passes
`--quick-input`. To override one for a single run without re-downloading the file, add the flag
to the `metadata.py` call inside the script; an explicit flag beats the file.

### Checklist
- Parent folder naming follows `Model_Experiment` (single underscore).
- Subdirectories are exactly `cntl` and `test`.
- If using the form: `input_args.json` sits at the top level of the FASTQ folder, next to
  `cntl`/`test`, under exactly that name.
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

