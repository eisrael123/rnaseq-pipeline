# RNA-Seq Pipeline Usage Guide 

## Table of Contents
- [General Information](#general-information)
- [Getting Started](#getting-started)
- [Running the Pipeline on Docker](#running-the-pipeline-on-docker)
- [Adding a New Species (Reference Files)](#adding-a-new-species)

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

### Input structure 
- Input parent directory name must be exactly in `Model_Experiment` format.
- Use exactly one underscore (`_`) in the parent directory name.
- You may use dashes (`-`) inside the model or experiment names.
- Parent directory must contain two condition subdirectories: `cntl` and `test`.
- FASTQ files are placed inside those two condition folders.

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

If you want to run a new species/reference not listed above, create a new reference directory under `referenceFiles` with the same file/folder pattern used by existing references, and use that directory name as `<species_name>` in `metadata.py`. **[Click here](#adding-a-new-species)** for full instructions on how to do so.

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
docker build -t rnaseqpipeline:latest .
```

### Pipeline Arguments 
#### For metadata.py: 
- `<fastq_root_dir>`: Filepath of Fastq Directory (the one that contains `cntl_*` and `test_*` subdirectories) 
- `<reference_dir>`: Filepath of `referenceFiles`. 
- `<species_name>`: e.g. hg38, mm39, etc. It must match an existing directory under `referenceFiles`.
- `<investigator_name>`: Name of investigator, should not contain spaces.
- `<PE | SE>`: use `PE` for paired-end data or `SE` for single-end data.
- `<results_dir>`: Filepath of where you want the output to exist. The folder must be empty.

#### For rnaseq.py: 
- `<metadata_file>`: The file path of generated metadata tsv file `output_dir/*.tsv`.
- `<reference_dir>`: Same argument as metadata.py.
- `<scripts_dir>`: The `./rnaseq_helper_scripts` folder in this directory.
- `<results_dir>`: Same argument as metadata.py. 


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

### Checklist
- Parent folder naming follows `Model_Experiment` (single underscore).
- Subdirectories are exactly `cntl` and `test`.
- Reference name exists in `referenceFiles`.
- Environment is `rnaseqpipeline` (native or container).
- Metadata and rnaseq are run with matching paths in the selected environment.

## Adding a New Species 

### Required files
`rnaseq.py` and all helper scripts in `/rnaseq_helper_scripts` resolve files from `referenceFiles/<SPECIES>/...` using the species directory name itself.  

Make sure these files exist and follow this exact naming:

```text
referenceFiles/<SPECIES>/
├── annotations/
│   ├── <SPECIES>.gtf
│   ├── <SPECIES>.bed
│   └── ERCC92_genes.txt
├── fasta/
│   ├── <SPECIES>.fa
│   ├── <SPECIES>.fa.fai
│   ├── <SPECIES>.transcriptome.fa
│   └── <SPECIES>.transcriptome.fa.fai
├── STAR/
│   ├── Genome
│   ├── SA
│   ├── SAindex
│   ├── chrLength.txt
│   ├── chrName.txt
│   ├── chrNameLength.txt
│   ├── chrStart.txt
│   ├── exonGeTrInfo.tab
│   ├── exonInfo.tab
│   ├── geneInfo.tab
│   ├── genomeParameters.txt
│   ├── sjdbInfo.txt
│   ├── sjdbList.fromGTF.out.tab
│   ├── sjdbList.out.tab
│   └── transcriptInfo.tab
├── kallisto/
│   ├── <SPECIES>.kallisto.idx
│   └── all_contaminants.kallisto.idx
├── biomart/
│   └── <SPECIES>.mart_export.txt
└── GSEA/
    └── *.gmt
```

### What your current annotations folders look like
To inspect the exact file naming pattern already in your references:

```bash
SPECIES=hg38
ls -1 "/Users/mac14/ethan/rnaseq/referenceFiles/${SPECIES}/annotations"
```

Replace `hg38` with any existing species directory (`mm39`, `hg38plusKSHV`, etc.) and mirror that naming style when creating a new one.

### Why generate references in Docker
Generate references inside the same Docker image/environment used for the pipeline so:
- tool versions (`STAR`, `kallisto`, `samtools`, `gffread`) match runtime behavior,
- indexes are created in a Linux container filesystem (same as pipeline execution),
- path assumptions stay consistent between generation and analysis.

---

### Interactive Docker workflow (step-by-step generation)
The command below opens an interactive shell with:
- host `referenceFiles` mounted read/write at `/data/referenceFiles`
- host source downloads mounted read-only at `/data/src`

```bash
docker run --rm -it \
  -v "/Users/mac14/ethan/rnaseq/referenceFiles:/data/referenceFiles" \
  -v "/path/to/raw_reference_downloads:/data/src:ro" \
  -w /work \
  rnaseqpipeline:latest \
  bash
```

Inside the container:

```bash
conda activate rnaseqpipeline
SPECIES="myNewSpecies"
REFROOT="/data/referenceFiles/${SPECIES}"
mkdir -p "${REFROOT}"/{annotations,biomart,fasta,GSEA,kallisto,STAR}
```

#### 1) Build merged genome FASTA and index
```bash
# Example: concatenate host + optional spike-in/virus FASTAs
cat /data/src/genome.primary.fa /data/src/spikeins_or_virus.fa > "${REFROOT}/fasta/${SPECIES}.fa"
samtools faidx "${REFROOT}/fasta/${SPECIES}.fa"
```

#### 2) Build merged annotation GTF and species BED
```bash
# Example: concatenate host + optional spike-in/virus GTFs
cat /data/src/genome.primary.gtf /data/src/spikeins_or_virus.gtf > "${REFROOT}/annotations/${SPECIES}.gtf"

# Convert GTF -> BED (requires gtf2bed from BEDOPS in PATH)
gtf2bed < "${REFROOT}/annotations/${SPECIES}.gtf" > "${REFROOT}/annotations/${SPECIES}.bed"
```

#### 3) Build transcriptome FASTA for kallisto
```bash
gffread \
  -g "${REFROOT}/fasta/${SPECIES}.fa" \
  "${REFROOT}/annotations/${SPECIES}.gtf" \
  -w "${REFROOT}/fasta/${SPECIES}.transcriptome.fa"
samtools faidx "${REFROOT}/fasta/${SPECIES}.transcriptome.fa"
```

#### 4) Build STAR index
```bash
STAR \
  --runThreadN 16 \
  --runMode genomeGenerate \
  --genomeDir "${REFROOT}/STAR" \
  --genomeFastaFiles "${REFROOT}/fasta/${SPECIES}.fa" \
  --sjdbGTFfile "${REFROOT}/annotations/${SPECIES}.gtf" \
  --sjdbOverhang 149
```

#### 5) Build kallisto indexes
```bash
kallisto index \
  -i "${REFROOT}/kallisto/${SPECIES}.kallisto.idx" \
  "${REFROOT}/fasta/${SPECIES}.transcriptome.fa"

kallisto index \
  -i "${REFROOT}/kallisto/all_contaminants.kallisto.idx" \
  "${REFROOT}/fasta/all_contaminants.fa"
```

#### 6) Add biomart and GSEA files
```bash
cp /data/src/biomart/"${SPECIES}.mart_export.txt" "${REFROOT}/biomart/"
cp /data/src/annotations/ERCC92_genes.txt "${REFROOT}/annotations/"
cp /data/src/GSEA/*.gmt "${REFROOT}/GSEA/"
```

#### 7) Quick completeness and naming check
```bash
test -f "${REFROOT}/annotations/${SPECIES}.gtf" && \
test -f "${REFROOT}/annotations/${SPECIES}.bed" && \
test -f "${REFROOT}/annotations/ERCC92_genes.txt" && \
test -f "${REFROOT}/fasta/${SPECIES}.fa" && \
test -f "${REFROOT}/fasta/${SPECIES}.fa.fai" && \
test -f "${REFROOT}/fasta/${SPECIES}.transcriptome.fa" && \
test -f "${REFROOT}/fasta/${SPECIES}.transcriptome.fa.fai" && \
test -f "${REFROOT}/STAR/Genome" && \
test -f "${REFROOT}/STAR/SA" && \
test -f "${REFROOT}/STAR/SAindex" && \
test -f "${REFROOT}/STAR/chrLength.txt" && \
test -f "${REFROOT}/STAR/chrName.txt" && \
test -f "${REFROOT}/STAR/chrNameLength.txt" && \
test -f "${REFROOT}/STAR/chrStart.txt" && \
test -f "${REFROOT}/STAR/exonGeTrInfo.tab" && \
test -f "${REFROOT}/STAR/exonInfo.tab" && \
test -f "${REFROOT}/STAR/geneInfo.tab" && \
test -f "${REFROOT}/STAR/genomeParameters.txt" && \
test -f "${REFROOT}/STAR/sjdbInfo.txt" && \
test -f "${REFROOT}/STAR/sjdbList.fromGTF.out.tab" && \
test -f "${REFROOT}/STAR/sjdbList.out.tab" && \
test -f "${REFROOT}/STAR/transcriptInfo.tab" && \
test -f "${REFROOT}/kallisto/${SPECIES}.kallisto.idx" && \
test -f "${REFROOT}/kallisto/all_contaminants.kallisto.idx" && \
test -f "${REFROOT}/biomart/${SPECIES}.mart_export.txt" && \
ls "${REFROOT}/GSEA/"*.gmt >/dev/null
echo "Reference check passed for ${SPECIES}"
```

You can also print the final structure:
```bash
ls -R "${REFROOT}"
```

---

### One-shot Docker command (generate all references for one species)
This does the same workflow in one command (non-interactive). Update host paths and source file names first.

```bash
./generate_reference_one_shot.sh
```

Edit variables at the top of `generate_reference_one_shot.sh` before running:
- `HOST_REFERENCEFILES_DIR`
- `HOST_RAW_REFERENCE_DOWNLOADS`
- `SPECIES`
- source filename variables (for FASTA/GTF/Biomart/ERCC)
