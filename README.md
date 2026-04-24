# RNA-Seq Pipeline Usage Guide 

## Table of Contents
- [General Information](#general-information)
- [Getting Started](#getting-started)
- [Running the Pipeline on Docker](#running-on-docker)
- [Adding a New Species (Reference Files)](#adding-a-new-species)

## General Information

### What this pipeline does
This RNA-seq pipeline processes raw sequencing data through differential analysis and reporting. Key analytical outputs include:

- Read quality assessment (`FastQC`, `fastp`)
- Genome alignment (`STAR`)
- Strandedness inference (`RSeQC`)
- Transcript abundance quantification (`kallisto`) and gene-level summarization
- Differential expression analysis (`DESeq2`)
- Optional ERCC-based normalization (if ERCC spike-ins are present)
- Transcript-level differential analysis (`Sleuth`)
- Gene set enrichment analysis (`GSEA`, MSigDB gene set collections)
- Alternative splicing analysis (`rMATS`)
- Final report generation (`report.html`) and organized output directories

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
The pipeline uses the reference directory path you pass in commands.
In this README, examples from inside the Docker Container use:

`/data/referenceFiles`

Supported species currently include:
- `hg38`
- `hg38plusAkataInverted`
- `hg38plusKSHV`
- `hg38plusKSHVALT`
- `mm39`
- `mm39plusMHV68`

Each reference scontains assets used by the pipeline (for example STAR index, annotation files, kallisto index, Biomart export, and GSEA gene sets). 

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
- `<scripts_dir>`: The `./scripts` folder in this directory.
- `<results_dir>`: Same argument as metadata.py. 

and will naturally exist inside the `/work` folder in the container. 

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
When you use `-v`, Docker maps folders from your computer to new paths inside the container.  
Use the **container paths** (right side of each `-v`) when calling `metadata.py` and `rnaseq.py` inside the container.

- `fastq_root_dir`: `/data/Model_Experiment`
- `reference_dir`: `/data/referenceFiles`
- `results_dir`: `/data/output`
- `scripts_dir`: `/work/scripts` (if helper scripts are in a `scripts/` subfolder), or `/work` (if helper scripts stay in repo root)

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
   python rnaseq.py /data/output/ethan_metadata_04232026_174629.tsv /data/referenceFiles /work/scripts_for_rnaseq /data/output
   ```
 

### 2) One-shot Docker command (non-interactive)
```bash
docker run --rm \
  -v "/path/to/referenceFiles/on/your/computer:/data/referenceFiles:ro" \
  -v "/path/to/Model_Experiment/on/your/computer:/data/Model_Experiment:ro" \
  -v "/path/to/output/folder/on/your/computer:/data/output" \
  -w /work \
  rnaseqpipeline:latest \
  bash -lc 'conda activate rnaseqpipeline && python metadata.py /data/Model_Experiment /data/referenceFiles <species_name> <investigator_name> <PE|SE> /data/output && python rnaseq.py /data/output/<investigator_name>_metadata_*.tsv /data/referenceFiles /work/scripts_for_rnaseq /data/output'
```

### Checklist
- Parent folder naming follows `Model_Experiment` (single underscore).
- Subdirectories are exactly `cntl` and `test`.
- Reference name exists in `referenceFiles`.
- Environment is `rnaseqpipeline` (native or container).
- Metadata and rnaseq are run with matching paths in the selected environment.

## Adding a New Species 
