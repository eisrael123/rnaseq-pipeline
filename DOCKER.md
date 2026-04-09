# Running the RNA-seq pipeline in Docker

This guide is for someone new to Docker who wants to build the image, attach reference data without putting it in the image, and run **`metadata.py`** followed by **`rnaseq.py`**.

## What Docker is doing here

- **Image**: A recipe + installed software (Python, conda, tools from `rnaseqpipeline.yml`, and your scripts copied into the image). It does **not** include your large reference genomes or indexes unless you add them on purpose.
- **Container**: A running instance of that image. When you **mount** (bind-mount) a folder from your computer into the container, programs inside the container see that data at a path you choose. Nothing is copied into the image at build time.

You will typically:

1. Build the image once.
2. Run a container with **volumes** that provide: reference files, input FASTQs, and an output directory.
3. Run **`metadata.py`** then **`rnaseq.py`** inside the container using paths that exist **inside** the container.

---

## Prerequisites

1. **Docker Desktop** (macOS) or Docker Engine installed and **running** (`docker ps` should work, not “Cannot connect to the Docker daemon”).
2. On macOS, ensure **File Sharing** in Docker Desktop includes the directories you will mount (for example your project folder and the folder that holds reference files).
3. This repository’s **`scripts/`** directory (contains `Dockerfile`, `rnaseqpipeline.yml`, and the Python/R scripts).

---

## 1. Build the image

From the directory that contains the `Dockerfile` (this repo’s `scripts/` folder):

```bash
cd /path/to/rnaseq-pipeline/scripts
docker build -t rnaseqpipeline:latest .
```

- **`-t rnaseqpipeline:latest`** assigns a name and tag so you can run `rnaseqpipeline:latest` later.
- The build reads `Dockerfile`, creates the conda env from `rnaseqpipeline.yml`, and copies the repo into `/work` in the image.

---

## 2. Paths the pipeline expects (reference layout)

Both **`metadata.py`** and **`rnaseq.py`** use a fixed base directory for references:

```text
/Applications/ngs/pipelines/rnaseq/referenceFiles
```

Under that base, for each genome build you use a **directory whose name equals the “reference” / species key** (for example `hg38plusAkataInverted`). The pipeline expects a layout like:

```text
referenceFiles/<reference_name>/
  STAR/                    # STAR genome index (genomeDir)
  STAR/chrNameLength.txt   # used when converting wiggle to bigWig
  annotations/
    <reference_name>.gtf
    <reference_name>.bed
  kallisto/
    <reference_name>.kallisto.idx
    all_contaminants.kallisto.idx   # optional, for contaminant kallisto runs
  biomart/
    <reference_name>.mart_export.txt
  GSEA/
    *.gmt                  # gene set files for GSEA
```

- **`metadata.py`** checks that **`referenceFiles/<reference_name>`** exists (it lists subfolders of `referenceFiles` to validate your chosen reference name).
- **`rnaseq.py`** reads STAR, GTF, BED, kallisto indexes, Biomart TSV, GSEA GMTs, and `chrNameLength.txt` from paths built from that same base + `<reference_name>`.

**Important:** The folder name **`reference_name`** must match **`Species`** in the metadata table (metadata sets `Species` from that argument).

---

## 3. Mounting reference files so paths line up

Because the code uses the absolute path `/Applications/ngs/pipelines/rnaseq/referenceFiles`, the straightforward approach **without changing Python** is:

**Mount your host reference directory at exactly that path inside the container** (read-only is enough for references):

```bash
-v "/path/on/your/computer/referenceFiles:/Applications/ngs/pipelines/rnaseq/referenceFiles:ro"
```

Replace `/path/on/your/computer/referenceFiles` with the real folder that contains your `hg38plusAkataInverted` (etc.) subdirectories. After this mount, the container sees the same tree the scripts expect.

---

## 4. Other hardcoded paths (needed for a *working* full pipeline)

The Python scripts assume additional absolute locations that **do not exist** in a stock Linux container unless you create them:

| Topic | What the code expects | Practical approach |
|--------|------------------------|---------------------|
| R and Perl helpers | Paths under `/Applications/ngs/pipelines/rnaseq/scripts/` (e.g. DESeq2/Sleuth R scripts, `reform_*.py`, Perl splice script) | **Recommended:** extend the `Dockerfile` with something like: `RUN mkdir -p /Applications/ngs/pipelines/rnaseq && ln -sfn /work /Applications/ngs/pipelines/rnaseq/scripts` so `/work` (where the image copies your repo) is also visible at the path the pipeline uses. |
| STAR temporary directory | A Mac-specific path under `/Users/mac14/Desktop/...` | **Recommended:** `RUN mkdir -p "/Users/mac14/Desktop/Misc_Desktop_Folders"` in the `Dockerfile` so STAR can write temp files inside the container. |

Without these, alignment or downstream R steps can fail even if references mount correctly.

### Conda environment check

`metadata.py` and `rnaseq.py` require the environment name **`rnaseqpipeline`** via `CONDA_DEFAULT_ENV`. The image puts that env on `PATH`, but **`CONDA_DEFAULT_ENV` may be unset** when you run commands.

**Recommended:** add to the `Dockerfile` (or pass when running):

```dockerfile
ENV CONDA_DEFAULT_ENV=rnaseqpipeline
```

Or on the command line:

```bash
-e CONDA_DEFAULT_ENV=rnaseqpipeline
```

---

## 5. FASTQs and outputs: paths must stay inside the container

The metadata TSV stores **absolute paths** to FASTQs (`Path Read 1`, `Path Read 2`). **`rnaseq.py` runs inside the container**, so those paths must be valid **there**.

**Rule:** Run **`metadata.py` inside the same container setup** (same volume mounts) as **`rnaseq.py`**, and use **container paths** for the parent directory and output directory.

**Example layout:**

- Host: `~/data/my_project/fastqs/` → mount as container: `/data/fastqs`
- Host: `~/data/my_project/output/` → mount as container: `/data/output`

Then run metadata with **`/data/fastqs/...`** (not `~/...` from the host) if execution is inside the container.

---

## 6. `metadata.py` usage

Intended usage:

```text
python metadata.py <parent_dir> <reference_name> <investigator> <PE|SE> <output_dir>
```

- **`parent_dir`**: Directory that **directly contains** FASTQ folders. Its **basename must be exactly `something_something`** (one underscore), e.g. `Mutu_siCNOT-Zta`.
- **`reference_name`**: Must match a subdirectory name under `referenceFiles` and becomes **`Species`** in the TSV.
- **`output_dir`**: Where the metadata TSV and log are written (use a mounted directory if you want files on the host).

---

## 7. `rnaseq.py` usage

```text
python rnaseq.py <metadata.tsv> <results_dir>
```

- **`metadata.tsv`**: Path to the file produced by `metadata.py` (inside the container).
- **`results_dir`**: Where pipeline outputs go; typically the **same** as `output_dir` from metadata, or a subfolder you prefer.

Run this from a working directory that allows the script to find any files it moves by relative name (the pipeline expects a consistent layout relative to `results_dir`).

---

## 8. Example: interactive shell, then run both steps

Adjust host paths to match your machine.

```bash
docker run --rm -it \
  -v "/Users/you/ethan/rnaseq/referenceFiles:/Applications/ngs/pipelines/rnaseq/referenceFiles:ro" \
  -v "/Users/you/data/Mutu_siCNOT-Zta:/data/Mutu_siCNOT-Zta:ro" \
  -v "/Users/you/data/output_run:/data/output" \
  -w /work \
  rnaseqpipeline:latest \
  bash
```

Inside the container:

```bash
conda activate rnaseqpipeline
python metadata.py /data/Mutu_siCNOT-Zta hg38plusAkataInverted myname PE /data/output
python rnaseq.py /data/output/myname_metadata_<timestamp>.tsv /data/output
exit
```

Use the **actual** metadata filename printed by `metadata.py`.

---

## 9. One-shot run (non-interactive)

Replace paths and the metadata filename:

```bash
docker run --rm \
  -e CONDA_DEFAULT_ENV=rnaseqpipeline \
  -v "/host/ref/referenceFiles:/Applications/ngs/pipelines/rnaseq/referenceFiles:ro" \
  -v "/host/fastq_root:/data/fastqs:ro" \
  -v "/host/output:/data/output" \
  -w /work \
  rnaseqpipeline:latest \
  bash -lc 'python metadata.py /data/fastqs/MyModel_MyExperiment refName investigator PE /data/output && python rnaseq.py /data/output/investigator_metadata_*.tsv /data/output'
```

Using a glob for the metadata file only works if exactly one file matches; otherwise pass the full filename.

---

## 10. Troubleshooting

| Symptom | What to check |
|---------|----------------|
| `Cannot connect to the Docker daemon` | Start Docker Desktop / Docker service. |
| Reference not found | Mount path must be exactly `.../referenceFiles` → `/Applications/ngs/pipelines/rnaseq/referenceFiles`; folder names must match `reference_name`. |
| `CONDA_DEFAULT_ENV` error | Set `ENV CONDA_DEFAULT_ENV=rnaseqpipeline` in the image or `-e CONDA_DEFAULT_ENV=rnaseqpipeline` at run. |
| FASTQs not found in `rnaseq.py` | Metadata was probably generated with host paths; regenerate metadata **inside** the container with the same `-v` mounts. |
| R/Perl “file not found” under `/Applications/ngs/...` | Add the `/work` → `/Applications/ngs/pipelines/rnaseq/scripts` symlink (or equivalent) as in section 4. |
| STAR fails on temp directory | Create the `Misc_Desktop_Folders` path in the image or change the temp path in code to `/tmp`. |

---

## 11. Summary checklist

1. Build: `docker build -t rnaseqpipeline:latest .` from `scripts/`.
2. Ensure **`CONDA_DEFAULT_ENV=rnaseqpipeline`** when running Python.
3. Mount references: host `referenceFiles` → `/Applications/ngs/pipelines/rnaseq/referenceFiles`.
4. Align reference **directory names** and **file names** with the layout in section 2.
5. Mount FASTQs and outputs; run **`metadata.py` and `rnaseq.py` in the container** using **container paths**.
6. Add Dockerfile fixes for **`/Applications/ngs/pipelines/rnaseq/scripts`** and STAR temp dir (section 4) so the full pipeline can run.
