#!/usr/bin/env python3

# Usage: rnaseq.py <reference_dir> <scripts_dir> <results_dir> [options]
#
# Reads and updates <results_dir>/metadata.tsv in place -- the file metadata.py must already
# have generated there. results_dir receives the canonical layout documented in README.md
# ("Outputs").

import subprocess
import os
import platform
from pathlib import Path
import pandas as pd
import sys
import resource
import time
soft, hard = 10000, 10000
resource.setrlimit(resource.RLIMIT_NOFILE, (soft, hard))
import shutil
import tempfile
import json
import matplotlib.pyplot as plt
from matplotlib import rcParams
import re
import argparse

# The helper modules live next to this file in both a checkout and the container image.
sys.path.insert(0, str(Path(__file__).resolve().parent / "rnaseq_helper_scripts"))

import bigwig
import checksums as checksums_module
import provenance
import strandedness as strandedness_module
import tables_de
import tables_expression
import tables_gsea
import tables_junction
import tables_qc
import tables_rmats
import validate_outputs
from annotation import Annotation
from outputs import RunContext, create_layout, stage_log_path
from schemas import PIPELINE_VERSION

REQUIRED_ENV = "rnaseqpipeline"

# Alignment settings shared by the paired- and single-end invocations, kept in one place so the
# manifest can record exactly what STAR was asked to do.
STAR_COMMON_ARGS = [
    "--outFilterMismatchNmax", "3",
    "--outSAMstrandField", "intronMotif",
    "--chimSegmentMin", "2",
    "--alignSJDBoverhangMin", "6",
    "--alignEndsType", "EndToEnd",
    "--alignIntronMax", "300000",
    "--outWigType", "wiggle",
    "--outWigNorm", "None",
    "--outSAMtype", "BAM", "SortedByCoordinate",
]
KALLISTO_COMMON_ARGS = ["-b", "100", "-t", "10"]
RMATS_COMMON_ARGS = ["--variable-read-length", "--nthread", "10"]
DESEQ2_ALPHA = 0.05

# Check if running in the right conda env
conda_env = os.environ.get("CONDA_DEFAULT_ENV")
if conda_env != REQUIRED_ENV:
    sys.stderr.write(
        f"\nERROR: This script must be run inside the '{REQUIRED_ENV}' conda environment.\n"
        f"Currently active environment: {conda_env or 'None'}\n"
        f"Please run:\n\n    conda activate {REQUIRED_ENV}\n\n"
    )
    sys.exit(1)

def sample_sort_key(name: str):
    # natural sort on trailing integer if present
    m = re.search(r"(\d+)(?!.*\d)", name)
    n = int(m.group(1)) if m else 10**9
    return (re.sub(r"\d+", "", name), n, name)

def ordered_sample_names(metadata):
    """Controls first, then tests, each in the order of metadata.tsv's ``replicate_index``.

    That index is what metadata.py baked into every sample_id (filename order), and
    deseq2_metadata.tsv follows it, so DESeq2's counts columns must too. Sorting by the trailing
    integer instead (sample_sort_key) picks up the lane in names like ``DDX5-2X-KO-1_L05`` and
    put replicate 2 (lane 4) before replicate 1 (lane 5), which DESeq2 rejects. sample_sort_key
    remains only as a fallback for metadata written before replicate_index existed.
    """
    cols = ["Sample name", "condition"]
    has_index = "replicate_index" in metadata.columns
    if has_index:
        cols.append("replicate_index")
    order = (
        metadata[cols]
        .drop_duplicates(subset=["Sample name"])
        .assign(_group=lambda d: d["condition"].str.strip().str.lower()
                .map(lambda x: 0 if x == "cntl" else 1))
    )
    if has_index:
        order = order.sort_values(by=["_group", "replicate_index"], kind="stable")
    else:
        order = order.sort_values(
            by=["_group", "Sample name"],
            key=lambda s: s.map(sample_sort_key) if s.name == "Sample name" else s,
        )
    return order["Sample name"].tolist()

def setup_logging(results_dir):
    log_file = stage_log_path(results_dir, "pipeline")
    sys.stdout = open(log_file, 'w')
    sys.stderr = sys.stdout

def log_total_time(start_time, results_dir):
    """
    Log the total processing time to logs/pipeline.log.

    Parameters:
    - start_time (float): The start time of the pipeline (from `time.time()`).
    - results_dir (str or Path): The directory where the results are stored.
    """
    # Calculate the total processing time
    end_time = time.time()
    total_time_seconds = end_time - start_time
    total_time_minutes = total_time_seconds / 60
    total_time_hours = total_time_minutes / 60

    # Format the processing time
    formatted_time = f"{int(total_time_hours)}h {int(total_time_minutes % 60)}m {int(total_time_seconds % 60)}s"

    # Define the path to the processing log file
    processing_log_file = stage_log_path(results_dir, "pipeline")

    # Write the processing time to the log file
    try:
        with open(processing_log_file, "a") as log_file:
            log_file.write(f"Total processing time:\t{formatted_time}\n")
        print(f"Total processing time logged: {formatted_time}")
    except Exception as e:
        print(f"Error writing total processing time to log file: {e}")

def run_fastqc(fastq1, fastq2, results_dir):
    """Run FastQC on FASTQ files."""
    os.makedirs(results_dir, exist_ok=True)
    if fastq2 is None:
        cmd = f"fastqc {fastq1} --threads 8 --outdir={results_dir}"
    else:
        cmd = f"fastqc {fastq1} {fastq2} --threads 8 --outdir={results_dir}"
    subprocess.run(cmd, shell=True, check=True)

def run_fastp_no_trimming(fastq1, fastq2, html_out, json_out):
    """Run fastp without trimming to analyze PCR duplication rates."""
    os.makedirs(os.path.dirname(html_out), exist_ok=True)
    os.makedirs(os.path.dirname(json_out), exist_ok=True)

    if fastq2 is None:
        cmd = f"fastp -w 10 -i {fastq1} --html {html_out} --json {json_out} --disable_trim_poly_g --disable_length_filtering"
    else:
        cmd = f"fastp -w 10 -i {fastq1} -I {fastq2} --detect_adapter_for_pe --html {html_out} --json {json_out} --disable_trim_poly_g --disable_length_filtering"
    subprocess.run(cmd, shell=True, check=True)

def run_star(sample_name, fastq1, fastq2, star_results_dir, species_name, is_paired_end=True):
    """Run STAR alignment."""
    os.makedirs(star_results_dir, exist_ok=True)
    decompression_command = "zcat" if platform.system() == "Linux" else "gzcat"
    genome_dir = REFERENCE_DIR / species_name / "STAR"
    gtf_file = REFERENCE_DIR / species_name / "annotations" / f"{species_name}.gtf"
    # macOS: legacy Desktop temp dir (ensure parent exists); other OS: system temp
    if platform.system() == "Darwin":
        _star_tmp_parent = Path("/Users/mac14/Desktop/Misc_Desktop_Folders")
        _star_tmp_parent.mkdir(parents=True, exist_ok=True)
        tmp_dir = str(_star_tmp_parent / f"rnaseqPipelineTmpDir_{sample_name}")
    else:
        tmp_dir = str(Path(tempfile.gettempdir()) / f"rnaseqPipelineTmpDir_{sample_name}")

    read_files = [fastq1, fastq2] if is_paired_end else [fastq1]
    star_cmd = [
        "STAR",
        "--genomeDir", str(genome_dir),
        "--readFilesIn", *read_files,
        "--readFilesCommand", decompression_command,
        "--outFileNamePrefix", f"{star_results_dir}/{sample_name}_",
        "--sjdbGTFfile", str(gtf_file),
        *STAR_COMMON_ARGS,
        f"--outTmpDir {tmp_dir}",
        "--limitBAMsortRAM 60000000000",
        "--runThreadN 10"
    ]

    try:
        subprocess.run(star_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running STAR for sample {sample_name}: {e}")
        sys.exit(1)
    finally:
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)  # Clean up the temporary directory after STAR run

def post_process_alignment(sample_name, star_results_dir):
    """
    Post-process the alignment by deleting the STARgenome directory,
    and indexing the sorted BAM file using samtools.

    Parameters:
    - sample_name (str): The name of the sample.
    - star_results_dir (str or Path): The output directory where the alignment files are stored.
    """
    star_results_dir = Path(star_results_dir)
    #unsorted_bam = star_results_dir / f"{sample_name}_Aligned.out.bam"
    sorted_bam = star_results_dir / f"{sample_name}_Aligned.sortedByCoord.out.bam"
    stargenome_dir = star_results_dir / f"{sample_name}__STARgenome"

    # Delete the unsorted BAM file
    #if unsorted_bam.exists():
        #unsorted_bam.unlink()
        #print(f"Deleted unsorted BAM file: {unsorted_bam}")

    # Delete the STARgenome directory
    if stargenome_dir.exists() and stargenome_dir.is_dir():
        shutil.rmtree(stargenome_dir)
        print(f"Deleted STARgenome directory: {stargenome_dir}")

    # Index the sorted BAM file using samtools
    cmd = f"samtools index {sorted_bam}"
    subprocess.run(cmd, shell=True, check=True)
    print(f"Indexed sorted BAM file: {sorted_bam}")

    return sorted_bam

def infer_strandedness(sorted_bam, species, results_dir):
    """Infer library strandedness with RSeQC and return the canonical library type.

    The raw RSeQC output is kept verbatim for the record; the interpretation is appended as a
    single machine-readable line. Every consumer reads it back through
    ``strandedness.parse_rseqc`` rather than string-matching this file.
    """
    rseqc_dir = results_dir / "rseqc"
    os.makedirs(rseqc_dir, exist_ok=True)

    bed_file = REFERENCE_DIR / species_name / "annotations" / f"{species_name}.bed"
    infer_command = ["infer_experiment.py", "-r", str(bed_file), "-i", str(sorted_bam), "-s", "5000000"]
    output_file = rseqc_dir / f"{sorted_bam.stem}_strandedness.txt"

    print(f"Running RSeQC infer_experiment.py: {' '.join(infer_command)}")
    result = subprocess.run(infer_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

    parsed = strandedness_module.parse_rseqc(result.stdout)
    with open(output_file, "w") as outfile:
        outfile.write(result.stdout)
        outfile.write(f"\n{strandedness_module.CANONICAL_LINE} {parsed['library_type']}\n")

    if parsed["library_type"] == strandedness_module.UNDETERMINED:
        raise RuntimeError(
            f"RSeQC could not determine strandedness for {sorted_bam.name}. Strand-resolved "
            f"bigWigs and rMATS both depend on it, so the run is stopping rather than "
            f"guessing. See {output_file}."
        )
    return parsed["library_type"]

def run_kallisto(fastq1, fastq2, sample_name, species_name, strandedness, results_dir, is_contaminant_run=False):
    """Run kallisto quant to measure abundance for species or contaminants."""
    # Set the correct output subdirectory with sampleName prefix
    output_subdir = f"{sample_name}_Foreign_sequences" if is_contaminant_run else f"{sample_name}_{species_name}"
    kallisto_dir = results_dir / "kallisto" / sample_name / output_subdir
    kallisto_dir.mkdir(parents=True, exist_ok=True)

    # Use the appropriate kallisto index based on species or contaminants.
    # On macOS (local runs), use the ORIGINAL idx that was built. On Docker, use the new idx.
    if is_contaminant_run:
        if platform.system() == "Darwin": #Mac
            kallisto_index = REFERENCE_DIR / species_name / "kallisto" / "all_contaminants_ORIGINAL.kallisto.idx"
            print(f"Kallisto contaminants: using index {kallisto_index}")
        else: #Docker
            kallisto_index = REFERENCE_DIR / species_name / "kallisto" / "all_contaminants.kallisto.idx"
            print(f"Kallisto contaminants: using index {kallisto_index}")
    else:
        if platform.system() == "Darwin": #Mac
            kallisto_index = REFERENCE_DIR / species_name / "kallisto" / f"{species_name}_ORIGINAL.kallisto.idx"
            print(f"Kallisto species: using index {kallisto_index}")
        else: #Docker
            kallisto_index = REFERENCE_DIR / species_name / "kallisto" / f"{species_name}.kallisto.idx"
            print(f"Kallisto species: using index {kallisto_index}")

    stranded_option = strandedness_module.kallisto_option(strandedness)

    # Configure the kallisto command based on PE/SE and strandedness
    kallisto_cmd = [
        "kallisto", "quant", "-i", str(kallisto_index),
        "-o", str(kallisto_dir), *KALLISTO_COMMON_ARGS
    ]
    if stranded_option:
        kallisto_cmd.append(stranded_option)
    # Only treat as paired-end if fastq2 is a non-empty string and not 'nan'
    if isinstance(fastq2, str) and fastq2 and fastq2.lower() != 'nan':  # PE data
        kallisto_cmd.extend([fastq1, fastq2])
    else:  # SE data
        kallisto_cmd.extend(["--single", "-l", "200", "-s", "20", fastq1])

    # Log and run the kallisto command
    log_file = kallisto_dir / f"{sample_name}.kallisto.log"
    err_file = kallisto_dir / f"{sample_name}.kallisto.err"
    with open(log_file, "w") as out, open(err_file, "w") as err:
        print(f"Running kallisto quant for {sample_name} {'contaminants' if is_contaminant_run else 'species'}")
        subprocess.run(kallisto_cmd, stdout=out, stderr=err, check=True)

def summarize_to_gene_level(sample_name, species, results_dir):
    """Summarize transcript-level abundance to gene-level abundance."""

    # Path to the mart_file
    mart_file = REFERENCE_DIR / species / "biomart" / f"{species}.mart_export.txt"

    # Load the transcript-to-gene mapping
    if not mart_file.exists():
        print(f"Mart file {mart_file} not found.")
        sys.exit(1)
    tx2gene = pd.read_csv(mart_file, sep='\t')

    # Check the column names in the mart file and assume first column is transcript, second is gene
    tx2gene_columns = tx2gene.columns.tolist()
    print(f"Mart file columns: {tx2gene_columns}")

    # Rename columns to standard names if needed
    if len(tx2gene_columns) >= 2:
        tx2gene.columns = ['transcript_id', 'gene'] + tx2gene_columns[2:]
    else:
        print(f"Error: Mart file should have at least 2 columns, found {len(tx2gene_columns)}")
        return

    # Path to the kallisto output directory for species
    kallisto_dir = results_dir / "kallisto" / sample_name / f"{sample_name}_{species}"

    # Path to the abundance.tsv file
    abundance_file = kallisto_dir / "abundance.tsv"
    if not abundance_file.exists():
        print(f"Abundance file {abundance_file} not found.")
        return

    # Read the abundance.tsv file
    abundance = pd.read_csv(abundance_file, sep='\t')
    print(f"Loaded abundance data with {len(abundance)} transcripts")

    # Show some EBV transcripts in abundance data for debugging
    ebv_transcripts = abundance[abundance['target_id'].str.contains('BVLF1_1|BVRF2_1|RPMS1_1', na=False)]
    if not ebv_transcripts.empty:
        print(f"Found {len(ebv_transcripts)} EBV transcripts in abundance data")
        print(f"Sample EBV transcripts: {ebv_transcripts['target_id'].head().tolist()}")
    else:
        print("No EBV transcripts found in abundance data")

    # Try direct matching first (transcript IDs in mart file match exactly)
    merged = pd.merge(abundance, tx2gene, how='left', left_on='target_id', right_on='transcript_id')

    # Check if direct matching worked for a few transcripts
    direct_match_count = merged['gene'].notna().sum()
    print(f"Direct matching resulted in {direct_match_count} successful mappings")

    # If direct matching fails for many transcripts, try with suffix removal
    if direct_match_count < len(abundance) * 0.5:  # If less than 50% matched
        print("Direct matching failed for many transcripts, trying with suffix removal...")

        # Create a mapping column by removing suffix from target_id for matching
        # This handles cases where abundance has "transcript_1" but mart has "transcript"
        abundance['transcript_for_mapping'] = abundance['target_id'].str.replace(r'_\d+$', '', regex=True)

        # Debug: show the mapping for some EBV transcripts
        ebv_mapping_debug = abundance[abundance['target_id'].str.contains('BVLF1_1|BVRF2_1|RPMS1_1', na=False)][['target_id', 'transcript_for_mapping']].head()
        if not ebv_mapping_debug.empty:
            print("EBV transcript mapping examples:")
            print(ebv_mapping_debug.to_string())

        # Merge abundance data with transcript-to-gene mapping using modified transcript IDs
        merged = pd.merge(abundance, tx2gene, how='left', left_on='transcript_for_mapping', right_on='transcript_id')
    else:
        print("Direct matching was successful, using original transcript IDs")
    print(f"After merge: {len(merged)} rows")

    # Check specifically for EBV transcripts after merge
    ebv_merged = merged[merged['target_id'].str.contains('BVLF1_1|BVRF2_1|RPMS1_1', na=False)]
    if not ebv_merged.empty:
        print(f"Found {len(ebv_merged)} EBV transcripts after merge")
        print(f"EBV genes mapped: {ebv_merged['gene'].dropna().unique().tolist()}")
        # Show details of EBV mapping
        print("EBV transcript mapping details:")
        mapping_cols = ['target_id', 'gene', 'est_counts', 'tpm']
        if 'transcript_for_mapping' in ebv_merged.columns:
            mapping_cols.insert(1, 'transcript_for_mapping')
        print(ebv_merged[mapping_cols].to_string())
    else:
        print("No EBV transcripts found after merge")

    # Check for transcripts without gene mapping
    if merged['gene'].isnull().any():
        unmapped_count = merged['gene'].isnull().sum()
        print(f"Warning: {unmapped_count} transcripts in {sample_name} do not have a gene mapping.")
        # Show some unmapped transcripts for debugging
        unmapped_transcripts = merged[merged['gene'].isnull()]['target_id'].head(10).tolist()
        print(f"Sample unmapped transcripts: {unmapped_transcripts}")

        # Check if any of the unmapped are EBV
        unmapped_ebv = merged[(merged['gene'].isnull()) & (merged['target_id'].str.contains('BVLF1_1|BVRF2_1|RPMS1_1', na=False))]
        if not unmapped_ebv.empty:
            print(f"WARNING: {len(unmapped_ebv)} EBV transcripts were not mapped!")
            debug_cols = ['target_id']
            if 'transcript_for_mapping' in unmapped_ebv.columns:
                debug_cols.append('transcript_for_mapping')
            print(unmapped_ebv[debug_cols].to_string())

    # Remove rows where gene mapping failed (keep only successfully mapped transcripts)
    mapped_data = merged.dropna(subset=['gene'])
    print(f"Successfully mapped {len(mapped_data)} transcripts to genes")

    # Group by gene and sum est_counts and tpm
    gene_counts = mapped_data.groupby('gene').agg({'est_counts': 'sum'}).reset_index()
    gene_tpm = mapped_data.groupby('gene').agg({'tpm': 'sum'}).reset_index()

    # Save gene-level est_counts
    gene_counts_file = kallisto_dir / "abundance_gene_est_counts.tsv"
    gene_counts.to_csv(gene_counts_file, sep='\t', index=False)
    print(f"Saved gene-level est_counts to {gene_counts_file} ({len(gene_counts)} genes)")

    # Save gene-level tpm
    gene_tpm_file = kallisto_dir / "abundance_gene_tpm.tsv"
    gene_tpm.to_csv(gene_tpm_file, sep='\t', index=False)
    print(f"Saved gene-level tpm to {gene_tpm_file} ({len(gene_tpm)} genes)")

    # Final check: verify EBV genes are in the output
    ebv_genes_in_counts = gene_counts[gene_counts['gene'].str.contains('BVLF1|BVRF2|RPMS1', na=False)]
    ebv_genes_in_tpm = gene_tpm[gene_tpm['gene'].str.contains('BVLF1|BVRF2|RPMS1', na=False)]

    if not ebv_genes_in_counts.empty:
        print(f"SUCCESS: Found {len(ebv_genes_in_counts)} EBV genes in final est_counts file")
        print(f"EBV genes in counts: {ebv_genes_in_counts['gene'].tolist()}")
    else:
        print("WARNING: No EBV genes found in final est_counts file")

    if not ebv_genes_in_tpm.empty:
        print(f"SUCCESS: Found {len(ebv_genes_in_tpm)} EBV genes in final tpm file")
        print(f"EBV genes in tpm: {ebv_genes_in_tpm['gene'].tolist()}")
    else:
        print("WARNING: No EBV genes found in final tpm file")

def generate_deseq2_metadata(metadata, results_dir):
    """Generate DESeq2 metadata file for differential expression analysis."""
    deseq2_dir = Path(results_dir) / "deseq2"
    deseq2_dir.mkdir(parents=True, exist_ok=True)
    deseq2_metadata_file = deseq2_dir / "deseq2_metadata.tsv"

    # Check for possible columns that can represent 'condition'
    condition_column = None
    for col in ['condition', 'Condition', 'group', 'Group', 'Test', 'test']:
        if col in metadata.columns:
            condition_column = col
            break

    if condition_column is None:
        print("Error: The metadata file must contain a 'condition' or similar column.")
        sys.exit(1)

    # Extract the relevant columns and drop duplicates
    deseq2_metadata = metadata[['Sample name', condition_column]].drop_duplicates().copy()
    deseq2_metadata.rename(columns={condition_column: 'condition'}, inplace=True)

    # Save the DESeq2 metadata file
    deseq2_metadata.to_csv(deseq2_metadata_file, sep='\t', index=False)
    print(f"DESeq2 metadata file saved to {deseq2_metadata_file}")

def generate_counts_matrix(samples, species, results_dir):
    """Generate a counts matrix for DESeq2 analysis."""
    deseq2_dir = Path(results_dir) / "deseq2"
    deseq2_dir.mkdir(parents=True, exist_ok=True)
    counts_matrix_file = deseq2_dir / "deseq2_counts_matrix.tsv"

    counts_list = []
    for sample_name in samples:
        # Path to the gene-level counts file
        counts_file = (Path(results_dir) / "kallisto" / sample_name / f"{sample_name}_{species}" / "abundance_gene_est_counts.tsv")
        if not counts_file.exists():
            print(f"Counts file {counts_file} not found for sample {sample_name}.")
            continue
        sample_counts = pd.read_csv(counts_file, sep='\t')
        sample_counts = sample_counts[['gene', 'est_counts']]
        sample_counts.rename(columns={'est_counts': sample_name}, inplace=True)
        counts_list.append(sample_counts)

    # Check if counts_list is empty
    if not counts_list:
        raise ValueError("No counts files found. Please check the input files and paths.")

    # Merge all counts data frames on 'gene'
    counts_matrix = counts_list[0]
    for df in counts_list[1:]:
        counts_matrix = pd.merge(counts_matrix, df, on='gene', how='outer')

    counts_matrix.fillna(0, inplace=True)
    counts_matrix.set_index('gene', inplace=True)
    counts_matrix = counts_matrix.astype(float)
    # Convert counts to integers by rounding
    counts_matrix = counts_matrix.round().astype(int)

    # Save the counts matrix to the deseq2 directory
    counts_matrix.to_csv(counts_matrix_file, sep='\t')
    print(f"Counts matrix saved to {counts_matrix_file}")

def run_deseq2_analysis(results_dir, species_name, reference_dir):
    """
    Execute the DESeq2 analysis using the deseq2_analysis_ercc.R script.

    Parameters:
    - results_dir (Path): Path object to the results directory containing input files.
    """
    # Define the path to the R script
    r_script = (scripts / "deseq2_analysis_ercc.R").resolve()
    if not r_script.exists():
        print("R script not found.")
        sys.exit(1)

    # Use absolute path: subprocess uses cwd=results_dir, so a relative R path would be wrong.
    cmd = [
        "Rscript",
        str(r_script),
        str(Path(results_dir).resolve()),
        species_name,
        str(Path(reference_dir).resolve()),
    ]
    print(f"Running DESeq2 analysis with command: {' '.join(cmd)}")  # Debug statement

    # Execute the R script
    try:
        # Run from results_dir so any implicit R outputs (e.g. Rplots.pdf) land there.
        subprocess.run(cmd, check=True, cwd=str(results_dir))
        print("DESeq2 analysis completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error running DESeq2 analysis: {e}")
        sys.exit(1)

def create_sleuth_metadata(metadata_file, sleuth_dir, species_name, results_dir):
    """Generate Sleuth metadata file."""

    # Construct the Sleuth directory path within the results directory
    sleuth_dir = Path(results_dir) / "sleuth"

    # Create the Sleuth directory if it doesn't exist
    sleuth_dir.mkdir(parents=True, exist_ok=True)

    # Define the output Sleuth metadata file path
    sleuth_metadata_file = sleuth_dir / "sleuth_metadata.tsv"

    # Read the metadata TSV file
    metadata = pd.read_csv(metadata_file, sep='\t')

    # Ensure the 'condition' column exists in the metadata
    if 'condition' not in metadata.columns:
        raise ValueError("The metadata file must contain a 'condition' column.")

    # Drop duplicate entries based on 'Sample name' and 'condition' columns
    metadata = metadata.drop_duplicates(subset=['Sample name', 'condition'])

    # Rename the 'Sample name' column to 'sample' as required by Sleuth
    metadata = metadata.rename(columns={'Sample name': 'sample'})

    # Recast the 'condition' column ('test'/'cntl') into Sleuth's own vocabulary
    # ('test'/'control'), case-insensitively.
    metadata['condition'] = metadata['condition'].apply(
        lambda x: 'control' if str(x).strip().lower() == 'cntl' else 'test'
    )

    # Update the 'path' column to point to the 'abundance.h5' files generated by Kallisto
    metadata['path'] = metadata['sample'].apply(
        lambda x: str(
            Path(sleuth_dir).parent / "kallisto" / str(x) / f"{x}_{species_name}" / "abundance.h5"
        )
    )

    # Select only the required columns for Sleuth metadata
    sleuth_metadata = metadata[['sample', 'condition', 'path']]

    # Verify that the 'condition' column contains only 'test' and 'control' levels
    if not set(sleuth_metadata['condition']).issubset({'test', 'control'}):
        raise ValueError("The 'condition' column must contain only 'test' and 'control' levels.")

    # Save the Sleuth metadata to a TSV file
    sleuth_metadata.to_csv(sleuth_metadata_file, sep='\t', index=False)
    print(f"Sleuth metadata file saved to {sleuth_metadata_file}")

def run_sleuth_analysis(results_dir, species_name, reference_dir):
    """Run Sleuth analysis using an R script."""
    r_script_path = (scripts / "sleuth_analysis_ercc.R").resolve()
    if not r_script_path.exists():
        print(f"R script not found: {r_script_path}")
        sys.exit(1)

    # Build the command to execute the R script
    command = [
        "Rscript",
        str(r_script_path),
        str(Path(results_dir).resolve()),
        species_name,
        str(Path(reference_dir).resolve()),
    ]
    print(f"Running Sleuth analysis with command: {' '.join(command)}")  # Debug statement

    # Execute the R script
    try:
        # Run from results_dir so any implicit R outputs (e.g. Rplots.pdf) land there.
        subprocess.run(command, check=True, cwd=str(results_dir))
        print("Sleuth analysis completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error running Sleuth analysis: {e}")
        sys.exit(1)

def create_gsea_input_file(metadata, results_dir, species_name):
    """Generate GSEA input file with collapsed TPMs from each sample."""

    # Define the path to kallisto within the specified output directory
    kallisto_dir = results_dir / "kallisto"

    # Define the GSEA directory
    gsea_dir = results_dir / "gsea"
    gsea_dir.mkdir(parents=True, exist_ok=True)
    print(f"GSEA directory created at '{gsea_dir}'.")

    # Load test and control samples from metadata
    control_samples = metadata[metadata['condition'].str.lower() == 'cntl']['Sample name'].tolist()
    test_samples = metadata[metadata['condition'].str.lower() == 'test']['Sample name'].tolist()

    # Combine control and test samples
    all_samples = control_samples + test_samples

    # Initialize a dictionary to store TPM dataframes for each sample
    sample_dfs = {}

    # Extract the test and control sample names for header ordering
    control_header = control_samples  # List of control sample names
    test_header = test_samples  # List of test sample names

    # Process each sample to load TPM data
    for sample in all_samples:
        sample_dir = kallisto_dir / str(sample) / f"{sample}_{species_name}"
        tpm_file = sample_dir / "abundance_gene_tpm.tsv"

        if tpm_file.exists():
            try:
                # Read TPM data; adjust 'gene' if your file has a different column name
                tpm_data = pd.read_csv(tpm_file, sep='\t', usecols=['gene', 'tpm'])

                # Rename 'tpm' column to the sample name
                tpm_data = tpm_data.rename(columns={'tpm': sample})

                # Set 'gene' as the index
                sample_dfs[sample] = tpm_data.set_index('gene')

                print(f"TPM data loaded for sample '{sample}'.")
            except Exception as e:
                print(f"Error reading TPM file '{tpm_file}': {e}")
        else:
            print(f"Warning: TPM file '{tpm_file}' not found.")

    # Concatenate all sample dataframes on 'gene' column
    if sample_dfs:
        try:
            merged_df = pd.concat(sample_dfs.values(), axis=1, join='outer').reset_index()
            print("Sample dataframes concatenated successfully.")
        except Exception as e:
            print(f"Error concatenating sample dataframes: {e}")
            sys.exit(1)

        # Verify the columns before renaming
        print("Columns before renaming:", merged_df.columns.tolist())

        # Fill NaN values with 0
        merged_df = merged_df.fillna(0)

        # Add 'DESCRIPTION' column with 'NA' values
        merged_df.insert(1, 'DESCRIPTION', 'NA')

        # Rename 'gene' column to 'NAME'
        if 'gene' in merged_df.columns:
            merged_df = merged_df.rename(columns={'gene': 'NAME'})
            print("Column 'gene' renamed to 'NAME'.")
        else:
            print("Error: 'gene' column not found in the merged DataFrame.")
            sys.exit(1)

        # Verify the columns after renaming
        print("Columns after renaming:", merged_df.columns.tolist())

        # Define the output file path within the GSEA directory
        gsea_input_file = gsea_dir / "test.txt"

        # Define the expected column order
        expected_columns = ['NAME', 'DESCRIPTION'] + control_header + test_header

        # Check if all expected columns are present
        missing_columns = [col for col in expected_columns if col not in merged_df.columns]
        if missing_columns:
            print(f"Error: Missing expected sample columns in merged DataFrame: {missing_columns}")
            sys.exit(1)

        # Reorder the columns
        merged_df = merged_df[expected_columns]
        print("Columns reordered successfully.")

        # Save the merged dataframe to the output file
        try:
            merged_df.to_csv(gsea_input_file, sep='\t', index=False)
            print(f"GSEA input file saved to '{gsea_input_file}'.")
        except Exception as e:
            print(f"Error saving GSEA input file '{gsea_input_file}': {e}")
            sys.exit(1)
    else:
        print("Error: No TPM data found for any samples.")
        sys.exit(1)

def create_gsea_cls_file(metadata, results_dir, species_name):
    """Generate GSEA class file (.cls) for input to GSEA."""

    # Define the GSEA directory
    gsea_dir = results_dir / "gsea"
    gsea_dir.mkdir(parents=True, exist_ok=True)

    # Extract control and test samples
    control_samples = metadata[metadata['condition'].str.lower() == 'cntl']['Sample name'].tolist()
    test_samples = metadata[metadata['condition'].str.lower() == 'test']['Sample name'].tolist()

    # Combine control and test samples
    all_samples = control_samples + test_samples

    # Prepare the content for the .cls file
    total_samples = len(all_samples)
    num_classes = 2
    class_pointers = 1  # Typically set to 1 unless multiple phenotype labels are present

    first_row = f"{total_samples} {num_classes} {class_pointers}"
    second_row = "# control test"  # GSEA expects a comment line starting with '#'
    third_row = " ".join(["0"] * len(control_samples) + ["1"] * len(test_samples))

    # Define the output file path within the GSEA directory
    gsea_cls_file = gsea_dir / "test.cls"

    # Write the content to the .cls file
    try:
        with open(gsea_cls_file, 'w') as f:
            f.write(f"{first_row}\n")
            f.write(f"{second_row}\n")
            f.write(f"{third_row}\n")
        print(f"GSEA class file saved to '{gsea_cls_file}'.")
    except Exception as e:
        print(f"Error saving GSEA class file '{gsea_cls_file}': {e}")
        sys.exit(1)

def run_gsea(results_dir, gsea_dir_name, species_name):
    """Run GSEA for each .gmt file in the referenceFiles/{species}/GSEA directory."""

    # Define paths
    gsea_dir = results_dir / gsea_dir_name
    gsea_dir.mkdir(parents=True, exist_ok=True)

    gsea_input_file = results_dir / "gsea" / "test.txt"
    gsea_cls_file = results_dir / "gsea" / "test.cls"
    gsea_results_dir = results_dir / "gsea"
    gsea_reference_dir = REFERENCE_DIR / species_name / "GSEA"
    gsea_cli_path = shutil.which("gsea-cli") or shutil.which("gsea-cli.sh")
    if not gsea_cli_path:
        raise FileNotFoundError("Neither gsea-cli nor gsea-cli.sh was found in PATH.")

    gmt_files = sorted(gsea_reference_dir.glob("*.gmt"))
    stems = {p.stem.lower() for p in gmt_files}
    # MSigDB ships GO both as the full C5 collection and as BP/MF/CC splits. The reference
    # tree has gene_ontology plus biological_process and molecular_function; the full set is
    # redundant, much larger, and is what OOM'd / crashed the smoke run after the splits had
    # already succeeded. Skip it when the splits are present.
    skip_full_go = (
        any(s.startswith("gene_ontology") for s in stems)
        and any(s.startswith("biological_process") for s in stems)
        and any(s.startswith("molecular_function") for s in stems)
    )

    # GSEA is a fat Java process; conda's default heap is too small for the largest GMTs.
    env = os.environ.copy()
    if "GSEA_JVM_HEAP" not in env and "JAVA_TOOL_OPTIONS" not in env:
        env["JAVA_TOOL_OPTIONS"] = "-Xmx16g"

    for gmt_file in gmt_files:
        if skip_full_go and gmt_file.stem.lower().startswith("gene_ontology"):
            print(f"Skipping {gmt_file.name}: covered by biological_process + "
                  f"molecular_function collections already in this reference set.")
            continue
        rpt_label = gmt_file.stem
        command = [
            gsea_cli_path, "GSEA",
            "-res", str(gsea_input_file),
            "-cls", str(gsea_cls_file),
            "-gmx", str(gmt_file),
            "-permute", "gene_set",
            "-collapse", "false",
            "-plot_top_x", "50",
            "-out", str(gsea_results_dir),
            "-rpt_label", rpt_label
        ]
        try:
            completed = subprocess.run(
                command, check=True, env=env, capture_output=True, text=True)
            if completed.stdout:
                print(completed.stdout, end="" if completed.stdout.endswith("\n") else "\n")
            print(f"GSEA analysis completed successfully for {gmt_file}")
        except subprocess.CalledProcessError as e:
            print(f"GSEA analysis failed for {gmt_file}: {e}")
            if e.stderr:
                print(e.stderr[-2000:])
            if e.stdout:
                print(e.stdout[-1000:])

def create_rnk_file(results_dir):
    """
    Create a .rnk file from DESeq2 results for GSEA pre-ranked analysis.

    Parameters:
    - results_dir (str or Path): Path to the results directory.
    """
    deseq2_results_file = results_dir / "deseq2" / "deseq2_results.csv"
    rnk_file = results_dir / "gsea" / "deseq2_ranked_genes.rnk"

    # Check if the DESeq2 results file exists
    if not deseq2_results_file.exists():
        print(f"Error: DESeq2 results file '{deseq2_results_file}' not found.")
        return None

    try:
        # Load the DESeq2 results with explicit column names
        column_names = ["gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
        deseq2_df = pd.read_csv(deseq2_results_file, sep=",", names=column_names, header=0)

        # Ensure required columns exist
        required_columns = {"gene", "log2FoldChange"}
        if not required_columns.issubset(deseq2_df.columns):
            print(f"Error: Missing required columns in DESeq2 results file: {required_columns - set(deseq2_df.columns)}")
            return None

        # Select 'gene' and 'log2FoldChange', drop any missing values
        rnk_df = deseq2_df[["gene", "log2FoldChange"]].dropna()

        # Sort by log2FoldChange (descending = upregulated at top)
        rnk_df = rnk_df.sort_values("log2FoldChange", ascending=False)

        # Save as tab-delimited .rnk file with NO header, NO index
        rnk_file.parent.mkdir(parents=True, exist_ok=True)
        rnk_df.to_csv(rnk_file, sep="\t", index=False, header=False)

        print(f"✅ .rnk file created successfully at '{rnk_file}'")
        return rnk_file
    except Exception as e:
        print(f"Error creating .rnk file: {e}")
        return None

def run_gsea_preranked(rnk_file, gsea_results_dir, gmt_file):
    """
    Run GSEA in pre-ranked mode using the .rnk file.

    Parameters:
    - rnk_file (str or Path): Path to the .rnk file.
    - gsea_results_dir (str or Path): Path to the GSEA results directory.
    - gmt_file (str or Path): Path to the .gmt file.
    """
    # Check if the .rnk file exists
    if not rnk_file.exists():
        raise FileNotFoundError(f".rnk file '{rnk_file}' not found.")

    # Define the output label for the pre-ranked analysis
    rpt_label = f"preranked_{gmt_file.stem}"

    # Ensure a GSEA CLI executable is available
    gsea_cli_path = shutil.which("gsea-cli") or shutil.which("gsea-cli.sh")
    if not gsea_cli_path:
        raise FileNotFoundError("Neither gsea-cli nor gsea-cli.sh was found in PATH. Please ensure GSEA is installed and accessible.")

    # Construct the GSEA pre-ranked command
    cmd = [
        gsea_cli_path, "GSEAPreranked",
        "-rnk", str(rnk_file),
        "-gmx", str(gmt_file),
        "-out", str(gsea_results_dir),
        "-rpt_label", rpt_label
    ]

    try:
        env = os.environ.copy()
        if "GSEA_JVM_HEAP" not in env and "JAVA_TOOL_OPTIONS" not in env:
            env["JAVA_TOOL_OPTIONS"] = "-Xmx16g"
        subprocess.run(cmd, check=True, env=env, capture_output=True, text=True)
        print(f"✅ GSEA pre-ranked analysis completed successfully for '{gmt_file}'")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running GSEA pre-ranked analysis for '{gmt_file}': {e}. Command: {' '.join(cmd)}")
        if e.stderr:
            print(e.stderr[-2000:])
        return False

def check_and_prepare_rmats(metadata, results_dir, sample_names, library_types):
    """Write the rMATS BAM lists.

    Replicate order in these files is the ``replicate_index`` recorded in
    splicing_event_replicate.tsv, so it follows metadata order rather than directory order.
    """
    distinct = set(library_types.values())
    if len(distinct) != 1:
        raise RuntimeError(
            f"samples do not share a library strandedness ({sorted(distinct)}); rMATS needs one "
            f"--libType for the whole comparison."
        )

    rmats_dir = results_dir / "rmats"
    rmats_dir.mkdir(parents=True, exist_ok=True)

    def ordered(condition_value):
        wanted = set(metadata[metadata['condition'].str.lower() == condition_value]['Sample name'])
        return [s for s in sample_names if s in wanted]

    def get_bam_paths(samples):
        bam_paths = []
        for sample in samples:
            bam_file = results_dir / "star" / sample / f"{sample}_Aligned.sortedByCoord.out.bam"
            if not bam_file.exists():
                raise FileNotFoundError(
                    f"BAM missing for {sample}: {bam_file}. rMATS replicate indices would no "
                    f"longer match metadata.tsv, so the run is stopping."
                )
            bam_paths.append(str(bam_file.resolve()))
        return bam_paths

    # --b1 is the test group; tables_rmats.GROUP_COLUMNS depends on that.
    test_bams_file = rmats_dir / "testBAMs.txt"
    test_bams_file.write_text(','.join(get_bam_paths(ordered('test'))))
    control_bams_file = rmats_dir / "controlBAMs.txt"
    control_bams_file.write_text(','.join(get_bam_paths(ordered('cntl'))))
    print(f"rMATS BAM lists written to '{rmats_dir}'.")

def run_rmats(metadata, results_dir, species_name, library_type):
    """Run rMATS analysis using generated BAM lists. Returns the arguments used."""

    # Define paths
    rmats_dir = results_dir / "rmats"
    test_bams_file = rmats_dir / "testBAMs.txt"
    control_bams_file = rmats_dir / "controlBAMs.txt"

    # Define GTF file path
    gtf_file = REFERENCE_DIR / species_name / "annotations" / f"{species_name}.gtf"
    if not gtf_file.exists():
        print(f"Error: GTF file '{gtf_file}' not found.")
        sys.exit(1)

    # Define output and temporary directories for rMATS
    rmats_output_dir = rmats_dir / "rmats_out"
    rmats_tmp_dir = rmats_dir / "rmats_tmp"
    rmats_output_dir.mkdir(parents=True, exist_ok=True)
    rmats_tmp_dir.mkdir(parents=True, exist_ok=True)

    # Determine if the dataset is paired-end or single-end
    if metadata['fastq_r2'].notnull().any():
        paired = True
        orientation = 'paired'
    else:
        paired = False
        orientation = 'single'
    print(f"Run type determined as '{orientation}'.")

    # Extract read lengths from fastp JSON files
    read_lengths = []
    all_samples = metadata['Sample name'].tolist()

    for sample in all_samples:
        json_file = results_dir / "fastp" / sample / f"{sample}.json"
        if not json_file.exists():
            print(f"Warning: JSON file '{json_file}' not found. Skipping read length extraction for this sample.")
            continue
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)
                read_length = data.get('read1_mean_length')
                if read_length:
                    read_lengths.append(int(read_length))
                else:
                    read_length = data.get('read1_before_filtering', {}).get('total_cycles')
                    if read_length:
                        read_lengths.append(int(read_length))
                        print(f"Info: Using 'total_cycles' as read length for sample '{sample}': {read_length}")
                    else:
                        print(f"Warning: Neither 'read1_mean_length' nor 'read1_before_filtering.total_cycles' found in '{json_file}'.")
        except Exception as e:
            print(f"Error reading '{json_file}': {e}")

    if not read_lengths:
        print("Error: No read lengths extracted from JSON files. Unable to set --readLength for rMATS.")
        sys.exit(1)

    if all(length == read_lengths[0] for length in read_lengths):
        read_length = read_lengths[0]
        print(f"All samples have the same read length: {read_length}.")
    else:
        read_length = int(sum(read_lengths) / len(read_lengths))
        print(f"Read lengths vary among samples. Using average read length: {read_length}.")

    lib_type = strandedness_module.rmats_lib_type(library_type)

    # Construct the rMATS command
    rmats_cmd = [
        "rmats.py",
        "--b1", str(test_bams_file),
        "--b2", str(control_bams_file),
        "--gtf", str(gtf_file),
        "-t", orientation,
        "--readLength", str(read_length),
        "--od", str(rmats_output_dir),
        "--tmp", str(rmats_tmp_dir),
        *RMATS_COMMON_ARGS,
        "--libType", lib_type
    ]

    print(f"Running rMATS with command: {' '.join(rmats_cmd)}")

    # Execute the rMATS command
    try:
        subprocess.run(rmats_cmd, check=True)
        print(f"rMATS analysis completed successfully. Results are in '{rmats_output_dir}'.")
    except subprocess.CalledProcessError as e:
        print(f"rMATS analysis failed: {e}")
        sys.exit(1)

    return " ".join(["-t", orientation, "--readLength", str(read_length),
                     *RMATS_COMMON_ARGS, "--libType", lib_type])

def generate_splicing_bar_chart(results_dir):
    """
    Generate bar charts of the splicing events using *.JCEC.txt files found in results_dir/rmats/rmats_out.

    Parameters:
    - results_dir (str or Path): The directory where the results are stored.
    """
    rmats_out_dir = Path(results_dir) / "rmats" / "rmats_out"

    # Load data
    a3ss = pd.read_csv(rmats_out_dir / "A3SS.MATS.JCEC.txt", sep='\t')
    a5ss = pd.read_csv(rmats_out_dir / "A5SS.MATS.JCEC.txt", sep='\t')
    mxe = pd.read_csv(rmats_out_dir / "MXE.MATS.JCEC.txt", sep='\t')
    ri = pd.read_csv(rmats_out_dir / "RI.MATS.JCEC.txt", sep='\t')
    se = pd.read_csv(rmats_out_dir / "SE.MATS.JCEC.txt", sep='\t')

    def plot_splicing_events(fdr_threshold, output_file):
        # Prepare dataframe with significant events
        a3ssSig = a3ss[a3ss['FDR'] <= fdr_threshold]
        a5ssSig = a5ss[a5ss['FDR'] <= fdr_threshold]
        mxeSig = mxe[mxe['FDR'] <= fdr_threshold]
        riSig = ri[ri['FDR'] <= fdr_threshold]
        seSig = se[se['FDR'] <= fdr_threshold]

        # Prepare dataframe with significant positive and negative events
        a3ssSigPos = a3ssSig[a3ssSig['IncLevelDifference'] > 0]
        a5ssSigPos = a5ssSig[a5ssSig['IncLevelDifference'] > 0]
        mxeSigPos = mxeSig[mxeSig['IncLevelDifference'] > 0]
        riSigPos = riSig[riSig['IncLevelDifference'] > 0]
        seSigPos = seSig[seSig['IncLevelDifference'] > 0]

        a3ssSigNeg = a3ssSig[a3ssSig['IncLevelDifference'] < 0]
        a5ssSigNeg = a5ssSig[a5ssSig['IncLevelDifference'] < 0]
        mxeSigNeg = mxeSig[mxeSig['IncLevelDifference'] < 0]
        riSigNeg = riSig[riSig['IncLevelDifference'] < 0]
        seSigNeg = seSig[seSig['IncLevelDifference'] < 0]

        # For positive events
        yp1 = len(a3ssSigPos)
        yp2 = len(a5ssSigPos)
        yp3 = len(mxeSigPos)
        yp4 = len(riSigPos)
        yp5 = len(seSigPos)

        yn1 = len(a3ssSigNeg)
        yn2 = len(a5ssSigNeg)
        yn3 = len(mxeSigNeg)
        yn4 = len(riSigNeg)
        yn5 = len(seSigNeg)

        # Prepare dataframe with number of altered events from all splice types
        x = ['A3SS', 'A5SS', 'MXE', 'RI', 'SE']
        y1 = [yp1, yp2, yp3, yp4, yp5]
        y2 = [yn1, yn2, yn3, yn4, yn5]

        # Plot bars in stacked manner
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.sans-serif'] = ['Arial']
        frameon = False
        plt.bar(x, y1, color='c', label='positive', edgecolor='k', linewidth=0.25)
        plt.bar(x, y2, bottom=y1, color='salmon', label='negative', edgecolor='k', linewidth=0.25)
        ax = plt.gca()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.ylabel("Number of altered splicing events")
        plt.legend(loc='upper left')
        plt.title(f"rMATS Alternative Splicing Events FDR {fdr_threshold}")
        plt.savefig(output_file, dpi=300)
        plt.close()

    # Generate bar charts for FDR 0.05 and FDR 0.0005
    plot_splicing_events(0.05, rmats_out_dir / "barChart_rmatsASE_fdr0.05.png")
    plot_splicing_events(0.0005, rmats_out_dir / "barChart_rmatsASE_fdr0.0005.png")

def collate_tpms(results_dir, sample_names, species_name):
    """
    Collate the TPMs for genes for each sample and write them to a TSV file.

    Parameters:
    - results_dir (str or Path): Path to the results directory.
    - sample_names (list of str): List of sample names.
    - species_name (str): Name of the species.
    """
    results_dir = Path(results_dir)
    collated_data = {}

    for sample in sample_names:
        sample_path = results_dir / "kallisto" / sample / f"{sample}_{species_name}" / "abundance_gene_tpm.tsv"
        if not sample_path.exists():
            print(f"Error: File '{sample_path}' not found.")
            continue

        print(f"Reading {sample_path}...")
        df = pd.read_table(sample_path, index_col=0, sep='\t')
        if 'gene' not in collated_data:
            collated_data['gene'] = df.index

        collated_data[sample] = df['tpm']

    # Create a DataFrame with the collated data
    collated_df = pd.DataFrame(collated_data)
    collated_df.index.name = 'gene'

    # Drop the duplicate 'gene' column if it exists
    if 'gene' in collated_df.columns:
        collated_df = collated_df.drop(columns=['gene'])

    output_path = results_dir / "deseq2" / "deseq2_collated_tpms.tsv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    collated_df.to_csv(output_path, sep='\t', index=True)
    print(f"deseq2_collated_tpms.tsv file created at '{output_path}'.")

def merge_deseq2_results(results_dir):
    """
    Merge the deseq2_collated_tpms.tsv file with the deseq2_results.csv file.

    Parameters:
    - results_dir (str or Path): Path to the results directory.
    """
    results_dir = Path(results_dir)
    collated_tpms_path = results_dir / "deseq2" / "deseq2_collated_tpms.tsv"
    deseq2_results_path = results_dir / "deseq2" / "deseq2_results.csv"
    metadata_path = results_dir / "deseq2" / "deseq2_metadata.tsv"
    output_path = results_dir / "deseq2" / "deseq2_results_genes.tsv"

    # Read the collated TPMs and DESeq2 results files
    collated_tpms_df = pd.read_csv(collated_tpms_path, sep='\t', index_col=0)
    deseq2_results_df = pd.read_csv(deseq2_results_path, sep=',', index_col=0)

    # Read the metadata file to determine control and test samples
    metadata_df = pd.read_csv(metadata_path, sep='\t')
    # deseq2_metadata.tsv labels controls 'cntl' (only the sleuth metadata is rewritten to
    # 'control'); matching 'control' alone silently dropped every control TPM column.
    control_samples = metadata_df[metadata_df['condition'].isin(['cntl', 'control'])]['Sample name'].tolist()
    test_samples = metadata_df[metadata_df['condition'] == 'test']['Sample name'].tolist()

    # Merge the dataframes on the gene index
    merged_df = collated_tpms_df.join(deseq2_results_df, how='inner')

    # Ensure all control and test samples are included in the column order
    sample_columns = control_samples + test_samples
    deseq2_columns = ['baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj']
    column_order = [col for col in sample_columns if col in merged_df.columns] + deseq2_columns

    # Reorder the columns
    merged_df = merged_df[column_order]

    # Append ' (TPMs)' to sample columns
    new_columns = []
    for col in merged_df.columns:
        if col in sample_columns:
            new_columns.append(f"{col} (TPMs)")
        else:
            new_columns.append(col)
    merged_df.columns = new_columns

    # Save the merged dataframe to a TSV file
    output_path.parent.mkdir(parents=True, exist_ok=True)
    merged_df.to_csv(output_path, sep='\t', index=True)
    print(f"deseq2_results_genes.tsv file created at '{output_path}'.")

def reform_deseq2_results(results_dir):
    """
    Reform the deseq2_results_genes.tsv file to an Excel file with specified formatting.

    Parameters:
    - results_dir (str or Path): Path to the results directory.
    """
    script_path = scripts / "reform_deseq2.py"
    try:
        subprocess.run(["python", str(script_path), str(results_dir)], check=True)
        print("Reformed DESeq2 results to Excel format.")
    except subprocess.CalledProcessError as e:
        print(f"Error reforming DESeq2 results: {e}")
        sys.exit(1)

def collate_transcript_tpms(results_dir, sample_names, species_name):
    """
    Collate the TPMs for transcripts for each sample and write them to a TSV file.

    Parameters:
    - results_dir (str or Path): Path to the results directory.
    - sample_names (list of str): List of sample names.
    - species_name (str): Name of the species.
    """
    results_dir = Path(results_dir)
    collated_data = {}

    for sample in sample_names:
        sample_path = results_dir / "kallisto" / sample / f"{sample}_{species_name}" / "abundance.tsv"
        if not sample_path.exists():
            print(f"Error: File '{sample_path}' not found.")
            continue

        print(f"Reading {sample_path}...")
        df = pd.read_table(sample_path, index_col=0, sep='\t')
        if 'target_id' not in collated_data:
            collated_data['target_id'] = df.index.to_series()

        collated_data[sample] = df['tpm']

    # Create a DataFrame with the collated data
    collated_df = pd.DataFrame(collated_data)
    collated_df.index.name = 'target_id'

    # Ensure the 'target_id' column is not duplicated
    if 'target_id' not in collated_df.columns:
        collated_df.insert(0, 'target_id', collated_data['target_id'])

    output_path = results_dir / "sleuth" / "sleuth_collated_transcripts_tpms.tsv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    collated_df.to_csv(output_path, sep='\t', index=False)
    print(f"sleuth_collated_transcripts_tpms.tsv file created at '{output_path}'.")

def merge_sleuth_wald_results(results_dir):
    """
    Merge the sleuth_collated_transcripts_tpms.tsv file with the sleuth_DETranscripts_wt.csv file.

    Parameters:
    - results_dir (str or Path): Path to the results directory.
    """
    sleuth_tpms_path = results_dir / "sleuth" / "sleuth_collated_transcripts_tpms.tsv"
    sleuth_results_path = results_dir / "sleuth" / "sleuth_DETranscripts_wt.csv"
    output_path = results_dir / "sleuth" / "sleuth_results_wald_transcripts.tsv"

    # Read the sleuth collated TPMs and sleuth results files
    sleuth_tpms_df = pd.read_csv(sleuth_tpms_path, sep='\t')
    sleuth_results_df = pd.read_csv(sleuth_results_path, sep=',')

    # Ensure the 'target_id' column is present in both dataframes
    if 'target_id' not in sleuth_tpms_df.columns or 'target_id' not in sleuth_results_df.columns:
        print("Error: 'target_id' column not found in one of the input files.")
        return

    # Drop the duplicated 'target_id' column if it exists in sleuth_tpms_df
    if 'target_id.1' in sleuth_tpms_df.columns:
        sleuth_tpms_df = sleuth_tpms_df.drop(columns=['target_id.1'])

    # Merge the dataframes on the target_id column
    merged_df = sleuth_tpms_df.merge(sleuth_results_df[['target_id', 'gene', 'b', 'qval']], on='target_id', how='inner')

    # Reorder columns as specified
    ordered_columns = ['target_id'] + sleuth_tpms_df.columns[1:].tolist() + ['gene', 'b', 'qval']
    merged_df = merged_df[ordered_columns]

    # Save the merged dataframe to a TSV file
    merged_df.to_csv(output_path, sep='\t', index=False)
    print(f"sleuth_results_wald_transcripts.tsv file created at '{output_path}'.")

def merge_sleuth_lrt_results(results_dir):
    """
    Merge the sleuth_collated_transcripts_tpms.tsv file with the sleuth_DETranscripts_lrt.csv file.

    Parameters:
    - results_dir (str or Path): Path to the results directory.
    """
    sleuth_tpms_path = results_dir / "sleuth" / "sleuth_collated_transcripts_tpms.tsv"
    sleuth_results_path = results_dir / "sleuth" / "sleuth_DETranscripts_lrt.csv"
    output_path = results_dir / "sleuth" / "sleuth_results_lrt_transcripts.tsv"

    # Read the sleuth collated TPMs and sleuth results files
    sleuth_tpms_df = pd.read_csv(sleuth_tpms_path, sep='\t')
    sleuth_results_df = pd.read_csv(sleuth_results_path, sep=',')

    # Ensure the 'target_id' column is present in both dataframes
    if 'target_id' not in sleuth_tpms_df.columns or 'target_id' not in sleuth_results_df.columns:
        print("Error: 'target_id' column not found in one of the input files.")
        return

    # Drop the duplicated 'target_id' column if it exists in sleuth_tpms_df
    if 'target_id.1' in sleuth_tpms_df.columns:
        sleuth_tpms_df = sleuth_tpms_df.drop(columns=['target_id.1'])

    # Merge the dataframes on the target_id column
    merged_df = sleuth_tpms_df.merge(sleuth_results_df[['target_id', 'gene', 'qval']], on='target_id', how='inner')

    # Reorder columns as specified
    ordered_columns = ['target_id'] + sleuth_tpms_df.columns[1:].tolist() + ['gene', 'qval']
    merged_df = merged_df[ordered_columns]

    # Save the merged dataframe to a TSV file
    merged_df.to_csv(output_path, sep='\t', index=False)
    print(f"sleuth_results_lrt_transcripts.tsv file created at '{output_path}'.")

def reformat_sleuth_results(results_dir):
    """
    Reformat the sleuth results files to Excel files with specified formatting.

    Parameters:
    - results_dir (str or Path): Path to the results directory.
    """
    script_path = scripts / "reform_sleuth.py"
    try:
        subprocess.run(["python", str(script_path), str(results_dir)], check=True)
        print("Reformed Sleuth results to Excel format.")
    except subprocess.CalledProcessError as e:
        print(f"Error reforming Sleuth results: {e}")
        sys.exit(1)

def check_mycoplasma_sequences(results_dir):
    """
    Check if samples are positive or negative for specific mycoplasma sequences using mycoplasma_check.py.

    Parameters:
    - results_dir (str or Path): Path to the results directory.
    """
    script_path = scripts / "mycoplasma_check.py"
    try:
        subprocess.run(["python", str(script_path), str(results_dir)], check=True)
        print("Mycoplasma sequences check completed.")
    except subprocess.CalledProcessError as e:
        print(f"Error checking mycoplasma sequences: {e}")
        sys.exit(1)

def prepare_volcano(results_dir):
    """
    Prepare volcano plots for DESeq2 and Sleuth results.

    Parameters:
    - results_dir (str or Path): Path to the results directory.
    """
    script_path = scripts / "volcano.py"
    try:
        subprocess.run(["python", str(script_path), str(results_dir)], check=True)
        print("Volcano plots prepared.")
    except subprocess.CalledProcessError as e:
        print(f"Error preparing volcano plots: {e}")
        sys.exit(1)

def convert_SJtab_to_bed(sample_name, results_dir):
    """Produce the IGV junction BEDs.

    junction.tsv is the canonical form of this data (Task 11); these files are kept only
    because they are what gets dragged into IGV, and they live under artifacts/ accordingly.
    """
    results_dir = Path(results_dir)
    star_dir = results_dir / "star" / sample_name
    sjtab_files = list(star_dir.glob("*.tab"))

    for sjtab_file in sjtab_files:
        cmd = f"perl {scripts / 'splice_junction_strand_separator.pl'} {sjtab_file}"
        subprocess.run(cmd, shell=True, check=True)
        print(f"Converted {sjtab_file} to BED format")

def remove_wig(results_dir):
    """
    Remove all files ending in '.wig' from results_dir/star/sampleName.

    Parameters:
    - results_dir (str or Path): The directory where the results are stored.
    """
    results_dir = Path(results_dir)
    star_dir = results_dir / "star"

    # Iterate over each sample directory in star_dir
    for sample_dir in star_dir.iterdir():
        if sample_dir.is_dir():
            for wig_file in sample_dir.glob("*.wig"):
                os.remove(wig_file)
                print(f"Removed file: {wig_file}")

def report_top_genes_and_transcripts(results_dir):
    """
    Report the top 10 genes with the lowest padj value and the top 10 transcripts with the lowest q-value.

    Parameters:
    - results_dir (str or Path): Path to the results directory.
    """
    results_dir = Path(results_dir)

    # Report top 10 genes with the lowest padj value
    deseq2_results_genes_path = results_dir / "deseq2" / "deseq2_results_genes.tsv"
    if deseq2_results_genes_path.exists():
        deseq2_results_genes_df = pd.read_csv(deseq2_results_genes_path, sep='\t')
        top_genes_df = deseq2_results_genes_df.nsmallest(10, 'padj')
        top_genes_path = results_dir / "topGenes.tsv"
        top_genes_df.to_csv(top_genes_path, sep='\t', index=False)
        print(f"Top 10 genes with the lowest padj value saved to {top_genes_path}")
    else:
        print(f"File {deseq2_results_genes_path} not found.")

    # Report top 10 transcripts with the lowest q-value
    sleuth_results_wald_transcripts_path = results_dir / "sleuth" / "sleuth_results_wald_transcripts.tsv"
    if sleuth_results_wald_transcripts_path.exists():
        sleuth_results_wald_transcripts_df = pd.read_csv(sleuth_results_wald_transcripts_path, sep='\t')
        top_transcripts_df = sleuth_results_wald_transcripts_df.nsmallest(10, 'qval')
        top_transcripts_path = results_dir / "topTranscripts.tsv"
        top_transcripts_df.to_csv(top_transcripts_path, sep='\t', index=False)
        print(f"Top 10 transcripts with the lowest q-value saved to {top_transcripts_path}")
    else:
        print(f"File {sleuth_results_wald_transcripts_path} not found.")

def report_top_gene_tpms(results_dir):
    """
    Inspect results_dir/deseq2/deseq2_collated_tpms.tsv and generate a tab-delimited file in results_dir
    called top_gene_tpms.tsv which for each sample lists the genes with the three highest TPM values.

    Parameters:
    - results_dir (str or Path): Path to the results directory.
    """
    results_dir = Path(results_dir)
    collated_tpms_path = results_dir / "deseq2" / "deseq2_collated_tpms.tsv"
    output_path = results_dir / "top_gene_tpms.tsv"

    if not collated_tpms_path.exists():
        print(f"File {collated_tpms_path} not found.")
        return

    collated_tpms_df = pd.read_csv(collated_tpms_path, sep='\t', index_col=0)

    top_genes_list = []
    for sample in collated_tpms_df.columns:
        top_genes = collated_tpms_df[sample].nlargest(3).reset_index()
        top_genes.columns = ['gene', f'{sample}_tpm']
        top_genes_list.append(top_genes)

    top_genes_df = pd.concat(top_genes_list, axis=1)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    top_genes_df.to_csv(output_path, sep='\t', index=False)
    print(f"Top genes with the highest TPM values saved to {output_path}")

def summarize_alignments(results_dir, sample_names):
    summary_file = results_dir / "alignmentSummary.tsv"
    summary_data = []

    # Define the headers for the summary file
    headers = [
        "SampleName",
        "Number of input reads",
        "Average input read length",
        "Average mapped length",
        "Uniquely mapped reads number",
        "Number of reads mapped to multiple loci",
        "Number of reads mapped to too many loci",
        "Mapped percentage"
    ]
    summary_data.append(headers)

    # Iterate over sample_names in caller-provided order
    for sample_name in sample_names:
        sample_dir = results_dir / "star" / sample_name
        log_file = sample_dir / f"{sample_name}_Log.final.out"

        if log_file.exists():
            with open(log_file, 'r') as f:
                log_content = f.read()

            # Extract relevant information using regex
            num_input_reads = re.search(r"Number of input reads\s+\|\s+(\d+)", log_content).group(1)
            avg_input_read_length = re.search(r"Average input read length\s+\|\s+(\d+)", log_content).group(1)
            avg_mapped_length = re.search(r"Average mapped length\s+\|\s+(\d+)", log_content).group(1)
            uniquely_mapped_reads_num = re.search(r"Uniquely mapped reads number\s+\|\s+(\d+)", log_content).group(1)
            reads_mapped_multiple_loci = re.search(r"Number of reads mapped to multiple loci\s+\|\s+(\d+)", log_content).group(1)
            reads_mapped_too_many_loci = re.search(r"Number of reads mapped to too many loci\s+\|\s+(\d+)", log_content).group(1)
            uniquely_mapped_percentage = float(re.search(r"Uniquely mapped reads %\s+\|\s+([\d\.]+)", log_content).group(1))
            multiple_mapped_percentage = float(re.search(r"% of reads mapped to multiple loci\s+\|\s+([\d\.]+)", log_content).group(1))
            mapped_percentage = uniquely_mapped_percentage + multiple_mapped_percentage

            # Append the extracted data to the summary_data list
            summary_data.append([
                sample_name,
                num_input_reads,
                avg_input_read_length,
                avg_mapped_length,
                uniquely_mapped_reads_num,
                reads_mapped_multiple_loci,
                reads_mapped_too_many_loci,
                f"{mapped_percentage:.2f}"
            ])

    # Write the summary data to the alignmentSummary.tsv file
    with open(summary_file, 'w') as f:
        for row in summary_data:
            f.write("\t".join(row) + "\n")

def summarize_transcript_coverage(results_dir, sample_names, species_name):
    output_file = results_dir / "transcriptCoverage.tsv"

    with open(output_file, 'w') as out_f:
        out_f.write("SampleName\tn_processed\tn_pseudoaligned\tp_pseudoaligned\n")

        for sample_name in sample_names:
            run_info_file = results_dir / f"kallisto/{sample_name}/{sample_name}_{species_name}/run_info.json"

            if run_info_file.exists():
                with open(run_info_file, 'r') as f:
                    run_info = json.load(f)
                    n_processed = run_info.get("n_processed", "NA")
                    n_pseudoaligned = run_info.get("n_pseudoaligned", "NA")
                    p_pseudoaligned = run_info.get("p_pseudoaligned", "NA")

                    out_f.write(f"{sample_name}\t{n_processed}\t{n_pseudoaligned}\t{p_pseudoaligned}\n")
            else:
                print(f"Warning: {run_info_file} does not exist")

def generate_results_html(results_dir, sample_names):
    results_file = results_dir / "reports" / "report.html"
    results_file.parent.mkdir(parents=True, exist_ok=True)

    # Paths to the required files
    mycoplasma_report_file = results_dir / "mycoplasma_report.tsv"
    top_genes_file = results_dir / "topGenes.tsv"
    alignment_summary_file = results_dir / "alignmentSummary.tsv"
    transcript_coverage_file = results_dir / "transcriptCoverage.tsv"

    # Read the mycoplasma report
    mycoplasma_report = []
    if mycoplasma_report_file.exists():
        with open(mycoplasma_report_file, 'r') as f:
            for line in f:
                mycoplasma_report.append(line.strip().split('\t'))

    # Read the top genes report
    top_genes_report = []
    if top_genes_file.exists():
        with open(top_genes_file, 'r') as f:
            for line in f:
                top_genes_report.append(line.strip().split('\t'))

    # Read the strandedness report
    strandedness_report = []
    for sample_name in sample_names:
        strandedness_file = strandedness_module.rseqc_path(results_dir, sample_name)
        if strandedness_file.exists():
            parsed = strandedness_module.parse_rseqc(strandedness_file.read_text())
            strandedness_report.append(parsed["library_type"])

    # Read the alignment summary
    alignment_summary = []
    if alignment_summary_file.exists():
        with open(alignment_summary_file, 'r') as f:
            for line in f:
                alignment_summary.append(line.strip().split('\t'))

    # Read the transcript coverage report
    transcript_coverage = []
    if transcript_coverage_file.exists():
        with open(transcript_coverage_file, 'r') as f:
            for line in f:
                transcript_coverage.append(line.strip().split('\t'))

    html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <title>RNA-Seq Results Summary</title>
        <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css" integrity="sha384-Gn5384xqQ1aoWXA+058RXPxPg6fy4IWvTNh0E263XmFcJlSAwiGgFAW/dAiS6JXm" crossorigin="anonymous">
        <script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/js/bootstrap.min.js" integrity="sha384-JZR6Spejh4U02d8jOt6vLEHfe/JQGiRRSQQxSfFWpi1MquVdAyjUar5+76PVCmYl" crossorigin="anonymous"></script>
        <style>
            body {{
                font-size: 12px; /* Decrease the font size */
                text-align: left; /* Left-align the text */
            }}
            .table th, .table td {{
                font-size: 12px; /* Decrease the font size for table cells */
                text-align: left; /* Left-align the table cells */
            }}
        </style>
    </head>
    <body>
        <div class="container">
            <div class="jumbotron p-3 p-md-5">
                <h1 class="jumbotron-heading">RNA-Seq Results Summary</h1>
            </div>
            <div class="container">
                <h2>Strandedness</h2>
                <table class="table">
                    <thead class="thead-dark">
                        <tr>
                            {"".join([f"<th>{sample_name}</th>" for sample_name in sample_names])}
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            {"".join([f"<td>{strandedness}</td>" for strandedness in strandedness_report])}
                        </tr>
                    </tbody>
                </table>
            </div>
            <div class="container">
                <h2>Mapping Statistics (from Kallisto for Abundance Measurement)</h2>
                <table class="table">
                    <thead class="thead-dark">
                        <tr>
                            {"".join([f"<th>{col}</th>" for col in transcript_coverage[0]])}
                        </tr>
                    </thead>
                    <tbody>
                        {"".join([f"<tr>{''.join([f'<td>{cell}</td>' for cell in row])}</tr>" for row in transcript_coverage[1:]])}
                    </tbody>
                </table>
            </div>
            <div class="container">
                <h2>Mapping Statistics (from STAR for Splicing Analysis)</h2>
                <table class="table">
                    <thead class="thead-dark">
                        <tr>
                            {"".join([f"<th>{col}</th>" for col in alignment_summary[0]])}
                        </tr>
                    </thead>
                    <tbody>
                        {"".join([f"<tr>{''.join([f'<td>{cell}</td>' for cell in row])}</tr>" for row in alignment_summary[1:]])}
                    </tbody>
                </table>
            </div>
            <div class="container">
                <h2>Mycoplasma Report</h2>
                <table class="table">
                    <thead class="thead-dark">
                        <tr>
                            {"".join([f"<th>{col}</th>" for col in mycoplasma_report[0]])}
                        </tr>
                    </thead>
                    <tbody>
                        {"".join([f"<tr>{''.join([f'<td>{cell}</td>' for cell in row])}</tr>" for row in mycoplasma_report[1:]])}
                    </tbody>
                </table>
            </div>
            <div class="container">
                <h2>Top DE Genes</h2>
                <table class="table">
                    <thead class="thead-dark">
                        <tr>
                            {"".join([f"<th>{col}</th>" for col in top_genes_report[0]])}
                        </tr>
                    </thead>
                    <tbody>
                        {"".join([f"<tr>{''.join([f'<td>{cell}</td>' for cell in row])}</tr>" for row in top_genes_report[1:]])}
                    </tbody>
                </table>
            </div>
        </div>
    </body>
    <script src="https://code.jquery.com/jquery-3.2.1.slim.min.js" integrity="sha384-KJ3o2DKtIkvYIK3UENzmM7KCkRr/rE9/Qpg6aAZGJwFDMVNA/GpGFF93hXpG5KkN" crossorigin="anonymous"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.12.9/umd/popper.min.js" integrity="sha384-ApNbgh9B+Y1QKtv3Rn7W3mgPxhU9K/ScQsAP7hUibX39j7fakFPskvXusvfa0b4Q" crossorigin="anonymous"></script>
    </html>
    """

    with open(results_file, 'w') as f:
        f.write(html_content)

def organize_artifacts(results_dir, sample_names):
    """Move raw tool output under artifacts/ and reports under reports/.

    Replaces the old 1_Preprocessing/2_Analysis/3_Results tree and the four post-hoc renaming
    passes. Filenames here are cosmetic: tables/ is the contract surface and the manifests are
    authoritative, so nothing downstream parses these paths. Raw rMATS and STAR output is moved,
    never deleted.
    """
    results_dir = Path(results_dir)
    artifacts = results_dir / "artifacts"

    def move_into(source, destination):
        source = Path(source)
        if not source.exists():
            return
        destination.mkdir(parents=True, exist_ok=True)
        target = destination / source.name
        if target.exists():
            shutil.rmtree(target) if target.is_dir() else target.unlink()
        shutil.move(str(source), str(target))

    def move_glob(pattern, destination, root=None):
        destination.mkdir(parents=True, exist_ok=True)
        for path in sorted((root or results_dir).glob(pattern)):
            move_into(path, destination)

    for stage in ("fastp", "fastqc", "rseqc"):
        move_into(results_dir / stage, artifacts / "qc")

    for sample_name in sample_names:
        star_dir = results_dir / "star" / sample_name
        if not star_dir.is_dir():
            continue
        move_glob(f"{sample_name}*.ba*", artifacts / "bam", root=star_dir)
        move_glob(f"{sample_name}*.bed", artifacts / "junction_bed", root=star_dir)
        move_glob(f"{sample_name}*.tab", artifacts / "star_logs", root=star_dir)
        move_glob(f"{sample_name}*.junction", artifacts / "star_logs", root=star_dir)
        move_glob(f"{sample_name}*.out", artifacts / "star_logs", root=star_dir)

    # Figures before the bulk moves, so they land in reports/ rather than artifacts/.
    for pattern in ("deseq2/*.png", "rmats/rmats_out/*.png"):
        move_glob(pattern, results_dir / "reports")

    move_into(results_dir / "kallisto", artifacts / "kallisto")
    for pattern, destination in (
        ("deseq2/*", artifacts / "deseq2"),
        ("sleuth/*", artifacts / "sleuth"),
        ("gsea/*", artifacts / "gsea"),
        ("rmats/rmats_out/*", artifacts / "rmats_raw"),
        ("rmats/rmats_tmp", artifacts / "rmats_raw"),
        ("rmats/*.txt", artifacts / "rmats_raw"),
    ):
        move_glob(pattern, destination)

    for name in ("topGenes.tsv", "topTranscripts.tsv", "top_gene_tpms.tsv",
                 "alignmentSummary.tsv", "transcriptCoverage.tsv", "mycoplasma_report.tsv",
                 "Rplots.pdf"):
        move_into(results_dir / name, results_dir / "reports")
    move_glob("*.png", results_dir / "reports")

    # Everything worth keeping has been moved; what is left is scratch.
    for scratch in ("star", "deseq2", "sleuth", "gsea", "rmats"):
        path = results_dir / scratch
        if path.is_dir():
            shutil.rmtree(path)

    # macOS Finder writes these onto network/local mounts while the run is in progress.
    # The validator rejects them, so scrub before checksums/validation rather than failing
    # a good analysis over desktop detritus.
    for junk in results_dir.rglob(".DS_Store"):
        junk.unlink(missing_ok=True)

def process_sample(sample_name, fastq1, fastq2):
    import datetime

    def log_processing_error(sample_name, step, error_msg):
        timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        with open(stage_log_path(results_dir, "errors"), 'a') as logf:
            logf.write(f"{timestamp}\t{sample_name}\t{step}\t{error_msg}\n")

    def check_fastq_integrity(fastq_path, max_records=100000):
        import gzip
        open_func = gzip.open if fastq_path.endswith('.gz') else open
        try:
            with open_func(fastq_path, 'rt') as f:
                line_num = 0
                while True:
                    header = f.readline()
                    seq = f.readline()
                    plus = f.readline()
                    qual = f.readline()
                    if not qual:
                        break
                    line_num += 4
                    if len(seq.strip()) != len(qual.strip()):
                        return f"ERROR: sequence and quality have different length in {fastq_path} at record starting with {header.strip()}"
                    if max_records and line_num >= max_records * 4:
                        break
        except Exception as e:
            return f"ERROR reading {fastq_path}: {e}"
        return None  # No error
    """Process a single sample from QC to Kallisto, with robust output validation."""
    import json
    results_path = Path(results_dir)
    fastp_results_dir = results_path / 'fastp' / sample_name
    fastqc_results_dir = results_path / 'fastqc' / sample_name
    star_results_dir = results_path / 'star' / sample_name
    genome_dir = REFERENCE_DIR / species_name / "STAR"
    gtf_file = REFERENCE_DIR / species_name / "annotations" / f"{species_name}.gtf"

    def file_nonempty(path):
        return path.exists() and path.stat().st_size > 0

    def fastp_valid():
        html = fastp_results_dir / f"{sample_name}.html"
        jsonf = fastp_results_dir / f"{sample_name}.json"
        return file_nonempty(html) and file_nonempty(jsonf)

    def fastqc_valid():
        # Dynamically determine expected FastQC output files based on FASTQ filenames
        import os
        fastq_files = []
        if isinstance(fastq1, str) and fastq1 and fastq1.lower() != 'nan':
            fastq_files.append(fastq1)
        if isinstance(fastq2, str) and fastq2 and fastq2.lower() != 'nan':
            fastq_files.append(fastq2)
        files = []
        for fq in fastq_files:
            fq_base = os.path.basename(fq)
            # Remove .gz if present
            if fq_base.endswith('.gz'):
                fq_base = fq_base[:-3]
            # Remove .fastq or .fq
            if fq_base.endswith('.fastq'):
                fq_base = fq_base[:-6]
            elif fq_base.endswith('.fq'):
                fq_base = fq_base[:-3]
            # FastQC output is <basename>_fastqc.html and <basename>_fastqc.zip
            files.append(fastqc_results_dir / f"{fq_base}_fastqc.html")
            files.append(fastqc_results_dir / f"{fq_base}_fastqc.zip")
        return all(file_nonempty(f) for f in files)

    def star_valid():
        # STAR: sorted BAM, BAM index, Log.final.out, SJ.out.tab
        bam = star_results_dir / f"{sample_name}_Aligned.sortedByCoord.out.bam"
        bai = star_results_dir / f"{sample_name}_Aligned.sortedByCoord.out.bam.bai"
        log = star_results_dir / f"{sample_name}_Log.final.out"
        sjtab = star_results_dir / f"{sample_name}_SJ.out.tab"
        return all(file_nonempty(f) for f in [bam, bai, log, sjtab])

    def kallisto_valid():
        # Kallisto: abundance_gene_tpm.tsv, run_info.json
        kallisto_dir = results_path / 'kallisto' / sample_name / f"{sample_name}_{species_name}"
        tpm = kallisto_dir / "abundance_gene_tpm.tsv"
        runinfo = kallisto_dir / "run_info.json"
        return file_nonempty(tpm) and file_nonempty(runinfo)

    # If all outputs are present and non-empty, print and return
    if fastp_valid() and fastqc_valid() and star_valid() and kallisto_valid():
        print(f"Sample {sample_name} has already been processed and all outputs are valid. Skipping.")
        return strandedness_module.read_library_type(results_dir, sample_name)

    print(f"Processing sample: {sample_name}")

    # Run fastp without trimming
    try:
        import subprocess, tempfile
        # Ensure fastp output directory exists
        fastp_results_dir.mkdir(parents=True, exist_ok=True)
        print(f"DEBUG: Created fastp output directory: {fastp_results_dir}")
        def run_fastp_no_trimming_checked(fastq1, fastq2, html_out, json_out, terminal_log_path):
            cmd = [
                "fastp", "-w", "10", "-i", fastq1
            ]
            # Only add -I if fastq2 is a valid string (not None, not empty, not nan/float)
            if fastq2 and isinstance(fastq2, str) and fastq2.lower() != 'nan':
                cmd += ["-I", fastq2]
            cmd += [
                "--detect_adapter_for_pe",
                "--html", str(html_out),
                "--json", str(json_out),
                "--disable_trim_poly_g",
                "--disable_length_filtering"
            ]
            print(f"DEBUG: Running fastp with command: {' '.join(str(x) for x in cmd)}")
            with open(terminal_log_path, 'w') as termlog:
                result = subprocess.run(cmd, stdout=termlog, stderr=subprocess.STDOUT, text=True)
            print(f"DEBUG: fastp finished with return code {result.returncode}")
            # Now check the log for errors
            with open(terminal_log_path, 'r') as termlog:
                log_content = termlog.read()
            if "ERROR: sequence and quality have different length" in log_content:
                raise RuntimeError(log_content)
            if result.returncode != 0:
                raise RuntimeError(f"fastp failed: {log_content}")
        terminal_log_path = stage_log_path(results_dir, f"fastp_{sample_name}")
        run_fastp_no_trimming_checked(
            fastq1, fastq2,
            html_out=fastp_results_dir / f"{sample_name}.html",
            json_out=fastp_results_dir / f"{sample_name}.json",
            terminal_log_path=terminal_log_path
        )
        print(f"DEBUG: fastp output should be at {fastp_results_dir / f'{sample_name}.html'} and {fastp_results_dir / f'{sample_name}.json'}")
    except Exception as e:
        log_processing_error(sample_name, "fastp", str(e))
        raise

    # Run FastQC
    try:
        run_fastqc(fastq1, fastq2, fastqc_results_dir)
    except Exception as e:
        log_processing_error(sample_name, "FastQC", str(e))
        raise

    # FASTQ integrity check before STAR
    # Only check integrity for non-empty fastq files
    # Only treat as paired-end if fastq2 is a non-empty string and not 'nan'
    if isinstance(fastq2, str) and fastq2 and fastq2.lower() != 'nan':
        fq_list = [fastq1, fastq2]
        label_list = ["read1", "read2"]
        is_paired_end = True
    else:
        fq_list = [fastq1]
        label_list = ["read1"]
        is_paired_end = False
    for fq, read_label in zip(fq_list, label_list):
        # Only check if fq is a valid string and not nan/float
        if isinstance(fq, str) and fq.lower() != 'nan':
            err = check_fastq_integrity(fq)
            if err:
                log_processing_error(sample_name, f"FASTQ integrity ({read_label})", err)
                raise RuntimeError(err)

    # Run STAR alignment (set is_paired_end correctly)
    try:
        run_star(sample_name, fastq1, fastq2 if is_paired_end else None, star_results_dir, species_name, is_paired_end=is_paired_end)
    except Exception as e:
        log_processing_error(sample_name, "STAR", str(e))
        raise

    # Post-process alignment
    try:
        sorted_bam = post_process_alignment(sample_name, star_results_dir)
    except Exception as e:
        log_processing_error(sample_name, "post_process_alignment", str(e))
        raise

    # Infer strandedness
    try:
        strandedness = infer_strandedness(sorted_bam, species_name, results_dir)
    except Exception as e:
        log_processing_error(sample_name, "infer_strandedness", str(e))
        raise

    # Run Kallisto quantification for species
    try:
        run_kallisto(fastq1, fastq2, sample_name, species_name, strandedness, results_dir)
    except Exception as e:
        log_processing_error(sample_name, "kallisto_species", str(e))
        raise

    # Summarize abundance to gene level
    try:
        summarize_to_gene_level(sample_name, species_name, results_dir)
    except Exception as e:
        log_processing_error(sample_name, "summarize_to_gene_level", str(e))
        raise

    # Run Kallisto quantification for contaminants
    try:
        run_kallisto(fastq1, fastq2, sample_name, species_name, strandedness, results_dir, is_contaminant_run=True)
    except Exception as e:
        log_processing_error(sample_name, "kallisto_contaminant", str(e))
        raise

    # After all steps, check again for output validity
    errors = []
    if not fastp_valid():
        errors.append("fastp output incomplete or missing")
    if not fastqc_valid():
        errors.append("FastQC output incomplete or missing")
    if not star_valid():
        errors.append("STAR output incomplete or missing")
    if not kallisto_valid():
        errors.append("Kallisto output incomplete or missing")
    # In SE mode, do not require fastq2 or paired-end outputs
    if not is_paired_end:
        # Optionally, could relax checks further if any are PE-specific
        pass
    if errors:
        err_msg = f"Sample {sample_name} failed output validation after processing: " + ", ".join(errors)
        log_processing_error(sample_name, "output_validation", err_msg)
        raise RuntimeError(err_msg)

    print(f"Finished processing sample: {sample_name}")
    return strandedness

def order_mycoplasma_report(results_dir, sample_names):
    """Force deterministic row order in the mycoplasma report."""
    path = Path(results_dir) / "mycoplasma_report.tsv"
    if not path.exists():
        return
    frame = pd.read_csv(path, sep="\t")
    if "Sample" not in frame.columns:
        return
    ordered = pd.Categorical(frame["Sample"], categories=sample_names, ordered=True)
    frame = frame.assign(_order=ordered).sort_values("_order", na_position="last")
    frame.drop(columns=["_order"]).to_csv(path, sep="\t", index=False)


def run_enrichment(metadata, results_dir, species_name):
    """Run GSEA, pre-ranked when a condition has fewer than three replicates.

    Returns the mode used so the manifest records which flavour of GSEA produced
    gene_set_enrichment.tsv.
    """
    gmt_dir = REFERENCE_DIR / species_name / "GSEA"
    gmt_files = sorted(gmt_dir.glob("*.gmt"))
    stems = {p.stem.lower() for p in gmt_files}
    skip_full_go = (
        any(s.startswith("gene_ontology") for s in stems)
        and any(s.startswith("biological_process") for s in stems)
        and any(s.startswith("molecular_function") for s in stems)
    )
    if skip_full_go:
        gmt_files = [p for p in gmt_files
                     if not p.stem.lower().startswith("gene_ontology")]

    condition_counts = metadata.groupby('condition')['Sample name'].nunique()
    if (condition_counts < 3).any():
        print("Less than 3 replicates detected for one or more conditions. Running GSEA "
              "pre-ranked analysis.")
        rnk_file = create_rnk_file(results_dir)
        if not rnk_file:
            raise RuntimeError("unable to create the .rnk file for GSEA pre-ranked analysis")
        gsea_results_dir = results_dir / "gsea"
        for gmt_file in gmt_files:
            run_gsea_preranked(rnk_file, gsea_results_dir, gmt_file)
        mode = "preranked"
    else:
        print("3 or more replicates detected for all conditions. Running standard GSEA "
              "workflow.")
        create_gsea_input_file(metadata, results_dir, species_name)
        create_gsea_cls_file(metadata, results_dir, species_name)
        run_gsea(results_dir, "gsea", species_name)
        mode = "standard"

    # GSEA creates an empty mmmDD directory in the working directory on some runs.
    import datetime
    month_day_folder = Path(datetime.datetime.now().strftime("%b%d").lower())
    if month_day_folder.is_dir() and not any(month_day_folder.iterdir()):
        shutil.rmtree(month_day_folder)
        print(f"Removed empty GSEA folder: {month_day_folder}")
    return mode


def write_human_readable_outputs(metadata_file, results_dir, sample_names, species_name):
    """Figures, spreadsheets and the HTML report.

    These are conveniences layered on top of the contract tables, so a failure here is reported
    but does not discard a completed run's results.
    """
    def optional(label, function, *arguments):
        try:
            function(*arguments)
        except Exception as error:  # noqa: BLE001 - cosmetic outputs must not abort the run
            print(f"WARNING: {label} failed: {type(error).__name__}: {error}")

    collate_tpms(results_dir, sample_names, species_name)
    optional("DESeq2 wide merge", merge_deseq2_results, results_dir)
    optional("DESeq2 XLSX", reform_deseq2_results, results_dir)
    optional("de_gene XLSX", tables_de.write_de_gene_xlsx, results_dir)
    optional("transcript TPM collation", collate_transcript_tpms, results_dir, sample_names,
             species_name)
    optional("Sleuth Wald merge", merge_sleuth_wald_results, results_dir)
    optional("Sleuth LRT merge", merge_sleuth_lrt_results, results_dir)
    optional("Sleuth XLSX", reformat_sleuth_results, results_dir)
    optional("mycoplasma check", check_mycoplasma_sequences, results_dir)
    optional("mycoplasma row ordering", order_mycoplasma_report, results_dir, sample_names)
    optional("volcano plot", prepare_volcano, results_dir)

    for sample_name in sample_names:
        optional(f"junction BED for {sample_name}", convert_SJtab_to_bed, sample_name,
                 results_dir)

    optional("top gene/transcript report", report_top_genes_and_transcripts, results_dir)
    optional("top TPM report", report_top_gene_tpms, results_dir)
    optional("alignment summary", summarize_alignments, results_dir, sample_names)
    optional("transcript coverage summary", summarize_transcript_coverage, results_dir,
             sample_names, species_name)
    optional("HTML report", generate_results_html, results_dir, sample_names)


def parse_cli(argv):
    parser = argparse.ArgumentParser(
        prog="rnaseq.py",
        description="Run the RNA-seq pipeline and emit the canonical warehouse tables.",
    )
    parser.add_argument("--version", action="version", version=PIPELINE_VERSION)
    parser.add_argument("reference_dir")
    parser.add_argument("scripts_dir", help="the rnaseq_helper_scripts directory")
    parser.add_argument("results_dir")
    parser.add_argument("--annotation-version",
                        help="override the annotation release recorded in the manifest")
    parser.add_argument("--normalization", default=bigwig.DEFAULT_NORMALIZATION,
                        choices=["cpm"],
                        help="bigWig normalization; must be uniform across all runs")
    parser.add_argument("--skip-transcript-expression", action="store_true",
                        help="write expression_transcript.tsv with headers only (it is the "
                             "largest table)")
    return parser.parse_args(argv)


def build_run_context(metadata, results_dir, run_id, genome_build):
    """Map internal sample names to warehouse sample_ids for every table writer."""
    missing = [c for c in ("sample_id", "experiment_id") if c not in metadata.columns]
    if missing:
        raise SystemExit(
            f"metadata file is missing {missing}. It was generated by an older metadata.py; "
            f"regenerate it so every table can carry a stable sample_id."
        )
    return RunContext(
        run_id=run_id,
        experiment_id=str(metadata["experiment_id"].iloc[0]),
        results_dir=Path(results_dir),
        genome_build=genome_build,
        sample_ids=dict(zip(metadata["Sample name"].astype(str),
                            metadata["sample_id"].astype(str))),
    )


def main():
    """Run every stage, then emit the tables, checksums, manifest and validation result."""
    args = parse_cli(sys.argv[1:])

    global REFERENCE_DIR, scripts, results_dir, species_name

    # Canonical paths so relative CLI args behave the same regardless of cwd / subprocess cwd.
    REFERENCE_DIR = Path(args.reference_dir).resolve()
    scripts = Path(args.scripts_dir).resolve()
    results_dir = Path(args.results_dir).resolve()
    create_layout(results_dir)
    setup_logging(results_dir)
    start_time = time.time()

    # The one metadata file: metadata.py must already have written it here. rnaseq.py reads it
    # and updates it in place (e.g. once RSeQC infers rseqc_measured_strandedness) rather than
    # copying from a separate input.
    metadata_file = results_dir / "metadata.tsv"
    if not metadata_file.is_file():
        print(f"Error: {metadata_file} not found. Run metadata.py against {results_dir} first.")
        sys.exit(1)

    metadata = pd.read_csv(metadata_file, sep='\t')
    required_columns = ['fastq_r1', 'fastq_r2', 'genome_build', 'Sample name', 'condition',
                        'experiment_id']
    missing_columns = [col for col in required_columns if col not in metadata.columns]
    if missing_columns:
        print(f"Error: Missing columns in metadata file: {', '.join(missing_columns)}")
        sys.exit(1)

    species_name = metadata['genome_build'].iloc[0]
    print(f"Species name: {species_name}")

    mart_file = REFERENCE_DIR / species_name / "biomart" / f"{species_name}.mart_export.txt"
    if not mart_file.exists():
        print(f"Mart file {mart_file} not found.")
        sys.exit(1)
    annotation = Annotation(mart_file)

    manifest = provenance.start_run(
        results_dir=results_dir,
        reference_dir=REFERENCE_DIR,
        scripts_dir=scripts,
        experiment_id=str(metadata['experiment_id'].iloc[0]),
        genome_build=species_name,
        investigator=str(metadata.get('investigator', pd.Series(['unknown'])).iloc[0]),
        library_layout="PE" if metadata['fastq_r2'].notnull().any() else "SE",
        parameters={
            "star_args": " ".join(STAR_COMMON_ARGS),
            "kallisto_args": " ".join(KALLISTO_COMMON_ARGS),
            "rmats_args": None,
            "deseq2_alpha": DESEQ2_ALPHA,
            "bigwig_normalization": args.normalization,
            "rmats_fdr_threshold": tables_rmats.RMATS_FDR_THRESHOLD,
            "rmats_inclusion_diff_threshold": tables_rmats.RMATS_INC_DIFF_THRESHOLD,
        },
        annotation_version_override=args.annotation_version,
    )
    manifest.set("gene_id_source", annotation.source)
    print(f"run_id: {manifest.run_id}")

    stage = "startup"
    try:
        ctx = build_run_context(metadata, results_dir, manifest.run_id, species_name)

        # rMATS replicate indices, the GSEA class file and the DESeq2 counts matrix all depend
        # on this ordering; see ordered_sample_names.
        sample_names = ordered_sample_names(metadata)
        grouped_metadata = metadata.groupby('Sample name')

        stage = "per_sample_processing"
        library_types = {}
        for sample_name in sample_names:
            group = grouped_metadata.get_group(sample_name)
            row = group.iloc[-1]
            fastq1 = row['fastq_r1']
            fastq2 = row.get('fastq_r2', None)
            if not fastq1:
                raise RuntimeError(f"Missing FASTQ files for sample {sample_name}.")
            if not (isinstance(fastq2, str) and fastq2 and fastq2.lower() != 'nan'):
                print(f"Single-end sample detected: {sample_name}")
                fastq2 = None
            library_types[sample_name] = process_sample(sample_name, fastq1, fastq2)

        stage = "qc_metrics"
        tables_qc.write_qc_metrics(ctx, results_dir, sample_names, species_name)
        tables_qc.update_metadata_strandedness(results_dir, library_types)
        library_type = next(iter(set(library_types.values())))

        stage = "deseq2"
        generate_deseq2_metadata(metadata, results_dir)
        generate_counts_matrix(sample_names, species_name, results_dir)
        run_deseq2_analysis(results_dir, species_name, REFERENCE_DIR)
        manifest.record_stage("de_gene_rows",
                              tables_de.write_de_gene(ctx, results_dir, annotation))

        stage = "expression"
        counts = tables_expression.write_expression_tables(
            ctx, results_dir, sample_names, species_name, annotation,
            write_transcripts=not args.skip_transcript_expression)
        for key, value in counts.items():
            manifest.record_stage(key, value)
        if args.skip_transcript_expression:
            manifest.record_stage("expression_transcript", "skipped")

        stage = "sleuth"
        sleuth_dir = results_dir / "sleuth"
        sleuth_dir.mkdir(parents=True, exist_ok=True)
        create_sleuth_metadata(metadata_file, sleuth_dir, species_name, results_dir)
        run_sleuth_analysis(sleuth_dir, species_name, REFERENCE_DIR)
        rows, status = tables_de.write_transcript_de(ctx, results_dir, annotation)
        manifest.record_stage("transcript_de_rows", rows)
        manifest.record_stage("sleuth", status)

        stage = "gsea"
        manifest.record_stage("gsea_mode", run_enrichment(metadata, results_dir, species_name))
        rows, status = tables_gsea.write_gene_set_enrichment(ctx, results_dir)
        manifest.record_stage("gene_set_enrichment_rows", rows)
        manifest.record_stage("gsea", status)

        stage = "rmats"
        check_and_prepare_rmats(metadata, results_dir, sample_names, library_types)
        manifest.data["parameters"]["rmats_args"] = run_rmats(
            metadata, results_dir, species_name, library_type)
        generate_splicing_bar_chart(results_dir)
        splicing = tables_rmats.write_splicing_tables(
            ctx, results_dir, results_dir / "rmats" / "rmats_out")
        for key, value in splicing.items():
            manifest.record_stage(f"splicing_{key}" if key != "status" else "rmats", value)

        stage = "junctions"
        junction_counts = tables_junction.write_junctions(ctx, results_dir, sample_names)
        manifest.record_stage("junction_rows", int(sum(junction_counts.values())))

        stage = "bigwig"
        chr_sizes = REFERENCE_DIR / species_name / "STAR" / "chrNameLength.txt"
        bigwig_rows = []
        for sample_name in sample_names:
            bigwig_rows.extend(bigwig.convert_sample_bigwigs(
                ctx, results_dir, sample_name, library_types[sample_name], chr_sizes,
                normalization=args.normalization))
        bigwig.write_bigwig_manifest(results_dir, bigwig_rows)
        remove_wig(results_dir)

        stage = "signal_over_gene"
        rows, status = bigwig.write_signal_over_gene(
            ctx, results_dir, REFERENCE_DIR, annotation, bigwig_rows)
        manifest.record_stage("signal_over_gene_rows", rows)
        manifest.record_stage("signal_over_gene", status)
        organism = str(metadata.get('organism', pd.Series(['unknown'])).iloc[0])
        manifest.record_stage("strand_check", bigwig.assert_strand_assignment(
            results_dir, organism, library_type, annotation))

        stage = "human_readable_outputs"
        write_human_readable_outputs(metadata_file, results_dir, sample_names, species_name)

        stage = "organize"
        log_total_time(start_time, results_dir)
        organize_artifacts(results_dir, sample_names)

        stage = "checksums"
        known_hashes = {row["file_path"]: row["sha256"] for row in bigwig_rows}
        checksums_module.write_checksums(results_dir, known=known_hashes)
    except Exception as error:
        manifest.set("error", f"{type(error).__name__}: {error}")
        manifest.write("failed", failed_stage=stage)
        print(f"RNA-seq pipeline FAILED during stage '{stage}': {error}")
        raise

    manifest.write("success")

    # Finder can drop .DS_Store between organize and here on macOS-mounted volumes.
    for junk in results_dir.rglob(".DS_Store"):
        junk.unlink(missing_ok=True)

    report = validate_outputs.validate(results_dir)
    for warning in report.warnings:
        print(f"validation warning: {warning}")
    for index, failure in enumerate(report.failures, start=1):
        print(f"validation failure {index}: {failure}")
    manifest.set("validation", "pass" if report.ok else "fail")
    manifest.write("success")

    print("RNA-seq pipeline completed successfully.")
    if not report.ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
