#!/usr/bin/env python3
# metadata.py
# metadata.py is same script as scriptb_bkup/metadata_620a.py

# Usage: metadata.py /PATH/TO/FASTQ_FILES <reference_name> <investigator_name> <PE|SE>

import pandas as pd
from pathlib import Path
import re
import sys
import os
from datetime import datetime
import logging
from tqdm import tqdm
import gzip
from itertools import chain

REQUIRED_ENV = "rnaseqpipeline"

# Check if running in the right conda env
conda_env = os.environ.get("CONDA_DEFAULT_ENV")
if conda_env != REQUIRED_ENV:
    sys.stderr.write(
        f"\nERROR: This script must be run inside the '{REQUIRED_ENV}' conda environment.\n"
        f"Currently active environment: {conda_env or 'None'}\n"
        f"Please run:\n\n    conda activate {REQUIRED_ENV}\n\n"
    )
    sys.exit(1)

REFERENCE_DIR = Path('/Applications/ngs/pipelines/rnaseq/referenceFiles')

def get_available_references():
    if not REFERENCE_DIR.exists() or not REFERENCE_DIR.is_dir():
        print("Error: 'referenceFiles' directory not found or not a directory.")
        sys.exit(1)
    return [d.name for d in REFERENCE_DIR.iterdir() if d.is_dir()]

def check_reference_availability(reference_name):
    """Check if the specified reference is available."""
    available_references = get_available_references()
    if reference_name not in available_references:
        print(f"Reference '{reference_name}' not available. Available references are: {', '.join(available_references)}")
        return False
    return True

def parse_fastq_files(parent_dir, experiment_type, reference_name):
    """
    Traverse the parent directory to find FASTQ files and extract metadata.

    Parameters:
        parent_dir (Path): Path to the parent directory containing FASTQ files.
        experiment_type (str): 'PE' for paired-end or 'SE' for single-end.
        reference_name (str): Reference genome name.

    Returns:
        List[Dict]: A list of dictionaries containing metadata for each sample.
    """
    import re
    from itertools import chain

    metadata_list = []

    # Define patterns for paired-end reads
    read1_patterns = [
        re.compile(r'(.*)_1\.fastq\.gz'),
        re.compile(r'(.*)_R1_001\.fastq\.gz'),
        re.compile(r'(.*)_1\.fq\.gz'),
        re.compile(r'(.*)_R1_001\.fq\.gz'),
        re.compile(r'(.*)_R1\.fastq\.gz'),
        re.compile(r'(.*)_R1\.fq\.gz')
    ]
    read2_patterns = [
        re.compile(r'(.*)_2\.fastq\.gz'),
        re.compile(r'(.*)_R2_001\.fastq\.gz'),
        re.compile(r'(.*)_2\.fq\.gz'),
        re.compile(r'(.*)_R2_001\.fq\.gz'),
        re.compile(r'(.*)_R2\.fastq\.gz'),
        re.compile(r'(.*)_R2\.fq\.gz')
    ]

    # Define patterns for single-end reads (*.fastq.gz and *.fq.gz)
    single_patterns = [
        re.compile(r'(.*)\.fastq\.gz'),
        re.compile(r'(.*)\.fq\.gz')
    ]

    # Iterate through subdirectories
    for subdir in parent_dir.iterdir():
        if subdir.is_dir():
            condition = 'test' if 'test' in subdir.name.lower() else 'control'
            control_flag = 'no' if condition == 'test' else 'yes'
            experiment_name = parent_dir.name

            if experiment_type == 'PE':
                paired_reads = {}
                for fastq_file in chain(subdir.rglob('*.fastq.gz'), subdir.rglob('*.fq.gz')):
                    for pattern in read1_patterns:
                        match1 = pattern.match(fastq_file.name)
                        if match1:
                            sample_id = match1.group(1)
                            paired_reads.setdefault(sample_id, {})['Read 1'] = str(fastq_file)
                            break
                    for pattern in read2_patterns:
                        match2 = pattern.match(fastq_file.name)
                        if match2:
                            sample_id = match2.group(1)
                            paired_reads.setdefault(sample_id, {})['Read 2'] = str(fastq_file)
                            break
                for sample_id, reads in paired_reads.items():
                    if 'Read 1' in reads and 'Read 2' in reads:
                        metadata_entry = {
                            'Path Read 1': reads['Read 1'],
                            'Path Read 2': reads['Read 2'],
                            'Species': reference_name,
                            'Sample name': sample_id,
                            'Condition': condition,
                            'Control?': control_flag,
                            'Experiment name': experiment_name
                        }
                        metadata_list.append(metadata_entry)
                    else:
                        logging.warning(f"Sample {sample_id} in {subdir} is missing paired reads. Skipping.")

            elif experiment_type == 'SE':
                for fastq_file in chain(subdir.rglob('*.fastq.gz'), subdir.rglob('*.fq.gz')):
                    matched = False
                    for pattern in single_patterns:
                        match = pattern.match(fastq_file.name)
                        if match:
                            sample_id = match.group(1)
                            metadata_entry = {
                                'Path Read 1': str(fastq_file),
                                'Path Read 2': '',
                                'Species': reference_name,
                                'Sample name': sample_id,
                                'Condition': condition,
                                'Control?': control_flag,
                                'Experiment name': experiment_name
                            }
                            metadata_list.append(metadata_entry)
                            matched = True
                            break
                    if not matched:
                        logging.warning(f"FASTQ file {fastq_file.name} does not match expected single-end pattern. Skipping.")

    return metadata_list

# Verify command-line arguments
if len(sys.argv) == 1:
    print("Usage: metadata.py <parent_dir> <reference_name> <investigator_name> <PE|SE>")
    print(f"Available references are: {', '.join(get_available_references())}")
    sys.exit(0)

if len(sys.argv) < 5:
    print("Usage: metadata.py <parent_dir> <reference_name> <investigator_name> <PE|SE>")
    print(f"Available references are: {', '.join(get_available_references())}")
    sys.exit(1)

# Extract command-line arguments
parent_dir = Path(sys.argv[1])
reference_name = sys.argv[2]
investigator = sys.argv[3]
experiment_type = sys.argv[4].upper()  # Expect PE or SE

# Check parent_dir name format: must be "model_experiment"
if len(parent_dir.name.split('_')) != 2:
    print('Oh no! The metadata file was not generated because the FASTQ files directory name should be in the format "model_experiment", \033[1;3mwith only ONE underscore\033[0m included in the directory name. Try again!')
    sys.exit(1)

if experiment_type not in ["PE", "SE"]:
    print("Error: Experiment type must be either 'PE' (paired-end) or 'SE' (single-end).")
    sys.exit(1)

# Check if the specified reference directory exists
if not check_reference_availability(reference_name):
    sys.exit(1)

# Set up logging
log_filename = f'metadata_log_{datetime.now().strftime("%Y%m%d_%H%M%S")}.txt'
logging.basicConfig(filename=log_filename, level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
logging.info("Starting metadata generation.")
logging.info(f"Parent Directory: {parent_dir}")
logging.info(f"Reference Name: {reference_name}")
logging.info(f"Investigator: {investigator}")
logging.info(f"Experiment Type: {experiment_type}")

# Parse FASTQ files and extract metadata
metadata_entries = parse_fastq_files(parent_dir, experiment_type, reference_name)

if not metadata_entries:
    logging.error("No valid samples found. Exiting.")
    print("Error: No valid samples found. Check the log for details.")
    sys.exit(1)

# Create DataFrame
metadata_df = pd.DataFrame(metadata_entries)

# Sort: controls first (lexicographically), then tests (lexicographically)
metadata_df['Condition_sort'] = metadata_df['Condition'].apply(lambda x: 0 if x == 'control' else 1)
metadata_df = metadata_df.sort_values(['Condition_sort', 'Sample name'], ascending=[True, True])
metadata_df = metadata_df.drop(columns=['Condition_sort'])

# Assign ConditionPrefix after sorting
def get_condition_prefix(row):
    subdir = Path(row['Path Read 1']).parts[-2]
    if subdir.lower().startswith('cntl'):
        return 'cntl'
    elif subdir.lower().startswith('test'):
        return 'test'
    else:
        return row['Condition']

metadata_df = metadata_df.copy()
metadata_df['ConditionPrefix'] = metadata_df.apply(get_condition_prefix, axis=1)
# Now assign replicate numbers in sorted order
metadata_df['ReplicateNum'] = metadata_df.groupby('ConditionPrefix').cumcount() + 1
metadata_df['ConditionReplicate'] = metadata_df['ConditionPrefix'] + metadata_df['ReplicateNum'].astype(str)

# Reorder columns as per requirement (include ConditionReplicate)
metadata_df = metadata_df[['Path Read 1', 'Path Read 2', 'Species', 'Sample name', 'Condition', 'Control?', 'ConditionReplicate', 'Experiment name']]

# Assign Species (all entries have the same species)
metadata_df['Species'] = reference_name

# Export to TSV
timestamp = datetime.now().strftime("%m%d%Y_%H%M%S")
metadata_filename = f'{investigator}_metadata_{timestamp}.tsv'
metadata_df.to_csv(metadata_filename, sep='\t', index=False)

logging.info(f"Metadata file generated: {metadata_filename}")
print(f"Metadata file successfully generated: {metadata_filename}")
