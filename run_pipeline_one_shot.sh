#!/usr/bin/env bash
set -euo pipefail

# One-shot metadata + rnaseq pipeline run in Docker.
# Edit these variables before running.

IMAGE="rnaseqpipeline:latest"
HOST_REFERENCE_DIR="/path/to/referenceFiles/on/your/computer"
HOST_FASTQ_ROOT_DIR="/path/to/Model_Experiment/on/your/computer"
HOST_OUTPUT_DIR="/path/to/output/folder/on/your/computer"

SPECIES_NAME="<species_name>"
INVESTIGATOR_NAME="<investigator_name>"
EXPERIMENT_TYPE="<PE|SE>"

# Container paths (do not change unless you also change mounted paths below).
CONTAINER_REFERENCE_DIR="/data/referenceFiles"
CONTAINER_FASTQ_ROOT_DIR="/data/Model_Experiment"
CONTAINER_OUTPUT_DIR="/data/output"
CONTAINER_SCRIPTS_DIR="/work/scripts_for_rnaseq"

docker run --rm \
  -v "${HOST_REFERENCE_DIR}:${CONTAINER_REFERENCE_DIR}:ro" \
  -v "${HOST_FASTQ_ROOT_DIR}:${CONTAINER_FASTQ_ROOT_DIR}:ro" \
  -v "${HOST_OUTPUT_DIR}:${CONTAINER_OUTPUT_DIR}" \
  -w /work \
  "${IMAGE}" \
  bash -lc "conda activate rnaseqpipeline && \
  python metadata.py ${CONTAINER_FASTQ_ROOT_DIR} ${CONTAINER_REFERENCE_DIR} ${SPECIES_NAME} ${INVESTIGATOR_NAME} ${EXPERIMENT_TYPE} ${CONTAINER_OUTPUT_DIR} && \
  python rnaseq.py ${CONTAINER_OUTPUT_DIR}/${INVESTIGATOR_NAME}_metadata_*.tsv ${CONTAINER_REFERENCE_DIR} ${CONTAINER_SCRIPTS_DIR} ${CONTAINER_OUTPUT_DIR}"
