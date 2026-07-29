#!/usr/bin/env bash
set -euo pipefail

# Edit these pipeline variables before running.
HOST_REFERENCE_DIR="/path/to/referenceFiles/on/your/computer"
HOST_FASTQ_ROOT_DIR="/path/to/Model_Experiment/on/your/computer"
HOST_OUTPUT_DIR="/path/to/output/folder/on/your/computer"

SPECIES_NAME="<species_name>"
INVESTIGATOR_NAME="<investigator_name>"
EXPERIMENT_TYPE="<PE|SE>"

# Experiment structure recorded in metadata.tsv. metadata.py validates these against the
# controlled vocabulary in rnaseq_helper_scripts/vocab.py and exits non-zero on an unknown
# value. Use NA where a field does not apply.
CELL_LINE="<Mutu|Akata|DG75|HepG2|Raji|SNU719|BCBL1|HEK293>"
PERTURBATION_TYPE="<transfection|siRNA|drug|BCR-crosslink|none>"
PERTURBATION_TARGET="<e.g. BMRF1, or NA>"
PERTURBATION_DOSE="<e.g. 100nM, or NA>"
TIMEPOINT_HOURS="<e.g. 24, or NA>"
SEQUENCING_RUN_DATE="<YYYY-MM-DD, or NA>"
NOTES=""

# Container paths (do not change unless you also change mounted paths below).
IMAGE="rnaseqpipeline:latest"
CONTAINER_REFERENCE_DIR="/data/referenceFiles"
CONTAINER_FASTQ_ROOT_DIR="/data/Model_Experiment"
CONTAINER_OUTPUT_DIR="/data/output"
CONTAINER_SCRIPTS_DIR="/work/rnaseq_helper_scripts"

# Recorded in run_manifest.json as docker_image_digest so a result can be traced to the exact
# image that produced it.
IMAGE_ID="$(docker image inspect --format '{{.Id}}' "${IMAGE}")"

docker run --rm \
  -v "${HOST_REFERENCE_DIR}:${CONTAINER_REFERENCE_DIR}:ro" \
  -v "${HOST_FASTQ_ROOT_DIR}:${CONTAINER_FASTQ_ROOT_DIR}:ro" \
  -v "${HOST_OUTPUT_DIR}:${CONTAINER_OUTPUT_DIR}" \
  -e "RNASEQ_DOCKER_IMAGE_DIGEST=${IMAGE_ID}" \
  -w /work \
  "${IMAGE}" \
  bash -lc "conda activate rnaseqpipeline && \
  python metadata.py ${CONTAINER_FASTQ_ROOT_DIR} ${CONTAINER_REFERENCE_DIR} ${SPECIES_NAME} ${INVESTIGATOR_NAME} ${EXPERIMENT_TYPE} ${CONTAINER_OUTPUT_DIR} \
    --non-interactive \
    --cell-line '${CELL_LINE}' \
    --perturbation-type '${PERTURBATION_TYPE}' \
    --perturbation-target '${PERTURBATION_TARGET}' \
    --perturbation-dose '${PERTURBATION_DOSE}' \
    --timepoint-hours '${TIMEPOINT_HOURS}' \
    --sequencing-run-date '${SEQUENCING_RUN_DATE}' \
    --notes '${NOTES}' && \
  python rnaseq.py ${CONTAINER_OUTPUT_DIR}/${INVESTIGATOR_NAME}_metadata_*.tsv ${CONTAINER_REFERENCE_DIR} ${CONTAINER_SCRIPTS_DIR} ${CONTAINER_OUTPUT_DIR}"
