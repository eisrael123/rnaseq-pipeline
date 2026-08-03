#!/usr/bin/env bash
set -euo pipefail

# ---------------------------------------------------------------------------
# Edit these two, and nothing else.
#
# Everything about the experiment itself -- cell line, perturbation, genome
# build, investigator, PE/SE, and the rest -- comes from the input_args.json
# that metadata_form.html produced, which must be sitting at the top level of
# HOST_FASTQ_ROOT_DIR alongside the cntl/ and test/ subdirectories.
# ---------------------------------------------------------------------------
HOST_FASTQ_ROOT_DIR="/path/to/Model_Experiment/on/your/computer"
HOST_OUTPUT_DIR="/path/to/output/folder/on/your/computer"

# Fixed on the Mac15. Change it only on a machine where referenceFiles lives elsewhere.
HOST_REFERENCE_DIR="/Applications/ngs/pipelines/docker-rnaseq/referenceFiles"

# Container paths (do not change unless you also change the mounted paths below).
IMAGE="rnaseqpipeline:latest"
CONTAINER_REFERENCE_DIR="/data/referenceFiles"
CONTAINER_OUTPUT_DIR="/data/output"
CONTAINER_SCRIPTS_DIR="/work/rnaseq_helper_scripts"

# metadata.py reads experiment_id from the *basename of the path it is given*, and that id is
# what sample_id and comparison_id are built from. Mounting to a fixed name would stamp every
# run with that name, so the mount point has to carry the real experiment directory name.
CONTAINER_FASTQ_ROOT_DIR="/data/$(basename "${HOST_FASTQ_ROOT_DIR}")"

# Fail here rather than three minutes into a container that is going to exit anyway.
QUICK_INPUT_FILE="${HOST_FASTQ_ROOT_DIR}/input_args.json"
if [[ ! -f "${QUICK_INPUT_FILE}" ]]; then
  echo "ERROR: no input_args.json found at ${QUICK_INPUT_FILE}" >&2
  echo "Fill out metadata_form.html and drag the file it downloads into" >&2
  echo "${HOST_FASTQ_ROOT_DIR}, next to the cntl/ and test/ subdirectories." >&2
  exit 1
fi

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
  python metadata.py \
    --root-fastq-dir ${CONTAINER_FASTQ_ROOT_DIR} \
    --output-dir ${CONTAINER_OUTPUT_DIR} \
    --quick-input \
    --non-interactive && \
  python rnaseq.py ${CONTAINER_OUTPUT_DIR}/*_metadata_*.tsv ${CONTAINER_REFERENCE_DIR} ${CONTAINER_SCRIPTS_DIR} ${CONTAINER_OUTPUT_DIR}"
