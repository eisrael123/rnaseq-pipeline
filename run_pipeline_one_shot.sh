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
#
# Either edit the two defaults below, or pass them in the environment (this is
# how ../tools/run_pipelines.sh calls this script):
#   HOST_FASTQ_ROOT_DIR=/local/Mutu_X_2026-01-01 HOST_OUTPUT_DIR=/local/pipeline_output/Mutu_X_2026-01-01 ./run_pipeline_one_shot.sh
# ---------------------------------------------------------------------------
HOST_FASTQ_ROOT_DIR="${HOST_FASTQ_ROOT_DIR:-/path/to/Model_Experiment/on/your/computer}"
HOST_OUTPUT_DIR="${HOST_OUTPUT_DIR:-/path/to/output/folder/on/your/computer}"

# Reference files: prefer the local copy on the Mac15 (fast local disk); fall back to the
# copy on the lab NAS. Set HOST_REFERENCE_DIR in the environment to force a specific one.
LOCAL_REFERENCE_DIR="/Applications/ngs/pipelines/docker-rnaseq/referenceFiles"
NAS_REFERENCE_DIR="/Volumes/FlemingtonLabMain1/VIRA/rnaseq/referenceFiles"
if [[ -z "${HOST_REFERENCE_DIR:-}" ]]; then
  if [[ -d "${LOCAL_REFERENCE_DIR}" ]]; then
    HOST_REFERENCE_DIR="${LOCAL_REFERENCE_DIR}"
  elif [[ -d "${NAS_REFERENCE_DIR}" ]]; then
    HOST_REFERENCE_DIR="${NAS_REFERENCE_DIR}"
  else
    echo "ERROR: referenceFiles not found at ${LOCAL_REFERENCE_DIR} or ${NAS_REFERENCE_DIR}" >&2
    exit 1
  fi
fi
echo "Using reference files: ${HOST_REFERENCE_DIR}"

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
  python rnaseq.py ${CONTAINER_REFERENCE_DIR} ${CONTAINER_SCRIPTS_DIR} ${CONTAINER_OUTPUT_DIR}"
