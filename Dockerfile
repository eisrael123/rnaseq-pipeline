# Set miniconda as the base image
FROM continuumio/miniconda3

WORKDIR /work

# Create conda environment from YAML 
COPY rnaseqpipeline.yml /tmp/rnaseqpipeline.yml
RUN conda env create -f /tmp/rnaseqpipeline.yml -n rnaseqpipeline \
  && conda clean -afy

# Use the pipeline env by default (pandas, STAR, etc. live here — not in conda "base").
# Interactive bash from miniconda images often runs `conda activate base`, which prepends
# /opt/conda/bin and shadows this PATH; disable that so `python` stays rnaseqpipeline's.
ENV CONDA_DEFAULT_ENV=rnaseqpipeline
ENV PATH="/opt/conda/envs/rnaseqpipeline/bin:${PATH}"

# Copy scripts into image
COPY . /work