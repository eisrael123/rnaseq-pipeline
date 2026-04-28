#!/usr/bin/env bash
set -euo pipefail

# One-shot reference file generation in Docker.
# Edit these variables before running.

IMAGE="rnaseqpipeline:latest"
HOST_REFERENCEFILES_DIR="/Users/mac14/ethan/rnaseq/referenceFiles"
HOST_RAW_REFERENCE_DOWNLOADS="/path/to/raw_reference_downloads"

SPECIES="myNewSpecies"
THREADS="16"
SJDB_OVERHANG="149"

# Example input source filenames under ${HOST_RAW_REFERENCE_DOWNLOADS}
PRIMARY_GENOME_FA="genome.primary.fa"
SPIKE_OR_VIRUS_FA="spikeins_or_virus.fa"
PRIMARY_GENOME_GTF="genome.primary.gtf"
SPIKE_OR_VIRUS_GTF="spikeins_or_virus.gtf"
BIOMART_FILE="${SPECIES}.mart_export.txt"
ERCC_FILE="ERCC92_genes.txt"

docker run --rm \
  -v "${HOST_REFERENCEFILES_DIR}:/data/referenceFiles" \
  -v "${HOST_RAW_REFERENCE_DOWNLOADS}:/data/src:ro" \
  -w /work \
  "${IMAGE}" \
  bash -lc '
    set -euo pipefail
    conda activate rnaseqpipeline
    REFROOT="/data/referenceFiles/'"${SPECIES}"'"
    mkdir -p "${REFROOT}"/{annotations,biomart,fasta,GSEA,kallisto,STAR}

    cat "/data/src/'"${PRIMARY_GENOME_FA}"'" "/data/src/'"${SPIKE_OR_VIRUS_FA}"'" > "${REFROOT}/fasta/'"${SPECIES}"'.fa"
    samtools faidx "${REFROOT}/fasta/'"${SPECIES}"'.fa"

    cat "/data/src/'"${PRIMARY_GENOME_GTF}"'" "/data/src/'"${SPIKE_OR_VIRUS_GTF}"'" > "${REFROOT}/annotations/'"${SPECIES}"'.gtf"
    gtf2bed < "${REFROOT}/annotations/'"${SPECIES}"'.gtf" > "${REFROOT}/annotations/'"${SPECIES}"'.bed"

    gffread -g "${REFROOT}/fasta/'"${SPECIES}"'.fa" "${REFROOT}/annotations/'"${SPECIES}"'.gtf" -w "${REFROOT}/fasta/'"${SPECIES}"'.transcriptome.fa"
    samtools faidx "${REFROOT}/fasta/'"${SPECIES}"'.transcriptome.fa"

    STAR --runThreadN '"${THREADS}"' --runMode genomeGenerate --genomeDir "${REFROOT}/STAR" --genomeFastaFiles "${REFROOT}/fasta/'"${SPECIES}"'.fa" --sjdbGTFfile "${REFROOT}/annotations/'"${SPECIES}"'.gtf" --sjdbOverhang '"${SJDB_OVERHANG}"'

    kallisto index -i "${REFROOT}/kallisto/'"${SPECIES}"'.kallisto.idx" "${REFROOT}/fasta/'"${SPECIES}"'.transcriptome.fa"
    kallisto index -i "${REFROOT}/kallisto/all_contaminants.kallisto.idx" "${REFROOT}/fasta/all_contaminants.fa"

    cp "/data/src/biomart/'"${BIOMART_FILE}"'" "${REFROOT}/biomart/"
    cp "/data/src/annotations/'"${ERCC_FILE}"'" "${REFROOT}/annotations/"
    cp /data/src/GSEA/*.gmt "${REFROOT}/GSEA/"

    test -f "${REFROOT}/annotations/'"${SPECIES}"'.gtf"
    test -f "${REFROOT}/annotations/'"${SPECIES}"'.bed"
    test -f "${REFROOT}/annotations/ERCC92_genes.txt"
    test -f "${REFROOT}/fasta/'"${SPECIES}"'.fa"
    test -f "${REFROOT}/fasta/'"${SPECIES}"'.fa.fai"
    test -f "${REFROOT}/fasta/'"${SPECIES}"'.transcriptome.fa"
    test -f "${REFROOT}/fasta/'"${SPECIES}"'.transcriptome.fa.fai"
    test -f "${REFROOT}/STAR/Genome"
    test -f "${REFROOT}/STAR/SA"
    test -f "${REFROOT}/STAR/SAindex"
    test -f "${REFROOT}/STAR/chrLength.txt"
    test -f "${REFROOT}/STAR/chrName.txt"
    test -f "${REFROOT}/STAR/chrNameLength.txt"
    test -f "${REFROOT}/STAR/chrStart.txt"
    test -f "${REFROOT}/STAR/exonGeTrInfo.tab"
    test -f "${REFROOT}/STAR/exonInfo.tab"
    test -f "${REFROOT}/STAR/geneInfo.tab"
    test -f "${REFROOT}/STAR/genomeParameters.txt"
    test -f "${REFROOT}/STAR/sjdbInfo.txt"
    test -f "${REFROOT}/STAR/sjdbList.fromGTF.out.tab"
    test -f "${REFROOT}/STAR/sjdbList.out.tab"
    test -f "${REFROOT}/STAR/transcriptInfo.tab"
    test -f "${REFROOT}/kallisto/'"${SPECIES}"'.kallisto.idx"
    test -f "${REFROOT}/kallisto/all_contaminants.kallisto.idx"
    test -f "${REFROOT}/biomart/'"${SPECIES}"'.mart_export.txt"
    ls "${REFROOT}/GSEA/"*.gmt >/dev/null

    echo "Reference generation complete: ${REFROOT}"
  '
