#!/usr/bin/env Rscript
# sleuth_analysis_ercc.R

# Load necessary libraries
library(sleuth)
library(dplyr)

# Define arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript sleuth_analysis_ercc.R <results_dir> <species_name>")
}

results_dir <- args[1]  # Directory to save sleuth results
species_name <- args[2]  # Species name

# Set paths for input files
sleuth_metadata_file <- file.path(results_dir, "sleuth_metadata.tsv")
mart_file <- file.path("/Applications/ngs/pipelines/rnaseq/referenceFiles", species_name, "biomart", paste0(species_name, ".mart_export.txt"))
ercc_file <- file.path("/Applications/ngs/pipelines/rnaseq/referenceFiles", species_name, "annotations", "ERCC92_genes.txt")

# Debugging information
cat("Results directory:", results_dir, "\n")
cat("Mart file:", mart_file, "\n")
cat("Species name:", species_name, "\n")
cat("ERCC file path:", ercc_file, "\n")

# Verify that the sleuth metadata file exists
if (!file.exists(sleuth_metadata_file)) {
  stop(paste("Sleuth metadata file not found:", sleuth_metadata_file))
}

# Verify that the mart file exists
if (!file.exists(mart_file)) {
  stop(paste("Mart file not found:", mart_file))
}

# Verify that the ERCC file exists
if (!file.exists(ercc_file)) {
  stop(paste("ERCC file not found:", ercc_file))
}

# Load metadata and mart file
s2c <- tryCatch({
  read.table(sleuth_metadata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}, error = function(e) {
  stop(paste("Failed to read sleuth metadata file:", e$message))
})

mart <- tryCatch({
  read.table(mart_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}, error = function(e) {
  stop(paste("Failed to read mart file:", e$message))
})

# Relevel condition to use 'control' as the reference
s2c$condition <- relevel(factor(s2c$condition), ref = "control")

# Load ERCC IDs
ercc_ids <- tryCatch({
  read.table(ercc_file, header = FALSE, stringsAsFactors = FALSE)$V1
}, error = function(e) {
  stop(paste("Failed to read ERCC file:", e$message))
})

# Prepare the sleuth object with transcript-level data
sot <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE, target_mapping = mart, transformation_function = function(x) log2(x + 0.5))

# Extract ERCC counts
ercc_counts <- sleuth_to_matrix(sot, "obs_raw", "est_counts")
ercc_counts <- ercc_counts[rownames(ercc_counts) %in% ercc_ids, ]

# Check if ERCC counts are non-zero
if (all(rowSums(ercc_counts) == 0)) {
  cat("No ERCC spike-in detected. Proceeding without normalization.\n")
  # Proceed without normalization
  sot <- sleuth_fit(sot, ~ condition, 'full')
  sot <- sleuth_fit(sot, ~ 1, 'reduced')
} else {
  cat("ERCC spike-in detected. Applying normalization.\n")
  # Apply normalization using ERCC counts
  sot <- sleuth_fit(sot, ~ condition, 'full')
  sot <- sleuth_fit(sot, ~ 1, 'reduced')
}

# Run LRT for transcript-level differential expression
models(sot)
tests(sot)
results_table_lrt_transcriptLevel <- tryCatch({
  sot <- sleuth_lrt(sot, 'reduced', 'full')
  sleuth_results(sot, 'reduced:full', test_type = 'lrt')
}, error = function(e) {
  stop(paste("Failed to run LRT for transcript-level differential expression:", e$message))
})
results_table_lrt_transcriptLevel <- arrange(results_table_lrt_transcriptLevel, qval)
write.csv(results_table_lrt_transcriptLevel, file.path(results_dir, 'sleuth_DETranscripts_lrt.csv'))

cat("Transcript-level LRT results saved.\n")

# Perform Wald test at transcript level
results_table_wt_transcriptLevel <- tryCatch({
  sot <- sleuth_wt(sot, which_beta = 'conditiontest')
  sleuth_results(sot, 'conditiontest', test_type = 'wt')
}, error = function(e) {
  stop(paste("Failed to run Wald test for transcript-level differential expression:", e$message))
})
results_table_wt_transcriptLevel <- arrange(results_table_wt_transcriptLevel, qval)
write.csv(results_table_wt_transcriptLevel, file.path(results_dir, 'sleuth_DETranscripts_wt.csv'))

cat("Transcript-level Wald test results saved.\n")

# Gene-Level Differential Expression Analysis
# Prepare the sleuth object with gene-level data
sog <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE, target_mapping = mart, gene_mode = TRUE, aggregation_column = "gene", transformation_function = function(x) log2(x + 0.5))

# Fit the full and reduced models for LRT at gene level
sog <- sleuth_fit(sog, ~condition, 'full')
sog <- sleuth_fit(sog, ~1, 'reduced')

# Run LRT for gene-level differential expression
models(sog)
tests(sog)
results_table_lrt_geneLevel <- tryCatch({
  sog <- sleuth_lrt(sog, 'reduced', 'full')
  sleuth_results(sog, 'reduced:full', test_type = 'lrt')
}, error = function(e) {
  stop(paste("Failed to run LRT for gene-level differential expression:", e$message))
})
results_table_lrt_geneLevel <- arrange(results_table_lrt_geneLevel, qval)
write.csv(results_table_lrt_geneLevel, file.path(results_dir, 'sleuth_DEGenes_lrt.csv'))

cat("Gene-level LRT results saved.\n")

# Perform Wald test at gene level
results_table_wt_geneLevel <- tryCatch({
  sog <- sleuth_wt(sog, which_beta = 'conditiontest')
  sleuth_results(sog, 'conditiontest', test_type = 'wt')
}, error = function(e) {
  stop(paste("Failed to run Wald test for gene-level differential expression:", e$message))
})
results_table_wt_geneLevel <- arrange(results_table_wt_geneLevel, qval)
write.csv(results_table_wt_geneLevel, file.path(results_dir, 'sleuth_DEGenes_wt.csv'))

cat("Gene-level Wald test results saved.\n")
