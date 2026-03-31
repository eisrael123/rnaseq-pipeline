#!/usr/bin/env Rscript
# deseq2_analysis_ercc.v5.R

# Assume 'results_dir' and 'species_name' are passed as command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 2){
  stop("Usage: Rscript deseq2_analysis_ercc.R <results_dir> <species_name>")
}
results_dir <- args[1]
species_name <- args[2]

# Correctly construct the paths to the input files
count_matrix_file <- file.path(results_dir, "deseq2", "deseq2_counts_matrix.tsv")
metadata_file <- file.path(results_dir, "deseq2", "deseq2_metadata.tsv")

# Verify that the input files exist
if(!file.exists(count_matrix_file)){
  stop(paste("Counts matrix file not found:", count_matrix_file))
}

if(!file.exists(metadata_file)){
  stop(paste("Metadata file not found:", metadata_file))
}

# Proceed with DESeq2 analysis
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(reshape2)

# Read in the count matrix and metadata with specified tab separator
counts <- tryCatch({
    as.matrix(read.table(count_matrix_file, header = TRUE, row.names = 1, sep = "\t"))
}, error = function(e) {
    stop(paste("Failed to read count matrix:", e$message))
})

metadata <- tryCatch({
    read.table(metadata_file, header = TRUE, row.names = 1, sep = "\t")
}, error = function(e) {
    stop(paste("Failed to read metadata file:", e$message))
})

# Debugging: Print dimensions of the data frames
print(paste("Counts matrix dimensions:", dim(counts)[1], "genes,", dim(counts)[2], "samples"))
print(paste("Metadata dimensions:", dim(metadata)[1], "samples,", dim(metadata)[2], "conditions"))

# Ensure that the sample names in metadata match the column names in counts
colnames(counts) <- sub("^X", "", colnames(counts))
# Convert dots back to hyphens to match metadata format
colnames(counts) <- gsub("\\.", "-", colnames(counts))

print("Column names in counts matrix:")
print(colnames(counts))

print("Row names in metadata:")
print(rownames(metadata))

## Also convert dots to hyphens in metadata row names to match counts
rownames(metadata) <- gsub("\\.", "-", rownames(metadata))

print("Row names in metadata after conversion:")
print(rownames(metadata))

if(!all(colnames(counts) %in% rownames(metadata))){
  stop("Sample names in counts matrix and metadata do not match even after conversion.")
}

# Path to the ERCC genes file
ercc_genes_file <- file.path("/Applications/ngs/pipelines/rnaseq/referenceFiles", species_name, "annotations", "ERCC92_genes.txt")

# Verify that the ERCC genes file exists
if(!file.exists(ercc_genes_file)){
  stop(paste("ERCC genes file not found:", ercc_genes_file))
}

# Read the ERCC genes
ercc_genes <- tryCatch({
    readLines(ercc_genes_file)
}, error = function(e) {
    stop(paste("Failed to read ERCC genes file:", e$message))
})

# Subset the count matrix to include only ERCC genes
ercc_counts <- counts[rownames(counts) %in% ercc_genes, ]

# Check if ERCC counts are zero across all samples
if(all(ercc_counts == 0)){
  print("ERCC counts are zero across all samples. Using default normalization method.")
  use_ercc <- FALSE
} else {
  use_ercc <- TRUE
}

# Convert design variables to factors
metadata$condition <- as.factor(metadata$condition)

# Check for samples with zero counts across all genes
zero_count_samples <- colnames(counts)[colSums(counts) == 0]
if(length(zero_count_samples) > 0){
  stop(paste("Samples with zero counts across all genes:", paste(zero_count_samples, collapse = ", ")))
}

# Proceed with DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ condition)

if(use_ercc){
  # Check for samples with NA size factors due to zero ERCC counts
  size_factors <- estimateSizeFactorsForMatrix(counts, controlGenes = rownames(counts) %in% ercc_genes)
  if(any(is.na(size_factors))){
    print("NA size factors detected. Using default normalization method.")
    use_ercc <- FALSE
  } else {
    # Plot ERCC counts before normalization
    ercc_counts_before <- counts[rownames(counts) %in% ercc_genes, ]
    ercc_counts_before_melt <- melt(ercc_counts_before)
    colnames(ercc_counts_before_melt) <- c("Gene", "Sample", "Count")

    p_before <- ggplot(ercc_counts_before_melt, aes(x = Sample, y = Count, color = Gene)) +
      geom_point() +
      theme_minimal() +
      ggtitle("ERCC Counts Before Normalization") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis text
            legend.text = element_text(size = 8))  # Decrease legend text size

    ggsave(file.path(results_dir, "deseq2", "ERCC_counts_before_normalization.png"), plot = p_before)

    # Estimate size factors using ERCC genes
    dds <- estimateSizeFactors(dds, controlGenes = rownames(dds) %in% ercc_genes)
  }
}

if(!use_ercc){
  # Use default normalization method
  dds <- estimateSizeFactors(dds)
}

# Get normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

# Subset to exclude ERCC genes
endogenous_counts <- normalized_counts[!rownames(normalized_counts) %in% ercc_genes, ]

# Write normalized counts to a file
write.table(endogenous_counts, file = file.path(results_dir, "deseq2", "deseq2_countsNormalized_matrix.tsv"), sep = "\t", quote = FALSE, col.names = NA)

if(use_ercc){
  # Plot ERCC counts after normalization
  ercc_counts_after <- normalized_counts[rownames(normalized_counts) %in% ercc_genes, ]
  ercc_counts_after_melt <- melt(ercc_counts_after)
  colnames(ercc_counts_after_melt) <- c("Gene", "Sample", "Count")

  p_after <- ggplot(ercc_counts_after_melt, aes(x = Sample, y = Count, color = Gene)) +
    geom_point() +
    theme_minimal() +
    ggtitle("ERCC Counts After Normalization") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis text
          legend.text = element_text(size = 8))  # Decrease legend text size

  ggsave(file.path(results_dir, "deseq2", "ERCC_counts_after_normalization.png"), plot = p_after)
}

# Continue with DESeq2 analysis
dds <- DESeq(dds)

# Save results
res <- results(dds)
write.csv(as.data.frame(res), file = file.path(results_dir, "deseq2", "deseq2_results.csv"))
print("DESeq2 results saved successfully.")

# Generate PCA plot
rld <- rlog(dds, blind = FALSE)
pcaData <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = name)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of RNA-seq Samples") +
  theme_minimal()

# Save the PCA plot
ggsave(filename = file.path(results_dir, "deseq2", "PCA_plot.png"))
print("PCA plot saved successfully.")