#!/usr/bin/env Rscript
# =============================================================================
# Script 01: Download and Prepare TCGA-BRCA RNA-seq Data
# =============================================================================
#
# Description:
#   Fully automated data acquisition from the GDC using TCGAbiolinks.
#   Downloads STAR-Counts RNA-seq data for TCGA-BRCA, extracts the count
#   matrix and sample metadata, and saves processed .rds files for
#   downstream analysis.
#
# Modes:
#   - "quick" (default): ~100 tumor + all available normal samples (~113)
#     Runs in < 10 minutes. Ideal for testing and GitHub demo.
#   - "full": All available tumor + normal samples (~1,100 tumor + ~113 normal)
#     Requires more time and disk space (~20 GB).
#
# Usage:
#   Rscript scripts/01_download_and_prepare.R            # quick-test mode
#   Rscript scripts/01_download_and_prepare.R full        # full mode
#
# Output:
#   data/counts_matrix_<mode>.rds   — integer count matrix (genes × samples)
#   data/col_data_<mode>.rds        — sample metadata (data.frame)
#   data/se_object_<mode>.rds       — full SummarizedExperiment object
#
# Author: DTU Bioinformatics MSc | 2026
# =============================================================================

# ── 0. Configuration ─────────────────────────────────────────────────────────

set.seed(123)

# Parse command-line argument for run mode
args <- commandArgs(trailingOnly = TRUE)
RUN_MODE <- ifelse(length(args) > 0 && tolower(args[1]) == "full", "full", "quick")

cat("============================================================\n")
cat("  TCGA-BRCA RNA-seq Data Download & Preparation Pipeline\n")
cat(sprintf("  Mode: %s\n", toupper(RUN_MODE)))
cat(sprintf("  Date: %s\n", Sys.time()))
cat("============================================================\n\n")

# Directories
PROJECT_DIR <- getwd()
DATA_DIR    <- file.path(PROJECT_DIR, "data")
GDC_DIR     <- file.path(DATA_DIR, "GDCdata")

# Create directories if they don't exist
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(GDC_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(PROJECT_DIR, "results"), showWarnings = FALSE, recursive = TRUE)

# Quick-test sample sizes
QUICK_TUMOR_N  <- 100
QUICK_NORMAL_N <- 100  # will use all available if fewer exist (BRCA has ~113)

# Output file paths
suffix          <- ifelse(RUN_MODE == "full", "full", "quick")
COUNTS_FILE     <- file.path(DATA_DIR, paste0("counts_matrix_", suffix, ".rds"))
COLDATA_FILE    <- file.path(DATA_DIR, paste0("col_data_", suffix, ".rds"))
SE_FILE         <- file.path(DATA_DIR, paste0("se_object_", suffix, ".rds"))

# ── 1. Check if data already exists ──────────────────────────────────────────

if (file.exists(COUNTS_FILE) && file.exists(COLDATA_FILE)) {
  cat("✓ Processed data files already exist:\n")
  cat(sprintf("  Counts: %s\n", COUNTS_FILE))
  cat(sprintf("  ColData: %s\n", COLDATA_FILE))
  
  counts <- readRDS(COUNTS_FILE)
  col_data <- readRDS(COLDATA_FILE)
  
  cat(sprintf("\n  Dimensions: %d genes × %d samples\n", nrow(counts), ncol(counts)))
  cat(sprintf("  Tumor samples:  %d\n", sum(col_data$sample_type == "Primary solid Tumor" |
                                              col_data$sample_type == "Primary Tumor")))
  cat(sprintf("  Normal samples: %d\n", sum(col_data$sample_type == "Solid Tissue Normal")))
  cat("\n→ Skipping download. Delete these files to re-download.\n")
  cat("============================================================\n")
  quit(save = "no", status = 0)
}

# ── 2. Load required packages ────────────────────────────────────────────────

cat("Loading required packages...\n")

required_packages <- c("TCGAbiolinks", "SummarizedExperiment", "dplyr", "stringr")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf(
      "Package '%s' is not installed. Please install it first:\n  BiocManager::install('%s')",
      pkg, pkg
    ))
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

cat("✓ All packages loaded successfully.\n\n")

# ── 3. Query the GDC ─────────────────────────────────────────────────────────

cat("Step 1/5: Querying the GDC for TCGA-BRCA STAR-Counts data...\n")

# Build query for all TCGA-BRCA RNA-Seq HTSeq-Counts (STAR workflow)
query <- tryCatch({
  GDCquery(
    project            = "TCGA-BRCA",
    data.category      = "Transcriptome Profiling",
    data.type          = "Gene Expression Quantification",
    workflow.type      = "STAR - Counts",
    sample.type        = c("Primary Tumor", "Solid Tissue Normal")
  )
}, error = function(e) {
  cat("\n✗ GDC query failed. Error message:\n")
  cat(paste0("  ", e$message, "\n"))
  cat("\nTroubleshooting:\n")
  cat("  1. Check your internet connection\n")
  cat("  2. The GDC API may be temporarily unavailable — try again later\n")
  cat("  3. If you need controlled-access data, set your GDC token:\n")
  cat("     Sys.setenv(GDC_TOKEN = 'your_token_here')\n")
  stop("GDC query failed. See above for details.")
})

# Extract results table to inspect sample types
results_table <- getResults(query)

cat(sprintf("✓ Query returned %d samples total.\n", nrow(results_table)))
cat("\n  Sample type breakdown:\n")
print(table(results_table$sample_type))
cat("\n")

# ── 4. Subset samples for quick-test mode ────────────────────────────────────

if (RUN_MODE == "quick") {
  cat(sprintf("Step 2/5: Subsetting for quick-test mode (%d tumor + up to %d normal)...\n",
              QUICK_TUMOR_N, QUICK_NORMAL_N))
  
  # Separate tumor and normal barcodes
  tumor_barcodes  <- results_table %>%
    filter(sample_type == "Primary Tumor") %>%
    pull(cases)
  
  normal_barcodes <- results_table %>%
    filter(sample_type == "Solid Tissue Normal") %>%
    pull(cases)
  
  cat(sprintf("  Available: %d tumor, %d normal\n",
              length(tumor_barcodes), length(normal_barcodes)))
  
  # Subset with random sampling (reproducible via set.seed)
  set.seed(123)
  selected_tumor  <- sample(tumor_barcodes, min(QUICK_TUMOR_N, length(tumor_barcodes)))
  selected_normal <- sample(normal_barcodes, min(QUICK_NORMAL_N, length(normal_barcodes)))
  
  selected_cases <- c(selected_tumor, selected_normal)
  
  cat(sprintf("  Selected:  %d tumor, %d normal\n",
              length(selected_tumor), length(selected_normal)))
  
  # Re-query with selected barcodes
  query <- GDCquery(
    project            = "TCGA-BRCA",
    data.category      = "Transcriptome Profiling",
    data.type          = "Gene Expression Quantification",
    workflow.type      = "STAR - Counts",
    sample.type        = c("Primary Tumor", "Solid Tissue Normal"),
    barcode            = selected_cases
  )
  
  results_table <- getResults(query)
  cat(sprintf("✓ Subset query returned %d samples.\n\n", nrow(results_table)))
} else {
  cat("Step 2/5: Full mode — using all available samples.\n\n")
}

# ── 5. Download data from GDC ────────────────────────────────────────────────

cat("Step 3/5: Downloading data from GDC (this may take a while)...\n")
cat(sprintf("  Download directory: %s\n", GDC_DIR))

download_success <- FALSE
max_retries <- 3

for (attempt in 1:max_retries) {
  tryCatch({
    GDCdownload(
      query,
      method    = "api",
      directory = GDC_DIR,
      files.per.chunk = 10
    )
    download_success <- TRUE
    break
  }, error = function(e) {
    if (grepl("already downloaded", e$message, ignore.case = TRUE) ||
        grepl("already been downloaded", e$message, ignore.case = TRUE)) {
      cat("✓ Data files already downloaded locally.\n")
      download_success <<- TRUE
    } else {
      cat(sprintf("\n⚠ Download attempt %d/%d failed: %s\n",
                  attempt, max_retries, e$message))
      if (attempt < max_retries) {
        cat("  Retrying in 10 seconds...\n")
        Sys.sleep(10)
      }
    }
  })
  if (download_success) break
}

if (!download_success) {
  stop("Failed to download data after ", max_retries, " attempts. ",
       "Check your internet connection and try again.")
}

cat("✓ Download complete.\n\n")

# ── 6. Prepare the SummarizedExperiment ──────────────────────────────────────

cat("Step 4/5: Preparing SummarizedExperiment and extracting count matrix...\n")

se <- GDCprepare(
  query,
  directory = GDC_DIR,
  save = FALSE,
  summarizedExperiment = TRUE
)

# Extract count matrix (unstranded counts)
counts_raw <- assay(se, "unstranded")

# Extract sample metadata
col_data <- as.data.frame(colData(se))

# Add a clean condition column for DESeq2
col_data$condition <- ifelse(
  col_data$sample_type %in% c("Primary Tumor", "Primary solid Tumor"),
  "Tumor",
  "Normal"
)
col_data$condition <- factor(col_data$condition, levels = c("Normal", "Tumor"))

cat(sprintf("  Count matrix dimensions: %d genes × %d samples\n",
            nrow(counts_raw), ncol(counts_raw)))
cat(sprintf("  Metadata columns: %d\n", ncol(col_data)))

# ── 7. Basic QC & verification ───────────────────────────────────────────────

cat("\nStep 5/5: Verification & QC...\n\n")

cat("╔══════════════════════════════════════════════╗\n")
cat("║         DATA VERIFICATION SUMMARY            ║\n")
cat("╠══════════════════════════════════════════════╣\n")

n_tumor  <- sum(col_data$condition == "Tumor")
n_normal <- sum(col_data$condition == "Normal")
n_genes  <- nrow(counts_raw)

cat(sprintf("║  Run mode:       %-25s ║\n", toupper(RUN_MODE)))
cat(sprintf("║  Tumor samples:  %-25d ║\n", n_tumor))
cat(sprintf("║  Normal samples: %-25d ║\n", n_normal))
cat(sprintf("║  Total genes:    %-25d ║\n", n_genes))
cat("╚══════════════════════════════════════════════╝\n\n")

# Check for known breast cancer genes
known_genes <- c("ESR1", "ERBB2", "TP53", "MKI67", "TOP2A", "BRCA1", "BRCA2",
                 "CDH1", "GATA3", "FOXA1", "PGR", "MMP9")

# Gene names are stored in rowData
gene_info <- as.data.frame(rowData(se))
gene_names <- gene_info$gene_name

found_genes <- known_genes[known_genes %in% gene_names]
missing_genes <- known_genes[!known_genes %in% gene_names]

cat(sprintf("Known breast cancer genes found: %d/%d\n", 
            length(found_genes), length(known_genes)))
if (length(found_genes) > 0) {
  cat(sprintf("  ✓ Found: %s\n", paste(found_genes, collapse = ", ")))
}
if (length(missing_genes) > 0) {
  cat(sprintf("  ✗ Missing: %s\n", paste(missing_genes, collapse = ", ")))
}

# Check for non-zero expression of key genes
cat("\nMean expression of key genes across all samples:\n")
for (gene in found_genes[1:min(5, length(found_genes))]) {
  gene_idx <- which(gene_names == gene)[1]
  if (!is.na(gene_idx)) {
    mean_expr <- mean(counts_raw[gene_idx, ], na.rm = TRUE)
    cat(sprintf("  %s: %.1f mean counts\n", gene, mean_expr))
  }
}

# ── 8. Save processed data ───────────────────────────────────────────────────

cat("\nSaving processed data...\n")

# Save count matrix
saveRDS(counts_raw, COUNTS_FILE)
cat(sprintf("  ✓ Counts matrix: %s (%.1f MB)\n",
            COUNTS_FILE, file.size(COUNTS_FILE) / 1e6))

# Save column data (metadata)
saveRDS(col_data, COLDATA_FILE)
cat(sprintf("  ✓ Sample metadata: %s (%.1f MB)\n",
            COLDATA_FILE, file.size(COLDATA_FILE) / 1e6))

# Save full SE object
saveRDS(se, SE_FILE)
cat(sprintf("  ✓ SummarizedExperiment: %s (%.1f MB)\n",
            SE_FILE, file.size(SE_FILE) / 1e6))

# ── Done ──────────────────────────────────────────────────────────────────────

cat("\n============================================================\n")
cat("  ✓ Data download and preparation complete!\n")
cat(sprintf("  Run mode: %s\n", toupper(RUN_MODE)))
cat(sprintf("  %d tumor + %d normal samples, %d genes\n", n_tumor, n_normal, n_genes))
cat("  Next: Run scripts/02_DESeq2_analysis.R\n")
cat("============================================================\n")

# Session info for reproducibility
cat("\n--- Session Info ---\n")
sessionInfo()
