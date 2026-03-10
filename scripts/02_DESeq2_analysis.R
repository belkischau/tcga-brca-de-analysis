#!/usr/bin/env Rscript
# =============================================================================
# Script 02: DESeq2 Differential Expression Analysis
# =============================================================================
#
# Description:
#   Performs differential gene expression analysis using DESeq2 on the
#   count matrix prepared by script 01. Compares Primary Tumor vs
#   Solid Tissue Normal.
#
# Input:
#   data/counts_matrix_<mode>.rds
#   data/col_data_<mode>.rds
#   data/se_object_<mode>.rds
#
# Output:
#   results/DEGs_full.csv         — all genes with DE statistics
#   results/DEGs_significant.csv  — filtered DEGs (|log2FC|>1 & padj<0.05)
#   results/deseq2_dds.rds        — fitted DESeqDataSet object
#   results/deseq2_results.rds    — DESeq2 results object
#   results/vst_counts.rds        — variance-stabilized counts (for viz)
#
# Usage:
#   Rscript scripts/02_DESeq2_analysis.R
#   Rscript scripts/02_DESeq2_analysis.R full
#
# Author: Belkis Chau
# =============================================================================

# ── 0. Configuration ─────────────────────────────────────────────────────────

set.seed(123)

args <- commandArgs(trailingOnly = TRUE)
RUN_MODE <- ifelse(length(args) > 0 && tolower(args[1]) == "full", "full", "quick")
suffix <- ifelse(RUN_MODE == "full", "full", "quick")

cat("============================================================\n")
cat("  DESeq2 Differential Expression Analysis\n")
cat(sprintf("  Mode: %s\n", toupper(RUN_MODE)))
cat(sprintf("  Date: %s\n", Sys.time()))
cat("============================================================\n\n")

PROJECT_DIR <- getwd()
DATA_DIR    <- file.path(PROJECT_DIR, "data")
RESULTS_DIR <- file.path(PROJECT_DIR, "results")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 1. Load packages ─────────────────────────────────────────────────────────

cat("Loading packages...\n")

required_packages <- c("DESeq2", "SummarizedExperiment", "dplyr", "tibble",
                        "readr", "stringr", "org.Hs.eg.db", "AnnotationDbi")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' not installed. Run: BiocManager::install('%s')", pkg, pkg))
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

cat("✓ Packages loaded.\n\n")

# ── 2. Load prepared data ────────────────────────────────────────────────────

cat("Loading prepared data...\n")

COUNTS_FILE <- file.path(DATA_DIR, paste0("counts_matrix_", suffix, ".rds"))
COLDATA_FILE <- file.path(DATA_DIR, paste0("col_data_", suffix, ".rds"))
SE_FILE <- file.path(DATA_DIR, paste0("se_object_", suffix, ".rds"))

if (!file.exists(COUNTS_FILE) || !file.exists(COLDATA_FILE)) {
  stop("Data files not found! Run scripts/01_download_and_prepare.R first.\n",
       "  Expected: ", COUNTS_FILE, "\n",
       "  Expected: ", COLDATA_FILE)
}

counts_raw <- readRDS(COUNTS_FILE)
col_data   <- readRDS(COLDATA_FILE)
se         <- readRDS(SE_FILE)

# Extract gene annotation from SE object
gene_info <- as.data.frame(rowData(se))

cat(sprintf("  Loaded count matrix: %d genes × %d samples\n",
            nrow(counts_raw), ncol(counts_raw)))
cat(sprintf("  Tumor: %d | Normal: %d\n",
            sum(col_data$condition == "Tumor"),
            sum(col_data$condition == "Normal")))
cat("✓ Data loaded.\n\n")

# ── 3. Sanity checks ─────────────────────────────────────────────────────────

cat("Running sanity checks...\n")

# Check that column names of count matrix match rownames of col_data
stopifnot("Column names of counts must match row names of col_data" =
            all(colnames(counts_raw) %in% rownames(col_data)))

# Align count matrix columns with col_data rows
counts_raw <- counts_raw[, rownames(col_data)]

# Ensure condition is a factor with Normal as reference
col_data$condition <- factor(col_data$condition, levels = c("Normal", "Tumor"))

# Check for minimum sample sizes
n_tumor  <- sum(col_data$condition == "Tumor")
n_normal <- sum(col_data$condition == "Normal")

if (n_normal < 3) {
  stop("Need at least 3 normal samples for DESeq2. Found: ", n_normal)
}
if (n_tumor < 3) {
  stop("Need at least 3 tumor samples for DESeq2. Found: ", n_tumor)
}

cat(sprintf("  ✓ %d tumor and %d normal samples — sufficient for DE analysis\n",
            n_tumor, n_normal))

# Remove genes with zero counts across all samples
nonzero_genes <- rowSums(counts_raw) > 0
cat(sprintf("  ✓ Removing %d genes with zero total counts (%d remain)\n",
            sum(!nonzero_genes), sum(nonzero_genes)))
counts_filtered <- counts_raw[nonzero_genes, ]

cat("✓ Sanity checks passed.\n\n")

# ── 4. Create DESeqDataSet ───────────────────────────────────────────────────

cat("Creating DESeqDataSet...\n")

# DESeq2 requires integer counts
storage.mode(counts_filtered) <- "integer"

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData   = col_data,
  design    = ~ condition
)

# Pre-filtering: keep genes with at least 10 counts across samples
# This speeds up computation and reduces multiple-testing burden
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]

cat(sprintf("  ✓ After pre-filtering: %d genes retained\n", nrow(dds)))
cat("✓ DESeqDataSet created.\n\n")

# ── 5. Run DESeq2 ────────────────────────────────────────────────────────────

cat("Running DESeq2 (this may take a few minutes)...\n")
t_start <- Sys.time()

dds <- DESeq(dds)

t_elapsed <- difftime(Sys.time(), t_start, units = "mins")
cat(sprintf("✓ DESeq2 complete in %.1f minutes.\n\n", as.numeric(t_elapsed)))

# ── 6. Extract results ───────────────────────────────────────────────────────

cat("Extracting results (Tumor vs Normal)...\n")

# Get results with shrunken log2 fold changes (apeglm shrinkage)
# Note: apeglm requires the coefficient name
resultsNames(dds)

res <- results(dds,
               contrast = c("condition", "Tumor", "Normal"),
               alpha = 0.05)

# Use lfcShrink for better log2FC estimates (ashr method for broad compatibility)
res_shrunk <- tryCatch({
  lfcShrink(dds,
            coef = "condition_Tumor_vs_Normal",
            type = "apeglm",
            res  = res)
}, error = function(e) {
  cat("  Note: apeglm shrinkage failed, falling back to 'normal' shrinkage.\n")
  lfcShrink(dds,
            coef = "condition_Tumor_vs_Normal",
            type = "normal",
            res  = res)
})

cat(sprintf("  Total genes tested: %d\n", nrow(res_shrunk)))
cat(sprintf("  Genes with padj < 0.05: %d\n",
            sum(res_shrunk$padj < 0.05, na.rm = TRUE)))
cat(sprintf("  DEGs (|log2FC| > 1 & padj < 0.05): %d\n",
            sum(abs(res_shrunk$log2FoldChange) > 1 & res_shrunk$padj < 0.05,
                na.rm = TRUE)))

# ── 7. Annotate results with gene symbols ────────────────────────────────────

cat("\nAnnotating results with gene symbols...\n")

# Convert results to data frame
res_df <- as.data.frame(res_shrunk) %>%
  tibble::rownames_to_column("gene_id")

# Map ENSEMBL IDs to gene symbols using the SE rowData
# gene_info has gene_id and gene_name columns from TCGAbiolinks
gene_map <- gene_info[, c("gene_id", "gene_name", "gene_type")]
# Ensure gene_id is character
gene_map$gene_id <- as.character(gene_map$gene_id)

# The rownames of DESeq2 results may be ENSEMBL IDs (with or without version)
# Let's check the format
sample_id <- res_df$gene_id[1]
cat(sprintf("  Sample gene ID format: %s\n", sample_id))

# Merge gene annotations
if ("gene_id" %in% colnames(gene_map)) {
  res_df <- res_df %>%
    left_join(gene_map, by = "gene_id") %>%
    distinct(gene_id, .keep_all = TRUE)
}

# If gene_name is still missing, try to map from org.Hs.eg.db
if (!"gene_name" %in% colnames(res_df) || all(is.na(res_df$gene_name))) {
  cat("  Mapping gene symbols from org.Hs.eg.db...\n")
  # Strip version numbers from ENSEMBL IDs
  ensembl_ids <- gsub("\\.\\d+$", "", res_df$gene_id)
  
  symbol_map <- tryCatch({
    AnnotationDbi::select(org.Hs.eg.db,
                          keys = ensembl_ids,
                          columns = "SYMBOL",
                          keytype = "ENSEMBL")
  }, error = function(e) {
    cat("  Warning: Gene symbol mapping failed.\n")
    data.frame(ENSEMBL = ensembl_ids, SYMBOL = NA_character_)
  })
  
  # Handle duplicates
  symbol_map <- symbol_map[!duplicated(symbol_map$ENSEMBL), ]
  res_df$ensembl_clean <- gsub("\\.\\d+$", "", res_df$gene_id)
  res_df <- res_df %>%
    left_join(symbol_map, by = c("ensembl_clean" = "ENSEMBL")) %>%
    mutate(gene_name = coalesce(gene_name, SYMBOL)) %>%
    dplyr::select(-SYMBOL, -ensembl_clean)
}

cat(sprintf("  ✓ %d / %d genes have gene symbols\n",
            sum(!is.na(res_df$gene_name)), nrow(res_df)))

# ── 8. Sort and classify DEGs ────────────────────────────────────────────────

cat("\nClassifying differentially expressed genes...\n")

# Define significance thresholds
LFC_THRESHOLD <- 1
PADJ_THRESHOLD <- 0.05

res_df <- res_df %>%
  mutate(
    DE_status = case_when(
      is.na(padj) ~ "Not Tested",
      padj >= PADJ_THRESHOLD ~ "Not Significant",
      abs(log2FoldChange) < LFC_THRESHOLD ~ "Significant (small LFC)",
      log2FoldChange > LFC_THRESHOLD ~ "Upregulated",
      log2FoldChange < -LFC_THRESHOLD ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  ) %>%
  arrange(padj)

# Summary table
de_summary <- table(res_df$DE_status)
cat("\n  DE Status Summary:\n")
print(de_summary)

n_up   <- sum(res_df$DE_status == "Upregulated", na.rm = TRUE)
n_down <- sum(res_df$DE_status == "Downregulated", na.rm = TRUE)
cat(sprintf("\n  Upregulated in Tumor: %d\n", n_up))
cat(sprintf("  Downregulated in Tumor: %d\n", n_down))

# ── 9. Verify known breast cancer genes ──────────────────────────────────────

cat("\n╔══════════════════════════════════════════════╗\n")
cat("║    KNOWN BREAST CANCER GENE VERIFICATION     ║\n")
cat("╠══════════════════════════════════════════════╣\n")

known_genes <- data.frame(
  gene = c("ESR1", "ERBB2", "TP53", "MKI67", "TOP2A", 
           "BRCA1", "CDH1", "GATA3", "FOXA1", "PGR"),
  expected_direction = c("variable", "up_in_tumor", "variable", "up_in_tumor", 
                          "up_in_tumor", "variable", "variable", "variable",
                          "variable", "variable"),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(known_genes)) {
  gene <- known_genes$gene[i]
  idx <- which(res_df$gene_name == gene)
  
  if (length(idx) > 0) {
    row <- res_df[idx[1], ]
    direction <- ifelse(row$log2FoldChange > 0, "UP", "DOWN")
    sig <- ifelse(!is.na(row$padj) && row$padj < PADJ_THRESHOLD, "***", "ns")
    cat(sprintf("║  %-8s  log2FC: %+7.2f  padj: %.2e  %s %s ║\n",
                gene,
                round(row$log2FoldChange, 2),
                ifelse(is.na(row$padj), NA, row$padj),
                direction, sig))
  } else {
    cat(sprintf("║  %-8s  NOT FOUND in results                  ║\n", gene))
  }
}
cat("╚══════════════════════════════════════════════╝\n")

# ── 10. Print top DEGs ───────────────────────────────────────────────────────

cat("\n── Top 10 Upregulated Genes (Tumor vs Normal) ──\n")
top_up <- res_df %>%
  filter(DE_status == "Upregulated") %>%
  arrange(desc(log2FoldChange)) %>%
  head(10) %>%
  dplyr::select(gene_name, log2FoldChange, padj, baseMean)

print(top_up, row.names = FALSE)

cat("\n── Top 10 Downregulated Genes (Tumor vs Normal) ──\n")
top_down <- res_df %>%
  filter(DE_status == "Downregulated") %>%
  arrange(log2FoldChange) %>%
  head(10) %>%
  dplyr::select(gene_name, log2FoldChange, padj, baseMean)

print(top_down, row.names = FALSE)

# ── 11. Compute variance-stabilized counts for visualization ─────────────────

cat("\nComputing variance-stabilized transformation (for visualization)...\n")
vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)

cat(sprintf("  ✓ VST matrix: %d genes × %d samples\n",
            nrow(vst_mat), ncol(vst_mat)))

# ── 12. Save results ─────────────────────────────────────────────────────────

cat("\nSaving results...\n")

# Full results
readr::write_csv(res_df, file.path(RESULTS_DIR, "DEGs_full.csv"))
cat(sprintf("  ✓ results/DEGs_full.csv (%d rows)\n", nrow(res_df)))

# Significant DEGs only
sig_df <- res_df %>%
  filter(DE_status %in% c("Upregulated", "Downregulated"))

readr::write_csv(sig_df, file.path(RESULTS_DIR, "DEGs_significant.csv"))
cat(sprintf("  ✓ results/DEGs_significant.csv (%d rows)\n", nrow(sig_df)))

# R objects for downstream scripts
saveRDS(dds, file.path(RESULTS_DIR, "deseq2_dds.rds"))
saveRDS(res_shrunk, file.path(RESULTS_DIR, "deseq2_results.rds"))
saveRDS(res_df, file.path(RESULTS_DIR, "deseq2_results_annotated.rds"))
saveRDS(vst_mat, file.path(RESULTS_DIR, "vst_counts.rds"))

cat("  ✓ R objects saved to results/\n")

# ── Done ──────────────────────────────────────────────────────────────────────

cat("\n============================================================\n")
cat("  ✓ DESeq2 analysis complete!\n")
cat(sprintf("  Total genes tested: %d\n", nrow(res_df)))
cat(sprintf("  Significant DEGs: %d (%d up, %d down)\n",
            n_up + n_down, n_up, n_down))
cat("  Next: Run scripts/03_pathway_enrichment.R\n")
cat("============================================================\n")

cat("\n--- Session Info ---\n")
sessionInfo()
