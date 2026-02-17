#!/usr/bin/env Rscript
# =============================================================================
# Setup renv for reproducible package management
# =============================================================================
#
# Run this script ONCE to initialize renv and install all required packages.
# After running, use renv::snapshot() to create/update the renv.lock file.
#
# Usage:
#   Rscript setup_renv.R
#
# =============================================================================

cat("============================================================\n")
cat("  Setting up renv for TCGA-BRCA DE Analysis\n")
cat("============================================================\n\n")

# ── 1. Initialize renv ───────────────────────────────────────────────────────

if (!requireNamespace("renv", quietly = TRUE)) {
  cat("Installing renv...\n")
  install.packages("renv", repos = "https://cloud.r-project.org")
}

# Initialize renv (if not already initialized)
if (!file.exists("renv.lock")) {
  cat("Initializing renv project...\n")
  renv::init(bare = TRUE)
} else {
  cat("renv already initialized (renv.lock exists).\n")
  cat("To restore packages: renv::restore()\n\n")
}

# ── 2. Install BiocManager ───────────────────────────────────────────────────

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# ── 3. Install all required packages ─────────────────────────────────────────

cat("\nInstalling required packages...\n")
cat("This may take 15-30 minutes on first run.\n\n")

# CRAN packages
cran_packages <- c(
  "dplyr", "ggplot2", "tibble", "readr", "stringr",
  "pheatmap", "RColorBrewer", "viridis", "ggrepel",
  "rmarkdown", "knitr"
)

# Bioconductor packages
bioc_packages <- c(
  "TCGAbiolinks",
  "DESeq2",
  "SummarizedExperiment",
  "clusterProfiler",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "EnhancedVolcano",
  "enrichplot",
  "DOSE",
  "maftools",
  "apeglm"
)

# Survival packages
survival_packages <- c("survival", "survminer")

# Install CRAN packages
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing %s...\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  } else {
    cat(sprintf("  ✓ %s already installed\n", pkg))
  }
}

# Install Bioconductor packages
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing %s (Bioconductor)...\n", pkg))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  } else {
    cat(sprintf("  ✓ %s already installed\n", pkg))
  }
}

# Install survival packages
for (pkg in survival_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing %s...\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  } else {
    cat(sprintf("  ✓ %s already installed\n", pkg))
  }
}

# ── 4. Snapshot ───────────────────────────────────────────────────────────────

cat("\nCreating renv snapshot (renv.lock)...\n")
renv::snapshot(prompt = FALSE)

cat("\n============================================================\n")
cat("  ✓ renv setup complete!\n")
cat("  All packages installed and locked in renv.lock.\n")
cat("  \n")
cat("  To reproduce this environment elsewhere:\n")
cat("    renv::restore()\n")
cat("  \n")
cat("  Next: Run the analysis pipeline:\n")
cat("    Rscript scripts/01_download_and_prepare.R\n")
cat("============================================================\n")
