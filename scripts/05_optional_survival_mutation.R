#!/usr/bin/env Rscript
# =============================================================================
# Script 05: Optional Extensions — Survival Analysis & Mutation Landscape
# =============================================================================
#
# Description:
#   1. Kaplan-Meier survival analysis for the top 3 DEGs
#   2. maftools oncoplot for TCGA-BRCA somatic mutations
#   3. Notes on Nextflow integration
#
# Input:
#   results/deseq2_results_annotated.rds
#   data/se_object_<mode>.rds
#   data/col_data_<mode>.rds
#
# Output:
#   results/KM_survival_<gene>.png    — Kaplan-Meier plots
#   results/oncoplot_BRCA.png         — Oncoplot of top mutated genes
#
# Usage:
#   Rscript scripts/05_optional_survival_mutation.R
#
# Note:
#   This script requires internet access for maftools to download
#   TCGA-BRCA MAF data. If maftools is not installed, it will skip
#   the mutation analysis gracefully.
#
# Author: Belkis Chau
# =============================================================================

# ── 0. Configuration ─────────────────────────────────────────────────────────

set.seed(123)

cat("============================================================\n")
cat("  Optional Extensions: Survival Analysis & Mutation Landscape\n")
cat(sprintf("  Date: %s\n", Sys.time()))
cat("============================================================\n\n")

PROJECT_DIR <- getwd()
RESULTS_DIR <- file.path(PROJECT_DIR, "results")
DATA_DIR    <- file.path(PROJECT_DIR, "data")

# ══════════════════════════════════════════════════════════════════════════════
# PART 1: KAPLAN-MEIER SURVIVAL ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════

cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  PART 1: Kaplan-Meier Survival Analysis\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# ── 1.1 Load packages ────────────────────────────────────────────────────────

survival_pkgs <- c("survival", "survminer", "ggplot2", "dplyr",
                     "SummarizedExperiment")
survival_available <- TRUE

for (pkg in survival_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  ⚠ Package '%s' not installed. Skipping survival analysis.\n", pkg))
    survival_available <- FALSE
    break
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

if (survival_available) {
  
  # ── 1.2 Load data ──────────────────────────────────────────────────────────
  
  cat("Loading data for survival analysis...\n")
  
  res_df <- readRDS(file.path(RESULTS_DIR, "deseq2_results_annotated.rds"))
  
  # Load SE object
  se_file <- if (file.exists(file.path(DATA_DIR, "se_object_quick.rds"))) {
    file.path(DATA_DIR, "se_object_quick.rds")
  } else {
    file.path(DATA_DIR, "se_object_full.rds")
  }
  se <- readRDS(se_file)
  
  # Load VST counts
  vst_mat <- readRDS(file.path(RESULTS_DIR, "vst_counts.rds"))
  
  col_data <- as.data.frame(colData(se))
  gene_info <- as.data.frame(rowData(se))
  
  # ── 1.3 Select top 3 DEGs for survival analysis ────────────────────────────
  
  # Only use tumor samples for survival analysis (normal tissue has no clinical follow-up)
  tumor_samples <- rownames(col_data)[col_data$sample_type %in% c("Primary Tumor", "Primary solid Tumor")]
  tumor_samples_in_vst <- intersect(tumor_samples, colnames(vst_mat))
  
  cat(sprintf("  Tumor samples with expression data: %d\n", length(tumor_samples_in_vst)))
  
  # Select top 3 most significant DEGs with gene symbols
  top3_genes <- res_df %>%
    filter(DE_status %in% c("Upregulated", "Downregulated")) %>%
    filter(!is.na(gene_name) & gene_name != "") %>%
    arrange(padj) %>%
    head(3)
  
  cat(sprintf("  Top 3 DEGs for survival: %s\n",
              paste(top3_genes$gene_name, collapse = ", ")))
  
  # ── 1.4 Prepare survival data ──────────────────────────────────────────────
  
  # Extract clinical survival information
  # TCGAbiolinks stores these in colData columns
  survival_cols <- c("days_to_death", "days_to_last_follow_up", "vital_status")
  
  has_survival <- all(survival_cols %in% colnames(col_data))
  
  if (has_survival) {
    cat("  ✓ Survival data found in colData.\n\n")
    
    # Prepare survival dataframe for tumor samples
    surv_data <- col_data[tumor_samples_in_vst, ]
    
    # Create overall survival time and event
    surv_data$OS_time <- as.numeric(ifelse(
      surv_data$vital_status == "Dead",
      surv_data$days_to_death,
      surv_data$days_to_last_follow_up
    ))
    
    surv_data$OS_event <- ifelse(surv_data$vital_status == "Dead", 1, 0)
    
    # Remove samples with missing survival data
    surv_data <- surv_data[!is.na(surv_data$OS_time) & surv_data$OS_time > 0, ]
    
    cat(sprintf("  Samples with valid survival data: %d\n", nrow(surv_data)))
    cat(sprintf("  Events (deaths): %d\n", sum(surv_data$OS_event)))
    
    if (nrow(surv_data) >= 20 && sum(surv_data$OS_event) >= 5) {
      
      # ── 1.5 Generate KM plots for each top DEG ─────────────────────────────
      
      for (i in 1:nrow(top3_genes)) {
        gene_name <- top3_genes$gene_name[i]
        gene_id   <- top3_genes$gene_id[i]
        
        cat(sprintf("\n  Generating KM plot for %s...\n", gene_name))
        
        # Check if gene exists in VST matrix
        if (!gene_id %in% rownames(vst_mat)) {
          cat(sprintf("    ⚠ %s not found in VST matrix. Skipping.\n", gene_name))
          next
        }
        
        # Get expression values for tumor samples with survival data
        common_samples <- intersect(rownames(surv_data), colnames(vst_mat))
        expr_vals <- vst_mat[gene_id, common_samples]
        
        # Stratify into High/Low expression groups (median split)
        median_expr <- median(expr_vals, na.rm = TRUE)
        surv_data_gene <- surv_data[common_samples, ]
        surv_data_gene$expression_group <- ifelse(
          expr_vals > median_expr, "High", "Low"
        )
        surv_data_gene$expression_group <- factor(
          surv_data_gene$expression_group, levels = c("Low", "High")
        )
        
        # Fit survival model
        surv_fit <- survfit(
          Surv(OS_time, OS_event) ~ expression_group,
          data = surv_data_gene
        )
        
        # Generate KM plot
        p_km <- ggsurvplot(
          surv_fit,
          data       = surv_data_gene,
          pval       = TRUE,
          pval.method = TRUE,
          risk.table = TRUE,
          risk.table.col = "strata",
          palette    = c("#2166AC", "#B2182B"),
          title      = sprintf("Overall Survival by %s Expression", gene_name),
          xlab       = "Time (days)",
          ylab       = "Survival Probability",
          legend.title = sprintf("%s Expression", gene_name),
          legend.labs = c("Low", "High"),
          ggtheme    = theme_minimal(base_size = 12),
          risk.table.height = 0.25
        )
        
        # Save
        km_file <- file.path(RESULTS_DIR, sprintf("KM_survival_%s.png", gene_name))
        
        png(km_file, width = 8, height = 7, units = "in", res = 300)
        print(p_km)
        dev.off()
        
        cat(sprintf("    ✓ Saved %s\n", basename(km_file)))
      }
      
    } else {
      cat("  ⚠ Insufficient survival data for meaningful KM analysis.\n")
      cat("    (Need ≥20 samples with ≥5 events)\n")
      cat("    This is expected in quick-test mode with limited samples.\n")
    }
    
  } else {
    cat("  ⚠ Survival columns not found in colData.\n")
    cat(sprintf("    Available columns: %s\n",
                paste(head(colnames(col_data), 20), collapse = ", ")))
  }
  
} # end if survival_available

# ══════════════════════════════════════════════════════════════════════════════
# PART 2: MUTATION LANDSCAPE WITH MAFTOOLS
# ══════════════════════════════════════════════════════════════════════════════

cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  PART 2: Mutation Landscape (maftools)\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

maftools_available <- requireNamespace("maftools", quietly = TRUE)

if (maftools_available) {
  suppressPackageStartupMessages(library(maftools))
  
  cat("Downloading TCGA-BRCA MAF data via maftools...\n")
  cat("(This uses pre-compiled TCGA MAF files from the GDC)\n\n")
  
  brca_maf <- tryCatch({
    tcgaLoad(study = "BRCA")
  }, error = function(e) {
    cat(sprintf("  ⚠ maftools::tcgaLoad failed: %s\n", e$message))
    cat("  Trying alternative: read.maf from TCGAbiolinks...\n")
    
    # Alternative: try GDCquery_Maf
    tryCatch({
      if (requireNamespace("TCGAbiolinks", quietly = TRUE)) {
        maf_data <- TCGAbiolinks::GDCquery_Maf(
          tumor = "BRCA",
          pipelines = "mutect2"
        )
        read.maf(maf = maf_data)
      } else {
        NULL
      }
    }, error = function(e2) {
      cat(sprintf("  ⚠ Alternative also failed: %s\n", e2$message))
      NULL
    })
  })
  
  if (!is.null(brca_maf)) {
    cat(sprintf("  ✓ Loaded MAF with %d variants from %d samples\n",
                nrow(brca_maf@data), length(unique(brca_maf@data$Tumor_Sample_Barcode))))
    
    # Summary plot
    cat("  Generating MAF summary plot...\n")
    png(file.path(RESULTS_DIR, "maf_summary.png"),
        width = 12, height = 8, units = "in", res = 300)
    plotmafSummary(maf = brca_maf, rmOutlier = TRUE, addStat = "median")
    dev.off()
    cat("  ✓ Saved results/maf_summary.png\n")
    
    # Oncoplot of top 20 mutated genes
    cat("  Generating oncoplot...\n")
    png(file.path(RESULTS_DIR, "oncoplot_BRCA.png"),
        width = 12, height = 8, units = "in", res = 300)
    oncoplot(maf = brca_maf, top = 20,
             fontSize = 0.6,
             titleFontSize = 1.0,
             legendFontSize = 0.8,
             annotationFontSize = 0.8)
    dev.off()
    cat("  ✓ Saved results/oncoplot_BRCA.png\n")
    
    # Check overlap with our top DEGs
    cat("\n  Checking overlap between top mutated & top DE genes...\n")
    top_mutated <- getGeneSummary(brca_maf)$Hugo_Symbol[1:50]
    
    sig_deg_symbols <- res_df %>%
      filter(DE_status %in% c("Upregulated", "Downregulated")) %>%
      filter(!is.na(gene_name)) %>%
      pull(gene_name)
    
    overlap <- intersect(top_mutated, sig_deg_symbols)
    cat(sprintf("  Overlap (top 50 mutated ∩ significant DEGs): %d genes\n",
                length(overlap)))
    if (length(overlap) > 0) {
      cat(sprintf("    %s\n", paste(head(overlap, 15), collapse = ", ")))
    }
    
  } else {
    cat("  ⚠ Could not load TCGA-BRCA MAF data. Skipping mutation analysis.\n")
  }
  
} else {
  cat("  ⚠ maftools not installed. Skipping mutation analysis.\n")
  cat("  Install with: BiocManager::install('maftools')\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# PART 3: NEXTFLOW INTEGRATION NOTES
# ══════════════════════════════════════════════════════════════════════════════

cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  PART 3: Notes on Nextflow Integration\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

nextflow_notes <- '
## Nextflow Pipeline Integration

This R pipeline can be converted to a Nextflow workflow for scalable,
reproducible execution on HPC clusters or cloud platforms.

### Proposed Nextflow Architecture:

```
workflow {
    // Step 1: Download and prepare data
    DOWNLOAD_DATA(params.run_mode)
    
    // Step 2: DESeq2 analysis (depends on Step 1)
    DESEQ2_ANALYSIS(DOWNLOAD_DATA.out.counts, DOWNLOAD_DATA.out.coldata)
    
    // Step 3: Pathway enrichment (depends on Step 2)
    PATHWAY_ENRICHMENT(DESEQ2_ANALYSIS.out.results)
    
    // Step 4: Visualization (depends on Steps 2 & 3)
    VISUALIZATION(DESEQ2_ANALYSIS.out.results, 
                  DESEQ2_ANALYSIS.out.vst,
                  PATHWAY_ENRICHMENT.out.go_results)
    
    // Step 5: Optional (depends on Step 2)
    SURVIVAL_ANALYSIS(DESEQ2_ANALYSIS.out.results, DOWNLOAD_DATA.out.se)
}
```

### Each process maps to an R script:

| Nextflow Process      | R Script                          | Container        |
|-----------------------|-----------------------------------|------------------|
| DOWNLOAD_DATA         | scripts/01_download_and_prepare.R | bioconductor/... |
| DESEQ2_ANALYSIS       | scripts/02_DESeq2_analysis.R      | bioconductor/... |
| PATHWAY_ENRICHMENT    | scripts/03_pathway_enrichment.R   | bioconductor/... |
| VISUALIZATION         | scripts/04_visualization.R        | bioconductor/... |
| SURVIVAL_ANALYSIS     | scripts/05_optional_survival_mutation.R | bioconductor/... |

### Key benefits of Nextflow wrapping:
1. **Reproducibility** via Docker/Singularity containers
2. **Scalability** — run on local, HPC (SLURM), or cloud (AWS Batch)
3. **Resume** — failed steps can restart without re-running completed ones
4. **Provenance** — full execution trace and reports
'

cat(nextflow_notes)

# Write notes to file for reference
writeLines(nextflow_notes, file.path(RESULTS_DIR, "nextflow_integration_notes.md"))
cat("  ✓ Notes saved to results/nextflow_integration_notes.md\n")

# ── Done ──────────────────────────────────────────────────────────────────────

cat("\n============================================================\n")
cat("  ✓ Optional extensions complete!\n")
cat("  Output files:\n")

opt_files <- list.files(RESULTS_DIR, pattern = "(KM_|oncoplot|maf_|nextflow)",
                         full.names = FALSE)
for (f in opt_files) {
  cat(sprintf("    • %s\n", f))
}

cat("\n  Pipeline complete! Render the report with:\n")
cat("    Rscript -e 'rmarkdown::render(\"reports/TCGA_BRCA_DE_Report.Rmd\")'\n")
cat("============================================================\n")

cat("\n--- Session Info ---\n")
sessionInfo()
