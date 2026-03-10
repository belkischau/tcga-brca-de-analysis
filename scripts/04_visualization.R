#!/usr/bin/env Rscript
# =============================================================================
# Script 04: Publication-Quality Visualizations
# =============================================================================
#
# Description:
#   Generates all publication-quality figures for the TCGA-BRCA DE analysis:
#   - Volcano plot (EnhancedVolcano)
#   - Heatmap of top 50 DEGs (pheatmap)
#   - PCA plot (ggplot2)
#   - MA plot
#   - DE summary bar plot
#   - Sample correlation heatmap
#
# Input:
#   results/deseq2_results_annotated.rds
#   results/deseq2_dds.rds
#   results/vst_counts.rds
#   data/col_data_<mode>.rds
#
# Output:
#   results/volcano_plot.png
#   results/top50_heatmap.png
#   results/PCA_plot.png
#   results/MA_plot.png
#   results/DE_summary_barplot.png
#   results/sample_correlation_heatmap.png
#
# Usage:
#   Rscript scripts/04_visualization.R
#
# Author: Belkis Chau
# =============================================================================

# ── 0. Configuration ─────────────────────────────────────────────────────────

set.seed(123)

cat("============================================================\n")
cat("  Publication-Quality Visualizations\n")
cat(sprintf("  Date: %s\n", Sys.time()))
cat("============================================================\n\n")

PROJECT_DIR <- getwd()
RESULTS_DIR <- file.path(PROJECT_DIR, "results")
DATA_DIR    <- file.path(PROJECT_DIR, "data")

# ── 1. Load packages ─────────────────────────────────────────────────────────

cat("Loading packages...\n")

required_packages <- c("EnhancedVolcano", "pheatmap", "DESeq2", "ggplot2",
                        "dplyr", "tibble", "RColorBrewer", "viridis",
                        "SummarizedExperiment", "ggrepel")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' not installed.", pkg))
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

cat("✓ Packages loaded.\n\n")

# ── 2. Load data ─────────────────────────────────────────────────────────────

cat("Loading data...\n")

res_df  <- readRDS(file.path(RESULTS_DIR, "deseq2_results_annotated.rds"))
dds     <- readRDS(file.path(RESULTS_DIR, "deseq2_dds.rds"))
vst_mat <- readRDS(file.path(RESULTS_DIR, "vst_counts.rds"))

# Try both quick and full mode col_data
col_data_file <- if (file.exists(file.path(DATA_DIR, "col_data_quick.rds"))) {
  file.path(DATA_DIR, "col_data_quick.rds")
} else {
  file.path(DATA_DIR, "col_data_full.rds")
}
col_data <- readRDS(col_data_file)

cat(sprintf("  DEGs: %d genes\n", nrow(res_df)))
cat(sprintf("  VST matrix: %d genes × %d samples\n", nrow(vst_mat), ncol(vst_mat)))
cat("✓ Data loaded.\n\n")

# ── 3. Volcano Plot ──────────────────────────────────────────────────────────

cat("Generating volcano plot...\n")

# Prepare data for EnhancedVolcano
volcano_data <- res_df %>%
  filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
  mutate(gene_label = ifelse(is.na(gene_name) | gene_name == "", gene_id, gene_name))

# Highlight known breast cancer genes
bc_genes_highlight <- c("ESR1", "ERBB2", "TP53", "MKI67", "TOP2A",
                         "BRCA1", "BRCA2", "CDH1", "GATA3", "FOXA1", "PGR")

# Create named vector for EnhancedVolcano — must use unique rownames
volcano_data$gene_label <- make.unique(volcano_data$gene_label)
rownames(volcano_data) <- volcano_data$gene_label

p_volcano <- EnhancedVolcano(
  volcano_data,
  lab             = volcano_data$gene_label,
  x               = "log2FoldChange",
  y               = "padj",
  title           = "TCGA-BRCA: Tumor vs Normal",
  subtitle        = "Differential Gene Expression (DESeq2)",
  pCutoff         = 0.05,
  FCcutoff        = 1.0,
  pointSize       = 1.5,
  labSize         = 3.0,
  labFace         = "italic",
  selectLab       = bc_genes_highlight,
  drawConnectors  = TRUE,
  widthConnectors = 0.5,
  colConnectors   = "grey50",
  maxoverlapsConnectors = 20,
  col             = c("grey80", "#2166AC", "#B2182B", "#D6604D"),
  legendPosition  = "right",
  legendLabSize   = 10,
  xlim            = c(-10, 10),
  ylim            = c(0, -log10(min(volcano_data$padj[volcano_data$padj > 0],
                                     na.rm = TRUE)) + 5)
)

ggsave(file.path(RESULTS_DIR, "volcano_plot.png"), p_volcano,
       width = 12, height = 8, dpi = 300)
cat("  ✓ Saved results/volcano_plot.png\n\n")

# ── 4. Heatmap of Top 50 DEGs ────────────────────────────────────────────────

cat("Generating heatmap of top 50 DEGs...\n")

# Select top 50 DEGs by adjusted p-value
top50 <- res_df %>%
  filter(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1) %>%
  arrange(padj) %>%
  head(50)

# Get VST values for these genes
top50_ids <- top50$gene_id
top50_ids_available <- top50_ids[top50_ids %in% rownames(vst_mat)]

if (length(top50_ids_available) > 0) {
  heatmap_mat <- vst_mat[top50_ids_available, ]
  
  # Use gene names as row labels
  gene_labels <- top50$gene_name[match(top50_ids_available, top50$gene_id)]
  gene_labels[is.na(gene_labels) | gene_labels == ""] <- top50_ids_available[is.na(gene_labels) | gene_labels == ""]
  rownames(heatmap_mat) <- make.unique(gene_labels)
  
  # Scale rows (z-score) for better visualization
  heatmap_mat_scaled <- t(scale(t(heatmap_mat)))
  
  # Annotation for columns (tumor vs normal)
  # Align sample order
  common_samples <- intersect(colnames(heatmap_mat), rownames(col_data))
  heatmap_mat_scaled <- heatmap_mat_scaled[, common_samples]
  
  annotation_col <- data.frame(
    Condition = col_data[common_samples, "condition"],
    row.names = common_samples
  )
  
  # Color scheme
  ann_colors <- list(
    Condition = c(Normal = "#2166AC", Tumor = "#B2182B")
  )
  
  # Clip extreme values for better color range
  heatmap_mat_clipped <- heatmap_mat_scaled
  heatmap_mat_clipped[heatmap_mat_clipped > 3] <- 3
  heatmap_mat_clipped[heatmap_mat_clipped < -3] <- -3
  
  png(file.path(RESULTS_DIR, "top50_heatmap.png"),
      width = 14, height = 12, units = "in", res = 300)
  
  pheatmap(
    heatmap_mat_clipped,
    annotation_col    = annotation_col,
    annotation_colors = ann_colors,
    show_colnames     = FALSE,
    show_rownames     = TRUE,
    fontsize_row      = 7,
    fontsize          = 10,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    clustering_method = "ward.D2",
    color             = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    main              = "Top 50 Differentially Expressed Genes\n(TCGA-BRCA: Tumor vs Normal)",
    border_color      = NA
  )
  
  dev.off()
  cat("  ✓ Saved results/top50_heatmap.png\n\n")
} else {
  cat("  ⚠ No matching genes found in VST matrix for heatmap.\n\n")
}

# ── 5. PCA Plot ──────────────────────────────────────────────────────────────

cat("Generating PCA plot...\n")

# Perform PCA on top 500 most variable genes
rv <- apply(vst_mat, 1, var)
top_var_genes <- head(order(rv, decreasing = TRUE), 500)
pca_mat <- t(vst_mat[top_var_genes, ])

pca_result <- prcomp(pca_mat, center = TRUE, scale. = FALSE)

# Variance explained
var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

# Build PCA data frame
pca_df <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Sample = rownames(pca_result$x)
)

# Match to condition
pca_df$Condition <- col_data[match(pca_df$Sample, rownames(col_data)), "condition"]

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c(Normal = "#2166AC", Tumor = "#B2182B")) +
  labs(
    title    = "PCA of TCGA-BRCA Samples",
    subtitle = "Top 500 most variable genes (VST-normalized)",
    x        = sprintf("PC1 (%.1f%% variance)", var_explained[1]),
    y        = sprintf("PC2 (%.1f%% variance)", var_explained[2]),
    color    = "Sample Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(RESULTS_DIR, "PCA_plot.png"), p_pca,
       width = 8, height = 6, dpi = 300)
cat("  ✓ Saved results/PCA_plot.png\n\n")

# ── 6. MA Plot ───────────────────────────────────────────────────────────────

cat("Generating MA plot...\n")

ma_df <- res_df %>%
  filter(!is.na(padj) & !is.na(log2FoldChange) & !is.na(baseMean)) %>%
  mutate(
    significant = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

p_ma <- ggplot(ma_df, aes(x = log10(baseMean + 1), y = log2FoldChange, color = significant)) +
  geom_point(size = 0.5, alpha = 0.4) +
  scale_color_manual(
    values = c(Up = "#B2182B", Down = "#2166AC", NS = "grey70"),
    name   = "DE Status"
  ) +
  geom_hline(yintercept = c(-1, 0, 1), linetype = c("dashed", "solid", "dashed"),
             color = c("blue", "black", "red"), alpha = 0.5) +
  labs(
    title    = "MA Plot — TCGA-BRCA Tumor vs Normal",
    x        = "log10(Mean Expression + 1)",
    y        = "log2 Fold Change"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

ggsave(file.path(RESULTS_DIR, "MA_plot.png"), p_ma,
       width = 8, height = 6, dpi = 300)
cat("  ✓ Saved results/MA_plot.png\n\n")

# ── 7. DE Summary Bar Plot ───────────────────────────────────────────────────

cat("Generating DE summary bar plot...\n")

de_counts <- res_df %>%
  filter(DE_status %in% c("Upregulated", "Downregulated")) %>%
  count(DE_status) %>%
  mutate(DE_status = factor(DE_status, levels = c("Upregulated", "Downregulated")))

p_bar <- ggplot(de_counts, aes(x = DE_status, y = n, fill = DE_status)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = format(n, big.mark = ",")), vjust = -0.5, size = 5) +
  scale_fill_manual(values = c(Upregulated = "#B2182B", Downregulated = "#2166AC")) +
  labs(
    title = "Number of Differentially Expressed Genes",
    subtitle = "|log2FC| > 1 & adjusted p-value < 0.05",
    x = NULL,
    y = "Number of DEGs"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  ylim(0, max(de_counts$n) * 1.15)

ggsave(file.path(RESULTS_DIR, "DE_summary_barplot.png"), p_bar,
       width = 6, height = 5, dpi = 300)
cat("  ✓ Saved results/DE_summary_barplot.png\n\n")

# ── 8. Sample Correlation Heatmap ────────────────────────────────────────────

cat("Generating sample correlation heatmap...\n")

# Use subset of samples for readability if > 50
sample_subset <- colnames(vst_mat)
if (length(sample_subset) > 60) {
  # Sample equal numbers of tumor and normal
  tumor_samples <- rownames(col_data)[col_data$condition == "Tumor"]
  normal_samples <- rownames(col_data)[col_data$condition == "Normal"]
  
  set.seed(123)
  n_show <- 30
  selected_tumor <- sample(intersect(tumor_samples, colnames(vst_mat)),
                            min(n_show, length(intersect(tumor_samples, colnames(vst_mat)))))
  selected_normal <- sample(intersect(normal_samples, colnames(vst_mat)),
                             min(n_show, length(intersect(normal_samples, colnames(vst_mat)))))
  sample_subset <- c(selected_tumor, selected_normal)
}

# Compute sample-sample correlation on top variable genes
cor_mat <- cor(vst_mat[top_var_genes, sample_subset], method = "spearman")

# Annotation
annotation_col_cor <- data.frame(
  Condition = col_data[sample_subset, "condition"],
  row.names = sample_subset
)

ann_colors_cor <- list(
  Condition = c(Normal = "#2166AC", Tumor = "#B2182B")
)

png(file.path(RESULTS_DIR, "sample_correlation_heatmap.png"),
    width = 10, height = 9, units = "in", res = 300)

pheatmap(
  cor_mat,
  annotation_col    = annotation_col_cor,
  annotation_row    = annotation_col_cor,
  annotation_colors = ann_colors_cor,
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  color             = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  clustering_method = "ward.D2",
  main              = "Sample-Sample Correlation (Spearman)\nTop 500 Variable Genes",
  border_color      = NA
)

dev.off()
cat("  ✓ Saved results/sample_correlation_heatmap.png\n\n")

# ── 9. Expression Box Plot for Key Genes ─────────────────────────────────────

cat("Generating expression box plots for key breast cancer genes...\n")

# Get SE object to access gene names
se_file <- if (file.exists(file.path(DATA_DIR, "se_object_quick.rds"))) {
  file.path(DATA_DIR, "se_object_quick.rds")
} else {
  file.path(DATA_DIR, "se_object_full.rds")
}

if (file.exists(se_file)) {
  se <- readRDS(se_file)
  gene_info <- as.data.frame(rowData(se))
  
  key_genes <- c("ESR1", "ERBB2", "MKI67", "TOP2A", "TP53", "FOXA1")
  
  # Find gene IDs for these genes
  key_gene_ids <- gene_info$gene_id[gene_info$gene_name %in% key_genes]
  key_gene_ids_in_vst <- key_gene_ids[key_gene_ids %in% rownames(vst_mat)]
  
  if (length(key_gene_ids_in_vst) > 0) {
    # Get expression data
    expr_long <- data.frame()
    for (gid in key_gene_ids_in_vst) {
      gname <- gene_info$gene_name[gene_info$gene_id == gid][1]
      expr_vals <- vst_mat[gid, ]
      conditions <- col_data[names(expr_vals), "condition"]
      
      df_tmp <- data.frame(
        Gene = gname,
        Expression = as.numeric(expr_vals),
        Condition = conditions,
        stringsAsFactors = FALSE
      )
      expr_long <- rbind(expr_long, df_tmp)
    }
    
    p_box <- ggplot(expr_long, aes(x = Condition, y = Expression, fill = Condition)) +
      geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
      facet_wrap(~ Gene, scales = "free_y", nrow = 2) +
      scale_fill_manual(values = c(Normal = "#2166AC", Tumor = "#B2182B")) +
      labs(
        title = "Expression of Key Breast Cancer Genes",
        subtitle = "VST-normalized expression (Tumor vs Normal)",
        y = "VST Expression",
        x = NULL
      ) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "none",
            strip.text = element_text(face = "bold.italic"))
    
    ggsave(file.path(RESULTS_DIR, "key_genes_boxplot.png"), p_box,
           width = 10, height = 6, dpi = 300)
    cat("  ✓ Saved results/key_genes_boxplot.png\n")
  }
}

# ── Done ──────────────────────────────────────────────────────────────────────

cat("\n============================================================\n")
cat("  ✓ All visualizations generated!\n")
cat("  Output files in results/:\n")

plot_files <- list.files(RESULTS_DIR, pattern = "\\.png$", full.names = FALSE)
for (f in plot_files) {
  cat(sprintf("    • %s\n", f))
}

cat("\n  Next: Run scripts/05_optional_survival_mutation.R (optional)\n")
cat("============================================================\n")

cat("\n--- Session Info ---\n")
sessionInfo()
