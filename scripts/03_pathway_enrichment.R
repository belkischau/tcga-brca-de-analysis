#!/usr/bin/env Rscript
# =============================================================================
# Script 03: Pathway Enrichment Analysis (GO & KEGG)
# =============================================================================
#
# Description:
#   Performs Gene Ontology (GO) and KEGG pathway enrichment analysis on
#   differentially expressed genes from the DESeq2 results. Uses
#   clusterProfiler for over-representation analysis (ORA) and gene set
#   enrichment analysis (GSEA).
#
# Input:
#   results/deseq2_results_annotated.rds
#
# Output:
#   results/GO_BP_enrichment.csv
#   results/GO_BP_enrichment.rds
#   results/KEGG_enrichment.csv
#   results/KEGG_enrichment.rds
#   results/GSEA_results.rds
#   results/GO_dotplot.png
#   results/KEGG_dotplot.png
#   results/GSEA_plot.png
#
# Usage:
#   Rscript scripts/03_pathway_enrichment.R
#
# Author: DTU Bioinformatics MSc | 2026
# =============================================================================

# ── 0. Configuration ─────────────────────────────────────────────────────────

set.seed(123)

cat("============================================================\n")
cat("  Pathway Enrichment Analysis (GO & KEGG)\n")
cat(sprintf("  Date: %s\n", Sys.time()))
cat("============================================================\n\n")

PROJECT_DIR <- getwd()
RESULTS_DIR <- file.path(PROJECT_DIR, "results")

# ── 1. Load packages ─────────────────────────────────────────────────────────

cat("Loading packages...\n")

required_packages <- c("clusterProfiler", "org.Hs.eg.db", "AnnotationDbi",
                        "enrichplot", "DOSE", "ggplot2", "dplyr", "stringr")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' not installed. Run: BiocManager::install('%s')", pkg, pkg))
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

cat("✓ Packages loaded.\n\n")

# ── 2. Load DESeq2 results ───────────────────────────────────────────────────

cat("Loading DESeq2 results...\n")

res_file <- file.path(RESULTS_DIR, "deseq2_results_annotated.rds")
if (!file.exists(res_file)) {
  stop("DESeq2 results not found! Run scripts/02_DESeq2_analysis.R first.\n",
       "  Expected: ", res_file)
}

res_df <- readRDS(res_file)

cat(sprintf("  Total genes: %d\n", nrow(res_df)))
cat(sprintf("  Significant DEGs: %d\n",
            sum(res_df$DE_status %in% c("Upregulated", "Downregulated"))))
cat("✓ Results loaded.\n\n")

# ── 3. Prepare gene lists ────────────────────────────────────────────────────

cat("Preparing gene lists for enrichment analysis...\n")

# Get significant DEGs with gene symbols
sig_genes <- res_df %>%
  filter(DE_status %in% c("Upregulated", "Downregulated")) %>%
  filter(!is.na(gene_name) & gene_name != "") %>%
  pull(gene_name) %>%
  unique()

# Separately get up- and down-regulated genes
up_genes <- res_df %>%
  filter(DE_status == "Upregulated", !is.na(gene_name), gene_name != "") %>%
  pull(gene_name) %>%
  unique()

down_genes <- res_df %>%
  filter(DE_status == "Downregulated", !is.na(gene_name), gene_name != "") %>%
  pull(gene_name) %>%
  unique()

cat(sprintf("  All significant DEGs: %d\n", length(sig_genes)))
cat(sprintf("  Upregulated: %d\n", length(up_genes)))
cat(sprintf("  Downregulated: %d\n", length(down_genes)))

# Convert gene symbols to Entrez IDs (required by KEGG)
symbol_to_entrez <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys      = sig_genes,
  columns   = "ENTREZID",
  keytype   = "SYMBOL"
)
symbol_to_entrez <- symbol_to_entrez[!is.na(symbol_to_entrez$ENTREZID), ]
symbol_to_entrez <- symbol_to_entrez[!duplicated(symbol_to_entrez$SYMBOL), ]

entrez_ids <- symbol_to_entrez$ENTREZID
cat(sprintf("  Mapped to Entrez IDs: %d\n", length(entrez_ids)))

# Background universe: all tested genes with gene symbols
universe_symbols <- res_df %>%
  filter(!is.na(gene_name) & gene_name != "") %>%
  pull(gene_name) %>%
  unique()

universe_entrez <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys      = universe_symbols,
  columns   = "ENTREZID",
  keytype   = "SYMBOL"
)
universe_entrez <- universe_entrez$ENTREZID[!is.na(universe_entrez$ENTREZID)]
universe_entrez <- unique(universe_entrez)

cat(sprintf("  Background universe: %d Entrez IDs\n", length(universe_entrez)))
cat("✓ Gene lists prepared.\n\n")

# ── 4. GO Biological Process Enrichment (ORA) ────────────────────────────────

cat("Running GO Biological Process enrichment...\n")

go_bp <- tryCatch({
  enrichGO(
    gene          = sig_genes,
    universe      = universe_symbols,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.10,
    readable      = TRUE,
    minGSSize     = 10,
    maxGSSize     = 500
  )
}, error = function(e) {
  cat(sprintf("  Warning: GO enrichment failed: %s\n", e$message))
  NULL
})

if (!is.null(go_bp) && nrow(as.data.frame(go_bp)) > 0) {
  go_df <- as.data.frame(go_bp)
  cat(sprintf("  ✓ Found %d enriched GO BP terms (padj < 0.05)\n", nrow(go_df)))
  
  # Show top 10 terms
  cat("\n  Top 10 GO Biological Process terms:\n")
  top_go <- head(go_df[, c("Description", "GeneRatio", "p.adjust", "Count")], 10)
  print(top_go, row.names = FALSE)
  
  # Save results
  readr::write_csv(go_df, file.path(RESULTS_DIR, "GO_BP_enrichment.csv"))
  saveRDS(go_bp, file.path(RESULTS_DIR, "GO_BP_enrichment.rds"))
  
  # Generate dotplot
  cat("\n  Generating GO dotplot...\n")
  p_go <- dotplot(go_bp, showCategory = 20, title = "GO Biological Process Enrichment") +
    theme(axis.text.y = element_text(size = 8))
  
  ggsave(file.path(RESULTS_DIR, "GO_dotplot.png"), p_go,
         width = 10, height = 8, dpi = 300)
  cat("  ✓ Saved results/GO_dotplot.png\n")
} else {
  cat("  ⚠ No significant GO terms found.\n")
}

cat("\n")

# ── 5. KEGG Pathway Enrichment (ORA) ─────────────────────────────────────────

cat("Running KEGG pathway enrichment...\n")

kegg_res <- tryCatch({
  enrichKEGG(
    gene          = entrez_ids,
    universe      = universe_entrez,
    organism      = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.10,
    minGSSize     = 10,
    maxGSSize     = 500
  )
}, error = function(e) {
  cat(sprintf("  Warning: KEGG enrichment failed: %s\n", e$message))
  cat("  This may be due to KEGG API issues. Continuing without KEGG results.\n")
  NULL
})

if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
  # Convert Entrez IDs to gene symbols in results for readability
  kegg_readable <- tryCatch({
    setReadable(kegg_res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  }, error = function(e) kegg_res)
  
  kegg_df <- as.data.frame(kegg_readable)
  cat(sprintf("  ✓ Found %d enriched KEGG pathways (padj < 0.05)\n", nrow(kegg_df)))
  
  # Show top 10 pathways
  cat("\n  Top 10 KEGG Pathways:\n")
  top_kegg <- head(kegg_df[, c("Description", "GeneRatio", "p.adjust", "Count")], 10)
  print(top_kegg, row.names = FALSE)
  
  # Save results
  readr::write_csv(kegg_df, file.path(RESULTS_DIR, "KEGG_enrichment.csv"))
  saveRDS(kegg_readable, file.path(RESULTS_DIR, "KEGG_enrichment.rds"))
  
  # Generate dotplot
  cat("\n  Generating KEGG dotplot...\n")
  p_kegg <- dotplot(kegg_readable, showCategory = 20,
                    title = "KEGG Pathway Enrichment") +
    theme(axis.text.y = element_text(size = 8))
  
  ggsave(file.path(RESULTS_DIR, "KEGG_dotplot.png"), p_kegg,
         width = 10, height = 8, dpi = 300)
  cat("  ✓ Saved results/KEGG_dotplot.png\n")
} else {
  cat("  ⚠ No significant KEGG pathways found.\n")
}

cat("\n")

# ── 6. Gene Set Enrichment Analysis (GSEA) ───────────────────────────────────

cat("Running Gene Set Enrichment Analysis (GSEA)...\n")

# Create ranked gene list (log2FC, named by Entrez ID)
# Map all genes (not just significant) to Entrez
all_genes_for_gsea <- res_df %>%
  dplyr::filter(!is.na(gene_name) & gene_name != "" & !is.na(log2FoldChange)) %>%
  dplyr::select(gene_name, log2FoldChange) %>%
  dplyr::distinct(gene_name, .keep_all = TRUE)

# Map to Entrez IDs
gsea_entrez_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = all_genes_for_gsea$gene_name,
  columns = "ENTREZID",
  keytype = "SYMBOL"
)
gsea_entrez_map <- gsea_entrez_map[!is.na(gsea_entrez_map$ENTREZID), ]
gsea_entrez_map <- gsea_entrez_map[!duplicated(gsea_entrez_map$SYMBOL), ]

# Merge and create named vector
gsea_input <- all_genes_for_gsea %>%
  inner_join(gsea_entrez_map, by = c("gene_name" = "SYMBOL"))

gene_list <- gsea_input$log2FoldChange
names(gene_list) <- gsea_input$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

cat(sprintf("  Ranked gene list: %d genes\n", length(gene_list)))

gsea_kegg <- tryCatch({
  gseKEGG(
    geneList      = gene_list,
    organism      = "hsa",
    minGSSize     = 15,
    maxGSSize     = 500,
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    verbose       = FALSE,
    seed          = TRUE
  )
}, error = function(e) {
  cat(sprintf("  Warning: GSEA KEGG failed: %s\n", e$message))
  NULL
})

if (!is.null(gsea_kegg) && nrow(as.data.frame(gsea_kegg)) > 0) {
  gsea_df <- as.data.frame(gsea_kegg)
  cat(sprintf("  ✓ Found %d enriched KEGG pathways via GSEA\n", nrow(gsea_df)))
  
  saveRDS(gsea_kegg, file.path(RESULTS_DIR, "GSEA_KEGG_results.rds"))
  
  # Generate GSEA plot for top pathways
  cat("  Generating GSEA plots...\n")
  
  # Ridgeplot
  p_ridge <- tryCatch({
    ridgeplot(gsea_kegg, showCategory = 15) +
      ggtitle("GSEA KEGG Pathway Ridge Plot") +
      theme(axis.text.y = element_text(size = 8))
  }, error = function(e) NULL)
  
  if (!is.null(p_ridge)) {
    ggsave(file.path(RESULTS_DIR, "GSEA_ridgeplot.png"), p_ridge,
           width = 10, height = 8, dpi = 300)
    cat("  ✓ Saved results/GSEA_ridgeplot.png\n")
  }
  
  # GSEA enrichment plot for top pathway
  p_gsea <- tryCatch({
    gseaplot2(gsea_kegg, geneSetID = 1:min(3, nrow(gsea_df)),
              pvalue_table = TRUE)
  }, error = function(e) NULL)
  
  if (!is.null(p_gsea)) {
    ggsave(file.path(RESULTS_DIR, "GSEA_plot.png"), p_gsea,
           width = 10, height = 6, dpi = 300)
    cat("  ✓ Saved results/GSEA_plot.png\n")
  }
} else {
  cat("  ⚠ No significant GSEA pathways found.\n")
}

# ── 7. GO Enrichment for Up and Down separately ──────────────────────────────

cat("\nRunning separate GO enrichment for up/down-regulated genes...\n")

# Upregulated
go_up <- tryCatch({
  enrichGO(
    gene          = up_genes,
    universe      = universe_symbols,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE,
    minGSSize     = 10,
    maxGSSize     = 500
  )
}, error = function(e) NULL)

if (!is.null(go_up) && nrow(as.data.frame(go_up)) > 0) {
  cat(sprintf("  ✓ Upregulated: %d enriched GO BP terms\n",
              nrow(as.data.frame(go_up))))
  saveRDS(go_up, file.path(RESULTS_DIR, "GO_BP_upregulated.rds"))
} else {
  cat("  ⚠ No enriched GO terms for upregulated genes.\n")
}

# Downregulated
go_down <- tryCatch({
  enrichGO(
    gene          = down_genes,
    universe      = universe_symbols,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE,
    minGSSize     = 10,
    maxGSSize     = 500
  )
}, error = function(e) NULL)

if (!is.null(go_down) && nrow(as.data.frame(go_down)) > 0) {
  cat(sprintf("  ✓ Downregulated: %d enriched GO BP terms\n",
              nrow(as.data.frame(go_down))))
  saveRDS(go_down, file.path(RESULTS_DIR, "GO_BP_downregulated.rds"))
} else {
  cat("  ⚠ No enriched GO terms for downregulated genes.\n")
}

# ── 8. Verify cancer-related pathways ────────────────────────────────────────

cat("\n╔══════════════════════════════════════════════╗\n")
cat("║     PATHWAY VERIFICATION — CANCER HALLMARKS  ║\n")
cat("╠══════════════════════════════════════════════╣\n")

cancer_keywords <- c("cell cycle", "DNA replication", "apoptosis",
                      "p53 signaling", "PI3K", "breast cancer",
                      "ECM", "focal adhesion", "angiogenesis",
                      "Wnt signaling", "MAPK", "mTOR")

# Check KEGG results
if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
  kegg_terms <- as.data.frame(kegg_readable)$Description
  for (kw in cancer_keywords) {
    found <- grep(kw, kegg_terms, ignore.case = TRUE, value = TRUE)
    if (length(found) > 0) {
      cat(sprintf("║  ✓ KEGG: %-38s ║\n", found[1]))
    }
  }
}

# Check GO results
if (!is.null(go_bp) && nrow(as.data.frame(go_bp)) > 0) {
  go_terms <- as.data.frame(go_bp)$Description
  go_cancer_words <- c("cell cycle", "cell division", "DNA repair",
                        "apoptot", "extracellular matrix", "angiogenesis",
                        "immune", "proliferat", "migration")
  for (kw in go_cancer_words) {
    found <- grep(kw, go_terms, ignore.case = TRUE, value = TRUE)
    if (length(found) > 0) {
      cat(sprintf("║  ✓ GO:   %-38s ║\n", substr(found[1], 1, 38)))
    }
  }
}

cat("╚══════════════════════════════════════════════╝\n")

# ── Done ──────────────────────────────────────────────────────────────────────

cat("\n============================================================\n")
cat("  ✓ Pathway enrichment analysis complete!\n")
cat("  Outputs saved to results/\n")
cat("  Next: Run scripts/04_visualization.R\n")
cat("============================================================\n")

cat("\n--- Session Info ---\n")
sessionInfo()
