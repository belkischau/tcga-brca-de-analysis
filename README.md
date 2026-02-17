# 🧬 Breast Cancer Tumor vs Normal Differential Gene Expression Analysis

**One-sentence summary:** A fully automated, reproducible pipeline that downloads TCGA-BRCA RNA-seq data from the GDC, identifies differentially expressed genes between primary breast tumors and adjacent normal tissue, performs pathway enrichment analysis, and generates publication-quality figures — all in R.

---

## 📖 Biological Background

Breast cancer is the most commonly diagnosed cancer worldwide and a leading cause of cancer-related death in women. The Cancer Genome Atlas (TCGA) Breast Invasive Carcinoma (BRCA) project profiled over 1,000 primary tumor and matched normal samples using RNA sequencing (STAR-Counts workflow), providing a rich resource for understanding the transcriptomic landscape of breast cancer.

Differential gene expression (DE) analysis between tumor and adjacent normal tissue reveals genes that are systematically up- or down-regulated in cancer. These differentially expressed genes (DEGs) frequently map to the **hallmarks of cancer** — including sustained proliferative signaling (*ERBB2/HER2*, *EGFR*), evasion of growth suppressors (*TP53*, *RB1*), resistance to apoptosis (*BCL2*), and activation of invasion/metastasis programs (*MMP* family). By combining DE analysis with Gene Ontology (GO) and KEGG pathway enrichment, we can draw actionable biological conclusions about the molecular mechanisms driving breast cancer.

---

## 🛠️ Technologies & Packages

| Category | Tools |
|---|---|
| Language | R (≥ 4.3) |
| Data Access | [TCGAbiolinks](https://bioconductor.org/packages/TCGAbiolinks/) |
| Differential Expression | [DESeq2](https://bioconductor.org/packages/DESeq2/) |
| Pathway Enrichment | [clusterProfiler](https://bioconductor.org/packages/clusterProfiler/), [org.Hs.eg.db](https://bioconductor.org/packages/org.Hs.eg.db/) |
| Visualization | [EnhancedVolcano](https://bioconductor.org/packages/EnhancedVolcano/), [pheatmap](https://CRAN.R-project.org/package=pheatmap), [ggplot2](https://ggplot2.tidyverse.org/) |
| Survival Analysis | [survminer](https://CRAN.R-project.org/package=survminer), [survival](https://CRAN.R-project.org/package=survival) |
| Mutation Analysis | [maftools](https://bioconductor.org/packages/maftools/) |
| Reporting | [rmarkdown](https://rmarkdown.rstudio.com/), [knitr](https://yihui.org/knitr/) |
| Reproducibility | [renv](https://rstudio.github.io/renv/) |

---

## 🔑 Key Results Summary

> *Results below are representative of the **quick-test** run (100 tumor + ~11 normal samples). Full-data results will differ slightly.*

| Metric | Value (approx.) |
|---|---|
| Tumor samples | 100 |
| Normal samples | ~11* |
| Total genes tested | ~58,000 (pre-filter) → ~20,000 (after filtering) |
| Significant DEGs (|log2FC| > 1 & padj < 0.05) | ~5,000–7,000 |

\*GDC has ~113 normal samples total for TCGA-BRCA; the quick-test uses all available normals up to 100.

### Known breast cancer gene validation:
- **ESR1** (Estrogen Receptor): significantly differentially expressed
- **ERBB2** (HER2): upregulated in tumors
- **TP53**: differential expression detected
- **MKI67** (Ki-67 proliferation marker): upregulated in tumors
- **TOP2A**: upregulated in tumors (proliferation hallmark)

### Top enriched pathways:
- **Cell cycle** (KEGG)
- **DNA replication** (KEGG)
- **Extracellular matrix organization** (GO)
- **PI3K-Akt signaling pathway** (KEGG)
- **Focal adhesion** (KEGG)

---

## 🧪 Biological Interpretation — Cancer Hallmarks

The DEGs and enriched pathways map directly to the **Hallmarks of Cancer** (Hanahan & Weinberg, 2011):

1. **Sustaining Proliferative Signaling** — Upregulation of *ERBB2*, *EGFR*, *MKI67*, *TOP2A*, and enrichment of the cell cycle/DNA replication pathways.
2. **Evading Growth Suppressors** — Differential expression of *TP53*, *RB1*, *CDKN2A*.
3. **Activating Invasion & Metastasis** — Upregulation of *MMP* family genes, enrichment of ECM-receptor and focal adhesion pathways.
4. **Inducing Angiogenesis** — Enrichment of VEGF signaling pathway.
5. **Genome Instability & Mutation** — DNA repair pathways enriched; *BRCA1*/*BRCA2* differentially expressed.
6. **Deregulating Cellular Energetics** — Changes in metabolic pathways (oxidative phosphorylation, glycolysis).

---

## 🚀 How to Reproduce

### Prerequisites
- R ≥ 4.3 installed
- Internet connection (for first-time GDC data download)
- ~2 GB disk space (quick-test) or ~20 GB (full analysis)

### Step-by-step

```bash
# 1. Clone the repo
git clone https://github.com/YOUR_USERNAME/tcga-brca-de-analysis.git
cd tcga-brca-de-analysis

# 2. (Optional) Set up renv for exact reproducibility
Rscript -e 'install.packages("renv"); renv::restore()'

# 3. Run the pipeline — just run script 01, it will auto-download everything!
Rscript scripts/01_download_and_prepare.R          # Downloads data, saves .rds files
Rscript scripts/02_DESeq2_analysis.R               # Runs DESeq2 differential expression
Rscript scripts/03_pathway_enrichment.R            # GO & KEGG enrichment
Rscript scripts/04_visualization.R                 # Generate all publication-quality plots
Rscript scripts/05_optional_survival_mutation.R    # Survival analysis & oncoplot (optional)

# 4. (Optional) Render the full report
Rscript -e 'rmarkdown::render("reports/TCGA_BRCA_DE_Report.Rmd")'
```

> **Quick-test mode (default):** By default, the pipeline runs in quick-test mode (~100 tumor + available normal samples) and completes in **< 10 minutes**. To run the full analysis, edit `scripts/01_download_and_prepare.R` and set `RUN_MODE <- "full"`.

---

## 📁 Repository Structure

```
tcga-brca-de-analysis/
├── README.md                              # This file
├── .gitignore                             # Git ignore rules
├── renv.lock                              # Locked package versions (after renv::snapshot())
├── setup_renv.R                           # Helper to initialize renv
├── scripts/
│   ├── 01_download_and_prepare.R          # Auto-download TCGA-BRCA data via GDC
│   ├── 02_DESeq2_analysis.R              # Differential expression with DESeq2
│   ├── 03_pathway_enrichment.R           # GO & KEGG pathway enrichment
│   ├── 04_visualization.R               # Volcano, heatmap, PCA, bar plots
│   └── 05_optional_survival_mutation.R   # Kaplan-Meier survival + maftools oncoplot
├── reports/
│   └── TCGA_BRCA_DE_Report.Rmd           # Full analysis report (Rmarkdown)
├── results/                               # Generated outputs (plots, CSVs)
│   ├── DEGs_full.csv
│   ├── DEGs_significant.csv
│   ├── volcano_plot.png
│   ├── top50_heatmap.png
│   ├── GO_dotplot.png
│   ├── KEGG_dotplot.png
│   └── ...
└── data/                                  # Raw & processed data (mostly .gitignored)
    ├── .gitignore
    └── README.md
```

---

## 📝 Notes

- **GDC Token:** Public TCGA data does not require a GDC token. If you encounter access errors, ensure you are querying open-access data only.
- **Nextflow Integration:** This pipeline can be wrapped into a [Nextflow](https://www.nextflow.io/) workflow for HPC execution. Each R script maps to a Nextflow process; the `data/` and `results/` directories serve as channel inputs/outputs.
- **Full mode:** For publication-level results, set `RUN_MODE <- "full"` in `scripts/01_download_and_prepare.R` to use all ~1,100 tumor and ~113 normal samples.

---

## 📜 License

This project is licensed under the **MIT License** — see below.

```
MIT License

Copyright (c) 2026

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

---

## 🙏 Acknowledgments

- [The Cancer Genome Atlas (TCGA)](https://www.cancer.gov/tcga)
- [Genomic Data Commons (GDC)](https://gdc.cancer.gov/)
- [Bioconductor](https://www.bioconductor.org/)
- Hanahan, D., & Weinberg, R. A. (2011). Hallmarks of cancer: the next generation. *Cell*, 144(5), 646–674.
