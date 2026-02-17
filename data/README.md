This directory stores raw and processed TCGA-BRCA data.

- Raw GDC downloads are placed here automatically by `scripts/01_download_and_prepare.R`.
- Processed `.rds` files (counts matrix, sample metadata) are saved here for use by downstream scripts.
- Raw data is **not** tracked by git (see `.gitignore`).
