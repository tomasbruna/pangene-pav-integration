# Evolutionary and methodological considerations when interpreting gene presence–absence variation in pangenomes — Analyses

This repository contains the scripts and notebooks used to generate the analyses for the [Evolutionary and methodological considerations when interpreting gene presence–absence variation in pangenomes](https://www.biorxiv.org/content/10.1101/2025.08.14.670405v1) paper.


All required datasets are automatically downloaded from [Zenodo](https://doi.org/10.5281/zenodo.16809323). Original data sources are listed in Supplementary Tables 1 and 2 of the paper.

## Repository structure

```
├── Cotton/                     # Cotton pangenome analysis notebook
│   └── cotton.ipynb
├── Soybean/                    # Soybean pangenome analysis notebook and metadata
│   ├── soybean.ipynb
│   └── soybean_md.txt
└── scripts/                    # Supporting scripts for analysis and figure generation
```

---

## Reproducing the analyses

1. **Clone the repository**
```bash
git clone https://github.com/username/pangene-pav-integration.git
cd pangene-pav-integration
```

2. **Work through the notebooks**
   - **Cotton**: `Cotton/cotton.ipynb`
   - **Soybean**: `Soybean/soybean.ipynb`

   The notebooks contain runnable code cells and instructions for steps that may need to be run externally (e.g., WDL workflows in `scripts/WDL/`).

3. **Generate R figures**
   After executing the **data preparation**, **GENESPACE**, and **PanKmer** cells in both notebooks, run:
```bash
Rscript scripts/generate_R_figures.R
```
   This will produce the R-based figures from the paper and save them in the `R_figures` output directory.

---
