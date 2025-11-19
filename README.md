# limma-DGE-Analysis-and-Integrated-Visualization-Pipeline-for-Breast-Cancer-
limma DGE Analysis and Integrated Visualization Pipeline for Breast Cancer Subtypes (GSE45827 Basal vs. Luminal A)
# 🧬 Microarray Differential Gene Expression (DGE) Pipeline: Breast Cancer Subtypes (GSE45827)

This R script automates a comprehensive bioinformatics pipeline for Differential Gene Expression (DGE) analysis of microarray data, specifically comparing **Basal-like (TNBC)** breast cancer samples against **Luminal A** samples using the publicly available **GSE45827** dataset from the Gene Expression Omnibus (GEO).

The pipeline utilizes the robust **`limma`** package for DGE and generates seven high-quality diagnostic and results plots, ensuring a complete and reproducible analysis.

## 🚀 Key Features

* **Automated Data Retrieval:** Downloads raw expression data and phenotype information directly from GEO using the **`GEOquery`** package.
* **Targeted Filtering:** Specifically filters the dataset to include only **Basal** and **Luminal A** samples.
* **Robust DGE Analysis:** Performs DGE using the **`limma`** package, which is optimized for microarray data.
* **Comprehensive QC & Visualization:** Generates a suite of seven plots covering quality control (QC) and final results visualization.
* **Reproducible Output:** All results (CSV and PNG plots) are saved to a dedicated output directory.

## 🔬 Analysis Overview

| Component | Method / Test | Purpose |
| :--- | :--- | :--- |
| **DGE** | `limma` (Linear Models for Microarray Data) | Identifies genes significantly up or down-regulated between Basal and Luminal A subtypes. |
| **Normalization** | Log2 Transformation (if needed) | Ensures data meets parametric assumptions for the `limma` model. |
| **Comparison** | Basal vs. Luminal A | Investigates expression differences between the most aggressive (Basal) and least aggressive (Luminal A) common breast cancer subtypes. |

---

## 🛠️ Prerequisites and Setup

### 📦 Packages

This script requires several Bioconductor and CRAN packages. The script includes an automated check and installation routine for all necessary libraries:
* **Bioconductor:** `GEOquery`, `limma`, `EnhancedVolcano`, `Biobase`, `genefilter`, `matrixStats`
* **CRAN:** `pheatmap`, `RColorBrewer`, `ggplot2`, `ggrepel`, `magrittr`, `dplyr`

### ⚙️ Execution

1.  **Download** the `limma DGE Analysis and Integrated Visualization Pipeline for Breast Cancer Subtypes (GSE45827 Basal vs. Luminal A).R` file.
2.  **Optional:** Modify the `output_path` variable at the beginning of the script to your desired saving location.
    ```R
    output_path <- "D:/DOWNLOADS" # Change this path
    ```
3.  **Execute** the script in your R environment:
    ```R
    source("limma DGE Analysis and Integrated Visualization Pipeline for Breast Cancer Subtypes (GSE45827 Basal vs. Luminal A).R")
    ```

---

## 📁 Output Files (7 Plots + 1 CSV)

All output files are saved to the specified `output_path`.

### Statistical Results

| Filename | Type | Description |
| :--- | :--- | :--- |
| `significant_genes_basal_vs_luminalA.csv` | CSV | Table containing all genes with an adjusted p-value (padj) $< 0.05$. Includes log2FoldChange, padj, and baseMean. |

### Visualization and QC Plots

| Filename | Type | Analysis Stage | Description |
| :--- | :--- | :--- | :--- |
| `pca_plot_GSE45827.png` | PNG | QC | **Principal Component Analysis (PCA)** plot for visual clustering of Basal vs. Luminal A samples. |
| `sample_distance_heatmap_GSE45827.png` | PNG | QC | **Sample-to-Sample Distance Heatmap** to assess batch effects and inter-group separation. |
| `heatmap_top50_variable_genes_GSE45827.png` | PNG | QC | **Heatmap of Top 50 Most Variable Genes** (independent of DGE results). |
| `heatmap_top50_significant_genes_GSE45827.png` | PNG | Results | **Heatmap of Top 50 Significant DEGs** (ranked by padj). |
| `enhanced_volcano_plot_GSE45827.png` | PNG | Results | **Volcano Plot** illustrating log2FC vs. $\log_{10}(P_{\text{adj}})$, highlighting significant DEGs. |
| `ma_plot_GSE45827.png` | PNG | Results | **MA Plot** (Mean-Difference Plot) showing log2FC (M) against Average Expression (A), colored by significance. |
| `top_gene_boxplot.png` | PNG | Results | **Boxplot** illustrating the expression distribution of the single most significant DEG across the two subtypes. |

---
*Developed using the limma DGE Analysis Pipeline for Bioinformatics.*
