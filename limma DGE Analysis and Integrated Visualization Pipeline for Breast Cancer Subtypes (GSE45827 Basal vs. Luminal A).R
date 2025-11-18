###################################################################
### Final DGE Analysis Script (Using limma for Microarray Data)
###
### GSE Title: Expression data from Breast cancer subtypes
### GSE ID: GSE45827
### Comparison: Basal-like (TNBC) vs. Luminal A
###################################################################

# --- 0. Setup: Define Output Path ---

output_path <- "D:/DOWNLOADS"
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}
cat("All output files will be saved to:", output_path, "\n")


# --- 0. Setup: Install and Load All Packages ---

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# List of Bioconductor packages (Functional packages removed)
bioc_packages <- c("GEOquery", "limma", "EnhancedVolcano", "Biobase", "genefilter", 
                   "matrixStats")

cat("Checking and installing Bioconductor packages...\n")
for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

# List of CRAN packages
cran_packages <- c("pheatmap", "RColorBrewer", "ggplot2", "ggrepel", "magrittr", "dplyr")

cat("Checking and installing CRAN packages...\n")
for (pkg in cran_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Load all required libraries
cat("Loading all required libraries...\n")
library(GEOquery)
library(limma)
library(Biobase) 
library(genefilter) 
library(matrixStats)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
library(dplyr)


###################################################################
### PART 1: Data Loading and Normalization (GSE45827)
###################################################################

cat("--- Starting Part 1: Data Loading ---\n")

# --- 1.1. Download Data from GEO ---
gse_id <- "GSE45827"
cat("Loading data for", gse_id, "from GEO...\n")
gse <- getGEO(gse_id, GSEMatrix = TRUE)
data <- gse[[1]]

# --- 1.2. Extract and Filter Data (Final, Corrected Version) ---
expression_data <- exprs(data)
pheno_data <- pData(data)

cat("Targeting column 'characteristics_ch1.1' for sample groups.\n")
subtype <- pheno_data$`characteristics_ch1.1`

correct_labels_to_keep <- c("tumor subtype: Basal", "tumor subtype: Luminal A")
keep_samples <- which(subtype %in% correct_labels_to_keep)

if (length(keep_samples) == 0) {
  stop("FATAL ERROR: Could not find 'tumor subtype: Basal' or 'tumor subtype: Luminal A'.")
}

expression_data <- expression_data[, keep_samples]
subtype <- subtype[keep_samples]

# Rename labels to Basal and LuminalA for clean analysis/plots
cat("Renaming labels for analysis...\n")
subtype[subtype == "tumor subtype: Basal"] <- "Basal"
subtype[subtype == "tumor subtype: Luminal A"] <- "LuminalA"

col_data <- data.frame(condition = factor(subtype))
col_data$condition <- relevel(col_data$condition, ref = "LuminalA")
design <- model.matrix(~condition, data = col_data) 

# --- 1.3. Data Transformation ---
if (max(expression_data) > 20) { 
  cat("Applying log2 transformation...\n")
  expression_data <- log2(expression_data)
} else {
  cat("Data appears to be already log-transformed.\n")
}
normalized_data <- expression_data 


###################################################################
### PART 2: limma DGE Analysis and Diagnostic QC Plots
###################################################################

cat("--- Starting Part 2: limma DGE Analysis and QC Plots ---\n")

# --- 2.1. limma DGE Analysis ---
cat("Running limma DGE analysis...\n")
fit <- lmFit(normalized_data, design) 
cont.matrix <- makeContrasts(BasalVsLuminalA = conditionBasal, levels = design) 
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2) 

# Extract results, renaming columns for plot compatibility
results_limma <- topTable(fit2, coef = "BasalVsLuminalA", number = Inf, adjust.method = "BH") %>%
  rename(log2FoldChange = logFC, padj = adj.P.Val, baseMean = AveExpr)

# --- 2.2. Get Annotation Data ---
annotation_df <- data.frame(Subtype = col_data$condition)
rownames(annotation_df) <- colnames(normalized_data)
annotation_colors <- list(Subtype = c(Basal = "red", LuminalA = "blue"))


# --- 2.3. PLOT 1: PCA Plot (Global Clustering) ---
cat("Generating and saving PCA Plot...\n")
pca_res <- prcomp(t(normalized_data))
pca_data <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2], 
                       condition = col_data$condition)

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  labs(title = "PCA of Basal vs. Luminal A Samples (GSE45827)") +
  scale_color_manual(values = annotation_colors$Subtype) +
  theme_minimal()
ggsave(file.path(output_path, "pca_plot_GSE45827.png"), plot = pca_plot)


# --- 2.4. PLOT 2: Sample-to-Sample Distance Heatmap ---
cat("Generating and saving Sample Distance Heatmap...\n")
sample_dist_matrix <- as.dist(1 - cor(normalized_data)) 
sample_dist_matrix_plot <- as.matrix(sample_dist_matrix)

png(file.path(output_path, "sample_distance_heatmap_GSE45827.png"), width=1000, height=1000)
pheatmap(sample_dist_matrix_plot,
         clustering_distance_rows = sample_dist_matrix,
         clustering_distance_cols = sample_dist_matrix,
         main = "Sample-to-Sample Distance Heatmap (GSE45827)",
         annotation_col = annotation_df,
         annotation_colors = annotation_colors,
         show_rownames = FALSE,
         show_colnames = FALSE)
dev.off()


# --- 2.5. PLOT 3: Top 50 Variable Genes Heatmap ---
cat("Generating and saving Top 50 Variable Gene Heatmap...\n")
top_50_genes_indices <- head(order(rowSds(normalized_data), decreasing = TRUE), 50)
heatmap_data <- normalized_data[top_50_genes_indices, ]

pheatmap(heatmap_data, 
         annotation_col = annotation_df, 
         scale = "row",
         show_rownames = FALSE,
         cluster_rows = TRUE, 
         main = "Heatmap of Top 50 Variable Genes (GSE45827)",
         annotation_colors = annotation_colors,
         filename = file.path(output_path, "heatmap_top50_variable_genes_GSE45827.png"))


###################################################################
### PART 3: limma Results and Visualization (Total 5 Results Plots)
###################################################################

cat("--- Starting Part 3: limma Results and Visualization ---\n")

# --- 3.1. Save Significant Gene Results (CSV) ---
significant_genes <- results_limma %>%
  filter(padj < 0.05) %>%
  select(log2FoldChange, padj, baseMean, everything()) 

cat("Number of significant DEGs found:", nrow(significant_genes), "\n")
cat("Saving significant gene list to significant_genes_basal_vs_luminalA.csv...\n")
write.csv(significant_genes, 
          file = file.path(output_path, "significant_genes_basal_vs_luminalA.csv"))


# --- 3.2. PLOT 4: Heatmap of Top 50 Significant Genes ---
cat("Generating and saving Top 50 Significant Gene Heatmap...\n")
if (nrow(significant_genes) > 0) {
  # Filter top 50 significant genes
  top_50_sig_genes <- significant_genes %>%
    arrange(padj) %>%
    head(50)
  # Get the expression data for these genes
  top_50_sig_heatmap_data <- normalized_data[rownames(top_50_sig_genes), ]
  
  pheatmap(top_50_sig_heatmap_data, 
           annotation_col = annotation_df, 
           scale = "row", 
           cluster_rows = TRUE,
           show_rownames = TRUE, 
           fontsize_row = 8,
           main = "Heatmap of Top 50 Significant DEGs (GSE45827)",
           annotation_colors = annotation_colors,
           filename = file.path(output_path, "heatmap_top50_significant_genes_GSE45827.png"))
} else {
  cat("Skipping Top 50 Significant Heatmap: Zero significant genes found.\n")
}


# --- 3.3. PLOT 5: Enhanced Volcano Plot ---
cat("Generating and saving Enhanced Volcano Plot...\n")
ev_plot <- EnhancedVolcano(results_limma,
                           lab = rownames(results_limma),
                           x = 'log2FoldChange',
                           y = 'padj',
                           title = 'Volcano Plot (Basal vs. Luminal A)',
                           subtitle = 'GSE45827',
                           pCutoff = 0.05,
                           FCcutoff = 1.0,
                           pointSize = 2.0,
                           labSize = 4.0,
                           legendPosition = 'bottom',
                           colAlpha = 0.5)
ggsave(file.path(output_path, "enhanced_volcano_plot_GSE45827.png"), plot = ev_plot, width = 12, height = 10)


# --- 3.4. PLOT 6: MA Plot ---
cat("Generating and saving MA Plot...\n")
results_limma$significant <- ifelse(results_limma$padj < 0.05, "Significant", "Not significant")

ma_plot <- ggplot(results_limma, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(alpha = 0.7, color = "gray") +
  geom_point(data = subset(results_limma, padj < 0.05), aes(color = significant), size = 2) +
  scale_color_manual(values = c("red")) +
  labs(title = "MA Plot (GSE45827)", x = "Average Expression", y = "log2 Fold Change") +
  theme_minimal()
ggsave(file.path(output_path, "ma_plot_GSE45827.png"), plot = ma_plot)


# --- 3.5. PLOT 7: Boxplot of Top Gene ---
cat("Generating and saving boxplot for the #1 significant gene...\n")
if (nrow(significant_genes) > 0) {
  top_gene_id <- rownames(significant_genes)[1]
  
  top_gene_exp <- data.frame(
    condition = col_data$condition,
    expression = normalized_data[top_gene_id, ]
  )
  
  top_gene_boxplot <- ggplot(top_gene_exp, aes(x = condition, y = expression, fill = condition)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.7) +
    scale_fill_manual(values = c("LuminalA" = "blue", "Basal" = "red")) +
    labs(title = paste("Expression of Top Gene:", top_gene_id),
         x = "Sample Group",
         y = "Log2 Expression Value") +
    theme_minimal()
  ggsave(file.path(output_path, "top_gene_boxplot.png"), plot = top_gene_boxplot)
} else {
  cat("Skipping Boxplot: Zero significant genes found.\n")
}


cat("\n#################################################\n")
cat("Analysis complete. All 7 plots and CSV saved to:", output_path, "\n")
cat("#################################################\n")
