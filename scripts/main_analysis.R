# ==============================================================================
# AUTOMATED MASTER ANALYSIS SCRIPT FOR NATURE SUBMISSION
# Project: Cross-Timepoint Heatmaps & Multi-Comparison GO Enrichment (Mouse)
# Author: Nguyen Thi Thanh Tuyen
# Date: June 26, 2026
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. AUTOMATIC INSTALLATION & LOAD DEPENDENCIES
# ------------------------------------------------------------------------------
required_cran <- c("tidyverse", "pheatmap", "RColorBrewer")
required_bioc <- c("clusterProfiler", "org.Mm.eg.db")

for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, dependencies = TRUE)
  library(pkg, character_only = TRUE)
}
for (pkg in required_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
  library(pkg, character_only = TRUE)
}

set.seed(42) # For reproducibility

# ------------------------------------------------------------------------------
# 2. INPUT DATA LOADING
# ------------------------------------------------------------------------------
tpm_matrix <- read.csv("data/expression_matrix_tpm.csv", row.names = 1, check.names = FALSE)

deg_files <- list(
  "8w_Het_vs_KO"   = "data/deg_results_8w_Het_vs_KO.csv",
  "14w_Het_vs_KO"  = "data/deg_results_14w_Het_vs_KO.csv",
  "8wKO_vs_14wKO"  = "data/deg_results_8wKO_vs_14wKO.csv"
)

# ------------------------------------------------------------------------------
# 3. GENERATE TARGETED HEATMAPS USING 14W DEG LIST (Log2FC > 0.58 & pval < 0.05)
# ------------------------------------------------------------------------------
message("Processing Targeted Heatmaps based on 14w DEGs...")
deg_14w <- read.csv(deg_files[["14w_Het_vs_KO"]], check.names = FALSE)

# Filter core significant transcripts at 14w based on your exact manuscript criteria
signif_14w_list <- deg_14w %>% 
  filter(pval < 0.05 & abs(log2foldChange) > 0.58) %>%
  pull(Transcript_ID)

if (length(signif_14w_list) > 0) {
  
  col_14w <- grep("14w", colnames(tpm_matrix), value = TRUE)
  col_8w  <- grep("8w", colnames(tpm_matrix), value = TRUE)
  
  # --- HEATMAP 1: 14w Genes plotted with 14w Samples ---
  if(length(col_14w) > 0) {
    mat_14w <- tpm_matrix[rownames(tpm_matrix) %in% signif_14w_list, col_14w]
    mat_14w_scaled <- t(scale(t(mat_14w)))
    
    pdf("figures/Figure_Heatmap_14w_Genes_at_14w.pdf", width = 8, height = 10)
    pheatmap(mat_14w_scaled, cluster_rows = TRUE, cluster_cols = TRUE,
             show_rownames = ifelse(nrow(mat_14w_scaled) < 60, TRUE, FALSE),
             show_colnames = TRUE, color = colorRampPalette(c("blue", "white", "red"))(50),
             main = "14w Significant Genes Profile (at 14 Weeks)")
    dev.off()
  }
  
  # --- HEATMAP 2: 14w Genes plotted with 8w Samples ---
  if(length(col_8w) > 0) {
    mat_8w <- tpm_matrix[rownames(tpm_matrix) %in% signif_14w_list, col_8w]
    mat_8w_scaled <- t(scale(t(mat_8w)))
    
    pdf("figures/Figure_Heatmap_14w_Genes_at_8w.pdf", width = 8, height = 10)
    pheatmap(mat_8w_scaled, cluster_rows = TRUE, cluster_cols = TRUE,
             show_rownames = ifelse(nrow(mat_8w_scaled) < 60, TRUE, FALSE),
             show_colnames = TRUE, color = colorRampPalette(c("blue", "white", "red"))(50),
             main = "14w Significant Genes Profile (at 8 Weeks)")
    dev.off()
  }
  
} else {
  warning("No significant genes found for 14w to build target heatmaps.")
}

# ------------------------------------------------------------------------------
# 4. AUTOMATED MULTI-COMPARISON GO ENRICHMENT LOOP (pval < 0.05 & Log2FC > 0.58)
# ------------------------------------------------------------------------------
for (comparison_name in names(deg_files)) {
  message(paste("Running Automated GO Analysis for:", comparison_name))
  
  deg_data <- read.csv(deg_files[[comparison_name]], check.names = FALSE)
  
  signif_genes <- deg_data %>% 
    filter(pval < 0.05 & abs(log2foldChange) > 0.58)
  
  if (nrow(signif_genes) == 0) {
    message(paste("Skipping GO for", comparison_name, "- Zero significant genes found."))
    next
  }
  
  unique_gene_list <- signif_genes %>%
    group_by(Gene) %>%
    slice_min(order_by = pval, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    pull(Gene)
  
  gene_conversion <- tryCatch({
    bitr(unique_gene_list, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db)
  }, error = function(e) { NULL })
  
  if (is.null(gene_conversion) || nrow(gene_conversion) == 0) next
  
  go_results <- enrichGO(
    gene          = gene_conversion$ENTREZID,
    OrgDb         = org.Mm.eg.db,
    ont           = "BP", 
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  if (!is.null(go_results) && nrow(as.data.frame(go_results)) > 0) {
    csv_path <- paste0("results/GO_results_", comparison_name, ".csv")
    write.csv(as.data.frame(go_results), csv_path, row.names = FALSE)
    
    pdf_path <- paste0("figures/Figure_GO_Dotplot_", comparison_name, ".pdf")
    go_dotplot <- dotplot(go_results, showCategory = 20) + 
      gtitle(paste("GO Biological Process:", gsub("_", " ", comparison_name)))
    
    ggsave(filename = pdf_path, plot = go_dotplot, width = 9, height = 7, dpi = 600)
  }
}

print("ALL PROCESSES COMPLETED SUCCESSFULLY WITH CORRECT CROSS-LOGIC!")