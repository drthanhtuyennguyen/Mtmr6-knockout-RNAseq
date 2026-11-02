# ==============================================================================
# Script Name : volcano_plots_analysis.R
# Project     : Mtmr6 RNA-seq Project (Het vs KO at 8w and 14w)
# Description : Publication-ready Volcano Plots for 8w and 14w Het vs KO
# Author      : Research Team
# ==============================================================================

# 1. Setup Environment & Packages ----------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})

# 2. Global Parameters & Directory Setup ---------------------------------------
p_adj_cut  <- 0.05
log2fc_cut <- 0.58

# Relative directory paths (Standard for publication code submission)
data_dir   <- "../data/"
output_dir <- "../figures/"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ==============================================================================
# PART 1: 8w Het vs KO Volcano Plot (No Legend)
# ==============================================================================

file_8w <- file.path(data_dir, "deg_results_8w_Het_vs_KO.csv")
deg_raw_8w <- read.csv(file_8w, stringsAsFactors = FALSE)

# Detect Beta/Log2FC and Adjusted P-value column names
x_col_8w <- names(deg_raw_8w)[grep("b$|beta|log2fc|foldchange|logfc", names(deg_raw_8w), ignore.case = TRUE)][1]
y_col_8w <- names(deg_raw_8w)[grep("qval|adj|padj|pvalue|p.value", names(deg_raw_8w), ignore.case = TRUE)][1]

deg_data_8w <- data.frame(
  transcript_id = as.character(deg_raw_8w[[1]]),
  gene_symbol   = as.character(deg_raw_8w[[5]]),
  x_val         = as.numeric(deg_raw_8w[[x_col_8w]]),
  y_val         = as.numeric(deg_raw_8w[[y_col_8w]]),
  stringsAsFactors = FALSE
)

# Fill missing gene symbols with transcript IDs
deg_data_8w$gene_symbol <- ifelse(
  is.na(deg_data_8w$gene_symbol) | deg_data_8w$gene_symbol == "" | deg_data_8w$gene_symbol == "NA",
  deg_data_8w$transcript_id,
  deg_data_8w$gene_symbol
)

# Aggregate to Gene level by selecting transcript with max effect size
deg_clean_8w <- deg_data_8w %>%
  filter(!is.na(x_val) & !is.na(y_val) & y_val > 0) %>%
  group_by(gene_symbol) %>%
  slice_max(order_by = abs(x_val), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    neg_log10_padj = -log10(y_val),
    Expression = case_when(
      x_val >= log2fc_cut & y_val < p_adj_cut ~ "Up-regulated",
      x_val <= -log2fc_cut & y_val < p_adj_cut ~ "Down-regulated",
      TRUE ~ "Not Significant"
    )
  )

# Strict filter for labels: only true DEGs at 8w
top_genes_8w <- deg_clean_8w %>%
  filter(y_val < p_adj_cut & abs(x_val) >= log2fc_cut) %>%
  arrange(y_val)

volcano_8w <- ggplot(deg_clean_8w, aes(x = x_val, y = neg_log10_padj)) +
  geom_point(aes(color = Expression), alpha = 0.65, size = 1.6) +
  scale_color_manual(
    values = c("Up-regulated" = "#E41A1C", "Down-regulated" = "#377EB8", "Not Significant" = "#B0B0B0")
  ) +
  geom_hline(yintercept = -log10(p_adj_cut), linetype = "dotted", color = "red", linewidth = 0.6) +
  geom_vline(xintercept = c(-log2fc_cut, log2fc_cut), linetype = "dotted", color = "red", linewidth = 0.6) +
  geom_text_repel(
    data              = top_genes_8w,
    aes(label         = gene_symbol),
    size              = 3.2,
    max.overlaps      = Inf,
    force             = 4,
    force_pull        = 1,
    box.padding       = 0.4,
    point.padding     = 0.3,
    min.segment.length = 0,
    segment.color     = "grey40",
    segment.size      = 0.4
  ) +
  coord_cartesian(ylim = c(0, 30)) +
  labs(
    x        = expression(log[2]~"Fold Change"),
    y        = expression(-log[10]~"(" * italic("Adjusted P-value") * ")"),
    title    = "8w Het vs KO",
    subtitle = expression(italic("Cutoffs: Adjusted P") < 0.05 ~ "|" ~ "|" * log[2] * "FC" * "|" > 0.58)
  ) +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    plot.title        = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle     = element_text(size = 10, hjust = 0.5, color = "grey30"),
    axis.title        = element_text(face = "bold", size = 11),
    axis.text         = element_text(size = 10, color = "black"),
    legend.position   = "none",
    panel.grid.major  = element_line(color = "grey95", linewidth = 0.2),
    panel.grid.minor  = element_blank()
  )

# Export 8w figures
ggsave(file.path(output_dir, "Volcano_8w_Het_vs_KO.pdf"), plot = volcano_8w, width = 6, height = 6, units = "in", dpi = 600)
ggsave(file.path(output_dir, "Volcano_8w_Het_vs_KO.tiff"), plot = volcano_8w, width = 6, height = 6, units = "in", dpi = 600, compression = "lzw")


# ==============================================================================
# PART 2: 14w Het vs KO Volcano Plot (No Legend, Custom Validation Targets)
# ==============================================================================

file_14w <- file.path(data_dir, "deg_results_14w_Het_vs_KO.csv")
deg_raw_14w <- read.csv(file_14w, stringsAsFactors = FALSE)

x_col_14w <- names(deg_raw_14w)[grep("b$|beta|log2fc|foldchange|logfc", names(deg_raw_14w), ignore.case = TRUE)][1]
y_col_14w <- names(deg_raw_14w)[grep("qval|adj|padj|pvalue|p.value", names(deg_raw_14w), ignore.case = TRUE)][1]

deg_data_14w <- data.frame(
  transcript_id = as.character(deg_raw_14w[[1]]),
  gene_symbol   = as.character(deg_raw_14w[[5]]),
  x_val         = as.numeric(deg_raw_14w[[x_col_14w]]),
  y_val         = as.numeric(deg_raw_14w[[y_col_14w]]),
  stringsAsFactors = FALSE
)

deg_data_14w$gene_symbol <- ifelse(
  is.na(deg_data_14w$gene_symbol) | deg_data_14w$gene_symbol == "" | deg_data_14w$gene_symbol == "NA",
  deg_data_14w$transcript_id,
  deg_data_14w$gene_symbol
)

deg_clean_14w <- deg_data_14w %>%
  filter(!is.na(x_val) & !is.na(y_val) & y_val > 0) %>%
  group_by(gene_symbol) %>%
  slice_max(order_by = abs(x_val), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    neg_log10_padj = -log10(y_val),
    Expression = case_when(
      x_val >= log2fc_cut & y_val < p_adj_cut ~ "Up-regulated",
      x_val <= -log2fc_cut & y_val < p_adj_cut ~ "Down-regulated",
      TRUE ~ "Not Significant"
    )
  )

# Specific validation targets
custom_targets <- c("Pea15a", "Espn", "Bik", "Defb36")

top_auto_14w <- deg_clean_14w %>%
  filter(y_val < p_adj_cut & abs(x_val) >= log2fc_cut) %>%
  arrange(y_val) %>%
  head(30)

target_df_14w <- deg_clean_14w %>%
  filter(tolower(gene_symbol) %in% tolower(custom_targets))

top_genes_14w <- bind_rows(top_auto_14w, target_df_14w) %>%
  distinct(gene_symbol, .keep_all = TRUE)

volcano_14w <- ggplot(deg_clean_14w, aes(x = x_val, y = neg_log10_padj)) +
  geom_point(aes(color = Expression), alpha = 0.65, size = 1.6) +
  scale_color_manual(
    values = c("Up-regulated" = "#E41A1C", "Down-regulated" = "#377EB8", "Not Significant" = "#B0B0B0")
  ) +
  geom_hline(yintercept = -log10(p_adj_cut), linetype = "dotted", color = "red", linewidth = 0.6) +
  geom_vline(xintercept = c(-log2fc_cut, log2fc_cut), linetype = "dotted", color = "red", linewidth = 0.6) +
  geom_text_repel(
    data              = top_genes_14w,
    aes(label         = gene_symbol),
    size              = 3.3,
    fontface          = ifelse(tolower(top_genes_14w$gene_symbol) %in% tolower(custom_targets), "bold.italic", "plain"),
    max.overlaps      = Inf,
    force             = 5,
    force_pull        = 1,
    box.padding       = 0.4,
    point.padding     = 0.3,
    min.segment.length = 0,
    segment.color     = "grey40",
    segment.size      = 0.4
  ) +
  coord_cartesian(ylim = c(0, 30)) +
  labs(
    x        = expression(log[2]~"Fold Change"),
    y        = expression(-log[10]~"(" * italic("Adjusted P-value") * ")"),
    title    = "14w Het vs KO",
    subtitle = expression(italic("Cutoffs: Adjusted P") < 0.05 ~ "|" ~ "|" * log[2] * "FC" * "|" > 0.58)
  ) +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    plot.title        = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle     = element_text(size = 10, hjust = 0.5, color = "grey30"),
    axis.title        = element_text(face = "bold", size = 11),
    axis.text         = element_text(size = 10, color = "black"),
    legend.position   = "none",
    panel.grid.major  = element_line(color = "grey95", linewidth = 0.2),
    panel.grid.minor  = element_blank()
  )

# Export 14w figures
ggsave(file.path(output_dir, "Volcano_14w_Het_vs_KO.pdf"), plot = volcano_14w, width = 6, height = 6, units = "in", dpi = 600)
ggsave(file.path(output_dir, "Volcano_14w_Het_vs_KO.tiff"), plot = volcano_14w, width = 6, height = 6, units = "in", dpi = 600, compression = "lzw")

message("All volcano plots generated successfully without legends!")
