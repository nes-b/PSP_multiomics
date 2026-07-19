# ---
# title: "PSP snRNA x CSF Integration"
# project: "snRNA+ATAC-seq in PSP frontal cortex"
# author: "Nils Briel, Center for Neuropathology, LMU Munich & Dept. Neurology, USZ Zurich"
# date: "2026-04-06"
# ---

library(limma)
library(stringr)
library(zoo)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggbeeswarm)
library(ggpubr)
library(readxl)
library(caret)
library(pheatmap)
library(stringr)
library(ComplexHeatmap)
library(tidyplots)
library(cowplot)

# intersect snRNA w/ CSF proteomics
csf_raw <- read_xlsx(
  "/Users/nils/Library/Mobile Documents/com~apple~CloudDocs/Research/01_multi_omics_tau/results/csf_integration/12014_2024_9507_MOESM4_ESM.xlsx",
  sheet = 1, col_names = FALSE
)

# data tidying
batch_row  <- csf_raw[1, ]
sample_row <- csf_raw[2, ]
data_block <- csf_raw[-c(1, 2), ]
orig_names <- as.character(unlist(sample_row))
uniq_names <- make.unique(orig_names)
colnames(data_block) <- uniq_names
gene_col_idx <- ncol(data_block)
expr_wide <- data_block[, -gene_col_idx]
expr_wide <- mutate_all(expr_wide, as.numeric)
gene_vec <- data_block[[gene_col_idx]] %>% as.character()
expr_mat <- t(as.matrix(expr_wide))
colnames(expr_mat) <- gene_vec
expr_mat <- t(expr_mat)
sample_ids <- uniq_names[-gene_col_idx]
colnames(expr_mat) <- sample_ids
rownames(expr_mat) <- gene_vec
all_sample_idx <- 1:length(sample_ids)
batch_vec <- na.locf(as.character(unlist(batch_row))[all_sample_idx], fromLast = FALSE)
diag_vec <- str_extract(orig_names[all_sample_idx], "(HC|PSP|PD|QC)$")
keep_idx <- which(!diag_vec %in% c("PD", "QC"))

csf_coldata <- tibble(
  sample = sample_ids[keep_idx],
  batch = factor(batch_vec[keep_idx]),
  diagnosis = factor(diag_vec[keep_idx])
)

expr_keep_idx <- keep_idx
expr_mat_subset <- expr_mat[, expr_keep_idx, drop = FALSE]

colnames(expr_mat_subset) <- csf_coldata$sample

valid_idx <- !is.na(csf_coldata$batch)
csf_coldata <- csf_coldata[valid_idx, ]
expr_mat_subset <- expr_mat_subset[, valid_idx, drop = FALSE]
design <- model.matrix(~ batch + diagnosis, data = csf_coldata)

# Diff abundance analysis
fit <- lmFit(expr_mat_subset, design)
fit <- eBayes(fit)
res <- topTable(fit, coef = "diagnosisPSP", number = Inf, sort.by = "none")
csf_limma <- as_tibble(res, rownames = "target") %>%
  dplyr::rename(log2FC = logFC) 
  # filter(!is.na(log2FC), P.Value < 0.05)

## snRNAseq DEGs
deg <- read_excel("/Users/nils/Library/Mobile Documents/com~apple~CloudDocs/Research/01_multi_omics_tau/doc/extended_data_tables.xlsx",
                  sheet = 3) %>%
  mutate(celltype   = str_remove_all(celltype, '"'),
         avg_log2FC = as.numeric(as.character(avg_log2FC))) 
  # filter(p_val <0.05)


celltype_cols <- c(
  "Exc-DLN" = "#39A43A",
  "Exc-ULN" = "#F0B44C",
  "Oligo" = "#578B9A",
  "Astro" = "#E15759",
  "OPC" = "#2F84A0FF",
  "Inh-Neu" = "#7D486EFF",
  "Endo-VLMC" = "#9C6B43",
  "Micro-PVM" = "#5C5C5C"
)

merged <- deg %>%
  filter(gene_name %in% csf_limma$ID) %>%
  inner_join(
    csf_limma %>% select(ID, log2FC_csf = log2FC, padj_csf = adj.P.Val, pval_csf = P.Value),
    by = c("gene_name" = "ID")
  ) %>%
  mutate(
    conc = case_when(
      avg_log2FC > 0 & log2FC_csf > 0 ~ "Up/Up",
      avg_log2FC < 0 & log2FC_csf < 0 ~ "Down/Down",
      TRUE ~ "Discordant"
    ),
    significance = case_when(
      p_val < 0.05 & pval_csf < 0.05 ~ "Significant",
      TRUE ~ "Non-significant"
    ),
    mean_fc = (avg_log2FC + log2FC_csf) / 2,
    diff_fc = log2FC_csf - avg_log2FC
  )

p1 <- ggplot(merged, aes(x = avg_log2FC, y = log2FC_csf, fill = conc, label = gene_name)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
  geom_hline(yintercept = 0, linetype = 3, color = "grey70") +
  geom_vline(xintercept = 0, linetype = 3, color = "grey70") +
  geom_point(aes(alpha = significance, color = significance), 
             size = 1.5, shape = 21) +
  geom_text_repel(max.overlaps = 50, size = 3, show.legend = FALSE, fontface = "italic") +
    geom_smooth(aes(fill = NULL, color = NULL), method = "lm", se = TRUE, 
              color = "black", linewidth = 0.5) +
  stat_cor(aes(color = NULL, fill = NULL), method = "pearson", 
           label.x.npc = "left", label.y.npc = "top", size = 3) +
  scale_fill_manual(values = c(
    "Up/Up" = "#D73027",
    "Down/Down" = "#4575B4",
    "Discordant" = "grey30"
  )) +
  scale_alpha_manual(values = c(
    "Significant" = 0.8,      
    "Non-significant" = 0.3   # More transparent
  )) +
  scale_color_manual(values = c(
    "Significant" = "black",
    "Non-significant" = "grey70"  # Lighter stroke
  )) +
  facet_wrap(~celltype, ncol = 2) +
  theme_bw() + theme(legend.position = "bottom") +
  labs(
    x = "Brain snRNAseq log2FC",
    y = "CSF proteomics log2FC",
    fill = "Concordance",
    alpha = NULL,
    color = NULL
  )


# Filter for concordant genes
concordant_genes <- merged %>%
  filter(p_val < 0.05 & pval_csf < 0.05, !conc == "Discordant", abs(mean_fc) >0.5) %>% 
  mutate(neg_log10_pval = -log(as.numeric(p_val)),
         gene_name = reorder(gene_name, mean_fc))

# Create dotplot
p2 <- ggplot(concordant_genes, aes(x = gene_name, y = celltype)) +
  geom_point(shape = 21, aes(fill = mean_fc, size = abs(neg_log10_pval))) +
  scale_fill_gradient2(name = "Mean log2(FC)",low = "blue",mid = "white",high = "red",midpoint = 0)  +
  scale_size_continuous(name = "-log10(P)", range = c(1, 6)) +  
  labs(x = "Concordant genes", y = "Cell types") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.text.y = element_text(face = "italic")) + 
  coord_flip()

p3 <- merged %>%
  group_by(celltype) %>% filter(p_val < 0.05 & pval_csf < 0.05,  abs(mean_fc) >0.5) %>%
  summarise(
    concordant = mean(conc != "Discordant", na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>% mutate(celltype = fct_reorder(celltype, -concordant)) %>% 
  ggplot(aes(x = celltype, y = concordant, fill = celltype)) +
  geom_col(width = 0.5, color = "grey4") +
  scale_fill_manual(values = celltype_cols, guide = "none") +
  theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  labs(y = "# concordant sig. genes / sig. genes", x = NULL) 

# combine
p23fig <- plot_grid(plotlist = list(p3, p2), labels = c("b", "c"), ncol = 1, align = "hv", rel_heights = c(0.6, 1))
fig <- plot_grid(plotlist = list(p1, p23fig), labels = c("a"), ncol = 2, align = "hv", rel_widths = c(1.5, 1))
ggsave("output_plot_grid_ext_.pdf", fig, width = 12, height = 10)
ggsave("output_plot_grid_ext_.png", fig, width = 12, height = 10)

sessionInfo()