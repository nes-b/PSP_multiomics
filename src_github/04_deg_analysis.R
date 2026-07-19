# ---
# title: "snRNA+ATAC-seq in PSP frontal cortex"
# author: "Nils Briel, Center for Neuropathology, LMU Munich & Dept. Neurology, USZ Zurich"
# date: "Created: 2024-03-13"
# ---

# 1 Setup Environment ####
# First, we need to load the necessary packages:
set.seed(1234)
library(Seurat)
options(Seurat.object.assay.version = "v5")
library(SeuratData)
library(Signac)

library(data.table)
library(Matrix)
library(parallel)
library(tidyverse)
library(RhpcBLASctl)
library(unix)
library(SCpubr)
source('../src/utils.R')
rpath <- "../results/psp/"

blas_set_num_threads(24)
future::plan("multicore", workers = 16)

options(future.globals.maxSize= 24000*1024^2) #24 000 MB
rlimit_as(150e9)

coembed <- readRDS("../results/psp/coembed.Rds")
combrna <- readRDS("../results/psp/combrna_psp_scl2.Rds") 

# filter out non-coding / not annotated genes:
unannotated_transcripts <- c(grep(c("^AL\\d+\\.\\d+"), rownames(combrna), value = TRUE),
                             grep(c("^AC\\d+\\.\\d+"), rownames(combrna), value = TRUE),
                             grep(c("^AD\\d+\\.\\d+"), rownames(combrna), value = TRUE),
                             grep(c("^AF\\d+\\.\\d+"), rownames(combrna), value = TRUE),
                             grep(c("^AJ\\d+\\.\\d+"), rownames(combrna), value = TRUE),
                             grep(c("^AP\\d+\\.\\d+"), rownames(combrna), value = TRUE))
combrna <- combrna[!rownames(combrna) %in% unannotated_transcripts, ]

# 2 GTSummary ####
library(dplyr)
library(gtsummary)
library(flextable)
library(conflicted)

sel_meta <- c("Pub_ID","Age","NPdiagnosis","pmi","gender","tech","brain.weight","Synuclein","LBs..McKeiths.",
              "LBs..Braak.","APOE","Braak.Braak..NFT.","Thal.Phase","TDP.43","FUS","nCount_RNA","nCount_peaks")
sel_f <- c("Pub_ID","NPdiagnosis", "gender", "tech", "Synuclein", "LBs..McKeiths.", "LBs..Braak.", "APOE", "TDP.43", "FUS")
sel_num <- base::setdiff(sel_meta, sel_f)

meta_df <- coembed@meta.data %>%
  dplyr::select(all_of(sel_meta)) %>%
  mutate(across(all_of(sel_f), as.factor),
         across(all_of(sel_num), as.numeric))

summary_df <- meta_df %>% 
  group_by(Pub_ID,tech) %>% 
  summarise_if(is.factor, ~names(which.max(table(.)))) %>%
  mutate(id_tech = paste0(Pub_ID,"_",tech))

# ugly fix of NA in combatac meta.data PSP NPdiagnosis column
summary_df$NPdiagnosis[summary_df$NPdiagnosis == "Ctrl" & grepl("^PSP", summary_df$Pub_ID)] <- "PSP"

summary_df_2 <- meta_df %>% 
  group_by(Pub_ID,tech) %>% add_count() %>%
  summarise_if(is.numeric, list(mean = mean)) %>%
  mutate(id_tech = paste0(Pub_ID,"_",tech)) 

meta_df <- left_join(summary_df,summary_df_2[,-c(1,2)], by= "id_tech") %>%
  distinct(id_tech, .keep_all = TRUE) %>%
  ungroup() %>%
  mutate(APOE = ifelse(Pub_ID %in% c("PSP1", "CBD2", "CBD3", "CBD7", "CBD8", "C1"), NA, APOE))

sel_num <- base::setdiff(colnames(meta_df), sel_f)
meta_df <- meta_df %>%
  mutate(across(all_of(sel_num), as.numeric),
         APOE = ifelse(APOE %in% c("Unknown","not specified","unassessed"), NA, APOE),
         FUS = ifelse(FUS %in% c("Unknown","not specified","unassessed"), NA, FUS))

conflicted::conflicts_prefer(gtsummary::select)
meta_df %>%
  select("NPdiagnosis","gender","Age_mean","pmi_mean","APOE","Braak.Braak..NFT._mean",
         "Thal.Phase_mean","TDP.43","FUS","nCount_peaks_mean","nCount_RNA_mean","n_mean","tech") %>%
  mutate(tech = paste0("Tech: sn", toupper(tech), "-seq"),
         NPdiagnosis = as.factor(NPdiagnosis),
         tech = as.factor(tech)) %>%
  tbl_strata(
    strata = tech,
    .tbl_fun =
      ~ .x %>%
      tbl_summary(by = NPdiagnosis, missing = "no",
                  statistic = list(all_continuous() ~ "{median}\n[{p25}, {p75}]",
                                   all_categorical() ~ "{n}")) %>%
      add_n(),
    .header = "**{strata}**, N = {n}") %>%
  bold_labels() %>%
  as_flex_table() %>%
  save_as_docx("metadata_table", path = paste0(rpath,"metadata_table.docx"))


###
Idents(combrna) <- "celltype"
pseudo <-
  AggregateExpression(
    combrna,
    assays = "RNA",
    return.seurat = T,
    group.by = c("NPdiagnosis","Pub_ID", "celltype")
  )
pseudo@meta.data <-
  pseudo@meta.data %>%
  as_tibble(rownames = "id_div") %>%
  separate(col = id_div, into = c("Condition", "Pub_ID", "celltype"), sep = "_", remove = FALSE) %>%
  dplyr::select(-"orig.ident") %>%
  as.data.frame() %>%
  mutate(celltype.NPdiagnosis = paste0(.$celltype,"_",.$Condition)) %>%
  `rownames<-`(.$id_div)

t_load <- readxl::read_excel("../data/tau_load.xlsx", sheet = "Sheet3") %>% 
  .[,c("Pub_ID","TA_load","NFT_load", "CB_load", "Threads_load")]
pseudo@meta.data <- left_join(pseudo@meta.data, t_load, by = 'Pub_ID') %>%
t_load <- left_join(t_load, subset(pseudo@meta.data, NPdiagnosis == "PSP")[,c("Pub_ID", "Age", "pmi", "Braak&Braak (NFT)", "Thal-Phase")], 
                    by = "Pub_ID") %>%
  column_to_rownames(., "Pub_ID") %>% data.matrix()

t_load2 <- data.frame(t_load) %>%
  mutate(Pub_ID = rownames(t_load)) %>% 
  pivot_longer(cols = c("TA_load","NFT_load", "CB_load", "Threads_load"),
               names_to = "NPtrait", values_to = "Value")

p1 <- ggpubr::ggscatter(data.frame(t_load), x = "CB_load", y = "TA_load", cor.coef = T, cor.method = "spearman", add = "reg.line") + 
  theme_bw()
p2 <- ggpubr::ggscatter(data.frame(t_load), x = "NFT_load", y = "Threads_load", cor.coef = T, cor.method = "spearman", add = "reg.line") + 
  theme_bw()
pdf("../../corr_nptrait_scatter.pdf", height = 4.5, width = 4)
p1|p2
dev.off()

ggplot(t_load2, aes(x = NPtrait, y = Value, fill = NPtrait)) +
  geom_bar(position = "dodge", stat = "identity", color="grey42") + facet_grid(.~Pub_ID) +
  theme_bw()
ggsave("../../barplot_nptrait.pdf",width = 8,height = 2)

pdf("../../demo_heatmap_case_nptrait.pdf", height = 4.5, width = 4)
library(ComplexHeatmap)
ComplexHeatmap::pheatmap(column_to_rownames(t_load, "Pub_ID"), 
         cluster_cols = T,
         cluster_rows = F,
         scale = "no",
         clustering_method = "complete",
         show_colnames = T, show_rownames = T,
         colorRampPalette(c("white", "#FEF12B","#FF6527", "red"))(4))
dev.off()


# 3 DEG/DAG/DTFME Analysis ####
library(ggrepel)
# We compute DEGs in pseudobulks of entire samples
# https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/pseudobulk_DESeq2_scrnaseq.md
# Create directories to save results if they don't already exist:
if (!dir.exists("../results/psp/DESeq2")) {
  dir.create("../results/psp/DESeq2")
}
if (!dir.exists("../results/psp/DESeq2/celltype")) {
  dir.create("../results/psp/DESeq2/celltype")
}
setwd("../results/psp/DESeq2/celltype/")

## 3.1 celltype: pseudobulk DEG ####
Idents(combrna) <- "celltype"
pseudo <-
  AggregateExpression(
    combrna,
    assays = "RNA",
    return.seurat = T,
    group.by = c("NPdiagnosis","Pub_ID", "celltype")
  )
pseudo@meta.data <-
  pseudo@meta.data %>%
  as_tibble(rownames = "id_div") %>%
  separate(col = id_div, into = c("Condition", "Pub_ID", "celltype"), sep = "_", remove = FALSE) %>%
  dplyr::select(-"orig.ident") %>%
  as.data.frame() %>%
  mutate(celltype.NPdiagnosis = paste0(.$celltype,"_",.$Condition)) %>%
  `rownames<-`(.$id_div)

Idents(pseudo) <- "celltype.NPdiagnosis"
plot_list <- list()
for (i in unique(pseudo$celltype)) {
  cat(i,"\n")
  try({
    ident1 <- paste0(i, "_PSP")
    ident2 <- paste0(i, "_Ctrl")
    condition.diffgenes <-
      FindMarkers(
        pseudo,
        ident.1 = ident1,
        ident.2 = ident2,
        pseudocount.use = 0.1,
        min.pct = 0.01,
        logfc.threshold = 0.1, 
        test.use = "DESeq2"
      )
   condition.diffgenes$gene <- rownames(condition.diffgenes)
   thresh_bh <- max(condition.diffgenes$p_val[condition.diffgenes$p_val_adj <= 0.1], na.rm = TRUE)
   
   plot_list[[i]] <- volcano(res_tbl = condition.diffgenes,
                             thresh_plot_1 = 0.1,
                             thresh_plot_2 = thresh_bh,
                             log2fc = 0.5,
                             top_x_to_plot = 25,
                             title = paste0(i, " in PSP vs. Ctrl"),
                             padj_col = "p_val",
                             log2FoldChange_col = "avg_log2FC",
                             ext_y = 1,
                             text_size = 3,
                             cols = c("#B53737", "#4774a8", "grey"))
   plot_list[[i]]
   ggsave(paste0("volc_",i,".pdf"), width = 9, height = 5)
   
    write.csv(condition.diffgenes, file = paste0(i, ".csv"))
  })
}

# Combine into panel
combined_plot <- cowplot::plot_grid(plotlist = plot_list, ncol = 2, labels = "AUTO")  # 4 rows x 2 cols for 8 plots
ggsave("volcano_panel_all_celltypes.pdf", combined_plot, width = 17, height = 18)

# bidirect_barplot
df <- lapply(list.files(pattern = "*.csv"), read.csv)
names(df) <-  str_remove(list.files(pattern = "*.csv"),".csv")
df <- imap(df, ~ mutate(.x, Category = .y)) %>% bind_rows(.)
write.csv(df, 'ct_combined_data.csv', row.names = FALSE)

bidirect_plot(df = df, ct = "celltype", multiple_p_vals = F, padj_tresh = 0.1) + scale_fill_manual(values = ct_cols_hyph)
ggsave("bidirect_barplot_celltype.pdf", width = 5, height = 3.5)

df %>% filter(p_val_adj < 0.1 & Category == "ct_combined_data") %>% nrow() #479
df %>% filter(p_val_adj < 0.1 & Category == "Exc-DLN") %>% nrow() #128
df %>% filter(p_val_adj < 0.1 & Category == "Exc-ULN") %>% nrow() #37
df %>% filter(p_val_adj < 0.1 & Category == "Astro") %>% nrow() #24
df %>% filter(p_val_adj < 0.1 & Category == "Oligo") %>% nrow() #28
df %>% filter(p_val_adj < 0.1 & Category == "OPC") %>% nrow() #11
df %>% filter(p_val_adj < 0.1 & Category == "Inh-Neu") %>% nrow() #7
df %>% filter(p_val_adj < 0.1 & Category == "Micro-PVM") %>% nrow() #1

# DESeq2
deseq_ls <- prep_deseq2(combrna = combrna, 
                        feature_sample = "Pub_ID",
                        feature_group = "NPdiagnosis",
                        feature_cell_identity = "celltype",
                        cores = 32)
p1 <- DESeq2::plotPCA(deseq_ls$rld, ntop = 1000, intgroup = "celltype") + 
  scale_color_manual(values = ct_cols) +  theme_bw() 
p2 <- DESeq2::plotPCA(deseq_ls$rld, ntop = 1000, intgroup = "NPdiagnosis") + 
  scale_color_manual(values = group_cols) +  theme_bw() 

plot_grid(plotlist = list(p1,p2), align = "hv")
ggsave("pca_celltype.pdf", width = 9, height = 6)

rld_mat <- assay(deseq_ls$rld)
rld_cor <- stats::cor(rld_mat)

# Plot heatmapdff
# Select cell type of interest
ann_col <- data.frame(row.names = deseq_ls$meta$clus_np_id, 
                      NPdiagnosis = deseq_ls$meta$NPdiagnosis,
                      Celltype = deseq_ls$meta$celltype)
pdf("cor_heatmap_celltype.pdf", height = 6, width = 8)
pheatmap(rld_cor, annotation_col = ann_col, 
         cluster_cols = T, 
         clustering_method = "ward.D2",
         annotation_colors = list(NPdiagnosis = group_cols,
                                  Celltype = ct_cols), 
         show_colnames = F, show_rownames = F,
         colorRampPalette(c("#4774a8", "white", "#B53737"))(50))
dev.off()

## 3.2 predicted.subclass: pseudobulk DEG ####
if (!dir.exists("../predicted.subclass")) {
  dir.create("../predicted.subclass")
}
setwd("../predicted.subclass/")
combrna$predicted.subclass <-
  stringr::str_replace_all(combrna$predicted.subclass,
                           pattern = "[ /]",
                           replacement = "_")
Idents(combrna) <- "predicted.subclass"
pseudo <-
  AggregateExpression(
    combrna,
    assays = "RNA",
    return.seurat = T,
    group.by = c("NPdiagnosis", "Pub_ID", "predicted.subclass")
  )
pseudo@meta.data <-
  pseudo@meta.data %>%
  as_tibble(rownames = "id_div") %>%
  separate(col = id_div, into = c("Condition", "Pub_ID", "predicted.subclass"), sep = "_", remove = FALSE) %>%
  dplyr::select(-"orig.ident") %>%
  as.data.frame() %>% 
  mutate(celltype.NPdiagnosis = paste0(.$predicted.subclass,"_",.$Condition)) %>%
  `rownames<-`(.$id_div)

Idents(pseudo) <- "celltype.NPdiagnosis"
for (i in unique(pseudo$predicted.subclass)) {
  cat(i,"\n")
  try({
    ident1 <- paste0(i, "_PSP")
    ident2 <- paste0(i, "_Ctrl")
    condition.diffgenes <-
      FindMarkers(
        pseudo,
        ident.1 = ident1,
        ident.2 = ident2,
        pseudocount.use = 0.1,
        min.pct = 0.01,
        logfc.threshold = 0.1,
        test.use = "DESeq2"
      )
    condition.diffgenes$gene <- rownames(condition.diffgenes)
    volcano(res_tbl = condition.diffgenes,
            padj_thresh = 0.05,
            padj_thresh_plot = 5e-5,
            log2fc = 0.5,
            top_x_to_plot = 15,
            title = paste0(i," in PSP vs. Ctrl"),
            padj_col = "p_val",
            log2FoldChange_col = "avg_log2FC",
            cols = c("#B53737","#4774a8","grey")
    )
    ggsave(paste0("volc_",i,".pdf"), width = 12, height = 5)
    
    write.csv(condition.diffgenes, file = paste0(i, ".csv"))
  })
}

# bidirect_barplot
df <- lapply(list.files(pattern = "*.csv"), read.csv)
names(df) <-  str_remove(list.files(pattern = "*.csv"),".csv")
df <- imap(df, ~ mutate(.x, Category = .y)) %>% bind_rows(.)

bidirect_plot(df = df, ct = "predicted.subclass")  + scale_fill_manual(values = sclass_cols)
ggsave("bidirect_barplot_subclass.pdf", width = 7, height = 5)

# DESeq2
deseq_ls <- prep_deseq2(combrna = combrna, 
                        feature_sample = "Pub_ID",
                        feature_group = "NPdiagnosis",
                        feature_cell_identity = "predicted.subclass",
                        cores = 32)

p1 <- DESeq2::plotPCA(deseq_ls$rld, ntop = 1000, intgroup = "celltype") + 
  scale_color_manual(values = sclass_cols) +  theme_bw() + coord_equal(1.5)
p2 <- DESeq2::plotPCA(deseq_ls$rld, ntop = 1000, intgroup = "NPdiagnosis") + 
  scale_color_manual(values = group_cols) +  theme_bw() + coord_equal(1.5)
plot_grid(plotlist = list(p1,p2), align = "hv")
ggsave("pca_subclass.pdf", width = 11, height = 9)

rld_mat <- assay(deseq_ls$rld)
rld_cor <- stats::cor(rld_mat)

# Plot heatmapdff
# Select cell type of interest
ann_col <- data.frame(row.names = deseq_ls$meta$clus_np_id, 
                      NPdiagnosis = deseq_ls$meta$NPdiagnosis,
                      Subclass = deseq_ls$meta$celltype)

names(sclass_cols)[13] <- "Sst_Chodl.NA.NA"

pdf("cor_heatmap_subclass.pdf", height = 6, width = 8)
pheatmap(rld_cor, annotation_col = ann_col, 
         cluster_cols = T, 
         clustering_method = "ward.D2",
         annotation_colors = list(NPdiagnosis = group_cols,
                                  Subclass = sclass_cols), 
         show_colnames = F, show_rownames = F,
         colorRampPalette(c("#4774a8", "white", "#B53737"))(50))
dev.off()

## 3.3 subcluster: pseudobulk DEG ####
if (!dir.exists("../subcluster")) {
  dir.create("../subcluster")
}
setwd("../subcluster/")
combrna$subcluster <- stringr::str_replace_all(combrna$subcluster,
                                               pattern = "[ /]",
                                               replacement = "_") %>% as.factor()
combrna$subcluster_hyph <-
  stringr::str_replace_all(combrna$subcluster,
                           pattern = "_",
                           replacement = "-") %>% as.factor()
Idents(combrna) <- "subcluster"
DefaultAssay(combrna) <- "RNA"

for (ident1 in levels(combrna$subcluster)) {
  cat(i,"\n")
  try({
  cls <- levels(combrna$subcluster)
  ct <- subset(combrna@meta.data, subcluster == ident1)$celltype[1]
  cls <- cls[grep(ct, cls)]
  ident2 <- cls[!cls == ident1]
  ident2 <- if (length(ident2) < 1) NULL else ident2
  
  condition.diffgenes <-
    FindMarkers(
      combrna,
      ident.1 = ident1,
      ident.2 = ident2,
      pseudocount.use = 0.1,
      min.pct = 0.25,
      logfc.threshold = 0.25,
      test.use = "MAST"
    )
  
  write.csv(condition.diffgenes, file = paste0(ident1, ".csv"))
  })
}

# plot gene expression over CB load
sbs <- subset(combrna, subset = celltype =="Oligo")
mat <- sbs[c("FKBP5","ZBTB16","MAPT","NC3C2"),] 
mat <- mat@assays$RNA$data %>% t %>% data.frame() 
mat$CB_load <- ifelse(is.na(sbs$CB_load), 0,  sbs$CB_load)
mat <- mat %>% pivot_longer(cols = c("FKBP5","ZBTB16","MAPT"), values_to = "SCT", names_to = "Gene")
mat$Gene <- as.factor(mat$Gene)
ggplot(mat, aes(x = CB_load, y = SCT, group = Gene, color = Gene)) + 
  stat_smooth(method = "loess") + theme_bw()
ggsave("fkbp5_zbtb16_mapt_over_CB_load.pdf")


# 4 Gene Ontology / KEGG ####
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(ggupset)
library(pathview)
library(stringi)


## 4.1 celltype ####
if (!dir.exists("../celltype_GO")) {
  dir.create("../celltype_GO")
}
setwd("../celltype_GO")
if (!dir.exists("plots")) {
  dir.create("plots")
}


for(i in list.files("../celltype/", pattern = ".csv")){
  cat(i)
    try(lollipop_go(i, lev="celltype"))
}

csv_files <- list.files(pattern = '\\joined.csv$')
data_list <- list()
for (file in csv_files) {
  df <- read.csv(file)
  df$celltype_id <- tools::file_path_sans_ext(file)
  data_list[[file]] <- df
}

## pack into funciton --> utils.R
combined_data <- bind_rows(data_list) %>%
  mutate(celltype_id = gsub("joined", "", .$celltype_id)) %>% 
  separate(col = GeneRatio,sep = "/",into = c("nom","denom"),remove = F,convert = T) %>%
  mutate(`Gene ratio` = .$nom / .$denom) %>%
  filter(celltype_id %in% c("Exc_DLN","Exc_ULN", "Astro", "Oligo", "OPC", "Inh_Neu","Micro_PVM")) %>%
  filter(p.adjust < 0.05) %>%
  group_by(celltype_id,ONTOLOGY,direct) %>% 
  slice_min(order_by = p.adjust, n = 3) %>% # adjust top vars
  ungroup() %>%
  mutate_at(vars(ONTOLOGY, ID, Description, geneID), as.factor) %>%
  mutate_at(vars(GeneRatio, BgRatio), as.factor) %>%
  mutate_at(vars(pvalue, p.adjust, qvalue, Count,`Gene ratio`), as.numeric) %>%
  mutate(`Gene ratio` = ifelse(direct == "up", `Gene ratio`, -`Gene ratio`),
         celltype_id = factor(celltype_id, levels = c("Exc_DLN","Exc_ULN", "Astro", "Oligo", "OPC", "Inh_Neu","Micro_PVM")),
         ONTOLOGY = factor(ONTOLOGY, levels = c("BP","CC","MF","KEGG"))) 

p <- ggplot(combined_data, aes(x = celltype_id, 
                     y = tidytext::reorder_within(Description, `Gene ratio`, list(ONTOLOGY)), 
                     fill = `Gene ratio`)) + 
  geom_tile(color = "grey4") + 
  facet_grid(ONTOLOGY~., drop = T, space = "free", scales = "free") +
  theme_test() +
  scale_fill_gradient2(
    low = "#4774a8", 
    mid = "white", 
    high = "#B53737", 
    midpoint = 0
  ) +
  labs(title = "",
       y = "",
       x = "") + 
  scale_y_discrete(limits=rev) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

p
ggsave("heatmap_cts_unselected.pdf", width = 7, height = 14)


## 4.2 subclass ####
if (!dir.exists("../subclass_GO")) {
  dir.create("../subclass_GO")
}
setwd("../subclass_GO")
if (!dir.exists("plots")) {
  dir.create("plots")
}

# generate GO for subcluster
conflicts_prefer(base::unname)
for(i in list.files("../predicted.subclass/", pattern = ".csv")){
  cat(i)
  try(lollipop_go(i, lev="predicted.subclass"))
}


## 4.3 subcluster ####
if (!dir.exists("../subcluster_GO")) {
  dir.create("../subcluster_GO")
}
setwd("../subcluster_GO")
if (!dir.exists("plots")) {
  dir.create("plots")
}

# generate GO for subclass
for(i in list.files("../subcluster/", pattern = ".csv")){
  cat(i)
  try(lollipop_go(i, lev="subcluster"))
}


# plot
fac_type <- plot_go_facet(df, facet_var = "celltype")
ggsave("plots/GO_subclusters_facet_type_astro.pdf",plot = fac_type$Astro, height = 7, width = 4.5)
ggsave("plots/GO_subclusters_facet_type_Endo_VLMC.pdf",plot = fac_type$Endo_VLMC, height = 8, width = 5.8)
ggsave("plots/GO_subclusters_facet_type_Exc_ULN.pdf",plot = fac_type$Exc_ULN, height = 3.5, width = 4)
ggsave("plots/GO_subclusters_facet_type_Exc_DLN.pdf",plot = fac_type$Exc_DLN, height = 10, width = 9)
ggsave("plots/GO_subclusters_facet_type_Inh_Neu.pdf",plot = fac_type$Inh_Neu, height = 5, width = 5)
ggsave("plots/GO_subclusters_facet_type_Micro_PVM.pdf",plot = fac_type$Micro_PVM, height = 3.5, width = 4)
ggsave("plots/GO_subclusters_facet_type_OPC.pdf",plot = fac_type$OPC, height = 4, width = 4.5)
ggsave("plots/GO_subclusters_facet_type_Oligo.pdf",plot = fac_type$Oligo, height = 6, width = 5.5)
fac_subclass <- plot_go_facet(df, facet_var = "subclass")
ggsave("plots/GO_subclusters_facet_subclass_lamp5.pdf",plot = fac_subclass$Lamp5, height = 4, width = 4.5)
ggsave("plots/GO_subclusters_facet_subclass_Pvalb.pdf",plot = fac_subclass$Pvalb, height = 4, width = 6)
ggsave("plots/GO_subclusters_facet_subclass_Sst.pdf",plot = fac_subclass$Sst, height = 6, width = 4.5)
ggsave("plots/GO_subclusters_facet_subclass_Vip.pdf",plot = fac_subclass$Vip, height = 6, width = 4.5)

# plot subclasses
ggsave("plots/go_dir_astro.pdf",plot = plot_ls$Astro, height = 10, width = 12)
ggsave("plots/go_dir_endo_vlmc.pdf",plot = plot_ls$Endo_VLMC, height = 9.5, width = 14)
ggsave("plots/go_dir_l23it.pdf",plot = plot_ls$`L2/3 IT`, height = 3.5, width = 12)
ggsave("plots/go_dir_l56np.pdf",plot = plot_ls$`L5/6 NP`, height = 3.5, width = 12)
ggsave("plots/go_dir_l5et.pdf",plot = plot_ls$`L5 ET`, height = 3.5, width = 12)
ggsave("plots/go_dir_l5it.pdf",plot = plot_ls$`L5 IT`, height = 20, width = 12)
ggsave("plots/go_dir_l6ct.pdf",plot = plot_ls$`L6 CT`, height = 3.5, width = 12)
ggsave("plots/go_dir_l6it.pdf",plot = plot_ls$`L6 IT`, height = 6, width = 12)
ggsave("plots/go_dir_l6itcar3.pdf",plot = plot_ls$`L6 IT Car3`, height = 3.5, width = 12)
ggsave("plots/go_dir_l6b.pdf",plot = plot_ls$L6b, height = 7, width = 12)
ggsave("plots/go_dir_mic.pdf",plot = plot_ls$`Micro-PVM`, height = 3.5, width = 12)
ggsave("plots/go_dir_opc.pdf",plot = plot_ls$OPC, height = 10, width = 12)
ggsave("plots/go_dir_oligo.pdf",plot = plot_ls$Oligo, height = 9, width = 12)
ggsave("plots/go_dir_pv.pdf",plot = plot_ls$Pvalb, height = 6, width = 12)
ggsave("plots/go_dir_sncg.pdf",plot = plot_ls$Sncg, height = 3.5, width = 12)
ggsave("plots/go_dir_sst.pdf",plot = plot_ls$Sst, height = 11, width = 12)
ggsave("plots/go_dir_sstchodl.pdf",plot = plot_ls$Sst_Chodl, height = 3.5, width = 12)
ggsave("plots/go_dir_vip.pdf",plot = plot_ls$Vip, height = 12, width = 12)


plot_ls <- list()
for(i in list.files()[grep(c("up_enrichGO.Rds"),list.files())]){
  egoup <- readRDS(i) %>% simplify()
  ## barplot of all 18 significant pathways
  plot_ls[[i]] <- barplot(egoup, x = "GeneRatio", font.size = 10,
                          showCategory = 10) + theme_bw() + coord_equal(0.05)
  pdf(paste0("plots/",i,"Barplot.pdf"), height = 5, width = 5)
  plot(plot_ls[[i]])
  ylab(label = "Gene ratio")
  dev.off()
  xx <- pairwise_termsim(egoup)
    ## Emapplot
  pdf(paste0("plots/",i,"Emapplot.pdf"), height = 16, width = 16)
  plot(emapplot(xx))
  dev.off()
}

cowplot::plot_grid(plotlist = plot_ls, ncol = 5, labels = names(plot_ls))


# 5 Motif Analysis ####
library(JASPAR2022) # if loading error: devtools::install_version("dbplyr", version = "2.3.4"), appeared on 2024-01-07
library(TFBSTools)
library(patchwork)
library(BiocParallel)
register(MulticoreParam(32))
set.seed(1234)

coembed <- readRDS("~/multiome_psp_cbd/results/psp/coembed.Rds")
combatac <- readRDS("~/multiome_psp_cbd/results/psp/combatac_psp.Rds")
DefaultAssay(combatac) <- "chromvar"
Idents(combatac) <- combatac$celltype

# annotate rownames
temp <- combatac@assays$chromvar@data
motifs_names <- sapply(row.names(temp), function(x) {name(getMatrixByID(JASPAR2022, ID = x))})
row.names(coembed@assays$chromvar@data) <- unname(motifs_names)

## 5.1 celltype
setwd("~/multiome_psp_cbd/results/psp/DESeq2/celltype")
if (!dir.exists("ct_motif")) {
  dir.create("ct_motif")
}
setwd("ct_motif")

combatac$NPdiagnosis <- ifelse(grepl("^P", combatac$Pub_ID), "PSP", "Ctrl")
combatac$celltype.NPdiagnosis = paste0(combatac$celltype,"_",combatac$NPdiagnosis)

Idents(combatac) <- combatac$celltype.NPdiagnosis
plot_list = list()
for (i in levels(combatac$celltype)) {
  cat(i,"\n")
  try({
    ident1 <- paste0(i, "_PSP")
    ident2 <- paste0(i, "_Ctrl")
    differential.activity <- FindMarkers(
      object = combatac,
      ident.1 = ident1,
      ident.2 = ident2,
      only.pos = F,
      test.use = 'LR',
      min.pct = 0.05,
      mean.fxn = rowMeans,
      fc.name = "avg_diff",
      latent.vars = "nCount_peaks"
    )
    
    differential.activity$gene <- paste0(names(motifs_names[rownames(differential.activity)]),"_",motifs_names[rownames(differential.activity)])
    differential.activity$log2FoldChange <- log2(differential.activity$avg_diff+1.5)
    differential.activity$log2FoldChange <- differential.activity$log2FoldChange - log2(1.5)
    thresh_bh <- max(differential.activity$p_val[differential.activity$p_val_adj <= 0.05], na.rm = TRUE)
    
    plot_list[[i]] <- volcano(res_tbl = differential.activity,
                              thresh_plot_1 = 0.05,
                              thresh_plot_2 = thresh_bh,
                              log2fc = 0.25,
                              top_x_to_plot = 25,
                              title = paste0(i, " TF in PSP vs. Ctrl"),
                              padj_col = "p_val",
                              log2FoldChange_col = "log2FoldChange",
                              ext_y = 1,
                              text_size = 3,
                              type = "gene",
                              cols = c("#B53737", "#4774a8", "grey"))
    plot_list[[i]]
    ggsave(paste0("volc_",i,".pdf"), width = 9, height = 5)
    
    write.csv(differential.activity, file = paste0(i, "_ct_tfme.csv"))
  })
}

# Combine into panel
combined_plot <- cowplot::plot_grid(plotlist = plot_list, ncol = 2, labels = "AUTO")  # 4 rows x 2 cols for 8 plots
ggsave("volcano_panel_all_celltypes_TFME.pdf", combined_plot, width = 17, height = 18)

# bidirect_barplot
df <- lapply(list.files(pattern = "*_ct_tfme.csv"), read.csv)
names(df) <-  str_remove(list.files(pattern = "*_ct_tfme.csv"),"_ct_tfme.csv")
df <- imap(df, ~ mutate(.x, Category = .y)) %>% bind_rows(.)
write.csv(df, 'ct_combined_data_tfme.csv', row.names = FALSE)

bidirect_plot(df = df, ct = "celltype", multiple_p_vals = F, log2FoldChange_col = "log2FoldChange", padj_tresh = 0.05) + 
  scale_fill_manual(values = ct_cols)
ggsave("bidirect_barplot_celltype_tfme.pdf", width = 7, height = 3.5)

df %>% filter(p_val_adj < 0.1) %>% nrow() #1686
df %>% filter(p_val_adj < 0.1 & Category == "Inh_Neu") %>% nrow() #400
df %>% filter(p_val_adj < 0.1 & Category == "Exc_ULN") %>% nrow() #275
df %>% filter(p_val_adj < 0.1 & Category == "Oligo") %>% nrow() #261
df %>% filter(p_val_adj < 0.1 & Category == "Exc_DLN") %>% nrow() #250
df %>% filter(p_val_adj < 0.1 & Category == "Astro") %>% nrow() #181
df %>% filter(p_val_adj < 0.1 & Category == "OPC") %>% nrow() #74
df %>% filter(p_val_adj < 0.1 & Category == "Micro-PVM") %>% nrow() #0

# DESeq2
deseq_ls <- prep_deseq2(combrna = combatac, 
                        feature_sample = "Pub_ID",
                        feature_group = "NPdiagnosis",
                        feature_cell_identity = "celltype",
                        cores = 32)
p1 <- DESeq2::plotPCA(deseq_ls$rld, ntop = 1000, intgroup = "celltype") + 
  scale_color_manual(values = ct_cols) + coord_equal(1.2) + theme_bw() 
p2 <- DESeq2::plotPCA(deseq_ls$rld, ntop = 1000, intgroup = "NPdiagnosis") + 
  scale_color_manual(values = group_cols) + coord_equal(1.2) + theme_bw() 
plot_grid(plotlist = list(p1,p2), align = "hv")
ggsave("pca_celltype.pdf", width = 9, height = 5)

rld_mat <- assay(deseq_ls$rld)
rld_cor <- cor(rld_mat)

# Plot heatmapdff
# Select cell type of interest
ann_col <- data.frame(row.names = deseq_ls$meta$clus_np_id, 
                      NPdiagnosis = deseq_ls$meta$NPdiagnosis,
                      Celltype = deseq_ls$meta$celltype)
pdf("cor_heatmap_celltype_tfme.pdf", height = 6, width = 8)
pheatmap(rld_cor, annotation_col = ann_col, 
         cluster_cols = T, 
         clustering_method = "ward.D2",
         annotation_colors = list(NPdiagnosis = group_cols,
                                  Celltype = ct_cols), 
         show_colnames = F, show_rownames = F,
         colorRampPalette(c("#4774a8", "white", "#B53737"))(50))
dev.off()

sessionInfo()
