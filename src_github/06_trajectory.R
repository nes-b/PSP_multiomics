# ---
# title: "snRNA+ATAC-seq in PSP frontal cortex"
# author: "Nils Briel, Center for Neuropathology, LMU Munich & Dept. Neurology, USZ Zurich"
# date: "Created: 2024-07-01"
# ---

# 1 Setup Environment ####
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
library(lme4)
library(viridis)
library(psupertime) #needs pandoc, BiocStyle, topGO to be installed
library(SeuratWrappers)
library(conflicted)
library(hdWGCNA)
library(flexplot)
library(PrettyCols)

source('../src/utils.R')
rpath <- "../results/psp/"

blas_set_num_threads(24)
future::plan("multicore", workers = 16)
options(future.globals.maxSize= 64000*1024^2) #24 000 MB
rlimit_as(150e9)

combrna <- readRDS("../results/psp/combrna_psp_scl2.Rds")

# 9 Trajectory Inference ####
if (!dir.exists("../results/psp/psupertime")) { dir.create("../results/psp/psupertime") }
setwd("../results/psp/psupertime")
set.seed(1234)

## 9.1 Oligo ####
cell_seu <- readRDS("../hdWGCNA/Oligosubset_hdwgcna_seuobj.Rds")
incols <- c("0"="#10abe8","1"="#a4e1f5","2"="#f5b295","3"="#e33b32")

### 9.1.1 Preprocess ####
# Set CB load to 0 (Ctrl) --> 4
cell_seu$CB_load <- ifelse(is.na(cell_seu$CB_load), 0, cell_seu$CB_load) %>% as.numeric()
Idents(cell_seu) <- cell_seu$Pub_ID

# cell_seu <- RunUMAP(cell_seu, dims = 1:30)
cell_seu <- IntegrateLayers(object = cell_seu, method = RPCAIntegration, normalization.method = "SCT", verbose = F)
conflicts_prefer(matrixStats::colSums2)
cell_seu <- FindNeighbors(cell_seu, reduction = "integrated.dr", dims = 1:30)
cell_seu <- FindClusters(cell_seu, resolution = 0.8)
cell_seu <- RunUMAP(cell_seu, dims = 1:30, reduction = "integrated.dr")
p1 <- do_DimPlot(cell_seu, reduction = "umap", group.by = c("subcluster"), raster = F, pt.size = 0.5, colors.use = sclus_cols) + 
  theme_void() + coord_fixed(1/1.2)
p2 <- do_DimPlot(cell_seu, reduction = "umap", group.by = c("CB_load"), raster = F, pt.size = 0.5, colors.use = incols) +
  theme_void() + coord_fixed(1/1.2)
p1|p2

# Area plot
data <- cell_seu@meta.data %>%
  group_by(CB_load, subcluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(CB_load) %>%
  mutate(percentage = n / sum(n),
         CB_load = factor(CB_load, levels = 0:3))
ggplot(data, aes(x = CB_load, y = percentage, fill = subcluster, group = subcluster)) +
  geom_area(alpha=0.75 , size=.5, colour="white") +
  scale_fill_brewer(palette = "RdBu",type = "qual", direction = -1) +
  labs(x = "CB Load", y = "Percentage") +
  theme_bw()
ggsave("proportion_npdiagnosis_subcluster_v2.pdf",width = 6, height=2)

sce <- Seurat::as.SingleCellExperiment(cell_seu, assay = "SCT")

# Filter the SCE object to remove unannotated transcripts
unannotated_transcripts <- c(grep(c("^AL\\d+\\.\\d+"), rownames(sce), value = TRUE), grep(c("^AC\\d+\\.\\d+"), rownames(sce), value = TRUE))
sce <- sce[!rownames(sce) %in% unannotated_transcripts, ]
sce <- sce[,sce$predicted.subclass == "Oligo"]

### 9.1.2 Psupertime ####
# prevent RAM overload
unix::rlimit_as(cur=100e9)
scepsu <- psupertime(sce, 
                     y = as.factor(sce$CB_load), 
                     sel_genes = "all",
                     assay_type = "logcounts",
                     scale = T,
                     smooth = T,
                     min_expression = 0.1)

### 9.1.3 Eval ####
plot_train_results(scepsu)
ggsave("olig_train_results.pdf", width = 7, height=3)
plot_labels_over_psupertime(scepsu, label_name='CB_load')
ggsave("olig_psup_density.pdf", width = 7, height=3)
psupertime::plot_identified_gene_coefficients(scepsu, n = 25)
ggsave("olig_psup_coeff_25.pdf", width = 5, height=4)

### make nebulosa plots!!!
plot_identified_genes_over_psupertime_modif(psuper_obj = scepsu, 
                                            label_name='CB_load', 
                                            size = 0.3,
                                            alpha = 0.5,
                                            se = T,
                                            sel = 2000,
                                            n_to_plot = 24, 
                                            plot_ratio = c(1.25))
ggsave("olig_psup_topgenes24.pdf", width = 11, height=5.2)
plot_predictions_against_classes(scepsu)
ggsave("olig_psup_predictions_against_classes.pdf", width = 6, height=5.2)
saveRDS(scepsu, "olig_sce.Rds")

### 9.1.4 Clustering and GO ####
unix::rlimit_as(cur=90e9)
go_list <- psup_go(scepsu, org_mapping = "org.Hs.eg.db")
saveRDS(go_list, "olig_go_list.Rds")

# heatmap of gene_clusters (go_list)
plot_dt = go_list$plot_dt[sample(2e6),]
ggplot(plot_dt) + 
  aes(x = cell_id, y = symbol, fill = value) +
  geom_tile() + 
  scale_fill_distiller(palette = "RdBu",limits = c(-3, 3)) +
  facet_grid(clust_label ~ ., scale = "free_y", space = "free_y") +
  theme_bw() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.text.y = element_text(angle = 0)) + 
  labs(x = "Cell", y = "Symbol",
       fill = "z-scored gene\nexpression")
ggsave("olig_heatmap_of_gene_cluster.png", width = 7, height=3.5)

plot_profiles_of_gene_clusters(go_list, label_name='CB Load', palette='RdBu')
ggsave("olig_profiles_of_gene_clusters.pdf", width = 5, height=6)


# plot specified genes over psupertime 
conflicted::conflicts_prefer(base::intersect)
toMatch <- unique (grep(paste(c("CDK5","TPKI","ERK","MEK","AMPK","DYRK","PKA",
                                "CAMK2","JNK","MAPK12","MAPK8","MAPT"), collapse="|"), 
                        rownames(sce), value=TRUE))
myelination <- c(
  "PLP1",  "MBP",  "MOBP",  "MAG",  "MOG",  "CD9",
  "CRYAB",   "CNP",  "NRG1",  "APP"
)

lipid <- c(
  "MID1IP1", "PLPP4", "AGMO", "NCSTN", "HACD2", 
  "ACADVL", "HMGCL", "LDLRAP1", "ELOVL1", "ABCA1", 
  "DHCR24", "HSD17B1", "FDFT1", "LBR", 
  "SQLE", "MSMO1", "DHCR7", "CYP51A1", "LSS"
)


plot_gene_psuper(scepsu, extra_genes = c(myelination,lipid), sel = 2000, plot_ratio = 2)
ggsave("olig_psup_sel_genes_myel_lipid_olig.pdf", width = 11, height=3.5)


### 9.1.5 Metadata Correlations ####
cell_seu$pseudotime <- scepsu$proj_dt$psuper
cell_seu$psuper_label <- scepsu$proj_dt$label_psuper
cell_seu$UMAP1 <- cell_seu@reductions$umap@cell.embeddings[,1]
cell_seu$UMAP2 <- cell_seu@reductions$umap@cell.embeddings[,2]
cell_seu$PC1 <- cell_seu@reductions$pca@cell.embeddings[,1]
cell_seu$PC2 <- cell_seu@reductions$pca@cell.embeddings[,2]
cell_seu$PC3 <- cell_seu@reductions$pca@cell.embeddings[,3]
cell_seu$PC4 <- cell_seu@reductions$pca@cell.embeddings[,4]

# 4-Plot
p2 <- FeaturePlot(cell_seu,features  = "CB_load", reduction = "integrated_dr", pt.size = 0.5, order = T) + 
  theme_test() +
  PrettyCols::scale_color_pretty_c(palette = "RedBlues", direction = -1)
p3 <- FeaturePlot(cell_seu,features  = "CB_load", reduction = "pca",dims = c(1,3), pt.size = 0.5, order = T) + 
  theme_test() +
  PrettyCols::scale_color_pretty_c(palette = "RedBlues", direction = -1)
p4 <- cell_seu@meta.data %>%
  ggplot(aes(y=subcluster, x=pseudotime)) +
  geom_jitter(aes(fill=pseudotime),shape = 21, size = 0.75, stroke = 0.1, color="grey89") +
  geom_boxplot(alpha = 0.75) +
  PrettyCols::scale_fill_pretty_c(palette = "RedBlues", direction = -1) +
  theme_test()
cowplot::plot_grid(plotlist =  list(p3,p2,p4), ncol = 1, align = "hv",
                   axis = "tblr") 
ggsave("olig_pseudotime_overview.pdf", width = 7, height=11)


#### 9.1.5.1 Pseudotime & CB load ####
library(ggpubr)
cell_seu@meta.data %>%
  group_by(Pub_ID, celltype) %>%
  summarise(
    CB_load = mean(CB_load),
    NFT_quant = mean(No_NFT_3.2e7mcsq),
    TA_quant = mean(No_TA_3.2e7mcsq),
    pseudotime = median(pseudotime)
  ) %>%
  ggscatter(
    x = "pseudotime",
    add = "reg.line",
    y = "CB_load",
    conf.int = T
  ) + geom_point(aes(x=pseudotime,y=CB_load, color = Pub_ID)) +
  stat_cor(method = "spearman") +
  theme_bw()
ggsave("pseudotime_cb_load_corr_opcoligo_overview.pdf", width = 5, height=4)


#### 9.1.5.1 Module Trajectory ####
p  <- plot_module_trajectory(
  cell_seu,
  ncol = 3,
  n_bins = 100,
  pseudotime_col = "pseudotime",
  harmonized = T,
  point_size = 0,
  line_size = 0.5,
  se = T,
) 
p
ggsave("pseudotime_me_opcoligo.pdf", width = 3, height=4,dpi = 500)

##
plot_gene_trajectory(
  seurat_obj = cell_seu,
  genes = c(scepsu$beta_dt$symbol[1:10]),
  pseudotime_col = "pseudotime",
  n_bins = 3,
  ncol = 5,
  point_size = 1,
  line_size = 1,
  se = TRUE,
  group_colors = T,
  patch = F
) + scale_fill_brewer(palette = "Set3") +
  scale_color_brewer(palette = "Set3")
ggsave("pseudotime_genes3bins_oligo.pdf", width = 4, height=4,dpi = 500)
plot_gene_trajectory(
  seurat_obj = cell_seu,
  genes = c(scepsu$beta_dt$symbol[1:10]),
  pseudotime_col = "pseudotime",
  n_bins = 25,
  ncol = 5,
  point_size = 1,
  line_size = 1,
  se = TRUE,
  group_colors = T,
  patch = F
)+ scale_fill_brewer(palette = "Set3") +
  scale_color_brewer(palette = "Set3")
ggsave("pseudotime_genes25bins_oligo.pdf", width = 4, height=4,dpi = 500)

mycols <- colorRampPalette(brewer.pal(8, "Set1"))(length(lipid))
plot_gene_trajectory(
  seurat_obj = cell_seu,
  genes = lipid,
  pseudotime_col = "pseudotime",
  n_bins = 3,
  ncol = 5,
  point_size = 1,
  line_size = 1,
  se = TRUE,
  group_colors = T,
  patch = F
)+   scale_fill_manual(values = mycols) +  
  scale_color_manual(values = mycols)
ggsave("pseudotime_lipid3bins_oligo.pdf", width = 4, height=4,dpi = 500)


mycols <- colorRampPalette(brewer.pal(8, "Set1"))(length(myelination))
plot_gene_trajectory(
  seurat_obj = cell_seu,
  genes = myelination,
  pseudotime_col = "pseudotime",
  n_bins = 3,
  ncol = 3,
  point_size = 1,
  line_size = 1,
  se = TRUE,
  group_colors = T,
  patch = F
)+   scale_fill_manual(values = mycols) +  
  scale_color_manual(values = mycols)
ggsave("pseudotime_myelin3bins_oligo.pdf", width = 4, height=4,dpi = 500)


#### 9.1.5.3 Explore Pseudoprogression - Oligo - Myelination ####


single <- cell_seu@meta.data %>%
  ggplot(aes(y=1, x=pseudotime, fill = NPdiagnosis)) +
  geom_jitter(aes(fill=NPdiagnosis), shape =21, alpha = 1, size = 1, color="grey20") +
  theme_test()  +  scale_fill_manual(values = group_cols)

comb <- cell_seu@meta.data %>%
  mutate(subcluster = fct_reorder(subcluster, pseudotime, .fun = median, .desc = T)) %>%
  ggplot(aes(y=subcluster, x=pseudotime, color = NPdiagnosis)) +
  geom_jitter(aes(fill=pseudotime),shape = 21, size = 1, stroke = 0.1, color="grey20") +
  geom_boxplot(alpha = 0.75) +
  PrettyCols::scale_fill_pretty_c(palette = "RedBlues", direction = -1) +
  scale_color_manual(values = group_cols) +
  theme_test()

props <- cell_seu@meta.data %>%
  mutate(pseudotime_bin = cut(pseudotime, breaks = 10)) %>%  # Bin pseudotime into 50 intervals
  group_by(pseudotime_bin, NPdiagnosis, subcluster) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(pseudotime_bin) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup() %>%
  mutate(bin_midpoint = as.numeric(sub("\\((.+),(.+)\\]", "\\1", pseudotime_bin)) + 
           (as.numeric(sub("\\((.+),(.+)\\]", "\\2", pseudotime_bin)) - 
              as.numeric(sub("\\((.+),(.+)\\]", "\\1", pseudotime_bin))) / 2) %>%
  ggplot(aes(x = bin_midpoint, y = proportion, fill = subcluster)) +
  geom_col(color = "grey20", size = 0.4) +
  theme_test() + 
  facet_grid(NPdiagnosis~.) +
  labs(
    x = "Pseudotime",
    y = "Proportion",
    color = "NP Diagnosis",
  ) +  scale_fill_manual(values = sclus_cols)

cowplot::plot_grid(
  plotlist =  list(single, comb, props),
  ncol = 1,
  align = "hv",
  rel_heights = c(0.6, 1, 1),
  axis = "tblr"
) 

ggsave("oligo_x_pseudotime_multipanel.pdf", width = 7, height = 8, dpi = 300)
##



SCpubr::do_BoxPlot(cell_seu[which(cell_seu@assays$RNA["ACADVL",]>0.1),],
                   feature = "ACADVL",
                   use_test = T,
                   use_silhouette  = T,
                   # min.cutoff = 0.1,
                   group.by = "subcluster",
                   comparisons = list(c("Oligo.1", "Oligo.2"),
                                      c("Oligo.1", "Oligo.3"),
                                      c("Oligo.1", "Oligo.4"))
) + theme_test()

SCpubr::do_GeyserPlot(sample = cell_seu,
                      features = "ACADVL",
                      slot = "data",
                      pt.size = 0.7,
                      min.cutoff = 0.1,
                      scale_type = "continuous",
                      group.by = "subcluster", split.by = "NPdiagnosis")

cowplot::plot_grid(plotlist =  list(p3,p2,p4), ncol = 1, align = "hv",
                   axis = "tblr") 

## 9.2 Inh Neu ####
# --> archived 20241103

## 9.3 Exc_ULN ####
# --> archived 20251019

## 9.4 Exc_DLN ####
# --> archived 20251019

## 9.5 All Neurons ####
cell_seu <- readRDS("../hdWGCNA/neursubset_hdwgcna_seuobj.Rds")
incols <- c("#033800","#99ce64","#9057c6","#5d2f89")

# Set NFT load to 0 (Ctrl) --> max(quant NFT)
cell_seu$NFT_load <- ifelse(is.na(cell_seu$No_NFT_3.2e7mcsq), 0, cell_seu$No_NFT_3.2e7mcsq)
# put labels into 0-4 scales
cell_seu@meta.data <- cell_seu@meta.data %>% mutate(
  NFT_load_dsct =  round(scales::rescale(.$NFT_load, to = c(0,4)),digits = 0)) %>%
  mutate(NFT_load_dsct = ifelse(NFT_load_dsct >= 3, NFT_load_dsct-1,NFT_load_dsct)) %>%
  mutate(NFT_load_dsct = as.factor(NFT_load_dsct))

Idents(cell_seu) <- cell_seu$Pub_ID

# Area plot
data <- cell_seu@meta.data %>%
  group_by(NFT_load_dsct, predicted.subclass) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(NFT_load_dsct) %>%
  mutate(percentage = n / sum(n))
ggplot(data, aes(x = NFT_load_dsct, y = percentage, fill = predicted.subclass, group = predicted.subclass)) +
  geom_area(alpha=0.75 , size=.5, colour="white") +
  scale_fill_manual(values = sclass_cols_orig) +
  labs(x = "NFT Load", y = "Percentage") +
  theme_bw()
ggsave("proportion_neur_nft_subcluster.pdf",width = 6, height=2.5)

### 9.5.2 Psupertime by Celltype ####
celltypes <- c("Exc_DLN", "Exc_ULN", "Inh_Neu")
scepsu_list <- list()
pseudotime_list <- list()
psuper_label_list <- list()

# prevent RAM overload
unix::rlimit_as(cur=1e12)

# Loop through each celltype
for (cell_type in celltypes) {
  
  cat("\n=== Processing celltype:", cell_type, "===\n")
  
  # Subset by celltype
  cell_seu_subset <- subset(cell_seu, celltype == cell_type)
  
  # Convert to SCE
  sce <- Seurat::as.SingleCellExperiment(cell_seu_subset, assay = "SCT")
  
  # Filter unannotated transcripts
  unannotated_transcripts <- c(
    grep("^AL\\d+\\.\\d+", rownames(sce), value = TRUE), 
    grep("^AC\\d+\\.\\d+", rownames(sce), value = TRUE)
  )
  sce <- sce[!rownames(sce) %in% unannotated_transcripts, ]
  
  # Run psupertime with celltype-specific variable features
  scepsu <- psupertime(sce, 
                       y = sce$NFT_load_dsct, 
                       sel_genes = "list",
                       gene_list = VariableFeatures(cell_seu_subset),
                       assay_type = "logcounts",
                       scale = T,
                       smooth = T,
                       min_expression = 0.1)
  
  # Store results
  scepsu_list[[cell_type]] <- scepsu
  pseudotime_list[[cell_type]] <- scepsu$proj_dt$psuper
  psuper_label_list[[cell_type]] <- scepsu$proj_dt$label_psuper
  
  # Save individual SCE object
  saveRDS(scepsu, paste0("neur_sce_", cell_type, ".Rds"))
  
  cat("Completed:", cell_type, "\n")
}

# Combine pseudotime and psuper_label back into cell_seu
cell_seu$pseudotime <- NA
cell_seu$psuper_label <- NA

for (cell_type in celltypes) {
  # Get cell IDs for this celltype
  cell_ids <- colnames(cell_seu)[cell_seu$celltype == cell_type]
  
  # Assign pseudotime and label
  cell_seu$pseudotime[cell_seu$celltype == cell_type] <- pseudotime_list[[cell_type]]
  cell_seu$psuper_label[cell_seu$celltype == cell_type] <- psuper_label_list[[cell_type]]
}

# Create celltype-specific subsets with pseudotime
cell_seu_exd <- cell_seu[,cell_seu$celltype == "Exc_DLN"]
cell_seu_exu <- cell_seu[,cell_seu$celltype == "Exc_ULN"]
cell_seu_inh <- cell_seu[,cell_seu$celltype == "Inh_Neu"]



### 9.5.3 Eval ####
# Generate evaluation plots for each celltype
for (cell_type in celltypes) {
  scepsu <- scepsu_list[[cell_type]]
  
  plot_train_results(scepsu)
  ggsave(paste0("neur_train_results_", cell_type, ".pdf"), width = 7, height=3)
  
  plot_labels_over_psupertime(scepsu, label_name='NFT_load') + 
    scale_fill_manual(values = incols) + 
    scale_color_manual(values = incols)
  ggsave(paste0("neur_psup_density_", cell_type, ".pdf"), width = 7, height=3)
  
  plot_identified_gene_coefficients(scepsu, n = 50)
  ggsave(paste0("neur_psup_coeff_", cell_type, ".pdf"), width = 9, height=4)
  
  plot_identified_genes_over_psupertime_modif(psuper_obj = scepsu, 
                                              label_name='NFT_load', 
                                              size = 0.3,
                                              alpha = 0.5,
                                              se = T,
                                              sel = 2000,
                                              n_to_plot = 40, 
                                              plot_ratio = c(1.25)) + 
    scale_fill_manual(values = incols)
  ggsave(paste0("neur_psup_topgenes_", cell_type, ".pdf"), width = 11, height=6)
  
  plot_predictions_against_classes(scepsu)
  ggsave(paste0("neur_predictions_", cell_type, ".pdf"), width = 7, height=5)
}


### 9.5.4 Module Trajectory ####
bins = 100
p1  <- plot_module_trajectory(
  cell_seu_exu,
  ncol = 3,
  n_bins = bins,
  pseudotime_col = "pseudotime",
  harmonized = T,
  point_size = 1,
  line_size = 1,
  se = T,
) + geom_smooth() + geom_vline(xintercept = c(50, 70), linetype = 2, alpha = 0.5) + 
scale_y_continuous(limits = c(-6,5)) 
p2  <- plot_module_trajectory(
  cell_seu_exd,
  ncol = 3,
  n_bins = bins,
  pseudotime_col = "pseudotime",
  harmonized = T,
  point_size = 1,
  line_size = 1,
  se = T,
) + geom_smooth() + geom_vline(xintercept = c(50, 70), linetype = 2, alpha = 0.5) + 
 scale_y_continuous(limits = c(-5,5))
p3  <- plot_module_trajectory(
  cell_seu_inh,
  ncol = 3,
  n_bins = bins,
  pseudotime_col = "pseudotime",
  harmonized = T,
  point_size = 1,
  line_size = 1,
  se = T,
) + geom_smooth() + geom_vline(xintercept = c(50, 70), linetype = 2, alpha = 0.5) + 
scale_y_continuous(limits = c(-5,5))

ggpubr::ggarrange(p1, p2, p3, ncol=3, nrow=1, common.legend = TRUE, legend="bottom")

ggsave("pseudotime_me_neur_integrated.pdf", width = 10, height=3.5,dpi = 500)


sessionInfo()