# ---
# title: "snRNA+ATAC-seq in PSP frontal cortex"
# author: "Nils Briel, Center for Neuropathology, LMU Munich & Dept. Neurology, USZ Zurich"
# date: "Created: 2024-07-01"
# ---

# 8 hdWGCNA ####
# First, we need to load the necessary packages:
set.seed(12345)
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
library(hdWGCNA)
library(umap)
library(igraph)
library(ggrepel)
source('../src/utils.R')
if (!dir.exists("../results/psp/hdWGCNA")) { dir.create("../results/psp/hdWGCNA") }
setwd("../results/psp/hdWGCNA")

# set envir
RhpcBLASctl::blas_set_num_threads(1)
rlimit_as(150e9)
theme_set(theme_test())
options(future.globals.maxSize= 64000*1024^2) #64 000 MB
WGCNA::enableWGCNAThreads(nThreads = 16)

# Load and set Seurat
combrna <- readRDS("../combrna_psp_scl2.Rds")
DefaultAssay(combrna) <- "SCT"
# filter out non-coding / not annotated genes:
unannotated_transcripts <- c(grep(c("^AL\\d+\\.\\d+"), rownames(combrna), value = TRUE),
                             grep(c("^AC\\d+\\.\\d+"), rownames(combrna), value = TRUE),
                             grep(c("^AD\\d+\\.\\d+"), rownames(combrna), value = TRUE),
                             grep(c("^AF\\d+\\.\\d+"), rownames(combrna), value = TRUE),
                             grep(c("^AJ\\d+\\.\\d+"), rownames(combrna), value = TRUE),
                             grep(c("^AP\\d+\\.\\d+"), rownames(combrna), value = TRUE))
combrna <- combrna[!rownames(combrna) %in% unannotated_transcripts, ]
combrna$celltype_neur <- case_when(combrna$celltype %in% c("Exc_DLN","Exc_ULN","Inh_Neu") ~ "neur",
                                   .default = as.character(combrna$celltype)) %>% as.factor()
seurat_obj <- SetupForWGCNA(combrna, gene_select = "fraction", fraction = 0.05, wgcna_name = "psp_hdwgcna")
seurat_obj_global_ds <- MetacellsByGroups(seurat_obj, group.by = c("celltype_neur", 'Pub_ID'),  
                                          reduction = 'umap', k = 25, max_shared = 10, 
                                          ident.group = 'celltype_neur')
seurat_obj_global_ds <- suppressMessages(suppressWarnings(NormalizeMetacells(seurat_obj_global_ds)))
saveRDS(seurat_obj_global_ds, "seurat_obj_global_ds_v3_celltype_neur.Rds")

# Loop through the cell types
conflicted::conflicts_prefer(dplyr::select)

for(cell_type in unique(seurat_obj_global_ds$celltype_neur)){
  message(paste0("hdWGCNA for ", cell_type,"\n"))
  
  if(cell_type == "neur"){
    seurat_obj <- custom_SetDatExpr_v5(seurat_obj_global_ds, group_name = cell_type, group.by='celltype_neur', assay = 'SCT', layer = 'data')
    seurat_obj$celltype <- seurat_obj$celltype_neur
  }else {
    seurat_obj <- custom_SetDatExpr_v5(seurat_obj_global_ds, group_name = cell_type, group.by='celltype', assay = 'SCT', layer = 'data')
  }
  
  seurat_obj <- TestSoftPowers(seurat_obj,powers = c(seq(1, 10, by = 1)), networkType = "signed", corFnc = "bicor")
  power <- GetPowerTable(seurat_obj) %>% dplyr::filter(SFT.R.sq > 0.8) %>% top_n(1) %>% .[1,1]
  hdWGCNA::PlotSoftPowers(seurat_obj)
  seurat_obj$subcluster_f <- factor(seurat_obj$subcluster, levels = names(sclus_cols))
  
  message("ConstructNetwork\n")
  conflicted::conflicts_prefer(WGCNA::cor)
  seurat_obj <- ConstructNetwork(seurat_obj,
                                 soft_power = power,
                                 overwrite_tom = T,
                                 minModuleSize = 50,
                                 tom_name = paste0(cell_type, collapse = "_"))
  
  pdf(paste0(cell_type, "dendro.pdf", collapse = "_"), height = 4, width = 7)
  PlotDendrogram(seurat_obj, main=paste(cell_type, "hdWGCNA Dendrogram"))
  dev.off()
  
  seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))
  seurat_obj <- ModEigengene(seurat_obj, exclude_grey = T,
                             vars.to.regress =  c('Age', "pmi"),
                             group.by.vars="Pub_ID")
  MEs <- GetMEs(seurat_obj, harmonized=T)
  seurat_obj <- ModuleConnectivity(seurat_obj, group.by = 'celltype',
                                   group_name = cell_type, harmonized = F)
  seurat_obj <- ResetModuleNames(seurat_obj,new_name = paste(paste0(cell_type,collapse = "_"), ".M", sep=""))
  hubs <- GetHubGenes(seurat_obj, n_hubs = 10)
  
  message("Plotting\n")
  library(conflicted)
  PlotKMEs(seurat_obj, n_hubs = 10, ncol=3, text_size = 2)
  ggsave(paste0(cell_type, "kme.pdf", collapse = "_"), width = 6, height=7)
  seurat_obj <- RunModuleUMAP(seurat_obj, n_hubs = 5)
  
  pdf(paste0(cell_type, "umap.pdf", collapse = "_"), width = 8, height = 7.5)
  ModuleUMAPPlot(seurat_obj,label_hubs=5)
  dev.off()
  
  # DMEs
  message("DMEs\n")
  mod_colors_df <- dplyr::select(GetModules(seurat_obj), c(module, color)) %>%
    dplyr::distinct(.) %>% dplyr::arrange(module)
  rownames(mod_colors_df) <- mod_colors_df$module
  group1 <- seurat_obj@meta.data %>% subset(celltype == cell_type & NPdiagnosis == "PSP") %>% rownames
  group2 <- seurat_obj@meta.data %>% subset(celltype == cell_type & NPdiagnosis == "Ctrl") %>% rownames
  DMEs <- FindDMEs(seurat_obj,
                   barcodes1 = group1,
                   barcodes2 = group2,
                   test.use='wilcox',
                   harmonized = F)
  
  PlotDMEsVolcano(seurat_obj, DMEs) + scale_fill_manual(values = mod_colors_df$color[-1]) + theme_test()
  ggsave(paste0(cell_type, "volcano.pdf", collapse = "_"), width = 6, height=5)
  PlotDMEsLollipop(seurat_obj, DMEs, wgcna_name="psp_hdwgcna", pvalue = "p_val_adj") + 
    scale_fill_manual(values = mod_colors_df$color[-1]) + theme_test()
  ggsave(paste0(cell_type, "lollipop.pdf", collapse = "_"), width = 4, height=5)
  
  # Find DMEs for subclusters
  cell_seu <- subset(seurat_obj, celltype == cell_type) 
  
  if(cell_type == "neur"){
    tmp <- left_join(cell_seu@meta.data, cell_identity_levels[,1:2], by = "predicted.subclass", relationship = "many-to-many")
    cell_seu$subcluster <- tmp$celltype.y
    clusters <- factor(unique(cell_seu$subcluster)) %>% 
      str_replace_all(., pattern = c("_"), "-") %>%
      as.factor() %>% sort()
  }else{
    clusters <- factor(unique(subset(cell_seu, celltype == cell_type)$subcluster)) %>% 
      str_replace_all(., pattern = c("_"), "-") %>%
      as.factor() %>% sort()
  }
  
  # Run DME loop for subclusters
  DMEs <- data.frame()
  
  if (length(clusters) > 1) {
      DMEs <- DME_loop(seurat_obj = cell_seu,
                       clusters = levels(cell_seu$subcluster))
    } else{
      DMEs <- FindDMEs(
        cell_seu,
        barcodes1 = subset(cell_seu@meta.data, NPdiagnosis == "PSP") %>% rownames,
        barcodes2 = subset(cell_seu@meta.data, NPdiagnosis == "Ctrl") %>% rownames,
        test.use = 'wilcox',
        pseudocount.use = 0.01
      )
      DMEs$cluster <- clusters
    }

  message("Compute GO")
  cell_seu <- compute_go_hdwgcna(seurat_obj = cell_seu, dbs = "ALL")
  
  df <- cell_seu@misc$psp_hdwgcna$enrichr_table %>%
    mutate_at(vars(ID, Description, geneID, module), as.factor) %>%
    mutate_at(vars(GeneRatio, BgRatio), as.factor) %>%
    mutate_at(vars(pvalue, p.adjust, qvalue, Count), as.numeric) %>% 
    group_by(module) %>%
    slice_min(order_by = p.adjust, n = 10) %>% # adjust top vars
    ungroup()  
  plot_go_facet(df, y = Description, x = module, direct = F) 
  ggsave(paste0(cell_type, "go_dotplot.pdf", collapse = "_"), width = nrow(df)/3.8+5, height=length(levels(df$module))/1.4+5)
  
  w_h <- data.frame(row.names = 1, module = c(levels(cell_seu$celltype), "OPCOligo", "pvalb", "neur"), 
                    #       Astro  Oligo  Exc_DLN  Exc_ULN  Endo_VLMC OPC Micro_PVM Inh_Neu OPCOligo pvalb neur 
                    width = c(11,  13,    12,      11,      12,       11,  11,       14,      12,     14, 11),
                    height = c(6,  9,     6,       6,       9,        6,   6,        6,        9,       5, 6))
  
  
  message("Saving.\n")
  saveRDS(cell_seu, paste0(cell_type, "subset_hdwgcna_seuobj.Rds", collapse = "_"))
  write_csv(hdWGCNA::GetModules(cell_seu), paste0(cell_type, "subset_hdwgcna_modules_genes.csv", collapse = "_"))
  
  message("Iterate over celltypes and construct complexHeatmaps\n")
  if(cell_type == "neur"){
    sccomp_res <- read_csv("../sccomp/sccomp_subcluster_structured_data.csv") %>% 
      dplyr::filter(celltype %in% c("Exc_DLN","Exc_ULN","Inh_Neu")) %>%
      mutate(celltype = "neur")
  }else {
    sccomp_res <- read_csv("../sccomp/sccomp_subcluster_structured_data.csv") %>% 
      dplyr::filter(celltype %in% cell_type)
  }
  
  
  if(length(clusters)>1){
    cell_seu$subcluster <- factor(cell_seu$subcluster)
    h2 <- plot_heatmap_mods(cell_seu, DMEs, GO = df, n_terms = 5, modus = "DME", cluster_columns = F,cluster_rows = F, 
                            sccomp_res = sccomp_res, 
                            color_pal = colorRamp2(c(-3, 0, 3), c("#4774a8", "white", "#B53737")),
                            w = w_h[cell_type,"width"]*1, h = w_h[cell_type,"height"]*1.8)
   pdf(paste0("compl_heat",cell_type,".pdf", collapse = "_"), width =  w_h[cell_type,"width"]*1.8,
        height =  w_h[cell_type,"height"]*1.8)
    draw(h2, row_title = "", row_title_gp = gpar(col = "red"), cluster_rows = T, 
         column_title = paste0(cell_type, collapse = " / "), column_title_gp = gpar(fontsize = 16, fontface= "bold"))
    dev.off()
  }else{
    pdf(paste0("compl_heat",cell_type,".pdf", collapse = "_"), width =  w_h[cell_type,"width"]*2.5,
        height =  w_h[cell_type,"height"]*1.8)
    draw(plot_heatmap_mods(cell_seu, DMEs, GO = df, n_terms = 5, modus = "DME",
                           sccomp_res = NULL,cluster_columns = F,
                           w = w_h[cell_type,"width"], h = w_h[cell_type,"height"]*2), row_title = "", row_title_gp = gpar(col = "red"), cluster_rows = T,
         column_title = paste0(cell_type), column_title_gp = gpar(fontsize = 16, fontface= "bold"))
    dev.off()
  }
  
  if(cell_type == "neur"){
    h2 <- plot_heatmap_mods(cell_seu, DMEs, GO = df, n_terms = 3, modus = "DME", cluster_columns = F, 
                            sccomp_res = sccomp_res,
                            color_pal = colorRamp2(c(-3, 0, 3), c("#8dc657", "white", "#9057c6")),
                            w = w_h[cell_type,"width"]*1, h = w_h[cell_type,"height"]*1.8)
    pdf(paste0("compl_heat_",cell_type,".pdf", collapse = "_"), width =  w_h[cell_type,"width"]*1.8,
        height =  w_h[cell_type,"height"]*2)
    draw(h2, row_title = "", row_title_gp = gpar(col = "red"), cluster_rows = T, 
         column_title = paste0(cell_type, collapse = " / "), column_title_gp = gpar(fontsize = 16, fontface= "bold"))
    dev.off()
  }
  
}


sessionInfo()
