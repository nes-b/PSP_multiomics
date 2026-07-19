# ---
# title: "snRNA-seq processing
# project: "snRNA+ATAC-seq in PSP frontal cortex"
# author: "Nils Briel, Center for Neuropathology, LMU Munich & Dept. Neurology, USZ Zurich"
# date: "Created: 2024-09-22"
# ---

## Setup Environment ####
# First, we need to load the necessary packages:
set.seed(1234)
library(Seurat)
options(Seurat.object.assay.version = "v5")
library(SeuratData)
devtools::install_github("satijalab/azimuth", "master")
library(Azimuth)
library(Signac)

library(scDblFinder)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
library(AUCell)
library(SummarizedExperiment)
library(GSEABase)

library(data.table)
library(Matrix)
library(parallel)
library(tidyverse)
library(RhpcBLASctl)
library(unix)
library(SCpubr)
library(parallel)
library(RColorBrewer)

blas_set_num_threads(1)
future::plan("multicore", workers = 16)
options(future.globals.maxSize= 64000*1024^2) #64 000 MB
rlimit_as(1e512)

source('../src/utils.R')
rpath <- "../results/psp/"
if(!dir.exists(rpath)){
  dir.create(rpath)
}

## snRNA-seq - Processing
# 1 QC & Metadata ####
meta <- readxl::read_excel("../data/case_metadata_coh1_2.xlsx", n_max = 24) %>% .[-5,-c(2,7,21,22,23,24,25)] %>% 
  dplyr::filter(NPdiagnosis %in% c("PSP","Ctrl")) %>%
  mutate(Sample = paste('num', .$Code_ID, sep=''))
sample_lib <- data.frame(Sample = meta$Sample, meta$Pub_ID) %>% 
  `rownames<-`(.$Sample) %>% 
  .[match(meta$Sample, .$Sample), ]

# Use lapply to read in the data for each sample and create a Seurat object
if(!file.exists("../results/psp/combrna_psp_scl2.Rds")){
  cat("File not found. Running snRNAseq processing pipeline.\n")
  seurat_list <- list()
  
## 1.1 + 1.2 QC 0 ####
  seurat_list <- mclapply(sample_lib$Sample, mc.cores = 8, function(i){
    print(i)
    data_dir <- paste0("../data/fastq_rna/", i, "/outs/filtered_feature_bc_matrix/")
    expression_matrix <- Read10X(data.dir = data_dir)
    seurat_object <- CreateSeuratObject(counts = expression_matrix, min.cells = 3, min.features = 200, project = "psp_rna")
    seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
    seurat_object[["percent.ribo1"]] <- PercentageFeatureSet(seurat_object, pattern = c("^RPS-"))
    seurat_object[["percent.ribo2"]] <- PercentageFeatureSet(seurat_object, pattern = c("^RPL-"))
    seurat_object[["percent.heme"]] <- PercentageFeatureSet(seurat_object, pattern = c("^HB[^(P)]"))
    seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5 & percent.heme < 0.5)
  
    # 1.2 Doublet scoring ####
    set.seed(123)
    sce = SingleCellExperiment(
      list(counts=expression_matrix[,colnames(seurat_object)]
           ),colData = seurat_object@meta.data
    ) 
    sce <- scDblFinder(sce, clusters=T)
    
    seurat_object <- NormalizeData(seurat_object)
    seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 5000)
    seurat_object[["scDblFinder_score"]] <- sce$scDblFinder.score
    seurat_object[["scDblFinder_class"]] <- sce$scDblFinder.class
    return(seurat_object) 
    rm(sce, data_dir, expression_matrix)
  })
  names(seurat_list) <- sample_lib$Sample
  
  library(conflicted) ### call here and not earlier
  conflicted::conflicts_prefer(Seurat::Assays)

  ## 1.3 Normalization, Clustering & UMAP inspection  ####
  set.seed(42)
  combrna <- merge(
    x = seurat_list[[1]],
    y = seurat_list[2:length(sample_lib$Sample)],
    add.cell.ids = sample_lib$meta.Pub_ID
  )
  
  combrna$Pub_ID <- rownames(combrna@meta.data) %>% 
    as.data.frame() %>% 
    separate('.', c("a", "b"), sep = '_') %>% 
    .[,1] %>% as.factor
  
  conflicted::conflicts_prefer(dplyr::rename)
  merged_df <- as.data.frame(combrna@meta.data) %>% 
    left_join(., meta, by = c('Pub_ID'), keep=F) %>% 
    rename(cohort = "...1") %>%
    data.frame(.)
  
  combrna@meta.data <- merged_df %>%
    `rownames<-`(colnames(combrna@assays$RNA))

  ## 1.4 Layers | SCTransform ####
  set.seed(42)
  combrna <- JoinLayers(combrna)
  combrna[["RNA"]] <- split(combrna[["RNA"]], f = combrna$Pub_ID)
  Idents(combrna) <- combrna$Pub_ID
  combrna <- SCTransform(combrna, vars.to.regress = c('nCount_RNA', 'Age', "percent.mt"), verbose = T)
  combrna <- RunPCA(combrna, npcs = 50, verbose= F)
  combrna <- RunUMAP(combrna, dims = 1:30)
  combrna <- IntegrateLayers(object = combrna, method = RPCAIntegration, normalization.method = "SCT", verbose = F)
  conflicts_prefer(matrixStats::colSums2)
  combrna <- FindNeighbors(combrna, reduction = "integrated.dr", dims = 1:30)
  combrna <- FindClusters(combrna, resolution = 0.8)
  combrna <- RunUMAP(combrna, dims = 1:30, reduction = "integrated.dr")
  saveRDS(combrna, "../results/psp/combrna_psp_unfiltered.Rds")
  
  # clean data
  rm(seurat_list, obj, meta, merged_df)
  gc()
  combrna$pmi <- as.numeric(combrna$pmi) 
  
  ## 1.5 QC 1: technical & biological ####
  ### 1.5.1 Diagnostic Plots ####
  p1 <- do_DimPlot(combrna, reduction = "umap", pt.size = 0.1,shuffle = T, group.by = 'NPdiagnosis',
             label = F, label.box = F, repel=F, colors.use = group_cols,  plot.title = "NPdiagnosis | QC 1") +
    theme_test() + coord_fixed(1/1.2)
  p2 <- do_DimPlot(combrna, reduction = "umap", pt.size = 0.1, shuffle = T, group.by = 'Pub_ID',
                   label = F, label.box = F, repel=F, plot.title = "ID | QC 1") +
    theme_test() + coord_fixed(1/1.2)
  p3 <- do_FeaturePlot(combrna, reduction = "umap", pt.size = 0.1, features = 'Age', raster.dpi = 512, 
                       legend.length = 1, legend.width = 7, plot.title = "Age | QC 1") + 
    theme_test() + coord_fixed(1/1.2)
  p4 <- do_FeaturePlot(combrna, reduction = "umap", pt.size = 0.1, features = 'pmi', raster.dpi = 512, 
                       legend.length = 1, legend.width = 7, plot.title = "PMI | QC 1") + 
    theme_test() + coord_fixed(1/1.2)
  p5 <- do_FeaturePlot(combrna, reduction = "umap", pt.size = 0.1, features = 'nCount_RNA', raster.dpi = 512, 
                       legend.length = 1, legend.width = 7, plot.title = "nCount RNA | QC 1") + 
    theme_test() + coord_fixed(1/1.2)
  p6 <- do_FeaturePlot(combrna, reduction = "umap", pt.size = 0.1, features = 'percent.mt', raster.dpi = 512, 
                       legend.length = 1, legend.width = 7, plot.title = "% mitochondrial reads | QC 1") + 
    theme_test() + coord_fixed(1/1.2)

  cowplot::plot_grid(plotlist = list(p1,p2,p3,p4, p5, p6), ncol = 2, align = "h")
  ggsave("QC1_diagnostics_prefilter.pdf", path = rpath, width = 10, height=14)
  
  ### 1.5.2 Dbl filter ####
  p1 <- do_FeaturePlot(subset(combrna, subset = scDblFinder_score < 0.7), reduction = "umap", features = 'scDblFinder_score',
                       pt.size = 0.05, raster.dpi = 512, legend.length = 1, legend.width = 7, order = T,
                       plot.title = "scDblFinder_score | doublet probability score filter: 0.7") + 
    theme_test()+ coord_fixed(1/1.2) + scale_color_gradientn(colours = c("grey99","#fad9d9","#fac0c0","red4"),  limits = c(0, 1))
  p2 <- do_FeaturePlot(subset(combrna, subset = scDblFinder_score < 0.5), reduction = "umap",features = 'scDblFinder_score',  
                       pt.size = 0.05, raster.dpi = 512, legend.length = 1, legend.width = 7, order = T,
                       plot.title = "scDblFinder_score | doublet probability score filter: 0.5") + 
    theme_test()+ coord_fixed(1/1.2) + scale_color_gradientn(colours = c("grey99","#fad9d9","#fac0c0","red4"),  limits = c(0, 1))
  p3 <- do_FeaturePlot(subset(combrna, subset = scDblFinder_score < 0.3), reduction = "umap", features = 'scDblFinder_score', 
                       pt.size = 0.05, raster.dpi = 512, legend.length = 1, legend.width = 7, order = T,
                       plot.title = "scDblFinder_score | doublet probability score filter: 0.3") + 
    theme_test()+ coord_fixed(1/1.2) + scale_color_gradientn(colours = c("grey99","#fad9d9","#fac0c0","red4"),  limits = c(0, 1))
  p4 <- ggplot(combrna@meta.data, aes(x = scDblFinder_score)) + 
    geom_histogram(binwidth = 0.01, fill = "grey80", color = "grey4") + 
    theme_test()
  cowplot::plot_grid(plotlist = list(p1,p2,p3,p4), ncol = 2, align = "h")
  ggsave("QC1_doublet_filter_steps.pdf", path = rpath, width = 10, height=11)
  
  ### 1.5.2 QC 1: Doublet & Mito filtering ####
  combrna <- subset(combrna, subset = scDblFinder_score < 0.5 & percent.mt < 3)

# 2 Cell Identity ####
  
  ### 2.1 Azimuth ####
  combrna <- JoinLayers(combrna, assay = "RNA")
  combrna <- RunAzimuth(combrna, reference = "humancortexref", assay = "RNA", umap.name = "umap.ref")
  DefaultAssay(combrna) <- "SCT"
  
  # define less granular celltypes
  ct_lib <- data.frame(predicted.subclass = c("Astro","Endo", "L2/3 IT", "L5 ET", "L5 IT", "L5/6 NP", "L6 CT","L6 IT","L6 IT Car3",
                                              "L6b","Lamp5","Micro-PVM","Oligo","OPC","Pvalb", "Sncg","Sst","Sst Chodl","Vip","VLMC"),
                       celltype = as.factor(c("Astro","Endo_VLMC", "Exc_ULN", "Exc_DLN", "Exc_DLN", "Exc_DLN", "Exc_DLN","Exc_DLN","Exc_DLN",
                                              "Exc_DLN","Inh_Neu","Micro_PVM","Oligo","OPC","Inh_Neu", "Inh_Neu",
                                              "Inh_Neu","Inh_Neu","Inh_Neu","Endo_VLMC")))
  
  combrna@meta.data <- left_join(combrna@meta.data, ct_lib, by="predicted.subclass")
  rownames(combrna@meta.data) = colnames(combrna@assays$RNA$data)
  saveRDS(combrna, "../results/psp/combrna_psp_libfiltered.Rds")
  
  ### 2.2 Compute sub-clusters ####
  RhpcBLASctl::blas_set_num_threads(1)
  if (!dir.exists("../results/psp/subclusters")) { dir.create("../results/psp/subclusters") }
  setwd("../results/psp/subclusters")
  conflicted::conflicts_prefer(MatrixGenerics::colSums2)
  unique_cell_types <- unique(combrna$predicted.subclass) %>% .[!. %in% c("Endo","Sst Chodl",NA)] 
  
  # to get the nearest subcluster in low-abundant cell types, where we don't compute subclusters
  city_neighborhood <- FindNeighbors(combrna, return.neighbor=TRUE)
  DefaultAssay(combrna) <- "RNA"
  cat("Finished set up. Computing subclusters...")
  
  for(i in unique_cell_types){
    # Runs scSHC::scSHC() under the hood
    subclust_celltype(combrna = combrna, cell_type = i, 
                      city_neighborhood = city_neighborhood,
                      cores = 6)
  }
  cat("Finished computing subclusters. Saving and merging...")
  saveRDS(combrna, "../results/psp/combrna_psp_scl_v2.Rds")
  
  ct <- janitor::make_clean_names(unique_cell_types)
  cl_ls <- list()
  for(i in 1:length(ct)){
    k = ct[[i]]
    cls_rds <- paste0(k, "_subclust.Rds")
    if(file.exists(paths = cls_rds)) {
      cl_ls[[i]] <- readRDS(cls_rds)
    }else{
      cls_rds <-
        paste0(unique_cell_types[i], "_subclust.Rds")
      cl_ls[[i]] <- readRDS(cls_rds)
    }
  }
  names(cl_ls) <- ct 
  ct_comb <- merge(x = cl_ls[[1]],
                   y = cl_ls[2:length(ct)])

  ct_comb$subcluster_scSHC_ct <- paste(ct_comb$predicted.subclass, ct_comb$subcluster_scSHC, sep=".")
  
  combrna$subcluster <- "0"
  combrna$subcluster[colnames(combrna)] <- ct_comb@meta.data[colnames(combrna), "subcluster_scSHC_ct"]
  combrna$subcluster <- str_replace(combrna$subcluster, pattern = c(" "), "_") 
  combrna$subcluster <- str_replace(combrna$subcluster, pattern = c(" "), "_") 
  combrna$subcluster <- str_replace(combrna$subcluster, pattern = c("/"), "_") 
  combrna$subcluster <- str_replace(combrna$subcluster, pattern = c("-"), "_") 
  combrna$subcluster <- str_replace(combrna$subcluster, pattern = c("_special_NN"), ".x") 
  combrna$subcluster <- factor(combrna$subcluster, levels = str_sort(levels(as.factor(combrna$subcluster)),numeric = TRUE))
  
  # clean data
  rm(ct_comb)
  gc()
  
  ### 2.3 QC 2: Cell type prediction ####
  # compare pre- and post Azimuth prediction score cut offs
  combrna_unfiltered <- combrna 
  combrna <- subset(combrna, subset = predicted.subclass.score > 0.6)
  ncol(combrna_unfiltered) # 88 883
  ncol(combrna) # 80 904
  
  p1 <- do_FeaturePlot(combrna_unfiltered, reduction = "umap", pt.size = 0.25, features = 'predicted.class.score', raster.dpi = 256, 
                       legend.length = 1, legend.width = 7, plot.title = 'Azimuth: Predicted.class.score | unfiltered') + 
    theme_test() + coord_fixed() + scale_color_gradientn( colours = c('red3',"grey70","grey99"),  limits = c(0, 1))
  p2 <- do_FeaturePlot(combrna, reduction = "umap", pt.size = 0.25, features = 'predicted.class.score', raster.dpi = 256, 
                       legend.length = 1, legend.width = 7, plot.title = 'Azimuth: Predicted.class.score | QC2') + 
    theme_test() + coord_fixed() + scale_color_gradientn( colours = c('red3',"grey70","grey99"),  limits = c(0, 1))
  
  p3 <- do_FeaturePlot(combrna_unfiltered, reduction = "umap", pt.size = 0.25, features = 'predicted.subclass.score', raster.dpi = 256, 
                       legend.length = 1, legend.width = 7, plot.title = 'Azimuth: Predicted.subclass.score | unfiltered') + 
    theme_test() + coord_fixed() + scale_color_gradientn( colours = c('red3',"grey70","grey99"),  limits = c(0, 1))
  p4 <- do_FeaturePlot(combrna, reduction = "umap", pt.size = 0.25, features = 'predicted.subclass.score', raster.dpi = 256, 
                       legend.length = 1, legend.width = 7, plot.title = 'Azimuth: Predicted.subclass.score | QC2') + 
    theme_test() + coord_fixed() + scale_color_gradientn( colours = c('red3',"grey70","grey99"),  limits = c(0, 1))
  cowplot::plot_grid(plotlist = list(p1,p2,p3,p4), ncol = 2, align = "h")
  ggsave("../QC2_azimuth_predictions.pdf", width = 10, height=10)
  
  combrna <- RunUMAP(combrna, dims = 1:30, reduction = "integrated.dr", reduction.name = "umap_v1")
  
  # clean data
  rm(combrna_unfiltered)
  gc()
  
  DefaultAssay(combrna) <- "SCT"
  Idents(combrna) <- combrna$subcluster
  saveRDS(cl_ls,"combrna_subclusters_v2.Rds")
  setwd("../")
  
  # 3 NP Trait Data Correlation ####
  cat("Adding NP trait data.\n")
  # Add Tau load data
  t_load <- readxl::read_excel("../../data/tau_load.xlsx", sheet = "Sheet3") %>% 
    .[,c("Pub_ID","TA_load","NFT_load", "CB_load", "Threads_load")]
  combrna@meta.data <- left_join(combrna@meta.data, t_load, by = 'Pub_ID') %>%
    `rownames<-`(rownames(combrna@meta.data))
  t_load <- left_join(t_load, subset(meta, NPdiagnosis == "PSP")[,c("Pub_ID", "Age", "pmi", "Braak&Braak (NFT)", "Thal-Phase")], 
                      by = "Pub_ID") %>%
    column_to_rownames(., "Pub_ID") %>% data.matrix()
  
  res <- Hmisc::rcorr(t_load, type="spearman") 
  res$P[which(is.na(res$P))] <- 0
  pdf("corrplot_tau_load.pdf",width = 6, height=5)
  corrplot::corrplot(res$r, type = "full", order = "AOE", method = "color",
           tl.col = "black", tl.srt = 45, p.mat = res$P, diag = F,insig = "pch")
  dev.off()
  
  saveRDS(combrna, "combrna_psp_scl2.Rds")
  setwd("../../src")
  
}else{
  combrna <- readRDS("../results/psp/combrna_psp_scl2.Rds") ## classic data set version
}

# plotting
cellsel <- sample(colnames(combrna), 2e4)
marker_genes <- c(
  # General excitatory neuron markers
  "SLC17A7", "GRIN1", "SNAP25", "SYT1",
  # Layer-specific excitatory neuron markers 
  "RELN", # Layer 1 neurons
  "CUX2", "SATB2", # L2/3 IT
  "RORB", # L4 IT
  "BCL11B", "FEZF2", # L5 IT
  "TLE4", # L6 IT
  # General inhibitory neuron markers
  "GAD1", "GAD2", "CALB1", "CALB2",
  # Specific inhibitory neuron subclass markers
  "SST", # Sst
  "PVALB", # Pvalb
  "VIP", # Vip
  "SNCG", # Sncg
  # Non-neuronal cell markers
  "AQP4", "GFAP", "C3", # Astrocytes 
  "CD74", "CSF1R", # Microglia
  "MBP", "MOBP", "PLP1", # Oligodendrocytes
  "PDGFRA", # OPCs
  "FLT1", "CLDN5", # Endothelial cells
  "VCAN" # Vascular and leptomeningeal cells
)
# subclusters
combrna$subcluster <- factor(combrna$subcluster, levels = names(sclus_cols))
DotPlot(combrna, features = marker_genes, group.by = "subcluster") + scale_y_discrete(limits=rev) + 
  theme_bw() + coord_fixed(1.2/1.2) + RotatedAxis() + 
  colorspace::scale_color_continuous_diverging()
ggsave("dotplot_subclus_markers_v2.pdf", path = rpath, width = 9, height=10)

# subclasses
combrna$predicted.subclass_f <- as.factor(combrna$predicted.subclass)
combrna$predicted.subclass <- factor(combrna$predicted.subclass, levels = names(sclass_cols_orig))
DoHeatmap(combrna, cells=cellsel, features = marker_genes, size = 3, 
          group.by = "predicted.subclass_f", group.colors = ct_cols, slot = "scale.data") + 
  colorspace::scale_fill_continuous_diverging(palette = "Blue-Red 3")
ggsave("heatmap_pred.subclass_markers_v2.png", path = rpath, width = 7, height=3)

DotPlot(combrna, features = marker_genes, group.by = "predicted.subclass_f") + 
  scale_y_discrete(limits=rev(c("L2/3 IT", "L5 ET", "L5 IT", "L5/6 NP", "L6 CT",
                                       "L6 IT", "L6 IT Car3", "L6b", "Lamp5", "Pvalb",
                                       "Sncg", "Sst", "Sst Chodl", "Vip", "Astro", "Micro-PVM",
                                       "Oligo", "OPC", "Endo_VLMC"))) + 
  theme_bw() + coord_fixed(1.2/1.2) + RotatedAxis() + 
  colorspace::scale_color_continuous_diverging() 
ggsave("dotplot_pred.subclass_markers_v2.pdf", path = rpath, width = 8, height=6)

# celltypes
combrna$celltype <- factor(as.factor(combrna$celltype), levels = names(ct_cols))
DoHeatmap(combrna, cells=cellsel, features = marker_genes, size = 3, 
          group.by = "celltype", group.colors = ct_cols, slot = "scale.data") + 
  colorspace::scale_fill_continuous_diverging(palette = "Blue-Red 3")
ggsave("heatmap_celltype_markers_v2.png", path = rpath, width = 7, height=3)

combrna$celltype <- factor(as.factor(combrna$celltype), levels = rev(names(ct_cols)))
DotPlot(combrna, features = marker_genes, group.by = "celltype") + theme_bw() + 
  coord_fixed(1.2/1.2) + RotatedAxis() + 
  colorspace::scale_color_continuous_diverging()
ggsave("dotplot_celltype_markers.pdf", path = rpath, width = 8, height=5)
combrna$celltype <- factor(combrna$celltype, levels = rev(levels(combrna$celltype)))

# technical
FeaturePlot(combrna, reduction = "umap.ref", features = c('predicted.subclass.score'), raster = F) + theme_test() + coord_fixed(1/1.2)+ 
  colorspace::scale_color_continuous_sequential(palette = "OrRd", limits = c(0,1.01)) 
ggsave("predicted.subclass.score_v2.pdf", path = rpath, width = 6, height=5)


p1 <- do_DimPlot(combrna, reduction = "umap_v1", pt.size = 0.5,shuffle = T, label = F, label.box = F, repel=F,colors.use = sclus_cols, group.by = 'subcluster')+  theme_test() + NoLegend() + coord_fixed(1/1.2)
p2 <- do_DimPlot(combrna, reduction = "umap_v1", pt.size = 0.5,shuffle = T, label = T, label.box = F, repel=T,colors.use = ct_cols, group.by = 'celltype')+  theme_test() + NoLegend() + coord_fixed(1/1.2)
p3 <- do_DimPlot(combrna, reduction = "umap_v1", pt.size = 0.5,shuffle = T, label = T, label.box = F, repel=T,colors.use = sclass_cols_orig, group.by = 'predicted.subclass')+  theme_test() + NoLegend() + coord_fixed(1/1.2)
p4 <- do_DimPlot(combrna, reduction = "umap_v1", pt.size = 0.5,shuffle = T, group.by = 'NPdiagnosis') + theme_test() + coord_fixed(1/1.2) + scale_color_manual(values= group_cols)
cowplot::plot_grid(plotlist = list(p1,p2,p3,p4), ncol = 2, align = "h")
ggsave("dimplot_clust_v2.pdf", path = rpath, width = 12, height=12, dpi = 600)
Idents(combrna) <- "predicted.subclass"

# Trait umaps
p1 <- FeaturePlot(combrna, reduction = "umap_v1", features = "TA_load", raster = T) + theme_test() + coord_fixed(1/1.2)+ 
  colorspace::scale_color_continuous_sequential(palette = "OrRd")
p2 <- FeaturePlot(combrna, reduction = "umap_v1", features = "CB_load", raster = T) + theme_test() + coord_fixed(1/1.2)+ 
  colorspace::scale_color_continuous_sequential(palette = "PurpOr")
p3 <- FeaturePlot(combrna, reduction = "umap_v1", features = "NFT_load", raster = T) + theme_test() + coord_fixed(1/1.2)+ 
  colorspace::scale_color_continuous_sequential(palette = "BuGn")
p4 <- FeaturePlot(combrna, reduction = "umap_v1", features = "Threads_load", raster = T) + theme_test() + coord_fixed(1/1.2)+ 
  colorspace::scale_color_continuous_sequential(palette = "Greens")
cowplot::plot_grid(plotlist = list(p1,p2,p3,p4), ncol = 2)
ggsave("dimplot_tau_load_v2.pdf", path = rpath, width = 12, height=12, dpi = 600)

# filter low number clusters
combrna@meta.data %>%
  dplyr::select(Pub_ID, NPdiagnosis, predicted.subclass, celltype) %>%
  group_by(Pub_ID, predicted.subclass,NPdiagnosis) %>%
  summarise(n_cells = n()) %>% 
  ggpubr::ggboxplot(x = "predicted.subclass", y = "n_cells") + coord_flip() +
  theme_bw() + geom_hline(yintercept = 25)+ geom_point(aes(color = Pub_ID)) 
ggsave("cluster_ncells.pdf", path = rpath, width = 4, height=7)

combrna <- subset(combrna, subset = predicted.subclass.score > 0.4)
ncol(combrna) #77 930

combrna$subcluster <- factor(combrna$subcluster, levels = names(sclus_cols))
prop <- SCpubr::do_BarPlot(combrna, group.by = "subcluster",split.by = 'Pub_ID',position = "fill",flip = T,
                   order.by = rev(levels(combrna$Pub_ID)), plot.title = "Subcluster Proportion by Donor",
                   colors.use = sclus_cols, legend.position = "right")
ggsave("proportion_predicted.subcluster_v2.pdf",plot = prop, path = rpath, width = 12, height=8)

combrna$predicted.subclass_f <- factor(as.factor(combrna$predicted.subclass), levels = names(sclass_cols_orig))
SCpubr::do_BarPlot(combrna, group.by = "predicted.subclass_f",split.by = 'Pub_ID',position = "fill",flip = F,legend.position = "right",
                   order.by = rev(levels(combrna$Pub_ID)),plot.title = "Individual-level Cell type Proportion",
                   colors.use = sclass_cols_orig)
ggsave("proportion_predicted.subclass_v2.pdf", path = rpath, width = 10, height=7)

SCpubr::do_BarPlot(combrna, group.by = "celltype",split.by = 'Pub_ID',position = "fill",flip = F,
                   order.by = rev(levels(combrna$Pub_ID)),plot.title = "Individual-level Cell type Proportion",
                   colors.use = ct_cols)
ggsave("proportion_celltype_v2.pdf", path = rpath, width = 6, height=9)

SCpubr::do_BarPlot(combrna, group.by = "predicted.subclass_f",split.by = 'NPdiagnosis',position = "fill",flip = F,
                  plot.title = "Group Comparison: Cell type Proportion",
                   colors.use = sclass_cols_orig)
ggsave("proportion_celltype_npdiagnosis.pdf", path = rpath, width = 3, height=9)


# Prep index for identity levels
temp <- data.frame(combrna@meta.data[,c("celltype", "predicted.subclass", "subcluster")]) %>%
  `rownames<-`(NULL) %>%
  unique(.) %>% 
  .[order(.$subcluster),] %>% 
  `rownames<-`(1:nrow(.))
saveRDS(temp, paste0("../cell_identity_levels.Rds"))

### Area plots ####
combrna$CB_load <- ifelse(is.na(combrna$CB_load), 0, combrna$CB_load)
combrna$NFT_load <- ifelse(is.na(combrna$No_NFT_3.2e7mcsq), 0, combrna$No_NFT_3.2e7mcsq)
combrna@meta.data <- combrna@meta.data %>% mutate(
  NFT_load_dsct =  round(scales::rescale(.$NFT_load, to = c(0,4)),digits = 0)) %>%
  mutate(NFT_load_dsct = ifelse(NFT_load_dsct >= 3, NFT_load_dsct-1,NFT_load_dsct)) %>%
  mutate(NFT_load_dsct = as.factor(NFT_load_dsct))
combrna$TA_load <- ifelse(is.na(combrna$TA_load), 0, combrna$TA_load)
combrna$predicted.subclass <- as.factor(combrna$predicted.subclass)

data_n <- combrna@meta.data %>%
  group_by(NFT_load_dsct, predicted.subclass) %>%
  summarise(n = n(), .groups = "keep") %>%
  group_by(NFT_load_dsct) %>%
  mutate(percentage = n / sum(n))

data_a <- combrna@meta.data %>%
  group_by(TA_load, predicted.subclass) %>%
  summarise(n = n(), .groups = "keep") %>%
  group_by(TA_load) %>%
  mutate(percentage = n / sum(n))

data_o <- combrna@meta.data %>%
  group_by(CB_load, predicted.subclass) %>%
  summarise(n = n(), .groups = "keep") %>%
  group_by(CB_load) %>%
  mutate(percentage = n / sum(n))

n <- ggplot(data_n, aes(x = NFT_load_dsct, y = percentage, fill = predicted.subclass, group = predicted.subclass)) +
  geom_area(aes(group =  predicted.subclass),alpha=0.75 , size=.5, colour="white") +
  scale_fill_manual(values = c(sclass_cols_orig)) +
  labs(x = "NFT Load", y = "Percentage") + 
  theme_bw() + theme(legend.position = "none")
t <- ggplot(data_a, aes(x = TA_load, y = percentage, fill = predicted.subclass, group = predicted.subclass)) +
  geom_area(aes(group =  predicted.subclass),alpha=0.75 , size=.5, colour="white") +
  scale_fill_manual(values = c(sclass_cols_orig)) +
  labs(x = "TA Load", y = "") +
  theme_bw() + theme(legend.position = "none")
o <- ggplot(data_o, aes(x = CB_load, y = percentage, fill = predicted.subclass, group = predicted.subclass)) +
  geom_area(aes(group =  predicted.subclass),alpha=0.75 , size=.5, colour="white") +
  scale_fill_manual(values = c(sclass_cols_orig)) +
  labs(x = "CB Load", y = "") +
  theme_bw()
n|t|o
ggsave("proportion_all_sclass.pdf",path = rpath, width = 15, height=4)


# Prep index file for somatic mutations inference
index_file <- colnames(Seurat::GetAssayData(combrna)) %>% as.data.frame() %>% 
  separate(col = '.', sep ='_', into = c('drop','barcode')) %>% 
  separate(col = 'barcode', sep ='-', into = c('bc','drop')) %>% .[,1] %>%
  data.frame(Index = ., Cell_type = gsub("[^[:alnum:]]", "_", combrna$predicted.subclass)) %>% unique  
write.table(index_file, ('../results/SComatic/index_file_subclasses_psp_rna2.csv'))

sessionInfo()


# 5 Cell type composition: sccomp ####
if (!dir.exists("../sccomp")) { dir.create("../sccomp") }
setwd("../sccomp")

library(sccomp)

combrna$NPdiagnosis <- as.factor(combrna$NPdiagnosis)
combrna$pmi <- as.numeric(combrna$pmi)

levels <- c("celltype","predicted.subclass","subcluster")

### 5.1 celltype ####
res <-
  combrna %>%
  sccomp_estimate( 
    formula_composition = ~ NPdiagnosis + pmi + Age, 
    formula_variability = ~ NPdiagnosis + pmi + Age,
    .sample =  Pub_ID, 
    .cell_group = celltype, 
    bimodal_mean_variability_association = TRUE,
    cores = 32,
    verbose = F,
    mcmc_seed = 100
  ) %>%
  sccomp_remove_outliers(verbose = F) %>%
  sccomp_test(test_composition_above_logit_fold_change = 0.1)
write_csv(res, "sccomp_celltype_restable.csv")
res_co <- res$count_data
names(res_co) <- res$celltype
saveRDS(res_co, "sccomp_celltype_restable_coutns.Rds")
plots <- plot(res)
plots$boxplot
ggsave("sccomp_celltype_boxplot.pdf", width = 6, height=7)
plots$credible_intervals_1D
ggsave(paste0("sccomp_celltype_credible_intervals_1D.pdf"), width = 6, height=5)
plots$credible_intervals_2D
ggsave(paste0("sccomp_celltype_credible_intervals_2D.pdf"), width = 12, height=7)

# celltypes
dat <- left_join(plots$credible_intervals_1D[[1]]$data, cell_identity_levels, by = "celltype") %>%
  mutate(celltype = factor(celltype, levels = names(ct_cols)))
write_csv(dat, "sccomp_celltype_structured_data.csv")

ggplot(dat, aes(x = celltype, y = effect, ymin = lower, ymax = upper, 
                colour = ifelse(FDR < 0.05, "True","False"))) +
  geom_hline(yintercept = 0.2, alpha = 0.25) +
  geom_hline(yintercept = -0.2, alpha = 0.25) +
  geom_point(position = "identity") +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  scale_color_manual(values = c("grey26","red")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  xlab("") +
  labs(colour="FDR < 0.05") 
ggsave(paste0("sccomp_celltype_credible_intervals_1D_modif2.pdf"), width = 5, height=3.5)

### 5.2 subclass ####
res_p.sc <-
  combrna %>%
  sccomp_estimate( 
    formula_composition = ~ NPdiagnosis + pmi + Age, 
    formula_variability = ~ NPdiagnosis + pmi + Age,
    .sample =  Pub_ID, 
    .cell_group = predicted.subclass, 
    bimodal_mean_variability_association = TRUE,
    cores = 32, 
    verbose = F,
    mcmc_seed = 100
  ) %>%
  sccomp_remove_outliers(verbose = F) %>%
  sccomp_test(test_composition_above_logit_fold_change = 0.1)
write_csv(res_p.sc, "sccomp_predicted.subclass_restable.csv")
res_co <- res_p.sc$count_data
names(res_co) <- res_p.sc$predicted.subclass
saveRDS(res_co, "sccomp_predicted.subclass_restable_coutns.Rds")
plots <- plot(res_p.sc)
plots$boxplot
ggsave("sccomp_predicted.subclass_boxplot.pdf", width = 6, height=7)
plots$credible_intervals_1D
ggsave(paste0("sccomp_predicted.subclass_credible_intervals_1D.pdf"), width = 6, height=5)
plots$credible_intervals_2D
ggsave(paste0("sccomp_predicted.subclass_credible_intervals_2D.pdf"), width = 12, height=7)

dat <- left_join(plots$credible_intervals_1D[[1]]$data, cell_identity_levels, by = "predicted.subclass") %>%
  mutate(predicted.subclass = factor(predicted.subclass, levels = names(sclass_cols_orig)))
write_csv(dat, "sccomp_predicted.subclass_structured_data.csv")

ggplot(dat, aes(x = predicted.subclass, y = effect, ymin = lower, ymax = upper, 
                colour = ifelse(FDR < 0.05, "True","False"))) +
  geom_hline(yintercept = 0.2, alpha = 0.25) +
  geom_hline(yintercept = -0.2, alpha = 0.25) +
  geom_point(position = "identity") +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  scale_color_manual(values = c("grey26","red")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  xlab("") +
  labs(colour="FDR < 0.05") 
ggsave(paste0("sccomp_celltype_credible_intervals_1D_modif2.pdf"), width = 5, height=3.5)

### 5.3 subcluster ####
res_scl <-
  combrna %>%
  sccomp_estimate( 
    formula_composition = ~ NPdiagnosis + pmi + Age, 
    formula_variability = ~ NPdiagnosis + pmi + Age,
    .sample =  Pub_ID, 
    .cell_group = subcluster, 
    bimodal_mean_variability_association = TRUE,
    cores = 32, 
    verbose = F,
    mcmc_seed = 99
  ) %>%
  sccomp_remove_outliers(verbose = F) %>%
  sccomp_test(test_composition_above_logit_fold_change = 0.2)
write_csv(res_scl, "sccomp_subcluster_restable_rev.csv")
res_co <- res_scl$count_data
names(res_co) <- res_scl$subcluster
saveRDS(res_co, "sccomp_subcluster_restable_counts_rev.Rds")
plots <- plot(res_scl)
plots$boxplot
ggsave("sccomp_subcluster_boxplot.pdf_rev", width = 12, height=10)
plots$credible_intervals_1D
ggsave(paste0("sccomp_subcluster_credible_intervals_1D_rev.pdf"), width = 6, height=10)
plots$credible_intervals_2D
ggsave(paste0("sccomp_subcluster_credible_intervals_2D_rev.pdf"), width = 12, height=7)


### 5.4 Modified Plotting ####
dat <- left_join(plots$credible_intervals_1D[[1]]$data, cell_identity_levels, by = "subcluster") %>%
  mutate(celltype = factor(celltype, levels = names(ct_cols)),
         subcluster = factor(subcluster, levels = names(sclus_cols)))
write_csv(dat, "sccomp_subcluster_structured_data.csv")

ggplot(dat, aes(x = subcluster, y = effect, ymin = lower, ymax = upper, 
                colour = ifelse(FDR < 0.05, "True","False"))) +
  geom_hline(yintercept = 0.2, alpha = 0.25) +
  geom_hline(yintercept = -0.2, alpha = 0.25) +
  geom_point(position = "identity") +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  scale_color_manual(values = c("grey26","red")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  xlab("") +
  labs(colour="FDR < 0.05") +
  facet_grid(.~celltype, drop = T, space = "free", scales = "free")
ggsave(paste0("sccomp_subcluster_credible_intervals_1D_modif2.pdf"), width = 10, height=3.5)

sessionInfo()