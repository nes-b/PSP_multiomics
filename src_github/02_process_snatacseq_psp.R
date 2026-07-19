# ---
# title: "snATAC-seq processing
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

blas_set_num_threads(12)
future::plan("multicore", workers = 16)
options(future.globals.maxSize= 64000*1024^2) #64 000 MB
rlimit_as(1e512)

source('../src/utils.R')
rpath <- "../results/psp/"
if(!dir.exists(rpath)){
  dir.create(rpath)
}

# snATAC-seq - Processing ####
# 1 QC & Metadata ####
meta <- readxl::read_excel("../data/case_metadata_coh1_2.xlsx", n_max = 24) %>% .[-5,-c(2,7,21,22,23,24,25)] %>% 
  dplyr::filter(NPdiagnosis  %in% c("PSP","Ctrl")) %>%
  mutate(Sample = paste('num', .$Code_ID, sep=''))
sample_lib <- data.frame(Sample = meta$Sample, meta$Pub_ID) %>% 
  `rownames<-`(.$Sample) %>% 
  .[match(meta$Sample, .$Sample), ]

if(!file.exists("../results/psp/combatac_psp.Rds")){
  
  ## 1.1 Data loading ####
  # Use lapply to read in the data for each sample and create a Seurat object
  if(!file.exists("../results/psp/combatac_wo_motifs_psp.Rds")){
    paths <- c(
      "~/multiome_psp_cbd/data/fastqs_atac/fastq/num102/outs/",
      "~/multiome_psp_cbd/data/fastqs_atac/fastq/num105/outs/",
      "~/multiome_psp_cbd/data/fastqs_atac/fastq/num112/outs/",
      "~/multiome_psp_cbd/data/fastqs_atac/fastq/num114/outs/",
      "~/multiome_psp_cbd/data/fastqs_atac/fastq/num115/outs/",
      "~/multiome_psp_cbd/data/fastqs_atac/Nils_fastq_Cohort2/num122/outs/",
      "~/multiome_psp_cbd/data/fastqs_atac/Nils_fastq_Cohort2/num123/outs/",
      "~/multiome_psp_cbd/data/fastqs_atac/Nils_fastq_Cohort2/num124/outs/",
      "~/multiome_psp_cbd/data/fastqs_atac/Nils_fastq_Cohort2/num125/outs/",
      "~/multiome_psp_cbd/data/fastqs_atac/fastq/numKontrolle_1_SFG/outs/",
      "~/multiome_psp_cbd/data/fastqs_atac/fastq/numKontrolle_2_SFG/outs/",
      "~/multiome_psp_cbd/data/fastqs_atac/fastq/numKontrolle_3_SFG/outs/",
      "~/multiome_psp_cbd/data/fastqs_atac/fastq/numKontrolle_4_SFG/outs/",
      "~/multiome_psp_cbd/data/fastqs_atac/fastq/numKontrolle_5_SFG/outs/",
      "~/multiome_psp_cbd/data/fastqs_atac/Nils_fastq_Cohort2/numC6/outs/",
      "~/multiome_psp_cbd/data/fastqs_atac/Nils_fastq_Cohort2/numC7/outs/",
      "~/multiome_psp_cbd/data/fastqs_atac/Nils_fastq_Cohort2/numC8/outs/"
    ) 
    
    # create common peak set
    bed_list <- lapply(paths, function(i) {
      print(i)
      peaks <- read.table(file = paste0(i,"peaks.bed"), col.names = c("chr", "start", "end"))
      peaks <- makeGRangesFromDataFrame(peaks)
    })
    combined.peaks <- Signac::reduce(x = c(bed_list[[1]], bed_list[[2]], bed_list[[3]], bed_list[[4]], bed_list[[5]], 
                                           bed_list[[6]], bed_list[[7]], bed_list[[8]], bed_list[[9]], bed_list[[10]],
                                           bed_list[[11]], bed_list[[12]], bed_list[[13]],bed_list[[14]],bed_list[[15]],
                                           bed_list[[16]] ,bed_list[[17]]))
    
        # Filter 
    peakwidths <- width(combined.peaks)
    combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
    samples <- sample_lib$Sample
    
    # Load the ATAC-seq data
    atac_list <- mclapply(1:17, mc.cores = 16, FUN = function(k) {
      print(samples[k])
      i <- paths[k]
      print(i)
      metadata <- read.csv(file = paste0(i, "singlecell.csv"),  header = TRUE,  row.names = 1)
      metadata <- metadata[metadata$passed_filters > 500 &  metadata$peak_region_fragments > 1000, ]
      metadata$Sample <- samples[k]
      frags <- CreateFragmentObject(path =paste0(i,"fragments.tsv.gz"),cells = rownames(metadata))
      counts <- FeatureMatrix(fragments = frags, features = combined.peaks, cells = rownames(metadata))
      assay <- CreateChromatinAssay(counts, fragments = frags, min.cells = 10, min.features = 200)
      seu_obj <- CreateSeuratObject(assay, assay = "peaks", meta.data=metadata)
      
      ## 1.2 Doublet filtering ####
      # snATAC-adjusted method 
      set.seed(123)
      sce <- scDblFinder(
        SingleCellExperiment(
          list(counts=counts)),
        aggregateFeatures=TRUE, 
        nfeatures=25, processing="normFeatures"
      )
      seu_obj[["scDblFinder_score"]] <- sce$scDblFinder.score
      seu_obj[["scDblFinder_class"]] <- sce$scDblFinder.class
      rm(sce)
      return(seu_obj)
    }
    )
    gc()
    
    ## Add meta.data and merge data
    names(atac_list) <- sample_lib$meta.Pub_ID
    
    library(conflicted) ### call here and not earlier
    conflicted::conflicts_prefer(Seurat::Assays)
    
    combatac <- merge(
      x = atac_list[[1]],
      y = atac_list[2:length(samples)],
      add.cell.ids = sample_lib$meta.Pub_ID
    )
    
    combatac$Pub_ID <- rownames(combatac@meta.data) %>% 
      as.data.frame() %>% 
      separate('.', c("a", "b"), sep = '_') %>% 
      .[,1] %>% as.factor
    detach("package:conflicted", unload = T)
    
    ## 1.3 QC ####
    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
    seqlevelsStyle(annotations) <- 'UCSC'
    Annotation(combatac) <- annotations
    combatac <- NucleosomeSignal(object = combatac)
    combatac <- TSSEnrichment(object = combatac, fast = F)
    combatac$pct_reads_in_peaks <- combatac$peak_region_fragments / combatac$passed_filters * 100
    combatac$blacklist_ratio <- combatac$blacklist_region_fragments / combatac$peak_region_fragments
    combatac$high.tss <- ifelse(combatac$TSS.enrichment > 2, 'High', 'Low')
    combatac$nucleosome_group <- ifelse(combatac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
    TSSPlot(combatac, group.by = 'high.tss') + theme_test() + NoLegend()
    ggsave("TSSPlot.pdf",path = rpath, width = 6, height=4)
    FragmentHistogram(object = combatac, group.by = 'nucleosome_group') + theme_test()
    ggsave("FragmentHistogram.pdf",path = rpath, width = 6, height=4)
    
    # QC filtering
    atac <- subset(
      x = combatac,
      subset = peak_region_fragments > 2000 &
        peak_region_fragments < 20000 &
        pct_reads_in_peaks > 10 &
        blacklist_ratio < 0.05 &
        nucleosome_signal < 4 &
        TSS.enrichment > 2 &
        scDblFinder_score < 0.5
    )
    
    # 2 Data Processing 
    ## 2.1 DimReduct & Normalization  ####
    atac <- RunTFIDF(atac)
    atac <- FindTopFeatures(atac, min.cutoff = 'q0')
    atac <- RunSVD(atac)
    atac <- RunUMAP(object = atac, reduction = 'lsi', dims = 2:30, reduction.name = "umap.atac", reduction.key = "atac_UMAP_")
    p1 <- DimPlot(atac, reduction = "umap.atac", raster = T, size = 1.5, group.by = 'Pub_ID') + theme_test() + coord_equal()
    rm(combatac, annotations, atac_list, bed_list, sample_lib)
    gc()
    
    ## 2.2 Batch Effekt Correction ####
    library(harmony)
    combatac <- RunHarmony(
      object = atac,
      group.by.vars = 'Pub_ID',
      reduction.use = 'lsi',
      assay.use = 'peaks',
      project.dim = FALSE
    )
    
    combatac <- RunUMAP(object = combatac, reduction = 'harmony', dims = 2:30, reduction.name = "umap.atac.ha")
    p2 <- DimPlot(combatac, reduction = "umap.atac.ha", group.by = 'Pub_ID', raster = T, size = 1.5) + theme_test() + coord_equal()
    p1 + p2
    ggsave("DimPlot_atac_harmony.pdf",path = rpath, width = 8, height=4)
    
    
    # 3 Gene Activity ####
    set.seed(42)
    blas_set_num_threads(32)
    combrna <- readRDS("../results/psp/combrna_psp_scl2.Rds")
    DefaultAssay(combrna) <- "RNA"
    combrna <- FindVariableFeatures(combrna, selection.method = "vst", nfeatures = 3000)
    gene.activities <- GeneActivity(combatac, assay = "peaks", features = VariableFeatures(combrna))
    
    # add the gene activity matrix to the Seurat object as a new assay and normalize it
    combatac[['GA']] <- CreateAssayObject(counts = gene.activities)
    DefaultAssay(combatac) <- "GA"
    combatac <- NormalizeData(combatac)
    combatac <- ScaleData(combatac, vars.to.regress = c('nCount_peaks',"mitochondrial"))
    saveRDS(combatac, "../results/psp/combatac_wo_motifs_psp.Rds")
  }else{
    combatac <- readRDS("../results/psp/combatac_wo_motifs_psp.Rds")
  }
  
  # 4 Motif Analysis ####
  library(JASPAR2022) # if loading error: devtools::install_version("dbplyr", version = "2.3.4"), appeared on 2024-01-07
  library(TFBSTools)
  library(patchwork)
  library(BiocParallel)
  register(MulticoreParam(32))
  set.seed(1234)
  
  # https://github.com/stuart-lab/signac/issues/780
  DefaultAssay(combatac) <- "peaks"
  gr <- granges(combatac)
  seq_keep <- seqnames(gr) %in% seqnames(BSgenome.Hsapiens.UCSC.hg38) 
  seq_keep <- as.vector(seq_keep)
  feat.keep <- GRangesToString(grange = gr[seq_keep])
  combatac[['peaks']] <- subset(combatac[["peaks"]], features = feat.keep)
  
  library(conflicted) 
  # Get a list of motif position frequency matrices from the JASPAR database
  pfm <- getMatrixSet(
    x = JASPAR2022, 
    opts = list(collection = "CORE", tax_group = 'vertebrates', species = 9606, all_versions = FALSE))
  
  # add motif information
  combatac <- AddMotifs(
    object = combatac,
    genome = BSgenome.Hsapiens.UCSC.hg38,
    pfm = pfm
  )
  
  conflicts_prefer(MatrixGenerics::colSums2)
  register(SerialParam())
  combatac <- RunChromVAR(
    object = combatac,
    genome = BSgenome.Hsapiens.UCSC.hg38
  )
  
  rm(gene.activities,gr,seq_keep,feat.keep,pfm)
  gc()

  # 5 Reference Integration & Celltype Annotation ####
  transfer.anchors <- FindTransferAnchors(reference = combrna, query = combatac,
                                          features = VariableFeatures(object = combrna),
                                          reference.assay = "RNA", query.assay = "GA",
                                          reduction = "cca")

  celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = combrna$predicted.subclass,
                                       weight.reduction = combatac[["lsi"]], dims = 2:30)

  combatac <- AddMetaData(combatac, metadata = celltype.predictions)

  # less granular celltypes
  ct_lib <- data.frame(predicted.id = c("Astro","Endo", "L2/3 IT", "L5 ET", "L5 IT", "L5/6 NP", "L6 CT","L6 IT","L6 IT Car3",
                                        "L6b","Lamp5","Micro-PVM","Oligo","OPC","Pvalb", "Sncg","Sst","Sst Chodl","Vip","VLMC"),
                       celltype = as.factor(c("Astro","Endo_VLMC", "Exc_ULN", "Exc_DLN", "Exc_DLN", "Exc_DLN", "Exc_DLN","Exc_DLN","Exc_DLN",
                                              "Exc_DLN","Inh_Neu","Micro_PVM","Oligo","OPC","Inh_Neu", "Inh_Neu",
                                              "Inh_Neu","Inh_Neu","Inh_Neu","Endo_VLMC"))
  )
  combatac@meta.data <- left_join(combatac@meta.data, ct_lib, by="predicted.id")
  rownames(combatac@meta.data) = colnames(combatac@assays$peaks)

  
  a1 <- do_DimPlot(combatac, group.by = "predicted.id", plot.title = "Predicted subclass",
                   pt.size = 1, raster.dpi = 400, shuffle = T, raster = T, colors.use = sclass_cols_orig,
                   reduction = "umap.atac.ha", label = T,  label.box = F, repel = T) + 
    theme_test() + coord_equal()
  a2 <- do_DimPlot(combatac, group.by = "celltype", plot.title = "Celltype",
                   pt.size = 1, raster.dpi = 400, shuffle = T, raster = T, colors.use = ct_cols,
                   reduction = "umap.atac.ha", label = T,  label.box = F, repel = T) + 
    theme_test() + coord_equal()
  a2|a1
  ggsave("DimPlot_annotated_atac_2025.pdf",path = rpath, width = 12, height=6)

  p1 <- FeaturePlot(combatac, features = "prediction.score.max",
                    reduction = "umap.atac.ha",
                    label = F) +
    colorspace::scale_color_continuous_sequential(palette = "OrRd", limits = c(0,1.01)) +
    theme_test() + coord_equal()
  
  # subset based on prediction score
  combatac <- subset(combatac, subset = prediction.score.max > 0.5)
  p2 <- FeaturePlot(combatac, features = "prediction.score.max",
                    reduction = "umap.atac.ha",
                    label = F) +
    colorspace::scale_color_continuous_sequential(palette = "OrRd", limits = c(0,1.01)) +
    theme_test() + coord_equal()
  p1|p2
  ggsave("prediction.score.max_atac.pdf",path = rpath, width = 8, height=4)
  
  # saving
  saveRDS(combatac, "../results/psp/combatac_psp.Rds")
}else{
  combatac <- readRDS("../results/psp/combatac_psp.Rds")
  combrna <- readRDS("../results/psp/combrna_psp_scl2.Rds")
}

# 6 Co-embedding snRNA- & snATAC-seq ####
set.seed(42)
blas_set_num_threads(16)
genes.use <- VariableFeatures(combrna)
refdata <- GetAssayData(combrna, assay = "RNA", layer = "data")[genes.use,]
DefaultAssay(combrna) <- "RNA"
transfer.anchors <- FindTransferAnchors(reference = combrna, query = combatac, 
                                        features = VariableFeatures(object = combrna),
                                        reference.assay = "RNA", query.assay = "GA", 
                                        reduction = "cca")
combatac[["RNA"]] <- TransferData(anchorset = transfer.anchors, refdata = refdata, 
                                  weight.reduction = combatac[["lsi"]], dims = 2:30)
coembed <- merge(x = combrna, y = combatac)
coembed <- ScaleData(coembed, do.scale = FALSE, features = genes.use)
coembed <- RunPCA(coembed, verbose = FALSE, features = genes.use)
coembed <- RunUMAP(coembed, dims = 1:30)
coembed$tech <- ifelse(is.na(coembed$predicted.id), "rna", "atac")
coembed$predicted.subclass <- ifelse(is.na(coembed$predicted.id), coembed$predicted.subclass, coembed$predicted.id) %>% as.factor()

p1 <- DimPlot(coembed, group.by = c("predicted.id"), reduction = "umap", cols = sclass_cols_orig, label = T, repel = T, raster = T, size = 1.5) + 
  theme_test() + NoLegend() + coord_fixed(1/1.2) + ggtitle("Predicted Cell Types in snATAC nuclei")
p2 <- DimPlot(coembed, group.by = c("celltype"), reduction = "umap", cols = sclass_cols_orig, label = T, repel = T, raster = T, size = 1.5) + 
  theme_test() + NoLegend() + coord_fixed(1/1.2) + ggtitle("Combined Cell Types")
p2.1 <- DimPlot(coembed, group.by = c("predicted.subclass"), reduction = "umap", cols = sclass_cols_orig, label = T, repel = T, raster = T, size = 1.5) + 
  theme_test() + NoLegend() + coord_fixed(1/1.2) + ggtitle("Combined Subclass")
p3 <- DimPlot(coembed, group.by = c("tech"), reduction = "umap", raster = T, size = 1.5,) +  
  theme_test() + coord_fixed(1/1.2) + ggtitle("Technology")
p4 <- DimPlot(coembed, group.by = c("NPdiagnosis"), reduction = "umap", cols = group_cols, raster = T, size = 1.5) + 
  theme_test() + coord_fixed(1/1.2) + ggtitle("Diagnosis")


# Predicted Cell Types in snATAC nuclei
p1 <- do_DimPlot(coembed, 
                 group.by = "predicted.id", 
                 plot.title = "Predicted Cell Types in snATAC nuclei",
                 pt.size = 1.5, 
                 raster.dpi = 256, 
                 shuffle = TRUE, 
                 raster = TRUE, 
                 colors.use = sclass_cols_orig,
                 reduction = "umap", 
                 label = TRUE, 
                 label.box = FALSE, 
                 repel = TRUE) +
  theme_test() + 
  coord_fixed(1/1.2)

# Combined Cell Types
p2 <- do_DimPlot(coembed, 
                 group.by = "celltype", 
                 plot.title = "Combined Cell Types",
                 pt.size = 1.5, 
                 raster.dpi = 256, 
                 shuffle = TRUE, 
                 raster = TRUE, 
                 colors.use = ct_cols,
                 reduction = "umap", 
                 label = TRUE, 
                 label.box = FALSE, 
                 repel = TRUE) +
  theme_test() + 
  coord_fixed(1/1.2)

# Combined Subclass
p2.1 <- do_DimPlot(coembed, 
                   group.by = "predicted.subclass", 
                   plot.title = "Combined Subclass",
                   pt.size = 1.5, 
                   raster.dpi = 256, 
                   shuffle = TRUE, 
                   raster = TRUE, 
                   colors.use = sclass_cols_orig,
                   reduction = "umap", 
                   label = TRUE, 
                   label.box = FALSE, 
                   repel = TRUE) +
  theme_test() + 
  coord_fixed(1/1.2)

# Technology
p3 <- do_DimPlot(coembed, 
                 group.by = "tech", 
                 plot.title = "Technology",
                 pt.size = 1.5, 
                 raster.dpi = 256, 
                 shuffle = TRUE, 
                 raster = TRUE,
                 reduction = "umap") +
  theme_test() + 
  coord_fixed(1/1.2)

# Diagnosis
p4 <- do_DimPlot(coembed, 
                 group.by = "NPdiagnosis", 
                 plot.title = "Diagnosis",
                 pt.size = 1.5, 
                 raster.dpi = 256, 
                 shuffle = TRUE, 
                 raster = TRUE, 
                 colors.use = group_cols,
                 reduction = "umap") +
  theme_test() + 
  coord_fixed(1/1.2)



p1 | p3
p2 | p2.1
p3 | p4
cowplot::plot_grid(plotlist = list(p1,p3,p2,p2.1,p3,p4), ncol = 2)
ggsave("DimPlot_coembed_3.pdf", path = rpath, width = 15, height=17)

# saving
saveRDS(coembed, "../results/psp/coembed.Rds")

sessionInfo()
