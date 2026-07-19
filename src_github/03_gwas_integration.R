# ---
## TITLE: Identification of brain cell populations with GWAS SNP enrichment with gchromVAR
## Author: "Nils Briel, Center for Neuropathology, LMU Munich & Dept. Neurology, USZ Zurich"
## Date: "2024-10-04"
# ---

#@ In this part, we determine population-specific GWAS SNP enrichment of brain cells in PSP and CBD frontal cortex samples.
#@ Therefore, we include data from two previously published and reproduced GWAS studies involving either PSP (https://www.niagads.org/summary-statistics-psp) or CBD individuals. 
  # download reference GWAS data from:
  # MONDO_0019037 --> PSP
  # MONDO_0022308 --> CBD
  # MONDO_0017276 --> FTD
  # MONDO_0004975 --> AD
  # MONDO_0005180 --> PD
  # MONDO_0004976 --> ALS

# Libraries
library(chromVAR)
library(gchromVAR)
library(SummarizedExperiment)
library(data.table)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg38)
library(data.table)
set.seed(123)
s
if(!dir.exists("../results/psp/gwas_integration")){
  dir.create("../results/psp/gwas_integration")
}
dir_gwas <- "../results/psp/gwas_integration"

#### 1 Prepare input ####
# 1.1 Prepare peaks
# First, we need to construct a SumamrizedExperiment object, fed with peak counts from previous peak calling in SnapATAC:
    
# take Peaks matrix from Seurat and prepare input with cluster assignment
x.sp <-  readRDS("../results/psp/combatac_psp.Rds")
temp <- rownames(x.sp@meta.data)
x.sp@meta.data <- left_join(x.sp@meta.data, sample_lib, by = "Sample")
rownames(x.sp@meta.data) <- temp
x.sp <- subset(x.sp, subset = meta.Pub_ID %in% grepl('PSP', x.sp$meta.Pub_ID))
x.sp_red <- x.sp[,grepl('PSP', x.sp$meta.Pub_ID)]
mat <- x.sp_red@assays$peaks

# exclude unmapped regions
newdata <- ! seqnames(mat@ranges) %in% c('s37d5', unique(as.vector(seqnames(mat@ranges[grepl('GL', seqnames(mat@ranges))]))))

# preaparation function
preparePeakMatFromXSP <- function(mat, 
                                  ct){
  # construct per cell type sums of count matrix
  x <- t(mat$data)
  l = list()
  for(i in levels(ct)){
    l[[i]] <- DelayedMatrixStats::colSums2(x, rows = ct==i)
    }
  x <- t(Reduce(cbind, l))
  rownames(x) <- levels(ct)
  compl_peaks <- which(Matrix::colSums(x) > 0)
  x <- x[,compl_peaks] %>% t()
  rowD <- mat@ranges %>% .[compl_peaks]
  
  # call functions from SummarizedExperiment to construct SE 
  SE <- SummarizedExperiment(assays = list(counts = x),
                             rowData = rowD, 
                             colData = DataFrame(names = colnames(x)))
  
  # ... and chromVAR to correct for hg38 GCbias
  SE <- suppressWarnings(addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg38))
  return(SE)
  print('SE object construciton successfull')
}

setwd(dir_gwas)
SE_tau <- preparePeakMatFromXSP(mat = mat, ct = x.sp_red$celltype)
SE_tau_scl <- preparePeakMatFromXSP(mat = mat, ct = as.factor(x.sp_red$predicted.id))
SE_tau@metadata$dis = 'PSP'
SE_tau_scl@metadata$dis = 'PSP'
saveRDS(SE_tau, 'SE_psp.Rds')
saveRDS(SE_tau_scl, 'SE_psp_scl.Rds')

# clean R.envir
rm(list=setdiff(ls(), c("SE_tau","SE_tau_scl")))
 
# 1.2 Prepare GWAS summary


# define function for input formatting - produce gchromVAR compatible GWAS summary input
 rearrGWASsum <- function(path, format, pvalue, beta = 0.1, chr=NULL, bp=NULL, psig=NULL){
  
  if(format == 'xlsx'){
    table <- readxl::read_xlsx(path) %>% 
      dplyr::select(chr, bp, psig)
  }else if(format == 'xls'){
    table <- readxl::read_xls(path) %>% 
      dplyr::select(chr, bp, psig)
  }else{
    table <- read_delim(file = path, 
                        "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
      dplyr::select('CHR_ID', 'CHR_POS', 'P-VALUE','REGION','OR or BETA','PVALUE_MLOG','MAPPED_GENE') %>%
      dplyr::filter(`P-VALUE` < pvalue & `OR or BETA` > beta) %>%
      separate(col = "MAPPED_GENE", sep = " - ",into = c("gene1","gene2")) %>%
      select(-`P-VALUE`)
  }
  colnames(table) = c('chr.1', 'location.bp.1', 'region','or','pval_mlog', 'gene1', "gene2")
  table$chr.1 = as.factor(paste0('chr', table$chr.1))
  table$loc2 = as.integer(as.integer(table$location.bp.1)+1)
  table$location.bp.1 = as.integer(table$location.bp.1)
  table <- table %>% dplyr::select('chr.1','location.bp.1','loc2','region','or','pval_mlog','gene1','gene2')     
  table$region = as.factor(table$region)
  colnames(table) = c('V1', 'V2', 'V3', 'V4', 'V5','V6','gene1','gene2')
  table <- table[complete.cases(table),]
  return(table)
}
  pval <- 0.01
    # apply function on PSP GWAS
    GWASsumPSP <-  rearrGWASsum(path = list.files()[grep('MONDO_0019037',list.files())], 
                                format='tsv', 
                                pvalue = pval)
    write.table(GWASsumPSP, 'GWASsumPSP.bed')
    
    # apply function on CBD GWAS
    GWASsumCBD <-  rearrGWASsum(list.files()[grep('MONDO_0022308',list.files())], 
                                format='tsv', 
                                pvalue = pval)
    write.table(GWASsumCBD, 'GWASsumCBD.bed')
    
    # apply function on FTD GWAS
    GWASsumFTD <-  rearrGWASsum(list.files()[grep('MONDO_0017276',list.files())], 
                                format='tsv', 
                                pvalue = pval)
    write.table(GWASsumFTD, 'GWASsumFTD.bed')
    
    # apply function on AD GWAS
    GWASsumAD <-  rearrGWASsum(list.files()[grep('MONDO_0004975',list.files())], 
                               format='tsv', 
                               pvalue = pval)
    write.table(GWASsumAD, 'GWASsumAD.bed')
    
    # apply function on PD GWAS
    GWASsumPD <-  rearrGWASsum(list.files()[grep('MONDO_0005180',list.files())], 
                               format='tsv', 
                               pvalue = pval)
    write.table(GWASsumPD, 'GWASsumPD.bed')

    # apply function on ALS GWAS
    GWASsumALS <-  rearrGWASsum(list.files()[grep('MONDO_0004976',list.files())], 
                                format='tsv', 
                                pvalue = pval)
    write.table(GWASsumALS, 'GWASsumALS.bed')
    
    # clean R.envir
    rm(list=setdiff(ls(), c("SE_tau")))
    files <- list.files(full.names = TRUE, pattern = ".bed$")
    head(read.table(files[6]))
    
    # Read in all the .bed files and combine them into a single dataframe
    all_gwas_data <- lapply(files, function(file) {
      data <- read_tsv(file)
      data$Disease <- gsub("GWASsum(.+)\\.bed", "\\1", basename(file))
      return(data)
    }) %>% bind_rows() %>%
      mutate(Disease = factor(Disease, levels = rev(c('PSP','CBD','FTD','AD','ALS','PD'))))
  
#### 2 Construct weight scores ####
    bedscores <- importBedScore(rowRanges(SE_tau), files, colidx = 5)
    bedscores_scl <- importBedScore(rowRanges(SE_tau_scl), files, colidx = 5)
    
#Compute weights deviation:
    wDEV <- computeWeightedDeviations(SE_tau, bedscores)
    zdf <- reshape2::melt(t(assays(wDEV)[["z"]]))
    zdf[,2] <- gsub("_PP001", "", zdf[,2])
    colnames(zdf) <- c("ct", "tr", "Zscore")
    zdf <- zdf[complete.cases(zdf),]
 
#### 3 Lineage-specific enrichment test ####
## Annotate your single cell lineages..
    Ast <- c('Astro')
    Mic <- c('Micro_PVM')
    Neu <- c('Exc_DLN','Exc_ULN','Inh_Neu')
    Oli <- c('Oligo' , 'OPC')
 
# annotate df lineage-specific (LS) enrichments:
zdf$LS  <-
  (zdf$tr == "GWASsumPSP" & zdf$ct %in% Ast) +
  (zdf$tr == "GWASsumPSP" & zdf$ct %in% Mic) + 
  (zdf$tr == "GWASsumPSP" & zdf$ct %in% Neu) +
  (zdf$tr == "GWASsumPSP" & zdf$ct %in% Oli) + 
  (zdf$tr == "GWASsumCBD" & zdf$ct %in% Ast) +
  (zdf$tr == "GWASsumCBD" & zdf$ct %in% Mic) + 
  (zdf$tr == "GWASsumCBD" & zdf$ct %in% Neu) +
  (zdf$tr == "GWASsumCBD" & zdf$ct %in% Oli) + 
  (zdf$tr == "GWASsumFTD" & zdf$ct %in% Ast) +
  (zdf$tr == "GWASsumFTD" & zdf$ct %in% Mic) + 
  (zdf$tr == "GWASsumFTD" & zdf$ct %in% Neu) +
  (zdf$tr == "GWASsumFTD" & zdf$ct %in% Oli) + 
  (zdf$tr == "GWASsumAD" & zdf$ct %in% Ast) +
  (zdf$tr == "GWASsumAD" & zdf$ct %in% Mic) + 
  (zdf$tr == "GWASsumAD" & zdf$ct %in% Neu) +
  (zdf$tr == "GWASsumAD" & zdf$ct %in% Oli) + 
  (zdf$tr == "GWASsumPD" & zdf$ct %in% Ast) +
  (zdf$tr == "GWASsumPD" & zdf$ct %in% Mic) + 
  (zdf$tr == "GWASsumPD" & zdf$ct %in% Neu) +
  (zdf$tr == "GWASsumPD" & zdf$ct %in% Oli) + 
  (zdf$tr == "GWASsumLBD" & zdf$ct %in% Ast) +
  (zdf$tr == "GWASsumLBD" & zdf$ct %in% Mic) + 
  (zdf$tr == "GWASsumLBD" & zdf$ct %in% Neu) +
  (zdf$tr == "GWASsumLBD" & zdf$ct %in% Oli) +
  (zdf$tr == "GWASsumALS" & zdf$ct %in% Ast) +
  (zdf$tr == "GWASsumALS" & zdf$ct %in% Mic) + 
  (zdf$tr == "GWASsumALS" & zdf$ct %in% Neu) +
  (zdf$tr == "GWASsumALS" & zdf$ct %in% Oli)
 
# Mann-Whitney Rank-sum statistic for relative enrichment
  set.seed(123)
  zdf$gchromVAR_pvalue <- pnorm(zdf$Zscore, lower.tail = FALSE)
  gchromVAR_ranksum <- sum(1:dim(zdf)[1]*zdf[order(zdf$gchromVAR_pvalue, decreasing = FALSE), "LS"])
  permuted <- sapply(1:10000, function(i) sum(1:dim(zdf)[1] * sample(zdf$LS, length(zdf$LS))))
  pnorm((mean(permuted) - gchromVAR_ranksum)/sd(permuted), lower.tail = FALSE)
   
#### 4 Extracting individual enrichment: ####
  zdf$adj.pvalue <- p.adjust(zdf$gchromVAR_pvalue,method = 'bonferroni')
  zdf$sign <- ifelse(zdf$adj.pvalue < .1, 'sign', 'n.s.') %>% as.factor()
  zdf$tr <- factor(zdf$tr, levels = rev(c('GWASsumPSP','GWASsumCBD','GWASsumFTD','GWASsumAD','GWASsumALS','GWASsumPD'))) 
  saveRDS(zdf, 'gchromVAR_z_df.Rds')

  gchrom <- ggplot(zdf, aes( x = ct, y = tr, fill = Zscore)) + 
    labs(fill = 'Z-score', y = '', x = '') + 
    colorspace::scale_fill_continuous_diverging(palette = 'Blue-Red 3') + 
    geom_tile(color = "grey22") + 
    geom_text(aes(label=signif(adj.pvalue, 3)), nudge_y = 0.15, size = 3, fontface=2) + 
    geom_text(aes(label=signif(gchromVAR_pvalue, 3)), nudge_y = -0.15, size = 3, fontface=1)+
    ggthemes::theme_tufte(base_family = 'Helvetica', base_size = 12) + 
    theme(legend.position = 'bottom', axis.text.x =element_text(angle = 30, hjust = 1))
  
#### 5 Calculate for Subclasses
  wDEV_scl <- computeWeightedDeviations(SE_tau_scl, bedscores_scl)
  zdf <- reshape2::melt(t(assays(wDEV_scl)[["z"]]))
  zdf[,2] <- gsub("_PP001", "", zdf[,2])
  colnames(zdf) <- c("ct", "tr", "Zscore")
  zdf <- zdf[complete.cases(zdf),]
  
    set.seed(123)
  zdf$gchromVAR_pvalue <- pnorm(zdf$Zscore, lower.tail = FALSE)
  gchromVAR_ranksum <- sum(1:dim(zdf)[1]*zdf[order(zdf$gchromVAR_pvalue, decreasing = FALSE), "LS"])
  permuted <- sapply(1:10000, function(i) sum(1:dim(zdf)[1] * sample(zdf$LS, length(zdf$LS))))
  pnorm((mean(permuted) - gchromVAR_ranksum)/sd(permuted), lower.tail = FALSE)
  
  #### 4 Extracting individual enrichment: ####
  zdf$adj.pvalue <- p.adjust(zdf$gchromVAR_pvalue,method = 'bonferroni')
  zdf$sign <- ifelse(zdf$adj.pvalue < .1, 'sign', 'n.s.') %>% as.factor()
  zdf$tr <- factor(zdf$tr, levels = rev(c('GWASsumPSP','GWASsumCBD','GWASsumFTD','GWASsumAD','GWASsumALS','GWASsumPD'))) 
  saveRDS(zdf, 'gchromVAR_z_df_scl.Rds')
  
  gchrom_scl <- ggplot(zdf, aes( x = ct, y = tr, fill = Zscore)) + 
    labs(fill = 'Z-score', y = '', x = '') + 
    colorspace::scale_fill_continuous_diverging(palette = 'Blue-Red 3') + 
    geom_tile(color = "grey22") + 
    geom_text(aes(label=signif(adj.pvalue, 3)), nudge_y = 0.15, size = 3, fontface=2) + 
    geom_text(aes(label=signif(gchromVAR_pvalue, 3)), nudge_y = -0.15, size = 3, fontface=1)+
    ggthemes::theme_tufte(base_family = 'Helvetica', base_size = 12) + 
    theme(legend.position = 'bottom', axis.text.x = element_text(angle = 30, hjust = 1))
  
  histogram <- ggplot(all_gwas_data, aes(x = Disease, fill = Disease)) +
    geom_bar(color = "grey22", width = 0.8) +
    coord_flip() +
    theme_test() +
    labs(title = "",x = "", y = "#GWAS Risk Vars") + 
    scale_fill_brewer(palette = "Set1", direction = -1)
  
  pdf(('heatmap_gchrom_scl.pdf'), width = 20, height = 5)
  plot( plot_grid(align = "h",plotlist = list(gchrom,gchrom_scl,histogram),ncol = 3, rel_widths = c(2,3.7,1)))
  dev.off()
  
  
sessionInfo()
   