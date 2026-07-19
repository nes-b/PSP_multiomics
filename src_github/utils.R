# ---
# title: Helper function script
# project: "snRNA+ATAC-seq in PSP frontal cortex"
# author: "Nils Briel, Center for Neuropathology, LMU Munich & Dept. Neurology, USZ Zurich"
# date: "Last modified: 2024-12-17"
# ---

# general libs
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

# define colors:
group_cols <- c("Ctrl" = "grey62", "PSP" = "red3") # "CBD" = "#7370ff",
group_cols_c <- c("Ctrl" = "grey62", "CBD" = "#f7b577")

ct_cols <- c("Exc_ULN" = '#f0c571',
             "Exc_DLN" = '#36b700',
             "Inh_Neu" = '#7e4974',
             "Astro" = '#e25759',
             "Micro_PVM" = "#c9c8c1",
             "Oligo" = '#598a9c',
             "OPC" = '#0b81a2',
             "Endo_VLMC" = "#9d2c00")
ct_cols_hyph <- ct_cols
names(ct_cols_hyph) <- names(ct_cols_hyph) %>% stringi::stri_replace_all_regex("_","-")


# subclass
sclass_cols <- c("L2_3_IT" = "#f2b774", 
                 "L5_ET" = "#f0c571", 
                 "L5_IT" = "#f2e574", 
                 "L5_6_NP" = '#afdb85', 
                 "L6_CT" = '#85db8b', 
                 "L6_IT" = '#2ab553',
                 "L6_IT_Car3" = "#36b700",
                 "L6b" = "#6fa374",
                 "Lamp5" = '#7e4974',
                 "Pvalb" = "#76a0e3",
                 "Sncg" = '#e65cd5',
                 "Sst" = '#9970ff',
                 "Sst_Chodl" = '#7370ff',
                 # "Sst_Chodl.NA.Na" = '#7370ff',
                 "Vip" = '#e695dc',
                 "Astro" = '#e25759',
                 "Micro-PVM" = "#c9c8c1",
                 "Oligo" = '#598a9c',
                 "OPC" = '#0b81a2',
                 "Endo_VLMC" = "#9d2c00",
                 "Endo" = "#9a4c00",
                 "VLMC" = "#9a6c00")

sclass_cols_hyph <- sclass_cols
names(sclass_cols_hyph) <- names(sclass_cols_hyph) %>% stringi::stri_replace_all_regex("_","-")

# subclass orig
sclass_cols_orig <- c("L2/3 IT" = "#f2b774", 
                      "L5 ET" = "#f0c571", 
                      "L5 IT" = "#f2e574", 
                      "L5/6 NP" = '#afdb85', 
                      "L6 CT" = '#85db8b', 
                      "L6 IT" = '#2ab553',
                      "L6 IT Car3" = "#36b700",
                      "L6b" = "#6fa374",
                      "Lamp5" = '#7e4974',
                      "Pvalb" = "#76a0e3",
                      "Sncg" = '#e65cd5',
                      "Sst" = '#9970ff',
                      # "Sst_Chodl" = '#7370ff',
                      "Sst_Chodl.NA.NA" = '#7370ff',
                      "Vip" = '#e695dc',
                      "Astro" = '#e25759',
                      "Micro-PVM" = "#c9c8c1",
                      "Oligo" = '#598a9c',
                      "OPC" = '#0b81a2',
                      "Endo_VLMC" = "#9d2c00",
                      "Endo" = "#9a4c00",
                      "VLMC" = "#9a6c00")


scl_cols <- c(
  "#e84d4d","#ba3434","#c24c4c","#CB181D","#9d2c00","#b05e3e","#c78165","#8f482c","#c2974c","#f7ca77","#f7b577",
  "#FDAE6B","#FD8D3C","#F16913","#D94801","#A63603","#7F2704","#a8f777","#67b550","#bede73","#ADDD8E","#78C679",
  "#41AB5D","#238443","#006837","#004529","#bffca7","#86d16f","#C7E9C0","#A1D99B","#74C476","#41AB5D","#238B45",
  "#006D2C","#00441B","#6fa374","#afdb85","#8a9ae3","#6375c9","#B3CDE3","#8C96C6","#c9c8c1","#7acbf5","#53adf5",
  "#5077d9","#535bf5","#9ECAE1","#6BAED6","#4292C6","#2171B5","#08519C","#c7a1ed","#a266de","#565494","#8a88c2", 
  "#BCBDDC","#9E9AC8","#807DBA","#6A51A3","#54278F","#d18ecf","#c38ed1","#D7B5D8","#DF65B0","#DD1C77","#980043", 
  "#9d2c00","#b05e3e","#c78165","#8f482c")

# alphabet
sclus_cols_alph <- c("Astro.1" = "#f77474","Astro.2" = "#fa5043", "Astro.3" = "#f54949",  "Astro.4" = "#cc0c12",
                     "Endo_VLMC.1" = "#9d2c00","Endo_VLMC.2" = "#b05e3e",  "Endo_VLMC.3" = "#c78165",  "Endo_VLMC.4" = "#8f482c",
                     "L2_3_IT.1" = "#c2974c",  "L5_6_NP.1" = "#f7ca77",  "L5_ET.1" = "#f7b577",  
                     "L5_IT.1" = "#006837",  "L5_IT.2" = "#004529",  "L5_IT.3" = "#bffca7",  "L5_IT.4" = "#41AB5D",  "L5_IT.5" = "#238443",  "L5_IT.6" = "#41AB53",
                     "L5_IT.7" = "#a8f777",  "L5_IT.8" = "#67b550",  "L5_IT.9" = "#bede73",  "L5_IT.10" = "#ADDD8E",  "L5_IT.11" = "#78C679",
                     "L6_CT.1" = "#86d16f",  
                     "L6_IT_Car3.1" = "#C7E9C0",  "L6_IT.1" = "#A1D99B",  "L6_IT.2" = "#74C476",
                     "L6b.1" = "#00441B",  "L6b.2" = "#6fa374",  "L6b.3" = "#afdb85",  
                     "Lamp5.1" = "#8a9ae3",  "Lamp5.2" = "#6375c9",  "Lamp5.3" = "#B3CDE3",  "Lamp5.4" = "#8C96C6",
                     "Micro-PVM.1" = "#c9c8c1",
                     "Oligo.1" = "#7acbf5",  "Oligo.2" = "#0aa8fa",  "Oligo.3" = "#1c56e8",  "Oligo.4" = "#0915eb",
                     "OPC.1" = "#6BAED6",  "OPC.2" = "#4292C6",   "OPC.3" = "#2171B5",  "OPC.4" = "#08519C",
                     "Pvalb.1" = "#c7a1ed",  "Pvalb.2" = "#a266de",  "Sncg.1" = "#565494",
                     "Sst_Chodl.1" = "#8a88c2",  "Sst.1" = "#BCBDDC",  "Sst.2" = "#9E9AC8",  "Sst.3" = "#807DBA",  "Sst.4" = "#6A51A3",  "Sst.5" = "#54278F", 
                     "Vip.1" = "#d18ecf",  "Vip.2" = "#c38ed1",  "Vip.3" = "#D7B5D8",  "Vip.4" = "#DF65B0",  "Vip.5" = "#DD1C77",  "Vip.6" = "#980043"
                     #"VLMC.1" = "#9d2c00",  "VLMC.2" = "#b05e3e",  "VLMC.3" = "#c78165",  "VLMC.4" = "#8f482c"
)

# cell classes
sclus_cols <- c("L2_3_IT.1" = "#c2974c",  "L5_6_NP.1" = "#f7ca77",  "L5_ET.1" = "#f7b577",  
                "L5_IT.1" = "#006837",  "L5_IT.2" = "#004529",  "L5_IT.3" = "#bffca7",  "L5_IT.4" = "#41AB5D",  "L5_IT.5" = "#238443",  "L5_IT.6" = "#41AB53",
                "L5_IT.7" = "#a8f777",  "L5_IT.8" = "#67b550",  "L5_IT.9" = "#bede73",  "L5_IT.10" = "#ADDD8E",  "L5_IT.11" = "#78C679",
                "L6_CT.1" = "#86d16f",  
                "L6_IT_Car3.1" = "#C7E9C0",  "L6_IT.1" = "#A1D99B",  "L6_IT.2" = "#74C476",
                "L6b.1" = "#00441B",  "L6b.2" = "#6fa374",  "L6b.3" = "#afdb85",  
                "Lamp5.1" = "#8a9ae3",  "Lamp5.2" = "#6375c9",  "Lamp5.3" = "#B3CDE3",  "Lamp5.4" = "#8C96C6",
                "Pvalb.1" = "#c7a1ed",  "Pvalb.2" = "#a266de",  "Sncg.1" = "#565494",
                "Sst_Chodl.1" = "#8a88c2",  "Sst.1" = "#BCBDDC",  "Sst.2" = "#9E9AC8",  "Sst.3" = "#807DBA",  "Sst.4" = "#6A51A3",  "Sst.5" = "#54278F", 
                "Vip.1" = "#d18ecf",  "Vip.2" = "#c38ed1",  "Vip.3" = "#D7B5D8",  "Vip.4" = "#DF65B0",  "Vip.5" = "#DD1C77",  "Vip.6" = "#980043",
                "Astro.1" = "#f77474","Astro.2" = "#fa5043", "Astro.3" = "#f54949",  "Astro.4" = "#cc0c12",
                "Micro-PVM.1" = "#c9c8c1",
                "Oligo.1" = "#7acbf5",  "Oligo.2" = "#53adf5",  "Oligo.3" = "#5077d9",  "Oligo.4" = "#535bf5",
                "OPC.1" = "#6BAED6",  "OPC.2" = "#4292C6",   "OPC.3" = "#2171B5",  "OPC.4" = "#08519C",
                "Endo_VLMC.1" = "#9d2c00","Endo_VLMC.2" = "#b05e3e",  "Endo_VLMC.3" = "#c78165",  "Endo_VLMC.4" = "#8f482c"
                # "VLMC.1" = "#9d2c00",  "VLMC.2" = "#b05e3e",  "VLMC.3" = "#c78165",  "VLMC.4" = "#8f482c"
)

sclus_cols_hyph <- sclus_cols
names(sclus_cols_hyph) <- names(sclus_cols_hyph) %>% stringi::stri_replace_all_regex("_","-")

###
if(file.exists("~/multiome_psp_cbd/results/psp/cell_identity_levels.Rds")){
  cell_identity_levels <- readRDS("~/multiome_psp_cbd/results/psp/cell_identity_levels.Rds")
}

###
get_cols <- function(level = "subcluster", subset = "Sst.", nam = T){
  if(file.exists("~/multiome_psp_cbd/results/psp/cell_identity_levels.Rds")){
    temp <- readRDS("~/multiome_psp_cbd/results/psp/cell_identity_levels.Rds")
  }
  ls <- list(
    celltype = ct_cols_nam,
    predicted.subclass = ct_cols_exp_clean,
    subcluster = scl_cols_nam_cl
  )
  sel <- temp[grep(subset, temp[[level]]), level]
  col <- suppressWarnings(ls[[level]][names(ls[[level]]) %in% sel])
  return(col)
}

cutpoints <- c(0, 0.001, 0.01, 0.05, 0.1, 1)
symbols <- c("***", "**", "*", "(*)", "")

# Create a function to convert p-values to significance stars
p_to_stars <- function(p_value) {
  stars <- symbols[findInterval(p_value, cutpoints)]
  return(stars)
}


# background gradient function
make_gradient <- function(deg = 45, n = 100, cols = blues9) {
  library(ggplot2) 
  library(grid)
  library(RColorBrewer)
  cols <- colorRampPalette(cols)(n + 1)
  rad <- deg / (180 / pi)
  mat <- matrix(
    data = rep(seq(0, 1, length.out = n) * cos(rad), n),
    byrow = TRUE,
    ncol = n
  ) +
    matrix(
      data = rep(seq(0, 1, length.out = n) * sin(rad), n),
      byrow = FALSE,
      ncol = n
    )
  mat <- mat - min(mat)
  mat <- mat / max(mat)
  mat <- 1 + mat * n
  mat <- matrix(data = cols[round(mat)], ncol = n)
  grid::rasterGrob(
    image = mat,
    width = unit(1, "npc"),
    height = unit(1, "npc"), 
    interpolate = TRUE
  )
}


plot_correlation_matrix <- function(data_matrix, method = "pearson") {
  library(Hmisc)
  library(corrplot)
  
  # Compute correlation matrix with p-values
  cor_matrix <- rcorr(as.matrix(data_matrix), type = method)
  
  # Extract correlation coefficients and p-values
  r <- cor_matrix$r
  p <- cor_matrix$P
  
  # Create a matrix to store significant correlations
  sig_correlations <- matrix(ifelse(p < 0.05, round(r, 2), ""), ncol = ncol(r))
  
  # Plot the correlation matrix
  corrplot(r, method = "color", type = "upper", order = "hclust", 
           tl.col = "black", tl.srt = 45, 
           p.mat = p, sig.level = 0.05, insig = "blank",
           addCoef.col = "black", number.cex = 0.7, 
           addCoefasPercent = FALSE, 
           cl.pos = "n")
  
  # Add significant correlations as text
  text(x = row(r)[upper.tri(r)],
       y = col(r)[upper.tri(r)],
       labels = sig_correlations[upper.tri(sig_correlations)],
       cex = 0.7, pos = 3)
  
  # Add color legend
  corrplot::corrplot.mixed(r, order = "hclust", tl.pos = "n", cl.pos = "b")
}

###### umap_plot

umap_plot <-   function(obj, 
                        main = NULL,
                        method = 'UMAPHarmony', 
                        point.size = 1.5, 
                        point.alpha = 0.8, 
                        point.color = NULL,
                        point.color.text = NULL,
                        text.add = NULL, 
                        text.size = 4, 
                        legend.add = F,
                        legend.pos = c("bottomleft", "bottom", "left", "topleft", "top", "topright", "right", "center"),
                        feature.mode = NULL,
                        scale_fill = NULL,
                        stroke = 0.1,
                        threeD = F,
                        hex.bin = F,
                        bin.size = 100,
                        ...) 
{
  if(!threeD){
    if(class(obj) == "ArchRProject"){
      data.use = as.data.frame(obj@embeddings[[method]][[1]]);  
      ncell = nrow(data.use);
      
      colnames(data.use) = c("UMAP1", "UMAP2")
      xlims = c(-max(abs(data.use[,1])) * 1.05, max(abs(data.use[,1])) * 1.2);
      ylims = c(-max(abs(data.use[,2])) * 1.05, max(abs(data.use[,2])) * 1.05);
      
      if(is.null(feature.mode)){
        data.use$col = factor(obj@cellColData[,point.color]);
      }else{
        data.use$col = feature.mode
      }
    }
    
    legend.pos = match.arg(legend.pos);
    if(!hex.bin){
      pt_plot <- ggplot(data.use, 
                        aes(x = UMAP1, 
                            y= UMAP2, 
                            fill = col),color = 'grey42') + 
        geom_point(size = point.size, 
                   alpha = point.alpha, 
                   shape=21, stroke = stroke) + 
        labs(title = main, 
             x = 'UMAP1', 
             y = 'UMAP2 ', 
             fill = point.color.text)  +
        theme_void(base_family="Helvetica") + theme(plot.title = element_text(hjust = 0.5))
      
      if(!is.null(feature.mode)){
        pt_plot <- pt_plot + scale_color_viridis_c()
      }
      if(!is.null(scale_fill)){
        pt_plot <- pt_plot + scale_fill_manual(values = scale_fill)
      }
    }else{
      pt_plot <- ggplot(data.use, 
                        aes(x = UMAP1, 
                            y = UMAP2, 
                            fill = as.factor(col)),color = 'grey42') + 
        geom_hex(bins = bin.size, color = 'grey22', size = stroke) + 
        labs(title = main, 
             x = 'UMAP1', 
             y = 'UMAP2 ', 
             fill = point.color.text)  +
        theme_void(base_family="Helvetica") + theme(plot.title = element_text(hjust = 0.5)) + 
        scale_fill_manual(values = scale_fill)
    }
    
    
    if(legend.add){
      pt_plot <- pt_plot + theme(legend.position=legend.pos) 
    }else{
      pt_plot <- pt_plot + theme(legend.position = "none") 
    }
    
    if(text.add){
      centroid <- aggregate(cbind(data.use$UMAP1, data.use$UMAP2) ~ col, data=data.use, FUN=mean)
      pt_plot <-  pt_plot + geom_label(data = centroid, mapping = aes(x=V1, y=V2, label=factor(col)), color = 'black', size = text.size) 
    }
    return(pt_plot)
    
    
    
  }else{
    if(class(obj) == 'SingleCellExperiment'){
      paste('sce-plot')
      ncell = ncol(obj)
      data.use = as.data.frame(reducedDims(obj)$UMAP)
      
      if(is.null(feature.mode)){
        data.use$col = as.factor(obj@metadata[[point.color]])
      }else{
        data.use$col = feature.mode
      }
      
    }else if(class(obj) == "snap"){
      data.use = as.data.frame(obj@embeddings[[method]][[1]]);  
      ncell = nrow(data.use);
      colnames(data.use) = c("UMAP1", "UMAP2", "UMAP3")
      xlims = c(-max(abs(data.use[,1])) * 1.05, max(abs(data.use[,1])) * 1.2);
      ylims = c(-max(abs(data.use[,2])) * 1.05, max(abs(data.use[,2])) * 1.05);
      zlims = c(-max(abs(data.use[,2])) * 1.05, max(abs(data.use[,2])) * 1.05);
      if(is.null(feature.mode)){
        data.use$col = factor(obj@metaData[,point.color]);
      }else{
        data.use$col = feature.mode
      }
    }
    legend.pos = match.arg(legend.pos);
    
    pt_plot <- ggplot(data.use, 
                      aes(x = UMAP1, 
                          y= UMAP2, 
                          z = UMAP3,
                          fill = col),color = 'grey42') + 
      geom_point(size = point.size, 
                 alpha = point.alpha, 
                 shape=21, stroke = stroke) + 
      labs(title = main, 
           x = 'UMAP1', 
           y = 'UMAP2',
           z = 'UMAP2',
           fill = point.color.text)  +
      theme_void(base_family="Helvetica") + theme(plot.title = element_text(hjust = 0.5))
    
    if(!is.null(feature.mode)){
      pt_plot <- pt_plot + scale_color_viridis_c()
    }
    if(!is.null(scale_fill)){
      pt_plot <- pt_plot + scale_fill_manual(values = scale_fill)
    }
    
    if(legend.add){
      pt_plot <- pt_plot + theme(legend.position=legend.pos) 
    }else{
      pt_plot <- pt_plot + theme(legend.position = "none") 
    }
    
    if(text.add){
      centroid <- aggregate(cbind(data.use$UMAP1, data.use$UMAP2, data.use$UMAP3) ~ col, data=data.use, FUN=mean)
      pt_plot <-  pt_plot + geom_label(data = centroid, mapping = aes(x=V1, y=V2, z=V3, label=factor(col)), color = 'black', size = text.size) 
    }
    return(pt_plot)
  }
}


###### ct_freqs_plot  

ct_freqs_plot <- function(obj = combrna){
  ct_cols = c("red", "lightgrey", "blue")
  
  ct_list_psp = list()
  meta_psp <- obj@meta.data[which(obj@meta.data$NPdiagnosis == 'PSP'),]
  lev <- levels(as.factor(meta_psp$Pub_ID))
  
  for(i in levels(as.factor(meta_psp$Pub_ID))){
    ct_prop <- subset(meta_psp, Pub_ID == i) %>% as.data.frame() %>% 
      dplyr::count(.$predicted.subclass) %>%     
      mutate(prop = prop.table(n))
    colnames(ct_prop)[1] = 'ct'
    ct_prop$ct <- as.character(ct_prop$ct)
    ct_prop <- ct_prop[order(ct_prop$ct),]
    ct_prop$sample = i
    ct_prop$group = 'PSP'
    #ct_list_psp[[paste0(i,'_ct_prop')]] = ct_prop
    
    bp <- ggplot(ct_prop, aes(x='', y=prop, fill=ct))+
      geom_bar(width = 2, color = 'black', stat = "identity") + labs(fill = '', x = i, y = '') + #scale_fill_manual(values = ct_cols)+ 
      theme_minimal(base_family="Arial", base_size = 14) +  theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 30))
    ct_list_psp[[i]] <- bp
  }
  ct_list_psp[['psp']] <- ggpubr::ggarrange(plotlist = ct_list_psp[lev],
                                            ncol = length(ct_list_psp), common.legend = T, legend = 'left')
  ct_list_psp <<- ct_list_psp
  
  ct_list_cbd = list()
  meta_cbd <- obj@meta.data[which(obj@meta.data$NPdiagnosis == 'CBD'),]
  lev <- levels(as.factor(meta_cbd$Pub_ID))
  
  for(i in levels(as.factor(meta_cbd$Pub_ID))){
    ct_prop <- subset(meta_cbd, Pub_ID == i) %>% as.data.frame() %>% 
      dplyr::count(.$predicted.subclass) %>%     
      mutate(prop = prop.table(n))
    colnames(ct_prop)[1] = 'ct'
    ct_prop$ct <- as.character(ct_prop$ct)
    ct_prop <- ct_prop[order(ct_prop$ct),]
    ct_prop$sample = i
    ct_prop$group = 'CBD'
    #ct_list_cbd[[paste0(i,'_ct_prop')]] = ct_prop
    
    bp <- ggplot(ct_prop, aes(x='', y=prop, fill=ct))+
      geom_bar(width = 2, color = 'black',stat = "identity") + labs(fill = '', x = i, y = '') + #scale_fill_manual(values = ct_cols)+ 
      theme_minimal(base_family="Arial", base_size = 14) +  theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 30))
    ct_list_cbd[[i]] <- bp
  }
  
  
  ct_list_cbd[['cbd']] <- ggpubr::ggarrange(plotlist = ct_list_cbd[lev],
                                            ncol = length(ct_list_cbd), common.legend = T, legend = 'none')
  ct_list_cbd <<- ct_list_cbd
  
  
  ct_list_ctrl = list()
  meta_ctrl <- subset(obj@meta.data, NPdiagnosis == 'Ctrl')
  lev <- levels(as.factor(meta_ctrl$Pub_ID))
  
  for(i in levels(as.factor(meta_ctrl$Pub_ID))){
    ct_prop <- subset(meta_ctrl, Pub_ID == i) %>% as.data.frame() %>%
      dplyr::count(.$predicted.subclass) %>%
      mutate(prop = prop.table(n))
    colnames(ct_prop)[1] = 'ct'
    ct_prop$ct <- as.character(ct_prop$ct)
    ct_prop <- ct_prop[order(ct_prop$ct),]
    ct_prop$sample = i
    ct_prop$group = 'Ctrl'
    c#t_list_ctrl[[paste0(i,'_ct_prop')]] = ct_prop
    
    bp <- ggplot(ct_prop, aes(x='', y=prop, fill=ct))+
      geom_bar(width = 2, color = 'black',stat = "identity") + labs(fill = '', x = i, y = '') + #scale_fill_manual(values = ct_cols)+ 
      theme_minimal(base_family="Arial", base_size = 14) +  theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 30))
    ct_list_ctrl[[i]] <- bp
  }
  
  
  ct_list_ctrl[['ctrl']] <- ggpubr::ggarrange(plotlist = ct_list_ctrl[lev], 
                                              ncol = length(ct_list_ctrl), common.legend = T, legend = 'none')
  ct_list_ctrl <<- ct_list_ctrl
  
  ## CT-colors:
  ct_final <- ggpubr::ggarrange(ct_list_psp[['psp']], ct_list_cbd[['cbd']], ct_list_ctrl[['ctrl']], labels = c('PSP', 'CBD', 'Ctrl'),
                                ncol = 3, common.legend = T, legend = 'none',  widths = c(1.3, 1.3, 1), heights = c(1,1,0.7))   
  
  return(ct_final)
}



###
ComputeModuleEigengene <- function (seurat_obj, cur_mod, modules, group.by.vars = NULL, 
                                    verbose = TRUE, vars.to.regress = NULL, scale.model.use = "linear", slot = "counts",
                                    pc_dim = 1, assay = NULL, wgcna_name = NULL, ...) 
{
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  if (is.null(assay)) {
    assay <- DefaultAssay(seurat_obj)
  }
  if (dim(seurat_obj@assays[[assay]]$data)[1] == 0) {
    stop(paste0("Normalized data slot not found in selected assay ", 
                assay))
  }
  cur_genes <- modules %>% subset(module == cur_mod) %>%  .$gene_name %>% 
    as.character()
  
  X_dat <- GetAssayData(seurat_obj, slot = slot, assay = assay)[cur_genes, 
  ]
  if (dim(seurat_obj@assays[[assay]]$counts)[1] == 0) {
    X <- X_dat
  }
  else {
    X <- GetAssayData(seurat_obj, slot = slot, assay = assay)[cur_genes, 
    ]
  }
  cur_seurat <- CreateSeuratObject(as.sparse(X), assay = assay, meta.data = seurat_obj@meta.data)
  cur_seurat <- SetAssayData(cur_seurat, slot = "data", new.data = X_dat, 
                             assay = assay)
  cur_seurat <- ScaleData(cur_seurat, features = rownames(cur_seurat), 
                          model.use = 'linear', vars.to.regress = NULL)
  
  cur_expr <- GetAssayData(cur_seurat, slot = "data")
  expr <- Matrix::t(cur_expr)
  averExpr <- Matrix::rowSums(expr)/ncol(expr)
  cur_pca <- Seurat::RunPCA(cur_seurat, features = cur_genes, 
                            reduction.key = paste0("pca", cur_mod), verbose = T)@reductions$pca
  pc <- cur_pca@cell.embeddings[, pc_dim]
  pc_loadings <- cur_pca@feature.loadings[, pc_dim]
  pca_cor <- cor(averExpr, pc)
 
  if (pca_cor < 0) {
    cur_pca@cell.embeddings[, pc_dim] <- -pc
    pc_loadings <- -pc_loadings
  }
  seurat_obj@reductions$ME <- Seurat::CreateDimReducObject(embeddings = cur_pca@cell.embeddings, 
                                                           assay = assay)
  seurat_obj <- SetMELoadings(seurat_obj, loadings = pc_loadings, 
                              harmonized = FALSE, wgcna_name = wgcna_name)
  seurat_obj
}



ModuleUMAP_custom <- function (seurat_obj, sample_edges = TRUE, edge_prop = 0.2, 
                               label_hubs = 5, edge.alpha = 0.25, vertex.label.cex = 0.5, 
                               label_genes = NULL, return_graph = FALSE, keep_grey_edges = TRUE, 
                               wgcna_name = NULL, ...) 
{
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  TOM <- GetTOM(seurat_obj, wgcna_name)
  modules <- GetModules(seurat_obj, wgcna_name)
  umap_df <- GetModuleUMAP(seurat_obj, wgcna_name)
  ####
  dict_tf <- atac@assays$chromvar@meta.features %>% mutate(gene = rownames(atac@assays$chromvar@meta.features))
  umap_df <- left_join(umap_df,
                       dict_tf,
                       by = c('gene')) %>% column_to_rownames('gene')
  
  umap_df$gene <- rownames(umap_df)
  ###
  mods <- levels(umap_df$module)
  mods <- mods[mods != "grey"]
  subset_TOM <- TOM[umap_df$gene, umap_df$gene[umap_df$hub == 
                                                 "hub"]]
  hub_list <- lapply(mods, function(cur_mod) {
    cur <- subset(modules, module == cur_mod)
    cur[, c("gene_name", paste0("kME_", cur_mod))] %>% top_n(label_hubs) %>% 
      .$gene_name
  })
  names(hub_list) <- mods
  hub_labels <- as.character(unlist(hub_list))
  print("hub labels")
  print(hub_labels)
  print(label_genes)
  if (is.null(label_genes)) {
    label_genes <- hub_labels
  }
  else {
    if (!any(label_genes %in% umap_df$gene)) {
      stop("Some genes in label_genes not found in the UMAP.")
    }
    label_genes <- unique(c(label_genes, hub_labels))
  }
  print(label_genes)
  selected_modules <- modules[umap_df$gene, ]
  selected_modules <- cbind(selected_modules, umap_df[, c("UMAP1", 
                                                          "UMAP2", "hub", "kME")])
  selected_modules$label <- ifelse(selected_modules$gene_name %in% 
                                     label_genes, selected_modules$gene_name, "")
  selected_modules$fontcolor <- ifelse(selected_modules$color == 
                                         "black", "gray50", "black")
  selected_modules$framecolor <- ifelse(selected_modules$gene_name %in% 
                                          label_genes, "black", selected_modules$color)
  edge_df <- subset_TOM %>% reshape2::melt()
  print(dim(edge_df))
  edge_df$color <- future.apply::future_sapply(1:nrow(edge_df), 
                                               function(i) {
                                                 gene1 = as.character(edge_df[i, "Var1"])
                                                 gene2 = as.character(edge_df[i, "Var2"])
                                                 col1 <- selected_modules[selected_modules$gene_name == 
                                                                            gene1, "color"]
                                                 col2 <- selected_modules[selected_modules$gene_name == 
                                                                            gene2, "color"]
                                                 if (col1 == col2) {
                                                   col = col1
                                                 }
                                                 else {
                                                   col = "grey90"
                                                 }
                                                 col
                                               })
  if (!keep_grey_edges) {
    edge_df <- edge_df %>% subset(color != "grey90")
  }
  groups <- unique(edge_df$color)
  if (sample_edges) {
    temp <- do.call(rbind, lapply(groups, function(cur_group) {
      cur_df <- edge_df %>% subset(color == cur_group)
      n_edges <- nrow(cur_df)
      cur_sample <- sample(1:n_edges, round(n_edges * 
                                              edge_prop))
      cur_df[cur_sample, ]
    }))
  }
  else {
    temp <- do.call(rbind, lapply(groups, function(cur_group) {
      cur_df <- edge_df %>% subset(color == cur_group)
      n_edges <- nrow(cur_df)
      cur_df %>% dplyr::top_n(round(n_edges * edge_prop), 
                              wt = value)
    }))
  }
  edge_df <- temp
  print(dim(edge_df))
  edge_df <- edge_df %>% group_by(color) %>% mutate(value = scale01(value))
  edge_df <- edge_df %>% arrange(value)
  edge_df <- rbind(subset(edge_df, color == "grey90"), subset(edge_df, 
                                                              color != "grey90"))
  edge_df$color_alpha <- ifelse(edge_df$color == "grey90", 
                                alpha(edge_df$color, alpha = edge_df$value/2), alpha(edge_df$color, 
                                                                                     alpha = edge_df$value))
  selected_modules <- rbind(subset(selected_modules, hub == 
                                     "other"), subset(selected_modules, hub != "other"))
  selected_modules <- rbind(subset(selected_modules, label == 
                                     ""), subset(selected_modules, label != ""))
  ####
  colnames(dict_tf)[2] <- 'Var1'
  edge_df<- left_join(edge_df,
                      dict_tf,
                      by = c('Var1'))
  colnames(dict_tf)[2] <- 'Var2'
  edge_df <- left_join(edge_df,
                       dict_tf,
                       by = c('Var2'))
  edge_df$Var1 <- as.factor(paste0(edge_df$Var2,':',edge_df$tf_name.x) )
  edge_df$Var2 <- as.factor(paste0(edge_df$Var2,':',edge_df$tf_name.y) )
  edge_df <- edge_df[,1:5]
  
  selected_modules$Var2 = rownames(selected_modules)
  selected_modules <- left_join(selected_modules,
                                dict_tf,
                                by = c('Var2')) 
  selected_modules$tf_name <- paste0(selected_modules$Var2,':',selected_modules$tf_name)
  selected_modules <- column_to_rownames(selected_modules,'tf_name')
  selected_modules$gene_name <- rownames(selected_modules) 
  selected_modules$label <- ifelse(selected_modules$gene_name %in% 
                                     label_genes, selected_modules$gene_name, "")
  
  x <- subset(edge_df, edge_df$Var1 %in% selected_modules$gene_name)
  ####
  g <- igraph::graph_from_data_frame(x, directed = FALSE, 
                                     vertices = selected_modules)
  print("making net")
  print(head(edge_df))
  print(head(selected_modules))
  if (return_graph) {
    return(g)
  }
  plot(g, layout = as.matrix(selected_modules[, c("UMAP1", 
                                                  "UMAP2")]), edge.color = adjustcolor(E(g)$color_alpha, 
                                                                                       alpha.f = edge.alpha), vertex.size = V(g)$kME * 3, edge.curved = 0, 
       edge.width = 0.5, vertex.color = V(g)$color, vertex.label = V(g)$label, 
       vertex.label.dist = 1.1, vertex.label.degree = -pi/4, 
       vertex.label.family = "Helvetica", vertex.label.font = 3, 
       vertex.label.color = V(g)$fontcolor, vertex.label.cex = 0, 
       vertex.frame.color = V(g)$framecolor, margin = 0)
}


### ModEigengene
ModEigengene <- function(seurat_obj = seurat_obj, 
                         vars.to.regress = c('Age', "pmi"),
                         exclude_grey = T,  
                         assay = DefaultAssay(seurat_obj), 
                         group.by.vars="Pub_ID", slot="counts")
{
  wgcna_name <- seurat_obj@misc$active_wgcna
  harmonized = !is.null(group.by.vars)
  me_list <- list()
  harmonized_me_list <- list()
  seurat_obj <- SetMELoadings(seurat_obj, loadings = c(""),
                              harmonized = TRUE, wgcna_name = wgcna_name)
  modules <- GetModules(seurat_obj, wgcna_name)
  projected <- FALSE
  mods <- levels(modules$module)
  mods_loop <- mods
  pc_dim = 1
  for (cur_mod in mods_loop) {
    print(cur_mod)
    seurat_obj <- ComputeModuleEigengene(seurat_obj = seurat_obj,
                                         cur_mod = cur_mod, modules = modules, group.by.vars = group.by.vars,
                                         vars.to.regress = NULL, scale.model.use = NULL, slot = slot,
                                         verbose = verbose, pc_dim = pc_dim, assay = assay,
                                         wgcna_name = wgcna_name)
    cur_me <- seurat_obj@reductions$ME@cell.embeddings[,
                                                       pc_dim]
    me_list[[cur_mod]] <- cur_me
  }
  me_df <- do.call(cbind, me_list)
  me_df <- WGCNA::orderMEs(me_df)
  seurat_obj <- SetMEs(seurat_obj, me_df, harmonized = FALSE, wgcna_name)
  MEs <- GetMEs(seurat_obj, harmonized, wgcna_name)
  modules$module <- factor(as.character(modules$module), levels = mods)
  seurat_obj <- SetModules(seurat_obj, modules, wgcna_name)
  seurat_obj@reductions$ME <- NULL
  seurat_obj@reductions$ME_harmony <- NULL
  return(seurat_obj)
}


### datExpr
custom_SetDatExpr <- function(
    seurat_obj = seurat_obj,
    group_name = "Astro",# the name of the group of interest in the group.by column
    group.by='predicted.subclass', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
    assay = 'RNA', # using RNA assay
    slot = 'scale.data' # using normalized data
)
{
  params <- GetWGCNAParams(seurat_obj, 'tau')
  genes_use <- GetWGCNAGenes(seurat_obj, 'tau')
  m_obj <- GetMetacellObject(seurat_obj, 'tau')
  seurat_meta <- seurat_obj@meta.data
  assay <- DefaultAssay(seurat_obj)
  seurat_meta <- seurat_meta %>% subset(get(group.by) %in% group_name)
  cells <- rownames(seurat_meta)
  datExpr <-BiocGenerics::as.data.frame(Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot)[genes_use, cells])
  datExpr <- BiocGenerics::as.data.frame(t(datExpr))
  gene_list <- genes_use[which(!is.na(genes_use[WGCNA::goodGenes(datExpr)]))]
  datExpr <- datExpr[,gene_list]
  seurat_obj <- SetWGCNAGenes(seurat_obj, gene_list, wgcna_name='tau')
  seurat_obj@misc[['tau']]$datExpr <- as.data.frame(datExpr)
  rm(group_name,group.by, assay, slot,params,genes_use,m_obj,seurat_meta,cells,datExpr,gene_list)
  gc()
  return(seurat_obj)
}

### for v5 Seurat
custom_SetDatExpr_v5 <- function(
    seurat_obj = seurat_obj_global_ds,
    group_name = "Astro",# the name of the group of interest in the group.by column
    group.by='celltype', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
    assay = 'SCT', # using RNA assay
    layer = 'data', # using normalized data
    wgcna_name = "psp_hdwgcna",use_metacells=T,
    multi.group.by = NULL, multi_group_name = NULL, return_seurat = TRUE){
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  params <- GetWGCNAParams(seurat_obj, wgcna_name)
  genes_use <- GetWGCNAGenes(seurat_obj, wgcna_name)
  m_obj <- GetMetacellObject(seurat_obj, wgcna_name)
  if (use_metacells & !is.null(m_obj)) {
    s_obj <- m_obj
  }  else {
    if (is.null(m_obj)) {
      warning("Metacell Seurat object not found. Using full Seurat object instead.")
    }
    s_obj <- seurat_obj
  }
  seurat_meta <- s_obj@meta.data
  if (is.null(assay)) {
    assay <- DefaultAssay(s_obj)
    warning(paste0("assay not specified, trying to use assay ",
                   assay))
  }
  if (!(assay %in% names(s_obj@assays))) {
    stop(
      "Assay not found. Check names(seurat_obj@assays) or names(GetMetacellObject(seurat_obj)@assays)"
    )
  }
  if (!is.null(group.by)) {
    if (!(group.by %in% colnames(s_obj@meta.data))) {
      m_cell_message <- ""
      if (use_metacells) {
        m_cell_message <- "metacell"
      }
      stop(
        paste0(
          group.by,
          " not found in the meta data of the ",
          m_cell_message,
          " Seurat object"
        )
      )
    }
    if (!all(group_name %in% s_obj@meta.data[[group.by]])) {
      groups_not_found <- group_name[!(group_name %in%
                                         s_obj@meta.data[[group.by]])]
      stop(paste0(
        "Some groups in group_name are not found in the seurat_obj: ",
        paste(groups_not_found, collapse = ", ")
      ))
    }
  }
  if (!is.null(group.by)) {
    seurat_meta <- seurat_meta %>% subset(get(group.by) %in%
                                            group_name)
  }
  if (!is.null(multi.group.by)) {
    seurat_meta <- seurat_meta %>% subset(get(multi.group.by) %in%
                                            multi_group_name)
  }
  cells <- rownames(seurat_meta)
  s_obj <- JoinLayers(s_obj, assay = assay, layers = layer)
  datExpr <-as.data.frame(Seurat::GetAssayData(s_obj, assay = assay,
                                       layer = layer)[genes_use, cells])
  datExpr <- as.data.frame(t(datExpr))
  if (return_seurat) {
    gene_list <- genes_use[WGCNA::goodGenes(datExpr)]
    datExpr <- datExpr[, gene_list]
    seurat_obj <- SetWGCNAGenes(seurat_obj, gene_list, wgcna_name)
    seurat_obj@misc[[wgcna_name]]$datExpr <- datExpr
    out <- seurat_obj
  }  else {
    out <- datExpr
  }
  out
}

###
custom_plot_GO <- function (seurat_obj, database, mods = "all", n_terms = 3, break_ties = TRUE, 
                            logscale = TRUE, wgcna_name = NULL, ...) 
{
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  modules <- GetModules(seurat_obj, wgcna_name)
  if (mods == "all") {
    mods <- levels(modules$module)
    mods <- mods[mods != "grey"]
  }
  enrichr_df <- GetEnrichrTable(seurat_obj, wgcna_name)
  mod_colors <- dplyr::select(modules, c(module, color)) %>% distinct
  enrichr_df$color <- mod_colors[match(enrichr_df$module, 
                                       mod_colors$module), "color"]
  wrapText <- function(x, len) {
    sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), 
           USE.NAMES = FALSE)
  }
  plot_df <- enrichr_df %>% subset(db == database & module %in% 
                                     mods) %>% group_by(module) %>% top_n(n_terms, wt = Combined.Score)
  if (break_ties) {
    plot_df <- do.call(rbind, lapply(plot_df %>% group_by(module) %>% 
                                       group_split, function(x) {
                                         x[sample(n_terms), ]
                                       }))
  }
  plot_df$Term <- wrapText(plot_df$Term, 45)
  plot_df$module <- factor(as.character(plot_df$module), levels = levels(modules$module))
  plot_df <- arrange(plot_df, module)
  plot_df$Term <- factor(as.character(plot_df$Term), levels = unique(as.character(plot_df$Term)))
  if (logscale) {
    plot_df$Combined.Score <- log(plot_df$Combined.Score)
    lab <- "Enrichment\nlog(combined score)"
    x <- 0.2
  }
  else {
    lab <- "Enrichment\n(combined score)"
    x <- 5
  }
  p <- plot_df %>% ggplot(aes(x = module, y = Term)) + geom_point(aes(size = Combined.Score), 
                                                                  color = plot_df$color) + RotatedAxis() + ylab("") + 
    xlab("") + labs(size = lab) + scale_y_discrete(limits = rev) + 
    ggtitle(database) + theme(plot.title = element_text(hjust = 0.5), 
                              axis.line.x = element_blank(), axis.line.y = element_blank(), 
                              panel.border = element_rect(colour = "black", fill = NA, 
                                                          size = 1))
  p
}

###
custom_plotKMEs <- function (seurat_obj, n_hubs = 10, text_size = 2, ncol = 5, 
                             plot_widths = c(3, 2), wgcna_name = NULL) 
{
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  modules <- GetModules(seurat_obj, wgcna_name) %>% subset(module != 
                                                             "grey")
  mods <- levels(modules$module)
  mods <- mods[mods != "grey"]
  mod_colors <- modules %>% subset(module %in% mods) %>% dplyr::select(c(module, 
                                                                         color)) %>% distinct
  hub_df <- do.call(rbind, lapply(mods, function(cur_mod) {
    print(cur_mod)
    cur <- subset(modules, module == cur_mod)
    cur <- cur[, c("gene_name", "module", paste0("kME_", 
                                                 cur_mod))]
    names(cur)[3] <- "kME"
    cur <- dplyr::arrange(cur, kME)
    top_genes <- cur %>% dplyr::top_n(n_hubs, wt = kME) %>% 
      .$gene_name
    cur$lab <- ifelse(cur$gene_name %in% top_genes, cur$gene_name, 
                      "")
    cur
  }))
  head(hub_df)
  plot_list <- lapply(mods, function(x) {
    print(x)
    cur_color <- subset(mod_colors, module == x) %>% .$color
    cur_df <- subset(hub_df, module == x)
    top_genes <- cur_df %>% dplyr::top_n(n_hubs, wt = kME) %>% 
      .$gene_name
    p <- cur_df %>% ggplot(aes(x = reorder(gene_name, kME), 
                               y = kME)) + geom_bar(stat = "identity", width = 1, 
                                                    color = cur_color, fill = cur_color) + ggtitle(x) + 
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
            plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), 
            axis.line.x = element_blank())
    p_anno <- ggplot() + annotate("label", x = 0, y = 0, 
                                  label = paste0(top_genes, collapse = "\n"), size = text_size, 
                                  fontface = "italic", label.size = 0) + theme_void()
    patch <- p + p_anno + plot_layout(widths = plot_widths)
    patch
  })
  wrap_plots(plot_list, ncol = ncol)
}

###
DME_loop <- function(seurat_obj = cell_seu, 
                     clusters, 
                     group_var = "NPdiagnosis", 
                     group_val = "PSP", 
                     wgcna_name = "psp_hdwgcna"
                     ) {
  ls <- list()
  # Loop through the clusters
  for(cur_cluster in clusters){
    
    # Identify barcodes for group1 and group2 in each cluster
    group1 <- subset(seurat_obj@meta.data, subcluster == cur_cluster & NPdiagnosis == group_val) %>% rownames
    group2 <- subset(seurat_obj@meta.data, subcluster == cur_cluster & NPdiagnosis == "Ctrl") %>% rownames
    
    if(length(group1) >5 && length(group2) >5){
     
        # Run the DME test
        cur_DMEs <- FindDMEs(
          seurat_obj,
          barcodes1 = group1,
          barcodes2 = group2,
          test.use='wilcox',
          pseudocount.use=0.01, # We can also change the pseudocount with this param
          wgcna_name = wgcna_name
        )
        
        # Add the cluster info to the table
        cur_DMEs$cluster <- cur_cluster
        
        # Append the table
        ls[[cur_cluster]] <- cur_DMEs
      }
  }
  DMEs <- do.call("rbind", ls)
  
  # Return the final DMEs data frame
  return(DMEs)
}


plot_heatmap_DME <- function(seurat_obj = cell_seu, 
                             DMEs, 
                             GO = NULL,
                             n_terms = 3,
                             w = c(1,0.9)){
  message("Depracated. Use plot_heatmap_mods with modus = 'DME' instead")
  modules <- GetModules(seurat_obj)
  mods <- levels(modules$module); mods <- mods[mods != 'grey'] %>% gsub("_","-",.)
  plot_df <- DMEs
  plot_df$module <- factor(plot_df$module, levels=rev(mods))
  maxval <- max(plot_df$avg_log2FC); minval <- min(plot_df$avg_log2FC)
  plot_df$avg_log2FC <- ifelse(plot_df$avg_log2FC > maxval, maxval, plot_df$avg_log2FC)
  plot_df$avg_log2FC <- ifelse(plot_df$avg_log2FC < minval, minval, plot_df$avg_log2FC)
  plot_df$Significance <- gtools::stars.pval(plot_df$p_val_adj)
  plot_df$textcolor <- ifelse(plot_df$avg_log2FC > 0.2, 'black', 'white')
  p <- plot_df %>% 
    ggplot(aes(x=str_sort(cluster), y= module, fill=avg_log2FC)) +
    geom_tile() + 
    geom_text(label=plot_df$Significance, color=plot_df$textcolor) + 
    scale_fill_gradient2(low='blue3', mid='white', high='red3') +
    RotatedAxis() +
    theme(
      panel.border = element_rect(fill=NA, color='black', size=1),
      axis.line.x = element_blank(),
      axis.line.y = element_blank(), 
      legend.position = "left",
      plot.margin=margin(0,0,0,0)
    ) + xlab('') + ylab('') 
  
  if(!is.null(GO)){
    GO <- GO %>% group_by(module) %>% dplyr::filter(p.adjust < .1) %>%
      slice_min(order_by = pvalue, n = n_terms) %>% # adjust top vars
      ungroup()
    
    GO$module <- factor(GO$module, levels = gtools::mixedsort(levels(GO$module)))

    p2 <- ggplot(GO, aes(x = Count, y = Description, fill = -log(p.adjust))) + 
      geom_col(color = "grey3") + 
      scale_fill_gradient2(low='white', high='red3')  + 
      ggplot2::facet_grid(module ~., space = "free_y", scales="free") + 
      scale_y_discrete(position = "right") +
      scale_x_continuous() +
      theme(
        panel.border = element_rect(fill=NA, color='white', size=1),
        strip.text.y = element_text(angle = 0)
      ) + xlab('Gene Count') + ylab('')
    
    return(
      cowplot::plot_grid(plotlist = list(p, p2), ncol = 2, 
                              align = "hv", axis = "b", 
                         rel_widths = w)
           )
    
  }else{
    
    return(p)
  }
  
}

###
plot_heatmap_mods <- function(seurat_obj = cell_seu, 
                              DMEs = NULL, 
                              GO = NULL, 
                              n_terms = 3, 
                              color_pal = colorRamp2(c(-3, 0, 3), c("blue", "white", "red3")),
                              w = 5,
                              h = 5,
                              modus = "DME",
                              annot = T,
                              cluster_rows = T,
                              cluster_columns = T,
                              sccomp_res = NULL,
                              annot_style = "anno_link") {
  library(dplyr)
  library(purrr)
  library(data.table)
  library(gtools)
  library(ComplexHeatmap)
  library(circlize)
  library(matrixTests)
  library(cowplot)
  # Common setup
  modules <- GetModules(seurat_obj)
  mods <- levels(modules$module); mods <- mods[mods != 'grey'] %>% gsub("-","_",.)
  clusters <- levels(as.factor(as.character(seurat_obj$subcluster))) %>% gsub("-","_",.)

    # DME modus / diff data prep required
  if (modus == "DME") {
    plot_df <- DMEs
    plot_df$module <- factor(gsub("-","_",plot_df$module), levels=rev(mods))
    maxval <- max(plot_df$avg_log2FC, na.rm = TRUE); minval <- min(plot_df$avg_log2FC, na.rm = TRUE)
    plot_df$avg_log2FC <- pmin(pmax(plot_df$avg_log2FC, minval), maxval)
    plot_df$Significance <- gtools::stars.pval(plot_df$p_val_adj)
    plot_df$textcolor <- ifelse(plot_df$avg_log2FC > 0.2, 'black', 'white')
    plot_df$cluster <- factor(gsub("-","_",plot_df$cluster), levels=clusters) #changed 2024-10-13

    
    # Prepare data for ComplexHeatmap
    heatmap_data <- matrix(plot_df$avg_log2FC, nrow = length(unique(plot_df$module)), 
                           ncol = length(unique(plot_df$cluster))) 
    rownames(heatmap_data) <- unique(plot_df$module)
    colnames(heatmap_data) <- unique(plot_df$cluster)
    heatmap_data <- heatmap_data[gtools::mixedsort(rownames(heatmap_data)),] 
    
    if(length(clusters)>1){
      heatmap_data <- sparrow::scale_rows(heatmap_data)
    }

    # Annotations for significance
    mat_pval <- matrix(plot_df$p_val_adj, nrow = length(unique(plot_df$module)), 
                       ncol = length(unique(plot_df$cluster)))
    rownames(mat_pval) <- unique(plot_df$module)
    colnames(mat_pval) <- unique(plot_df$cluster)
    
    ego <- GO %>% group_by(module) %>% dplyr::filter(p.adjust < .1) %>%
      slice_min(order_by = pvalue, n = n_terms) %>% # adjust top vars
      ungroup()
    
    ego$module <- factor(ego$module, levels = gtools::mixedsort(levels(ego$module)))
    
    # scl modus / diff data prep required
  } else if (modus == "scl") {
    
    ME_cluster <- sapply(unique(clusters), function(c) {
      cells <- base::rownames(subset(seurat_obj@meta.data, subcluster == c))
      seurat_obj@misc$psp_hdwgcna$MEs[cells, str_sort(mods)]
    }) %>% .[, clusters]
    
    heatmap_data <- sapply(unique(clusters), function(c) {
      cells <- rownames(subset(seurat_obj@meta.data, subcluster == c))
      colSums(seurat_obj@misc$psp_hdwgcna$MEs[cells, str_sort(mods)])
    }) %>% sparrow::scale_rows()
    heatmap_data <- heatmap_data[gtools::mixedsort(rownames(heatmap_data)),] 
    
    mat_pval <- matrix(0, nrow = length(rownames(ME_cluster)), ncol = length(colnames(ME_cluster)))
    rownames(mat_pval) <- rownames(ME_cluster)
    colnames(mat_pval) <- colnames(ME_cluster)
    
    # Iterate over each row in ME_cluster
    for(m in rownames(ME_cluster)){
      for (i in colnames(ME_cluster)) {
        j = colnames(ME_cluster)[!colnames(ME_cluster) == i]
        # Perform the Wilcoxon test comparing row i against row j
        test_result <- wilcox.test(
          x = ME_cluster[m,][[i]],
          y = unlist(ME_cluster[m,][j])
        )
        
        # Store the p-value
        mat_pval[m, i] <- test_result$p.value
      }
      
      mat_pval[m, ] <- p.adjust(mat_pval[m, ], method = "bonf")
    }

    ego <-  seurat_obj@misc$psp_hdwgcna$enrichr_table %>%
      group_by(module) %>%
      dplyr::filter(p.adjust < 0.1) %>% 
      slice_min(order_by = pvalue, n = n_terms, na_rm = T) %>%
      ungroup() 
   
  }
  
  hubs <- GetHubGenes(seurat_obj, n_hubs = 5) %>% mutate(module = factor(module))
  hubs <- split(hubs$gene_name, hubs$module)%>% 
    sapply(., function(x) paste(x, collapse = "/"),
           simplify = FALSE)
  text <- split(ego$Description, ego$module) %>% 
    sapply(., function(x) paste(x, collapse = " | "),
           simplify = FALSE) 
  
  # Add list element, for layout in complexHeatmap
  if(!is.null(base::setdiff(unique(mods),names(text)))){
    missing <- base::setdiff(unique(mods),names(text))
    for(i in missing){
      text[[i]] = "*No GO enrichment*"  
    }
  }
  
  text <- mapply(function(x, y) paste(x, y, sep = ": "), hubs[names(hubs)], text[names(hubs)], SIMPLIFY = F)
 
  cell_fun <- function(j, i, x, y, width, height, fill) {
    if (mat_pval[i, j] < 0.001) {
      grid.text("***", x, y, just = "center")
    } else if (mat_pval[i, j] < 0.01) {
      grid.text("**", x, y, just = "center")
    } else if (mat_pval[i, j] < 0.05) {
      grid.text("*", x, y, just = "center")
    }
  }
  
  module_df <- data.frame(
    module = factor(unique(seurat_obj@misc$psp_hdwgcna$module_umap$module)),
    color = unique(seurat_obj@misc$psp_hdwgcna$module_umap$color)
  ) %>% deframe()
  
  # Create the left & top annotation
  left_ann <- rowAnnotation(
    Module = names(module_df), 
    col = list(Module = module_df
    ),
    gp = gpar(col = "grey42")
  )
  
  top_ann <- NULL
  
  if(!is.null(sccomp_res)){
    continuous_values <- setNames(sccomp_res$effect, sccomp_res$subcluster)
    if(length(clusters)>1){
      continuous_values <- continuous_values[colnames(heatmap_data)]
    }
    # Create the top annotation
    top_ann <- HeatmapAnnotation(
      subcluster_abundance_psp = continuous_values,
      col = list(subcluster_abundance_psp = colorRamp2(c(min(continuous_values), 0, max(continuous_values)),
                                                    c("black", "white", "red3"))),
      gp = gpar(col = "grey42"),
      subcluster_abundance_psp_bar = anno_barplot(
        continuous_values, gp = gpar(fill = "grey")
      ), annotation_name_side = "right"
    )
  }
  
  if(annot){
    ann <- rowAnnotation(
      textbox = anno_textbox(
        mods,
        text,by = annot_style,
        gp = gpar(
          col = "black",
          fontsize = 8),
        background_gp = gpar(fill = "grey92", 
                             col = "grey4"),
        word_wrap = TRUE, 
        add_new_line = TRUE,
      ))
    
    # Create the heatmap with the left annotation
    if(!is.null(nrow(heatmap_data))){
    heat <- Heatmap(
      heatmap_data,
      name = "data",
      col = color_pal,
      border = TRUE,
      cluster_rows = cluster_rows,
      cluster_columns = cluster_columns,
      column_names_gp = gpar(fontsize = 10),
      row_names_gp = gpar(fontsize = 10),
      cell_fun = cell_fun,
      rect_gp = gpar(col = "grey42"),
      width = ncol(heatmap_data) * unit(w, "mm"),
      height = nrow(heatmap_data) * unit(h, "mm"),
      column_title = modus,
      heatmap_legend_param = list(title = "ME\nDifferential ME"),
      right_annotation = ann,
      left_annotation = left_ann,
      top_annotation = top_ann)
    heat
    }else{
      heat <- Heatmap(
        heatmap_data,
        name = "data",
        col = color_pal,
        border = TRUE,
        cluster_rows = cluster_rows,
        cluster_columns = cluster_columns,
        column_names_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 10),
        cell_fun = cell_fun,
        rect_gp = gpar(col = "grey42"),
        width = unit(w, "mm"),
        height = length(heatmap_data) * unit(h, "mm"),
        column_title = modus,
        heatmap_legend_param = list(title = "ME\nDifferential ME"),
        right_annotation = ann,
        left_annotation = left_ann)
      heat
    }
  }else{
    heat <- Heatmap(
      heatmap_data,
      name = "data",
      col = color_pal,
      border = TRUE,
      cluster_rows = cluster_rows,
      cluster_columns = cluster_columns,
      column_names_gp = gpar(fontsize = 10),
      row_names_gp = gpar(fontsize = 10),
      right_annotation = NULL,
      cell_fun = cell_fun,
      rect_gp =  gpar(col = "grey42"),
      width = ncol(heatmap_data)*unit(w, "mm"),
      height = nrow(heatmap_data)*unit(h, "mm"),
      column_title = modus,
      heatmap_legend_param = list(title = "ME\nDifferential ME"),
      left_annotation = left_ann,
      top_annotation = top_ann)
  }
  return(suppressMessages(heat))
}

###
compute_go_hdwgcna <- function(seurat_obj, 
                               dbs = c("BP"), 
                               me_cutoff = 0.05, 
                               wgcna_name = NULL,  
                               pAdjustMethod = "BH",
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.1,
                               max_genes = 5000) 
{
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(DOSE)
  library(ggupset)
  library(pathview)
  
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  
  modules <- GetModules(seurat_obj, wgcna_name)
  mods <- levels(modules$module)
  mods <- mods[mods != "grey"]
  combined_output <- data.frame()
  
  for (i in 1:length(mods)) {
    cur_mod <- mods[i]
    if (max_genes != Inf) {
      cur_info <- subset(modules, module == cur_mod)
      nam <- paste0("kME_", cur_mod)
      cur_info <- cur_info[, c("gene_name", nam)]
      cur_genes <- cur_info[which(cur_info[[nam]] > me_cutoff),] %>% .$gene_name %>% 
        as.character
    }  else {
      cur_genes <- subset(modules, module == cur_mod) %>% 
        .$gene_name %>% as.character
    }
    
    # Run GO enrichment analysis
    ego <- enrichGO(
      gene = cur_genes,
      universe =  rownames(modules),
      keyType = "SYMBOL",
      OrgDb = org.Hs.eg.db,
      ont = dbs,
      pAdjustMethod = pAdjustMethod,
      pvalueCutoff = pvalueCutoff,
      qvalueCutoff = qvalueCutoff,
      readable = T
    )
    
    # Simplify the enrichment results to reduce redundancy
    if (nrow(ego@result)> 1) {
      ego_simplified <- clusterProfiler::simplify(ego)
      cluster_summary <- data.frame(ego)
      
      saveRDS(ego, file = paste0("hdwgcna_enrichGO",cur_mod,".Rds"))
      write.csv(
        cluster_summary,
        file = paste0("hdwgcna_enrichGO",cur_mod,".csv"),
        row.names = F
      )
    
    
    if (nrow(cluster_summary) > 1) {
      cluster_summary$module <- cur_mod
      combined_output <- rbind(combined_output, cluster_summary)
    }
  }
  seurat_obj <- SetEnrichrTable(seurat_obj, 
                                combined_output, 
                                wgcna_name)
  }
  return(seurat_obj)
}


####
plot_module_trajectory <-
  function(seurat_obj = cell_seu,
           pseudotime_col = "pseudotime",
           n_bins = 20,
           harmonized = TRUE,
           ncol = 4,
           point_size = 1,
           line_size = 1,
           se = TRUE,
           group_colors = NULL,
           wgcna_name = NULL,
           patch = F)
    
{
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  if (!is.null(group_colors)) {
    if (length(group_colors) != length(pseudotime_col)) {
      stop("group_colors should contain the same number of elements as pseudotime_col")
    }
    group_cp <- group_colors
    names(group_cp) <- pseudotime_col
  }
  for (ps_col in pseudotime_col) {
    seurat_obj <- BinPseudotime(seurat_obj, 
                                pseudotime_col = ps_col, 
                                n_bins = n_bins)
  }
  MEs <- GetMEs(seurat_obj, harmonized = harmonized, wgcna_name = wgcna_name)
  modules <- GetModules(seurat_obj, wgcna_name = wgcna_name)
  mods <- levels(modules$module)
  mods <- mods[mods != "grey"]
  module_colors <- modules %>% dplyr::select(c(module, color)) %>% 
    dplyr::distinct()
  mod_colors <- module_colors$color
  names(mod_colors) <- module_colors$module
  if(!any(mods %in% colnames(seurat_obj@meta.data))){
    plot_df <- cbind(seurat_obj@meta.data, MEs[rownames(seurat_obj@meta.data),])
  }else{
    plot_df <- seurat_obj@meta.data
  }
  # print(head(plot_df))
  avg_list <- list()
  for (ps_col in pseudotime_col) {
    ps_bin_name <- paste0(ps_col, "_bins_", n_bins)
    
    avg_scores <- plot_df %>% 
      dplyr::group_by(get(ps_bin_name)) %>% 
      dplyr::select(all_of(mods)) %>% 
      summarise_all(mean) %>%
      # Calculate relative deviations from first bin
      mutate(across(all_of(mods), ~ . - dplyr::first(.)))  # Key modification
    
    colnames(avg_scores)[1] <- "bin"
    avg_df <- reshape2::melt(avg_scores)
    avg_df$bin <- as.numeric(avg_df$bin)
    avg_df$group <- ps_col
    avg_list[[ps_col]] <- avg_df
  }
  plot_df <- do.call(rbind, avg_list)
  if (length(pseudotime_col) == 1) {
    p <- ggplot(plot_df, aes(x = bin, y = value, color = variable, 
                             fill = variable)) + 
      # geom_point(size = point_size, aes(fill = variable), color = "grey42", shape = 21) + 
      geom_smooth(size = line_size, se = se) + 
      geom_hline(yintercept = 0, 
                 linetype = "dashed", color = "grey22") + 
      scale_color_manual(values = mod_colors) + 
      scale_fill_manual(values = mod_colors) + 
      xlab("Pseudotime [bin]") + 
      ylab("Delta(Module Eigengene)") + 
      theme_test() +
      theme(
        # axis.ticks.x = element_blank(),
        # axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        panel.border = element_rect(
          size = 1,
          fill = NA,
          color = "black"
        )
      ) +
      NoLegend()
  }else {
    p <- ggplot(plot_df, aes(x = as.numeric(bin), y = value, 
                             fill = group)) + 
      geom_smooth(se = se, size = line_size) + 
      geom_hline(yintercept = 0, linetype = "dashed", 
                 color = "grey") + xlab("Pseudotime") + ylab("Module Eigengene") + 
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
            axis.line.x = element_blank(), axis.line.y = element_blank(), 
            panel.border = element_rect(size = 1, fill = NA, 
                                        color = "black")) + labs(color = "Trajectory")
    if (!is.null(group_colors)) {
      p <- p + scale_color_manual(values = group_cp)
    }
  }
  
  if(patch){
    patch <- p + facet_wrap(~variable, ncol = ncol, scales = "free")
    return(patch)
    
  }else{
    
    return(p)
  }
  
  }


###
plot_gene_trajectory <- function(seurat_obj = cell_seu,
                                 genes = NULL,
                                 pseudotime_col = "pseudotime",
                                 n_bins = 20,
                                 ncol = 4,
                                 point_size = 1,
                                 line_size = 1,
                                 se = TRUE,
                                 group_colors = NULL,
                                 patch = FALSE) {
  
  if (is.null(genes)) {
    stop("Please provide a list of genes to plot.")
  }
  
  for (ps_col in pseudotime_col) {
    seurat_obj <- BinPseudotime(seurat_obj, 
                                pseudotime_col = ps_col, 
                                n_bins = n_bins)
  }
  
  # Extract expression data for the specified genes
  gene_expr <- FetchData(seurat_obj, vars = genes)
  
  # Add pseudotime bins to the expression data
  plot_df <- cbind(seurat_obj@meta.data, gene_expr)
  
  avg_list <- list()
  for (ps_col in pseudotime_col) {
    ps_bin_name <- paste0(ps_col, "_bins_", n_bins)
    
    avg_scores <- plot_df %>%
      dplyr::group_by(get(ps_bin_name)) %>%
      dplyr::select(all_of(genes)) %>%
      summarise_all(mean) %>%
      # Calculate relative deviations from the first bin
      mutate(across(all_of(genes), ~ . - dplyr::first(.))) 
    
    colnames(avg_scores)[1] <- "bin"
    avg_df <- reshape2::melt(avg_scores)
    avg_df$bin <- as.numeric(avg_df$bin)
    avg_df$group <- ps_col
    avg_list[[ps_col]] <- avg_df
  }
  
  plot_df <- do.call(rbind, avg_list)
  
  if (length(pseudotime_col) == 1) {
    p <- ggplot(plot_df, aes(x = bin, y = value, color = variable, 
                             fill = variable)) + 
      geom_smooth(size = line_size, se = se) + 
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey22") + 
      xlab("Pseudotime [bin]") + 
      ylab("Delta(Expression)") + 
      theme_test() +
      theme(
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        panel.border = element_rect(
          size = 1,
          fill = NA,
          color = "black"
        )
      )
  } else {
    p <- ggplot(plot_df, aes(x = as.numeric(bin), y = value, 
                             fill = group)) + 
      geom_smooth(se = se, size = line_size) + 
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey") + 
      xlab("Pseudotime") + ylab("Gene Expression") + 
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
            axis.line.x = element_blank(), axis.line.y = element_blank(), 
            panel.border = element_rect(size = 1, fill = NA, color = "black")) +
      labs(color = "Trajectory")
    
    if (!is.null(group_colors)) {
      p <- p + scale_color_manual(values = group_colors)
    }
  }
  
  if (patch) {
    patch <- p + facet_wrap(~variable, ncol = ncol, scales = "free")
    return(patch)
    
  } else {
    return(p)
  }
}



####
get_dds_resultsAvsB <- function(clustx, 
                                A, 
                                B, 
                                padj_cutoff = 0.05,
                                save_dir = "DESeq2/pairwise/results/") {
  print(clustx)
  idx <- which(names(counts_ls) == clustx)
  cluster_counts <- counts_ls[[idx]]
  cluster_metadata <- metadata_ls[[idx]]
  contrast <- paste(c("group_id", A, "vs", B), collapse = "_")
  
  # if ( all(colnames(cluster_counts) != rownames(cluster_metadata))) 
  #   print("ERROR: sample names in counts matrix columns and metadata rows do not match!")
  # }
  clustx <- janitor::make_clean_names(clustx)
  
  dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                colData = cluster_metadata, 
                                design = ~ group_id)
  rld <- vst(dds, blind = TRUE)
  
  DESeq2::plotPCA(rld, intgroup = "group_id")
  if (!dir.exists(save_dir)) { dir.create(save_dir) }
  ggsave(paste0(save_dir, clustx,"_",contrast, "_specific_PCAplot.png"))
  
  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)
  
  png(paste0(save_dir, clustx,"_",contrast, "_specific_heatmap.png"),
      height = 6, width = 7.5, units = "in", res = 300)
  pheatmap(rld_cor)
  dev.off()
  
  dds <- DESeq(dds, parallel = T, 
               BPPARAM = BiocParallel::MulticoreParam(32))
  
  png(paste0(save_dir, clustx,"_",contrast, "_dispersion_plot.png"),
      height = 5, width = 6, units = "in", res = 300)
  plotDispEsts(dds)
  dev.off()
  
  res <- results(dds, contrast = c("group_id", A, B), alpha = 0.05)
  # res <- lfcShrink(dds, coef = contrast, res = res)
  
  res_tbl <- res %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble()
  
  
  write.csv(res_tbl,
            paste0(save_dir, clustx, "_", contrast, "_all_genes.csv"),
            quote = FALSE, 
            row.names = FALSE)
  
  sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
    dplyr::arrange(padj)
  
  write.csv(sig_res,
            paste0(save_dir, clustx, "_", contrast, "_signif_genes.csv"),
            quote = FALSE, 
            row.names = FALSE)
  
  normalized_counts <- counts(dds, normalized = TRUE)
  
  top20_sig_genes <- sig_res %>%
    dplyr::arrange(padj) %>%
    dplyr::pull(gene) %>%
    head(n = 20)
  
  top20_sig_counts <- normalized_counts[rownames(normalized_counts) %in% top20_sig_genes, ]
  
  top20_sig_df <- data.frame(top20_sig_counts) %>% janitor::clean_names()
  top20_sig_df$gene <- rownames(top20_sig_counts)
  
  if(length(top20_sig_df$gene) >= 10){  
    top20_sig_df <- melt(setDT(top20_sig_df), 
                         id.vars = c("gene"),
                         variable.name = "clus_id") %>% 
      data.frame()
    
    top20_sig_df$clus_id <- gsub("\\.", "_", top20_sig_df$clus_id)
    
    df_sce <- as.data.frame(colData(dds))
    df_sce$clus_id <- janitor::make_clean_names(df_sce$clus_id)
    
    top20_sig_df <- plyr::join(top20_sig_df, df_sce,
                               by = "clus_id")
    
    ggplot(top20_sig_df, aes(y = value, x = group_id, col = group_id)) +
      geom_jitter(height = 0, width = 0.15) +
      scale_y_continuous(trans = 'log10') +
      ylab("log10 of normalized expression level") +
      xlab("condition") +
      ggtitle("Top 20 Significant DE Genes") +
      theme(plot.title = element_text(hjust = 0.5)) +
      facet_wrap(~ gene)
    
    ggsave(paste0(save_dir, clustx, "_", contrast, "_top20_DE_genes.png"))
  }
}

####

get_dds_LRTresults <- function(clustx,
                               multi_comparison=F,
                               dis_entities=NULL,
                               save_dir = "DESeq2/lrt/results/"){
  print(clustx) 
  idx <- which(names(counts_ls) == clustx)
  cluster_counts <- counts_ls[[idx]]
  cluster_metadata <- metadata_ls[[idx]]
  clustx <- janitor::make_clean_names(clustx)
  
  dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                colData = cluster_metadata, 
                                design = ~ group_id)
  dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ 1,
                   parallel = T, 
                   BPPARAM = BiocParallel::MulticoreParam(32))
  
  if(multi_comparison){
    if (!dir.exists(save_dir)) { dir.create(save_dir) }
    for(gp in dis_entities){
      dir <- paste0(save_dir,gp)
      
      if(!dir.exists(dir)){dir.create(dir)}
      res_LRT <- results(dds_lrt, contrast = c("group_id",gp,"Ctrl"))
      res_LRT_tb <- res_LRT %>%
        data.frame() %>%
        rownames_to_column(var = "gene") %>% 
        mutate(comparison = paste0(gp,"_vs_Ctrl")) %>%
        as_tibble()
      
      write.csv(res_LRT_tb,
                paste0(dir,"/", clustx, "_LRT_all_genes.csv"),
                quote = FALSE, 
                row.names = FALSE)
      sigLRT_genes <- res_LRT_tb %>% 
        dplyr::filter(padj < 0.05)
      write.csv(sigLRT_genes,
                paste0(dir,"/", clustx, "_LRT_signif_genes.csv"),
                quote = FALSE, 
                row.names = FALSE)
      
      if(length(sigLRT_genes$gene) >= 15){
        rld <- rlog(dds_lrt, blind = TRUE)
        rld_mat <- assay(rld)
        rld_cor <- cor(rld_mat)
        cluster_rlog <- rld_mat[sigLRT_genes$gene, ]
        rownames(cluster_metadata) <- cluster_metadata$clus_id
        cluster_meta_sig <- cluster_metadata[which(rownames(cluster_metadata) %in% colnames(cluster_rlog)), ]
        cluster_meta_sig$np <- cluster_meta_sig$group_id
        cluster_groups <- suppressMessages(degPatterns(cluster_rlog, metadata = cluster_meta_sig,
                                                       time = "group_id", col = 'np'))
        cluster_groups$plot+scale_color_manual(values = group_cols) + theme_bw() + labs(title=paste0(gp,": ",clustx))
        ggsave(paste0(dir,"/", clustx, "_LRT_DEgene_groups.png"))
        write.csv(cluster_groups$df,
                  paste0(dir,"/", clustx, "_LRT_DEgene_groups.csv"),
                  quote = FALSE, 
                  row.names = FALSE)
        saveRDS(cluster_groups,paste0(dir,"/",clustx, "_LRT_DEgene_groups.rds"))
      }
      save(dds_lrt, res_LRT, sigLRT_genes, 
           file = paste0(dir,"/", clustx, "_all_LRTresults.Rdata"))
    }
  }else{
    
    res_LRT <- results(dds_lrt)
    res_LRT_tb <- res_LRT %>%
      data.frame() %>%
      rownames_to_column(var = "gene") %>% 
      as_tibble()
    
    if (!dir.exists(save_dir)) { dir.create(save_dir) }
    write.csv(res_LRT_tb,
              paste0(save_dir, clustx, "_LRT_all_genes.csv"),
              quote = FALSE, 
              row.names = FALSE)
    sigLRT_genes <- res_LRT_tb %>% 
      dplyr::filter(padj < 0.05)
    write.csv(sigLRT_genes,
              paste0(save_dir, clustx, "_LRT_signif_genes.csv"),
              quote = FALSE, 
              row.names = FALSE)
    
    if(length(sigLRT_genes$gene) >= 15){
      rld <- rlog(dds_lrt, blind = TRUE)
      rld_mat <- assay(rld)
      rld_cor <- cor(rld_mat)
      cluster_rlog <- rld_mat[sigLRT_genes$gene, ]
      rownames(cluster_metadata) <- cluster_metadata$clus_id
      cluster_meta_sig <- cluster_metadata[which(rownames(cluster_metadata) %in% colnames(cluster_rlog)), ]
      cluster_meta_sig$np <- cluster_meta_sig$group_id
      cluster_groups <- degPatterns(cluster_rlog, metadata = cluster_meta_sig,
                                    time = "group_id", col = 'np')
      cluster_groups$plot+scale_color_manual(values = group_cols) + theme_bw()
      ggsave(paste0(save_dir, clustx, "_LRT_DEgene_groups.png"))
      write.csv(cluster_groups$df,
                paste0(save_dir, clustx, "_LRT_DEgene_groups.csv"),
                quote = FALSE, 
                row.names = FALSE)
      saveRDS(cluster_groups, paste0(save_dir, clustx, "_LRT_DEgene_groups.rds"))
    }
    save(dds_lrt, cluster_groups, res_LRT, sigLRT_genes, 
         file = paste0(save_dir, clustx, "_all_LRTresults.Rdata"))
  }
}


####
volcano_from_DESeq <- function(res_tbl = res_tbl,
                               padj_tresh = 0.05,
                               log2fc = 0.5,
                               top_x_to_plot = 15,
                               title = "GE in Tau vs Ctrl")
                               
{
  # Ensure the specified columns exist in the data frame
  if (!all(c(padj_col, log2FoldChange_col) %in% names(res_tbl))) {
    stop("Specified columns not found in the data frame")
  }

  res_table_thres <- res_tbl[!is.na(res_tbl[[padj_col]]), ] %>% 
    mutate(threshold = !!sym(padj_col) < padj_thresh & abs(!!sym(log2FoldChange_col)) >= log2fc)
  
  min_padj <- min(log10(res_table_thres[[padj_col]]))
  
  res_table_thres <- res_table_thres %>%
    mutate(gene_type = case_when(
      !!sym(log2FoldChange_col) >= 1 & !!sym(padj_col) <= padj_thresh ~ "up",
      !!sym(log2FoldChange_col) <= -1 & !!sym(padj_col) <= padj_thresh ~ "down",
      TRUE ~ "ns"
    ))   
  cols <- c("up" = cols[1], "down" = cols[2], "ns" = cols[3]) 
  sizes <- c("up" = 2, "down" = 2, "ns" = 1) 
  alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)
  sig_il_genes <- subset(res_table_thres, gene_type %in% c('up','down')) %>% as.data.frame %>%
    mutate(selector = abs(log2FoldChange) * -log10(padj)) %>%
    arrange(selector, .by_group = TRUE) %>%
    top_n(top_x_to_plot)
  p1 <- ggplot(data = res_table_thres,
               aes(x = log2FoldChange,
                   y = -log10(padj))) + 
    geom_point(aes(fill = gene_type), 
               color='grey42',
               alpha = 0.5, 
               shape = 21,
               size = 2) + 
    geom_hline(yintercept = -log10(padj_tresh),
               linetype = "dashed") + 
    geom_vline(xintercept = c(log2(0.5), log2(2)),
               linetype = "dashed") +
    geom_label_repel(data = sig_il_genes,  
                     aes(label = gene),
                     force = 2,max.overlaps = 50,
                     nudge_y = 5, nudge_x = 0) +
    scale_fill_manual(values = cols) + 
    scale_x_continuous(breaks = c(seq(-10, 10, 2)),     
                       limits = c(-8, 8))  +
    labs(title=title,
         x = "log2(fold change)",
         y = "-log10(adjusted P-value)",
         colour = "Expression \nchange") +
    theme_bw()  
  return(p1)
}


volcano <- function(res_tbl,
                    thresh_plot_1 = 0.1,
                    thresh_plot_2 = 5e-5,
                    ext_y = 10,
                    log2fc = 0.2,
                    top_x_to_plot = 15,
                    text_size = 3,
                    title = "GE in Tau vs Ctrl",
                    padj_col = "padj",
                    log2FoldChange_col = "log2FoldChange",
                    cols = c("#ffad73","#26b3ff","grey"),
                    type = "gene"
)
{
  # Ensure the specified columns exist in the data frame
  if (!all(c(padj_col, log2FoldChange_col) %in% names(res_tbl))) {
    stop("Specified columns not found in the data frame")
  }
  
  res_table_thres <- res_tbl[!is.na(res_tbl[[padj_col]]), ] %>% 
    mutate(threshold = !!sym(padj_col) < thresh_plot_1 & abs(!!sym(log2FoldChange_col)) >= log2fc) %>%
    complete()
  
  min_padj <- is.finite(-log10(res_table_thres[[padj_col]]))
  min_padj <- max(-log10(res_table_thres[[padj_col]])[min_padj], na.rm = T)
  
  max_fc <- max(abs(res_table_thres[[log2FoldChange_col]]))
  
  res_table_thres <- res_table_thres %>%
    mutate(gene_type = case_when(
      !!sym(log2FoldChange_col) >= log2fc & !!sym(padj_col) <= thresh_plot_1 ~ "up",
      !!sym(log2FoldChange_col) <= -log2fc & !!sym(padj_col) <= thresh_plot_1 ~ "down",
      TRUE ~ "ns"
    ))   
  
  cols <- c("up" = cols[1], "down" = cols[2], "ns" = cols[3]) 
  sizes <- c("up" = 2, "down" = 2, "ns" = 1) 
  alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)
  
  sig_il_genes <- subset(res_table_thres, gene_type %in% c('up','down')) %>% as.data.frame() %>%
    mutate(selector = abs(!!sym(log2FoldChange_col)) * -log10(!!sym(padj_col))) %>%
    arrange(desc(selector))
  
  x <- subset(sig_il_genes, gene_type == "up") %>%
    arrange(desc(selector)) %>% 
    head(top_x_to_plot)
  y <- subset(sig_il_genes, gene_type == "down") %>%
    arrange(desc(selector)) %>% 
    head(top_x_to_plot)
  
  sig_il_genes <- rbind(x,y)
  
  p1 <- ggplot(data = res_table_thres,
               aes(x = !!sym(log2FoldChange_col),
                   y = -log10(!!sym(padj_col)))) + 
    geom_point(aes(fill = gene_type), 
               color = 'grey22',
               alpha = 0.5, 
               shape = 21,
               size = 2) + 
    geom_hline(yintercept = -log10(thresh_plot_1),
               linetype = "dashed") + 
    geom_hline(yintercept = -log10(thresh_plot_2),
               linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = c(-(log2fc), log2fc),
               linetype = "dashed") +
    geom_text_repel(data = subset(sig_il_genes, gene_type == 'up'),
                    aes(label = gene, 
                        segment.square = T,
                        segment.inflect = F),
                    force = 0.5, 
                    size = text_size,
                    fontface = "italic", 
                    max.overlaps = Inf,
                    nudge_x = max_fc - x[[log2FoldChange_col]] + 0.3,  # 
                    segment.size = 0.2,
                    segment.curvature = -0.01,
                    direction = "y") +
    geom_text_repel(data = subset(sig_il_genes, gene_type == 'down'),
                    aes(label = gene,
                        segment.square = T,
                        segment.inflect = F),
                    force = 0.5,
                    size = text_size,
                    fontface = "italic", 
                    max.overlaps = Inf,
                    nudge_x = -max_fc - y[[log2FoldChange_col]] - 0.3,  # 
                    segment.size = 0.2,
                    segment.curvature = 0.01,
                    direction = "y") +
    
    scale_fill_manual(values = cols) + 
    scale_x_continuous(breaks = c(seq(-10, 10, 0.5)),     
                       limits = c(-max_fc-0.5,max_fc+0.5)) +
    scale_y_continuous(breaks = c(seq(0, 1000, 10)),
                       limits = c(0,min_padj+ext_y)) +
    labs(title = title,
         x = "log2(fold change)",
         y = "-log10(P-value)",
         colour = "Expression \nchange") +
    theme_bw()  +
    theme(legend.position = "none")
  
  return(p1)
}


bidirect_plot <- function(df = df, 
                          ct = "celltype", 
                          multiple_p_vals = F,
                          col_categ = "Category", 
                          padj_tresh = 0.1,
                          log2FoldChange_col="avg_log2FC"){
 
  # Bar Diagram
  get_category_counts <- function(data, p_value_col, cutoff, label) {
    data %>%
      dplyr::filter(!!sym(p_value_col) < cutoff) %>%
      mutate(Regulation = ifelse(!!sym(log2FoldChange_col) > 0, "Up", "Down")) %>%
      dplyr::count(Category, Regulation) %>%
      mutate(n = ifelse(Regulation == "Down", -n, n),
             Cutoff = label)
  }
  
  if(multiple_p_vals){
    combined_counts <- get_category_counts(df, "p_val_adj", 0.05, "Adjusted p < 0.05")
    category_counts1 <- get_category_counts(df, "p_val_adj", 0.1, "Adjusted p < 0.1")
    category_counts2 <- get_category_counts(df, "p_val_adj", 0.2, "Adjusted p < 0.2")
    
    combined_counts <- bind_rows(combined_counts, category_counts1, category_counts2)
  }else{
    combined_counts <- get_category_counts(df, "p_val_adj", padj_tresh, "Adjusted p < 0.1")  
  }

  total_counts <- combined_counts %>%
    group_by(Category) %>%
    summarize(total = sum(abs(n)))
  
  combined_counts$Category <- factor(combined_counts$Category, 
                                     levels = total_counts$Category[order(total_counts$total, decreasing = FALSE)])
  
  p <- ggplot(combined_counts, aes(x = n, y = Category, fill = Category, group = Cutoff)) + 
    geom_col(position = position_dodge(width = 0.7), color = "grey7") + 
    geom_text(aes(label = abs(n), x = n), 
              position = position_dodge(width = 0.7),
              hjust = ifelse(combined_counts$n > 0, -0.3, 1.3),
              size = 3.5) +
    scale_x_continuous(labels = abs, n.breaks = 20) + 
    labs(title = paste0("Number of DEGs by ", ct), 
         x = "Number of DEGs", 
         y = "Cell Types") + 
    theme_bw() +
    theme(legend.position = "bottom") 
  
 
  if(col_categ == "Cutoff"){
    p <- p + scale_fill_manual(values =c("#4774a8", "#B53737"))
    }
                    
 
 return(p)
}

####
prep_deseq2 <-
  function(combrna,
           feature_cell_identity = "celltype",
           feature_sample = "Pub_ID",
           feature_group = "NPdiagnosis",
           cores = 16, 
           single_ct = F) {
    suppressPackageStartupMessages(library(SingleCellExperiment))
    suppressPackageStartupMessages(library(Seurat))
    suppressPackageStartupMessages(library(SeuratObject))
    suppressPackageStartupMessages(library(DESeq2))
    suppressPackageStartupMessages(library(tidyverse))
    
    # Prepare data slots
   if(combrna@active.assay == "chromvar"){
     counts <- combrna@assays$chromvar$data
   }else{
     counts <- combrna@assays$RNA$data 
   }
    metadata <- combrna@meta.data %>% as.data.frame()
    metadata$cluster_id <- c(unlist(combrna[[feature_cell_identity]])) #%>% stringr::str_replace_all(pattern = "_",replacement = "-")
    metadata$temp2 <- c(unlist(combrna[[feature_sample]]))
    metadata$clus_ct_id <- paste0(metadata$cluster_id, "__", metadata$temp2)
    
    # Create single cell experiment object
    sce <- SingleCellExperiment(assays = list(counts = counts),
                                colData = metadata)
    aggr_counts <- Matrix.utils::aggregate.Matrix(t(counts(sce)),
                                                  groupings = sce$clus_ct_id,
                                                  fun = "sum")
    
    # Extract sample-level variables
    meta <- c(rownames(aggr_counts)) %>% 
      # stringi::stri_replace_last_regex("_", "-") %>% stringi::stri_split_regex(".")
      as.data.frame() %>%
      separate(
        col = '.',
        into = c('celltype', 'feature_sample'),
        sep = "__"
      )
    meta$clus_np_id <- rownames(aggr_counts)
    meta <- lapply(meta, as.factor) %>% as.data.frame()
    if(combrna@active.assay == "RNA"){
      meta <- left_join(meta, unique(metadata[,c(feature_sample,feature_group,
                                                 "No_TA_3.2e7mcsq", "TA_load", "No_NFT_3.2e7mcsq", "NFT_load", "CB_load", "Threads_load")]), 
                        by = c("feature_sample" = colnames(metadata)[colnames(metadata) == feature_sample])) 
      input_counts <- round(t(aggr_counts))
    }else{
      meta$NPdiagnosis <- ifelse(grepl("^P", meta$feature_sample), "PSP", "Ctrl")
      input_counts <- scales::rescale(as.matrix(round(t(aggr_counts))), to = c(0,10000)) %>% round
    }
    
    # Create DESeqDataSet object
    if(single_ct){
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = input_counts,
                                            colData = meta,
                                            design = ~ NPdiagnosis)
    }else{
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = input_counts,
                                  colData = meta,
                                  design = ~ NPdiagnosis + celltype)
    }
    dds <- DESeq(dds, BPPARAM = cores)
  
  # Extract results and perform PCA plot
  res <- results(dds)
  rld <- vst(dds, blind = TRUE, nsub = 100) # changed 1000 -> 100 for tfme
  pca_plot <- DESeq2::plotPCA(rld, ntop = 10000, intgroup = "celltype") + 
    ggrepel::geom_label_repel(label = meta$feature_sample) + 
    theme_bw() + 
    coord_equal(1)
  
  # Return list of results
  return(list(dds = dds, res = res, rld = rld, meta = meta, pca_plot = pca_plot))
}


### 
subclust_celltype <- function(combrna,
                              cell_type,
                              city_neighborhood = city_neighborhood,
                              cores = 4) {
  cat("Processing cell type:", cell_type, "\n")
  RhpcBLASctl::blas_set_num_threads(1)
  
  # Subset the Seurat object "combrna" into the current cell type
  
  # Special case: combine VLMC and Endo cell types
  if (cell_type %in% c("Endo","VLMC")){
    current_cell_type <- subset(combrna, predicted.subclass %in% c("Endo", "VLMC"))
    current_cell_type$predicted.subclass <- "Endo_VLMC"
  # }else{  
    current_cell_type <- subset(combrna, predicted.subclass == cell_type)
  # }
  DefaultAssay(current_cell_type) <- "RNA"
  # Step 1: Remove genes expressed in fewer than 15 nuclei and non-coding and non-annotated genes
  current_cell_type <- current_cell_type[rowSums(current_cell_type@assays$RNA$counts > 0) >= 15, ] 
  current_cell_type <- current_cell_type[!grepl("^(AC\\d{3}|AL\\d{3}|AP\\d{3}|LINC\\d{3})", rownames(current_cell_type@assays$RNA$counts)), ]
  
  # Step 2: cluster robustness 
  clusters <- scSHC::scSHC(current_cell_type@assays$RNA$counts, 
                           batch = current_cell_type$Pub_ID, 
                           cores = cores)
  saveRDS(clusters,paste0(janitor::make_clean_names(cell_type),"_clusters.Rds"))
  current_cell_type$subcluster_scSHC <- clusters[[1]]
  message("Clusters computed.")
  
  if(ncol(current_cell_type) > 250){
    message("Path A\n")
    
    # # Iteratively run the normalization and clustering pipeline
    current_cell_type <- SCTransform(current_cell_type, variable.features.n = ifelse(cell_type %in% c("Endo", "VLMC"), 2000, 4000))
    current_cell_type <- RunPCA(current_cell_type, npcs = 50)
    current_cell_type <- FindNeighbors(current_cell_type, dims = 1:50)
    current_cell_type <- FindClusters(current_cell_type, resolution = 0.5, algorithm = 4, method = "igraph")
    current_cell_type <- RunUMAP(current_cell_type, dims = 1:50, min.dist = 0.1)
    current_cell_type$celltype.subclusters <- paste0(cell_type,".",current_cell_type$seurat_clusters)
    
    # Join clusters with less than three contributor cases
    for(sc in levels(as.factor(current_cell_type$subcluster_scSHC))){
      df_temp <- subset(current_cell_type@meta.data, subcluster_scSHC == sc)
      if(sum(table(df_temp$Pub_ID) > 5) < 3){
        top <-TopNeighbors(city_neighborhood@neighbors$SCT.nn, cell=colnames(current_cell_type), n = 10)
        new_subclass <- subset(current_cell_type@meta.data[top,], !subcluster_scSHC == sc)$subcluster_scSHC %>% 
          table(.) %>% which.max(.) %>% names()
        current_cell_type$subcluster_scSHC[which(current_cell_type$subcluster_scSHC == sc)] <- new_subclass
      }
    }
    current_cell_type$subcluster_scSHC <- as.factor(current_cell_type$subcluster_scSHC) %>% as.numeric()
    clOrder<-c(paste(cell_type, 1:length(levels(as.factor(current_cell_type$subcluster_scSHC))), sep="."))
    current_cell_type$subcluster_scSHC <- as.factor(current_cell_type$subcluster_scSHC)
    
    saveRDS(object = current_cell_type, file = paste0(janitor::make_clean_names(cell_type),"_subclust.Rds"))
    return(current_cell_type)
    
  }else{
    if(ncol(current_cell_type) > 25){
      message("Path B: low-abundant cell types\n")
      # # Special case: low abundant subclasses:
      top <-TopNeighbors(city_neighborhood@neighbors$SCT.nn, cell=colnames(current_cell_type),n=10)
      new_subclass <- subset(combrna@meta.data[top,], !predicted.subclass == cell_type)$predicted.subclass %>%
        table(.) %>% which.max(.) %>% names()
      current_cell_type$subcluster_scSHC <- paste0(new_subclass, "_special_NN")

      current_cell_type <- SCTransform(current_cell_type, variable.features.n = ifelse(cell_type %in% c("Endo", "VLMC"), 2000, 4000))
      current_cell_type <- RunPCA(current_cell_type, npcs = 30)
      current_cell_type <- FindNeighbors(current_cell_type, dims = 1:30)
      current_cell_type <- RunUMAP(current_cell_type, dims = 1:30, min.dist = 0.1)
      current_cell_type$celltype.subclusters <- paste0(cell_type,".1")
      saveRDS(object = current_cell_type, file = paste0(janitor::make_clean_names(cell_type),"_subclust.Rds"))
      return(current_cell_type)
      
    }else{
      message(cell_type, " omitted. Has only ", ncol(current_cell_type), " cells.")
    }
  }
  cat("\n")
}

###
compute_go <-
  function(de_results,
           direction = "up",
           universe = all_genes,
           mode = "GO",
           ont = "BP",
           pAdjustMethod = "BH",
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.1,
           readable = T,
           write_table = T) {
    library(org.Hs.eg.db)
    # Ensure that the DE results dataframe has the expected format
    if (!("p_val_adj" %in% colnames(de_results)) ||
        !("avg_log2FC" %in% colnames(de_results))) {
      stop("The input data frame does not have the expected format.")
    }
    
    # Filter significant genes based on adjusted p-value
    if (direction == "up") {
      sig_genes <-
        de_results[de_results$p_val < pvalueCutoff &
                     de_results$avg_log2FC > 0,]
    } else{
      sig_genes <-
        de_results[de_results$p_val < pvalueCutoff &
                     de_results$avg_log2FC < 0,]
    }
    
    if(mode == "GO"){
      # Run GO enrichment analysis
      ego <- enrichGO(
        gene = rownames(sig_genes),
        universe =  rownames(de_results),
        keyType = "SYMBOL",
        OrgDb = org.Hs.eg.db,
        ont = ont,
        pAdjustMethod = pAdjustMethod,
        pvalueCutoff = pvalueCutoff,
        qvalueCutoff = qvalueCutoff,
        readable = readable
      )   
    }else if(mode == "KEGG"){
      gene_id <- mapIds(org.Hs.eg.db, keys = rownames(sig_genes),
                        column = "ENTREZID", keytype = "SYMBOL")
      gene_id_univ <- mapIds(org.Hs.eg.db, keys = rownames(de_results),
                        column = "ENTREZID", keytype = "SYMBOL")
      ego <- enrichKEGG(
        gene = unname(gene_id),
        universe =  unname(gene_id_univ),
        organism = "hsa",
        keyType = "ncbi-geneid",
        pAdjustMethod = pAdjustMethod,
        pvalueCutoff = qvalueCutoff,
        minGSSize = 5,
        maxGSSize = 500,
        use_internal_data = FALSE
      )
    }
    
    # Simplify the enrichment results to reduce redundancy
    cluster_summary <- data.frame(ego)
    tryCatch({
      cluster_summary$celltype_id <- as.factor(de_results$celltype_id[1])
      
      if (write_table) {
        saveRDS(ego, file = paste0(de_results$celltype_id[1], "_", direction, "_enrich",mode,".Rds"))
        write.csv(
          cluster_summary,
          file = paste0(de_results$celltype_id[1], "_", direction, "_enrich",mode,".csv"),
          row.names = F
        )
      }
    }, error = function(e) {
      message("Combination: ",
              de_results$celltype_id[1],
              " & ",
              direction,
              paste0(" not giving any significant ",mode," enrichment."))
      return(cluster_summary)
    })
    
    return(cluster_summary)
  }

###
reorder_description <-
  function(description, celltype_id_order) {
    factor(description, levels = unique(description[order(celltype_id_order, decreasing = TRUE)]))
  }

###
plot_go_facet <-
  function(df,
           facet_var = NULL,
           ncol = 3,
           x = NULL,
           y = NULL,
           rel = F,
           col_lib = scl_cols_nam_cl,
           direct = "up") {
    library(stringi)
    if(!is.null(facet_var)){
      # Get unique levels of the facetting variable
      facet_levels <- unique(df[[facet_var]]) %>% as.factor() %>% levels
      
      # Create an empty list to store the individual plots
      plot_list <- list()
      
      # Loop through each facet level and create a plot for each subset
      for (level in facet_levels) {
        # Subset the data by the current facet level
        subset_df <- df[df[[facet_var]] == level, ]
        col_sub <- col_lib[unique(subset_df$celltype_id)]
        
        # Create the plot for the current subset
        h <- length(unique(df$Description))
        w <- length(unique(df$celltype_id))
        
        plot <- ggplot(
          subset_df,
          aes(
            x = stri_sort(celltype_id),
            y = reorder_description(Description, celltype_id_order),
            size = -log10(p.adjust),
            fill = celltype_id,
            shape = direct
          ),
          group = class
        ) +
          geom_point(alpha = 0.7, color = "grey4") +
          scale_size_continuous(range = c(1, 4)) +
          scale_shape_manual(values = c(21)) + #25,24 if down/up
          scale_fill_manual(values = col_sub, guide = "none") +
          labs(
            x = "Cell identity",
            y = "",
            size = "-log10(Adjusted P-value)",
            fill = "Cell Type"
          ) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),
                legend.position = "bottom")
        
        # Add the plot to the plot list
        plot_list[[level]] <- plot
      }
      
      if (!rel) {
        return(plot_list)
        
      } else{
        heights <- sapply(plot_list, function(plot) {
          length(unique(plot$data$Description))
        })
        heights <- heights / max(heights)
        
        # Calculate the relative widths based on the number of items in the celltype_id column
        widths <- sapply(plot_list, function(plot) {
          length(unique(plot$data$celltype_id))
        })
        widths <- widths / max(widths)
        
        # Combine the plots using plot_grid()
        combined_plot <- cowplot::plot_grid(
          plotlist = plot_list,
          ncol = ncol,
          rel_heights = heights,
          rel_widths = widths,align = "hv"
        )
        return(combined_plot)
      }
    }else{
      plot <- ggplot(
        df,
        aes(
          x = module,
          y = reorder_description(df$Description, df$module),
          size = -log10(p.adjust),
          fill = Count
        )) +
        geom_point(alpha = 0.7, shape = 21, color = "grey4") +
        colorspace::scale_fill_continuous_sequential(palette = "OrRd", p1 = 0.5) +
        labs(
          x = "Cell identity",
          y = "",
          size = "-log10(Adjusted P-value)",
          fill = "Gene Count"
        ) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "right") + coord_fixed()
      
      return(plot)
    }
    
  }


# define function
lollipop_go <- function(ct, lev, dir = "../") {
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(DOSE)
  library(ggupset)
  library(pathview)
  de_results <-
    read_csv(paste0(dir,lev,"/", ct)) %>% 
    column_to_rownames(colnames(.)[1]) %>% 
    mutate(celltype_id = ct)
  
  k_up <- compute_go(de_results, direction = "up", mode = "KEGG", pvalueCutoff = 0.05, write_table = F) %>% mutate(direct = "up", ONTOLOGY = "KEGG")
  k_down <- compute_go(de_results, direction = "down", mode = "KEGG", pvalueCutoff = 0.05, write_table = F) %>% mutate(direct = "down", ONTOLOGY = "KEGG")
  GO_up <- compute_go(de_results, direction = "up", mode = "GO", ont = "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.1, write_table = F) %>% mutate(direct = "up")
  GO_down <- compute_go(de_results, direction = "down", mode = "GO", ont = "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.1, write_table = F) %>% mutate(direct = "down")
  
  sel = c("ONTOLOGY","ID","Description", "GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","celltype_id","direct")
  out <- rbind(k_up %>% dplyr::select(any_of(sel)),
               k_down %>% dplyr::select(any_of(sel)),
               GO_down %>% dplyr::select(any_of(sel)),
               GO_up %>% dplyr::select(any_of(sel))) %>% 
    data.frame() %>%
    mutate(celltype_id = stringr::str_replace_all(.$celltype_id,
                                                  pattern = "[ -]",
                                                  replacement = "_")) %>% 
    mutate(celltype_id = gsub("\\.csv.*", "", celltype_id)) 
  write.csv(out, paste0(out$celltype_id[1],"joined.csv"))
  
  col_grep <- ls <- list(celltype = ct_cols, predicted.subclass = sclass_cols, subcluster = sclus_cols)[[lev]]
  
  # filter, wrangle
  dfd <- out %>% group_by(ONTOLOGY,direct) %>% 
    slice_min(order_by = p.adjust, n = 5) %>% # adjust top vars
    ungroup() %>%
    mutate_at(vars(ONTOLOGY, ID, Description, geneID), as.factor) %>%
    mutate_at(vars(GeneRatio, BgRatio), as.factor) %>%
    mutate_at(vars(pvalue, p.adjust, qvalue, Count), as.numeric) %>% 
    left_join(., unique(data.frame(class = combrna$predicted.class, 
                                   celltype = combrna$celltype, 
                                   subclass = combrna$predicted.subclass,
                                   subcluster = combrna$subcluster,
                                   row.names = NULL)), 
              by = c("celltype_id"=lev)) %>% 
    mutate(celltype_id_order = factor(celltype_id, levels = names(col_grep)),
           Count_dir = ifelse(direct == "up", Count, -Count))
  
  p <- ggplot(dfd, aes(x = Count_dir, 
                       y = tidytext::reorder_within(Description, Count_dir, list(ONTOLOGY)), 
                       size = -log(pvalue), 
                       fill = celltype_id)) + 
    geom_vline(xintercept = 0, alpha = 0.6) +
    geom_segment(aes(y = tidytext::reorder_within(Description, Count_dir,  list(ONTOLOGY)), 
                     xend=0, xend=Count_dir), size = 0.5, color="grey4") +
    geom_point(stat = "identity", shape = 21) + facet_grid(ONTOLOGY~., drop = T, space = "free", scales = "free") +
    theme_bw() + scale_fill_manual(values = col_grep)  +
    labs(title = "",
         y = "",
         x = "Gene Count, DEG DOWN vs. UP") + 
    scale_y_discrete(limits=rev) + 
    theme(strip.text.y = element_text(angle = 0))
  ggsave(paste0("plots/go_dir_",dfd$celltype_id[1],".pdf"),
         plot = p, height = 6, width = 8)
  return(p)
}


## modif psupertime plot
make_col_vals <- function(y_labels, palette='RdBu') {
  n_labels 	= length(levels(y_labels))
  require(RColorBrewer)
  max_col 	= 11
  if (n_labels==1) {
    col_vals 	= brewer.pal(3, palette)
    col_vals 	= col_vals[1]
  } else if (n_labels==2) {
    col_vals 	= brewer.pal(3, palette)
    col_vals 	= col_vals[-2]
  } else if (n_labels<=max_col) {
    col_vals 	= brewer.pal(n_labels, palette)
  } else {
    col_pal 	= brewer.pal(max_col, palette)
    col_vals 	= colorRampPalette(col_pal)(n_labels)
  }
  col_vals 	= rev(col_vals)
  
  return(col_vals)
}

plot_identified_genes_over_psupertime_modif <- function (psuper_obj = scepsu, 
                                                         label_name = "Ordered labels", 
                                                         n_to_plot = 10, 
                                                         size = 0.75,
                                                         alpha = 0.5,
                                                         se = F,
                                                         sel = 7000,
                                                         palette = "RdBu", 
                                                         plot_ratio = 1.25,
                                                         density = F,
                                                         hex = 50) {
  require(scales)
  sel <- sample(row.names(psuper_obj$proj_dt),sel) %>% as.numeric()
  proj_dt = psuper_obj$proj_dt[sel,]
  beta_dt = psuper_obj$beta_dt
  x_data = psuper_obj$x_data[sel,]
  ps_params = psuper_obj$ps_params
  beta_nzero = beta_dt[abs_beta > 0]
  n_nzero = nrow(beta_nzero)
  top_genes = as.character(beta_nzero[1:min(n_to_plot, nrow(beta_nzero))]$symbol)
  plot_wide = cbind(proj_dt, data.frame(x_data[, top_genes, drop = FALSE]))
  plot_dt = melt.data.table(plot_wide, id = c("cell_id", "psuper", 
                                              "label_input", "label_psuper"), measure = top_genes, 
                            variable.name = "symbol")
  plot_dt[, `:=`(symbol, factor(symbol, levels = top_genes))]
  col_vals = make_col_vals(plot_dt$label_input, palette)
  n_genes = length(top_genes)
  ncol = ceiling(sqrt(n_genes * plot_ratio))
  nrow = ceiling(n_genes/ncol)
  
  if(!density){
  g = ggplot(plot_dt) + 
    aes(x = psuper, y = value) + 
    geom_point(size = size, aes(fill = label_input), shape = 21, color = "grey42", alpha = alpha) + 
    geom_smooth(se = se, 
                colour = "black") + 
    scale_fill_manual(values = col_vals) + 
    scale_shape_manual(values = c(1, 16)) + 
    scale_x_continuous(breaks = pretty_breaks()) + 
    scale_y_continuous(breaks = pretty_breaks()) + 
    facet_wrap(~symbol, 
               scales = "free_y", nrow = nrow, ncol = ncol) + theme_bw() + 
    theme(axis.text.x = element_blank()) + labs(x = "psupertime", 
                                                y = "z-scored log2 expression", colour = label_name)
  }else{
    g = ggplot(plot_dt) + 
      aes(x = psuper, y = value) +
      geom_hex(bins = hex, linewidth = size) +
      # scale_fill_continuous(type = "viridis") +
      ggpubr::gradient_fill("YlOrRd") +
      geom_smooth(se = se, 
                  colour = "black") + 
      scale_shape_manual(values = c(1, 16)) + 
      scale_x_continuous(breaks = pretty_breaks()) + 
      scale_y_continuous(breaks = pretty_breaks()) + 
      facet_wrap(~symbol, 
                 scales = "free_y", nrow = nrow, ncol = ncol) + theme_bw() + 
      theme(axis.text.x = element_blank()) + labs(x = "psupertime", 
                                                  y = "z-scored log2 expression", colour = label_name)
  }
  
  
  
  
  
  return(g)
}


plot_gene_psuper <- function(psuper_obj = scepsu, 
                             extra_genes = "MAPT",
                             sel = 3000,
                             size = 0.75,
                             alpha = 0.5,
                             se = F,
                             label_name = "Ordered labels", 
                             palette = "RdBu", 
                             plot_ratio = 1.25) 
{
  sel <- sample(row.names(psuper_obj$proj_dt),sel) %>% as.numeric()
  
  proj_dt = psuper_obj$proj_dt[sel,]
  beta_dt = psuper_obj$beta_dt
  x_data = psuper_obj$x_data[sel,]
  ps_params = psuper_obj$ps_params
  extra_genes = intersect(extra_genes, colnames(x_data))
  if (length(extra_genes) == 0) {
    warning("genes not found; did not plot")
    return()
  }
  plot_wide = cbind(proj_dt, data.frame(x_data[, extra_genes, 
                                               drop = FALSE]))
  plot_dt = melt.data.table(plot_wide, id = c("psuper", "label_input", 
                                              "label_psuper"), measure = extra_genes, variable.name = "symbol")
  plot_dt[, `:=`(symbol, factor(symbol, levels = extra_genes))]
  col_vals = make_col_vals(plot_dt$label_input, palette)
  n_genes = length(extra_genes)
  ncol = ceiling(sqrt(n_genes * plot_ratio))
  nrow = ceiling(n_genes/ncol)
  g = ggplot(plot_dt) + 
    aes(x = psuper, y = value) + 
    geom_point(size = size, aes(fill = label_input), shape = 21, color = "grey42", alpha = alpha) + 
    geom_smooth(se = se, 
                colour = "black") + 
    scale_fill_manual(values = col_vals) + 
    scale_shape_manual(values = c(1, 16)) + 
    scale_x_continuous(breaks = pretty_breaks()) + 
    scale_y_continuous(breaks = pretty_breaks()) + 
    facet_wrap(~symbol, 
               scales = "free_y", nrow = nrow, ncol = ncol) + theme_bw() + 
    theme(axis.text.x = element_blank()) + labs(x = "psupertime", 
                                                y = "z-scored log2 expression", colour = label_name)
  return(g)
}


calc_clusters_dt <- function(hclust_obj, x_data, proj_dt, k=5) {
  # make thing
  clusters_dt 	= data.table( h_clust=cutree(hclust_obj, k=k), symbol=colnames(x_data))
  # add clustering
  clusters_dt[, N:=.N, by=h_clust ]
  
  # order by correlation with psupertime
  temp_dt 			= data.table(melt(as.matrix(x_data), varnames=c('cell_id', 'symbol')))
  temp_dt 			= clusters_dt[ temp_dt, on='symbol' ]
  means_dt 			= temp_dt[, list(mean=mean(value)), by=list(cell_id, h_clust) ]
  means_dt 			= proj_dt[ means_dt, on='cell_id' ]
  corrs_dt 			= means_dt[, list( cor=cor(mean, psuper) ), by=h_clust]
  setorder(corrs_dt, cor)
  corrs_dt[, clust := 1:.N ]
  corrs_dt[, clust := factor(clust)]
  
  # add clusters ordered by size back in
  clusters_dt			= corrs_dt[ clusters_dt, on='h_clust' ]
  clusters_dt[, clust_label := factor(sprintf('%02d (%d genes)', clust, N)) ]
  clusters_dt[, h_clust := NULL ]
  setorder(clusters_dt, clust, symbol)
  
  return(clusters_dt)
}


do_topgo_for_cluster <- function (clusters_dt, sig_cutoff, org_mapping) 
{
  all_clusters = unique(clusters_dt[N >= sig_cutoff]$clust)
  go_results = data.table()
  message(sprintf("calculating GO enrichments for %d clusters:", 
                  length(all_clusters)))
  for (c in all_clusters) {
    message(".", appendLF = FALSE)
    gene_list = subset(clusters_dt,clust == c)$symbol

    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(enrichplot)
    library(DOSE)
    library(ggupset)
    library(pathview)

    # Run GO enrichment analysis
     ego <- enrichGO(
        gene = gene_list,
        universe =  clusters_dt$symbol,
        keyType = "SYMBOL",
        OrgDb = org.Hs.eg.db,
        ont = c("ALL"), 
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.1,
        readable = T
      )
      
     temp_results <- data.frame(ego)
        
    go_results = rbind(go_results, temp_results)
  }
  message("")
  return(go_results)
}


make_plot_dt <- function(x_data, hclust_obj, proj_dt, clusters_dt) {
  # plot
  plot_dt 			= data.table(melt(as.matrix(x_data), varnames=c('cell_id', 'symbol')))
  
  # nice ordering
  symbol_order 		= colnames(x_data)[hclust_obj$order]
  plot_dt[, symbol 	:= factor(symbol, levels=symbol_order)]
  plot_dt[, cell_id 	:= factor(cell_id, levels=proj_dt$cell_id)]
  
  # put this into plotting 
  plot_dt 			= clusters_dt[ plot_dt, on='symbol' ]
  plot_dt 			= proj_dt[ plot_dt, on='cell_id' ]
  
  return(plot_dt)
}


psup_go <- function(psuper_obj, org_mapping = "org.Hs.eg.db", k = 5, sig_cutoff = 1) 
{
  if (!requireNamespace("topGO", quietly = TRUE)) {
    message("topGO not installed; not doing GO analysis")
    return()
  }
  if (!requireNamespace("fastcluster", quietly = TRUE)) {
    message("fastcluster not installed; not doing GO analysis")
    return()
  }
  glmnet_best = psuper_obj$glmnet_best
  best_lambdas = psuper_obj$best_lambdas
  proj_dt = copy(psuper_obj$proj_dt)
  x_data = copy(psuper_obj$x_data)
  beta_dt = psuper_obj$beta_dt
  cuts_dt = psuper_obj$cuts_dt
  rownames(x_data) = sprintf("cell_%04d", 1:nrow(x_data))
  set(proj_dt, i = NULL, "cell_id", rownames(x_data))
  setorder(proj_dt, psuper)
  message("clustering genes")
  hclust_obj = fastcluster::hclust(dist(t(x_data)), method = "complete")
  clusters_dt = calc_clusters_dt(hclust_obj, x_data, proj_dt, 
                                  k)
  go_results = do_topgo_for_cluster(clusters_dt, sig_cutoff, 
                                     org_mapping)
  plot_dt = make_plot_dt(x_data, hclust_obj, proj_dt, clusters_dt)
  go_list = list(clusters_dt = clusters_dt, 
                 go_results = go_results, 
                 plot_dt = plot_dt, 
                 cuts_dt = copy(psuper_obj$cuts_dt))
  
  # order: https://github.com/wmacnair/psupertime/issues/16 ##credits!!
  x <- go_list$plot_dt
  x_ord <-  with(x, x[order(label_psuper),])
  x_ord$cell_id <- factor(x_ord$cell_id, levels = unique(x_ord$cell_id))
  
  go_list$plot_dt <- x_ord
  return(go_list)
}

# define function for input formatting - produce gchromVAR compatible GWAS summary input
loadGWASsum <- function(path, format, pvalue = 0.05){
      table <- read_delim(file = path, 
                        "\t", escape_double = FALSE, trim_ws = F) %>% 
      dplyr::select('CHR_ID', 'CHR_POS', 'P-VALUE','REGION','OR or BETA','PVALUE_MLOG','MAPPED_GENE') %>%
      dplyr::filter(`P-VALUE` < pvalue) %>%
      separate(col = "MAPPED_GENE", sep = " - ",into = c("gene1","gene2")) %>%
        separate(col = "gene1", sep = ", ", into = paste0("gene",1:4)) %>% 
        pivot_longer(cols = paste0("gene",1:4),names_to = "gene_type",
                     values_to = "symbol", values_drop_na = T) %>% 
        dplyr::select(-gene_type) %>%
        unique() 
      
  return(table)
}

