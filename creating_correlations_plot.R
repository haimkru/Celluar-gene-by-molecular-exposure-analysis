#creating_correlations_plot
#Source of code: is mainly almost 100% marie.
#/labs/smontgom/grps_smontgom/hkrupkin/GxE/Code/creating_correlations_plot
#author: Him Krupkin
#date: 05/17/2026
#descrption: this is to create the correlations plot


.libPaths(c("/scg/apps/software/r/4.3.3/lib", .libPaths()))

.libPaths("/home/mahuynh/R/x86_64-pc-linux-gnu-library/4.2")

library(WGCNA)
library(ComplexHeatmap)
library(isva)
library(magrittr)
library(tidyverse)
options(stringsAsFactors = FALSE);
library(readxl)
library(ggplot2)
setwd("/labs/smontgom/grps_smontgom/mahuynh/metabolites")

##################################
# We define the timepoint/datatime
##################################
datatime = 1

if (datatime == 1){
  ######################
  # Read Exam for metabolites
  ######################
  datatime = 1
  
  datatype = "metabolites"
  norm_method = "VST" # Options : VST, ranknorm
  what_to_regress <- "batch_and_ancestry" # Options : batch, batch_and_ancestry, batch_ancestry_sex_age
  num_rows = 541 
  
  if (datatype == "metabolites"){
    if (datatime == 1){
      num_modules = 52
      softPower = 8
      deep_split = 2
      MEDissThres = 0.15
      minModuleSize = 20
    }else{
      num_modules = 52
      softPower = 8
      deep_split = 2
      MEDissThres = 0.15
      minModuleSize = 20
    }
    data_path = paste("MESA/data/adjusted_data/Final_Version_Exam_", datatime,"_", norm_method, "_", datatype, "_", num_rows, "samples_", what_to_regress, ".rds",sep="")
    datExpr <- readRDS(file = data_path)
    base_dir = paste0("MESA/data/WGCNA/Final_Exam_", datatime, "_", datatype, "_norm_", norm_method, "_samples_", num_rows)
    path = paste0(base_dir, "/", 
                  num_modules, "modules_", "soft_thres_", softPower, "_deep_split_", deep_split, "_me_", MEDissThres, 
                  "_min_size_", minModuleSize, ".RData", sep = "")
    load(path)
  } else {
  }
  
  rownames <- c("ID", rownames(datExpr))
  metadf <- read.csv(paste0("/oak/stanford/groups/smontgom/mahuynh/metabolites/MESA/data/raw_data/Metabolites_Exam_", datatime,".txt") ,sep="\t")
  metadf <- metadf[, rownames]
  
  
  MEs1 <- MEs
  moduleColors1 <- moduleColors
  moduleLabels1 <- moduleLabels
  
  ######################
  # Read Exam for proteins 
  ######################
  data_path="/labs/smontgom/grps_smontgom/mahuynh/metabolites/MESA/data/adjusted_data/Final_Version_Exam_1_VST_proteins_541samples_batch_and_ancestry.rds"
  data_path
  datExpr <- readRDS(file = data_path)
  path = paste0("/labs/smontgom/grps_smontgom/mahuynh/metabolites/MESA/data/WGCNA/Final_Exam_1_proteins_norm_VST_samples_541", "/","18modules_soft_thres_7_deep_split_3_me_0.25_min_size_20_4_removed_pcs.RData",sep="")
  path
  load(path)
  
  MEs2 <- MEs
  moduleColors2 <- moduleColors
  moduleLabels2 <- moduleLabels
}

# Define number of samples
nSamples <- nrow(datExpr)
print(nSamples)

# Correlate Eigenfeatures
MEsCor = cor(MEs1, MEs2, use = "p")
cat("Dimensions of MEsCor:\n")
cat(dim(MEsCor), "\n")

# Calculate the p-values
MEsPValue = corPvalueStudent(MEsCor, nSamples)
cat("Dimensions of MEsPValue:\n")
cat(dim(MEsPValue), "\n")
cat(min(MEsPValue))



# Assuming MEsCor is your object
dimensions <- dim(MEsCor)
cat("Dimensions of MEsCor:", dimensions, "\n")

cat("Number of Modules:", length(names(MEs1)))
cat("Number of Modules:", length(names(MEs2)))

rownames(MEsCor) <- names(MEs1)
colnames(MEsCor) <- names(MEs2)

MEsCor


saveRDS(MEsCor,
        "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Code/heatmaping_habibi.rds")

# Remove "ME" from row names and column names
rownames(MEsCor) <- gsub("^ME", "", rownames(MEsCor))
colnames(MEsCor) <- gsub("^ME", "", colnames(MEsCor))

# Also do the same for MEsPValue since it's used in cell_fun
rownames(MEsPValue) <- gsub("^ME", "", rownames(MEsPValue))
colnames(MEsPValue) <- gsub("^ME", "", colnames(MEsPValue))


png("/oak/stanford/groups/smontgom/hkrupkin/GxE/plots/plots_for_figures/heatmap_plot_correlation_between_MEMs.png",
    width = 7, height = 6, units = "in", res = 300)

Heatmap(
  MEsCor,
  name = "Pearson Correlation",
  column_title = "Proteomic MEMs",
  row_title = "Metabolomic MEMs",
  column_title_gp = gpar(fontsize = 14),
  row_title_gp = gpar(fontsize = 14),
  row_names_gp = gpar(fontsize = 7), 
  column_names_gp = gpar(fontsize = 8),
  row_names_side = "left",
  row_dend_side = "right",
  column_names_side = "top",
  column_dend_side = "bottom", 
  cell_fun = function(j, i, x, y, width, height, fill) {
    p_value <- MEsPValue[i, j]
    if (p_value < 0.1) {
      grid.text(sprintf("%.2f", p_value), x, y, gp = gpar(fontsize = 4, col = "white"))
    }
  }
)


dev.off()

###################################################################



