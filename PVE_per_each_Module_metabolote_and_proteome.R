#PVE_per_each_Module_metabolote_and_proteome
#based on : /labs/smontgom/grps_smontgom/hkrupkin/GxE/Code/PVE_per_each_Module_metabolote_and_proteome
#/labs/smontgom/grps_smontgom/hkrupkin/GxE/Code/PVE_per_each_Module.R
#author: Haim Krupkin
#date: 04/26/2026

.libPaths(c("/scg/apps/software/r/4.3.3/lib", .libPaths()))
#.libPaths("/home/mahuynh/R/x86_64-pc-linux-gnu-library/4.1")
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
setwd("/labs/smontgom/grps_smontgom/mahuynh/metabolites/")
library("ggplot2")
library(dplyr)
library(gridExtra)


################################################################
# Load the previous parts
################################################################
#We load the expression and trait data saved in the first part


#Proteins_18modules_soft_thres_7_deep_split_3_me_0.25_min_size_20_4_removed_pcs.txt

data_path="/labs/smontgom/grps_smontgom/mahuynh/metabolites/MESA/data/adjusted_data/Final_Version_Exam_1_VST_proteins_541samples_batch_and_ancestry.rds"
data_path
datExpr <- readRDS(file = data_path)
path = paste0("/labs/smontgom/grps_smontgom/mahuynh/metabolites/MESA/data/WGCNA/Final_Exam_1_proteins_norm_VST_samples_541", "/","18modules_soft_thres_7_deep_split_3_me_0.25_min_size_20_4_removed_pcs.RData",sep="")
path
load(path)


################################################################
# What variance is explained by the metabolome and proteome?
################################################################
#getting total varaince epxlainiend across the entire ome
# 1. Total number of genes (including grey)
total_genes <- length(moduleColors)

# 2. Unique modules excluding grey
unique_mods <- setdiff(unique(moduleColors), "grey")

# 3. Sum of per-gene R² across all modules
sum_r2 <- sum(sapply(unique_mods, function(m) {
  
  # Indices of genes in this module
  mod_genes_idx <- which(moduleColors == m)
  
  # This module's eigengene
  me <- MEs[, paste0("ME", m)]
  
  # Per-gene R² against the module eigengene
  # (identical to what propVarExplained computes internally, 
  #  but retained at gene level rather than averaged)
  r2_per_gene <- cor(datExpr[, mod_genes_idx], me, use = "p")^2
  
  # Sum R² across all genes in this module
  sum(r2_per_gene)
  
}))

# 4. Divide by total genes (including grey) to get total PVE
total_pve <- sum_r2 / total_genes

cat("Total variance explained by all modules:", round(total_pve * 100, 2), "%\n")
###############################################################################
# First, we want to do a histogram of PVE

PVE <- propVarExplained(datExpr, moduleColors, MEs)
#We sort the PVE
sortedPVE <- sort(PVE, decreasing = TRUE)

print(mean(sortedPVE))
#We compute the number of genes per eigenmodule
colorCounts <- table(moduleColors)

print(length(colorCounts))
print(length(sortedPVE))

# First, create a data frame with your data
df <- data.frame(Module = names(sortedPVE), 
                 PVE = unlist(sortedPVE), 
                 ModuleSize = colorCounts[order(PVE, decreasing = TRUE)])

df <- dplyr::rename(df, ModuleColors = ModuleSize.moduleColors, Size = ModuleSize.Freq)

print(mean(df$PVE, na.rm = TRUE))

saveRDS(df,"/labs/smontgom/grps_smontgom/hkrupkin/GxE/plots/plots_for_figures/PVE_Exam1_protein_Colors.rds")
df<-readRDS("/labs/smontgom/grps_smontgom/hkrupkin/GxE/plots/plots_for_figures/PVE_Exam1_protein_Colors.rds")
png(paste0("/labs/smontgom/grps_smontgom/hkrupkin/GxE/plots/plots_for_figures/PVE_Exam1_protein_Colors.png"),
    width = 8, height = 5, units = "in", res = 900)
p1 <- ggplot(df, aes(x = reorder(ModuleColors, -PVE), y = PVE, fill = Size)) +
  geom_bar(stat = "identity") +
  theme_bw()+
  # Add trans = "log10" for a logarithmic scale
  scale_fill_gradient(low = "lightblue", high = "darkblue", trans = "log10") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Proteomic Module") +
  ylab("PVE in Each Module") +
  labs(fill = "Module size") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = margin(1, 1, 0.2, 1, "cm"))

p1

dev.off()

###combined proteomic and metabolomic####
df_proteome<-readRDS("/labs/smontgom/grps_smontgom/hkrupkin/GxE/plots/plots_for_figures/PVE_Exam1_protein_Colors.rds")
df_proteome$ome<-"Proteome"
df_metabolome<-readRDS("/labs/smontgom/grps_smontgom/hkrupkin/GxE/plots/plots_for_figures/PVE_Exam1_Metabolite_Colors.rds")
df_metabolome$ome<-"Metabolome"
both_omes<-rbind(df_proteome,df_metabolome)
p_both <- ggplot(both_omes, aes(x = reorder(ModuleColors, -PVE), y = PVE, fill = Size)) +
  geom_bar(stat = "identity") +
  theme_bw()+
  facet_wrap(~ome,ncol=1,scales="free_x")+
  # Add trans = "log10" for a logarithmic scale
  scale_fill_gradient(low = "lightblue", high = "darkblue", trans = "log10") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Module") +
  ylab("PVE in Each Module") +
  labs(fill = "Module size") +
  theme(
    strip.text = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 12, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    plot.title = element_text(size = 18, hjust = 0.5)
  )

p_both
png(paste0("/labs/smontgom/grps_smontgom/hkrupkin/GxE/plots/plots_for_figures/PVE_Exam1_Metabolome_and_proteome.png"),
    width = 8, height = 8, units = "in", res = 900)
p_both
dev.off()

###V2 renaming y axis name and addint label of total variance explained.
# Create annotation data frame
both_omes <- both_omes %>%
  group_by(ome) %>%
  arrange(desc(PVE)) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  mutate(unique_id = paste(ome, ModuleColors, sep = "___")) %>%
  mutate(unique_id = factor(unique_id, levels = unique_id[order(ome, rank)]))

# Create annotation data frame
annotations <- data.frame(
  ome = c("Metabolome", "Proteome"),
  label = c("Total PVE = 37.1%", "Total PVE = 10%")
)

p_both <- ggplot(both_omes, aes(x = unique_id, y = PVE, fill = Size)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  facet_wrap(~ome, ncol = 1, scales = "free") +
  scale_x_discrete(labels = function(x) gsub(".*___", "", x)) +
  scale_fill_gradient(low = "lightblue", high = "darkblue", trans = "log10") +
  geom_text(
    data = annotations,
    aes(x = Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = 1.1,
    vjust = 1.5,
    size = 5,
    fontface = "bold"
  ) +
  coord_cartesian(clip = "off") +
  xlab("Module") +
  ylab("PVE within Each Module") +
  labs(fill = "Module size") +
  theme(
    strip.text = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 12, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    plot.title = element_text(size = 18, hjust = 0.5)
  )

p_both
png(paste0("/labs/smontgom/grps_smontgom/hkrupkin/GxE/plots/plots_for_figures/PVE_Exam1_Metabolome_and_proteomeV2.png"),
    width = 8, height = 8, units = "in", res = 900)
p_both
dev.off()



