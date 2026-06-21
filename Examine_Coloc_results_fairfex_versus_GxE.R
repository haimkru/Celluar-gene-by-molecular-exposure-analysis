###Examine_Coloc_results_fairfex_versus_GxE####
#/oak/stanford/groups/smontgom/hkrupkin/Code/Examine_Coloc_results_fairfex_versus_GxE
#author: Haim Krupkin
#date: 03/17/2026

###description:#####
#this script is meant to examine the results from the colonization between FairFex and GxE

###packages######
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(patchwork)
#################

###analysis#####
colocs<-fread("/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/coloc_results/MEtan/MEtan_coloc_results.csv")
colocs%>%
  dplyr::filter(PP.H4>0.5)%>%
  View()
####
colocs<-fread("/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/coloc_results_big/MEtan_coloc_results.csv")
#colocs%>%
#  dplyr::filter(PP.H4>0.5)%>%
#  View()
###Explaining_results######

colocs%>%
  dplyr::filter(PP.H4>0.5)%>%
  View()

classified_colocs_sig<-colocs%>%
  dplyr::filter(PP.H4>0.5)%>%
  group_by(lead_snp, signal) %>%
  mutate(
    has_naive    = any(condition == "Naive"),
    has_response = any(condition != "Naive")
  ) %>%
  ungroup() %>%
  mutate(
    variant_class = case_when(
      has_naive & has_response  ~ "naive_and_response",
      has_naive & !has_response ~ "naive_only",
      !has_naive & has_response ~ "response_only",
      .default = "unclassified"
    ),
    classification = case_when(
      signal == "interaction" & variant_class == "naive_and_response" ~ "interaction_naive_and_response",
      signal == "interaction" & variant_class == "naive_only"         ~ "interaction_naive_only",
      signal == "interaction" & variant_class == "response_only"      ~ "interaction_response_only",
      signal != "interaction" & variant_class == "naive_and_response" ~ "main_effect_naive_and_response",
      signal != "interaction" & variant_class == "naive_only"         ~ "main_effect_naive_only",
      signal != "interaction" & variant_class == "response_only"      ~ "main_effect_response_only",
      .default = "unclassified"
    )
  )

###lets_load_all_of_them####
coloc_file_paths <- list.files(
  "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/coloc_results_big/",
  pattern = ".csv",
  full.names = TRUE
)

coloc_file_paths <- setdiff(
  coloc_file_paths,
  "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/coloc_results_big//MEtan_coloc_results.csv"
)

# Get file info (includes modification time)
file_info <- file.info(coloc_file_paths)

# Sort by modification time (newest → oldest)
files <- coloc_file_paths[order(file_info$mtime, decreasing = TRUE)]

files

### read + add true path column ####
dfs <- lapply(files, function(f) {
  print(f)
  dt <- fread(f, header = TRUE)
  dt[, true_file_path := normalizePath(f, winslash = "/", mustWork = FALSE)]
  dt
})

combined <- rbindlist(dfs, use.names = TRUE, fill = TRUE)

###saving intermediate files####
write.csv(combined,
          "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/all_coloc_fairfex_GxE_results.csv")

###classifitng the colocs####
classified_combined<-combined%>%
  dplyr::filter(PP.H4>0.5)%>%
  group_by(true_file_path,lead_snp, signal) %>%
  mutate(
    has_naive    = any(condition == "Naive"),
    has_response = any(condition != "Naive")
  ) %>%
  ungroup() %>%
  mutate(
    variant_class = case_when(
      has_naive & has_response  ~ "naive_and_response",
      has_naive & !has_response ~ "naive_only",
      !has_naive & has_response ~ "response_only",
      .default = "unclassified"
    ),
    classification = case_when(
      signal == "interaction" & variant_class == "naive_and_response" ~ "interaction_naive_and_response",
      signal == "interaction" & variant_class == "naive_only"         ~ "interaction_naive_only",
      signal == "interaction" & variant_class == "response_only"      ~ "interaction_response_only",
      signal != "interaction" & variant_class == "naive_and_response" ~ "main_effect_naive_and_response",
      signal != "interaction" & variant_class == "naive_only"         ~ "main_effect_naive_only",
      signal != "interaction" & variant_class == "response_only"      ~ "main_effect_response_only",
      .default = "unclassified"
    )
  )
classified_combined%>%
  View()




classified_combined%>%
  dplyr::filter(classification=="interaction_response_only")%>%
  dplyr::filter(grepl("ENSG",phenotype_id))%>%
  pull(phenotype_id)%>%
  unique()%>%
  cat(sep = "\n")


###plotting the summary of the general trends of colonization with Fairfex####
#locuszoomr
library(locuscomparer)
#ENSG00000234439.1
#chr21_20701702_T_C_b38
#/oak/stanford/groups/smontgom/hkrupkin/GxE/Data/coloc_results_big/nk_expression_MEtan_trans_coloc_results.csv
#
relevant_nk_coloc<-fread("/oak/stanford/groups/smontgom/hkrupkin/GxE/Data/coloc_results_big/nk_expression_MEtan_trans_coloc_results.csv")
relevant_nk_GxE<-fread("/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/for_developments/output_bigs/nk_expression_MEtan_trans_qtl.tsv.gz")
fairfex_results<-fread("/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/coloc_results_big/tmp/fx_hg38.csv")
##
target_lead_snp <- "chr21_20701702_T_C_b38"
target_condition <- "LPS24"
fairfax_pval_col <- "LPS24.p.value" # Parameterized column name

# ---------------------------------------------------------
# 2. Extract Window Info from your Coloc Results
# ---------------------------------------------------------
# Get the specific row for this colocalization
coloc_hit <- relevant_nk_coloc[lead_snp == target_lead_snp & condition == target_condition]

target_gene_ensg <- coloc_hit$phenotype_id[1]
window_start <- coloc_hit$window_start[1]
window_end <- coloc_hit$window_end[1]

# ---------------------------------------------------------
# 3. Clean and Filter GxE Data
# ---------------------------------------------------------
# Filter for the target gene
gxe_sub <- relevant_nk_GxE[phenotype_id == target_gene_ensg]

# Extract position from variant_id (chr_pos_ref_alt_b38)
gxe_sub[, pos_hg38 := as.numeric(sapply(strsplit(variant_id, "_"), `[`, 2))]

# Filter to the genomic window
gxe_sub <- gxe_sub[pos_hg38 >= window_start & pos_hg38 <= window_end]

# ---------------------------------------------------------
# 4. Clean and Filter Fairfax Data
# ---------------------------------------------------------
# Filter Fairfax by hg38 position (Note: Fairfax uses HGNC symbols, GxE uses ENSG, 
# so merging by position is the safest and most accurate method).
fairfax_sub <- fairfex_results[pos_hg38 >= window_start & pos_hg38 <= window_end]

# Extract the dynamically chosen p-value column
fairfax_sub[, fairfax_p := get(fairfax_pval_col)]

# ---------------------------------------------------------
# 5. Merge Datasets
# ---------------------------------------------------------
# Merge by hg38 position
plot_data <- merge(
  gxe_sub[, .(pos_hg38, variant_id, pval_GxE = pval)], # Assuming 'pval' is the interaction p-value
  fairfax_sub[, .(pos_hg38, Gene, pval_Fairfax = fairfax_p)],
  by = "pos_hg38",
  all=T
)
plot_data
# Calculate -log10 P-values
plot_data[, logP_GxE := -log10(pval_GxE)]
plot_data[, logP_Fairfax := -log10(pval_Fairfax)]

# Identify the lead SNP for coloring
plot_data[, is_lead := ifelse(variant_id == target_lead_snp,
                              "Lead SNP", "Other SNPs")]

# ---------------------------------------------------------
# 6. Build the Plots
# ---------------------------------------------------------
# Color palette
snp_colors <- c("Lead SNP" = "red", "Other SNPs" = "grey50")
snp_sizes <- c("Lead SNP" = 3, "Other SNPs" = 1.5)

# Plot A: GxE Stack
p_gxe <- ggplot(plot_data, aes(x = pos_hg38 / 1e6, y = logP_GxE, color = is_lead, size = is_lead)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = snp_colors) +
  scale_size_manual(values = snp_sizes) +
  labs(title = paste("GxE eQTL:", target_gene_ensg), y = expression(-log[10](P[GxE]))) +
  theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank())

# Plot B: Fairfax Stack
p_fairfax <- plot_data%>%
  dplyr::select(logP_Fairfax,is_lead,pos_hg38)%>%
  ggplot( aes(x = pos_hg38 / 1e6, y = logP_Fairfax))+#, color = is_lead, size = is_lead)) +
  geom_point(alpha = 0.8)
  #scale_color_manual(values = snp_colors) +
  #scale_size_manual(values = snp_sizes) +
  #labs(title = paste("Fairfax QTL:", plot_data$Gene[1], paste0("(", target_condition, ")")), 
  #     x = "Position (hg38, Mb)", y = expression(-log[10](P[Fairfax]))) +
  #theme_classic() +
  #theme(legend.position = "none")
p_fairfax
# Plot C: LocusCompare Scatter (No R-squared)
p_scatter <- ggplot(plot_data, aes(x = logP_GxE, y = logP_Fairfax, color = is_lead, size = is_lead)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = snp_colors) +
  scale_size_manual(values = snp_sizes) +
  labs(title = "Colocalization Scatter", 
       x = expression(-log[10](P[GxE])), 
       y = expression(-log[10](P[Fairfax]))) +
  theme_classic() +
  theme(legend.position = "bottom", legend.title = element_blank())

# ---------------------------------------------------------
# 7. Combine using Patchwork
# ---------------------------------------------------------
# Layout: Stacked loci on the left, scatter plot on the right
final_plot <- (p_gxe / p_fairfax) | p_scatter
final_plot <- final_plot + plot_annotation(
  title = paste("Colocalization Locus:", target_lead_snp),
  subtitle = paste("Condition:", target_condition)
)

print(final_plot)



































