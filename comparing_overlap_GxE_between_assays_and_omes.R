###comparing_overlap_GxE_between_assays_and_omes####
#/labs/smontgom/grps_smontgom/hkrupkin/GxE/Code/comparing_overlap_GxE_between_assays_and_omes
#author: Haim Krupkin
#date: 05/11/2026

##Description:
#ok things we are going to do is:
#for the same cell type, we are going to check for overlaps In GxE hits by assay and by proteome,
#I am trying to think whats the best way to visualize this
#cause there are many cell types.
#options include 
#upset plot one cell type as an example
#
.libPaths(c("/scg/apps/software/r/4.3.3/lib", .libPaths()))
###packages#########
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(stringr)
#install.packages("UpSetR",lib = "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Code/")
library(UpSetR,lib.loc="/labs/smontgom/grps_smontgom/hkrupkin/GxE/Code/")

###packages####
library(ggplot2)
library(ggpubr)
library(dplyr)
library(data.table)
library(readr)
library(stringr)
library(dplyr)
library(ggplot2)
library(stringr)
library(stringr)
library(cowplot)
library(tidyverse)
library(purrr)

install.packages("data.table","/labs/smontgom/grps_smontgom/hkrupkin/GxE/Code/")
library(data.table)


###to get gene symbols for hubs#####



# 1. Define the GTF path
gtf_path <- "/oak/stanford/groups/smontgom/dnachun/projects/topmed/gencode.v39.annotation.gtf.gz"

# 2. Use a system command to extract the mapping (FAST)
# We only care about lines where the 3rd column is 'gene'
# We then extract gene_id and gene_name from the attributes
cmd <- paste0("zcat ", gtf_path, " | awk '$3 == \"gene\" {print $0}' | sed -E 's/.*gene_id \"([^ \"]+)\".*gene_name \"([^ \"]+)\".*/\\1\\t\\2/'")

# 3. Read into R
gene_map <- fread(cmd, col.names = c("gene_id", "gene_name"), header = FALSE)

# 4. Remove duplicates (GTF often has many entries per ID)
gene_map <- unique(gene_map)

# 5. Fast Lookup function
get_symbol <- function(ensembl_id, map_dt) {
  # This works even if your ID has a version (e.g., .1) and the map doesn't, or vice-versa
  # Strip version for matching if necessary
  symbol <- map_dt[gene_id == ensembl_id, gene_name]
  
  if (length(symbol) == 0) {
    # Try matching without version number
    clean_id <- gsub("\\..*", "", ensembl_id)
    symbol <- map_dt[gsub("\\..*", "", gene_id) == clean_id, gene_name]
  }
  
  return(if(length(symbol) > 0) symbol[1] else NA)
}
#####


options(scipen = 999)

all_hits<-fread("/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/all_cell_types_sig_hits.tsv")

#head(all_hits)
all_hits$unique_terms<-paste0(all_hits$term,all_hits$Omics,"_")

gxe_hubs <- all_hits %>%
  filter(effect_type == "interaction_effect") %>%
  group_by(phenotype_id, cell_type, Omics, assay) %>%
  summarise(
    n_unique_envs = n_distinct((unique_terms)),
    env_list = paste(unique((unique_terms)), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_unique_envs))

thresholds_by_assay <- gxe_hubs %>%
  group_by(assay) %>%
  summarise(
    threshold = quantile(n_unique_envs, 0.99),
    .groups = "drop"
  )

print(thresholds_by_assay)

# 4. Join thresholds back and define hubs
gxe_hubs <- gxe_hubs %>%
  left_join(thresholds_by_assay, by = "assay") %>%
  mutate(is_hub = n_unique_envs >= threshold)

# 5. Merge with gene map
expression_hubs <- merge(
  gxe_hubs, 
  gene_map, 
  by.x = "phenotype_id", 
  by.y = "gene_id", 
  all.x = TRUE
)

# 6. Plot with per-assay threshold lines
plot_number_of_terms_per_phenotype <- expression_hubs %>%
  ggplot(aes(x = n_unique_envs)) +
  geom_histogram(binwidth = 1, fill = "#005088", color = "white") +
  geom_vline(
    data = thresholds_by_assay,
    aes(xintercept = threshold),
    linetype = "dashed", color = "red"
  ) +
  geom_text(
    data = thresholds_by_assay,
    aes(x = threshold + 0.5, y = Inf, label = paste0("Top 1% = ", round(threshold, 1))),
    color = "red", hjust = 0, vjust = 1.5, size = 3
  ) +
  scale_y_continuous(transform = "log10") +
  theme_bw() +
  labs(
    x = "N Associated MEMs",
    y = "N Cell Type-Molecular Phenotype Pairs"
  ) +
  facet_wrap(~assay, ncol = 1, scales = "free_x")

ggsave(plot_number_of_terms_per_phenotype,
       width = 4,
       height = 4.5,
       dpi = 900,
       file = "/labs/smontgom/grps_smontgom/hkrupkin/GxE/plots/plots_for_figures/plot_number_of_terms_per_phenotype.png")


###global reactivitiy####
all_hits$unique_terms<-paste0(all_hits$term,all_hits$Omics,"_")
global_reactivity <- all_hits[effect_type == "interaction_effect" & cell_type!="Bulk", .(
  n_envs       = n_distinct(unique_terms),            # Total unique environmental modules
  n_cell_types = n_distinct(cell_type),       # Breadth across tissues/cells
  n_modalities = n_distinct(Omics),           # Sensitivity to both Prot and Metab
  cell_list    = paste(unique(cell_type), collapse = ", "),
  top_env      = names(which.max(table(term))) # The environment it likes most
), by = .(phenotype_id,assay)]

####



####




how_universal_is_p<-global_reactivity%>%
  ggplot(aes(x=n_envs))+
  geom_histogram(binwidth = 1)+
  theme_bw()+
  facet_wrap(~assay,scale="free")+
  labs(x="Number of Modules",y="Number of GxE Hits")
ggsave(plot=how_universal_is_p,
       filename="",
       width=5,
       height=3,
       dpi=900)

# 2. The Fixed Plot
plot_cell_breadth <- ggplot(global_reactivity, aes(x = n_cell_types, y = n_envs)) +
  
  # Keep the boxes grouped by n_cell_types
  geom_boxplot(aes(group = n_cell_types, fill = assay), 
               alpha = 0.7, outlier.shape = NA) +
  # Spearman Correlation (Numeric X-axis)
  stat_cor(aes(label = paste(after_stat(r.label), 
                             ifelse(after_stat(p) < 0.005, "p < 0.005", after_stat(p.label)), 
                             sep = "~`, `~")),
           method = "spearman", 
           label.x.npc = "left", 
           label.y.npc = "top", 
           size = 3) +
  facet_wrap(~assay, scales = "free_y") +
  
  # Formatting
  scale_fill_manual(values = c("Expression" = "#005088", "Methylation" = "#880050")) +
  scale_x_continuous(breaks = sort(unique(global_reactivity$n_cell_types))) +
  
  labs(
    x = "Number of Cell Types with GxE",
    y = "N Modules",
  ) +
  theme_bw()+
  theme(legend.position = "none")

print(plot_cell_breadth)
ggsave(plot_cell_breadth,
       width=4,
       height=3,
       dpi=900,
       file="/labs/smontgom/grps_smontgom/hkrupkin/GxE/plots/plots_for_figures/number_of_cell_types_to_unique_modules.png")
###
library(data.table)
library(dplyr)


# 1. Aggregate at the Gene + Assay + Omics level
global_reactivity_raw <- all_hits[effect_type == "interaction_effect" & cell_type != "Bulk", .(
  n_envs_by_ome = n_distinct(unique_terms), 
  n_cell_types  = n_distinct(cell_type)      
), by = .(phenotype_id, assay, Omics)]

# 2. Reshape wide, but include 'assay' in the grouping so it isn't collapsed
global_reactivity_wide <- dcast(
  global_reactivity_raw, 
  phenotype_id + assay ~ Omics, # Keeps assay separate
  value.var = c("n_envs_by_ome", "n_cell_types"),
  fill = 0 
)

# 3. Final data cleaning for the plot
global_reactivity_plot <- global_reactivity_wide %>%
  mutate(
    n_eigen_proteins    = n_envs_by_ome_Proteome,
    n_eigen_metabolites = n_envs_by_ome_Metabolome,
    max_cell_breadth    = pmax(n_cell_types_Proteome, n_cell_types_Metabolome)
  )

plot_ome_crosstalk_faceted <- global_reactivity_plot%>%
  mutate(max_cell_breadth=as.factor(max_cell_breadth))%>%
  ggplot( 
    aes(x = n_eigen_proteins, 
        y = n_eigen_metabolites, 
        color = max_cell_breadth)) + 
  
  # Jitter to prevent integer points from stacking
  geom_jitter(alpha = 0.7, size = 2.5, width = 0.2, height = 0.2) +
  
  # Trendlines calculated independently for Expression and Methylation
  geom_smooth(
    aes(group = factor(max_cell_breadth)),
    method = "lm",
    se = FALSE
  ) +  
  # Color gradient mapping cellular breadth
  scale_color_viridis_c(option = "plasma", 
                        name = "N Cell Types\nwith GxE",
                        breaks = sort(unique(global_reactivity_plot$max_cell_breadth))) +
  
  # Facet by the upstream regulatory layer (Expression vs Methylation)
  facet_wrap(~assay) +
  
  labs(
    x = "Number of Unique Eigen-Protein Envs",
    y = "Number of Unique Eigen-Metabolite Envs",
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.justification = "center"
  )+
  stat_cor(method="spearman")

print(plot_ome_crosstalk_faceted)


####making sure this is true statistically###

library(data.table)

# 1. Split the data by assay layer
expression_data  <- global_reactivity_plot[assay == "Expression"]
methylation_data <- global_reactivity_plot[assay == "Methylation"]

# 2. Run the GLM for the Expression Layer (Dynamic Cascade)
# We include an interaction term (*) to see if cellular breadth amplifies the protein-metabolite link
expression_model <- glm(n_eigen_metabolites ~ n_eigen_proteins * max_cell_breadth, 
                        data = expression_data, 
                        family = poisson(link = "log"))

cat("--- EXPRESSION MODEL RESULTS ---\n")
summary(expression_model)


# 3. Run the GLM for the Methylation Layer (Homeostatic Brake)
methylation_model <- glm(n_eigen_metabolites ~ n_eigen_proteins * max_cell_breadth, 
                         data = methylation_data, 
                         family = poisson(link = "log"))

cat("\n--- METHYLATION MODEL RESULTS ---\n")
summary(methylation_model)

# Bi-directional Proof for Expression
exp_model_flipped <- glm(n_eigen_proteins ~ n_eigen_metabolites * max_cell_breadth, 
                         data = expression_data, family = poisson(link = "log"))
summary(exp_model_flipped)
# Bi-directional Proof for Methylation
meth_model_flipped <- glm(n_eigen_proteins ~ n_eigen_metabolites * max_cell_breadth, 
                          data = methylation_data, family = poisson(link = "log"))
summary(meth_model_flipped)
###examining the exponential theory####
install.packages("ggeffects",lib="/labs/smontgom/grps_smontgom/hkrupkin/GxE/Code/")
library(ggeffects,lib.loc="/labs/smontgom/grps_smontgom/hkrupkin/GxE/Code/")
library(ggplot2)
library(patchwork) # To bind the two plots together side-by-side

# 1. Generate Exponential Predictions for Expression
# We look at cellular breadth across representative numbers of eigen-proteins (e.g., 1, 3, and 5)
pred_exp <- predict_response(expression_model, terms = c("max_cell_breadth", "n_eigen_proteins [1, 3, 5]"))

plot_exp <- ggplot(pred_exp, aes(x = x, y = predicted, color = group, fill = group)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, color = NA) +
  scale_color_viridis_d(option = "mako", end = 0.8, name = "Eigen-Proteins") +
  scale_fill_viridis_d(option = "mako", end = 0.8, name = "Eigen-Proteins") +
  labs(
    x = "Cellular Breadth (Number of Cell Types)",
    y = "Predicted Number of Eigen-Metabolites",
    title = "Expression: Multi-Omic Cascade",
    subtitle = "Cellular breadth multiplicatively amplifies cross-ome impact"
  ) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"), legend.position = "right")


# 2. Generate Exponential Predictions for Methylation
pred_meth <- predict_response(methylation_model, terms = c("max_cell_breadth", "n_eigen_proteins [1, 3, 5]"))

plot_meth <- ggplot(pred_meth, aes(x = x, y = predicted, color = group, fill = group)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, color = NA) +
  scale_color_viridis_d(option = "flare", end = 0.8, name = "Eigen-Proteins") +
  scale_fill_viridis_d(option = "flare", end = 0.8, name = "Eigen-Proteins") +
  labs(
    x = "Cellular Breadth (Number of Cell Types)",
    y = "Predicted Number of Eigen-Metabolites",
    title = "Methylation: The Threshold Breakthrough",
    subtitle = "High cellular breadth shatters the epigenetic constraint"
  ) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"), legend.position = "right")


# 3. Combine them into a single publication-ready figure
exponential_proof_plot <- plot_exp + plot_meth + 
  plot_layout(guides = "collect") + 
  plot_annotation(
    title = "Evidence of Exponential Multi-Omic Breakthroughs",
    subtitle = "Curves represent multiplicative marginal effects directly from the Poisson GLMs",
    theme = theme(plot.title = element_text(face = "bold", size = 16))
  )

print(exponential_proof_plot)

plot_ome_crosstalk_faceted <- global_reactivity_plot %>%
  mutate(max_cell_breadth = factor(max_cell_breadth, 
                                   levels = sort(unique(max_cell_breadth)))) %>%
  ggplot(aes(x = n_eigen_proteins,
             y = n_eigen_metabolites,
             color = max_cell_breadth,
             group = max_cell_breadth)) +
  geom_jitter(alpha = 0.7, size = 2.5, width = 0.2, height = 0.2) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_viridis_d(option = "plasma",
                        name = "N Cell Types\nwith GxE") +
  facet_wrap(~assay) +
  stat_cor(method = "spearman") +
  labs(x = "Number of Unique Eigen-Protein Envs",
       y = "Number of Unique Eigen-Metabolite Envs") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.justification = "center")

print(plot_ome_crosstalk_faceted)

####fully self contained interaction plot####
library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)

library(ggplot2)
library(dplyr)
library(patchwork)

# ==============================================================================
# 1. GENERATE THE TEXT ANNOTATIONS FROM YOUR GLMS
# ==============================================================================
get_clean_stats_text <- function(model) {
  s <- summary(model)$coefficients
  
  # Extract exact Estimates and P-values
  est_prot  <- s["n_eigen_proteins", "Estimate"]
  p_prot    <- s["n_eigen_proteins", "Pr(>|z|)"]
  
  est_cell  <- s["max_cell_breadth", "Estimate"]
  p_cell    <- s["max_cell_breadth", "Pr(>|z|)"]
  
  est_inter <- s["n_eigen_proteins:max_cell_breadth", "Estimate"]
  p_inter   <- s["n_eigen_proteins:max_cell_breadth", "Pr(>|z|)"]
  
  fmt_p <- function(p) if(p < 0.001) "p < 0.001" else paste0("p = ", round(p, 3))
  
  paste0(
    "Poisson GLM Effects (Log Scale):\n",
    "• n_eigen_proteins: β = ", round(est_prot, 3), " (", fmt_p(p_prot), ")\n",
    "• max_cell_breadth: β = ", round(est_cell, 3), " (", fmt_p(p_cell), ")\n",
    "• Interaction Term: β = ", round(est_inter, 3), " (", fmt_p(p_inter), ")"
  )
}

exp_annot_text  <- get_clean_stats_text(expression_model)
meth_annot_text <- get_clean_stats_text(methylation_model)

# Standardize the grouping factor
plot_data_clean <- global_reactivity_plot %>%
  mutate(max_cell_breadth = factor(max_cell_breadth, levels = sort(unique(max_cell_breadth))))

# ==============================================================================
# 2. PANEL A: EXPRESSION EMPIRICAL DISTRIBUTION + GLM STATS
# ==============================================================================
plot_final_exp <- plot_data_clean %>% filter(assay == "Expression") %>%
  ggplot(aes(x = n_eigen_proteins, y = n_eigen_metabolites, color = max_cell_breadth, group = max_cell_breadth)) +
  geom_jitter(alpha = 0.5, size = 2.2, width = 0.2, height = 0.2) +
  geom_smooth(method = "lm", se = T, linewidth = 1.0) +
  scale_color_viridis_d(option = "plasma", end = 0.85, name = "N Cell Types") +
  # Dynamically place the text box inside the plot panel
  annotate("text", x = -Inf, y = Inf, label = exp_annot_text, 
           hjust = -0.05, vjust = 1.15, size = 2, fontface = "plain", color = "black") +
  labs(x = "Number of EigenProtein", 
       y = "Number of EigenMetabolite") +
  theme_bw()

# ==============================================================================
# 3. PANEL B: METHYLATION EMPIRICAL DISTRIBUTION + GLM STATS
# ==============================================================================
plot_final_meth <- plot_data_clean %>% filter(assay == "Methylation") %>%
  ggplot(aes(x = n_eigen_proteins, y = n_eigen_metabolites, color = max_cell_breadth, group = max_cell_breadth)) +
  geom_jitter(alpha = 0.5, size = 2.2, width = 0.2, height = 0.2) +
  geom_smooth(method = "lm", se = T, linewidth = 1.0) +
  scale_color_viridis_d(option = "plasma", end = 0.85, name = "N Cell Types") +
  # Dynamically place the text box inside the plot panel
  annotate("text", x = -Inf, y = Inf, label = meth_annot_text, 
           hjust = -0.05, vjust = 1.15, size = 2, fontface = "plain", color = "black") +
  xlim(0,18)+
  ylim(0,52)+
  labs(x = "Number of EigenProtein", 
       y = "Number of EigenMetabolite") +
  theme_bw() 

# ==============================================================================
# 4. ASSEMBLE SIDE-BY-SIDE MASTER FIGURE
# ==============================================================================
master_two_panel_plot <- plot_final_exp + plot_final_meth + 
  plot_layout(guides = "collect") + 
  plot_annotation(
    theme = theme(
      plot.subtitle = element_text(italic = TRUE, size = 10, color = "gray20")
    )
  )

print(master_two_panel_plot)

ggsave(plot=master_two_panel_plot,
       file="/labs/smontgom/grps_smontgom/hkrupkin/GxE/plots/plots_for_figures/assocation_between_number_of_metabolites_to_proteins.png",
       width=7,
       height = 4,
       dpi=900)
###examining hub genes####
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

# 1. Define the ultra-strict 0.1% threshold
top_01_cutoff <- 0.999

# 2. Score and filter the absolute heavyweights
heavyweight_genes <- global_reactivity_plot %>%
  group_by(assay) %>% 
  mutate(
    rank_proteins    = percent_rank(n_eigen_proteins),
    rank_metabolites = percent_rank(n_eigen_metabolites),
    rank_breadth     = percent_rank(max_cell_breadth)
  ) %>%
  ungroup()

elite_cohort <- heavyweight_genes %>%
  filter(
    rank_proteins    >= top_01_cutoff,
    rank_metabolites >= top_01_cutoff,
    rank_breadth     >= top_01_cutoff
  ) %>%
  # Strip version dots (e.g., ENSG00000182261.5 -> ENSG00000182261)
  mutate(clean_ensembl = gsub("\\..*", "", phenotype_id))

# Extract clean unique gene vectors for the enrichment tool
expression_hub_ids  <- elite_cohort %>% filter(assay == "Expression")  %>% pull(clean_ensembl) %>% unique()
methylation_hub_ids <- elite_cohort %>% filter(assay == "Methylation") %>% pull(clean_ensembl) %>% unique()library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

# 1. Define the ultra-strict 0.1% threshold
top_01_cutoff <- 0.99

# 2. Score and filter the absolute heavyweights
heavyweight_genes <- global_reactivity_plot %>%
  group_by(assay) %>% 
  mutate(
    rank_proteins    = percent_rank(n_eigen_proteins),
    rank_metabolites = percent_rank(n_eigen_metabolites),
    rank_breadth     = percent_rank(max_cell_breadth)
  ) %>%
  ungroup()

elite_cohort <- heavyweight_genes %>%
  filter(
    rank_proteins    >= top_01_cutoff,
    rank_metabolites >= top_01_cutoff,
    rank_breadth     >= 0.6
  ) %>%
  # Strip version dots (e.g., ENSG00000182261.5 -> ENSG00000182261)
  mutate(clean_ensembl = gsub("\\..*", "", phenotype_id))

elite_cohort%>%
  merge(gene_map,
               by.x="phenotype_id",
               by.y="gene_id")%>%
  View()

# Extract clean unique gene vectors for the enrichment tool
expression_hub_ids  <- elite_cohort %>% filter(assay == "Expression")  %>% pull(clean_ensembl) %>% unique()
expression_hub_ids

cat(expression_hub_ids, sep = "\n")
# 3. Run enrichment for the Expression Master Hubs
go_enrich_exp <- enrichGO(
  gene          = expression_hub_ids,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP",           # Biological Process
  pAdjustMethod = "BH",           # Benjamini-Hochberg FDR correction
  pvalueCutoff  = 0.05,
  minGSSize = 3,
  qvalueCutoff  = 0.05,
  readable      = TRUE            # Converts Ensembl IDs back to readable Gene Symbols in the output
)

View(go_enrich_exp@result)
