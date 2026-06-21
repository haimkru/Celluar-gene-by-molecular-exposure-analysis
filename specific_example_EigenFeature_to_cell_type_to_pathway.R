###specific_example_EigenFeature_to_cell_type_to_pathway####
#pathway: /oak/stanford/groups/smontgom/hkrupkin/GxE/Code/specific_example_EigenFeature_to_cell_type_to_pathway
#author: Haim Krupkin
#date: 01/30/2026
#updated on 05/01/2026

#description: the purpose of this file is to have an example GxE story
#pseudocode:
#1.we are going to : 1. identify an eigenFeature with a high vriance explain of a phenotype
#using the data from code: /oak/stanford/groups/smontgom/hkrupkin/GxE/Code/correlating_phenotypes_with_eigenFeatures_simple.R
#2. we are going to drill down on the gene pathways associated with it.
#using teh data from code:
#

###packages####
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(stringr)
library(data.table)

library(ggplot2)
library(ggpubr)
library(dplyr)
library(data.table)
library(readr)
library(stringr)

library(tidyr)
library(dynamicTreeCut)
library(colorspace)
library(dendextend )
library(ape)
library(broom)
library(clusterProfiler)
library(org.Hs.eg.db)  # For human genes; replace with the correct database for your organism
library(AnnotationDbi)
library(ggrepel)
library(pbmcapply)
library(parallel)
library(doParallel)


###function defintion###
run_go_enrichment <- function(df,
                              gene_col = "phenotype_id",
                              orgdb = org.Hs.eg.db,
                              ont = "BP",
                              pAdjustMethod = "BH",
                              pvalueCutoff = 1,
                              qvalueCutoff = 1) {
  
  library(dplyr)
  library(clusterProfiler)
  
  # 1. Clean gene IDs (remove version numbers like ENSG... .1)
  df <- df %>%
    mutate(phenotype_id_noversion = sub("\\..*$", "", .data[[gene_col]]))
  
  # 2. Extract unique gene list
  gene_list <- df %>%
    pull(phenotype_id_noversion) %>%
    unique()
  
  # 3. GO enrichment
  go_enrichment <- enrichGO(
    gene          = gene_list,
    OrgDb         = orgdb,
    keyType       = "ENSEMBL",
    ont           = ont,
    pAdjustMethod = pAdjustMethod,
    pvalueCutoff  = pvalueCutoff,
    qvalueCutoff  = qvalueCutoff
  )
  
  # 4. Extract results safely
  result <- as.data.frame(go_enrichment)
  
  # 5. Compute Fold Enrichment
  result$FoldEnrichment <- with(
    result,
    (as.numeric(sub("/.*", "", GeneRatio)) /
       as.numeric(sub(".*/", "", GeneRatio))) /
      (as.numeric(sub("/.*", "", BgRatio)) /
         as.numeric(sub(".*/", "", BgRatio)))
  )
  
  return(result)
}
###optimized for speed_version
run_go_enrichment_optimized <- function(gene_list,
                              gene_col = "phenotype_id",
                              orgdb = org.Hs.eg.db,
                              ont = "BP",
                              pAdjustMethod = "BH",
                              pvalueCutoff = 1,
                              qvalueCutoff = 1) {
  

  
  # 3. GO enrichment
  go_enrichment <- enrichGO(
    gene          = gene_list,
    OrgDb         = orgdb,
    keyType       = "ENSEMBL",
    ont           = ont,
    pAdjustMethod = pAdjustMethod,
    pvalueCutoff  = pvalueCutoff,
    qvalueCutoff  = qvalueCutoff
  )
  
  # 4. Extract results safely
  result <- as.data.frame(go_enrichment)
  
  # 5. Compute Fold Enrichment
  result$FoldEnrichment <- with(
    result,
    (as.numeric(sub("/.*", "", GeneRatio)) /
       as.numeric(sub(".*/", "", GeneRatio))) /
      (as.numeric(sub("/.*", "", BgRatio)) /
         as.numeric(sub(".*/", "", BgRatio)))
  )
  
  return(result)
}





###analysis####
#loading_eigenFeature to Phenotype Assocatiosn
term_table<-read.table("/oak/stanford/groups/smontgom/hkrupkin/GxE/Data/Outputs/single_lm_models_term_results.tsv")
#loading eigenFeature Genes-pathway associations
#expression_pathways<-fread("/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/Outputs/Expression_GO_Enrichemnts_for_all_by_genes_from_Each_cell_type_eigenFeature_ome.tsv",
#                                       sep="\t")

expression_pathways<-fread("/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/Outputs/Expression_GO_Enrichemnts_for_all_by_genes_from_Each_cell_type_eigenFeature_omeV3.tsv",
                           sep="\t")
####################################################

###followup on a specific example G versus GxE hits####
#loading all the hits####
all_hits<-fread("/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/all_cell_types_sig_hits.tsv")
all_hits$term<-str_replace(str_replace(all_hits$term,pattern = "ME",""),"_Eigen","")
all_hits$phenotype_id_noversion <- sub("\\..*$", "", all_hits$phenotype_id)

###here we define what are the GxE versus G we are comparing
#getting the GxE data
Expresion_hits_interaction<-all_hits%>%
  dplyr::filter(assay=="Expression")%>%
  dplyr::filter(effect_type=="interaction_effect")%>%
  dplyr::filter(cell_type=="CD4⁺ T" & term=="purple"  & Omics=="Proteome")
head(Expresion_hits_interaction)

#getting the G data
Expresion_hits_main_effect<-all_hits%>%
  dplyr::filter(assay=="Expression")%>%
  dplyr::filter(effect_type=="main_effect")%>%
  dplyr::filter(cell_type=="CD4⁺ T"  & Omics=="Proteome")
head(Expresion_hits_main_effect)

#running the GO ontology on both
main_pathways<-run_go_enrichment(Expresion_hits_main_effect)
interaction_pathways<-run_go_enrichment(Expresion_hits_interaction)

main_pathways$effect_type<-"main_effect"
interaction_pathways$effect_type<-"interaction_effect"
str(main_pathways)
str(interaction_pathways)

###plotting the differences###
main_and_interaction<-rbind(main_pathways,interaction_pathways)

####
main_and_interaction%>%
  as.data.frame() %>%
  dplyr::select(ID, Description, qvalue, effect_type) %>%
  tidyr::pivot_wider(names_from = effect_type, values_from = qvalue) %>%
  mutate(label_path=ifelse(main_effect<0.05 | interaction_effect<0.05,Description,NA ))%>%
  ggplot(aes(x = -log10(main_effect), y = -log10(interaction_effect))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
  theme_bw()+
  geom_text_repel(aes(label=label_path))+
  labs(x="GxE Genes GO Pathway Enrichment (-log10(p.adjusted))",y="GxE Genes GO Pathway Enrichment (-log10(p.adjusted))")
###testing in mono#
Expresion_hits_interaction<-all_hits%>%
  dplyr::filter(assay=="Expression")%>%
  dplyr::filter(effect_type=="interaction_effect")%>%
  dplyr::filter(cell_type=="Mono" & term=="purple"  & Omics=="Proteome")
head(Expresion_hits_interaction)

#getting the G data
Expresion_hits_main_effect<-all_hits%>%
  dplyr::filter(assay=="Expression")%>%
  dplyr::filter(effect_type=="main_effect")%>%
  dplyr::filter(cell_type=="Mono"  & Omics=="Proteome")
head(Expresion_hits_main_effect)

#running the GO ontology on both
main_pathways<-run_go_enrichment(Expresion_hits_main_effect)
interaction_pathways<-run_go_enrichment(Expresion_hits_interaction)

main_pathways$effect_type<-"main_effect"
interaction_pathways$effect_type<-"interaction_effect"
str(main_pathways)
str(interaction_pathways)

###plotting the differences###
main_and_interaction<-rbind(main_pathways,interaction_pathways)

####
main_and_interaction%>%
  as.data.frame() %>%
  dplyr::select(ID, Description, qvalue, effect_type) %>%
  tidyr::pivot_wider(names_from = effect_type, values_from = qvalue) %>%
  mutate(label_path=ifelse(main_effect<0.05 | interaction_effect<0.05,Description,NA ))%>%
  ggplot(aes(x = -log10(main_effect), y = -log10(interaction_effect))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
  theme_bw()+
  geom_text_repel(aes(label=label_path))+
  labs(x="GxE Genes GO Pathway Enrichment (-log10(p.adjusted))",y="GxE Genes GO Pathway Enrichment (-log10(p.adjusted))")




#######

### 1. Setup the Workload ####
# Filter to only Proteome and Expression up front
proteome_hits <- all_hits %>%
  dplyr::filter(assay == "Expression" & Omics == "Proteome")

# Identify all unique combinations of Cell Type and Module (term)
task_grid <- proteome_hits %>%
  dplyr::select(cell_type, term) %>%
  distinct()

print(paste("Processing", nrow(task_grid), "combinations in parallel..."))

### 2. Parallelized Execution ####
# Set number of cores (usually leave 1-2 free for the system)
num_cores <- 4

# Path for log file
log_file <- "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/Outputs/go_progress_log.txt"
writeLines(paste("Job started at:", Sys.time()), log_file)


cat("Found", nrow(task_grid), "unique cell type & module combinations to process.\n")

### 4. Parallel Execution with Progress Bar ###
# pbmclapply works exactly like mclapply but gives you a live progress bar
results_list <- pbmclapply(1:nrow(task_grid), function(i) {
  
  current_cell <- task_grid$cell_type[i]
  current_term <- task_grid$term[i]
  #current_cell <- task_grid$cell_type[42]
  #current_term <- task_grid$term[42]
  
  

  # Get GxE genes
  genes_int <- proteome_hits %>%
    dplyr::filter(cell_type == current_cell, term == current_term, effect_type == "interaction_effect") %>%
    pull(phenotype_id_noversion)%>%
    unique()
  
  # Get Main Effect genes
  genes_main <- proteome_hits %>%
    dplyr::filter(cell_type == current_cell, effect_type == "main_effect") %>%
    pull(phenotype_id_noversion)%>%
    unique()
  
  # Run GO enrichment
  res_int <- run_go_enrichment_optimized(genes_int)
  res_main <- run_go_enrichment_optimized(genes_main)
  
  # Labeling results correctly
  res_int$effect_type <- "interaction_effect"
  res_int$cell_type <- current_cell
  res_int$term <- current_term
  res_main$effect_type <- "main_effect"
  res_main$cell_type <- current_cell
  res_main$term <- current_term
  
  return(rbind(res_int, res_main))
  
  }
  )
  
}, mc.cores = num_cores)

### 5. Combine and Save ###
# Bind all the lists together into one master data frame
master_results <- rbindlist(results_list, fill = TRUE)
master_results%>%
  dplyr::filter(term=="purple" & Omics=="Proteome")%>%
  View()
master_results%>%
  dplyr::filter(term=="yellow" & Omics=="Proteome")%>%
  View()

master_results%>%
  dplyr::filter(term=="tan" & Omics=="Metabolome")%>%
  View()
###########################
  
