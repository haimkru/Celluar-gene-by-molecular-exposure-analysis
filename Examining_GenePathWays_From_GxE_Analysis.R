###Examining_GenePathWays_From_GxE_Analysis####
#/oak/stanford/groups/smontgom/hkrupkin/Code/Examining_GenePathWays_From_GxE_Analysis
#author: Haim Krupkin
#date: 01/15/2026
#update date: 04/30/2026 - we corrected the GxE analysis to be better, more correct, using the correct refrence file and the correct column...
#description:This script is meant to take all the signficant hits from GxE.
#and looking if the genes from the gene expression or the regions from the methylations are all from specific pathways



###packages###
library(ggplot2)
library(ggpubr)
library(dplyr)
library(data.table)
library(readr)
library(stringr)

library(tidyr)
library(data.table)
library(stringr)
library(dynamicTreeCut)
library(colorspace)
library(dendextend )
library(ape)
library(broom)
library(clusterProfiler)
library(org.Hs.eg.db)  # For human genes; replace with the correct database for your organism
library(AnnotationDbi)
library(ggrepel)


#install.packages("ggvenn",lib="/labs/smontgom/grps_smontgom/hkrupkin/Code/packages")
###analysis####
#####loading data####
#
#all_hitsV2<-fread("/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/all_hits_methylation_and_expression_with_methylation_genesV2.tsv")

all_hits<-fread("/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/all_cell_types_sig_hits.tsv")



all_hits$term<-str_replace(str_replace(all_hits$term,pattern = "ME",""),"_Eigen","")
###Expression GSEA analysis#
#we are going to start with the GSEA for the expression first
#this was updated on 04/30/2026 because it did not have the effect type distinction.... orignally it was:
#Expresion_hits<-all_hits%>%dplyr::filter(assay=="Expression")
Expresion_hits<-all_hits%>%dplyr::filter(assay=="Expression")%>%
  dplyr::filter(effect_type=="interaction_effect")
dim(Expresion_hits)
dim(all_hits)
# we are going to do it sepretly per cell type and per eigenFeature
#testing is done initally on one;, but the loop can and should be runnable on all.
###example####
#this required an update on 04/30/2026, before we had:
#Expresion_hits<-Expresion_hits%>%dplyr::filter(cell_type=="Mono" & term=="lightgreen" & omics=="Proteome")
#now we have:
Expresion_hits<-Expresion_hits%>%dplyr::filter(cell_type=="Mono" & term=="lightgreen" & Omics=="Proteome")
dim(Expresion_hits)
Expresion_hits$phenotype_id_noversion <- sub("\\..*$", "", Expresion_hits$phenotype_id)

all_genes_associated_with_eigenFeature<-Expresion_hits%>%
  pull(phenotype_id_noversion)%>%unique()
all_genes_associated_with_eigenFeature
dim(all_genes_associated_with_eigenFeature)
go_enrichment <- enrichGO(gene = all_genes_associated_with_eigenFeature,
                                                  OrgDb        = org.Hs.eg.db,
                                                  keyType       = 'ENSEMBL',
                                                  ont          = "BP",      # Ontology: "BP" (Biological Process), "MF" (Molecular Function), or "CC" (Cellular Component)
                                                  pAdjustMethod = "BH",     # Adjust p-values using Benjamini-Hochberg method
                                                  pvalueCutoff  = 1,
                                                  qvalueCutoff  = 1)
result_go_enrichment<-go_enrichment@result
result_go_enrichment$FoldEnrichment <- with(
  result_go_enrichment,
  (as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio))) /
    (as.numeric(sub("/.*", "", BgRatio)) / as.numeric(sub(".*/", "", BgRatio)))
)
####
#result_go_enrichment%>%
#  dplyr::filter(p.adjust<0.05)%>%
#  View()

plot_Volvano_result_go_enrichment<-result_go_enrichment%>%
  mutate(labeling_bp=ifelse(p.adjust<0.05,Description,NA))%>%
  ggplot(aes(x=log2(FoldEnrichment),y=-log10(p.adjust)))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05),
             color="red",
             linetype = "dashed")+
  geom_text_repel(
    aes(label = labeling_bp ),
    force = 3,               # Increase the repulsion force to spread labels further apart
    box.padding = 0.5,       # Add padding around each label box
    point.padding = 0.3,     # Add padding between points and labels
    max.overlaps = 10,      # Allow as many labels as possible to be shown
    min.segment.length = 0,  # Ensures connecting lines are drawn even for close labels
    segment.size = 0.2,
    size=4,# Makes the line segment thinner for a cleaner look
    segment.color = "grey50" # Color of the connecting lines (optional)
  )+
  theme_minimal() +                    # Optional: gives a clean theme
  theme(plot.title = element_text(hjust = 0.5, size = 14))+
  labs(x = "log2(FoldEnrichment)",
       y = "-log10(Adjusted P-Value)",
       title="Nearest Genes GO Analysis Biological Process")+
  theme(
    axis.title = element_text(size = 24), # Axis titles
    axis.text = element_text(size = 24)   # Axis text
  )
plot_Volvano_result_go_enrichment
plot_Volvano_result_go_enrichment

ggsave(plot=plot_Volvano_result_go_enrichment,
       "/labs/smontgom/grps_smontgom/hkrupkin/GxE/plots/Example_plot_Volvano_result_go_enrichment_monocytes_lightgreen.png",
       width=6,
       height=6)




###running it on all the cell types and all the possible 
#this is the version prior to 04/30/2026:
#Expresion_hits<-all_hits%>%dplyr::filter(assay=="Expression")
#POST 04/30/2026 WE CORRECTED:
Expresion_hits<-all_hits%>%dplyr::filter(assay=="Expression")%>%dplyr::filter(effect_type=="interaction_effect")

Expresion_hits$phenotype_id_noversion <- sub("\\..*$", "", Expresion_hits$phenotype_id)
possible_cell_types<-unique(Expresion_hits$cell_type)
list_of_eigenFeatures<-unique(Expresion_hits$term)
possible_mixes<-length(possible_cell_types)*length(list_of_eigenFeatures)
print(paste0("we have : ",length(possible_cell_types), " cell types "))
Expresion_hits$ome_and_eigenFeature<-paste0(Expresion_hits$term,"_",Expresion_hits$Omics)
list_of_eigenFeature_ome<-unique(Expresion_hits$ome_and_eigenFeature)
print(paste0("we have : ", length(list_of_eigenFeature_ome) ," Possible list_of_eigenFeature_ome"))


possible_mixes<-length(possible_cell_types)*length(list_of_eigenFeature_ome)
print(possible_mixes)
what_mix_are_we_on=0
all_pathways_enrichment=data.frame()
hs_db <- org.Hs.eg.db
list_of_go_enrichemtns<-list()
index=0
for (cell_type_processed in possible_cell_types){
  print("we are working on cell type :")
  print(cell_type_processed)
  for (eigenFeature_ome in list_of_eigenFeature_ome){
    print("we are on eigenFeature : ")
    print(eigenFeature_ome)
    print("we are on mix number: ")
    print(what_mix_are_we_on)
    what_mix_are_we_on=what_mix_are_we_on+1
    
    all_genes_associated_with_eigenFeature<-Expresion_hits%>%
      dplyr::filter(cell_type==cell_type_processed & ome_and_eigenFeature==eigenFeature_ome)%>%
      pull(phenotype_id_noversion)%>%unique()
    
    go_enrichment <- enrichGO(gene = all_genes_associated_with_eigenFeature,
                              OrgDb        = hs_db,
                              keyType       = 'ENSEMBL',
                              ont          = "BP",      # Ontology: "BP" (Biological Process), "MF" (Molecular Function), or "CC" (Cellular Component)
                              pAdjustMethod = "BH",     # Adjust p-values using Benjamini-Hochberg method
                              pvalueCutoff  = 1,
                              qvalueCutoff  = 1)
    result_go_enrichment<-go_enrichment@result
    result_go_enrichment$FoldEnrichment <- with(
      result_go_enrichment,
      (as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio))) /
        (as.numeric(sub("/.*", "", BgRatio)) / as.numeric(sub(".*/", "", BgRatio)))
    )
    
    result_go_enrichment$eigenFeature_ome<-eigenFeature_ome
    result_go_enrichment$cell_type_processed<-cell_type_processed
    list_of_go_enrichemtns[[what_mix_are_we_on]]<-result_go_enrichment
    
    #the breaking is for testing if this will work :)
    #if (what_mix_are_we_on>3){
    #  double_break=T
    #  break
    #}
  }
  #if (double_break==T){
  #  break
  #}
}
all_pathways_enrichment <- rbindlist(list_of_go_enrichemtns, use.names = TRUE, fill = TRUE)

#saving the big file
#readr::write_tsv(all_pathways_enrichment,
 #                "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/Outputs/Expression_GO_Enrichemnts_for_all_by_genes_from_Each_cell_type_eigenFeature_omeV2.tsv")
readr::write_tsv(all_pathways_enrichment,
                 "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/Outputs/Expression_GO_Enrichemnts_for_all_by_genes_from_Each_cell_type_eigenFeature_omeV3.tsv")

###geting_pathway_enrichments_for_methylations#####
#PRE 04/30/2026
Meth_hits<-all_hits%>%dplyr::filter(assay=="Methylation")
#AFTER 04/30/2026
Meth_hits<-all_hits%>%dplyr::filter(assay=="Methylation")%>%dplyr::filter(effect_type=="interaction_effect")

Meth_hits$phenotype_id_noversion <- sub("\\..*$", "", Meth_hits$phenotype_id)
possible_cell_types<-unique(Meth_hits$cell_type)
list_of_eigenFeatures<-unique(Meth_hits$term)
possible_mixes<-length(possible_cell_types)*length(list_of_eigenFeatures)
print(paste0("we have : ",length(possible_cell_types), " cell types "))
Meth_hits$ome_and_eigenFeature<-paste0(Meth_hits$term,"_",Meth_hits$Omics)
list_of_eigenFeature_ome<-unique(Meth_hits$ome_and_eigenFeature)
print(paste0("we have : ", length(list_of_eigenFeature_ome) ," Possible list_of_eigenFeature_ome"))


possible_mixes<-length(possible_cell_types)*length(list_of_eigenFeature_ome)
print(possible_mixes)
what_mix_are_we_on=0
all_pathways_enrichment_met=data.frame()
hs_db <- org.Hs.eg.db

for (cell_type_processed in possible_cell_types){
  print("we are working on cell type :")
  print(cell_type_processed)
  for (eigenFeature_ome in list_of_eigenFeature_ome){
    print("we are on eigenFeature : ")
    print(eigenFeature_ome)
    print("we are on mix number: ")
    print(what_mix_are_we_on)
    what_mix_are_we_on=what_mix_are_we_on+1
    
    all_genes_associated_with_eigenFeature<-Meth_hits%>%
      #dplyr::filter(cell_type==cell_type_processed & ome_and_eigenFeature==eigenFeature_ome)%>%
      dplyr::filter(cell_type==cell_type_processed & ome_and_eigenFeature==eigenFeature_ome & annotation!="Intergenic")%>%
      pull(gene_id)%>%unique()
    
    go_enrichment <- enrichGO(gene = all_genes_associated_with_eigenFeature,
                              OrgDb        = hs_db,
                              keyType       = 'ENSEMBL',
                              ont          = "BP",      # Ontology: "BP" (Biological Process), "MF" (Molecular Function), or "CC" (Cellular Component)
                              pAdjustMethod = "BH",     # Adjust p-values using Benjamini-Hochberg method
                              pvalueCutoff  = 1,
                              qvalueCutoff  = 1)
    result_go_enrichment<-go_enrichment@result
    result_go_enrichment$FoldEnrichment <- with(
      result_go_enrichment,
      (as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio))) /
        (as.numeric(sub("/.*", "", BgRatio)) / as.numeric(sub(".*/", "", BgRatio)))
    )
    
    result_go_enrichment$eigenFeature_ome<-eigenFeature_ome
    result_go_enrichment$cell_type_processed<-cell_type_processed
    all_pathways_enrichment_met<-rbind(all_pathways_enrichment_met,result_go_enrichment)
    
    #the breaking is for testing if this will work :)
    #if (what_mix_are_we_on>3){
    #  double_break=T
    #  break
    #}
  }
  #if (double_break==T){
  #  break
  #}
}

#all_pathways_enrichment_met%>%
#  dplyr::filter(qvalue<0.05)%>%
#  View()

readr::write_tsv(all_pathways_enrichment_met,
                 "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/Outputs/meth_gene_enrichment_hits.tsv")

