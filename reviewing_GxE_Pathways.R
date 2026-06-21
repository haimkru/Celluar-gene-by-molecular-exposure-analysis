###reviewing_GxE_Pathways###
#author: Haim Krupkin
#date: 01/16/2026
#path: /oak/stanford/groups/smontgom/hkrupkin/Code/reviewing_GxE_Pathways
#description: This script is meant to load GxE pathway hits from the script /oak/stanford/groups/smontgom/hkrupkin/Code/Examining_GenePathWays_From_GxE_Analysis.R
#and examine them

###Packages###
library(stringr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)
###Analysis####
#loading methylation data
meth_pathways<-fread("/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/Outputs/meth_gene_enrichment_hits.tsv")
unique(meth_pathways$eigenFeature_ome)
unique(meth_pathways$cell_type_processed)
meth_pathways%>%
  dplyr::filter(p.adjust<0.05 & Count>3)%>%
  View()

#loading expression data
expression_pathways<-fread("/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/Outputs/Expression_GO_Enrichemnts_for_all_by_genes_from_Each_cell_type_eigenFeature_ome.tsv")
unique(expression_pathways$eigenFeature_ome)
unique(expression_pathways$cell_type_processed)
expression_pathways%>%
  dplyr::filter(p.adjust<0.05)%>%# & Count>5)%>%
  View()
