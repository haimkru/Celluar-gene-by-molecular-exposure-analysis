####Manually_checkig_gencode_versions####
#author: Haim Krupkin
#date: 12/12/2025
#description: this code is meant for manually checking gencode versions.

#packages###
library(dplyr)
library(rtracklayer)
####

#b cell
file_path_expression_standardized<-"/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/expression/bcell_expression_exam1_standardized.tsv"

#bulk RNA seq
#file_path_expression_standardized<-"/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/expression/bulk_expression_exam1_standardized.tsv"
#number of genes:

expression_standardized<-read.table(header=T,file_path_expression_standardized)
gencode_versions<-list.files(full.names = T,
           path="/labs/smontgom/grps_smontgom/hkrupkin/GxE/Code/gencode_files")
gencode_versions
overlap_df=data.frame()
for (file in rev(gencode_versions[1:10])){
  #file<-gencode_versions[1]
  print("working on file : ")
  print(file)
  gencode_file<-readGFF(file)
  gencode_genes <- filter(gencode_file, type == "gene")
  gencode_genes_positions <- dplyr::select(gencode_genes, seqid, start, end, gene_id)
  number_of_diffrences<-setdiff(expression_standardized$Name,gencode_genes_positions$gene_id)
  print("number of diffrences is : ")
  diffrneces=length(number_of_diffrences)
  print(diffrneces)
  print("number of genes originally: ")
  number_of_genes_originally=length(expression_standardized$Name)
  print(number_of_genes_originally)
  percent_match=(number_of_genes_originally-diffrneces)/number_of_genes_originally
  print(percent_match)
  overlap_df_for_this_gencode=data.frame(number_of_genes_originally=number_of_genes_originally,
                                         diffrneces=diffrneces,
                                         percent_match=percent_match,
                                         file=file)
  overlap_df=rbind(overlap_df,overlap_df_for_this_gencode)
}





