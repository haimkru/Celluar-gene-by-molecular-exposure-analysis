###Creating_TensorQTL_Appropriate_FilesV2###
#author: Haim Krupkin
#date : 12/11/2025
#description : This R script is meant to create all the files required to run tensorQTL for GxE analysis on MESA data 

###Pseudocode####
#this codes has part a which is for creation of normalized expression files and part b which is for the covariate files
#part a
#1. we are gonna create a function, that corrects the expresion standardized values and saves them in the correct format
#2. then we are gonna loop through the function on each file per cell type
#3. Then we aare gonna save the correct format TensorQTL data

#part b
#1. we load the genotype pc
#2. we then have a loop , that then will go through each cell type and create a specific 

#part c
#we create compatiable interaction files of the eigen metabolites and eigen proteins
###notes####

#!!!!!!!!!!THERE IS AN ISSUE WITH MAPPING REGIONS, I SUSPECT GENCODE VERSIONS!!!!!!!!!!!#

# 4 inputs for tensorQTL
# - .bed/.bim/.fam plink files for genotypes - already generated, 1 file for all runs - done
# - .bed for standardized methylation with gene/CpG positions - 6 files for expresion, 7 for methylation - do now expression
# - .txt (tab-separated) file with expression/methylation PCs, genotype PCs 
#-  .txt eigenmetabolites/eigenproteins 
# for methylation, need to
# - QC/filter out probes for variable sites
# - Standardize
# - Run PCA
####relies on which code files#####
#1. /oak/stanford/groups/smontgom/dnachun/projects/topmed/subset_data.R

###packages###
library(dplyr)
library(tidyverse)
library(rtracklayer)
library(data.table)
library(readr)
library(data.table)

###Analysis ####
###Parameters####
input_directory=""
path_output_directory="/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/Ready_for_running_inputs"



###making output directories###
path_expression_outputs=paste0(path_output_directory,
                               "/standardized_phenotypesV2/")
print("making directory if needed for expression:")
print(path_expression_outputs)
dir.create(path_expression_outputs)

path_saving_covariates<-paste0(path_output_directory,"/covariatesV2/")
print('we are saving the covriate files to : ')
print(path_saving_covariates)
dir.create(path_saving_covariates)

####making standardized expression fie with gene positions####
#####loading files that are generally needed to the analysis#####
######################         
#loading the gene code positosn as positions are required by tensorQTL
###
gene_code_positions_gencode39_easy<-read.table(header=T,
                                               "/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/expression/gencode_genes_positions.tsv")
gencode_genes_positions<-gene_code_positions_gencode39_easy

###
#gencode_coordinates <- readGFF("/oak/stanford/groups/smontgom/dnachun/projects/topmed/gencode.v39.annotation.gtf.gz")
#gencode_genes <- filter(gencode_coordinates, type == "gene")
#gencode_genes_positions <- select(gencode_genes, seqid, start, end, gene_id)

#Might require later fixing, due to incosistent gencode version.
#currently we simply remove the version from the ensembl id and merge like that...
#gencode_genes_positions<-gencode_genes_positions %>%
#mutate(gene_id_version_stripped = str_remove(gene_id, "\\.\\d+$"))%>%
#dplyr::select(-gene_id)
######################         

#loading the sample key so that we will have a universal id that is the same and will be used by tensorQTL, we decided to use the genomic one.
ids <- readr::read_lines("/oak/stanford/groups/smontgom/dnachun/projects/topmed/intersection_samples.txt")
sample_key <- read_tsv("/oak/stanford/groups/smontgom/shared/topmed/MESA_may2019/MESA_TOPMed_WideID_20190517.txt")
#this is so that we have a df we can use to reference and correct the expression values#
key_identifier_genotype_to_RNA<-sample_key%>%
  dplyr::filter(nwdid1 %in%ids )%>%
  dplyr::select(nwdid1,tor_id11 )


##oranigse everything according to this order.
sorted_ids_key<-sort(key_identifier_genotype_to_RNA$nwdid1)


#a function to have the correct colnames changed.
rename_expression_columns <- function(expr_df, key_df,
                                      from_col = "tor_id11",
                                      to_col   = "nwdid1") {
  if (!all(c(from_col, to_col) %in% colnames(key_df))) {
    stop("Key dataframe must contain columns: ", from_col, " and ", to_col)
  }
  
  id_map <- setNames(key_df[[to_col]], key_df[[from_col]])
  
  # Rename only columns that appear in the map
  colnames(expr_df) <- ifelse(
    colnames(expr_df) %in% names(id_map),
    id_map[colnames(expr_df)],
    colnames(expr_df)
  )
  
  return(expr_df)
}
sort_df_by_rownames <- function(df, row_order_vector) {
  # Ensure the row names exist in df
  row_order_vector <- intersect(row_order_vector, rownames(df))
  
  # Reorder rows
  df_ordered <- df[row_order_vector, , drop = FALSE]
  
  return(df_ordered)
}

###here we create a tensorQTL compatiable format for the exxpression data per cell type.

create_tensorQTL_compatiable_standarzied_expression_file<-function(file_path_expression_standardized){
  #if we are developing code, its benefitial to put this out of note and the row below it into a note.
  #file_path_expression_standardized<-per_cell_type_expression_standardized_files[1]
  
  file_path_expression_standardized<-file_path_expression_standardized #this is importatn to enable running of multiple files
  
  
  print("working on file: ")
  print(file_path_expression_standardized)
  cell_type<-gsub("_expression_exam1_standardized.tsv","",basename(file_path_expression_standardized))
  print("corresponding to cell type : ")
  print(cell_type)
  
  expression_standardized<-fread(header=T,sep="\t",file_path_expression_standardized)
  #here we change to correct column ids for genotype data#
  expression_standardized_genotype_ided <-rename_expression_columns(expression_standardized, key_identifier_genotype_to_RNA,
                                                                    from_col = "tor_id11",
                                                                    to_col   = "nwdid1")
  ########################
  #Might require later fixing, due to incosistent gencode version.
  #currently we simply remove the version from the ensembl id and merge like that...
  #expression_standardized_genotype_ided<-expression_standardized_genotype_ided%>%
  #  mutate(Name_version_stripped = str_remove(Name, "\\.\\d+$"))%>%
  #  dplyr::select(-Name)
  
  #we are still missign 252 genes, we need to discuss this with Dan
  setdiff(expression_standardized_genotype_ided$Name,
          gencode_genes_positions$gene_id)
  #merging
  genecode39_regions_with_expression_standardized_genotyped_ided<-gencode_genes_positions%>%
    merge(expression_standardized_genotype_ided,
          by.x="gene_id",
          by.y="Name")
  
  tensorQTL_format_expression_standardied<-genecode39_regions_with_expression_standardized_genotyped_ided%>%
    dplyr::mutate(phenotype_id=gene_id)%>%
    mutate(chr = gsub("chr", "", seqid))%>%
    dplyr::select(chr,start,end,phenotype_id,everything())
  
  tensorQTL_format_expression_standardied<-tensorQTL_format_expression_standardied%>%
    dplyr::select(-seqid)
  ######################
  
  file_path_output_for_tensor_QTL_dataframes<-paste0(path_expression_outputs,
                                                     cell_type,
                                                     ".bed")
  print("we are outputing to path:")
  print(file_path_output_for_tensor_QTL_dataframes)
  tensorQTL_format_expression_standardied<-tensorQTL_format_expression_standardied%>%dplyr::select(-gene_id)
  
  #backUp_column_names<-colnames(tensorQTL_format_expression_standardied)
  colnames(tensorQTL_format_expression_standardied)[1] <- "#chr"
  print(colnames(tensorQTL_format_expression_standardied))
  cols_to_order <- grep("^NWD", names(tensorQTL_format_expression_standardied), value = TRUE)
  cols_to_order
  # Keep other columns in original order
  other_cols <- setdiff(names(tensorQTL_format_expression_standardied), cols_to_order)
  other_cols
  # Reorder df
  tensorQTL_format_expression_standardied_ordered <- tensorQTL_format_expression_standardied[, c(other_cols, intersect(sorted_ids_key, cols_to_order))]
  chr_levels <- c(as.character(1:22), "X", "Y", "M")
  
  tensorQTL_format_expression_standardied_ordered <- tensorQTL_format_expression_standardied_ordered %>%
    mutate(`#chr` = factor(`#chr`, levels = chr_levels)) %>%
    arrange(`#chr`, start, end)
  print(tensorQTL_format_expression_standardied_ordered[1:10,1:10])
  
  print(file_path_output_for_tensor_QTL_dataframes)
  if (basename(file_path_output_for_tensor_QTL_dataframes)=="b_cell.bed"){
    print("this is b cell we are correct the cell type")
    file_path_output_for_tensor_QTL_dataframes=gsub("b_cell","bcell",file_path_output_for_tensor_QTL_dataframes)
  }
  file_path_output_for_tensor_QTL_dataframes=gsub(".bed","_expression.bed",file_path_output_for_tensor_QTL_dataframes)
  print(file_path_output_for_tensor_QTL_dataframes)
  
  tensorQTL_format_expression_standardied_ordered<-tensorQTL_format_expression_standardied_ordered%>%
    dplyr::filter(`#chr`!="M")
  tensorQTL_format_expression_standardied_ordered<-tensorQTL_format_expression_standardied_ordered%>%
    dplyr::filter(`#chr`!="Y")
  print("the chromsomes we have are : ")
  print(table(tensorQTL_format_expression_standardied_ordered$`#chr`))
  
  readr::write_tsv(tensorQTL_format_expression_standardied_ordered,
                   file_path_output_for_tensor_QTL_dataframes)
  
}

#loading all the files of expression per cell type
per_cell_type_expression_standardized_files<-list.files(path="/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/expression/",
                                                        pattern="*standardized.tsv",
                                                        full.names = T)
print(per_cell_type_expression_standardized_files)
#running loop to run this per file
file_number=0
for (expression_file_path in per_cell_type_expression_standardized_files){
  print("we are in file number : ")
  print(file_number)
  file_number=file_number+1
  print("out of total")
  print(length(per_cell_type_expression_standardized_files))
  print("starting function on file : ")
  print(expression_file_path)
  create_tensorQTL_compatiable_standarzied_expression_file(expression_file_path)
  
}

s###Covariate Files####
####loading genotype data####
genotype_pcs<-read_table("/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/genotypes/ancestry_pcs.tsv")

#this is the function to the all the merging#
create_tensorQTL_compatiable_covirate_file<-function(file_path_expression_pcs_per_cell_type){
  
  #extracting the correct file to save to later
  #what_eigen_analysis<-"Metabolite"
  #file_path_expression_pcs_per_cell_type<-"/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/expression/bcell_pcs.tsv"
  #file_path_eigen_value<-"/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/TensorQTL_related/Metabolite_52modules_soft_thres_8_deep_split_2_me_0.15_min_size_20.txt"
  #we need to get cell type
  #file_path_expression_pcs_per_cell_type<-"/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/expression//bulk_exam1_pcs.tsv"
  cell_type_covariates<-gsub("_pcs.tsv","",basename(file_path_expression_pcs_per_cell_type))
  cell_type_covariates
  print("corresponding to cell type : ")
  print(cell_type_covariates)
  
  
  #making expression data mergeable#
  expression_pcs <- fread(header=T,sep="\t",
                          file_path_expression_pcs_per_cell_type
  )
  expression_pcs_with_genotype_id<-expression_pcs%>%
    merge(key_identifier_genotype_to_RNA,
          by.x="torid11",
          by.y="tor_id11")
  expression_pcs_with_genotype_id<-expression_pcs_with_genotype_id%>%
    dplyr::select(nwdid1 ,everything())%>%
    dplyr::select(-torid11)
  colnames(expression_pcs_with_genotype_id)
  
  colnames(expression_pcs_with_genotype_id)[2:length(colnames(expression_pcs_with_genotype_id))]<-paste0(colnames(expression_pcs_with_genotype_id)[2:length(colnames(expression_pcs_with_genotype_id))],
                                                                                                         "_Expression_PCs")
  colnames(expression_pcs_with_genotype_id)
  #making eigenvalue mergable
  expression_pcs_with_genotype_id
  
  #final merges
  data_frame_covariates<-merge(genotype_pcs,
                               expression_pcs_with_genotype_id,
                               by.x="nwdid",
                               by.y="nwdid1")
  
  #saving the file
  
  print("we are saving the covirates file to : ")
  if (cell_type_covariates=="bulk_exam1"){
    cell_type_covariates="bulk"
  }
  path_saving_covirates<-paste0(path_saving_covariates,
                                cell_type_covariates,
                                "_expression_covariates.txt")
  print(path_saving_covirates)
  
  
  data_frame_covariates_ordered<-data_frame_covariates%>%
    mutate(nwdid =  factor(nwdid  , levels = sorted_ids_key)) %>%
    arrange(nwdid  )
  data_frame_covariates_ordered<-data_frame_covariates_ordered%>%dplyr::mutate(sample=nwdid)%>%dplyr::select(-nwdid)
  data_frame_covariates_ordered<-data_frame_covariates_ordered%>%dplyr::select(sample,everything())
  print("number of na samples:")
  print(sum(is.na(data_frame_covariates_ordered$sample)))
  
  print(data_frame_covariates_ordered[1:10,])
  
  ids<-data_frame_covariates_ordered$sample
  transposed_covariates<-as.data.frame(t(data_frame_covariates_ordered)[-1,])
  colnames(transposed_covariates)<-ids
  print(transposed_covariates[1:10,])
  transposed_covariates$variable<-row.names(transposed_covariates) 
  transposed_covariates<-transposed_covariates%>%dplyr::select(variable,everything())
  readr::write_tsv(transposed_covariates,
                   path_saving_covirates)
  
}

#here we are going to run over each file

list_of_expression_pc_files<-list.files(path="/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/expression/",
                                        pattern="*pcs.tsv",
                                        full.names = T)
list_of_expression_pc_files


#running loop to run this per file
file_number=0
for (expression_pcs_file in list_of_expression_pc_files){
  print("we are in file number : ")
  print(file_number)
  file_number=file_number+1
  print("out of total")
  print(length(list_of_expression_pc_files))
  print("starting function on file : ")
  print(expression_pcs_file)
  #this is for the protein eignvalues
  print('saving merged protein')
  create_tensorQTL_compatiable_covirate_file(expression_pcs_file)
  
  
}

###writing eigen values sepretly###
#writing protein data
expression_eigen_dir="/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/Ready_for_running_inputs/interaction_termsV2"
dir.create(expression_eigen_dir)


eigen_value<-read.table("/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/TensorQTL_related/Proteins_18modules_soft_thres_7_deep_split_3_me_0.25_min_size_20_4_removed_pcs.txt",
                        header=T)
eigen_value$sample<-rownames(eigen_value)
#we are now going to merge the expression with the eigen value, then with the pcs, then we are going to save
colnames(eigen_value)[1:length(colnames(eigen_value))-1]<-paste0(colnames(eigen_value)[1:length(colnames(eigen_value))-1],
                                                                 "_Eigen_",
                                                                 "Protein")
eigen_value<-eigen_value%>%dplyr::select(sample,everything())
eigen_value<-eigen_value%>%
  mutate(sample = factor(sample   , levels = sorted_ids_key)) %>%
  arrange(sample   )
print(eigen_value[1:10,1:10])
readr::write_tsv(eigen_value,
                 paste0(expression_eigen_dir,"/","Eigen_Proteins_Expression_interaction.txt"))

#writing metabolite data
eigen_value<-read.table("/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/TensorQTL_related/52modules_soft_thres_8_deep_split_2_me_0.15_min_size_20.txt",
                        header=T)
eigen_value$sample<-rownames(eigen_value)
#we are now going to merge the expression with the eigen value, then with the pcs, then we are going to save
colnames(eigen_value)[1:length(colnames(eigen_value))-1]<-paste0(colnames(eigen_value)[1:length(colnames(eigen_value))-1],
                                                                 "_Eigen_",
                                                                 "Metabolite")
eigen_value<-eigen_value%>%dplyr::select(sample,everything())
eigen_value<-eigen_value%>%
  mutate(sample = factor(sample   , levels = sorted_ids_key)) %>%
  arrange(sample   )
print(eigen_value[1:10,1:10])

colnames(eigen_value)
colnames(eigen_value[,1:19])
colnames(eigen_value[,c(1,20:37)])
colnames(eigen_value[,c(1,38:53)])

readr::write_tsv(eigen_value,
                 paste0(expression_eigen_dir,"/","Eigen_Metabolite_Expression_interaction.txt"))
readr::write_tsv(eigen_value[,1:19],
                 paste0(expression_eigen_dir,"/","Eigen_Metabolite_Expression_interaction_chunk_1.txt"))
readr::write_tsv(eigen_value[,c(1,20:37)],
                 paste0(expression_eigen_dir,"/","Eigen_Metabolite_Expression_interaction_chunk_2.txt"))
readr::write_tsv(eigen_value[,c(1,38:53)],
                 paste0(expression_eigen_dir,"/","Eigen_Metabolite_Expression_interaction_chunk3.txt"))



###creating methyhlation covariates#####
###getting key identifier for methylation pcs####

#this is so that we have a df we can use to reference and correct the expression values#
ids <- readr::read_lines("/oak/stanford/groups/smontgom/dnachun/projects/topmed/intersection_samples.txt")


key_identifier_genotype_to_methylation<-sample_key%>%
  dplyr::filter(nwdid1 %in%ids )%>%
  dplyr::select(nwdid1,toe_id11)
print(key_identifier_genotype_to_methylation)

path_methylation_outputs=paste0(path_output_directory,
                                "/covariatesV2/")
print("making directory if needed for expression:")
print(path_methylation_outputs)
#dir.create(path_methylationn_outputs)

create_tensorQTL_compatiable_covirate_file_methlation<-function(methylation_pcs_file){
  #methylation_pcs_file<-list_of_methylation_pc_files[1]
  #print("woring on file : ")
  #print(methylation_pcs_file)
  #extracting the correct file to save to later
  #we need to get cell type
  cell_type_covariates_methylation<-gsub("_methylation_pcs.tsv","",basename(methylation_pcs_file))
  cell_type_covariates_methylation
  print("corresponding to cell type : ")
  print(cell_type_covariates_methylation)
  
  
  #making expression data mergeable#
  methylation_pcs <- fread(header=T,sep="\t",
                           methylation_pcs_file
  )
  methylation_pcs_with_genotype_id<-methylation_pcs%>%
    merge(key_identifier_genotype_to_methylation,
          by.x="toeid11",
          by.y="toe_id11")
  methylation_pcs_with_genotype_id<-methylation_pcs_with_genotype_id%>%
    dplyr::select(nwdid1 ,everything())%>%
    dplyr::select(-toeid11)
  colnames(methylation_pcs_with_genotype_id)
  
  colnames(methylation_pcs_with_genotype_id)[2:length(colnames(methylation_pcs_with_genotype_id))]<-paste0(colnames(methylation_pcs_with_genotype_id)[2:length(colnames(methylation_pcs_with_genotype_id))],
                                                                                                           "_Methylation_PCs")
  colnames(methylation_pcs_with_genotype_id)
  #making eigenvalue mergable
  
  
  #final merges
  data_frame_covariates_methylation<-merge(genotype_pcs,
                                           methylation_pcs_with_genotype_id,
                                           by.x="nwdid",
                                           by.y="nwdid1")
  
  #saving the file
  
  print("we are saving the covirates file to : ")
  path_methylation_outputs="/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/Ready_for_running_inputs/covariatesV2"
  path_saving_covirates_methylation<-paste0(path_methylation_outputs,
                                            "/",
                                            cell_type_covariates_methylation,
                                            "_methylation_covariates.txt")
  


  
  
  data_frame_covariates_methylation<-data_frame_covariates_methylation%>%
    mutate(nwdid =  factor(nwdid  , levels = sorted_ids_key)) %>%
    arrange(nwdid  )
  ###
  
  print(data_frame_covariates_methylation[1:10,1:10])
  print("we are saving to : ")
  print(path_saving_covirates_methylation)
  print("dim of methaltion covariates")
  print(dim(data_frame_covariates_methylation))
  
  print(data_frame_covariates_methylation[1:10,])
  
  ids<-data_frame_covariates_methylation$nwdid 
  transposed_covariates<-as.data.frame(t(data_frame_covariates_methylation)[-1,])
  colnames(transposed_covariates)<-ids
  print(transposed_covariates[1:10,])
  transposed_covariates$variable<-row.names(transposed_covariates) 
  transposed_covariates<-transposed_covariates%>%dplyr::select(variable,everything())
  print(transposed_covariates[1:10,])
  
  readr::write_tsv(transposed_covariates,
                   path_saving_covirates_methylation)
}

methylation_pcs_path<-"/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/methylation"

list_of_methylation_pc_files<-list.files(path=methylation_pcs_path,
                                         pattern="*_methylation_pcs.tsv",
                                         full.names = T)
list_of_methylation_pc_files


#running loop to run this per file
file_number_methylation=0
for (methylation_pcs_file in list_of_methylation_pc_files){
  print("we are in file number : ")
  print(file_number_methylation)
  file_number_methylation=file_number_methylation+1
  print("out of total")
  print(length(list_of_methylation_pc_files))
  print("starting function on file : ")
  print(methylation_pcs_file)
  #this is for the protein eignvalues
  print('saving merged protein')
  create_tensorQTL_compatiable_covirate_file_methlation(methylation_pcs_file)
  
  
}

###saving_standardized_methylation####
path_methylation_outputs="/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/Ready_for_running_inputs/standardized_phenotypesV2"
print("we are creating the path : ")
print(path_methylation_outputs)
dir.create(path_methylation_outputs)



rename_expression_columns <- function(expr_df, key_df,
                                      from_col = "tor_id11",
                                      to_col   = "nwdid1") {
  if (!all(c(from_col, to_col) %in% colnames(key_df))) {
    stop("Key dataframe must contain columns: ", from_col, " and ", to_col)
  }
  
  id_map <- setNames(key_df[[to_col]], key_df[[from_col]])
  
  # Rename only columns that appear in the map
  colnames(expr_df) <- ifelse(
    colnames(expr_df) %in% names(id_map),
    id_map[colnames(expr_df)],
    colnames(expr_df)
  )
  
  return(expr_df)
}
positions_probes_path<-"/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/methylation/epic_manifest.tsv"
positions_probes<-read_delim(positions_probes_path)


create_tensorQTL_compatiable_standarzied_methylation_file<-function(file_path_methylation_standardized){
  
  #file_path_methylation_standardized<-per_cell_type_methylation_standardized_files[2] #this is importatn to enable running of multiple files
  
  
  print("working on file: ")
  print(file_path_methylation_standardized)
  cell_type<-gsub("_norm.tsv","",basename(file_path_methylation_standardized))
  print("corresponding to cell type : ")
  print(cell_type)
  if (cell_type=="beta_var"){
    cell_type="bulk"
  }
  cell_type<-gsub("tensor_methylation_","",basename(cell_type))
  print("after adjusting the cell type name the cell type is : ")
  print(cell_type)
  
  print("loading file")
  methylation_standardized<-fread(header=T,sep="\t",file_path_methylation_standardized)
  #here we change to correct column ids for genotype data#
  print("renaming columns")
  methylation_standardized_genotype_ided <-rename_expression_columns(methylation_standardized,
                                                                     key_identifier_genotype_to_methylation,
                                                                     from_col = "toe_id11",
                                                                     to_col   = "nwdid1")
  
  setdiff(methylation_standardized_genotype_ided$Probe_ID,
          positions_probes$Probe_ID)
  #merging
  print("merging with probe positions")
  
  probe_regions_standardized_genotyped_ided<-positions_probes%>%
    merge(methylation_standardized_genotype_ided,
          by.x="Probe_ID",
          by.y="Probe_ID")
  print("making correct columns for tensor QTL compatiability")
  
  tensorQTL_format_methylation_standardied<-probe_regions_standardized_genotyped_ided%>%
    dplyr::mutate(phenotype_id=Probe_ID)%>%
    mutate(chr = gsub("chr", "", seqnames ))%>%
    dplyr::select(-Probe_ID,-width ,-strand,-nextBase,-seqnames)%>%
    dplyr::select(chr,start,end,phenotype_id,everything())
  
  ######################
  path_methylation_outputs="/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/Ready_for_running_inputs/standardized_phenotypesV2"
  file_path_output_for_tensor_QTL_dataframes_methylation<-paste0(path_methylation_outputs,"/",
                                                                 cell_type,
                                                                 ".bed")
  print("we are outputing to path:")
  print(file_path_output_for_tensor_QTL_dataframes_methylation)
  print(tensorQTL_format_methylation_standardied[1:10,1:10])
  print(file_path_output_for_tensor_QTL_dataframes_methylation)
  #print("Saving this big file....")
  
  
  ##

  #backUp_column_names<-colnames(tensorQTL_format_expression_standardied)
  colnames(tensorQTL_format_methylation_standardied)[1] <- "#chr"
  print(colnames(tensorQTL_format_methylation_standardied))
  cols_to_order <- grep("^NWD", names(tensorQTL_format_methylation_standardied), value = TRUE)
  cols_to_order
  # Keep other columns in original order
  other_cols <- setdiff(names(tensorQTL_format_methylation_standardied), cols_to_order)
  other_cols
  # Reorder df
  tensorQTL_format_methylation_standardied <- tensorQTL_format_methylation_standardied[, c(other_cols, intersect(sorted_ids_key, cols_to_order))]
  colnames(tensorQTL_format_methylation_standardied)
  chr_levels <- c(as.character(1:22), "X", "Y", "M")
  
  tensorQTL_format_methylation_standardied_ordered <- tensorQTL_format_methylation_standardied %>%
    mutate(`#chr` = factor(`#chr`, levels = chr_levels)) %>%
    arrange(`#chr`, start, end)
  
  tensorQTL_format_methylation_standardied_ordered<-tensorQTL_format_methylation_standardied_ordered%>%
    dplyr::filter(`#chr`!="M")
  tensorQTL_format_methylation_standardied_ordered<-tensorQTL_format_methylation_standardied_ordered%>%
    dplyr::filter(`#chr`!="Y")
  print(tensorQTL_format_methylation_standardied_ordered[1:10,1:10])
  print("the chromsomes we have are : ")
  print(table(tensorQTL_format_methylation_standardied_ordered$`#chr`))

  file_path_output_for_tensor_QTL_dataframes_methylation=gsub(".bed","_methylation.bed",file_path_output_for_tensor_QTL_dataframes_methylation)
  print(file_path_output_for_tensor_QTL_dataframes_methylation)
  
  readr::write_tsv(tensorQTL_format_methylation_standardied_ordered,
                   file_path_output_for_tensor_QTL_dataframes_methylation)
  print("finished for this file : ) ")

  print("-----------------------------------------------------------------------")
}

#loading all the files of expression per cell type
per_cell_type_methylation_standardized_files<-list.files(path="/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/methylation",
                                                         pattern="*_norm.tsv",
                                                         full.names = T)
print(per_cell_type_methylation_standardized_files)

#running loop to run this per file
file_number_methylation=0
for (methylation_standardized_file_path in per_cell_type_methylation_standardized_files){
  print("we are in file number : ")
  print(file_number_methylation)
  file_number_methylation=file_number_methylation+1
  print("out of total")
  print(length(per_cell_type_methylation_standardized_files))
  print("starting function on file : ")
  print(methylation_standardized_file_path)
  create_tensorQTL_compatiable_standarzied_methylation_file(methylation_standardized_file_path)
  
}


#saving_eigen_value_for_methylation_number_of_samples
methylation_eigen_dir="/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/Ready_for_running_inputs/interaction_termsV2"
dir.create(methylation_eigen_dir)


eigen_value<-read.table("/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/TensorQTL_related/Proteins_18modules_soft_thres_7_deep_split_3_me_0.25_min_size_20_4_removed_pcs.txt",
                        header=T)
eigen_value$sample<-rownames(eigen_value)
#we are now going to merge the expression with the eigen value, then with the pcs, then we are going to save
colnames(eigen_value)[1:length(colnames(eigen_value))-1]<-paste0(colnames(eigen_value)[1:length(colnames(eigen_value))-1],
                                                                 "_Eigen_",
                                                                 "Protein")
eigen_value<-eigen_value%>%dplyr::select(sample,everything())
eigen_value<-eigen_value%>%
  mutate(sample =  factor(sample   , levels = sorted_ids_key)) %>%
  arrange(sample   )
eigen_value_methylation_samples<-eigen_value%>%dplyr::filter(sample %in% (data_frame_covariates_methylation$nwdid))
dim(eigen_value_methylation_samples)

print(eigen_value_methylation_samples[1:10,1:10])
readr::write_tsv(eigen_value_methylation_samples,
                 paste0(methylation_eigen_dir,"/","Eigen_Proteins_Methylation_interaction.txt"))

#writing metabolite data
eigen_value<-read.table("/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/TensorQTL_related/52modules_soft_thres_8_deep_split_2_me_0.15_min_size_20.txt",
                        header=T)
eigen_value$sample<-rownames(eigen_value)
#we are now going to merge the expression with the eigen value, then with the pcs, then we are going to save
colnames(eigen_value)[1:length(colnames(eigen_value))-1]<-paste0(colnames(eigen_value)[1:length(colnames(eigen_value))-1],
                                                                 "_Eigen_",
                                                                 "Metabolite")
eigen_value<-eigen_value%>%dplyr::select(sample,everything())
eigen_value<-eigen_value%>%
  mutate(sample =  factor(sample   , levels = sorted_ids_key)) %>%
  arrange(sample   )
dim(eigen_value)
eigen_value_methylation_samples<-eigen_value%>%dplyr::filter(sample %in% (data_frame_covariates_methylation$nwdid))
dim(eigen_value_methylation_samples)

print(eigen_value_methylation_samples[1:10,1:10])

colnames(eigen_value)
colnames(eigen_value[,1:19])
colnames(eigen_value[,c(1,20:37)])
colnames(eigen_value[,c(1,38:53)])



readr::write_tsv(eigen_value_methylation_samples,
                 paste0(methylation_eigen_dir,"/","Eigen_Metabolite_Methylation_interaction.txt"))
readr::write_tsv(eigen_value[,1:19],
                 paste0(expression_eigen_dir,"/","Eigen_Metabolite_Methylation_interaction_chunk_1.txt"))
readr::write_tsv(eigen_value[,c(1,20:37)],
                 paste0(expression_eigen_dir,"/","Eigen_Metabolite_Methylation_interaction_chunk_2.txt"))
readr::write_tsv(eigen_value[,c(1,38:53)],
                 paste0(expression_eigen_dir,"/","Eigen_Metabolite_Methylation_interaction_chunk3.txt"))



