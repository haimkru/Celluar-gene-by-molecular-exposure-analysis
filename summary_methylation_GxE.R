###summary_methylation_GxE####
#/oak/stanford/groups/smontgom/hkrupkin/GxE/Code/summary_methylation_GxE
#/oak/stanford/groups/smontgom/hkrupkin/GxE/Code/predictions_of_hits_methylation.R
#Author:Haim Krupkin
#Date: 03/20/2026
#Date of correcitons: 04/26/2026
#further update of corrections on 05/11/2026 - specifically the purpose is to gett the total varaince explained acorss the 
#metabolome and proteome

#Description: 
#This script is meant to average to get the per cell type averages across samples

###package#####
library(stringr)
library(WGCNA)
###Analysis#####
#wrong expression files from a long time ago
#expression_per_cell_type<-list.files("/oak/stanford/groups/smontgom/dnachun/projects/topmed",
#           pattern="_exam1.tsv",
#           full.names = T)
#correct expression files, corrected on 04/27/2026
expression_per_cell_type<-list.files("/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/Ready_for_running_inputs/standardized_phenotypesV2/",
                                     pattern="_methylation.bed",
                                     full.names = T)

files <- sort(expression_per_cell_type)
files
# sanity checks
cat("N files:", length(files), "\n")
print(head(basename(files), 5))
print(tail(basename(files), 5))

### read + add true path column ####
first <- files[1]
cell_Type_first <- str_replace_all(basename(first), "_methylation.bed", "")

dt_first <- fread(first, sep = "\t", header = TRUE)
colnames(dt_first)[1]<-"chr"
dt_first<-dt_first%>%
  dplyr::select(-start,-end,-chr)
head(dt_first)
dt_avg_first <- dt_first[, .(
  phenotype_id = phenotype_id,
  tmp = rowMeans(.SD, na.rm = TRUE)
), .SDcols = patterns("^NWD")]
setnames(dt_avg_first, "tmp", cell_Type_first)

first_cell_type <- dt_avg_first

# Loop through remaining files
for (i in seq_along(files)[-1]) {
  f <- files[i]
  
  print("Current number of rows:")
  print(nrow(first_cell_type))
  print(f)
  
  cell_Type <- str_replace_all(basename(f), "_methylation.bed", "")
  
  dt <- fread(f, sep = "\t", header = TRUE)
  colnames(dt)[1]<-"chr"
  
  dt_avg <- dt[, .(
    phenotype_id = phenotype_id,
    tmp = rowMeans(.SD, na.rm = TRUE)
  ), .SDcols = patterns("^NWD")]
  
  setnames(dt_avg, "tmp", cell_Type)
  
  # ✅ FIX: merge into accumulating object
  first_cell_type <- merge(
    first_cell_type,
    dt_avg,
    by = "phenotype_id",
    all = TRUE   # optional: keeps all genes across files
  )
}

first_cell_type

#fread("/oak/stanford/groups/smontgom/dnachun/projects/topmed/nk_expr_exam1.tsv")

write.csv(first_cell_type,
          "/oak/stanford/groups/smontgom/hkrupkin/GxE/Code/average_methylation_in_diffrent_cell_types.csv")
head(first_cell_type)











