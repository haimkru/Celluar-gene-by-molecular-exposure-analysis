###summary_tpms_for_transcriptome_GxE####
#/oak/stanford/groups/smontgom/hkrupkin/GxE/Code/summary_tpms_for_transcriptome_GxE
#/oak/stanford/groups/smontgom/hkrupkin/GxE/Code/predictions_of_hits
#Author:Haim Krupkin
#Date: 03/20/2026
#Date of correcitons: 04/26/2026

#Description: 
#This script is meant to average to get the per cell type averages across samples

###package#####
library(stringr)
###Analysis#####
#wrong expression files from a long time ago
#expression_per_cell_type<-list.files("/oak/stanford/groups/smontgom/dnachun/projects/topmed",
#           pattern="_exam1.tsv",
#           full.names = T)
#correct expression files, corrected on 04/27/2026
expression_per_cell_type<-list.files("/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/expression",
           pattern="_exam1_standardized.tsv",
           full.names = T)

files <- sort(expression_per_cell_type)
files
# sanity checks
cat("N files:", length(files), "\n")
print(head(basename(files), 5))
print(tail(basename(files), 5))

### read + add true path column ####
first <- files[1]
cell_Type_first <- str_replace_all(basename(first), "_expr_exam1.tsv", "")

dt_first <- fread(first, sep = "\t", header = TRUE)

dt_avg_first <- dt_first[, .(Name, tmp = rowMeans(.SD, na.rm = TRUE)), .SDcols = !"Name"]
setnames(dt_avg_first, "tmp", cell_Type_first)

first_cell_type <- dt_avg_first

# Loop through remaining files
for (i in seq_along(files)[-1]) {
  f <- files[i]
  
  print("Current number of rows:")
  print(nrow(first_cell_type))
  print(f)
  
  cell_Type <- str_replace_all(basename(f), "_expression_exam1_standardized.tsv", "")
  
  dt <- fread(f, sep = "\t", header = TRUE)
  
  dt_avg <- dt[, .(Name, tmp = rowMeans(.SD, na.rm = TRUE)), .SDcols = !"Name"]
  setnames(dt_avg, "tmp", cell_Type)
  
  # ✅ FIX: merge into accumulating object
  first_cell_type <- merge(
    first_cell_type,
    dt_avg,
    by = "Name",
    all = TRUE   # optional: keeps all genes across files
  )
}

first_cell_type

#fread("/oak/stanford/groups/smontgom/dnachun/projects/topmed/nk_expr_exam1.tsv")

write.csv(first_cell_type,
          "/oak/stanford/groups/smontgom/hkrupkin/GxE/Code/average_expression_in_diffrent_cell_types.csv")
head(first_cell_type)











