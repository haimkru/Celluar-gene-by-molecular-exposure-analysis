###Single_example_good_example_plot################
#/labs/smontgom/grps_smontgom/hkrupkin/GxE/Code/single_variant_example_tensorQTL/Single_example_good_example_plot
#author: Haim Krupkin using AI
#date: 04/26/2026
#description: we are reverse engineerign tensorQTL to be able to replicate a single interaction from the GxE analysis we did
##to get singel exmaple you need to extract all the data relvent to that isngel example
#expression in cell type
#covarites
#genotypes
#and the valeus of all the exposures
#to get the genoypes you use plink, specifically you use:
#plink --bfile /oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/genotypes/mesa_1331samples.maf01.biallelic.intersect --extract /oak/stanford/groups/smontgom/dnachun/projects/topmed/t
#ensorqtl_inputs/genotypes/variant_of_intrest.txt --recode A  --out /oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/genotypes/variant_of_intrest.txt

# ============================================================
# Guided Lab: tensorQTL-style Interaction Pipeline in R
# ============================================================
# Purpose:
# - This file is meant to be run section-by-section in an editor.
# - It is deliberately repetitive and explicit.
# - After most sections, there is a "STOP AND INSPECT" note.
#
# Recommended workflow:
# 1) Open this file in your editor.
# 2) Run one section at a time.
# 3) After each section, inspect the printed objects before moving on.
#
# This guided lab uses the same logic as:
# - interaction_summary_and_plot.py
# - interaction_summary_and_plot_walkthrough.R
#
# But here the emphasis is learning the pipeline, not compactness.
###lets start by picking a good variant, with high maf and a good gene and not main effect.

###packags####

library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(data.table)
library(stringr)
####loading gene annoations
library(data.table)

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
INPUT_FILE <- "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/all_cell_types_sig_hits.tsv"


dt <- fread(INPUT_FILE, sep="\t", header=TRUE)
dt_backup<-dt
dt_interaction_only <- dt %>%
  group_by(cell_type, assay, Omics,term,phenotype_id) %>%
  filter(all(effect_type == "interaction_effect")) %>%
  ungroup()


dt_interaction_only_across_all_terms <- dt %>%
  group_by(cell_type, assay, Omics,phenotype_id) %>%
  filter(all(effect_type == "interaction_effect")) %>%
  ungroup()

dt_interaction_only_across_all_terms_and_cell_types <- dt %>%
  group_by( assay, Omics,phenotype_id) %>%
  filter(all(effect_type == "interaction_effect")) %>%
  ungroup()

dim(dt_interaction_only_across_all_terms)

dt_main_only <- dt %>%
  group_by(cell_type, assay, Omics,phenotype_id) %>%
  filter(all(effect_type == "main_effect")) %>%
  ungroup()
#dt<-dt_backup
dt_main_only_Across_cell_types <- dt %>%
  group_by( assay, Omics,phenotype_id) %>%
  filter(all(effect_type == "main_effect")) %>%
  ungroup()

dim(dt_main_only_Across_cell_types)
dt_both <- dt %>%
  group_by(cell_type, assay, Omics, phenotype_id) %>%
  filter(all(c("main_effect", "interaction_effect") %in% effect_type)) %>%
  ungroup()

genes_both<-unique(dt_both$phenotype_id)
genes_mainOnly<-unique(dt_main_only$phenotype_id)
genes_interactionOnly_across_all_terms<-unique(dt_interaction_only_across_all_terms$phenotype_id)
genes_mainOnly_across_cell_types<-unique(dt_main_only_Across_cell_types$phenotype_id)
genes_interactionOnly_across_cell_types<-unique(dt_interaction_only_across_all_terms_and_cell_types$phenotype_id)

intersect(genes_interactionOnly_across_all_terms,genes_mainOnly)[1:5]

intersect(genes_mainOnly_across_cell_types,genes_interactionOnly_across_cell_types)

count_per_gene<-dt %>%
  count(cell_type, assay, Omics, term, phenotype_id, effect_type) %>%
  arrange(desc(n))
interaction_counts<-count_per_gene%>%
  dplyr::filter(effect_type=="interaction_effect")

main_counts<-count_per_gene%>%
  dplyr::filter(effect_type=="main_effect")

unique_variant_ids<-dt%>%
  dplyr::select(variant_id)%>%
  distinct()

####
head(dt_interaction_only)
head(dt_main_only)

####

nrow(dt_both)+nrow(dt_main_only)+nrow(dt_interaction_only)

allele_frequencies<-fread("/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/genotypes/variant_frequencies.frq")

coloc_hits<-fread("/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/coloc_interaction_with_GWASesV2/candidates_pph4_under_0.8.csv")
dim(coloc_hits)
head(coloc_hits)

############################################################
###deciphering what is what in the colocs#####
# 1. Define a helper function to clean up naming messiness
clean_strings <- function(x) {
  x %>%
    # Standardize Cell Types: "B Cells" -> "B cell", "CD4T" -> "CD4 T" etc.
    str_replace_all("B Cells", "B cell") %>%
    str_replace_all("NK Cells", "NK") %>%
    # Standardize Modules: Remove "_Eigen" suffix if present
    str_remove_all("_Eigen") %>%
    # Trim whitespace just in case
    str_trim()
}

# 2. Re-generate the keys with cleaning applied
prepare_for_join_v2 <- function(df, gene_col, term_col) {
  df %>%
    mutate(
      gene_join  = str_remove(!!sym(gene_col), "\\..*"),
      # Apply the cleaning to both cell_type and term
      cell_clean = clean_strings(cell_type),
      term_clean = clean_strings(!!sym(term_col)),
      assay_join = str_to_title(assay),
      join_key   = paste(cell_clean, gene_join, term_clean, assay_join, sep = "_")
    )
}

# Process the Coloc hits
coloc_ready <- coloc_hits %>%
  mutate(assay = ifelse(grepl("expression", gxe_name), "Expression", "Methylation")) %>%
  prepare_for_join_v2("Gene", "module")

# Process your category keys
interaction_only_keys <- prepare_for_join_v2(dt_interaction_only_across_all_terms, "phenotype_id", "term")$join_key
main_only_keys        <- prepare_for_join_v2(dt_main_only_Across_cell_types, "phenotype_id", "term")$join_key

# 3. Re-assign Categories
coloc_classified <- coloc_ready %>%
  mutate(effect_category = case_when(
    join_key %in% interaction_only_keys ~ "Interaction Only",
    join_key %in% main_only_keys        ~ "Main Only",
    # If it's a hit but NOT only Interaction and NOT only Main, it's Both
    TRUE                                ~ "Both"
  ))

# 4. Check alignment again
cat("Sample Coloc Key: ", coloc_ready$join_key[1], "\n")
cat("Sample Interaction Key: ", interaction_only_keys[grepl("B cell", interaction_only_keys)][1], "\n")

# View updated counts
print(table(coloc_classified$effect_category))
#########
library(data.table)
library(stringr)

small_coloc_classified<-coloc_classified%>%
  mutate(term=module)%>%
  dplyr::select(gwas_short,Gene,term,cell_type,PP.H4,effect_category)
head(small_coloc_classified)
# Ensure everything is a data.table
setDT(coloc_hits)
setDT(small_coloc_classified)
setDT(dt)
setDT(allele_frequencies)

# --- 1. Prepare Coloc Reference ---

# A. Extract just the metadata we need and hardcode Metabolome
coloc_hits_meta <- coloc_hits[, .(Gene, term = module, cell_type, gxe_name, source_file, Omics = "Metabolome")]

# B. Deduplicate it BEFORE the merge so we have exactly 1 row per join key
# This strictly prevents the vecseq Cartesian explosion error
coloc_hits_unique <- unique(coloc_hits_meta, by = c("Gene", "term", "cell_type"))

# C. Now the merge is safe (one-to-one or many-to-one)
coloc_meta <- merge(
  small_coloc_classified, 
  coloc_hits_unique, 
  by = c("Gene", "term", "cell_type"),
  all.x = TRUE
)

# Fast string standardizer function
std_str <- function(x) {
  str_to_lower(str_trim(str_remove_all(str_replace_all(str_replace_all(x, "B Cells", "B cell"), "NK Cells", "NK"), "_Eigen")))
}

# Apply to coloc and get unique combinations
coloc_for_sampling <- unique(coloc_meta[, .(
  gene_join   = str_remove(Gene, "\\..*"),
  cell_join   = std_str(cell_type),
  term_join   = std_str(term),
  assay_join  = str_to_lower(ifelse(grepl("expression", gxe_name), "Expression", "Methylation")),
  omics_join  = str_to_lower(Omics), # Uses our hardcoded "Metabolome" string
  effect_category
)])

valid_genes <- unique(coloc_for_sampling$gene_join)


# --- 2. Aggressive Prefilter of Allele Frequencies ---
high_maf <- allele_frequencies[MAF > 0.3, .(variant_id = SNP, MAF)]


# --- 3. Aggressive Prefilter of DT (The Speed Optimization) ---

# A. Prefilter by Gene/Phenotype without modifying DT directly:
all_phenos <- unique(dt$phenotype_id)
keep_phenos <- all_phenos[str_remove(all_phenos, "\\..*") %in% valid_genes]

dt_sub <- dt[phenotype_id %in% keep_phenos]

# B. Prefilter by MAF:
dt_sub <- dt_sub[variant_id %in% high_maf$variant_id]


# --- 4. Format and Join (Running on a tiny fraction of the data) ---

# Format remaining rows using update-in-place `:=`
dt_sub[, `:=`(
  gene_join  = str_remove(phenotype_id, "\\..*"),
  cell_join  = std_str(cell_type),
  term_join  = std_str(term),
  assay_join = str_to_lower(assay),
  omics_join = str_to_lower(Omics)
)]

# Inner join with Coloc to assign categories (instantly drops mismatches)
dt_matched <- merge(
  dt_sub, 
  coloc_for_sampling, 
  by = c("gene_join", "cell_join", "term_join", "assay_join", "omics_join"), 
  all = FALSE  
)

# Bring MAF back in
dt_final <- merge(dt_matched, high_maf, by = "variant_id", all = FALSE)
dim(dt_final)

# --- 5. Stratified Sampling ---
set.seed(42)

# Native data.table sampling by group 
final_variants_100 <- dt_final[, .SD[sample(.N, min(.N, 50))], by = effect_category]

# --- 6. Verification ---
print("Final Variant Counts per Category:")
#still missing main effect only so we are going to load that into the plan ya habibi
dt_main <- dt_main_only_Across_cell_types[
  dt_main_only_Across_cell_types$variant_id %in% high_maf$variant_id,
]
dim(dt_main)
#btw generally working in DT is faster
setDT(dt_main)
# Strip gene IDs and drop anything that colocalized
dt_main[, temp_gene := str_remove(phenotype_id, "\\..*")]
dt_main_pure <- dt_main[!(temp_gene %in% valid_genes)]
dim(dt_main_pure)
# Sample 50
set.seed(42)
sampled_main <- dt_main_pure[sample(.N, min(.N, 50))]

# Format to match the other table seamlessly
sampled_main[, `:=`(
  gene_join  = temp_gene,
  cell_join  = std_str(cell_type),
  term_join  = std_str(term),
  assay_join = str_to_lower(assay),
  omics_join = str_to_lower(Omics),
  effect_category = "Main Only"
)]
sampled_main[, temp_gene := NULL]

# Bring MAF back in
sampled_main <- merge(sampled_main, high_maf, by = "variant_id", all.x = TRUE)


# --- 5. Combine and Verify ---
# Intersect ensures we only bind matching columns
common_cols <- intersect(names(final_variants_100), names(sampled_main))
common_cols
final_variants_150 <- rbind(final_variants_100[, ..common_cols], sampled_main[, ..common_cols])

print("Final Variant Counts per Category:")
print(table(final_variants_150$effect_category))

write.csv(final_variants_150,"/labs/smontgom/grps_smontgom/hkrupkin/GxE/Code/single_variant_example_tensorQTL/example_150_variants.csv")

###creating the genotype file or these varaints####
bfile_path <- "/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/genotypes/mesa_1331samples.maf01.biallelic.intersect"

output_dir <- "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Code/variant_followup"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

variant_list_file <- file.path(output_dir, "target_150_variants.txt")
out_prefix        <- file.path(output_dir, "genotypes_150_variants")

# 2. Export the Variant List
variant_list_file
fwrite(
  final_variants_150[, .(variant_id)], 
  file = variant_list_file, 
  col.names = FALSE, 
  quote = FALSE
)

# 3. Construct the Bash Command (Now with module load)
# Using '&&' ensures PLINK only runs if the module loads successfully
plink_cmd <- paste(
  "module load plink && plink",
  "--bfile", bfile_path,
  "--extract", variant_list_file,
  "--recode A",
  "--out", out_prefix
)

cat("Running bash command:\n", plink_cmd, "\n\n")

# 4. Execute from inside R
# in reality you need to run it inside the bash system... wow
#system(plink_cmd)

# 5. Load the Genotypes Back into R
geno_raw_file <- paste0(out_prefix, ".raw")

if (file.exists(geno_raw_file)) {
  genotypes_dt <- fread(geno_raw_file)
  cat("Successfully loaded", nrow(genotypes_dt), "samples and", ncol(genotypes_dt) - 6, "variants.\n")
} else {
  stop("PLINK output not found. Check the console for PLINK errors.")
}


#########################################################################

################################################################################
# Script: Extract Specific Phenotype Rows from Large BED Files
# Purpose: Prevents file collisions by naming outputs: Gene_Cell_Assay.txt
################################################################################
library(data.table)

# --- 1. SET PATHS ---
pheno_dir <- "/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/Ready_for_running_inputs/standardized_phenotypesV2"
out_dir   <- "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Code/single_variant_example_tensorQTL"

# --- 2. PREPARE MAPPING DATA WITH TRANSLATION ---
targets <- unique(final_variants_150[, .(phenotype_id, cell_join, assay_join)])

# TRANSLATION LAYER: Fix the naming quirks
targets[, cell_clean := tolower(cell_join)]           # Lowercase everything
targets[, cell_clean := gsub(" ", "", cell_clean)]    # Remove spaces ("b cell" -> "bcell")
targets[, cell_clean := gsub("[⁺+]", "", cell_clean)] # Remove plus signs ("cd8⁺t" -> "cd8t")
targets[, cell_clean := gsub("neutro", "neu", cell_clean)] # Shorten "neutro" to "neu"

# Now construct the paths
targets[, bed_path := file.path(pheno_dir, paste0(cell_clean, "_", assay_join, ".bed"))]

# --- 3. LOOP AND EXTRACT ---
unique_beds <- unique(targets$bed_path)

for (bed in unique_beds) {
  
  if (!file.exists(bed)) {
    # Helpful debugging: print what it was looking for vs what exists
    message("!! Still missing: ", basename(bed))
    next
  }
  
  current_subset <- targets[bed_path == bed]
  cat("-> Extracting", nrow(current_subset), "genes from:", basename(bed), "\n")
  
  # Get header
  header_line <- system(paste("head -n 1", shQuote(bed)), intern = TRUE)
  
  for (i in 1:nrow(current_subset)) {
    gene  <- current_subset$phenotype_id[i]
    cell  <- current_subset$cell_clean[i]
    assay <- current_subset$assay_join[i]
    
    file_name <- paste0(gene, "_", cell, "_", assay, ".txt")
    out_path  <- file.path(out_dir, file_name)
    
    writeLines(header_line, out_path)
    
    awk_cmd <- paste0(
      "awk -v gene=", shQuote(gene), 
      " '$4 == gene' ", shQuote(bed), 
      " >> ", shQuote(out_path)
    )
    print("ruuning command: ")
    print(awk_cmd)
    system(awk_cmd)
  }
}
##########
###getting data for a specific correct row####
# --- 1. SETUP PATHS ---
base_dir   <- "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Code/single_variant_example_tensorQTL"
covar_dir  <- "/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/Ready_for_running_inputs/covariatesV2"
setwd(base_dir)

# --- 2. SELECT YOUR TARGET ROW ---
# Let's say we want to analyze the first row
nrow(final_variants_150)
for (index_row in c(1:nrow(final_variants_150))){
#for (index_row in c(100:nrow(final_variants_150))){
    
  print("working on row number : ")
  print(index_row)
  
row_idx <- index_row

#View(target)

target  <- final_variants_150[row_idx]
target
# Extract Metadata
v_id    <- target$variant_id     # e.g., "chr3_135122610_C_T_b38"
p_id    <- target$phenotype_id   # e.g., "ENSG00000241905.1"
assay   <- target$assay_join     # e.g., "expression"
term_id <- target$term           # e.g., "MEtan_Eigen"

# Clean Cell Name (The translation layer we built)
cell_clean <- tolower(target$cell_join)
cell_clean <- gsub(" ", "", cell_clean)
cell_clean <- gsub("[⁺+]", "", cell_clean)
cell_clean <- gsub("neutro", "neu", cell_clean)

# --- 3. CONSTRUCT FILE PATHS ---
# Expression file created in the previous sub-task
expression_path <- file.path(base_dir, paste0(p_id, "_", cell_clean, "_", assay, ".txt"))

# Covariate file from the OAK directory
covariates_path <- file.path(covar_dir, paste0(cell_clean, "_", assay, "_covariates.txt"))

# Genotypes (The master .raw file from PLINK)
genotypes_path  <- file.path(base_dir, "genotypes_150_variants.raw")

# --- 4. DATA EXTRACTION ---

# A. Load Expression (Just the row we extracted)
expr_dt <- fread(expression_path)
# Transpose to get samples as rows for the model
# Row 1 is metadata, so we take everything from column 5 onwards
expr_s  <- as.numeric(expr_dt[1, 5:ncol(expr_dt), with=FALSE])
names(expr_s) <- colnames(expr_dt)[5:ncol(expr_dt)]

# B. Load Genotype using your extract_variant_series function
geno_raw <- fread(geno_raw_file)
var_res  <- extract_variant_series(variant_df = geno_raw, 
                                   variant_id = v_id)
geno_s   <- var_res$series

# C. Load Covariates
covar_dt <- fread(covariates_path)
# Ensure covariates are transposed: samples as rows, covariates as columns
# Usually in TensorQTL, covars are samples in columns. Let's fix that:
t_covar  <- as.data.frame(t(covar_dt[, -1, with=FALSE]))
colnames(t_covar) <- covar_dt[[1]]
t_covar$IID <- rownames(t_covar)


###########

# STOP AND INSPECT:
# - Are these the file names you actually want?
# - Is plot_environment a shorthand like 'tan' or an exact column name?


# ============================================================
# SECTION 2: Lightweight helpers
# ============================================================

# STOP AND INSPECT:
# - These are generic helpers: path resolution, reading, rownames, numeric conversion.
# - Nothing tensorQTL-specific has happened yet.


# ============================================================
# SECTION 3: Parse expression into one sample-indexed vector
# ============================================================
expr_s <- expr_s #still loading expression data
phenotype_id_resolved <- p_id #getting gene id

cat("Expression vector length:", length(expr_s), "\n")
print(head(expr_s))

geno_s <- var_res$series #continue from above
variant_id_resolved <- var_res$variant_id #continue from above 

cat("Genotype vector length:", length(geno_s), "\n")
cat("Resolved variant ID:", variant_id_resolved, "\n")
print(table(geno_s, useNA = "ifany"))


extract_interaction_df <- function(interaction_df, sample_hint) {
  df <- maybe_set_rownames(interaction_df)
  idx_overlap <- sum(rownames(df) %in% sample_hint)
  col_overlap <- sum(colnames(df) %in% sample_hint)
  
  if (col_overlap > idx_overlap && col_overlap >= 3) {
    df <- as.data.frame(t(df), stringsAsFactors = FALSE, check.names = FALSE)
    log_msg("  - Interaction table transposed to samples x environments")
  }
  
  num_df <- data.frame(lapply(df, function(x) suppressWarnings(as.numeric(x))),
                       row.names = rownames(df), check.names = FALSE)
  keep <- colSums(!is.na(num_df)) > 0
  num_df <- num_df[, keep, drop = FALSE]
  if (ncol(num_df) == 0) stop("No numeric interaction/environment columns found")
  
  log_msg(sprintf("  - Interaction parsed with %d environment(s)", ncol(num_df)))
  num_df
}

extract_covariates_df <- function(cov_df, sample_hint) {
  df <- maybe_set_rownames(cov_df)
  idx_overlap <- sum(rownames(df) %in% sample_hint)
  col_overlap <- sum(colnames(df) %in% sample_hint)
  
  if (col_overlap > idx_overlap && col_overlap >= 3) {
    df <- as.data.frame(t(df), stringsAsFactors = FALSE, check.names = FALSE)
    log_msg("  - Covariate table transposed to samples x covariates")
  }
  
  num_df <- data.frame(lapply(df, function(x) suppressWarnings(as.numeric(x))),
                       row.names = rownames(df), check.names = FALSE)
  keep <- colSums(!is.na(num_df)) > 0
  num_df <- num_df[, keep, drop = FALSE]
  if (ncol(num_df) == 0) return(data.frame(row.names = rownames(df)))
  
  log_msg(sprintf("  - Covariates parsed with %d covariate(s)", ncol(num_df)))
  num_df
}

cov_raw <- t_covar #we load the covariate data

resolve_input_path <- function(path) {
  if (file.exists(path)) return(path)
  candidates <- c(paste0(path, ".txt"), paste0(path, ".txt.raw"), paste0(path, ".raw"))
  for (candidate in candidates) {
    if (file.exists(candidate)) return(candidate)
  }
  stop(sprintf("Input file not found: %s", path))
}
extract_interaction_df <- function(interaction_df, sample_hint) {
  df <- maybe_set_rownames(interaction_df)
  idx_overlap <- sum(rownames(df) %in% sample_hint)
  col_overlap <- sum(colnames(df) %in% sample_hint)
  
  if (col_overlap > idx_overlap && col_overlap >= 3) {
    df <- as.data.frame(t(df), stringsAsFactors = FALSE, check.names = FALSE)
    log_msg("  - Interaction table transposed to samples x environments")
  }
  
  num_df <- data.frame(lapply(df, function(x) suppressWarnings(as.numeric(x))),
                       row.names = rownames(df), check.names = FALSE)
  keep <- colSums(!is.na(num_df)) > 0
  num_df <- num_df[, keep, drop = FALSE]
  if (ncol(num_df) == 0) stop("No numeric interaction/environment columns found")
  
  log_msg(sprintf("  - Interaction parsed with %d environment(s)", ncol(num_df)))
  num_df
}

extract_covariates_df <- function(cov_df, sample_hint) {
  df <- maybe_set_rownames(cov_df)
  idx_overlap <- sum(rownames(df) %in% sample_hint)
  col_overlap <- sum(colnames(df) %in% sample_hint)
  
  if (col_overlap > idx_overlap && col_overlap >= 3) {
    df <- as.data.frame(t(df), stringsAsFactors = FALSE, check.names = FALSE)
    log_msg("  - Covariate table transposed to samples x covariates")
  }
  
  num_df <- data.frame(lapply(df, function(x) suppressWarnings(as.numeric(x))),
                       row.names = rownames(df), check.names = FALSE)
  keep <- colSums(!is.na(num_df)) > 0
  num_df <- num_df[, keep, drop = FALSE]
  if (ncol(num_df) == 0) return(data.frame(row.names = rownames(df)))
  
  log_msg(sprintf("  - Covariates parsed with %d covariate(s)", ncol(num_df)))
  num_df
}
read_table <- function(path) {
  if (!file.exists(path)) stop(sprintf("Input file not found: %s", path))
  tryCatch(
    read.table(path, header = TRUE, sep = "", quote = "", comment.char = "",
               check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) {
      read.table(path, header = TRUE, sep = "\t", quote = "", comment.char = "",
                 check.names = FALSE, stringsAsFactors = FALSE)
    }
  )
}

find_interaction_file <- function(dt, interaction_dir) {
  
  omics_map <- c(
    "metabolome" = "Metabolite",
    "proteins"   = "Proteins",
    "proteome"   = "Proteins"
  )
  
  assay_map <- c(
    "expression"  = "Expression",
    "methylation" = "Methylation"
  )
  
  # omics that have no chunking
  no_chunk_omics <- c("proteins", "proteome")
  
  dt <- copy(dt)
  
  dt[, chunk := regmatches(path, regexpr("chunk_?\\d+", path))]
  dt[chunk == "", chunk := NA]
  
  dt[, interaction_file := ifelse(
    omics_join %in% no_chunk_omics,
    paste0("Eigen_", omics_map[omics_join], "_", assay_map[assay_join], "_interaction.txt"),
    paste0("Eigen_", omics_map[omics_join], "_", assay_map[assay_join], "_interaction_", chunk, ".txt")
  )]
  
  dt[, interaction_path := file.path(interaction_dir, interaction_file)]
  dt[, file_exists := file.exists(interaction_path)]
  
  if (any(!dt$file_exists)) {
    warning("Some interaction files not found:\n",
            paste(dt[file_exists == FALSE, interaction_file], collapse = "\n"))
  }
  
  return(dt[, .(interaction_file, interaction_path, file_exists)])
}
interaction_file<-find_interaction_file(
  dt = target,
  interaction_dir = "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Code/single_variant_example_tensorQTL"
)
interaction_path<-interaction_file$interaction_path
print("interaction path is :")
print(interaction_path)

#interaction_path <- "Eigen_Metabolite_Expression_interaction_chunk_1.txt"
int_raw <- read_table(resolve_input_path(interaction_path)) #we load the enviromenta data
interaction_df <- extract_interaction_df(int_raw, names(expr_s)) #we make enviromental data into a pretty data frame
cov_df <- extract_covariates_df(cov_raw, names(expr_s)) #we make the covariate data into a pretty data frame as well

cat("Interaction shape:", nrow(interaction_df), "x", ncol(interaction_df), "\n")
cat("Covariate shape:", nrow(cov_df), "x", ncol(cov_df), "\n")
print(colnames(interaction_df)[1:min(5, ncol(interaction_df))])

# STOP AND INSPECT:
# - interaction_df should be samples x environments.
interaction_df[1:5,1:5]
# - cov_df should be samples x covariates.
#we did lose the covaraite names which is sad but not super important.
ncol(cov_df)
cov_df[1:5,1:5]
# - rownames of both should be sample IDs.

# ============================================================
# SECTION 6: Align shared samples
# ============================================================

align_on_samples <- function(expr_s, geno_s, interaction_df, cov_df = NULL) {
  common <- Reduce(intersect, list(
    names(expr_s),
    names(geno_s),
    rownames(interaction_df),
    if (!is.null(cov_df) && ncol(cov_df) > 0) rownames(cov_df) else names(expr_s)
  ))
  
  if (length(common) < 10) stop(sprintf("Too few overlapping samples (%d)", length(common)))
  common <- sort(common)
  
  list(
    expr_s = expr_s[common],
    geno_s = geno_s[common],
    interaction_df = interaction_df[common, , drop = FALSE],
    cov_df = if (!is.null(cov_df) && ncol(cov_df) > 0) cov_df[common, , drop = FALSE] else NULL
  )
}

aligned <- align_on_samples(expr_s, geno_s, interaction_df, cov_df) # we make sure they are all the same order and on the same data
expr_s <- aligned$expr_s
geno_s <- aligned$geno_s
interaction_df <- aligned$interaction_df
cov_df <- aligned$cov_df
expr_s[1:5]
geno_s[1:5]
interaction_df[1:5,1:5]
cov_df[1:5,1:5]

cat("Aligned sample count:", length(expr_s), "\n")
cat("First few aligned sample IDs:\n")
print(head(names(expr_s)))

# STOP AND INSPECT:
# - From here on, every model matrix is using exactly the same samples in the same order.
# - This alignment step is essential; without it, coefficients are meaningless.


# ============================================================
# SECTION 7: Understand residualization
# ============================================================
#i made all impute mean commands comments, as we dont need, and cannot, impute genotype here.

impute_mean <- function(genotypes_mat, missing = -9) {
  out <- genotypes_mat
  for (r in seq_len(nrow(out))) {
    miss <- out[r, ] == missing
    if (any(miss, na.rm = TRUE)) {
      out[r, miss] <- mean(out[r, !miss], na.rm = TRUE)
    }
  }
  out
}
#i think this is just QR reguilirization
#we take the rows, center them across samples, project onto covariate space, and substract that... so we get what the covariates cannot explain
#
create_residualizer <- function(C) {
  C_center <- scale(C, center = TRUE, scale = FALSE)
  Q <- qr.Q(qr(C_center))
  dof <- nrow(C) - 2 - ncol(C)
  list(Q = Q, dof = dof)
}

residualize_rows <- function(M, residualizer, center = TRUE) {
  M <- as.matrix(M)
  M_center <- M - rowMeans(M)
  if (center) {
    M_center - (M_center %*% residualizer$Q) %*% t(residualizer$Q)
  } else {
    M - (M_center %*% residualizer$Q) %*% t(residualizer$Q)
  }
}

# Build raw model pieces for one variant / one phenotype.
phenotype_mat <- matrix(as.numeric(expr_s), nrow = 1)
genotype_mat <- matrix(as.numeric(geno_s), nrow = 1)
interaction_mat <- as.matrix(interaction_df)
#genotype_mat <- impute_mean(genotype_mat) #we are not imputing as we only have one matrix...

#this centers our data
g0 <- genotype_mat - rowMeans(genotype_mat) #this only centers
i0 <- scale(interaction_mat, center = TRUE, scale = FALSE) #this actually just centers, no scaling.
p0 <- phenotype_mat - rowMeans(phenotype_mat) #this also just centers

cat("Raw centered shapes:\n")
cat("g0:", dim(g0), "\n")
cat("i0:", dim(i0), "\n")
cat("p0:", dim(p0), "\n")

if (!is.null(cov_df) && ncol(cov_df) > 0) {
  R <- create_residualizer(as.matrix(cov_df))
  p0_resid <- residualize_rows(p0, R, center = FALSE) #here we residualize the expression
  g0_resid <- residualize_rows(g0, R, center = FALSE) #here we residulize the genotpy
  i0_resid <- t(residualize_rows(t(i0), R, center = FALSE)) #here we residualize the enviroment (our interaction term)
  cat("Residualizer degrees of freedom base:", R$dof, "\n")
} else {
  R <- NULL
  p0_resid <- p0
  g0_resid <- g0
  i0_resid <- i0
}

# STOP AND INSPECT:
# - Residualization means: remove the part explained by covariates.
# - Q spans the covariate space, and Q Q^T is the projection operator.
# - After this step, phenotype/genotype/environment are orthogonal to covariates.


# ============================================================
# SECTION 8: Build the full interaction design matrix by hand
# ============================================================
plot_environment<-term_id
requested_env <- tolower(plot_environment)
env_matches <- colnames(interaction_df)[grepl(requested_env, tolower(colnames(interaction_df)), fixed = TRUE)]
env_matches
if (length(env_matches) == 0) stop(sprintf("No environment matched '%s'", plot_environment))
env_name_manual <- env_matches[1]

ni <- ncol(interaction_mat)
k_env <- match(env_name_manual, colnames(interaction_df))

# Build all g:i terms
ngi <- array(NA_real_, dim = c(1, ncol(genotype_mat), ni))
for (k in seq_len(ni)) {
  tmp <- genotype_mat * matrix(interaction_mat[, k], nrow = 1)
  ngi[1, , k] <- tmp - rowMeans(tmp)
}

if (!is.null(R)) {
  for (k in seq_len(ni)) {
    ngi[1, , k] <- residualize_rows(matrix(ngi[1, , k], nrow = 1, ncol = ncol(genotype_mat)), R, center = FALSE)
  }
}

X <- cbind(
  g0_resid[1, ],
  i0_resid,
  do.call(cbind, lapply(seq_len(ni), function(k) matrix(ngi[1, , k], nrow = ncol(genotype_mat), ncol = 1)))
)

colnames(X) <- c(
  "g",
  paste0("i:", colnames(interaction_df)),
  paste0("g:i:", colnames(interaction_df))
)

y <- as.numeric(p0_resid[1, ])

cat("Design matrix shape:", dim(X), "\n")
cat("First 10 coefficient names:\n")
print(colnames(X)[1:min(10, ncol(X))])
dim(X)
X
# STOP AND INSPECT:
# - X has one row per sample.
# - X has columns for g, every i, and every g:i.
# - If there are ni environments, then number of coefficients is 1 + 2*ni.


# ============================================================
# SECTION 9: Solve OLS by matrix algebra
# ============================================================

XtX_inv <- solve(crossprod(X))
beta <- XtX_inv %*% crossprod(X, y)
resid <- y - as.numeric(X %*% beta)
resid
if (!is.null(R)) {
  dof <- R$dof - 2 * ni
} else {
  dof <- length(y) - 2 - 2 * ni
}

###so they do turn the p.value from the t.stat, as seen here.

rss <- sum(resid^2)
se <- sqrt(diag(XtX_inv) * rss / dof)
tstat <- as.numeric(beta) / se
pval <- 2 * pt(-abs(tstat), df = dof)
pval
#this is a msall dataframe to see the genomic pvalues and such####

manual_df <- data.frame(
  term = colnames(X),
  beta = as.numeric(beta),
  se = se,
  tstat = tstat,
  pval = pval,
  row.names = NULL,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

print((manual_df))
# STOP AND INSPECT:
# - This is the heart of the pipeline.
# - Everything else is parsing, alignment, exporting, and plotting.
# - The interaction p-value for a chosen environment comes from one row here.


# ============================================================
# SECTION 10: Wrap this into a reusable summary function
# ============================================================

get_allele_stats <- function(genotype_mat) {
  n2 <- 2 * ncol(genotype_mat)
  af <- rowSums(genotype_mat) / n2
  ma_samples <- integer(length(af))
  ma_count <- integer(length(af))
  for (i in seq_along(af)) {
    g <- genotype_mat[i, ]
    if (af[i] <= 0.5) {
      ma_samples[i] <- sum(g > 0.5)
      ma_count[i] <- as.integer(sum(g[g > 0.5]))
    } else {
      ma_samples[i] <- sum(g < 1.5)
      ma_count[i] <- as.integer(n2 - sum(g[g > 0.5]))
    }
  }
  list(af = af, ma_samples = ma_samples, ma_count = ma_count)
}

compute_interaction_stats <- function(expr_s, geno_s, interaction_df, cov_df, variant_id, phenotype_id) {
  phenotype_mat <- matrix(as.numeric(expr_s), nrow = 1)
  genotype_mat <- matrix(as.numeric(geno_s), nrow = 1)
  interaction_mat <- as.matrix(interaction_df)
  #genotype_mat <- impute_mean(genotype_mat)
  
  g0 <- genotype_mat - rowMeans(genotype_mat)
  i0 <- scale(interaction_mat, center = TRUE, scale = FALSE)
  p0 <- phenotype_mat - rowMeans(phenotype_mat)
  
  gi0 <- array(NA_real_, dim = c(1, ncol(genotype_mat), ncol(interaction_mat)))
  for (k in seq_len(ncol(interaction_mat))) {
    tmp <- genotype_mat * matrix(interaction_mat[, k], nrow = 1)
    gi0[1, , k] <- tmp - rowMeans(tmp)
  }
  
  if (!is.null(cov_df) && ncol(cov_df) > 0) {
    R <- create_residualizer(as.matrix(cov_df))
    p0 <- residualize_rows(p0, R, center = FALSE)
    g0 <- residualize_rows(g0, R, center = FALSE)
    i0 <- t(residualize_rows(t(i0), R, center = FALSE))
    for (k in seq_len(ncol(interaction_mat))) {
      gi0[1, , k] <- residualize_rows(matrix(gi0[1, , k], nrow = 1, ncol = ncol(genotype_mat)), R, center = FALSE)
    }
    dof <- R$dof - 2 * ncol(interaction_mat)
  } else {
    dof <- ncol(phenotype_mat) - 2 - 2 * ncol(interaction_mat)
  }
  
  X <- cbind(
    g0[1, ],
    i0,
    do.call(cbind, lapply(seq_len(ncol(interaction_mat)), function(k) matrix(gi0[1, , k], nrow = ncol(genotype_mat), ncol = 1)))
  )
  y <- as.numeric(p0[1, ])
  XtX_inv <- solve(crossprod(X))
  beta <- XtX_inv %*% crossprod(X, y)
  resid <- y - as.numeric(X %*% beta)
  rss <- sum(resid^2)
  se <- sqrt(diag(XtX_inv) * rss / dof)
  tstat <- as.numeric(beta) / se
  
  allele <- get_allele_stats(genotype_mat)
  envs <- colnames(interaction_df)
  ni <- length(envs)
  
  row <- list(
    phenotype_id = phenotype_id,
    variant_id = variant_id,
    af = allele$af[1],
    ma_samples = allele$ma_samples[1],
    ma_count = allele$ma_count[1],
    pval_g = 2 * pt(-abs(tstat[1]), df = dof),
    b_g = as.numeric(beta[1]),
    b_g_se = se[1]
  )
  
  for (i in seq_len(ni)) {
    row[[paste0("pval_", envs[i])]] <- 2 * pt(-abs(tstat[1 + i]), df = dof)
    row[[paste0("b_", envs[i])]] <- as.numeric(beta[1 + i])
    row[[paste0("b_", envs[i], "_se")]] <- se[1 + i]
  }
  
  for (i in seq_len(ni)) {
    idx <- 1 + ni + i
    row[[paste0("pval_g-", envs[i])]] <- 2 * pt(-abs(tstat[idx]), df = dof)
    row[[paste0("b_g-", envs[i])]] <- as.numeric(beta[idx])
    row[[paste0("b_g-", envs[i], "_se")]] <- se[idx]
  }
  
  as.data.frame(row, check.names = FALSE, stringsAsFactors = FALSE)
}

stats_df <- compute_interaction_stats(expr_s, geno_s, interaction_df, cov_df, variant_id_resolved, phenotype_id_resolved)
#
stats_df
#write.table(stats_df, stats_out, sep = "\t", quote = FALSE, row.names = FALSE)
print(stats_df[, 1:min(12, ncol(stats_df)), drop = FALSE])

# STOP AND INSPECT:
# - stats_df is the final one-row summary table.
# - It contains one set of coefficients/p-values for every environment.
# - Compare its selected environment values to manual_df above.


# ============================================================
# SECTION 11: Build residualized export used for plotting
# ============================================================

build_residualized_regression_dataframe <- function(expr_s, geno_s, interaction_df, cov_df = NULL) {
  phenotype_mat <- matrix(as.numeric(expr_s), nrow = 1)
  genotype_mat <- matrix(as.numeric(geno_s), nrow = 1)
  interaction_mat <- as.matrix(interaction_df)
  #genotype_mat <- impute_mean(genotype_mat)
  
  g0 <- genotype_mat - rowMeans(genotype_mat)
  i0 <- scale(interaction_mat, center = TRUE, scale = FALSE)
  p0 <- phenotype_mat - rowMeans(phenotype_mat)
  
  gi0 <- array(NA_real_, dim = c(1, ncol(genotype_mat), ncol(interaction_mat)))
  for (k in seq_len(ncol(interaction_mat))) {
    tmp <- genotype_mat * matrix(interaction_mat[, k], nrow = 1)
    gi0[1, , k] <- tmp - rowMeans(tmp)
  }
  
  if (!is.null(cov_df) && ncol(cov_df) > 0) {
    R <- create_residualizer(as.matrix(cov_df))
    p0 <- residualize_rows(p0, R, center = FALSE)
    g0 <- residualize_rows(g0, R, center = FALSE)
    i0 <- t(residualize_rows(t(i0), R, center = FALSE))
    for (k in seq_len(ncol(interaction_mat))) {
      gi0[1, , k] <- residualize_rows(matrix(gi0[1, , k], nrow = 1, ncol = ncol(genotype_mat)), R, center = FALSE)
    }
  }
  
  out <- data.frame(
    sample_id = names(expr_s),
    phenotype_centered_residualized = as.numeric(p0[1, ]),
    genotype_centered_residualized = as.numeric(g0[1, ]),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  for (k in seq_len(ncol(interaction_df))) {
    env <- colnames(interaction_df)[k]
    out[[paste0(env, "_centered_residualized")]] <- as.numeric(i0[, k])
    out[[paste0("g_x_", env, "_centered_residualized")]] <- as.numeric(gi0[1, , k])
  }
  
  out
}

residualized_df <- build_residualized_regression_dataframe(expr_s, geno_s, interaction_df, cov_df)
#write.table(residualized_df, residualized_out, sep = "\t", quote = FALSE, row.names = FALSE)
print(head(residualized_df[, 1:min(8, ncol(residualized_df)), drop = FALSE]))

# STOP AND INSPECT:
# - This file contains the exact transformed variables used by the model.
# - Plotting these is more faithful than plotting raw expression/raw environment.


# ============================================================
# SECTION 12: Compare to reference
# ============================================================

compare_with_reference <- function(stats_df, reference_path, phenotype_id, variant_id, max_ratio_diff_ok = 0.005) {
  ref <- read.csv(reference_path)
  keep <- as.character(ref$phenotype_id) == phenotype_id & as.character(ref$variant_id) == variant_id
  if (!any(keep)) stop("Reference row not found")
  
  ref_row <- ref[which(keep)[1], , drop = FALSE]
  shared <- intersect(colnames(stats_df), colnames(ref_row))
  shared <- setdiff(shared, c("phenotype_id", "variant_id"))
  
  pieces <- list()
  for (nm in shared) {
    a <- suppressWarnings(as.numeric(stats_df[[nm]][1]))
    b <- suppressWarnings(as.numeric(ref_row[[nm]][1]))
    if (is.na(a) || is.na(b)) next
    
    abs_diff <- abs(a - b)
    rel_diff_ratio <- if (b == 0) if (a == 0) 0 else Inf else abs_diff / abs(b)
    
    pieces[[length(pieces) + 1]] <- data.frame(
      column = nm,
      computed = a,
      reference = b,
      abs_diff = abs_diff,
      rel_diff_ratio = rel_diff_ratio,
      within_ratio_diff_limit = rel_diff_ratio <= max_ratio_diff_ok,
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, pieces)
}


# STOP AND INSPECT:
# - This section answers: how close are we to file_check?
# - The practical acceptance rule used here is the relative ratio threshold.
# ============================================================
# SECTION 13: Plot residualized interaction
# ============================================================

resolve_environment_name <- function(requested, columns) {
  cols <- as.character(columns)
  if (requested %in% cols) return(requested)
  req <- tolower(requested)
  exact <- cols[tolower(cols) == req]
  if (length(exact) == 1) return(exact)
  partial <- cols[grepl(req, tolower(cols), fixed = TRUE)]
  if (length(partial) == 1) return(partial)
  stop(sprintf("Could not uniquely resolve environment: %s", requested))
}

fit_plot_model <- function(y, g, env) {
  X <- cbind(1, g, env, g * env)
  XtX_inv <- solve(crossprod(X))
  beta <- XtX_inv %*% crossprod(X, y)
  fitted <- as.numeric(X %*% beta)
  resid <- y - fitted
  dof <- nrow(X) - ncol(X)
  sigma2 <- sum(resid^2) / dof
  list(beta = as.numeric(beta), fitted = fitted, XtX_inv = XtX_inv, sigma2 = sigma2)
}

make_interaction_plot <- function(
    residualized_df,
    geno_s,
    env_name,
    gene_name,
    out_path,
    g_pval = NA_real_,
    e_pval = NA_real_,
    gxe_pval = NA_real_,
    conf_level = 0.95,
    phenotype_type
) {
  ###manual setting for troubleshooting manually####
  #residualized_df=residualized_df
  #geno_s=geno_s
  #env_name=env_name
  #gene_name=phenotype_id_resolved
  #out_path=plot_out
  #g_pval = g_pval
  #e_pval = e_pval
  #gxe_pval = gxe_pval
  #conf_level = 0.95
  #phenotype_type=phenotype_type
  #phenotype_type
  ##done manual annotation of stuff
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Install with install.packages('ggplot2').")
  }
  
  y <- residualized_df$phenotype_centered_residualized
  y
  env <- residualized_df[[paste0(env_name, "_centered_residualized")]]
  env
  g <- as.numeric(geno_s)
  
  fit <- fit_plot_model(y, g, env)
  beta <- fit$beta
  XtX_inv <- fit$XtX_inv
  sigma2 <- fit$sigma2
  g_cls <- pmax(0, pmin(2, round(g)))
  cols <- c("0" = "#1f77b4", "1" = "#ff7f0e", "2" = "#2ca02c")
  
  plot_df <- data.frame(
    env = as.numeric(env),
    y = as.numeric(y),
    genotype = factor(g_cls, levels = c(0, 1, 2))
  )
  
  env_grid <- seq(min(env, na.rm = TRUE), max(env, na.rm = TRUE), length.out = 120)
  alpha <- 1 - conf_level
  tcrit <- qt(1 - alpha / 2, df = max(1, length(y) - 4))
  line_df <- do.call(rbind, lapply(0:2, function(gv) {
    xmat <- cbind(1, gv, env_grid, gv * env_grid)
    yhat <- as.numeric(xmat %*% beta)
    se_fit <- sqrt(rowSums((xmat %*% XtX_inv) * xmat) * sigma2)
    data.frame(
      env = env_grid,
      yhat = yhat,
      lower = yhat - tcrit * se_fit,
      upper = yhat + tcrit * se_fit,
      genotype = factor(gv, levels = c(0, 1, 2))
    )
  }))
  
  gene_name
  # Example Usage:
  if (phenotype_type=="Expression"){
    target_id <- gene_name
    my_symbol <- get_symbol(target_id, gene_map)
    cat("The symbol for", target_id, "is", my_symbol, "\n")
    ome="Expression"
  }else {
    my_symbol<-gene_name
    ome="Methylation"
    
  }

  env_name
  # 1. Detect if it's Metabolic or Proteomic based on the string contents
  mod_type <- ifelse(grepl("Metabolite", env_name, ignore.case = TRUE), "Metabolic",
                     ifelse(grepl("Proteome", env_name, ignore.case = TRUE), "Proteomic", "Unknown"))
  
  # 2. Extract the raw color (grabs everything after 'ME' and before the first '_')
  color_raw <- sub("^ME([^_]+)_.*", "\\1", env_name)
  
  # 3. Capitalize the first letter of the extracted color
  color_cap <- paste0(toupper(substr(color_raw, 1, 1)), substring(color_raw, 2))
  
  # 4. Construct the final logical name
  env_name_better <- paste("Residualized", mod_type, "Module", color_cap)
  env_name_better
  
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = env, y = y, color = genotype)) +
    ggplot2::geom_point(alpha = 0.65, size = 1.7) +
    ggplot2::geom_ribbon(
      data = line_df,
      ggplot2::aes(x = env, ymin = lower, ymax = upper, fill = genotype, group = genotype),
      alpha = 0.18,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_line(
      data = line_df,
      ggplot2::aes(x = env, y = yhat, color = genotype, group = genotype),
      linewidth = 1.05,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = c("0" = "#1f77b4", "1" = "#ff7f0e", "2" = "#2ca02c"),
      breaks = c("0", "1", "2"),
      labels = c("C/C", "C/T", "T/T"),
      name = "Genotype"
    ) +
    ggplot2::scale_fill_manual(
      values = c("0" = "#1f77b4", "1" = "#ff7f0e", "2" = "#2ca02c"),
      breaks = c("0", "1", "2"),
      labels = c("C/C", "C/T", "T/T"),
      guide = "none"
    ) +
    ggplot2::labs(
      x = paste0(env_name_better),
      y = paste0("Residualized ", my_symbol," ", ome)
    ) +
    ggplot2::theme_bw()
  
  if (is.finite(g_pval) || is.finite(e_pval) || is.finite(gxe_pval)) {
    pval_label <- sprintf(
      "G:P.value=%s \nE:P.value=%s \nGxE:P.value=%s",
      if (is.finite(g_pval)) formatC(g_pval, format = "e", digits = 3) else "NA",
      if (is.finite(e_pval)) formatC(e_pval, format = "e", digits = 3) else "NA",
      if (is.finite(gxe_pval)) formatC(gxe_pval, format = "e", digits = 3) else "NA"
    )
    p <- p + ggplot2::annotate(
      "text",
      x = min(env, na.rm = TRUE),
      y = max(y, na.rm = TRUE),
      label = pval_label,
      hjust = 0,
      vjust = 1,
      size = 3.3
    )
  }
  
  p
}
plot_environment
if (plot_environment=="g"){
  env_name=str_replace_all(names(residualized_df)[4],"_centered_residualized","")
} else{
  env_name <- resolve_environment_name(plot_environment, colnames(interaction_df))
  
}
pval_col <- paste0("pval_g-", env_name)
interaction_pval <- if (pval_col %in% colnames(stats_df)) as.numeric(stats_df[[pval_col]][1]) else NA_real_
g_pval <- if ("pval_g" %in% colnames(stats_df)) as.numeric(stats_df[["pval_g"]][1]) else NA_real_
e_pval_col <- paste0("pval_", env_name)
e_pval <- if (e_pval_col %in% colnames(stats_df)) as.numeric(stats_df[[e_pval_col]][1]) else NA_real_
gxe_pval <- interaction_pval

dir_output<-"/labs/smontgom/grps_smontgom/hkrupkin/GxE/plots/single_variant_examples"
phenotype_type<-ifelse(grepl("ENSG",phenotype_id_resolved),"Expression","Methylation")
phenotype_type
path_output_file=paste0(dir_output,"/",cell_clean,"_",env_name,"_",target$effect_type,"_",phenotype_type,
                        "_",var_res$variant_id,".png")
path_output_file

interaction_plot <- make_interaction_plot(
  residualized_df,
  geno_s,
  env_name,
  phenotype_id_resolved,
  plot_out,
  g_pval = g_pval,
  e_pval = e_pval,
  gxe_pval = gxe_pval,
  conf_level = 0.95,
  phenotype_type
)
print(interaction_plot)

print("we are saving to : ")
print(path_output_file)
ggsave(plot=interaction_plot,
       file=path_output_file,
       width=5,
       height=5,
       dpi=900)

# STOP AND INSPECT:
# - The lines come from the residualized model-scale variables.
# - The interaction p-value displayed should match stats_df.

}
