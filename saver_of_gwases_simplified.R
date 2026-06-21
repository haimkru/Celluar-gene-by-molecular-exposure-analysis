### GWAS_PreProcessor.R ###
library(data.table)
library(GenomicRanges)
args <- commandArgs(trailingOnly = TRUE)
library(stringr)
library(dplyr)

# --- CONFIG (Change these for each GWAS) ---
if (length(args) < 1) {
  stop("Usage: script.R <GWAS_path> <GxE_path> <output_file_directory>")
}

# Assign arguments
GWAS_RAW <- args[1]
#GWAS_NAME <- "LDL"
#GWAS_RAW  <- "/oak/stanford/groups/smontgom/shared/gwas_summary_stats/barbeira_gtex_imputed/imputed_gwas_hg38_1.1/imputed_MAGNETIC_LDL.C.txt.gz"

GWAS_NAME<-str_replace_all(basename(GWAS_RAW),"txt.gz","")
OUT_DIR   <- "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/processed_rds_gwas/"
#dir.create(OUT_DIR)
# 1. Load data
print(paste("Processing GWAS:", GWAS_NAME))
gwas <- fread(GWAS_RAW)

# 2. Identify Windows around GWAS peaks
print("Finding GWAS peaks...")
# Standard GWAS significance is 5e-8
peaks <- gwas[pvalue < 5e-8]
ranges <- reduce(GRanges(
  seqnames = peaks$chromosome,
  ranges = IRanges(start = peaks$position - 500000, 
                   end = peaks$position + 500000)
))

# 3. Keep all SNPs in these active regions
print("Filtering to active loci...")
gwas_gr <- GRanges(gwas$chromosome, IRanges(gwas$position, gwas$position))
gwas_hits <- findOverlaps(gwas_gr, ranges)
gwas_active <- gwas[queryHits(gwas_hits)]

# 4. Save
print("Indexing and saving...")
# Keying by chr/pos makes regional lookups instant in the coloc step
setkey(gwas_active, chromosome, position) 
# Also key by ID for fast merging
setkey(gwas_active, panel_variant_id) 
path_saving_gwas<-paste0(OUT_DIR, GWAS_NAME, "_active_loci.rds")
path_saving_gwas
saveRDS(gwas_active,path_saving_gwas )
print(paste("Done! GWAS shrunk to", nrow(gwas_active), "rows."))
