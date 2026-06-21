### GxE_PreProcessor.R ###
library(data.table)
library(stringr)
library(GenomicRanges)

# --- CONFIG ---
# Using 'my_args' to avoid conflict with the built-in 'args()' function
my_args <- commandArgs(trailingOnly = TRUE)

# Debug: Print the arguments to the log so you can see what Bash passed
print("Arguments received:")
print(my_args)

if (length(my_args) < 1) {
  stop("Error: No GxE file path provided to the script.")
}

GxE_path <- my_args[1]
GxE_name <- str_replace_all(basename(GxE_path), "_trans_qtl.tsv.gz", "")
OUT_DIR  <- "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/processed_rds_GxEs_simplfiied/"

if(!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# 1. Load data
print(paste("Processing GxE:", GxE_name))
if (!file.exists(GxE_path)) stop(paste("File not found:", GxE_path))

GxE <- fread(GxE_path)

# --- DYNAMIC COLUMN CHECK ---
p_col <- if ("pval_ME" %in% names(GxE)) "pval_ME" else "pval"
print(paste("Using column for filtering:", p_col))

# 2. Identify Windows
print("Finding peaks and defining windows...")
# Use get() to handle the dynamic column name string
peaks <- GxE[get(p_col) < 5e-5] 

if (nrow(peaks) == 0) {
  print("No peaks found below threshold. Saving empty result to prevent array failure.")
  # Create an empty version of the expected output format if necessary
  saveRDS(GxE[0], paste0(OUT_DIR, GxE_name, "_active_loci.rds"))
  quit(save = "no", status = 0)
}

# Now split strings ONLY for the peaks
peaks[, c("chr", "pos") := tstrsplit(variant_id, "_", keep=1:2)]
peaks[, pos := as.numeric(pos)]

# Create merged 1MB windows
ranges <- reduce(GRanges(
  seqnames = peaks$chr,
  ranges = IRanges(start = peaks$pos - 500000, end = peaks$pos + 500000)
))

# 3. Filter the main table
print("Filtering rows to active loci...")
GxE[, c("chr", "pos") := tstrsplit(variant_id, "_", keep=1:2)]
GxE[, pos := as.numeric(pos)]

ranges_dt <- as.data.table(ranges)
setnames(ranges_dt, c("seqnames", "start", "end"), c("chr", "win_start", "win_end"))
ranges_dt[, chr := as.character(chr)]

# Set keys for foverlaps
setkey(ranges_dt, chr, win_start, win_end)
GxE[, pos_end := pos]
setkey(GxE, chr, pos, pos_end)

# Run the join
GxE_active <- foverlaps(
  GxE, 
  ranges_dt, 
  by.x = c("chr", "pos", "pos_end"), 
  by.y = c("chr", "win_start", "win_end"), 
  nomatch = 0
)

# Clean up
GxE_active[, pos_end := NULL]
GxE_active[, c("win_start", "win_end") := NULL]

# 4. Save
print("Indexing and saving...")
setkey(GxE_active, variant_id)

path_saving <- paste0(OUT_DIR, GxE_name, "_active_loci.rds")
saveRDS(GxE_active, path_saving)
print(paste("Successfully saved to:", path_saving))