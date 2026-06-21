#author: haim krupkin
#date updated: 05/15/2026


###packages#####
.libPaths(c("/scg/apps/software/r/4.3.3/lib", .libPaths()))
library(data.table)
library(stringr)
library(GenomicRanges)

###analysis####
#/labs/smontgom/grps_smontgom/hkrupkin/GxE/Code/saver_of_GxE_simplifiedV2.R


args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript combined_GxE_to_per_gene_active_loci.R <GxE_input.{tsv.gz|tsv|rds}> [output_dir] [peak_pval_threshold]")
}

GxE_path <- args[1]
#before 05/15/2026
#GxE_path <- "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/for_developments/output_bigs/nk_expression_MEcyan_trans_qtl.tsv.gz"
#after 05/15/2026
#GxE_path <- "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/for_developments/output_bigsV2/nk_metabolome_chunk_1_expression_MEcyan_trans_qtl.tsv.gz"
#GxE_path <- "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/for_developments/output_bigsV2/tmp_metabolome_chunk_1_expression_MEcyan_trans_qtl.tsv.gz"

#before 05/15/2026 
#output_dir <- if (length(args) >= 2) args[2] else "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/saving_for_coloc_susieV2/"
#after 05/15/2026
output_dir <- if (length(args) >= 2) args[2] else "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/saving_for_coloc_susieV3/"



PEAK_PVAL_THRESHOLD <- if (length(args) >= 3) as.numeric(args[3]) else 5e-5
WINDOW_HALF_SIZE <- 500000

if (is.na(PEAK_PVAL_THRESHOLD)) {
  stop("peak_pval_threshold must be numeric")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("[START] Combined GxE preprocessing")
message("[INPUT] ", GxE_path)
message("[OUTPUT_DIR] ", output_dir)
message("[THRESHOLD] ", PEAK_PVAL_THRESHOLD)

if (!file.exists(GxE_path)) {
  stop("Input file not found: ", GxE_path)
}

input_name <- basename(GxE_path)
GxE_name <- input_name
GxE_name <- str_replace(GxE_name, "_trans_qtl\\.tsv\\.gz$", "")
GxE_name <- str_replace(GxE_name, "\\.tsv\\.gz$", "")
GxE_name <- str_replace(GxE_name, "\\.tsv$", "")
GxE_name <- str_replace(GxE_name, "\\.rds$", "")

if (grepl("\\.rds$", GxE_path, ignore.case = TRUE)) {
  message("[LOAD] Reading RDS input")
  GxE <- as.data.table(readRDS(GxE_path))
} else {
  message("[LOAD] Reading TSV input")
  GxE <- fread(GxE_path)
}

message("[LOAD] Rows loaded: ", nrow(GxE))
message("[LOAD] Columns loaded: ", ncol(GxE))

required_cols <- c("variant_id", "phenotype_id")
missing_required <- setdiff(required_cols, names(GxE))
if (length(missing_required) > 0) {
  stop("Missing required columns: ", paste(missing_required, collapse = ", "))
}

peak_cols <- intersect(c("pval_ME", "pval", "pval_g"), names(GxE))
if (length(peak_cols) == 0) {
  stop("None of the supported p-value columns were found: pval_ME, pval, pval_g")
}
message("[FILTER] Supported p-value columns: ", paste(peak_cols, collapse = ", "))

peak_mask <- Reduce(`|`, lapply(peak_cols, function(col_name) {
  vals <- GxE[[col_name]]
  !is.na(vals) & vals < PEAK_PVAL_THRESHOLD
}))

peaks <- copy(GxE[peak_mask])
message("[FILTER] Rows passing threshold in any supported column: ", nrow(peaks))

if (nrow(peaks) == 0) {
  marker_path <- file.path(output_dir, paste0(GxE_name, "_NO_PEAKS.rds"))
  saveRDS(GxE[0], marker_path)
  message("[SKIP] No qualifying peaks found. Saved empty marker: ", marker_path)
  quit(save = "no", status = 0)
}

message("[PARSE] Parsing variant positions for full table")
GxE[, c("chr", "pos") := tstrsplit(variant_id, "_", keep = 1:2)]
GxE[, pos := as.numeric(pos)]
GxE <- GxE[!is.na(chr) & !is.na(pos)]
message("[PARSE] Rows with valid chr/pos in full table: ", nrow(GxE))

message("[PARSE] Parsing variant positions for peak table")
peaks[, c("chr", "pos") := tstrsplit(variant_id, "_", keep = 1:2)]
peaks[, pos := as.numeric(pos)]
peaks <- peaks[!is.na(chr) & !is.na(pos)]
message("[PARSE] Rows with valid chr/pos in peak table: ", nrow(peaks))

list_of_genes <- unique(peaks$phenotype_id)
message("[GENES] Genes with qualifying peaks: ", length(list_of_genes))

saved_count <- 0L
skipped_count <- 0L

list_of_genes[1:5]

setkey(GxE, phenotype_id)
setkey(peaks, phenotype_id)

for (idx in seq_along(list_of_genes)) { #forward order
#for (idx in rev(seq_along(list_of_genes))) { #reverse order  
#for (idx in rev(6694:length(list_of_genes))) { #reverse order whats left
  gene <- list_of_genes[[idx]]
  message("[GENE_START] ", idx, "/", length(list_of_genes), " :: ", gene)
  
  gene <- list_of_genes[[idx]]
  message("[GENE_START] ", idx, "/", length(list_of_genes), " :: ", gene)
  
  # FAST SUBSET: Uses binary search via the .() syntax
  gene_dt <- copy(GxE[.(gene)])
  gene_peaks <- peaks[.(gene)]
  
  if (nrow(gene_dt) == 0) {
    message("[GENE_SKIP] No rows found in full table for gene: ", gene)
    skipped_count <- skipped_count + 1L
    next
  }
  
  if (nrow(gene_peaks) == 0) {
    message("[GENE_SKIP] No qualifying peaks found for gene: ", gene)
    skipped_count <- skipped_count + 1L
    next
  }
  
  ranges <- reduce(GRanges(
    seqnames = gene_peaks$chr,
    ranges = IRanges(
      start = pmax(1, gene_peaks$pos - WINDOW_HALF_SIZE),
      end = gene_peaks$pos + WINDOW_HALF_SIZE
    )
  ))
  
  ranges_dt <- as.data.table(ranges)
  setnames(ranges_dt, c("seqnames", "start", "end"), c("chr", "win_start", "win_end"))
  ranges_dt[, chr := as.character(chr)]
  
  gene_dt[, pos_end := pos]
  setkey(ranges_dt, chr, win_start, win_end)
  setkey(gene_dt, chr, pos, pos_end)
  
  gene_active <- foverlaps(
    gene_dt,
    ranges_dt,
    by.x = c("chr", "pos", "pos_end"),
    by.y = c("chr", "win_start", "win_end"),
    nomatch = 0
  )
  
  gene_active[, pos_end := NULL]
  gene_active[, c("win_start", "win_end") := NULL]
  
  if (nrow(gene_active) == 0) {
    message("[GENE_SKIP] No rows remained after active-locus overlap for gene: ", gene)
    skipped_count <- skipped_count + 1L
    next
  }
  
  setkey(gene_active, variant_id)
  path_saving <- file.path(output_dir, paste0(GxE_name, "_", gene, "_active_loci.rds"))
  saveRDS(gene_active, path_saving)
  saved_count <- saved_count + 1L
  
  message("[GENE_SAVE] Saved ", nrow(gene_active), " rows to: ", path_saving)
}

message("[DONE] Genes saved: ", saved_count)
message("[DONE] Genes skipped: ", skipped_count)
message("[DONE] Combined GxE preprocessing finished")
