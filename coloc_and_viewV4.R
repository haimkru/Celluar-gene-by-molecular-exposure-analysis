###############################################################################
# coloc_GxE_and_view.R
# Author: Haim Krupkin
# Date: 04/12/2026  |  Fixed: 05/22/2026
#
# FIX LOG (05/22/2026):
#   1. eqtl_dir: removed erroneous "http://" prefix
#   2. Removed dplyr pipe — replaced with data.table syntax (dplyr not loaded)
#   3. resolve_gxe_rds_path: completely rewritten to parse filename column correctly
#      and construct path without list.files() or nonexistent seed_row$chunk_id
#   4. resolve_coloc_gene_rds_path: rewritten — no more list.files() on 14M-file dir
#   5. resolve_eqtl_rds_path: rewritten — constructs path directly, no list.files()
#   6. plot_coloc_hit: parses all 5 prefix parts from gxe_path (was only parsing 2),
#      passes file_prefix to resolve_eqtl_rds_path instead of cell/modality alone
#
# Filename convention (derived from actual files on disk):
#   saved RDS: {cell}_{omics}_{chunk}_{N}_{assay}_{term}_{gene}_active_loci.rds
#   e.g.  bcell_metabolome_chunk_1_expression_MEdarkorange_ENSG00000268903.1_active_loci.rds
#
#   filename col in all_hits TSV: {omics}_{chunk}_{N}_{cell}_{assay}_sig.tsv
#   e.g.  metabolome_chunk_1_bcell_expression_sig.tsv
#
#   term col: MEdarkorange_Eigen  → strip _Eigen → MEdarkorange
#             g                  → map to       → main_effect
###############################################################################

###############################################################################
# CONFIG
###############################################################################

SINGLE_MIN_PP_H4     <- 0.50
SINGLE_MIN_P_INT     <- 1e-5
SINGLE_MIN_P_GWAS    <- 1e-5
COLOC_REQUIRE_TERM_MATCH <- FALSE

gwas_dir     <- "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/processed_rds_gwas"
output_dir   <- "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/coloc_interaction_with_GWASesV3"
all_hits_tsv <- "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/all_cell_types_sig_hits.tsv"
# FIX 1: was "http://labs/smontgom/..." — corrected to absolute path
eqtl_dir     <- "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/saving_for_coloc_susieV3/"

INTERACTION_PVAL_THRESH <- 1e-5
GWAS_PVAL_THRESH        <- 1e-5
WINDOW_SIZE             <- 200000
GWAS_NEAR_WINDOW        <- 100000
MHC_CHR                 <- "chr6"
MHC_START               <- 25726063
MHC_END                 <- 33400644
MIN_VARIANTS            <- 50
SAMPLE_N_GXE            <- 541

CC_GWAS <- list()

###############################################################################
# PACKAGES
###############################################################################

library(data.table)
library(GenomicRanges)
library(coloc)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
term_arg <- NA_character_
if (length(args) > 0) {
  term_hits <- grep("^--term=", args, value = TRUE)
  if (length(term_hits) > 0) {
    term_arg <- sub("^--term=", "", term_hits[1])
    term_arg <- trimws(term_arg)
    if (!nzchar(term_arg)) term_arg <- NA_character_
  }
}

coloc_results <- data.table()
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

###############################################################################
# HELPERS
###############################################################################

parse_variant_chr_pos <- function(variant_ids) {
  m <- str_match(variant_ids, "^(chr[^_]+)_([0-9]+)_")
  data.table(
    variant_id = variant_ids,
    chromosome = m[, 2],
    position   = as.numeric(m[, 3])
  )
}

is_in_mhc <- function(chromosome, position) {
  chromosome == MHC_CHR & !is.na(position) &
    position >= MHC_START & position <= MHC_END
}

# ---------------------------------------------------------------------------
# INTERNAL HELPER: parse filename column into an ordered prefix string
#
# filename col: "metabolome_chunk_1_bcell_expression_sig.tsv"
# parts after stripping _sig.tsv: metabolome | chunk | 1 | bcell | expression
# RDS prefix order:               bcell_metabolome_chunk_1_expression
# ---------------------------------------------------------------------------
.parse_filename_prefix <- function(filename) {
  fn_base <- sub("_sig\\.tsv$", "", as.character(filename))
  parts   <- strsplit(fn_base, "_")[[1]]
  # Guard: need at least 5 tokens (omics, "chunk", N, cell, assay)
  if (length(parts) < 5) {
    warning("Unexpected filename format — cannot parse prefix: ", filename)
    return(NA_character_)
  }
  omics <- parts[1]
  chunk <- paste(parts[2], parts[3], sep = "_")  # "chunk_1"
  cell  <- parts[4]
  assay <- tolower(parts[5])
  paste(cell, omics, chunk, assay, sep = "_")    # "bcell_metabolome_chunk_1_expression"
}

# ---------------------------------------------------------------------------
# INTERNAL HELPER: map raw term value to the token used in RDS filenames
# ---------------------------------------------------------------------------
.clean_term <- function(term_raw) {
  term_raw <- as.character(term_raw)
  if (term_raw == "g") return("main_effect")
  sub("_Eigen$", "", term_raw)
}

# ---------------------------------------------------------------------------
# FIX 3: resolve_gxe_rds_path
#
# Constructs the path to the per-gene active-loci RDS from seed_row columns:
#   filename, term, phenotype_id
# No list.files(), no nonexistent seed_row$chunk_id.
# ---------------------------------------------------------------------------
resolve_gxe_rds_path <- function(seed_row) {
  prefix     <- .parse_filename_prefix(seed_row$filename)
  if (is.na(prefix)) return(NA_character_)
  
  term_clean <- .clean_term(seed_row$term)
  gene       <- as.character(seed_row$phenotype_id)
  
  fname <- paste0(prefix, "_", term_clean, "_", gene, "_active_loci.rds")
  path  <- file.path(eqtl_dir, fname)
  cat("      [resolve_gxe] Constructed path:", path, "\n")
  return(path)
}

# ---------------------------------------------------------------------------
# FIX 4: resolve_coloc_gene_rds_path
#
# Was scanning the entire 14M-file directory with list.files() on every call.
# Now constructs the path directly — same convention as resolve_gxe_rds_path.
# Falls back to main_effect if the interaction file is absent.
# ---------------------------------------------------------------------------
resolve_coloc_gene_rds_path <- function(seed_row, gxe_path_hint = NA_character_) {
  # Primary: use the hint path if it already exists on disk
  if (!is.na(gxe_path_hint) && file.exists(gxe_path_hint)) {
    return(gxe_path_hint)
  }
  
  # Construct interaction path from seed columns
  primary <- resolve_gxe_rds_path(seed_row)
  if (!is.na(primary) && file.exists(primary)) return(primary)
  
  # Fallback: main_effect file for the same prefix + gene
  if (!is.na(seed_row$filename)) {
    prefix <- .parse_filename_prefix(seed_row$filename)
    gene   <- as.character(seed_row$phenotype_id)
    me_fname <- paste0(prefix, "_main_effect_", gene, "_active_loci.rds")
    me_path  <- file.path(eqtl_dir, me_fname)
    if (file.exists(me_path)) {
      cat("      [resolve_coloc] Falling back to main_effect RDS:", me_path, "\n")
      return(me_path)
    }
  }
  
  cat("      [resolve_coloc] No RDS found for gene:", seed_row$phenotype_id, "\n")
  return(NA_character_)
}

# ---------------------------------------------------------------------------
# FIX 5: resolve_eqtl_rds_path
#
# Was scanning the entire directory with list.files() for every plot call.
# New signature: takes file_prefix (e.g. "bcell_metabolome_chunk_1_expression")
# and gene_id, constructs the main_effect path directly.
# ---------------------------------------------------------------------------
resolve_eqtl_rds_path <- function(file_prefix, gene_id) {
  fname <- paste0(file_prefix, "_main_effect_", gene_id, "_active_loci.rds")
  path  <- file.path(eqtl_dir, fname)
  if (file.exists(path)) {
    cat("      [resolve_eqtl] Found ME RDS:", path, "\n")
    return(path)
  }
  cat("      [resolve_eqtl] No main_effect RDS found for gene:", gene_id,
      "prefix:", file_prefix, "\n")
  return(NA_character_)
}

###############################################################################
# COLOC RUNNER
###############################################################################

run_coloc_seed <- function(gwas, gwas_basename, seed_row, gxe_path,
                           gwas_is_cc = FALSE) {
  
  append_coloc_result <- function(status, msg = NA, extra = list()) {
    row <- data.table(
      gene       = as.character(seed_row$phenotype_id),
      variant    = as.character(seed_row$variant_id),
      gwas       = gwas_basename,
      gxe_rds    = as.character(gxe_path),
      status     = status,
      msg        = ifelse(is.na(msg), "", as.character(msg)),
      PP.H0      = ifelse(!is.null(extra$PP.H0),      extra$PP.H0,      NA),
      PP.H4      = ifelse(!is.null(extra$PP.H4),      extra$PP.H4,      NA),
      nsnps      = ifelse(!is.null(extra$nsnps),      extra$nsnps,      NA),
      Lead_P_int = ifelse(!is.null(extra$Lead_P_int), extra$Lead_P_int, NA),
      Lead_P_GWAS= ifelse(!is.null(extra$Lead_P_GWAS),extra$Lead_P_GWAS,NA)
    )
    assign("coloc_results",
           rbind(get("coloc_results", envir = .GlobalEnv), row),
           envir = .GlobalEnv)
  }
  
  cat("\n--- Running seed coloc for ", seed_row$phenotype_id,
      " @ ", seed_row$variant_id, " ---\n", sep = "")
  
  if (is.na(gxe_path) || !file.exists(gxe_path)) {
    cat("  GxE RDS not found:", gxe_path, "\n")
    append_coloc_result("fail", "GxE RDS file not found")
    return(NULL)
  }
  
  gxe_dt      <- as.data.table(readRDS(gxe_path))
  needed      <- c("variant_id", "phenotype_id", "pval", "b_int", "b_int_se")
  missing_cols <- setdiff(needed, names(gxe_dt))
  if (length(missing_cols) > 0) {
    cat("  MISSING COLUMNS in GxE RDS:", paste(missing_cols, collapse = ", "), "\n")
    append_coloc_result("fail", paste("Missing columns:", paste(missing_cols, collapse = ",")))
    rm(gxe_dt); gc(verbose = FALSE)
    return(NULL)
  }
  
  lead_variant <- seed_row$variant_id
  lead_gene    <- seed_row$phenotype_id
  lead_chr     <- seed_row$chromosome
  lead_pos     <- seed_row$position
  
  region_start <- max(1, lead_pos - WINDOW_SIZE)
  region_end   <- lead_pos + WINDOW_SIZE
  
  cat("  Region:", lead_chr, region_start, region_end, "\n")
  
  gene_gxe <- gxe_dt[phenotype_id == lead_gene]
  if (nrow(gene_gxe) == 0) {
    cat("  Gene not found in selected GxE RDS.\n")
    append_coloc_result("fail", "Gene not found in GxE RDS")
    rm(gxe_dt); gc(verbose = FALSE)
    return(NULL)
  }
  
  seed_in_rds <- gene_gxe[variant_id == lead_variant][1]
  if (nrow(seed_in_rds) == 0) {
    cat("  Seed variant not found in selected GxE RDS; skipping.\n")
    append_coloc_result("fail", "Seed variant not found in GxE RDS")
    rm(gxe_dt); gc(verbose = FALSE)
    return(NULL)
  }
  if (is.na(seed_in_rds$pval) || seed_in_rds$pval > SINGLE_MIN_P_INT) {
    cat(sprintf("  Seed not interaction-significant (p_int=%.2e > %.2e); skipping.\n",
                seed_in_rds$pval, SINGLE_MIN_P_INT))
    append_coloc_result("fail",
                        sprintf("Seed not interaction-significant (p_int=%.2e)", seed_in_rds$pval))
    rm(gxe_dt); gc(verbose = FALSE)
    return(NULL)
  }
  
  gwas_region <- gwas[chromosome == lead_chr &
                        position %between% c(region_start, region_end)]
  if (nrow(gwas_region) < MIN_VARIANTS) {
    cat("  Not enough GWAS variants in window (", nrow(gwas_region), ").\n", sep = "")
    append_coloc_result("fail",
                        sprintf("Not enough GWAS variants in window (%d)", nrow(gwas_region)))
    rm(gxe_dt); gc(verbose = FALSE)
    return(NULL)
  }
  
  gwas_lead <- gwas_region[!is.na(pvalue)][order(pvalue)][1]
  if (nrow(gwas_lead) == 0)
    gwas_lead <- data.table(panel_variant_id = NA_character_,
                            pvalue = NA_real_, position = NA_real_)
  
  merged <- gene_gxe[gwas_region, on = .(variant_id = panel_variant_id), nomatch = NULL]
  if (nrow(merged) < MIN_VARIANTS) {
    cat("  Not enough overlapping variants for coloc (", nrow(merged), ").\n", sep = "")
    append_coloc_result("fail",
                        sprintf("Not enough overlapping variants for coloc (%d)", nrow(merged)))
    rm(gxe_dt); gc(verbose = FALSE)
    return(NULL)
  }
  
  merged <- merged[!is.na(b_int) & !is.na(b_int_se) & !is.na(pvalue) &
                     !is.na(frequency) & !is.na(variant_id)]
  merged[pvalue <= 0, pvalue := 1e-300]
  
  if (nrow(merged) < MIN_VARIANTS) {
    cat("  Not enough non-missing overlapping variants (", nrow(merged), ").\n", sep = "")
    append_coloc_result("fail",
                        sprintf("Not enough non-missing overlapping variants (%d)", nrow(merged)))
    rm(gxe_dt); gc(verbose = FALSE)
    return(NULL)
  }
  
  merged[, maf_gwas := ifelse(frequency > 0.5, 1 - frequency, frequency)]
  
  # Check for af column (used as MAF for GxE dataset)
  if (!"af" %in% names(merged)) {
    cat("  WARNING: 'af' column missing from merged table; falling back to maf_gwas.\n")
    merged[, af := maf_gwas]
  }
  
  dataset_gxe <- list(
    beta    = merged$b_int,
    varbeta = merged$b_int_se^2,
    type    = "quant",
    snp     = merged$variant_id,
    N       = SAMPLE_N_GXE,
    MAF     = merged$af
  )
  
  if (gwas_is_cc) {
    dataset_gwas <- list(
      pvalues = merged$pvalue,
      N       = merged$sample_size[1],
      MAF     = merged$maf_gwas,
      type    = "cc",
      s       = CC_GWAS[[gwas_basename]],
      snp     = merged$variant_id
    )
  } else {
    dataset_gwas <- list(
      pvalues = merged$pvalue,
      N       = merged$sample_size[1],
      MAF     = merged$maf_gwas,
      type    = "quant",
      snp     = merged$variant_id
    )
  }
  
  coloc_res <- tryCatch({
    junk <- capture.output(
      coloc_call <- coloc.abf(
        dataset1 = dataset_gxe,
        dataset2 = dataset_gwas,
        p1  = 1e-4,
        p2  = 1e-4,
        p12 = 5e-6
      )
    )
    coloc_call
  }, error = function(e) {
    cat("  coloc error:", e$message, "\n")
    append_coloc_result("fail", paste("coloc error:", e$message))
    NULL
  })
  
  if (is.null(coloc_res)) {
    rm(gxe_dt); gc(verbose = FALSE)
    return(NULL)
  }
  
  lead_row <- merged[variant_id == lead_variant & !is.na(pval) & !is.na(pvalue)]
  if (nrow(lead_row) == 0) {
    candidates <- merged[!is.na(pval) & !is.na(pvalue)]
    if (nrow(candidates) > 0) {
      lead_row <- candidates[which.min(pval)]
    } else {
      seed_gwas <- gwas[panel_variant_id == lead_variant & !is.na(pvalue), pvalue][1]
      lead_row  <- data.table(
        variant_id = lead_variant,
        pval       = seed_in_rds$pval,
        pvalue     = seed_gwas,
        position   = lead_pos
      )
    }
  } else {
    lead_row <- lead_row[1]
  }
  
  out <- data.table(
    Region      = paste(lead_chr, region_start, region_end, sep = "_"),
    Gene        = lead_gene,
    nsnps       = coloc_res$summary["nsnps"],
    PP.H0       = coloc_res$summary["PP.H0.abf"],
    PP.H1       = coloc_res$summary["PP.H1.abf"],
    PP.H2       = coloc_res$summary["PP.H2.abf"],
    PP.H3       = coloc_res$summary["PP.H3.abf"],
    PP.H4       = coloc_res$summary["PP.H4.abf"],
    Lead_SNP    = lead_row$variant_id,
    Lead_P_int  = lead_row$pval,
    Lead_P_GWAS = lead_row$pvalue,
    Lead_pos    = lead_row$position,
    GWAS_Lead_SNP  = gwas_lead$panel_variant_id,
    GWAS_Lead_P    = gwas_lead$pvalue,
    GWAS_Lead_pos  = gwas_lead$position,
    Lead_SNP_matches_GWAS_lead = as.logical(
      !is.na(lead_row$variant_id) && !is.na(gwas_lead$panel_variant_id) &&
        lead_row$variant_id == gwas_lead$panel_variant_id),
    Seed_term  = if ("term" %in% names(seed_row)) seed_row$term else NA_character_,
    gwas_name  = gwas_basename,
    gxe_name   = str_replace(basename(gxe_path), "\\.rds$", "")
  )
  
  cat(sprintf("  coloc done: PP.H0=%.3f | PP.H4=%.3f | nsnps=%d | p_int=%.1e | p_gwas=%.1e\n",
              out$PP.H0, out$PP.H4, out$nsnps, out$Lead_P_int, out$Lead_P_GWAS))
  append_coloc_result("success", NA,
                      list(PP.H0 = out$PP.H0, PP.H4 = out$PP.H4, nsnps = out$nsnps,
                           Lead_P_int = out$Lead_P_int, Lead_P_GWAS = out$Lead_P_GWAS))
  rm(gxe_dt); gc(verbose = FALSE)
  return(list(summary = out, raw_coloc = coloc_res))
}

###############################################################################
# STACKED PLOT FUNCTION
###############################################################################

plot_coloc_hit <- function(hit, gwas, gxe_path, plot_dir = output_dir) {
  
  library(ggplot2)
  
  cat("\n┌──────────────────────────────────────────────────────────┐\n")
  cat("│  COLOC HIT REPORT                                       │\n")
  cat("├──────────────────────────────────────────────────────────┤\n")
  cat("  Region:      ", hit$Region, "\n")
  cat("  Gene:        ", hit$Gene, "\n")
  cat("  GWAS:        ", hit$gwas_name, "\n")
  cat("  GxE:         ", hit$gxe_name, "\n")
  cat("  PP.H4:       ", sprintf("%.4f", hit$PP.H4), "\n")
  cat("  PP.H3:       ", sprintf("%.4f", hit$PP.H3), "\n")
  cat("  Lead SNP:    ", hit$Lead_SNP, "\n")
  cat("  Lead p(int): ", sprintf("%.2e", hit$Lead_P_int), "\n")
  cat("  Lead p(GWAS):", sprintf("%.2e", hit$Lead_P_GWAS), "\n")
  
  parts   <- strsplit(hit$Region, "_")[[1]]
  chr_r   <- parts[1]
  start_r <- as.numeric(parts[2])
  end_r   <- as.numeric(parts[3])
  gene_r  <- hit$Gene
  
  gwas_sub <- gwas[chromosome == chr_r & position >= start_r & position <= end_r]
  gwas_sub[pvalue <= 0 | is.na(pvalue), pvalue := 1e-300]
  gwas_sub[, neglog10p := -log10(pvalue)]
  cat("    GWAS variants:", nrow(gwas_sub), "\n")
  
  cat("  Loading GxE for plotting...\n")
  gxe_dt <- as.data.table(readRDS(gxe_path))
  
  avail_genes <- unique(gxe_dt$phenotype_id)
  if (!(gene_r %in% avail_genes)) {
    base_g  <- sub("\\..*", "", gene_r)
    partial <- grep(base_g, avail_genes, value = TRUE)
    if (length(partial) > 0) {
      cat("    Gene fallback:", partial[1], "\n")
      gene_r <- partial[1]
    }
  }
  
  gxe_gene <- gxe_dt[phenotype_id == gene_r]
  gxe_gene <- merge(
    gxe_gene,
    gwas[, .(panel_variant_id, chromosome, position)],
    by.x = "variant_id", by.y = "panel_variant_id",
    all.x = FALSE, suffixes = c("", ".gwas")
  )
  gxe_sub <- gxe_gene[chromosome == chr_r & position >= start_r & position <= end_r]
  gxe_sub[, neglog10p := -log10(pval)]
  cat("    GxE variants:", nrow(gxe_sub), "\n")
  cat("    Interaction min p:", sprintf("%.2e", min(gxe_sub$pval)), "\n")
  
  pval_cols <- grep("^pval", names(gxe_sub), value = TRUE)
  for (pc in pval_cols) {
    cat(sprintf("      %-12s min=%.2e  #<1e-5=%d\n",
                pc, min(gxe_sub[[pc]], na.rm = TRUE),
                sum(gxe_sub[[pc]] < 1e-5, na.rm = TRUE)))
  }
  rm(gxe_dt); gc(verbose = FALSE)
  
  # FIX 6: Parse all 5 prefix parts from gxe_path so we can build the ME path correctly.
  # gxe_path basename example:
  #   bcell_metabolome_chunk_1_expression_MEdarkorange_ENSG00000268903.1_active_loci.rds
  # After stripping .rds and trailing _active_loci:
  #   bcell_metabolome_chunk_1_expression_MEdarkorange_ENSG00000268903.1
  # parts_gxe: [1]=bcell [2]=metabolome [3]=chunk [4]=1 [5]=expression [6]=term [7..]=gene
  gxe_base   <- sub("\\.rds$", "", basename(gxe_path))
  gxe_base   <- sub("_active_loci$", "", gxe_base)
  parts_gxe  <- strsplit(gxe_base, "_")[[1]]
  
  cell_type   <- parts_gxe[1]   # e.g. "bcell"
  modality    <- parts_gxe[2]   # e.g. "metabolome"
  # file_prefix = cell_omics_chunk_N_assay — first 5 tokens
  file_prefix <- paste(parts_gxe[1:5], collapse = "_")
  # e.g. "bcell_metabolome_chunk_1_expression"
  
  preferred_term <- if ("Seed_term" %in% names(hit)) hit$Seed_term else NA_character_
  
  # FIX 5 in action: resolve_eqtl_rds_path now takes file_prefix + gene_id
  me_rds <- resolve_eqtl_rds_path(file_prefix, gene_r)
  me_ok  <- !is.na(me_rds) && file.exists(me_rds)
  me_sub <- NULL
  
  if (me_ok) {
    cat("  Loading main effect from:", me_rds, "\n")
    me_dt   <- as.data.table(readRDS(me_rds))
    me_gene <- me_dt[phenotype_id == gene_r]
    if (nrow(me_gene) == 0) {
      gene_base <- sub("\\..*", "", gene_r)
      me_gene   <- me_dt[grepl(paste0("^", gene_base, "(\\..*)?$"), phenotype_id)]
    }
    me_gene <- merge(
      me_gene,
      gwas[, .(panel_variant_id, chromosome, position)],
      by.x = "variant_id", by.y = "panel_variant_id",
      all.x = FALSE, suffixes = c("", ".gwas")
    )
    cat(sprintf("    [ME panel] overlap: %d variants\n", nrow(me_gene)))
    me_sub  <- me_gene[chromosome == chr_r & position >= start_r & position <= end_r]
    me_pcol <- intersect(c("pval", "pval_ME", "pvalue"), names(me_sub))[1]
    if (length(me_pcol) == 0 || is.na(me_pcol)) {
      cat("    Main effect file has no p-value column. Skipping ME panel.\n")
      me_sub <- NULL
    } else {
      me_sub[, neglog10p := -log10(get(me_pcol))]
      cat("    ME variants:", nrow(me_sub), "| col:", me_pcol, "\n")
    }
    rm(me_dt); gc(verbose = FALSE)
  } else {
    cat("  No main_effect RDS found for gene:", gene_r,
        "prefix:", file_prefix, "\n")
  }
  
  lead_pos <- hit$Lead_pos
  
  has_me <- me_ok && !is.null(me_sub) && nrow(me_sub) > 0
  
  ymax_vals <- c(gwas_sub$neglog10p, gxe_sub$neglog10p)
  if (has_me) ymax_vals <- c(ymax_vals, me_sub$neglog10p)
  shared_ymax <- max(ymax_vals, na.rm = TRUE) * 1.05
  
  th <- theme_bw() +
    theme(panel.grid.minor = element_blank(),
          plot.title  = element_text(size = 10, face = "bold"),
          axis.title  = element_text(size = 9),
          axis.text   = element_text(size = 8),
          plot.margin = margin(5, 10, 2, 10))
  
  p1 <- ggplot(gwas_sub, aes(position / 1e6, neglog10p)) +
    geom_point(size = 0.8, alpha = 0.6, color = "grey40") +
    geom_point(data = gwas_sub[abs(position - lead_pos) < 1],
               size = 2.5, color = "red", shape = 18) +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed",
               color = "red", linewidth = 0.3) +
    coord_cartesian(ylim = c(0, shared_ymax)) +
    labs(title = paste0("GWAS: ", hit$gwas_name),
         y = expression(-log[10](p)), x = NULL) +
    th + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  p2 <- ggplot(gxe_sub, aes(position / 1e6, neglog10p)) +
    geom_point(size = 0.8, alpha = 0.6, color = "steelblue") +
    geom_point(data = gxe_sub[abs(position - lead_pos) < 1],
               size = 2.5, color = "red", shape = 18) +
    geom_hline(yintercept = -log10(1e-5), linetype = "dashed",
               color = "blue", linewidth = 0.3) +
    coord_cartesian(ylim = c(0, shared_ymax)) +
    labs(title = paste0("GxE Interaction: ", cell_type, " ", modality,
                        " - ", gene_r, "  [pval = b_int]"),
         y = expression(-log[10](p[interaction])), x = NULL) +
    th
  
  panels  <- list(p1, p2)
  heights <- c(1, 1.15)
  
  if (has_me) {
    p2 <- p2 + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    panels[[2]] <- p2
    
    p3 <- ggplot(me_sub, aes(position / 1e6, neglog10p)) +
      geom_point(size = 0.8, alpha = 0.6, color = "darkorange") +
      geom_point(data = me_sub[abs(position - lead_pos) < 1],
                 size = 2.5, color = "red", shape = 18) +
      geom_hline(yintercept = -log10(1e-5), linetype = "dashed",
                 color = "blue", linewidth = 0.3) +
      coord_cartesian(ylim = c(0, shared_ymax)) +
      labs(title = paste0("eQTL Panel: ", cell_type, " ", modality, " - ", gene_r,
                          ifelse(!is.na(preferred_term) && nzchar(preferred_term),
                                 paste0(" [", preferred_term, "]"), "")),
           y = expression(-log[10](p[eQTL])),
           x = paste0("Position on ", chr_r, " (Mb)")) +
      th
    
    panels  <- list(p1, p2, p3)
    heights <- c(1, 1, 1.15)
  } else {
    p2 <- p2 + labs(x = paste0("Position on ", chr_r, " (Mb)"))
    panels[[2]] <- p2
  }
  
  title_str <- paste0(cell_type, " ", modality, " x ", hit$gwas_name,
                      "  |  ", gene_r,
                      "  |  PP.H4=", round(hit$PP.H4, 3),
                      "  |  p_int=", sprintf("%.1e", hit$Lead_P_int))
  
  stacked <- cowplot::plot_grid(plotlist = panels, ncol = 1, align = "v",
                                rel_heights = heights)
  final   <- cowplot::ggdraw() +
    cowplot::draw_label(title_str, x = 0.5, y = 0.99, vjust = 1,
                        size = 10, fontface = "bold") +
    cowplot::draw_plot(stacked, y = 0, height = 0.97)
  
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  fname <- paste0("coloc_interaction_", cell_type, "_", modality, "_",
                  hit$gwas_name, "_", gene_r, ".pdf")
  fpath <- file.path(plot_dir, fname)
  ggsave(fpath, final,
         width  = 8,
         height = ifelse(length(panels) == 3, 9, 6.5),
         dpi    = 300)
  cat("  Saved:", fpath, "\n")
  cat("└──────────────────────────────────────────────────────────┘\n")
  return(invisible(final))
}

###############################################################################
# MAIN EXECUTION
###############################################################################

cat("\n========================================\n")
cat("  Fast mode: seed from all GxE hits TSV\n")
cat("  Criteria: PP.H4 >", SINGLE_MIN_PP_H4,
    "| p_int <", SINGLE_MIN_P_INT,
    "| p_GWAS <", SINGLE_MIN_P_GWAS, "\n")
cat("========================================\n")

if (!file.exists(all_hits_tsv)) stop("all_hits_tsv not found: ", all_hits_tsv)

cat("Loading all GxE hits table...\n")
all_hits <- fread(all_hits_tsv, sep = "\t", header = TRUE)
cat("  Total rows:", nrow(all_hits), "\n")

needed_hits_cols <- c("phenotype_id", "variant_id", "pval", "filename", "effect_type")
missing_hits_cols <- setdiff(needed_hits_cols, names(all_hits))
if (length(missing_hits_cols) > 0)
  stop("Missing required columns in all_hits_tsv: ",
       paste(missing_hits_cols, collapse = ", "))

SINGLE_MIN_P_INT<-1e-5

all_hits <- all_hits[effect_type == "interaction_effect" & pval <= SINGLE_MIN_P_INT]

dim(all_hits)
#SINGLE_MIN_P_INT<-1e-8

all_hits <- all_hits[effect_type == "interaction_effect" & pval <= SINGLE_MIN_P_INT]
dim(all_hits)


all_hits <- unique(all_hits, by = c("phenotype_id", "variant_id", "filename"))
# FIX 2: was all_hits %>% dplyr::filter(Omics=="Metabolome") — dplyr not loaded
#all_hits <- all_hits[Omics == "Metabolome"]

if (!is.na(term_arg)) {
  if (!("term" %in% names(all_hits)))
    stop("--term was provided but column 'term' is missing in all_hits_tsv")
  all_hits <- all_hits[term == term_arg]
  cat("  Term filter active:", term_arg,
      "| rows after term filter:", nrow(all_hits), "\n")
}



gwas_files <- list.files(gwas_dir, pattern = "\\.rds$", full.names = TRUE)
cat("GWAS files:", length(gwas_files), "\n")

best_hit      <- NULL
best_gwas     <- NULL
best_gxe_path <- NA_character_
count_skipped <- 0

double_break  <- FALSE
running_count <- 0
failed_rows   <- list()
count_failed_rows <- 0
analyses_outputs  <- list()
raw_coloc_outputs <- list()

for (gi in seq_along(gwas_files)) {
  if (double_break) break
  gf      <- gwas_files[gi]
  gwas_bn <- str_replace(basename(gf), "\\.rds$", "")
  cat("\n══ GWAS [", gi, "/", length(gwas_files), "]:", gwas_bn, " ══\n")
  
  gwas       <- as.data.table(readRDS(gf))
  setkey(gwas, panel_variant_id)
  gwas_is_cc <- gwas_bn %in% names(CC_GWAS)
  
  gwas_near <- gwas[pvalue < SINGLE_MIN_P_GWAS & !is.na(chromosome) & !is.na(position),
                    .(chromosome, gwas_pos = position, gwas_pvalue = pvalue)]
  cat("  GWAS rows:", nrow(gwas),
      "| GWAS variants with p<", SINGLE_MIN_P_GWAS, ":", nrow(gwas_near), "\n", sep = "")
  if (nrow(gwas_near) == 0) {
    rm(gwas); gc(verbose = FALSE)
    cat("  Skip: no GWAS variants pass near-locus p-value gate.\n")
    next
  }
  
  gene_hits <- copy(all_hits)
  parsed    <- parse_variant_chr_pos(gene_hits$variant_id)
  gene_hits[, chromosome := parsed$chromosome]
  gene_hits[, position   := parsed$position]
  gene_hits <- gene_hits[!is.na(chromosome) & !is.na(position)]
  before_mhc <- nrow(gene_hits)
  gene_hits  <- gene_hits[!is_in_mhc(chromosome, position)]
  cat("  Gene-level interaction seeds:", nrow(gene_hits),
      "(MHC excluded:", before_mhc - nrow(gene_hits), ")\n")
  if (nrow(gene_hits) == 0) {
    rm(gwas, gwas_near, gene_hits); gc(verbose = FALSE)
    cat("  Skip: no non-MHC interaction seeds.\n")
    next
  }
  
  gwas_near[, `:=`(near_start = pmax(1, gwas_pos - GWAS_NEAR_WINDOW),
                   near_end   = gwas_pos + GWAS_NEAR_WINDOW)]
  near_pairs   <- gene_hits[gwas_near,
                            on = .(chromosome,
                                   position >= near_start,
                                   position <= near_end),
                            nomatch = 0L, allow.cartesian = TRUE]
  gwas_p_col   <- if ("i.gwas_pvalue" %in% names(near_pairs)) "i.gwas_pvalue" else "gwas_pvalue"
  near_summary <- near_pairs[, .(gwas_pvalue = min(get(gwas_p_col), na.rm = TRUE)),
                             by = .(phenotype_id, variant_id)]
  overlap      <- merge(gene_hits, near_summary,
                        by = c("phenotype_id", "variant_id"), all = FALSE)
  overlap      <- unique(overlap)
  cat("  Seeds with >=1 GWAS variant (p<", SINGLE_MIN_P_GWAS,
      ") within ±", GWAS_NEAR_WINDOW, "bp:", nrow(overlap), "\n", sep = "")
  if (nrow(overlap) == 0) {
    rm(gwas, gwas_near, gene_hits, overlap, near_pairs, near_summary); gc(verbose = FALSE)
    cat("  Skip: 0 qualifying seeds near GWAS signal.\n")
    next
  }
  
  dup_cols <- which(duplicated(names(overlap)))
  if (length(dup_cols) > 0) overlap <- overlap[, !dup_cols, with = FALSE]
  names(overlap) <- make.unique(names(overlap))
  
  overlap[, region_start := pmax(1, position - WINDOW_SIZE)]
  overlap[, region_end   := position + WINDOW_SIZE]
  setorder(overlap, phenotype_id, chromosome, position)
  
  overlap[, window_id := {
    max_end_so_far <- cummax(data.table::shift(region_end, n = 1, fill = -1))
    new_window     <- region_start > max_end_so_far
    as.integer(cumsum(new_window))
  }, by = .(phenotype_id, cell_type, term, chromosome)]
  
  overlap[, global_window_id := .GRP,
          by = .(phenotype_id, cell_type, term, chromosome, window_id)]
  
  setorder(overlap, window_id, phenotype_id, pval, gwas_pvalue)
  unique_windows <- unique(overlap[!is.na(global_window_id),
                                   .(global_window_id, phenotype_id)])
  cat("  Trying", nrow(unique_windows),
      "non-overlapping region-gene windows for this GWAS...\n")
  
  for (ri in seq_len(nrow(unique_windows))) {
    running_count <- running_count + 1
    cat("  [window", ri, "/", nrow(unique_windows),
        "| total:", running_count, "]\n")
    
    win       <- unique_windows[ri]
    win_seeds <- overlap[global_window_id == win$global_window_id]
    seed      <- win_seeds[1]
    
    cat(sprintf("    Region-Gene-Window [%d/%d]: %s:%d-%d | %s | %s | p_int=%s | p_gwas=%s\n",
                ri, nrow(unique_windows),
                seed$chromosome, seed$region_start, seed$region_end,
                seed$phenotype_id, seed$variant_id,
                sprintf("%.1e", seed$pval),
                sprintf("%.1e", seed$gwas_pvalue)))
    
    gxe_hint_path <- resolve_gxe_rds_path(seed)
    
    res_wrapper <- tryCatch(
      run_coloc_seed(gwas, gwas_bn, seed, gxe_hint_path, gwas_is_cc = gwas_is_cc),
      error = function(e) {
        message("Failed for seed: ", seed$phenotype_id, " -> ", e$message)
        count_skipped    <<- count_skipped + 1
        seed$reason_for_fail <<- e$message
        failed_rows[[length(failed_rows) + 1]] <<- seed
        return(NULL)
      }
    )
    
    if (is.null(res_wrapper)) next
    
    res            <- res_wrapper$summary
    raw_coloc_data <- res_wrapper$raw_coloc
    if (is.null(res)) next
    
    list_name <- paste(res$Gene, res$gwas_name, res$Seed_term, res$Lead_SNP, sep = "_")
    raw_coloc_outputs[[list_name]]              <- raw_coloc_data
    analyses_outputs[[length(analyses_outputs) + 1]] <- res
    
    if (res$PP.H4 > SINGLE_MIN_PP_H4) {
      cat("      -> Accepted hit (PP.H4 above threshold).\n")
      best_hit      <- res
      best_gwas     <- gwas
      best_gxe_path <- gxe_hint_path
    } else {
      cat(sprintf("      -> Rejected: PP.H4=%.4f < %.4f\n",
                  res$PP.H4, SINGLE_MIN_PP_H4))
    }
  }
  
  rm(gwas, gwas_near, overlap, near_pairs, near_summary); gc(verbose = FALSE)
}

if (is.null(best_hit)) {
  cat("\nNo seed overlap produced a valid coloc result.\n")
} else {
  cat("\n*** PLOTTING TOP FAST HIT ***\n")
  cat("  ", best_hit$gxe_name, " x ", best_hit$gwas_name, "\n")
  cat("  Gene:", best_hit$Gene, "| PP.H4:", round(best_hit$PP.H4, 4),
      "| p_int:", sprintf("%.1e", best_hit$Lead_P_int),
      "| p_GWAS:", sprintf("%.1e", best_hit$Lead_P_GWAS), "\n")
  plot_coloc_hit(best_hit, best_gwas, best_gxe_path)
}

cat("count_skipped:", count_skipped, "\n")

if (length(analyses_outputs) > 0) {
  all_results_dt  <- data.table::rbindlist(analyses_outputs, fill = TRUE)
  term_tag        <- if (is.na(term_arg)) "all_terms" else gsub("[^A-Za-z0-9._-]", "_", term_arg)
  all_results_rds <- file.path(output_dir, paste0("all_coloc_seed_outputs_V2", term_tag, ".rds"))
  all_results_tsv <- file.path(output_dir, paste0("all_coloc_seed_outputs_V2", term_tag, ".tsv"))
  raw_data_path   <- file.path(output_dir, paste0("all_coloc_seed_outputs_V2", term_tag, "_raw_data_colocs_.rds"))
  saveRDS(all_results_dt,    all_results_rds)
  saveRDS(raw_coloc_outputs, raw_data_path)
  data.table::fwrite(all_results_dt, all_results_tsv, sep = "\t")
  cat("Saved all successful seed outputs:\n  ", all_results_rds,
      "\n  ", all_results_tsv, "\n", sep = "")
} else {
  cat("No successful seed outputs to save.\n")
}

cat("\n=== SCRIPT FINISHED ===\n")