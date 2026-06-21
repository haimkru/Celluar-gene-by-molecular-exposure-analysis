#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

msg <- function(...) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("[%s] %s", ts, paste0(..., collapse = "")))
}

# ---------------------------
# Hardcoded configuration
# ---------------------------
fairfax_file <- "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/other_responseQTLs/fairfex_all_response_qtls.csv"
gxe_sig_hits_file <- "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/all_cell_types_sig_hits.tsv"
out_prefix <- "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/fisher_fairfex/fairfax_gxe_gene"

# Preferred: a local 2-column mapping file with Ensembl ID and gene symbol.
# If missing, script will try org.Hs.eg.db fallback.
ensembl_symbol_map_file <- "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/fisher_fairfex/ensembl_to_symbol.tsv"

fairfax_fdr_threshold <- 0.01
fairfax_naive_fdr_threshold <- 0.01
fisher_alternative <- "greater"

# Universe is derived from genes actually tested in Fairfax (computed below from fairfax_all_ensembl).
# Do not set manually.

gxe_pval_threshold <- 1e-5
gxe_pval_col <- "pval"

# Optional: file with one column of tested Ensembl gene IDs (header allowed).
# If provided, both Fairfax and GxE sets are restricted to this universe.
tested_gene_universe_file <- ""

# Monocyte cell type grep pattern — check msg output and adjust if needed.
monocyte_pattern <- "mono"

# Fairfax columns expected in the file header.
fairfax_gene_col <- "Gene"
fairfax_naive_fdr_col <- "Naive.FDR"
fairfax_non_naive_fdr_cols <- c("LPS2.FDR", "LPS24.FDR", "IFN.FDR")

# ---------------------------
# Validation
# ---------------------------
if (is.na(fairfax_fdr_threshold) || fairfax_fdr_threshold <= 0 || fairfax_fdr_threshold > 1) {
  stop("fairfax_fdr_threshold must be in (0,1].")
}
if (is.na(fairfax_naive_fdr_threshold) || fairfax_naive_fdr_threshold <= 0 || fairfax_naive_fdr_threshold > 1) {
  stop("fairfax_naive_fdr_threshold must be in (0,1].")
}
if (!fisher_alternative %in% c("greater", "two.sided", "less")) {
  stop("fisher_alternative must be one of greater/two.sided/less")
}
if (is.na(gxe_pval_threshold) || gxe_pval_threshold <= 0 || gxe_pval_threshold > 1) {
  stop("gxe_pval_threshold must be in (0,1].")
}
if (!file.exists(fairfax_file)) stop("Missing file: ", fairfax_file)
if (!file.exists(gxe_sig_hits_file)) stop("Missing file: ", gxe_sig_hits_file)

out_dir <- dirname(out_prefix)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

msg("Starting gene-level pipeline")
msg("fairfax_file: ", fairfax_file)
msg("gxe_sig_hits_file: ", gxe_sig_hits_file)
msg("ensembl_symbol_map_file: ", ensembl_symbol_map_file)
msg("out_prefix: ", out_prefix)
msg("fairfax_fdr_threshold (treatment): ", fairfax_fdr_threshold)
msg("fairfax_naive_fdr_threshold (naive): ", fairfax_naive_fdr_threshold)
msg("gxe_pval_threshold: ", gxe_pval_threshold)
msg("gxe_pval_col: ", gxe_pval_col)
if (nzchar(tested_gene_universe_file)) {
  msg("tested_gene_universe_file: ", tested_gene_universe_file)
}

pipeline_start <- proc.time()[["elapsed"]]

# ---------------------------
# Fairfax gene hit definition
# ---------------------------
stage_start <- proc.time()[["elapsed"]]

fairfax_all_fdr_cols <- c(fairfax_naive_fdr_col, fairfax_non_naive_fdr_cols)
fairfax <- fread(
  fairfax_file,
  header = TRUE,
  select = c(fairfax_gene_col, fairfax_all_fdr_cols),
  colClasses = c(
    setNames("character", fairfax_gene_col),
    setNames(rep("numeric", length(fairfax_all_fdr_cols)), fairfax_all_fdr_cols)
  ),
  showProgress = TRUE
)

required_fairfax_cols <- c(fairfax_gene_col, fairfax_all_fdr_cols)
missing_fairfax_cols <- setdiff(required_fairfax_cols, names(fairfax))
if (length(missing_fairfax_cols) > 0) {
  stop("Missing required Fairfax columns: ", paste(missing_fairfax_cols, collapse = ", "))
}

# Treat missing FDR values as non-significant.
for (col_nm in fairfax_all_fdr_cols) {
  fairfax[is.na(get(col_nm)), (col_nm) := 1]
}

fairfax[, gene_symbol := toupper(trimws(get(fairfax_gene_col)))]
fairfax[gene_symbol %chin% c("", "NA", "NAN"), gene_symbol := NA_character_]

fairfax_all_genes <- unique(fairfax[!is.na(gene_symbol), gene_symbol])

# Aggregate to per-gene minimum FDR across all rows (SNPs) before thresholding.
# This prevents inflating hits by counting genes with many SNPs tested.
fairfax_gene_level <- fairfax[!is.na(gene_symbol), lapply(.SD, min, na.rm = TRUE),
                              by = gene_symbol, .SDcols = fairfax_all_fdr_cols]
# Replace Inf (from all-NA columns after min) with 1.
for (col_nm in fairfax_all_fdr_cols) {
  fairfax_gene_level[is.infinite(get(col_nm)), (col_nm) := 1]
}

msg("Fairfax rows: ", format(nrow(fairfax), big.mark = ","),
    " → unique genes after per-gene min-FDR aggregation: ",
    format(nrow(fairfax_gene_level), big.mark = ","))

fairfax_non_naive_sig_genes <- unique(
  fairfax_gene_level[
    get(fairfax_non_naive_fdr_cols[1]) <= fairfax_fdr_threshold |
      get(fairfax_non_naive_fdr_cols[2]) <= fairfax_fdr_threshold |
      get(fairfax_non_naive_fdr_cols[3]) <= fairfax_fdr_threshold,
    gene_symbol
  ]
)

# Naive-significant genes: also significant in unstimulated monocytes.
# Uses a more lenient threshold to aggressively remove constitutive eQTLs.
fairfax_naive_sig_genes <- unique(
  fairfax_gene_level[
    get(fairfax_naive_fdr_col) <= fairfax_naive_fdr_threshold,
    gene_symbol
  ]
)

# Interaction-hit genes: significant in stimulated BUT NOT naive.
# These are response/interaction eQTLs — condition-specific effects.
fairfax_hit_genes <- setdiff(fairfax_non_naive_sig_genes, fairfax_naive_sig_genes)

msg("Fairfax all genes: ", format(length(fairfax_all_genes), big.mark = ","))
msg("Fairfax non-naive significant genes: ", format(length(fairfax_non_naive_sig_genes), big.mark = ","))
msg("Fairfax naive-only significant genes: ", format(length(fairfax_naive_sig_genes), big.mark = ","))
msg("Fairfax interaction-hit genes (stimulated minus naive): ", format(length(fairfax_hit_genes), big.mark = ","))
msg("Stage time Fairfax genes: ", round(proc.time()[["elapsed"]] - stage_start, 2), " sec")

# ---------------------------
# GxE expression interaction genes (Ensembl)
# ---------------------------
stage_start <- proc.time()[["elapsed"]]

current_gxe_meta <- list(
  source_file = gxe_sig_hits_file,
  pval_col = gxe_pval_col,
  pval_threshold = gxe_pval_threshold
)

reuse_gxe_in_memory <- FALSE
if (exists("gxe_hits", inherits = FALSE)) {
  required_cols <- c("phenotype_id", "Omics", "assay", "cell_type", "term", "ensembl_gene_id")
  has_required_cols <- all(required_cols %in% names(gxe_hits))
  cached_meta <- attr(gxe_hits, "gxe_meta")
  meta_matches <- is.list(cached_meta) &&
    identical(cached_meta$source_file, current_gxe_meta$source_file) &&
    identical(cached_meta$pval_col, current_gxe_meta$pval_col) &&
    isTRUE(all.equal(cached_meta$pval_threshold, current_gxe_meta$pval_threshold))
  
  reuse_gxe_in_memory <- has_required_cols && meta_matches
  if (reuse_gxe_in_memory) {
    msg("Using existing in-memory object: gxe_hits (metadata matched)")
  } else {
    msg("Existing in-memory gxe_hits is stale/incompatible; rebuilding from file")
  }
}

if (!reuse_gxe_in_memory) {
  gxe_hits <- fread(
    gxe_sig_hits_file,
    select = c("phenotype_id", "term", gxe_pval_col, "effect_type", "Omics", "assay", "cell_type"),
    colClasses = c(
      "phenotype_id" = "character",
      "term" = "character",
      setNames("numeric", gxe_pval_col),
      "effect_type" = "character",
      "Omics" = "character",
      "assay" = "character",
      "cell_type" = "character"
    ),
    sep = "\t",
    showProgress = TRUE
  )
  
  gxe_hits <- gxe_hits[
    effect_type == "interaction_effect" &
      tolower(trimws(assay)) == "expression" &
      !is.na(get(gxe_pval_col)) &
      get(gxe_pval_col) <= gxe_pval_threshold,
    .(phenotype_id, Omics, assay, cell_type, term)
  ]
  
  if (nrow(gxe_hits) == 0) {
    stop("No GxE rows remain after filtering to interaction_effect + Expression + pval threshold")
  }
  
  gxe_hits[, ensembl_gene_id := sub("\\..*$", "", trimws(phenotype_id))]
  gxe_hits <- gxe_hits[!is.na(ensembl_gene_id) & ensembl_gene_id != ""]
  gxe_hits <- unique(gxe_hits, by = c("Omics", "assay", "cell_type", "term", "ensembl_gene_id"))
  attr(gxe_hits, "gxe_meta") <- current_gxe_meta
  msg("Loaded and filtered GxE from file")
}

if (nrow(gxe_hits) == 0) {
  stop("No GxE rows available after cache/load/filtering")
}

# Print cell types so user can verify monocyte_pattern matches correctly.
msg("GxE cell types present: ", paste(sort(unique(gxe_hits$cell_type)), collapse = ", "))
msg("GxE expression interaction rows (dedup by stratum+gene): ", format(nrow(gxe_hits), big.mark = ","))

# Never mutate the reusable in-memory cache object in downstream steps.
gxe_hits_work <- copy(gxe_hits)

# ---------------------------
# Ensembl -> gene symbol mapping
# ---------------------------
stage_map <- proc.time()[["elapsed"]]

map_dt <- NULL
if (file.exists(ensembl_symbol_map_file)) {
  msg("Using mapping file: ", ensembl_symbol_map_file)
  map_raw <- fread(ensembl_symbol_map_file, showProgress = TRUE)
  map_names_lc <- tolower(names(map_raw))
  
  ens_idx <- which(map_names_lc %in% c("ensembl_gene_id", "ensembl", "gene_id", "phenotype_id"))
  sym_idx <- which(map_names_lc %in% c("gene_symbol", "symbol", "gene_name", "hgnc_symbol"))
  
  if (length(ens_idx) == 0 || length(sym_idx) == 0) {
    stop("Mapping file must contain Ensembl and gene symbol columns")
  }
  
  ens_col <- names(map_raw)[ens_idx[1]]
  sym_col <- names(map_raw)[sym_idx[1]]
  
  map_dt <- map_raw[, .(
    ensembl_gene_id = sub("\\..*$", "", trimws(as.character(get(ens_col)))),
    gene_symbol = toupper(trimws(as.character(get(sym_col))))
  )]
} else {
  msg("Mapping file not found. Trying org.Hs.eg.db fallback with Fairfax SYMBOL->ENSEMBL mapping")
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) || !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("No mapping file and org.Hs.eg.db not available. Provide ensembl_symbol_map_file.")
  }
  
  fairfax_syms <- unique(fairfax_all_genes)
  ens_ids <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = fairfax_syms,
    keytype = "SYMBOL",
    column = "ENSEMBL",
    multiVals = "first"
  )
  map_dt <- data.table(
    gene_symbol = toupper(trimws(as.character(names(ens_ids)))),
    ensembl_gene_id = sub("\\..*$", "", trimws(as.character(unname(ens_ids))))
  )
}

map_dt <- map_dt[
  !is.na(ensembl_gene_id) & ensembl_gene_id != "" &
    !is.na(gene_symbol) & gene_symbol != "" & gene_symbol != "NA",
  .(ensembl_gene_id, gene_symbol)
]
map_dt <- unique(map_dt)

msg("Mapped Ensembl<->symbol pairs (unique pairs): ", format(nrow(map_dt), big.mark = ","))

# Build symbol->Ensembl dictionary for Fairfax-to-Ensembl harmonization.
symbol_to_ens <- unique(map_dt[, .(gene_symbol, ensembl_gene_id)])
fairfax_map <- symbol_to_ens[data.table(gene_symbol = fairfax_all_genes), on = "gene_symbol"]
fairfax_all_ensembl <- unique(fairfax_map[!is.na(ensembl_gene_id), ensembl_gene_id])

fairfax_non_naive_map <- symbol_to_ens[data.table(gene_symbol = fairfax_non_naive_sig_genes), on = "gene_symbol"]
fairfax_non_naive_ensembl <- unique(fairfax_non_naive_map[!is.na(ensembl_gene_id), ensembl_gene_id])

fairfax_hit_map <- symbol_to_ens[data.table(gene_symbol = fairfax_hit_genes), on = "gene_symbol"]
fairfax_hit_ensembl <- unique(fairfax_hit_map[!is.na(ensembl_gene_id), ensembl_gene_id])

fairfax_unmapped_symbols_n <- uniqueN(data.table(gene_symbol = fairfax_all_genes)[!symbol_to_ens, on = "gene_symbol", gene_symbol])
msg("Fairfax symbols mapped to Ensembl: ", format(length(fairfax_all_ensembl), big.mark = ","))
msg("Fairfax symbols not mapped to Ensembl: ", format(fairfax_unmapped_symbols_n, big.mark = ","))

# ---------------------------
# Derive analysis universe from tested genes (not from significant-hit overlap).
# Default: Fairfax tested Ensembl genes.
# If tested_gene_universe_file is provided: intersect Fairfax tested genes with that file.
# ---------------------------
gxe_all_genes <- unique(gxe_hits_work$ensembl_gene_id)
msg("GxE expression interaction Ensembl genes: ", format(length(gxe_all_genes), big.mark = ","))

analysis_universe_ids <- unique(fairfax_all_ensembl)
if (nzchar(tested_gene_universe_file)) {
  if (!file.exists(tested_gene_universe_file)) {
    stop("Missing tested_gene_universe_file: ", tested_gene_universe_file)
  }
  tested_u <- fread(tested_gene_universe_file, showProgress = TRUE)
  if (ncol(tested_u) < 1) stop("tested_gene_universe_file has no columns")
  tested_ids <- unique(sub("\\..*$", "", trimws(as.character(tested_u[[1]]))))
  tested_ids <- tested_ids[!is.na(tested_ids) & tested_ids != ""]
  if (length(tested_ids) == 0) stop("tested_gene_universe_file has no usable IDs")
  analysis_universe_ids <- intersect(analysis_universe_ids, tested_ids)
  msg("Universe source: Fairfax tested genes ∩ tested_gene_universe_file")
} else {
  msg("Universe source: Fairfax tested genes (no tested_gene_universe_file provided)")
}

analysis_universe_n <- length(analysis_universe_ids)
msg("Analysis universe size (Ensembl): ", format(analysis_universe_n, big.mark = ","))
if (analysis_universe_n == 0) {
  stop("Universe is empty after tested-gene restriction. Check identifier harmonization.")
}
msg("Stage time mapping: ", round(proc.time()[["elapsed"]] - stage_map, 2), " sec")

# ---------------------------
# Gene-level Fisher with Fairfax-derived universe
# ---------------------------
stage_start <- proc.time()[["elapsed"]]

fairfax_set_all <- unique(fairfax_hit_ensembl)
fairfax_set <- intersect(fairfax_set_all, analysis_universe_ids)
msg("Fairfax hit Ensembl genes (all): ", format(length(fairfax_set_all), big.mark = ","))
msg("Fairfax hit Ensembl genes (in universe): ", format(length(fairfax_set), big.mark = ","))

if (nrow(gxe_hits_work) == 0) {
  stop("No GxE expression genes remain after Ensembl harmonization")
}

# Restrict GxE genes to the analysis universe so Fisher denominator is valid.
gxe_unique_pre <- length(unique(gxe_hits_work$ensembl_gene_id))
gxe_hits_work <- gxe_hits_work[ensembl_gene_id %chin% analysis_universe_ids]
gxe_unique_post <- length(unique(gxe_hits_work$ensembl_gene_id))
gxe_dropped_n <- gxe_unique_pre - gxe_unique_post
if (gxe_dropped_n > 0) {
  msg(
    "Dropped GxE genes not in analysis universe: ",
    format(gxe_dropped_n, big.mark = ","),
    " (kept ", format(gxe_unique_post, big.mark = ","), ")"
  )
}
if (nrow(gxe_hits_work) == 0) {
  stop("No GxE genes remain after restricting to analysis universe")
}

setkey(gxe_hits_work, Omics, assay, cell_type, term, ensembl_gene_id)
gxe_hits_work <- unique(gxe_hits_work)

gxe_global_set <- unique(gxe_hits_work$ensembl_gene_id)

# Diagnostic overlaps.
ov_all       <- length(intersect(gxe_global_set, fairfax_all_ensembl))
ov_non_naive <- length(intersect(gxe_global_set, fairfax_non_naive_ensembl))
ov_interaction <- length(intersect(gxe_global_set, fairfax_set))
msg("Diagnostic overlap counts (GxE vs Fairfax): all=", ov_all, ", non_naive=", ov_non_naive, ", interaction=", ov_interaction)

union_n <- length(unique(c(gxe_global_set, fairfax_set)))
msg("Global set sizes: |Fairfax|=", length(fairfax_set), ", |GxE|=", length(gxe_global_set), ", |Union|=", union_n, ", N=", analysis_universe_n)
if (union_n > analysis_universe_n) {
  stop(
    "Inconsistent universe after universe restriction: union(Fairfax,GxE)=", union_n,
    " exceeds N=", analysis_universe_n,
    ". Check identifier harmonization in Ensembl IDs and tested-gene universe inputs."
  )
}

# ---------------------------
# Per-term Fisher (all cell types)
# ---------------------------
gxe_counts <- gxe_hits_work[, .(
  gxe_hits_n = .N,
  a_overlap  = sum(ensembl_gene_id %chin% fairfax_set)
), by = .(Omics, assay, cell_type, term)]

res <- vector("list", nrow(gxe_counts))
for (j in seq_len(nrow(gxe_counts))) {
  om    <- gxe_counts$Omics[j]
  asy   <- gxe_counts$assay[j]
  ct    <- gxe_counts$cell_type[j]
  trm   <- gxe_counts$term[j]
  a     <- gxe_counts$a_overlap[j]
  gxe_n <- gxe_counts$gxe_hits_n[j]
  
  b <- gxe_n - a
  c <- length(fairfax_set) - a
  d <- analysis_universe_n - (a + b + c)
  
  if (d < 0) stop("Negative d for term ", trm, ". Check universe restriction and identifier harmonization.")
  
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE), alternative = fisher_alternative)
  
  res[[j]] <- data.table(
    Omics = om, assay = asy, cell_type = ct, term = trm,
    a_overlap = a, b_gxe_only = b, c_fairfax_only = c, d_neither = d,
    gxe_hits_n = gxe_n,
    fairfax_hits_n = length(fairfax_set),
    shared_universe_n = analysis_universe_n,
    odds_ratio = unname(ft$estimate),
    p_value = ft$p.value,
    conf_low = ft$conf.int[1],
    conf_high = ft$conf.int[2]
  )
  
  if (j %% 25 == 0 || j == nrow(gxe_counts)) {
    msg("Processed gene-level terms: ", j, "/", nrow(gxe_counts))
  }
}

res <- rbindlist(res, fill = TRUE)
res[, p_adj_bh := p.adjust(p_value, method = "BH")]
setorder(res, p_value)

# ---------------------------
# Global Fisher helper
# ---------------------------
run_global_fisher <- function(gxe_gene_set, fairfax_set, universe_n, label) {
  a <- length(intersect(gxe_gene_set, fairfax_set))
  b <- length(gxe_gene_set) - a
  c <- length(fairfax_set) - a
  d <- universe_n - (a + b + c)
  if (d < 0) stop("Negative d in global overlap test (", label, "). Check universe.")
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE), alternative = fisher_alternative)
  data.table(
    label = label,
    a_overlap = a, b_gxe_only = b, c_fairfax_only = c, d_neither = d,
    gxe_hits_n = length(gxe_gene_set),
    fairfax_hits_n = length(fairfax_set),
    shared_universe_n = universe_n,
    odds_ratio = unname(ft$estimate),
    p_value = ft$p.value,
    conf_low = ft$conf.int[1],
    conf_high = ft$conf.int[2],
    alternative = fisher_alternative
  )
}

# Global — all cell types
global_all <- run_global_fisher(gxe_global_set, fairfax_set, analysis_universe_n, "all_cell_types")

# Global — monocyte only
gxe_mono_set <- unique(
  gxe_hits_work[grepl(monocyte_pattern, cell_type, ignore.case = TRUE), ensembl_gene_id]
)
msg("GxE monocyte-matching cell types: ",
    paste(sort(unique(gxe_hits_work[grepl(monocyte_pattern, cell_type, ignore.case = TRUE), cell_type])), collapse = ", "))
msg("GxE monocyte genes: ", format(length(gxe_mono_set), big.mark = ","))

global_mono <- run_global_fisher(gxe_mono_set, fairfax_set, analysis_universe_n, "monocyte_only")

global_overlap <- rbind(global_all, global_mono)
msg("Global Fisher (all cell types): OR=", round(global_all$odds_ratio, 3), "  p=", signif(global_all$p_value, 3))
msg("Global Fisher (monocyte only):  OR=", round(global_mono$odds_ratio, 3), "  p=", signif(global_mono$p_value, 3))
if (global_all$d_neither == 0 || global_mono$d_neither == 0) {
  msg("WARNING: d_neither is 0 for at least one global test; Fisher result may be uninformative due to saturated universe.")
}

gene_overlap_set <- intersect(gxe_global_set, fairfax_set)

msg("Stage time gene-level Fisher: ", round(proc.time()[["elapsed"]] - stage_start, 2), " sec")

# ---------------------------
# Outputs
# ---------------------------
out_table          <- paste0(out_prefix, ".fisher_by_term.tsv")
out_global         <- paste0(out_prefix, ".global_overlap_fisher.tsv")
out_overlap_genes  <- paste0(out_prefix, ".global_overlap_genes.tsv.gz")
out_observed_genes <- paste0(out_prefix, ".observed_gene_ids.tsv.gz")
out_counts         <- paste0(out_prefix, ".counts_summary.tsv")

fwrite(res,                                                              out_table,          sep = "\t")
fwrite(global_overlap,                                                   out_global,         sep = "\t")
fwrite(data.table(gene_id = gene_overlap_set),                           out_overlap_genes,  sep = "\t")
fwrite(data.table(gene_id = unique(c(fairfax_all_ensembl, gxe_all_genes))), out_observed_genes, sep = "\t")

counts_summary <- data.table(
  metric = c(
    "fairfax_all_genes_n",
    "fairfax_non_naive_sig_genes_n",
    "fairfax_interaction_hit_genes_n",
    "fairfax_mapped_to_ensembl_n",
    "fairfax_unmapped_symbols_n",
    "fairfax_non_naive_mapped_ensembl_n",
    "fairfax_interaction_hit_ensembl_n",
    "gxe_expression_interaction_ensembl_n",
    "gxe_monocyte_ensembl_n",
    "diagnostic_overlap_all_n",
    "diagnostic_overlap_non_naive_n",
    "diagnostic_overlap_interaction_n",
    "analysis_universe_n",
    "fairfax_hit_ensembl_n",
    "gxe_global_genes_n",
    "fairfax_gxe_union_genes_n",
    "global_overlap_n_all_celltypes",
    "global_overlap_fisher_p_all_celltypes",
    "global_overlap_fisher_or_all_celltypes",
    "global_overlap_n_mono",
    "global_overlap_fisher_p_mono",
    "global_overlap_fisher_or_mono",
    "fisher_tests_n"
  ),
  value = c(
    length(fairfax_all_genes),
    length(fairfax_non_naive_sig_genes),
    length(fairfax_hit_genes),
    length(fairfax_all_ensembl),
    fairfax_unmapped_symbols_n,
    length(fairfax_non_naive_ensembl),
    length(fairfax_hit_ensembl),
    length(gxe_all_genes),
    length(gxe_mono_set),
    ov_all,
    ov_non_naive,
    ov_interaction,
    analysis_universe_n,
    length(fairfax_set),
    length(gxe_global_set),
    union_n,
    global_all$a_overlap,
    global_all$p_value,
    global_all$odds_ratio,
    global_mono$a_overlap,
    global_mono$p_value,
    global_mono$odds_ratio,
    nrow(res)
  )
)

fwrite(counts_summary, out_counts, sep = "\t")

msg("Wrote: ", out_table)
msg("Wrote: ", out_global)
msg("Wrote: ", out_overlap_genes)
msg("Wrote: ", out_observed_genes)
msg("Wrote: ", out_counts)
msg("Total pipeline time: ", round(proc.time()[["elapsed"]] - pipeline_start, 2), " sec")
msg("Done")
