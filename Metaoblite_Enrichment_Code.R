###Metaoblite_Enrichment_Code####
#/oak/stanford/groups/smontgom/hkrupkin/Code/Metaoblite_Enrichment_Code
#author: haim krupkin
#date: 04/25/2026

###description#####
#this script This file is meant to check what annotations/pathways each module of metabolites is associated with



#!/usr/bin/env Rscript
# =============================================================================
# HMDB Multi-Type Metabolite Enrichment Analysis — ALL MODULES
# =============================================================================

# ── 0. Configuration ──────────────────────────────────────────────────────────
ANNOTATION_FILE <- "/labs/smontgom/grps_smontgom/mahuynh/metabolites/MESA/data/annotation_metabolites_met_ID/annotation_table2.csv"
MODULES_DIR     <- "/oak/stanford/groups/smontgom/mahuynh/metabolites/MESA/data/associated_features/Exam 1 - metabolites/Modules"
RESULTS_DIR     <- "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Code/METABOLITE_PATHWAY_ERNICHEMTN/outputV2"

USE_CSV_FILES <- TRUE
CSV_DIR       <- "/oak/stanford/groups/smontgom/hkrupkin/GxE/Code/METABOLITE_PATHWAY_ERNICHEMTN/hmdb_enrichment_package"
HMDB_XML      <- "hmdb_metabolites.xml/hmdb_metabolites.xml"

MIN_PATH_SIZE <- 3
MAX_PATH_SIZE <- 500
NOMINAL_ALPHA <- 0.05
FDR_ALPHA     <- 0.05

# ── 1. Load packages ───────────────────────────────────────────────────────────
message("=== Step 1: Loading packages ===")
library(XML)

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 2. Load HMDB annotations ──────────────────────────────────────────────────
message("=== Step 2: Loading HMDB annotations ===")

if (USE_CSV_FILES) {
  message("  Loading pre-built CSV annotation files...")
  required_csvs <- c("hmdb_pathways.csv","hmdb_diseases.csv","hmdb_biospecimen.csv",
                     "hmdb_tissue.csv","hmdb_taxonomy.csv","hmdb_names.csv")
  missing <- required_csvs[!file.exists(file.path(CSV_DIR, required_csvs))]
  if (length(missing) > 0)
    stop("Missing CSV files in CSV_DIR='", CSV_DIR, "': ", paste(missing, collapse=", "))
  
  df_pw  <- read.csv(file.path(CSV_DIR, "hmdb_pathways.csv"),    stringsAsFactors=FALSE)
  df_dis <- read.csv(file.path(CSV_DIR, "hmdb_diseases.csv"),    stringsAsFactors=FALSE)
  df_bio <- read.csv(file.path(CSV_DIR, "hmdb_biospecimen.csv"), stringsAsFactors=FALSE)
  df_tis <- read.csv(file.path(CSV_DIR, "hmdb_tissue.csv"),      stringsAsFactors=FALSE)
  df_tax <- read.csv(file.path(CSV_DIR, "hmdb_taxonomy.csv"),    stringsAsFactors=FALSE)
  df_nm  <- read.csv(file.path(CSV_DIR, "hmdb_names.csv"),       stringsAsFactors=FALSE)
  
  to_list <- function(df, key_col, val_col) {
    df <- df[nzchar(trimws(df[[val_col]])), ]
    tapply(df[[val_col]], df[[key_col]], function(x) unique(trimws(x)), simplify=FALSE)
  }
  
  pathways_map    <- to_list(df_pw,  "accession", "pathway")
  diseases_map    <- to_list(df_dis, "accession", "disease")
  biospecimen_map <- to_list(df_bio, "accession", "biospecimen")
  tissue_map      <- to_list(df_tis, "accession", "tissue")
  names_map       <- setNames(df_nm$name, df_nm$accession)
  
  taxonomy_map <- setNames(
    lapply(seq_len(nrow(df_tax)), function(i) {
      r   <- df_tax[i, ]
      out <- list()
      for (col in c("kingdom","super_class","class","sub_class","direct_parent")) {
        if (nzchar(trimws(r[[col]]))) out[[col]] <- trimws(r[[col]])
      }
      out
    }),
    df_tax$accession
  )
  taxonomy_map <- taxonomy_map[sapply(taxonomy_map, length) > 0]
  
  message(sprintf("  Loaded: %d pathway, %d disease, %d biospecimen, %d tissue, %d taxonomy entries",
                  length(pathways_map), length(diseases_map),
                  length(biospecimen_map), length(tissue_map), length(taxonomy_map)))
  
} else {
  message("  Parsing ", HMDB_XML, " via SAX (this takes ~6 minutes)...")
  if (!file.exists(HMDB_XML))
    stop("HMDB XML not found: ", HMDB_XML)
  
  pathways_map <- diseases_map <- biospecimen_map <- tissue_map <-
    taxonomy_map <- list()
  names_map <- character(0)
  
  env <- new.env(parent = emptyenv())
  env$acc <- NULL; env$cname <- NULL
  env$acc_done <- env$name_done <- env$in_met <- FALSE
  env$in_pathways <- env$in_pathway <- env$in_diseases <- env$in_disease <- FALSE
  env$in_biospec <- env$in_tissue <- FALSE
  env$pw_name <- env$dis_name <- env$bio_name <- env$tis_name <- NULL
  env$tax <- list(); env$cur_tag <- NULL; env$n <- 0L; env$t0 <- proc.time()["elapsed"]
  
  strip_ns    <- function(name) sub("^.*:", "", name)
  collapse_pw <- function(name) {
    if (is.null(name) || !nzchar(name)) return(name)
    trimws(gsub("\\s+[A-Z]{1,4}\\([0-9:ZEzea,/()\\s-]+\\).*$", "", name))
  }
  
  startElement <- function(name, attrs) {
    tag <- strip_ns(name); env$cur_tag <- tag
    if (tag == "metabolite") {
      env$in_met <- TRUE; env$acc <- NULL; env$cname <- NULL
      env$acc_done <- env$name_done <- FALSE
      env$in_pathways <- env$in_pathway <- env$in_diseases <- env$in_disease <- FALSE
      env$in_biospec  <- env$in_tissue  <- FALSE
      env$pw_name <- env$dis_name <- env$bio_name <- env$tis_name <- NULL
      env$tax <- list()
    } else if (tag=="pathways"              && env$in_met)      { env$in_pathways <- TRUE
    } else if (tag=="pathway"               && env$in_pathways) { env$in_pathway  <- TRUE; env$pw_name <- NULL
    } else if (tag=="diseases"              && env$in_met)      { env$in_diseases <- TRUE
    } else if (tag=="disease"               && env$in_diseases) { env$in_disease  <- TRUE; env$dis_name <- NULL
    } else if (tag=="biospecimen_locations" && env$in_met)      { env$in_biospec  <- TRUE
    } else if (tag=="tissue_locations"      && env$in_met)      { env$in_tissue   <- TRUE }
  }
  
  text_handler <- function(text, ...) {
    text <- trimws(text)
    if (!nzchar(text) || !env$in_met) return()
    tag <- env$cur_tag
    if      (tag=="accession" && !env$acc_done  && !env$in_pathways && !env$in_diseases) {
      env$acc <- text; env$acc_done <- TRUE
    } else if (tag=="name" && !env$name_done && !env$in_pathways && !env$in_diseases &&
               !env$in_biospec && !env$in_tissue) {
      env$cname <- text; env$name_done <- TRUE
    } else if (tag=="name"        && env$in_pathway) { env$pw_name  <- text
    } else if (tag=="name"        && env$in_disease) { env$dis_name <- text
    } else if (tag=="biospecimen" && env$in_biospec)  { env$bio_name <- text
    } else if (tag=="tissue"      && env$in_tissue)   { env$tis_name <- text
    } else if (tag %in% c("kingdom","super_class","class","sub_class","direct_parent") &&
               env$in_met && !env$in_pathways && !env$in_diseases) {
      env$tax[[tag]] <- text }
  }
  
  endElement <- function(name) {
    tag <- strip_ns(name)
    if (tag=="pathway" && env$in_pathway) {
      if (!is.null(env$acc) && !is.null(env$pw_name)) {
        pw <- collapse_pw(env$pw_name)
        if (!is.null(pw) && nzchar(pw))
          pathways_map[[env$acc]] <<- unique(c(pathways_map[[env$acc]], pw))
      }
      env$in_pathway <- FALSE; env$pw_name <- NULL
    } else if (tag=="pathways")  { env$in_pathways <- FALSE
    } else if (tag=="disease" && env$in_disease) {
      if (!is.null(env$acc) && !is.null(env$dis_name))
        diseases_map[[env$acc]] <<- unique(c(diseases_map[[env$acc]], env$dis_name))
      env$in_disease <- FALSE; env$dis_name <- NULL
    } else if (tag=="diseases")  { env$in_diseases <- FALSE
    } else if (tag=="biospecimen" && env$in_biospec && !is.null(env$bio_name)) {
      if (!is.null(env$acc))
        biospecimen_map[[env$acc]] <<- unique(c(biospecimen_map[[env$acc]], env$bio_name))
      env$bio_name <- NULL
    } else if (tag=="biospecimen_locations") { env$in_biospec <- FALSE
    } else if (tag=="tissue" && env$in_tissue && !is.null(env$tis_name)) {
      if (!is.null(env$acc))
        tissue_map[[env$acc]] <<- unique(c(tissue_map[[env$acc]], env$tis_name))
      env$tis_name <- NULL
    } else if (tag=="tissue_locations") { env$in_tissue <- FALSE
    } else if (tag=="metabolite") {
      if (!is.null(env$acc)) {
        names_map[env$acc] <<- if (!is.null(env$cname)) env$cname else ""
        if (length(env$tax) > 0) taxonomy_map[[env$acc]] <<- env$tax
      }
      env$in_met <- FALSE; env$n <- env$n + 1L
      if (env$n %% 50000L == 0L)
        message(sprintf("  Parsed %d metabolites in %.0fs...",
                        env$n, proc.time()["elapsed"] - env$t0))
    }
    env$cur_tag <- NULL
  }
  
  xmlEventParse(HMDB_XML,
                handlers = list(startElement=startElement, text=text_handler, endElement=endElement),
                ignoreBlanks = TRUE)
  message(sprintf("  Done: %d metabolites in %.0fs", env$n, proc.time()["elapsed"]-env$t0))
}

# ── 3. Load annotation table & build background ───────────────────────────────
message("=== Step 3: Loading annotation table ===")
annot <- read.csv(ANNOTATION_FILE, stringsAsFactors = FALSE)

bg      <- annot[!is.na(annot$Level) & annot$Level == 1 &
                   !is.na(annot$Compound.name) & !is.na(annot$HMDB.ID), ]
bg      <- bg[!duplicated(bg$Compound.name), ]
bg_hmdb <- bg$HMDB.ID
N       <- length(bg_hmdb)
message(sprintf("  Background: %d unique Level-1 compounds with HMDB IDs", N))

# ── 4. Build annotation sets (once against full background) ───────────────────
message("=== Step 4: Building annotation sets ===")

build_sets <- function(data_list, bg_set) {
  sets <- list()
  for (acc in names(data_list)) {
    if (!(acc %in% bg_set)) next
    for (term in data_list[[acc]]) {
      term <- trimws(term)
      if (!nzchar(term)) next
      sets[[term]] <- unique(c(sets[[term]], acc))
    }
  }
  sets
}

build_taxonomy_sets <- function(taxonomy_map, bg_set, level) {
  sets <- list()
  for (acc in names(taxonomy_map)) {
    if (!(acc %in% bg_set)) next
    val <- taxonomy_map[[acc]][[level]]
    if (is.null(val) || length(val) == 0) next
    val <- trimws(val[1])
    if (!nzchar(val)) next
    sets[[val]] <- unique(c(sets[[val]], acc))
  }
  sets
}

all_annotation_types <- list(
  "HMDB_Pathway"            = build_sets(pathways_map,    bg_hmdb),
  "HMDB_Disease"            = build_sets(diseases_map,    bg_hmdb),
  "HMDB_Biospecimen"        = build_sets(biospecimen_map, bg_hmdb),
  "HMDB_Tissue"             = build_sets(tissue_map,      bg_hmdb),
  "ClassyFire_Kingdom"      = build_taxonomy_sets(taxonomy_map, bg_hmdb, "kingdom"),
  "ClassyFire_SuperClass"   = build_taxonomy_sets(taxonomy_map, bg_hmdb, "super_class"),
  "ClassyFire_Class"        = build_taxonomy_sets(taxonomy_map, bg_hmdb, "class"),
  "ClassyFire_SubClass"     = build_taxonomy_sets(taxonomy_map, bg_hmdb, "sub_class"),
  "ClassyFire_DirectParent" = build_taxonomy_sets(taxonomy_map, bg_hmdb, "direct_parent")
)

# ── 5. Define enrichment function ─────────────────────────────────────────────
run_enrichment <- function(mod_set, bg_set, annotation_sets, atype_label) {
  N_bg  <- length(bg_set)
  n_mod <- length(mod_set)
  rows  <- list()
  for (term in names(annotation_sets)) {
    bg_members <- annotation_sets[[term]]
    K  <- length(bg_members)
    if (K < MIN_PATH_SIZE || K > MAX_PATH_SIZE) next
    hits <- intersect(mod_set, bg_members)
    k    <- length(hits)
    mat  <- matrix(c(k, n_mod-k, K-k, N_bg-n_mod-K+k), nrow=2, byrow=TRUE)
    ft   <- fisher.test(mat, alternative="greater")
    a <- k; b <- n_mod-k; cc <- K-k; d <- N_bg-n_mod-K+k
    or_val <- if (b==0||cc==0) { if(k>0) Inf else NA_real_ } else (a*d)/(b*cc)
    hit_names <- names_map[hits]
    hit_names[is.na(hit_names)|!nzchar(hit_names)] <- hits[is.na(hit_names)|!nzchar(hit_names)]
    rows[[length(rows)+1L]] <- data.frame(
      annotation_type  = atype_label,
      term             = term,
      term_size_bg     = K,
      module_hits      = k,
      hit_fraction     = paste0(k,"/",K),
      hit_compounds    = paste(sort(hit_names), collapse="; "),
      odds_ratio       = round(or_val, 3),
      p_value          = ft$p.value,
      stringsAsFactors = FALSE)
  }
  if (length(rows) == 0L) return(data.frame())
  df             <- do.call(rbind, rows)
  df$fdr         <- p.adjust(df$p_value, method="BH")
  df$sig_nominal <- df$p_value < NOMINAL_ALPHA
  df$sig_fdr     <- df$fdr     < FDR_ALPHA
  df$neg_log10_p <- -log10(pmax(df$p_value, 1e-300))
  df[order(df$p_value), ]
}

# ── 6. Loop over all modules ───────────────────────────────────────────────────
message("=== Step 6: Running enrichment across all modules ===")

module_files <- list.files(MODULES_DIR, pattern = "_members\\.txt$", full.names = TRUE)
message(sprintf("  Found %d module files", length(module_files)))

all_modules_results <- list()

for (mf in module_files) {
  
  module_name <- sub("_members\\.txt$", "", basename(mf))
  message(sprintf("\n--- Module: %s ---", module_name))
  
  module_qi <- trimws(readLines(mf))
  module_qi <- module_qi[nzchar(module_qi)]
  
  mod_raw  <- annot[!is.na(annot$Level) & annot$Level == 1 & annot$name %in% module_qi, ]
  mod_raw  <- mod_raw[!is.na(mod_raw$Compound.name) & !is.na(mod_raw$HMDB.ID), ]
  mod      <- mod_raw[!duplicated(mod_raw$Compound.name), ]
  mod_hmdb <- mod$HMDB.ID
  n_mod    <- length(mod_hmdb)
  
  message(sprintf("  %d QI features -> %d unique annotated compounds", length(module_qi), n_mod))
  
  if (n_mod == 0) {
    message("  Skipping: no annotated compounds.")
    next
  }
  
  not_in_bg <- mod_hmdb[!(mod_hmdb %in% bg_hmdb)]
  if (length(not_in_bg) > 0) {
    warning(module_name, ": ", length(not_in_bg), " compounds not in background — excluded.")
    mod_hmdb <- mod_hmdb[mod_hmdb %in% bg_hmdb]
    n_mod    <- length(mod_hmdb)
  }
  
  # Coverage summary
  res_df <- data.frame()
  for (atype in names(all_annotation_types)) {
    n_terms   <- length(all_annotation_types[[atype]])
    n_covered <- sum(mod_hmdb %in% unlist(all_annotation_types[[atype]]))
    res_df    <- rbind(res_df, data.frame(
      annotation_type = atype,
      n_terms         = n_terms,
      n_covered       = n_covered,
      n_mod           = n_mod,
      stringsAsFactors = FALSE))
  }
  
  # Enrichment
  results_list <- list()
  for (atype in names(all_annotation_types)) {
    df <- run_enrichment(mod_hmdb, bg_hmdb, all_annotation_types[[atype]], atype)
    results_list[[atype]] <- df
    if (nrow(df) > 0)
      message(sprintf("  %-30s %4d tested | %3d p<0.05 | %3d FDR<0.05",
                      atype, nrow(df), sum(df$sig_nominal), sum(df$sig_fdr)))
  }
  
  results_all <- do.call(rbind, results_list)
  rownames(results_all) <- NULL
  
  if (nrow(results_all) == 0) {
    message("  No terms passed size filters — skipping.")
    next
  }
  
  results_all <- merge(results_all, res_df, by = "annotation_type")
  results_all$module <- module_name
  
  write.csv(results_all,
            file = file.path(RESULTS_DIR, paste0("enrichment_", module_name, "_all.csv")),
            row.names = FALSE)
  
  sig <- results_all[results_all$sig_nominal, ]
  write.csv(sig,
            file = file.path(RESULTS_DIR, paste0("enrichment_", module_name, "_sig.csv")),
            row.names = FALSE)
  
  all_modules_results[[module_name]] <- results_all
  message(sprintf("  -> Done: %s", module_name))
}

# ── 7. Combine & save master table ────────────────────────────────────────────
message("\n=== Step 7: Combining all modules into master table ===")

if (length(all_modules_results) == 0) {
  message("  WARNING: No results generated for any module.")
} else {
  master_df <- do.call(rbind, all_modules_results)
  rownames(master_df) <- NULL
  
  col_order <- c("module", setdiff(names(master_df), "module"))
  master_df <- master_df[, col_order]
  
  out_path <- file.path(RESULTS_DIR, "enrichment_ALL_MODULES.csv") 
  master_df$proportion_of_hits<-master_df$module_hits/master_df$term_size_bg
  write.csv(master_df, file = out_path, row.names = FALSE)
  out_path

}

message("\n=== Done! ===")