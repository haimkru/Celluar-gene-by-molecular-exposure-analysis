###predictions_of_hits_methylation#####
#/oak/stanford/groups/smontgom/hkrupkin/GxE/Code/predictions_of_hits_methylation
#author: Haim Krupkin
#date: 03/21/2026
#date of improvment: 04/27/2026
#description: The idea is to prove a prediction model is able to use transfer learning to predict if
#a gene will be a GxE in another cell type based on its expression and the characteristics of the
#module\eigenFactor and of the gene.


#/oak/stanford/groups/smontgom/hkrupkin/GxE/Code/predictions_of_hits
##to calcualte average expression in a cell type we use:
#



################################################################################
# CEG INTERACTION MODEL: ZERO-SHOT CELL TYPE PREDICTION
# ============================================================================
# Framework: Cell-Environment-Gene (CEG) triplet classification
# Design:    Leave-One-Cell-Type-Out (LOCTO) -- model is never shown cell labels,
#            only expression vectors. Zero-shot prediction on unseen cell types.
# Metric:    AUPRC (primary), AUROC (secondary), Lift over random baseline
#
# Expression input: CIBERSORT deconvolution output (wide format)
#   * Columns: Name (Ensembl ID with version), bcell, cd4t, cd8t, mono
#   * Values:  deconvolution coefficients -- can be negative (below-noise)
#   * Transform: pmax(score, 0) -> log1p  [negatives treated as 0]
#   * is_expressed: raw score > 0
#
# DATA LEAKAGE FIX (vs prior version):
#   module_stats, gene_stats, positive_pairs, ALL_GENES, ALL_MODULES, negative
#   sampling, and all feature merges are now recomputed INSIDE the LOCTO loop
#   using only train_hits (all_hits filtered to cell_code != test_cell).
#   The test cell's hits no longer leak into gene_resp_rank or log_module_potency.
################################################################################


# ==============================================================================
# SECTION 0: INSTALL AND LOAD DEPENDENCIES
# ==============================================================================
cran_pkgs <- c("data.table", "xgboost", "PRROC", "ggplot2",
               "stringr", "scales", "patchwork")
library(dplyr)
new_pkgs <- cran_pkgs[!cran_pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs)) {
  message("Installing: ", paste(new_pkgs, collapse = ", "))
  install.packages(new_pkgs, repos = "https://cloud.r-project.org")
}
invisible(lapply(cran_pkgs, library, character.only = TRUE))


# ==============================================================================
# SECTION 1: LOAD AND PREPARE all_hits
# ==============================================================================
all_hits <- fread("/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/all_cell_types_sig_hits.tsv")
b_all_hits<-all_hits

#all_hits <- all_hits %>% dplyr::filter(Omics == "Metabolome")%>%dplyr::filter()
all_hits <- all_hits %>%
  #before 04/27/2026 this was not only interaction, and is thus kinda wrong. right? yeah...
  dplyr::filter(effect_type=="interaction_effect")%>%
  dplyr::select(phenotype_id, term, cell_type, Omics, assay) %>% distinct()

head(all_hits)

#all_hits <- all_hits %>% dplyr::filter(cell_type != "nk")
stopifnot("all_hits must be loaded in your session" = exists("all_hits"))
setDT(all_hits)

message("all_hits: ", nrow(all_hits), " rows")
message("Cell types found: ", paste(unique(all_hits$cell_type), collapse = " | "))

CELL_TYPE_MAP <- c(
  "B Cells" = "bcell", "B cells" = "bcell",
  "b cells" = "bcell", "bcell"   = "bcell",
  
  "CD4 T Cells" = "cd4", "CD4 T cells" = "cd4",
  "cd4 t cells" = "cd4", "CD4" = "cd4",
  "CD4⁺ T" = "cd4",   # ✅ fix
  
  "CD8 T Cells" = "cd8", "CD8 T cells" = "cd8",
  "cd8 t cells" = "cd8", "CD8" = "cd8",
  "CD8⁺ T" = "cd8",   # ✅ fix
  
  "Monocytes" = "mono", "monocytes" = "mono",
  "Mono" = "mono", "mono" = "mono",
  
  "NK" = "nk"  # ✅ optional: include NK if you want it
)

all_hits[, cell_code := CELL_TYPE_MAP[cell_type]]
bad_cells <- unique(all_hits$cell_type[is.na(all_hits$cell_code)])
if (length(bad_cells))
  warning("Unmapped cell types: ", paste(bad_cells, collapse = ", "))

all_hits[, gene_id := sub("\\..*", "", phenotype_id)]
all_hits[, module  := sub("_Eigen$", "", term)]

CELL_CODES <- sort(unique(na.omit(all_hits$cell_code)))
message("Standardized cell codes: ", paste(CELL_CODES, collapse = ", "))


# ==============================================================================
# SECTION 2: LOAD CIBERSORT EXPRESSION DATA
# ==============================================================================
message("\n>>> Loading CIBERSORT expression data...")
#expression file is created by:
#/oak/stanford/groups/smontgom/hkrupkin/GxE/Code/summary_tpms_for_transcriptome_GxE.R
first_cell_type <- fread("/oak/stanford/groups/smontgom/hkrupkin/GxE/Code/average_methylation_in_diffrent_cell_types.csv")
cibersort_wide  <- copy(first_cell_type)
setDT(cibersort_wide)
GENE_NAME_COL <- "Name"
CIBER_COL_MAP <- c("bcell" = "bcell", "cd4t" = "cd4", "cd8t" = "cd8", "mono" = "mono","nk"="nk")

missing_cols <- setdiff(names(CIBER_COL_MAP), names(cibersort_wide))


message("CIBERSORT columns confirmed: ", paste(names(CIBER_COL_MAP), collapse = ", "))
message("Genes in CIBERSORT table   : ", nrow(cibersort_wide))


# ==============================================================================
# SECTION 3: RESHAPE AND TRANSFORM CIBERSORT DATA
# ==============================================================================
GENE_NAME_COL<-"phenotype_id"
expr_long <- melt(cibersort_wide, id.vars = GENE_NAME_COL,
                  measure.vars  = names(CIBER_COL_MAP),
                  variable.name = "ciber_col", value.name = "raw_score")
setDT(expr_long)
expr_long[, cell_code := CIBER_COL_MAP[as.character(ciber_col)]]
expr_long[, ciber_col := NULL]
expr_long[, gene_id   := sub("\\..*", "", get(GENE_NAME_COL))]
expr_long[, (GENE_NAME_COL) := NULL]
expr_long[, `:=`(is_expressed = as.integer(raw_score > 0),
                 log_expr     = log1p(pmax(raw_score, 0)))]

expr_long<-expr_long%>%
  dplyr::filter(!is.na(is_expressed))

message("\nExpression coverage after winsorize+log1p:")
print(expr_long[, .(n_genes       = uniqueN(gene_id),
                    pct_expressed = round(mean(is_expressed) * 100, 1)),
                by = cell_code][order(cell_code)])

# Safe outside loop: derived from expression data, not from hits
cohort_expr      <- expr_long[, .(gene_id, cell_code, raw_score, log_expr, is_expressed)]
cell_median_expr <- cohort_expr[, .(median_log_expr = median(log_expr, na.rm = TRUE)),
                                by = cell_code]

missing_cells <- CELL_CODES[!CELL_CODES %in% unique(cohort_expr$cell_code)]
if (length(missing_cells))
  warning("Cell types missing from CIBERSORT: ", paste(missing_cells, collapse = ", "),
          "\n  -> These cells get log_expr=0, is_expressed=0.")

# Feature names are fixed across all folds
FEATURE_COLS <- c(
  "log_expr",           # Gene expression in this cell [KEY ZERO-SHOT FEATURE]
  "is_expressed",       # Binary: CIBERSORT coeff > 0
  "gene_resp_rank",     # How broadly responsive this gene is (0-1 percentile)
  "gene_n_modules",     # Distinct modules the gene hits in training
  "gene_n_cells",       # Cell types in which gene is significant in training
  "log_module_potency", # Log(total hits) for this module across training
  "module_n_genes",     # Distinct genes this module perturbs in training
  "module_n_cells",     # Cell types in which this module is active in training
  "expr_x_potency"      # Interaction: log_expr x log_module_potency
)


# ==============================================================================
# SECTIONS 4+5: LOCTO LOOP
# Everything derived from all_hits is recomputed here on train_hits only.
# The test cell's hits never touch any feature seen during training/evaluation.
# ==============================================================================
message("\n>>> LOCTO evaluation -- ", length(CELL_CODES), " folds\n")
locto_results        <- vector("list", length(CELL_CODES))
names(locto_results) <- CELL_CODES

for (i in seq_along(CELL_CODES)) {
  
  test_cell  <- CELL_CODES[i]
  message(sprintf("[%d/%d] Hold-out: %-8s", i, length(CELL_CODES), test_cell))
  
  # -- 4a. Split ---------------------------------------------------------------
  train_hits <- all_hits[!is.na(cell_code) & cell_code != test_cell]
  test_hits  <- all_hits[!is.na(cell_code) & cell_code == test_cell]
  
  # -- 4b. Module stats from training data only --------------------------------
  module_stats <- train_hits[, .(module_potency = .N,
                                 module_n_genes  = uniqueN(gene_id),
                                 module_n_cells  = uniqueN(cell_code)),
                             by = module]
  module_stats[, log_module_potency := log1p(module_potency)]
  
  # -- 4c. Gene stats from training data only ----------------------------------
  gene_stats <- train_hits[, .(gene_total_hits = .N,
                               gene_n_modules  = uniqueN(module),
                               gene_n_cells    = uniqueN(cell_code)),
                           by = gene_id]
  gene_stats[, gene_resp_rank := rank(gene_total_hits, ties.method = "average") / .N]
  
  # -- 4d. Positive pairs ------------------------------------------------------
  train_pos <- unique(train_hits[, .(gene_id, module, cell_code, hit = 1L)])
  test_pos  <- unique(test_hits[,  .(gene_id, module, cell_code, hit = 1L)])
  message(sprintf("  Train positives: %d | Test positives: %d",
                  nrow(train_pos), nrow(test_pos)))
  
  if (nrow(test_pos) == 0) {
    warning("  No positives for ", test_cell, " -- skipping.")
    next
  }
  
  # -- 4e. Negative sampling from training universe only -----------------------
  TRAIN_GENES   <- unique(train_hits$gene_id)
  TRAIN_MODULES <- unique(train_hits$module)
  
  .gen_neg <- function(pos_dt, tgt_cell, g_univ, m_univ, ratio = 5L) {
    pos_keys <- pos_dt[cell_code == tgt_cell, paste(gene_id, module)]
    n_need   <- length(pos_keys) * ratio
    neg      <- character(0)
    iter     <- 0L
    while (length(neg) < n_need && iter < 25L) {
      iter  <- iter + 1L
      n_try <- (n_need - length(neg)) * 4L
      cands <- paste(sample(g_univ, n_try, replace = TRUE),
                     sample(m_univ, n_try, replace = TRUE))
      neg   <- unique(c(neg, cands[!cands %in% c(pos_keys, neg)]))
    }
    neg   <- neg[seq_len(min(n_need, length(neg)))]
    parts <- strsplit(neg, " ")
    data.table(gene_id   = vapply(parts, `[`, character(1), 1L),
               module    = vapply(parts, `[`, character(1), 2L),
               cell_code = tgt_cell, hit = 0L)
  }
  
  set.seed(2024L + i)
  train_neg <- rbindlist(lapply(unique(train_pos$cell_code), .gen_neg,
                                pos_dt = train_pos,
                                g_univ = TRAIN_GENES, m_univ = TRAIN_MODULES))
  test_neg  <- .gen_neg(test_pos, test_cell, TRAIN_GENES, TRAIN_MODULES)
  
  train_pairs <- rbind(train_pos, train_neg)
  test_pairs  <- rbind(test_pos,  test_neg)
  message(sprintf("  Train: %d pos + %d neg | Test: %d pos + %d neg",
                  nrow(train_pos), nrow(train_neg),
                  nrow(test_pos),  nrow(test_neg)))
  
  # -- 4f. Attach expression (safe outside loop, no hit leakage) ---------------
  .attach_expr <- function(dt) {
    dt <- merge(dt, cohort_expr[, .(gene_id, cell_code, log_expr, is_expressed)],
                by = c("gene_id", "cell_code"), all.x = TRUE)
    dt <- merge(dt, cell_median_expr, by = "cell_code", all.x = TRUE)
    dt[is.na(log_expr),     log_expr     := fifelse(is.na(median_log_expr), 0, median_log_expr)]
    dt[is.na(is_expressed), is_expressed := 0L]
    dt[, median_log_expr := NULL]
    dt
  }
  train_pairs <- .attach_expr(train_pairs)
  test_pairs  <- .attach_expr(test_pairs)
  
  # -- 4g. Attach training-derived statistics ----------------------------------
  .attach_stats <- function(dt) {
    dt <- merge(dt, gene_stats,   by = "gene_id", all.x = TRUE)
    dt <- merge(dt, module_stats, by = "module",  all.x = TRUE)
    for (sc in c("gene_total_hits","gene_n_modules","gene_n_cells","gene_resp_rank",
                 "module_potency","module_n_genes","module_n_cells","log_module_potency"))
      dt[is.na(get(sc)), (sc) := 0]
    dt[, expr_x_potency := log_expr * log_module_potency]
    dt
  }
  train_pairs <- .attach_stats(train_pairs)
  test_pairs  <- .attach_stats(test_pairs)
  
  # -- 5. Train XGBoost --------------------------------------------------------
  # Early stopping leakage fix:
  #   We carve a validation set (20%) out of train_pairs ONLY, stratified by
  #   cell type so every training cell contributes to the val set proportionally.
  #   dtest is never shown to the training loop -- it is used for evaluation only.
  set.seed(42L + i)
  val_idx     <- train_pairs[, .I[sample(.N, max(1L, round(.N * 0.20)))],
                             by = cell_code]$V1
  val_pairs   <- train_pairs[val_idx]
  fit_pairs   <- train_pairs[-val_idx]
  
  X_fit   <- as.matrix(fit_pairs[,  ..FEATURE_COLS])
  y_fit   <- fit_pairs$hit
  X_val   <- as.matrix(val_pairs[,  ..FEATURE_COLS])
  y_val   <- val_pairs$hit
  X_test  <- as.matrix(test_pairs[, ..FEATURE_COLS])
  y_test  <- test_pairs$hit
  
  n_pos <- sum(y_fit == 1); n_neg <- sum(y_fit == 0)
  spw   <- n_neg / max(n_pos, 1)
  message(sprintf("  Fit: %d pos / %d neg | Val: %d rows | spw=%.1f",
                  n_pos, n_neg, nrow(val_pairs), spw))
  
  dfit  <- xgb.DMatrix(X_fit,  label = y_fit)
  dval  <- xgb.DMatrix(X_val,  label = y_val)
  dtest <- xgb.DMatrix(X_test, label = y_test)
  
  model <- xgb.train(
    params = list(objective = "binary:logistic", eval_metric = c("aucpr","logloss"),
                  eta = 0.05, max_depth = 5, subsample = 0.80,
                  colsample_bytree = 0.80, min_child_weight = 10,
                  gamma = 0.10, lambda = 1.0, scale_pos_weight = spw, seed = 42L),
    data                  = dfit,
    nrounds               = 500,
    watchlist             = list(train = dfit, val = dval),  # dtest never seen here
    verbose               = 0,
    early_stopping_rounds = 40,
    print_every_n         = 100
  )
  
  preds   <- predict(model, dtest)
  pr_obj  <- pr.curve(scores.class0 = preds[y_test == 1],
                      scores.class1 = preds[y_test == 0], curve = TRUE)
  roc_obj <- roc.curve(scores.class0 = preds[y_test == 1],
                       scores.class1 = preds[y_test == 0], curve = TRUE)
  prevalence <- mean(y_test)
  lift       <- pr_obj$auc.integral / prevalence
  
  message(sprintf("  AUPRC=%.4f | Baseline=%.4f | Lift=%.2fx | AUROC=%.4f | Rounds=%d",
                  pr_obj$auc.integral, prevalence, lift,
                  roc_obj$auc, model$best_iteration))
  
  locto_results[[test_cell]] <- list(
    cell_type  = test_cell, model = model, feature_cols = FEATURE_COLS,
    preds = preds, labels = y_test,
    pr_curve = pr_obj, roc_curve = roc_obj,
    auprc = pr_obj$auc.integral, auroc = roc_obj$auc,
    prevalence = prevalence, lift = lift,
    best_round = model$best_iteration,
    n_train_pos = n_pos, n_train_neg = n_neg,
    n_test_pos  = sum(y_test == 1), n_test_neg = sum(y_test == 0)
  )
}

locto_results <- Filter(Negate(is.null), locto_results)


# ==============================================================================
# SECTION 6: SUMMARY TABLE
# ==============================================================================
summary_dt <- rbindlist(lapply(locto_results, function(r) {
  data.table(cell_type        = r$cell_type,
             n_test_pos       = r$n_test_pos,
             n_test_total     = r$n_test_pos + r$n_test_neg,
             prevalence       = round(r$prevalence, 5),
             auprc            = round(r$auprc, 4),
             auroc            = round(r$auroc, 4),
             lift_over_random = round(r$lift, 2),
             best_round       = r$best_round)
}))

message("\n", strrep("=", 70))
message("LOCTO EVALUATION SUMMARY")
message(strrep("=", 70))
print(summary_dt)
message(strrep("-", 70))
message(sprintf("Mean AUPRC across folds : %.4f", mean(summary_dt$auprc)))
message(sprintf("Mean lift over baseline : %.2fx", mean(summary_dt$lift_over_random)))
message(strrep("=", 70))


# ==============================================================================
# SECTION 7: VISUALIZATION
# ==============================================================================
message("\n>>> Generating visualizations...")
pr_df <- rbindlist(lapply(locto_results, function(r) {
  cd <- as.data.table(r$pr_curve$curve)
  setnames(cd, c("recall", "precision", "threshold"))
  cd[, cell_type := r$cell_type][, label := sprintf("%s (%.3f)", r$cell_type, r$auprc)]
}))


imp_all <- rbindlist(lapply(locto_results, function(r) {
  imp <- as.data.table(xgb.importance(model = r$model))
  imp[, cell_type := r$cell_type]; imp
}))
imp_avg <- imp_all[, .(mean_gain = mean(Gain), sd_gain = sd(Gain)),
                   by = Feature][order(-mean_gain)]

summary_dt<-summary_dt%>%
  #dplyr::filter(cell_type!="nk")%>%
  mutate(cell_type=ifelse(cell_type=="mono","Mono",cell_type))%>%
  mutate(cell_type=ifelse(cell_type=="nk","NK",cell_type))%>%
  
  mutate(cell_type=ifelse(cell_type=="cd8","CD8+ T",cell_type))%>%
  mutate(cell_type=ifelse(cell_type=="cd4","CD4+ T",cell_type))%>%
  mutate(cell_type=ifelse(cell_type=="bcell","B Cells",cell_type))


saveRDS(summary_dt,"/labs/smontgom/grps_smontgom/hkrupkin/GxE/plots/plots_for_figures/lift_over_simple_modelsV2_methylation.rds")
summary_dt<-readRDS("/labs/smontgom/grps_smontgom/hkrupkin/GxE/plots/plots_for_figures/lift_over_simple_modelsV2_methylation.rds")

p4 <- ggplot(summary_dt, aes(x = reorder(cell_type, lift_over_random),
                             y = lift_over_random, fill = cell_type)) +
  geom_col(alpha = 0.85, width = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "firebrick", linewidth = 0.9) +
  geom_text(aes(label = paste0(lift_over_random, "x")), hjust = 1, size = 4.2) +
  scale_fill_brewer(palette = "Dark2", guide = "none") + coord_flip(clip = "off") +
  labs(
    x = "Held-Out Cell Type", y = "Lift") + theme_bw()
p4
ggsave(plot=p4,
       file="/labs/smontgom/grps_smontgom/hkrupkin/GxE/plots/plots_for_figures/lift_over_simple_modelsV2_methylation.png",
       width=3,
       height=3,
       dpi=300
       
)
