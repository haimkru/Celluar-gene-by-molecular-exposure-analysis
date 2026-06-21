###################
#/labs/smontgom/grps_smontgom/hkrupkin/GxE/Code/single_variant_example_tensorQTL/code_to_replicate_in_R
#author: Haim Krupkin using AI
#date: 04/09/2026
#description: we are reverse engineerign tensorQTL to be able to replicate a single interaction from the GxE analysis we did
##to get singel exmaple you need to extract all the data relvent to that isngel example
#expression in cell type
#covarites
#genotypes
#and the valeus of all the exposures
#to get the genoypes you use plink, specifically you use:
#plink --bfile /oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/genotypes/mesa_1331samples.maf01.biallelic.intersect --extract /oak/stanford/groups/smontgom/dnachun/projects/topmed/t
#ensorqtl_inputs/genotypes/variant_of_intrest.txt --recode A  --out /oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/genotypes/variant_of_intrest.txt


setwd("/labs/smontgom/grps_smontgom/hkrupkin/GxE/Code/single_variant_example_tensorQTL")
#!/usr/bin/env Rscript
#!/usr/bin/env Rscript
#!/usr/bin/env Rscript

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


# ============================================================
# SECTION 1: Setup
# ============================================================

covariates_path <- "nk_expression_covariates.txt"
variant_path <- "variant_of_intrest.txt.raw"
interaction_path <- "Eigen_Metabolite_Expression_interaction_chunk_1.txt"
expression_path <- "ENSG00000273100.1.txt"
reference_path <- "file_check.txt"

variant_id <- "chr6_169284690_C_T_b38"
phenotype_id <- "ENSG00000273100.1"
plot_environment <- "tan"

max_ratio_diff_ok <- 0.005

stats_out <- "interaction_summary_stats.tsv"
comparison_out <- "interaction_summary_reference_comparison.tsv"
residualized_out <- "interaction_residualized_inputs.tsv"
plot_out <- "interaction_plot.png"

cat("Working directory:", getwd(), "\n")

# STOP AND INSPECT:
# - Are these the file names you actually want?
# - Is plot_environment a shorthand like 'tan' or an exact column name?


# ============================================================
# SECTION 2: Lightweight helpers
# ============================================================

log_msg <- function(msg) cat(msg, "\n")

resolve_input_path <- function(path) {
  if (file.exists(path)) return(path)
  candidates <- c(paste0(path, ".txt"), paste0(path, ".txt.raw"), paste0(path, ".raw"))
  for (candidate in candidates) {
    if (file.exists(candidate)) return(candidate)
  }
  stop(sprintf("Input file not found: %s", path))
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

maybe_set_rownames <- function(df) {
  if (ncol(df) < 2) return(df)
  first <- tolower(colnames(df)[1])
  id_like <- c("sample", "sample_id", "iid", "id", "s", "phenotype_id", "variant_id", "gene_id")
  if (first %in% id_like) {
    rownames(df) <- as.character(df[[1]])
    df <- df[, -1, drop = FALSE]
  }
  df
}

coerce_numeric <- function(x, name) {
  out <- suppressWarnings(as.numeric(x))
  if (all(is.na(out))) stop(sprintf("All values are non-numeric in: %s", name))
  out
}

canonical_variant_id <- function(v) {
  sub("^(.*_[^_]+)_[ACGT]$", "\\1", as.character(v), perl = TRUE)
}

escape_regex <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}

# STOP AND INSPECT:
# - These are generic helpers: path resolution, reading, rownames, numeric conversion.
# - Nothing tensorQTL-specific has happened yet.


# ============================================================
# SECTION 3: Parse expression into one sample-indexed vector
# ============================================================

extract_expression_series <- function(expression_df, phenotype_id = NULL) {
  df <- maybe_set_rownames(expression_df)
  num_cols <- vapply(df, is.numeric, logical(1))
  
  if (ncol(df) == 1 && any(num_cols)) {
    pid <- if (!is.null(phenotype_id)) phenotype_id else colnames(df)[1]
    s <- coerce_numeric(df[[1]], "expression")
    names(s) <- rownames(df)
    log_msg("  - Expression parsed as samples x 1 column")
    return(list(series = s, phenotype_id = pid))
  }
  
  if (!is.null(phenotype_id) && phenotype_id %in% rownames(df)) {
    s <- coerce_numeric(df[phenotype_id, ], paste("expression row", phenotype_id))
    names(s) <- colnames(df)
    log_msg(sprintf("  - Expression parsed from row phenotype_id=%s", phenotype_id))
    return(list(series = s, phenotype_id = phenotype_id))
  }
  
  if (!is.null(phenotype_id) && phenotype_id %in% colnames(df)) {
    s <- coerce_numeric(df[[phenotype_id]], paste("expression column", phenotype_id))
    names(s) <- rownames(df)
    log_msg(sprintf("  - Expression parsed from column phenotype_id=%s", phenotype_id))
    return(list(series = s, phenotype_id = phenotype_id))
  }
  
  if (nrow(df) == 1) {
    pid <- if (!is.null(phenotype_id)) phenotype_id else rownames(df)[1]
    s <- coerce_numeric(df[1, ], "expression single-row")
    names(s) <- colnames(df)
    log_msg("  - Expression parsed as 1 phenotype row x samples")
    return(list(series = s, phenotype_id = pid))
  }
  
  stop("Could not infer expression layout")
}

expr_raw <- read_table(resolve_input_path(expression_path)) #loading data
expr_res <- extract_expression_series(expr_raw, phenotype_id) #loading more data
expr_s <- expr_res$series #still loading expression data
phenotype_id_resolved <- expr_res$phenotype_id #getting gene id

cat("Expression raw shape:", nrow(expr_raw), "x", ncol(expr_raw), "\n")
cat("Expression vector length:", length(expr_s), "\n")
print(head(expr_s))

# STOP AND INSPECT:
# - expr_s should be one numeric vector.
# - names(expr_s) should be sample IDs.
# - phenotype_id_resolved should be the phenotype you expect.


# ============================================================
# SECTION 4: Parse genotype into one sample-indexed vector
# ============================================================

extract_variant_series <- function(variant_df, variant_id = NULL, sample_hint = NULL) {
  df <- maybe_set_rownames(variant_df)
  cols_lower <- tolower(colnames(df))
  
  variant_col <- colnames(df)[match(TRUE, cols_lower %in% c("variant_id", "variant"), nomatch = 0)]
  sample_col <- colnames(df)[match(TRUE, cols_lower %in% c("sample", "sample_id", "iid"), nomatch = 0)]
  geno_col <- colnames(df)[match(TRUE, cols_lower %in% c("genotype", "dosage", "gt"), nomatch = 0)]
  
  if (length(variant_col) == 1 && length(sample_col) == 1 && length(geno_col) == 1) {
    if (is.null(variant_id)) {
      uniq <- unique(as.character(df[[variant_col]]))
      if (length(uniq) != 1) stop("Variant file has multiple variants; set variant_id")
      variant_id <- uniq[1]
    }
    sdf <- df[as.character(df[[variant_col]]) == variant_id, , drop = FALSE]
    s <- coerce_numeric(sdf[[geno_col]], paste("genotype", variant_id))
    names(s) <- as.character(sdf[[sample_col]])
    log_msg("  - Variant parsed from long format")
    return(list(series = s, variant_id = variant_id))
  }
  
  plink_sample_col <- colnames(df)[match(TRUE, cols_lower %in% c("iid", "sample", "sample_id", "id"), nomatch = 0)]
  if (length(plink_sample_col) == 0 && "IID" %in% colnames(df)) plink_sample_col <- "IID"
  if (length(plink_sample_col) == 0 && "FID" %in% colnames(df)) plink_sample_col <- "FID"
  
  if (length(plink_sample_col) == 1) {
    meta_cols <- c("fid", "iid", "pat", "mat", "sex", "phenotype", "pheno", "id", "sample", "sample_id")
    candidate_cols <- colnames(df)[!(tolower(colnames(df)) %in% meta_cols)]
    chosen_col <- NULL
    
    if (!is.null(variant_id)) {
      if (variant_id %in% colnames(df)) {
        chosen_col <- variant_id
      } else {
        by_prefix <- candidate_cols[grepl(paste0("^", escape_regex(variant_id), "_"), candidate_cols)]
        if (length(by_prefix) == 1) chosen_col <- by_prefix
        if (length(by_prefix) > 1) stop(sprintf("Ambiguous variant match: %s", paste(by_prefix, collapse = ", ")))
        if (is.null(chosen_col)) {
          by_contains <- candidate_cols[grepl(escape_regex(variant_id), candidate_cols)]
          if (length(by_contains) == 1) chosen_col <- by_contains
          if (length(by_contains) > 1) stop(sprintf("Ambiguous variant match: %s", paste(by_contains, collapse = ", ")))
        }
      }
    } else if (length(candidate_cols) == 1) {
      chosen_col <- candidate_cols
    }
    
    if (!is.null(chosen_col)) {
      s <- coerce_numeric(df[[chosen_col]], paste("genotype column", chosen_col))
      names(s) <- as.character(df[[plink_sample_col]])
      out_variant_id <- if (!is.null(variant_id)) variant_id else canonical_variant_id(chosen_col)
      log_msg(sprintf("  - Variant parsed from PLINK raw style column=%s", chosen_col))
      return(list(series = s, variant_id = out_variant_id))
    }
  }
  
  if (!is.null(variant_id) && variant_id %in% rownames(df)) {
    s <- coerce_numeric(df[variant_id, ], paste("variant row", variant_id))
    names(s) <- colnames(df)
    log_msg(sprintf("  - Variant parsed from row variant_id=%s", variant_id))
    return(list(series = s, variant_id = variant_id))
  }
  
  stop("Could not infer variant layout")
}

var_raw <- read_table(resolve_input_path(variant_path)) #this helps load the varaint data and parse it into the correct format
var_res <- extract_variant_series(var_raw, variant_id, names(expr_s)) #continue from above
geno_s <- var_res$series #continue from above
variant_id_resolved <- var_res$variant_id #continue from above 

cat("Variant raw shape:", nrow(var_raw), "x", ncol(var_raw), "\n")
cat("Genotype vector length:", length(geno_s), "\n")
cat("Resolved variant ID:", variant_id_resolved, "\n")
print(table(geno_s, useNA = "ifany"))

# STOP AND INSPECT:
# - geno_s should usually look like dosages 0/1/2, maybe with -9 for missing.
# - names(geno_s) should be sample IDs.
# - variant_id_resolved should match the reference row naming convention.


# ============================================================
# SECTION 5: Parse interaction and covariate matrices
# ============================================================

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

int_raw <- read_table(resolve_input_path(interaction_path)) #we load the enviromenta data
cov_raw <- read_table(resolve_input_path(covariates_path)) #we load the covariate data

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
dim(expr_s)
dim(geno_s)
dim(interaction_df)
dim(cov_df)
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

requested_env <- tolower(plot_environment)
env_matches <- colnames(interaction_df)[grepl(requested_env, tolower(colnames(interaction_df)), fixed = TRUE)]
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
write.table(stats_df, stats_out, sep = "\t", quote = FALSE, row.names = FALSE)
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
write.table(residualized_df, residualized_out, sep = "\t", quote = FALSE, row.names = FALSE)
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

if (file.exists("/labs/smontgom/grps_smontgom/hkrupkin/GxE/Code/single_variant_example_tensorQTL/file_check.txt")) {
  comp_df <- compare_with_reference(stats_df, reference_path, phenotype_id_resolved, variant_id_resolved, max_ratio_diff_ok)
  write.table(comp_df, comparison_out, sep = "\t", quote = FALSE, row.names = FALSE)
  print(head(comp_df[order(-comp_df$abs_diff), ], 15))
  cat("Ratio pass count:", sum(comp_df$within_ratio_diff_limit), "/", nrow(comp_df), "\n")
} else {
  comp_df <- NULL
  cat("Reference file not found; comparison skipped.\n")
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
    conf_level = 0.95
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Install with install.packages('ggplot2').")
  }
  
  y <- residualized_df$phenotype_centered_residualized
  env <- residualized_df[[paste0(env_name, "_centered_residualized")]]
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
      x = paste0("Residualized ", env_name),
      y = paste0("Residualized ", gene_name, " Expression")
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
  
  if (!is.null(out_path) && nzchar(out_path)) {
    ggplot2::ggsave(filename = out_path, plot = p, width = 8, height = 5.6, dpi = 150)
  }
  
  p
}

env_name <- resolve_environment_name(plot_environment, colnames(interaction_df))
pval_col <- paste0("pval_g-", env_name)
interaction_pval <- if (pval_col %in% colnames(stats_df)) as.numeric(stats_df[[pval_col]][1]) else NA_real_
g_pval <- if ("pval_g" %in% colnames(stats_df)) as.numeric(stats_df[["pval_g"]][1]) else NA_real_
e_pval_col <- paste0("pval_", env_name)
e_pval <- if (e_pval_col %in% colnames(stats_df)) as.numeric(stats_df[[e_pval_col]][1]) else NA_real_
gxe_pval <- interaction_pval

interaction_plot <- make_interaction_plot(
  residualized_df,
  geno_s,
  env_name,
  phenotype_id_resolved,
  plot_out,
  g_pval = g_pval,
  e_pval = e_pval,
  gxe_pval = gxe_pval,
  conf_level = 0.95
)
print(interaction_plot)
cat("Plot written:", plot_out, "\n")

ggsave(plot=interaction_plot,
       file="/labs/smontgom/grps_smontgom/hkrupkin/GxE/plots/plots_for_figures/single_exanmple_chr6_169284690_C_T_b38_T_tan_nk_cell_ESNG00000273100.1.png",
       width=4,
       height=4,
       dpi=900)

# STOP AND INSPECT:
# - The lines come from the residualized model-scale variables.
# - The interaction p-value displayed should match stats_df.


# ============================================================
# SECTION 14: Final recap
# ============================================================

cat("\nOutputs generated:\n")
cat("-", stats_out, "\n")
cat("-", residualized_out, "\n")
if (!is.null(comp_df)) cat("-", comparison_out, "\n")
cat("-", plot_out, "\n")

cat("\nWhat this taught you:\n")
cat("1. How raw tables become aligned model vectors/matrices.\n")
cat("2. How QR residualization removes covariate effects.\n")
cat("3. How the design matrix for y ~ g + i + g:i is assembled.\n")
cat("4. How beta, se, t, and p are computed from matrix algebra.\n")
cat("5. How the final plot relates to the residualized model inputs.\n")

#

