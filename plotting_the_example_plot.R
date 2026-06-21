################################################################################
# Replicated Single Interaction GxE Pipeline (Faithful Mirror of Expansion Pack)
# Target: ENSG00000180660.8 | chr13_34629200_C_G_b38 | MEbisque4 (B cell)
################################################################################

library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(data.table)

# ============================================================
# SECTION 1: Paths & Manual Override Injection
# ============================================================
base_dir   <- "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Code/single_variant_example_tensorQTL"
covar_dir  <- "/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/Ready_for_running_inputs/covariatesV2"
pheno_dir  <- "/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/Ready_for_running_inputs/standardized_phenotypesV2"
setwd(base_dir)

# Build target exactly as structured in data.table rows
#target <- data.table(
#  variant_id   = "chr13_34629200_C_G_b38",
#  phenotype_id = "ENSG00000180660.8",
#  cell_join    = "B cell",          
#  assay_join   = "expression",
#  omics_join   = "metabolome",      
#  term         = "MEbisque4",       
#  path         = "chunk_1",
#  effect_type  = "interaction_effect"
#)

#target <- data.table(
#  variant_id   = "chr1_148438355_G_A_b38",
#  phenotype_id = "ENSG00000272824.1",
#  cell_join    = "B cell",          
#  assay_join   = "expression",
#  omics_join   = "metabolome",      
#  term         = "MEdarkred",       
#  path         = "chunk_2",
#  effect_type  = "interaction_effect"
#)



# 1. Create a target variant file
#echo "chr13_34629200_C_G_b38" > /labs/smontgom/grps_smontgom/hkrupkin/GxE/Code/single_variant_example_tensorQTL/single_target_variant.txt

# 2. Extract dosage data using PLINK
#module load plink && plink \
#--bfile /oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/genotypes/mesa_1331samples.maf01.biallelic.intersect \
#--extract /labs/smontgom/grps_smontgom/hkrupkin/GxE/Code/single_variant_example_tensorQTL/single_target_variant.txt \
#--recode A \
#--out /labs/smontgom/grps_smontgom/hkrupkin/GxE/Code/single_variant_example_tensorQTL/single_variant

#bcell_metabolome_chunk_3_expression_MEivory_trans_qtl.tsv.gz    conv=-0.0491    ENSG00000108474.17      chr17_16319648_G_A_b38  0.06561922      65      71      8.710393183412991e-13  -
#  0.32864434      0.044719834     0.961166328064871       -0.016751472    0.3438629       7.129686781328784e-09   -6.697511       1.1363434
target <- data.table(
  variant_id   = "chr1_148438355_G_A_b38",
  phenotype_id = "ENSG00000272824.1",
  cell_join    = "B cell",          
  assay_join   = "expression",
  omics_join   = "metabolome",      
  term         = "MEdarkred",       
  path         = "chunk_2",
  effect_type  = "interaction_effect"
)



v_id       <- target$variant_id     
p_id       <- target$phenotype_id   
assay      <- target$assay_join     
term_id    <- target$term           

# Clean Cell Name (using your exact translation rules)
cell_clean <- tolower(target$cell_join)
cell_clean <- gsub(" ", "", cell_clean)
cell_clean <- gsub("[⁺+]", "", cell_clean)
cell_clean <- gsub("neutro", "neu", cell_clean)

geno_raw_file <- file.path(base_dir, "single_variant.raw")

# Load gene annotations exactly as specified
gene_map <- fread("/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/gene_map.csv", header=TRUE)
colnames(gene_map) <- c("n", "n2", "gene_id", "gene_name")

# ============================================================
# SECTION 2: Exact Helper Functions From Your Code
# ============================================================
get_symbol <- function(ensembl_id, map_dt) {
  symbol <- map_dt[gene_id == ensembl_id, gene_name]
  if (length(symbol) == 0) {
    clean_id <- gsub("\\..*", "", ensembl_id)
    symbol <- map_dt[gsub("\\..*", "", gene_id) == clean_id, gene_name]
  }
  return(if(length(symbol) > 0) symbol[1] else NA)
}

extract_variant_series <- function(variant_df, variant_id) {
  col_idx <- grep(paste0("^", variant_id), colnames(variant_df))
  if (length(col_idx) == 0) stop(paste("Could not find variant column for", variant_id))
  return(list(
    series = setNames(as.numeric(variant_df[[col_idx[1]]]), variant_df$IID),
    variant_id = variant_id
  ))
}

maybe_set_rownames <- function(df) {
  if (!is.null(df[[1]]) && is.character(df[[1]])) {
    rownames(df) <- df[[1]]
    df <- df[, -1, drop = FALSE]
  }
  df
}

extract_interaction_df <- function(interaction_df, sample_hint) {
  df <- maybe_set_rownames(interaction_df)
  idx_overlap <- sum(rownames(df) %in% sample_hint)
  col_overlap <- sum(colnames(df) %in% sample_hint)
  
  if (col_overlap > idx_overlap && col_overlap >= 3) {
    df <- as.data.frame(t(df), stringsAsFactors = FALSE, check.names = FALSE)
    message("  - Interaction table transposed to samples x environments")
  }
  
  num_df <- data.frame(lapply(df, function(x) suppressWarnings(as.numeric(x))),
                       row.names = rownames(df), check.names = FALSE)
  keep <- colSums(!is.na(num_df)) > 0
  num_df <- num_df[, keep, drop = FALSE]
  if (ncol(num_df) == 0) stop("No numeric interaction/environment columns found")
  
  message(sprintf("  - Interaction parsed with %d environment(s)", ncol(num_df)))
  num_df
}

extract_covariates_df <- function(cov_df, sample_hint) {
  df <- maybe_set_rownames(cov_df)
  idx_overlap <- sum(rownames(df) %in% sample_hint)
  col_overlap <- sum(colnames(df) %in% sample_hint)
  
  if (col_overlap > idx_overlap && col_overlap >= 3) {
    df <- as.data.frame(t(df), stringsAsFactors = FALSE, check.names = FALSE)
    message("  - Covariate table transposed to samples x covariates")
  }
  
  num_df <- data.frame(lapply(df, function(x) suppressWarnings(as.numeric(x))),
                       row.names = rownames(df), check.names = FALSE)
  keep <- colSums(!is.na(num_df)) > 0
  num_df <- num_df[, keep, drop = FALSE]
  if (ncol(num_df) == 0) return(data.frame(row.names = rownames(df)))
  
  message(sprintf("  - Covariates parsed with %d covariate(s)", ncol(num_df)))
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
  omics_map <- c("metabolome" = "Metabolite", "proteins" = "Proteins", "proteome" = "Proteins")
  assay_map <- c("expression" = "Expression", "methylation" = "Methylation")
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
  
  return(dt[, .(interaction_file, interaction_path, file_exists)])
}

# ============================================================
# SECTION 3: Stream Phenotypes & Data Loading
# ============================================================
expression_path <- file.path(base_dir, paste0(p_id, "_", cell_clean, "_", assay, ".txt"))
bed_source      <- file.path(pheno_dir, paste0(cell_clean, "_", assay, ".bed"))

if (!file.exists(expression_path)) {
  header_line <- system(paste("head -n 1", shQuote(bed_source)), intern = TRUE)
  writeLines(header_line, expression_path)
  awk_cmd     <- paste0("awk -v gene=", shQuote(p_id), " '$4 == gene' ", shQuote(bed_source), " >> ", shQuote(expression_path))
  system(awk_cmd)
}

expr_dt <- fread(expression_path)
expr_s  <- as.numeric(expr_dt[1, 5:ncol(expr_dt), with=FALSE])
names(expr_s) <- colnames(expr_dt)[5:ncol(expr_dt)]

geno_raw <- fread(geno_raw_file)
var_res  <- extract_variant_series(variant_df = geno_raw, variant_id = v_id)
geno_s   <- var_res$series

covariates_path <- file.path(covar_dir, paste0(cell_clean, "_", assay, "_covariates.txt"))
covar_dt <- fread(covariates_path)
t_covar  <- as.data.frame(t(covar_dt[, -1, with=FALSE]))
colnames(t_covar) <- covar_dt[[1]]
t_covar$IID <- rownames(t_covar)
cov_raw <- t_covar

interaction_file <- find_interaction_file(dt = target, interaction_dir = base_dir)
interaction_path <- interaction_file$interaction_path

int_raw        <- read_table(interaction_path)
interaction_df <- extract_interaction_df(int_raw, names(expr_s))
cov_df         <- extract_covariates_df(cov_raw, names(expr_s))

# ============================================================
# SECTION 4: Sample Alignment & Residualization Framework
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

aligned        <- align_on_samples(expr_s, geno_s, interaction_df, cov_df)
expr_s         <- aligned$expr_s
geno_s         <- aligned$geno_s
interaction_df <- aligned$interaction_df
cov_df         <- aligned$cov_df

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

# ============================================================
# SECTION 5: Exact Original Stats & Plot Framework Functions
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

build_residualized_regression_dataframe <- function(expr_s, geno_s, interaction_df, cov_df = NULL) {
  phenotype_mat <- matrix(as.numeric(expr_s), nrow = 1)
  genotype_mat <- matrix(as.numeric(geno_s), nrow = 1)
  interaction_mat <- as.matrix(interaction_df)
  
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
    residualized_df, geno_s, env_name, gene_name, out_path,
    g_pval = NA_real_, e_pval = NA_real_, gxe_pval = NA_real_,
    conf_level = 0.95, phenotype_type, variant_id = v_id
) {
  y   <- residualized_df$phenotype_centered_residualized
  env <- residualized_df[[paste0(env_name, "_centered_residualized")]]
  g   <- as.numeric(geno_s)
  
  fit     <- fit_plot_model(y, g, env)
  beta    <- fit$beta
  XtX_inv <- fit$XtX_inv
  sigma2  <- fit$sigma2
  g_cls   <- pmax(0, pmin(2, round(g)))
  
  # ============================================================
  # DYNAMIC GENOTYPE LABEL PARSING
  # ============================================================
  parts <- unlist(strsplit(variant_id, "_"))
  if (length(parts) >= 4) {
    ref_allele <- parts[3] # e.g., "C" or "G"
    alt_allele <- parts[4] # e.g., "G" or "A"
    
    # PLINK --recode A counts the alternative allele (0, 1, 2 copies of ALT)
    geno_labels <- c(
      "0" = paste0(ref_allele, "/", ref_allele),
      "1" = paste0(ref_allele, "/", alt_allele),
      "2" = paste0(alt_allele, "/", alt_allele)
    )
  } else {
    # Fallback to numbers if the string format is unexpected
    geno_labels <- c("0" = "0", "1" = "1", "2" = "2")
  }
  # ============================================================
  
  plot_df <- data.frame(
    env      = as.numeric(env),
    y        = as.numeric(y),
    genotype = factor(g_cls, levels = c(0, 1, 2))
  )
  
  env_grid <- seq(min(env, na.rm = TRUE), max(env, na.rm = TRUE), length.out = 120)
  alpha    <- 1 - conf_level
  tcrit    <- qt(1 - alpha / 2, df = max(1, length(y) - 4))
  
  line_df <- do.call(rbind, lapply(0:2, function(gv) {
    xmat   <- cbind(1, gv, env_grid, gv * env_grid)
    yhat   <- as.numeric(xmat %*% beta)
    se_fit <- sqrt(rowSums((xmat %*% XtX_inv) * xmat) * sigma2)
    data.frame(
      env      = env_grid,
      yhat     = yhat,
      lower    = yhat - tcrit * se_fit,
      upper    = yhat + tcrit * se_fit,
      genotype = factor(gv, levels = c(0, 1, 2))
    )
  }))
  
  if (phenotype_type == "Expression"){
    target_id <- gene_name
    my_symbol <- get_symbol(target_id, gene_map)
    cat("The symbol for", target_id, "is", my_symbol, "\n")
    ome       <- "Expression"
  } else {
    my_symbol <- gene_name
    ome       <- "Methylation"
  }
  
  mod_type  <- ifelse(grepl("Metabolite", env_name, ignore.case = TRUE), "Metabolic",
                      ifelse(grepl("Proteome", env_name, ignore.case = TRUE), "Proteomic", "Unknown"))
  color_raw <- sub("^ME([^_]+)_.*", "\\1", env_name)
  color_cap <- paste0(toupper(substr(color_raw, 1, 1)), substring(color_raw, 2))
  env_name_better <- paste("Residualized", mod_type, "Module", color_cap)
  
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
      labels = geno_labels, # Dynamically assigned labels
      name = "Genotype"
    ) +
    ggplot2::scale_fill_manual(
      values = c("0" = "#1f77b4", "1" = "#ff7f0e", "2" = "#2ca02c"),
      breaks = c("0", "1", "2"),
      labels = geno_labels, # Dynamically assigned labels
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
      x     = min(env, na.rm = TRUE),
      y     = max(y, na.rm = TRUE),
      label = pval_label,
      hjust = 0,
      vjust = 1,
      size  = 3.3
    )
  }
  return(p)
}

# ============================================================
# SECTION 6: Processing & Execution For Specific Target Row
# ============================================================
stats_df        <- compute_interaction_stats(expr_s, geno_s, interaction_df, cov_df, v_id, p_id)
residualized_df <- build_residualized_regression_dataframe(expr_s, geno_s, interaction_df, cov_df)

plot_environment <- term_id
if (plot_environment == "g") {
  env_name <- str_replace_all(names(residualized_df)[4], "_centered_residualized", "")
} else {
  resolve_environment_name <- function(requested, columns) {
    cols <- as.character(columns)
    if (requested %in% cols) return(requested)
    req   <- tolower(requested)
    exact <- cols[tolower(cols) == req]
    if (length(exact) == 1) return(exact)
    partial <- cols[grepl(req, tolower(cols), fixed = TRUE)]
    if (length(partial) == 1) return(partial)
    stop(sprintf("Could not uniquely resolve environment: %s", requested))
  }
  env_name <- resolve_environment_name(plot_environment, colnames(interaction_df))
}

pval_col         <- paste0("pval_g-", env_name)
interaction_pval <- if (pval_col %in% colnames(stats_df)) as.numeric(stats_df[[pval_col]][1]) else NA_real_
g_pval           <- if ("pval_g" %in% colnames(stats_df)) as.numeric(stats_df[["pval_g"]][1]) else NA_real_
e_pval_col       <- paste0("pval_", env_name)
e_pval           <- if (e_pval_col %in% colnames(stats_df)) as.numeric(stats_df[[e_pval_col]][1]) else NA_real_
gxe_pval         <- interaction_pval

dir_output     <- "/labs/smontgom/grps_smontgom/hkrupkin/GxE/plots/single_variant_examples"
phenotype_type <- ifelse(grepl("ENSG", p_id), "Expression", "Methylation")
path_output_file <- paste0(dir_output, "/", cell_clean, "_", env_name, "_", target$effect_type, "_", phenotype_type, "_", v_id, ".png")

interaction_plot <- make_interaction_plot(
  residualized_df = residualized_df,
  geno_s          = geno_s,
  env_name        = env_name,
  gene_name       = p_id,
  out_path        = path_output_file,
  g_pval          = g_pval,
  e_pval          = e_pval,
  gxe_pval        = gxe_pval,
  conf_level      = 0.95,
  phenotype_type  = phenotype_type
)

# Render and Save Plot Exactly to Specified Location
print(interaction_plot)
ggsave(plot   = interaction_plot,
       file   = path_output_file,
       width  = 5,
       height = 5,
       dpi    = 900)

message("Success! Plot matched and exported exactly to: ", path_output_file)
