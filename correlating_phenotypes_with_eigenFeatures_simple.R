###correlating_phenotypes_with_eigenFeatures_simple###
#/oak/stanford/groups/smontgom/hkrupkin/GxE/Code/correlating_phenotypes_with_eigenFeatures_simple
#author: Haim Krupkin
#DATE: 01/26/2026

#Description: This script is meant to examine the correlaiton between phenotypes and eigenFeatures in GxE analysis#

###packages#####
.libPaths(c("/scg/apps/software/r/4.3.3/lib", .libPaths()))

library(dplyr)
library(broom)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(stringr)
library(glmnet)
library(purrr)
library(stringr)
library(tibble)
library(tidyverse)
###analysis####
#loading phenotype data
phenotype_data<-read.csv("/oak/stanford/groups/smontgom/mahuynh/metabolites/MESA/data/phenotype/phenotype_MESA.csv")
phenotype_data_exam1<-read.csv("/oak/stanford/groups/smontgom/mahuynh/metabolites/MESA/data/phenotype/phenotype_exam_1.csv")

##subject_id_to_sample_id
ID_Table<-phenotype_data_exam1%>%
  dplyr::select(ID,SUBJECT_ID)

phenotype_data_wID<-merge(ID_Table,
                          phenotype_data,
                          by.x="SUBJECT_ID",
                          by.y="SUBJECT_ID")
###
no_unique_info <- names(phenotype_data_wID)[
  sapply(phenotype_data_wID, function(x) {
    all(is.na(x)) || length(unique(x[!is.na(x)])) <= 1
  })
]
phenotype_data_wID <- phenotype_data_wID[, !names(phenotype_data_wID) %in% no_unique_info]



#loading eigenFeature Data
eigen_value_protein<-read.table("/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/TensorQTL_related/Proteins_18modules_soft_thres_7_deep_split_3_me_0.25_min_size_20_4_removed_pcs.txt",
                                header=T)
colnames(eigen_value_protein)<-paste0("Protein_",colnames(eigen_value_protein))
eigen_value_protein$sample<-rownames(eigen_value_protein)
eigen_value_metabolite<-read.table("/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/TensorQTL_related/52modules_soft_thres_8_deep_split_2_me_0.15_min_size_20.txt",
                                   header=T)
colnames(eigen_value_metabolite)<-paste0("Metabolite",colnames(eigen_value_metabolite))
eigen_value_metabolite$sample<-rownames(eigen_value_metabolite)

all_eigen_values<-merge(eigen_value_metabolite,
                        eigen_value_protein,
                        by.x="sample",
                        by.y="sample",
                        all=T)
####
lall_eigen_values<-all_eigen_values%>%
  pivot_longer(cols=colnames(all_eigen_values)[2:length(colnames(all_eigen_values))],
               names_to="Eigen_Feature",
               values_to="Value")
age_cols <- grep("^age_at_.*_\\d+$", names(phenotype_data_wID), value = TRUE)

## dataframe of just those columns (optionally keep IDs too)
age_df <- phenotype_data_wID[, age_cols, drop = FALSE]
###
id_pheno <- "sample"  # in all_eigen_values
id_pred  <- "ID"      # in phenotype_data_wID
global_covariates <- c("annotated_sex_1", "race_us_1", "geographic_site_1")
age_prefix <- "age_at_"
dat <- phenotype_data_wID %>%
  inner_join(all_eigen_values, by = setNames(id_pheno, id_pred))



exclude_cols <- c(id_pred, "SUBJECT_ID", "unique_subject_key", "topmed_abbreviation")
cn <- setdiff(names(dat), exclude_cols)
age_cols <- cn[str_starts(cn, fixed(age_prefix))]
phenotype_cols <- setdiff(names(all_eigen_values), id_pheno)

pred_candidates <- setdiff(cn, c(age_cols, phenotype_cols))

pairs <- tibble(
  predictor = pred_candidates,
  age_col   = paste0(age_prefix, pred_candidates)
) %>%
  filter(age_col %in% names(dat))
####
pairs$global_covariates<-paste(global_covariates,
                               collapse = " + ")
pairs$per_row_formula<-formulas<-paste0(pairs$predictor," + ",pairs$age_col,
                                                        sep = " + ")
pairs$per_model_formula<-paste0(pairs$per_row_formula,pairs$global_covariates)
list_of_lists_of_formulas<-c()
for (phenotype in phenotype_cols){
  #phenotype<-phenotype_cols[[1]]
  print(phenotype)
  list_of_formulas<-paste0(phenotype,"~",unique(pairs$per_model_formula))
  list_of_formulas
  list_of_lists_of_formulas<-c(list_of_lists_of_formulas,list_of_formulas)

}
#####
head(list_of_lists_of_formulas)
head(dat)
dat_filtered<-dat%>%
  dplyr::filter(!is.na(annotated_sex_1) | !is.na(race_us_1) | !is.na(geographic_site_1))
#

fit_one <- function(form, data) {
  form_chr <- paste(deparse(form), collapse = " ")  # <- key fix
  
  m <- lm(form, data = data, na.action = na.omit)
  
  g <- broom::glance(m)
  df2 <- m$df.residual
  model_p <- pf(g$statistic, g$df, df2, lower.tail = FALSE)
  
  td <- broom::tidy(m) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      partial_r2 = (statistic^2) / (statistic^2 + df2),
      formula = form_chr
    ) %>%
    select(formula, term, estimate, std.error, statistic, p.value, partial_r2)
  
  model_row <- tibble(
    formula = form_chr,
    n = g$nobs,
    r2 = g$r.squared,
    adj_r2 = g$adj.r.squared,
    f = g$statistic,
    model_p = model_p,
    aic = g$AIC,
    bic = g$BIC
  )
  
  list(model = model_row, terms = td)
}

forms <- lapply(list_of_lists_of_formulas, as.formula)
forms
res <- map(forms, fit_one, data = dat)
###scaled_version###
fit_one <- function(form, data) {
  form_chr <- paste(deparse(form), collapse = " ")
  
  # identify variables used in this formula
  vars <- all.vars(form)
  
  # work on a copy; standardize continuous predictors only
  d <- data
  pred_vars <- setdiff(vars, as.character(form[[2]]))  # RHS vars only
  
  cont_preds <- pred_vars[
    vapply(d[pred_vars], function(x) is.numeric(x) && length(unique(x[!is.na(x)])) > 2, logical(1))
  ]
  
  if (length(cont_preds) > 0) {
    d[cont_preds] <- lapply(d[cont_preds], function(x) as.numeric(scale(x)))
  }
  
  m <- lm(form, data = d, na.action = na.omit)
  
  g <- broom::glance(m)
  df2 <- m$df.residual
  model_p <- pf(g$statistic, g$df, df2, lower.tail = FALSE)
  
  td <- broom::tidy(m) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      partial_r2 = (statistic^2) / (statistic^2 + df2),
      formula = form_chr
    ) %>%
    select(formula, term, estimate, std.error, statistic, p.value, partial_r2)
  
  model_row <- tibble(
    formula = form_chr,
    n = g$nobs,
    r2 = g$r.squared,
    adj_r2 = g$adj.r.squared,
    f = g$statistic,
    model_p = model_p,
    aic = g$AIC,
    bic = g$BIC
  )
  
  list(model = model_row, terms = td)
}

forms <- lapply(list_of_lists_of_formulas, as.formula)
res <- purrr::map(forms, fit_one, data = dat)


###



model_table <- bind_rows(map(res, "model"))
term_table  <- bind_rows(map(res, "terms"))
###



write.table(term_table,
            "/oak/stanford/groups/smontgom/hkrupkin/GxE/Data/Outputs/single_lm_models_term_results.tsv")
write.table(model_table,
            "/oak/stanford/groups/smontgom/hkrupkin/GxE/Data/Outputs/single_lm_models_per_model_results.tsv")

#term_table<-fread("/oak/stanford/groups/smontgom/hkrupkin/GxE/Data/Outputs/single_lm_models_term_results.tsv")
model_table<-fread("/oak/stanford/groups/smontgom/hkrupkin/GxE/Data/Outputs/single_lm_models_per_model_results.tsv")
#######
term_table<-read.table("/oak/stanford/groups/smontgom/hkrupkin/GxE/Data/Outputs/single_lm_models_term_results.tsv")
term_table$formula[1:10]

x <- term_table$formula

formula_cols <- tibble(formula = x) %>%
  mutate(formula = str_squish(formula)) %>%          # remove weird spacing/newlines
  separate(formula, into = c("predicted", "rhs"), sep = "\\s*~\\s*", remove = FALSE) %>%
  mutate(vars = str_split(rhs, "\\s*\\+\\s*")) %>%
  mutate(n_vars = map_int(vars, length)) %>%
  unnest_wider(vars, names_sep = "", names_repair = "minimal") %>%
  dplyr::rename(
    variable1 = vars1,
    variable2 = vars2,
    variable3 = vars3,
    variable4 = vars4,
    variable5 = vars5
  ) %>%
  dplyr::select(predicted, variable1:variable5, formula, rhs, n_vars)

# quick check (should be all 5 if truly fixed)
term_table$formula<-str_replace_all(term_table$formula," ","")
formula_cols$formula<-str_replace_all(formula_cols$formula," ","")

vterm_table<-merge(term_table,
                   formula_cols%>%distinct(),
      by.x="formula",
      by.y="formula")
nrow(vterm_table)
only_key_variables<-vterm_table%>%
  dplyr::filter(!(term %in% c("annotated_sex_1male","race_us_1Other","race_us_1White","geographic_site_1MESA_JHU","geographic_site_1MESA_UMN","geographic_site_1MESA_WFU")) & !(str_detect(term,"age_at_")))
only_key_variables

###plotting our results####
only_key_variables<-only_key_variables%>%
  dplyr::mutate(term=str_replace_all(term,"_1",""))

# Will display correlations and their p-values
df <- only_key_variables %>%
  mutate(predicted = as.character(predicted),
         term      = as.character(term))

# make wide matrices
unique(df$term)

df$term_short <- dplyr::case_when(
  df$term == "antihypertensive_meds"     ~ "Anti-HTN meds",
  df$term == "basophil_ncnc_bld"         ~ "Basophils",
  df$term == "bmi_baseline"              ~ "BMI",
  df$term == "bp_diastolic"              ~ "DBP",
  df$term == "bp_systolic"               ~ "SBP",
  df$term == "cac_score"                 ~ "CAC score",
  df$term == "cac_volume"                ~ "CAC volume",
  df$term == "carotid_plaque"            ~ "Carotid plaque",
  df$term == "carotid_stenosis"          ~ "Carotid stenosis",
  df$term == "cd40"                      ~ "CD40",
  df$term == "cimt"                      ~ "CIMT",
  df$term == "cimt_2"                    ~ "CIMT (2)",
  df$term == "crp"                       ~ "CRP",
  df$term == "current_smoker_baseline"   ~ "Current smoker",
  df$term == "eosinophil_ncnc_bld"       ~ "Eosinophils",
  df$term == "eselectin"                 ~ "E-selectin",
  df$term == "ever_smoker_baseline"      ~ "Ever smoker",
  df$term == "fasting_lipids"            ~ "Fasting lipids",
  df$term == "hdl"                       ~ "HDL-C",
  df$term == "height_baseline"           ~ "Height",
  df$term == "hematocrit_vfr_bld"        ~ "Hematocrit",
  df$term == "hemoglobin_mcnc_bld"       ~ "Hemoglobin",
  df$term == "icam1"                     ~ "ICAM-1",
  df$term == "il10"                      ~ "IL-10",
  df$term == "il6"                       ~ "IL-6",
  df$term == "ldl"                       ~ "LDL-C",
  df$term == "lipid_lowering_medication" ~ "Lipid meds",
  df$term == "lppla2_act"                ~ "Lp-PLA2 act",
  df$term == "lppla2_mass"               ~ "Lp-PLA2 mass",
  df$term == "lymphocyte_ncnc_bld"       ~ "Lymphocytes",
  df$term == "mch_entmass_rbc"           ~ "MCH",
  df$term == "mchc_mcnc_rbc"             ~ "MCHC",
  df$term == "mcv_entvol_rbc"            ~ "MCV",
  df$term == "mmp9"                      ~ "MMP-9",
  df$term == "monocyte_ncnc_bld"         ~ "Monocytes",
  df$term == "neutrophil_ncnc_bld"       ~ "Neutrophils",
  df$term == "platelet_ncnc_bld"         ~ "Platelets",
  df$term == "rbc_ncnc_bld"              ~ "RBC",
  df$term == "sleep_duration"            ~ "Sleep duration",
  df$term == "tnfa"                      ~ "TNF-a",
  df$term == "tnfa_r1"                   ~ "TNFR1",
  df$term == "total_cholesterol"         ~ "Total chol",
  df$term == "triglycerides"             ~ "Triglycerides",
  df$term == "wbc_ncnc_bld"              ~ "WBC",
  df$term == "weight_baseline"           ~ "Weight",
  TRUE                                   ~ df$term
)
r2_wide <- df %>%
  dplyr::select(predicted, term_short, partial_r2) %>%
  pivot_wider(names_from = term_short, values_from = partial_r2)

p_wide <- df %>%
  dplyr::select(predicted, term_short, p.value) %>%
  pivot_wider(names_from = term_short, values_from = p.value)

# convert to matrices with matching dimnames
r2_wide <- df %>%
  dplyr::select(predicted, term_short, partial_r2) %>%
  pivot_wider(names_from = term_short, values_from = partial_r2)

r2_mat <- as.matrix(r2_wide[,-1])
rownames(r2_mat) <- r2_wide$predicted

# r2_mat assumed rows=predicted, cols=term, values in [0,1]

# logistic color transform (keeps ordering but spreads mid-range)
logistic01 <- function(x, mid = 0.10, k = 25) {
  1 / (1 + exp(-k * (x - mid)))
}

r2_colmat <- logistic01(r2_mat, mid = 0.10, k = 25)  # tune mid/k if you want

# show text only if r2 >= 0.10
labels_mat <- matrix("", nrow = nrow(r2_mat), ncol = ncol(r2_mat),
                     dimnames = dimnames(r2_mat))
idx <- !is.na(r2_mat) & r2_mat >= 0.01
labels_mat[idx] <- sprintf("%.2f", r2_mat[idx])

# color palette + breaks for the transformed values (0..1)
cols <- colorRampPalette(c("blue", "white", "red"))(100)
bk   <- seq(0, 1, length.out = length(cols) + 1)

colnames(r2_colmat)

pdf("/oak/stanford/groups/smontgom/hkrupkin/GxE/plots/Module-Trait Relationships_all_simple_pheatmap.png",
    width = 12, height = 10)

pheatmap(
  mat = r2_colmat[
    apply(r2_colmat, 1, max, na.rm = TRUE) > 0.1,
    apply(r2_colmat, 2, max, na.rm = TRUE) > 0.1
  ],
  color = cols,
  breaks = bk,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize = 6,
  fontsize_row = 6,
  fontsize_col = 6,
  display_numbers = labels_mat[
    apply(r2_colmat, 1, max, na.rm = TRUE) > 0.1,
    apply(r2_colmat, 2, max, na.rm = TRUE) > 0.1
  ],
  number_color = "black",
  border_color = NA,
  main = "Module-trait relationships (R^2)"
)

dev.off()


###for figure making####
saveRDS(r2_colmat,"/labs/smontgom/grps_smontgom/hkrupkin/GxE/plots/plots_for_figures/heatmap_data_for_clinical_traits_subsample.rds")
saveRDS(cols,"/labs/smontgom/grps_smontgom/hkrupkin/GxE/plots/plots_for_figures/colros_for_heatmap_for_clinical_traits_subsample.rds")
saveRDS(bk,"/labs/smontgom/grps_smontgom/hkrupkin/GxE/plots/plots_for_figures/breaks_for_heatmap_data_for_clinical_traits_subsample.rds")
saveRDS(labels_mat,"/labs/smontgom/grps_smontgom/hkrupkin/GxE/plots/plots_for_figures/labels_mat_for_heatmap_data_for_clinical_traits_subsample.rds")

color = cols
breaks = bk
# 1. Capture the clustering from the original data
new_rownames <- gsub("ME", "_", gsub("_ME", "_", rownames(r2_colmat)))
rownames(r2_colmat) <- new_rownames
rownames(labels_mat) <- new_rownames

# 2. Define your desired subset with the EXACT names from your matrix
rows_to_show <- c("Protein_purple", "Protein_yellow", "Metabolite_tan",
                  "Metabolite_blue", "Metabolite_darkmagenta",
                  "Metabolite_brown4", "Metabolite_turquoise",
                  "Metabolite_white", "Metabolite_violet",
                  "Metabolite_greenyellow", "Metabolite_darkorange",
                  "Metabolite_darkorange2", "Metabolite_skyblue",
                  "Metabolite_cyan", "Metabolite_lightgreen",
                  "Metabolite_salmon", "Metabolite_grey60", # FIX: Changed from "Metabolite_grey"
                  "Protein_blue", "Protein_lightcyan",
                  "Protein_red", "Protein_cyan",
                  "Metabolite_darkgrey")

cols_to_show <- c("LDL-C", "Total chol", "HDL-C", "Triglycerides", "TNFR1", # FIX: Corrected "Total Chol" and "HDL"
                  "BMI", "Weight", "MMP-9", "IL-6", "CRP") # FIX: Corrected "IL6"

# 3. Capture clustering from the full (now renamed) data
##############################################################################################################
library(dendextend)

# Get all row/col labels from the full dendrogram
all_row_labels <- labels(as.dendrogram(p_full$tree_row))
all_col_labels <- labels(as.dendrogram(p_full$tree_col))

# Determine which labels to REMOVE (prune)
rows_to_remove <- setdiff(all_row_labels, rows_to_show)
cols_to_remove <- setdiff(all_col_labels, cols_to_show)

# Prune the dendrograms
row_dend <- as.dendrogram(p_full$tree_row)
for (r in rows_to_remove) {
  row_dend <- prune(row_dend, r)
}

col_dend <- as.dendrogram(p_full$tree_col)
for (c in cols_to_remove) {
  col_dend <- prune(col_dend, c)
}

# Now use in pheatmap




##############################################################################################################


# 4. Draw the new heatmap with the subset
png("/labs/smontgom/grps_smontgom/hkrupkin/GxE/plots/plots_for_figures/Module_Trait_Relationships_subsetV2.png",
    width = 5, height = 4,unit="in",
    res=600)
pheatmap(
  mat = r2_colmat[rows_to_show, cols_to_show],
  color = cols,
  breaks = bk,
  #cluster_rows = as.hclust(row_dend),
  #cluster_cols = as.hclust(col_dend),
  fontsize = 6,
  fontsize_row = 6,
  fontsize_col = 6,
  display_numbers = labels_mat[rows_to_show, cols_to_show],
  number_color = "black",
  border_color = NA,
  main = "Module-trait relationships (R^2) - Subset"
)
dev.off()
#######################


#######################



###actual paper quality figures####
# row annotation: blue = protein, gold = metabolite (based on rowname prefix)
row_type <- dplyr::case_when(
  grepl("^Protein",    rownames(r2_colmat), ignore.case = TRUE) ~ "protein",
  grepl("^Metabolite", rownames(r2_colmat), ignore.case = TRUE) ~ "metabolite",
  TRUE ~ "other"
)

ann_row <- data.frame(EigenFeature = row_type)
rownames(ann_row) <- rownames(r2_colmat)

ann_colors <- list(EigenFeature = c(protein = "blue", metabolite = "gold"))
ann_colors$type
pheatmap(
  mat = r2_colmat,
  color = cols,
  breaks = bk,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize = 6,
  fontsize_row = 6,
  fontsize_col = 6,
  angle_col = 45,
  annotation_row = ann_row,
  annotation_colors = ann_colors,
  display_numbers = labels_mat,
  number_color = "black",
  border_color = NA,
  main = "Module-trait relationships (R^2)"
)

###
r2_flip     <- r2_colmat[rev(rownames(r2_colmat)), rev(colnames(r2_colmat)), drop = FALSE]
labels_flip <- labels_mat[rev(rownames(labels_mat)), rev(colnames(labels_mat)), drop = FALSE]

# rebuild row annotation to match flipped rows
row_type <- dplyr::case_when(
  grepl("^Protein",    rownames(r2_flip), ignore.case = TRUE) ~ "protein",
  grepl("^Metabolite", rownames(r2_flip), ignore.case = TRUE) ~ "metabolite",
  TRUE ~ "other"
)
ann_row <- data.frame(type = row_type)
rownames(ann_row) <- rownames(r2_flip)
pdf("/oak/stanford/groups/smontgom/hkrupkin/GxE/plots/Module-Trait Relationships_all_simple_pheatmap.pdf",
    width = 12, height = 10)

pheatmap::pheatmap(
  mat = r2_flip,
  color = cols,
  breaks = bk,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize = 6,
  fontsize_row = 6,
  fontsize_col = 6,
  angle_col = 45,
  
  annotation_row = ann_row,
  annotation_colors = list(type = c(protein="blue", metabolite="gold")),
  display_numbers = labels_flip,
  number_color = "black",
  border_color = NA,
  main = "Module-trait relationships (R^2)"
)
dev.off()
###pheatmaping estimates###
# convert to matrices with matching dimnames
e_wide <- df %>%
  mutate(estimate = if_else(p.value > 0.05, 0, estimate)) %>%   # non-sig -> 0
  select(predicted, term_short, estimate) %>%
  pivot_wider(names_from = term_short, values_from = estimate, values_fill = 0)

# drop term_short columns that are all 0
keep_cols <- names(e_wide)[-1][colSums(abs(as.matrix(e_wide[,-1])) > 0, na.rm = TRUE) > 0]
signed_logistic <- function(x, k = 10) {
  sign(x) * (1 / (1 + exp(-k * abs(x))) - 0.5) * 2
}

## =========================
## Flip matrix + labels
## =========================
e_mat_flip <- e_mat[
  rev(rownames(e_mat)),
  rev(colnames(e_mat)),
  drop = FALSE
]

labels_flip <- labels_mat[
  rev(rownames(labels_mat)),
  rev(colnames(labels_mat)),
  drop = FALSE
]

## =========================
## Apply logistic transform
## =========================
e_mat_log <- signed_logistic(e_mat_flip, k = 50)

## =========================
## Build row annotation
## =========================
row_type <- dplyr::case_when(
  grepl("^Protein", rownames(e_mat_log), ignore.case = TRUE) ~ "protein",
  grepl("^Metabolite", rownames(e_mat_log), ignore.case = TRUE) ~ "metabolite",
  TRUE ~ "other"
)

ann_row <- data.frame(type = row_type)
rownames(ann_row) <- rownames(e_mat_log)

## =========================
## Color scale (diverging)
## =========================
lim <- max(abs(e_mat_log), na.rm = TRUE)

bk <- seq(-lim, lim, length.out = 101)

cols <- colorRampPalette(c(
  "#2166AC",  # deep blue
  "white",
  "#B2182B"   # deep red
))(length(bk) - 1)

## =========================
## Plot
## =========================
pdf(
  "/oak/stanford/groups/smontgom/hkrupkin/GxE/plots/Module-Trait_Relationships_eMat_logistic.pdf",
  width = 12,
  height = 10
)

pheatmap(
  mat = e_mat,
  color = cols,
  breaks = bk,
  
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  
  fontsize = 6,
  fontsize_row = 6,
  fontsize_col = 6,
  angle_col = "90",
  
  annotation_row = ann_row,
  annotation_colors = list(
    type = c(
      protein = "blue",
      metabolite = "gold"
    )
  ),
  
  number_color = "black",
  
  border_color = NA,
  
  main = "Module–trait relationships (signed logistic scale; 0 = white)"
)

dev.off()

###additional attempt###
library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)

# Step 1: Rename terms to short names
only_key_variables$term_short <- dplyr::case_when(
  only_key_variables$term == "antihypertensive_meds"     ~ "Anti-HTN meds",
  only_key_variables$term == "basophil_ncnc_bld"         ~ "Basophils",
  only_key_variables$term == "bmi_baseline"              ~ "BMI",
  only_key_variables$term == "bp_diastolic"              ~ "DBP",
  only_key_variables$term == "bp_systolic"               ~ "SBP",
  only_key_variables$term == "cac_score"                 ~ "CAC score",
  only_key_variables$term == "cac_volume"                ~ "CAC volume",
  only_key_variables$term == "carotid_plaque"            ~ "Carotid plaque",
  only_key_variables$term == "carotid_stenosis"          ~ "Carotid stenosis",
  only_key_variables$term == "cd40"                      ~ "CD40",
  only_key_variables$term == "cimt"                      ~ "CIMT",
  only_key_variables$term == "cimt_2"                    ~ "CIMT (2)",
  only_key_variables$term == "crp"                       ~ "CRP",
  only_key_variables$term == "current_smoker_baseline"   ~ "Current smoker",
  only_key_variables$term == "eosinophil_ncnc_bld"       ~ "Eosinophils",
  only_key_variables$term == "eselectin"                 ~ "E-selectin",
  only_key_variables$term == "ever_smoker_baseline"      ~ "Ever smoker",
  only_key_variables$term == "fasting_lipids"            ~ "Fasting lipids",
  only_key_variables$term == "hdl"                       ~ "HDL-C",
  only_key_variables$term == "height_baseline"           ~ "Height",
  only_key_variables$term == "hematocrit_vfr_bld"        ~ "Hematocrit",
  only_key_variables$term == "hemoglobin_mcnc_bld"       ~ "Hemoglobin",
  only_key_variables$term == "icam1"                     ~ "ICAM-1",
  only_key_variables$term == "il10"                      ~ "IL-10",
  only_key_variables$term == "il6"                       ~ "IL-6",
  only_key_variables$term == "ldl"                       ~ "LDL-C",
  only_key_variables$term == "lipid_lowering_medication" ~ "Lipid meds",
  only_key_variables$term == "lppla2_act"                ~ "Lp-PLA2 act",
  only_key_variables$term == "lppla2_mass"               ~ "Lp-PLA2 mass",
  only_key_variables$term == "lymphocyte_ncnc_bld"       ~ "Lymphocytes",
  only_key_variables$term == "mch_entmass_rbc"           ~ "MCH",
  only_key_variables$term == "mchc_mcnc_rbc"             ~ "MCHC",
  only_key_variables$term == "mcv_entvol_rbc"            ~ "MCV",
  only_key_variables$term == "mmp9"                      ~ "MMP-9",
  only_key_variables$term == "monocyte_ncnc_bld"         ~ "Monocytes",
  only_key_variables$term == "neutrophil_ncnc_bld"       ~ "Neutrophils",
  only_key_variables$term == "platelet_ncnc_bld"         ~ "Platelets",
  only_key_variables$term == "rbc_ncnc_bld"              ~ "RBC",
  only_key_variables$term == "sleep_duration"            ~ "Sleep duration",
  only_key_variables$term == "tnfa"                      ~ "TNF-a",
  only_key_variables$term == "tnfa_r1"                   ~ "TNFR1",
  only_key_variables$term == "total_cholesterol"         ~ "Total chol",
  only_key_variables$term == "triglycerides"             ~ "Triglycerides",
  only_key_variables$term == "wbc_ncnc_bld"              ~ "WBC",
  only_key_variables$term == "weight_baseline"           ~ "Weight",
  TRUE                                                   ~ only_key_variables$term
)

# Step 2: Clean module names (keep prefix, remove "ME")
only_key_variables$module_clean <- gsub("ME", "_", only_key_variables$predicted)
only_key_variables$module_clean <- gsub("__", "_", only_key_variables$module_clean)

# Step 3: Build full R² matrix (rows = traits, cols = modules)
heatmap_data <- only_key_variables %>%
  dplyr::select(module_clean, term_short, partial_r2) %>%
  distinct()

# Check for duplicates that would cause pivot issues
dupes <- heatmap_data %>%
  group_by(module_clean, term_short) %>%
  filter(n() > 1)
if (nrow(dupes) > 0) {
  cat("WARNING: Found duplicate module-trait pairs. Taking the mean.\n")
  heatmap_data <- heatmap_data %>%
    group_by(module_clean, term_short) %>%
    summarise(partial_r2 = mean(partial_r2, na.rm = TRUE), .groups = "drop")
}

r2_mat <- heatmap_data %>%
  pivot_wider(names_from = module_clean, values_from = partial_r2) %>%
  column_to_rownames("term_short") %>%
  as.matrix()

# Replace NAs with 0
r2_mat[is.na(r2_mat)] <- 0

cat("Full matrix dimensions:", dim(r2_mat), "\n")
cat("Full matrix range:", range(r2_mat), "\n")

# Step 4: Cluster on the FULL matrix
p_full <- pheatmap(
  mat = r2_mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  silent = TRUE
)

# Step 5: Iteratively filter for subset (partial R² > 0.1)
mat_filtered <- r2_mat
prev_dim <- c(0, 0)

while (!identical(dim(mat_filtered), prev_dim)) {
  prev_dim <- dim(mat_filtered)
  
  # Keep rows (traits) where at least one module has R² > 0.1
  rows_keep <- which(apply(mat_filtered, 1, function(x) any(x > 0.1, na.rm = TRUE)))
  mat_filtered <- mat_filtered[rows_keep, , drop = FALSE]
  
  # Keep columns (modules) where at least one trait has R² > 0.1
  cols_keep <- which(apply(mat_filtered, 2, function(x) any(x > 0.1, na.rm = TRUE)))
  mat_filtered <- mat_filtered[, cols_keep, drop = FALSE]
}

cat("Filtered matrix dimensions:", dim(mat_filtered), "\n")
cat("Filtered matrix range:", range(mat_filtered), "\n")

# Verify filtering worked
cat("\nMax per column (all should have at least one > 0.1):\n")
print(sort(apply(mat_filtered, 2, max, na.rm = TRUE)))
cat("\nMax per row (all should have at least one > 0.1):\n")
print(sort(apply(mat_filtered, 1, max, na.rm = TRUE)))

# Step 6: Get the full clustering order, apply to filtered subset
full_row_order <- rownames(r2_mat)[p_full$tree_row$order]
full_col_order <- colnames(r2_mat)[p_full$tree_col$order]

ordered_rows <- full_row_order[full_row_order %in% rownames(mat_filtered)]
ordered_cols <- full_col_order[full_col_order %in% colnames(mat_filtered)]

# Reorder filtered matrix by full clustering order
mat_plot <- mat_filtered[ordered_rows, ordered_cols]

# Step 7: Define color scheme
max_val <- max(mat_plot, na.rm = TRUE)
bk <- seq(0, max_val, length.out = 100)
colors <- colorRampPalette(c("white", "yellow", "orange", "red"))(length(bk) - 1)

# Step 8: Plot
png("/labs/smontgom/grps_smontgom/hkrupkin/GxE/plots/plots_for_figures/Module_Trait_Relationships_filtered_0.1.png",
    width = 7, height = 6, unit = "in", res = 600)

print(pheatmap(
  mat = mat_plot,
  color = colors,
  breaks = bk,
  fontsize = 6,
  fontsize_row = 7,
  fontsize_col = 7,
  display_numbers = TRUE,
  number_format = "%.2f",
  number_color = "black",
  border_color = NA,
  angle_col = 90,
  main = "Module-Trait Relationships Filtered Partial R² > 0.1"
))

dev.off()

