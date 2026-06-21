###explainabiltiy_of_xgboosted_GxE_predictions###############
#/labs/smontgom/grps_smontgom/hkrupkin/GxE/Code/explainabiltiy_of_xgboosted_GxE_predictions
#author: Haim Krupkin
#date: 05/15/2026

###description
#

###package####
.libPaths(c("/scg/apps/software/r/4.3.3/lib", .libPaths()))

library(xgboost)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(ggpubr)

#
#install.packages("SHAPforxgboost",loc="/labs/smontgom/grps_smontgom/hkrupkin/Code/")
library(SHAPforxgboost,lib.loc="/labs/smontgom/grps_smontgom/hkrupkin/Code/") # Great wrapper library for R visualizations

###analysis#####

# 1. Load your saved workspace
locto_data <- readRDS("/labs/smontgom/grps_smontgom/hkrupkin/GxE/plots/plots_for_figures/locto_xgboost_models_with_data.rds")

# ── loop over all cell types or pick one ──────────────────────────────────
cell_types   <- names(locto_data)          # e.g. c("mono","cd4","cd8",...)
target_cells <- cell_types                 # change to c("mono") for a single cell

output_dir <- "/labs/smontgom/grps_smontgom/hkrupkin/GxE/plots/shap_ggplots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

###############################################################
# helper — build tidy SHAP long table from SHAPforxgboost
###############################################################
build_shap_long <- function(cell, data_list) {
  model  <- data_list[[cell]]$model
  X_test <- data_list[[cell]]$X_test
  
  # SHAPforxgboost::shap.prep returns a tidy data.frame
  shap_long <- shap.prep(xgb_model = model, X_train = X_test)
  shap_long$cell_type <- cell
  shap_long
}

###############################################################
# 2.  Plot A — SHAP Beeswarm summary  (one per cell type)
###############################################################

plot_shap_beeswarm <- function(shap_long, cell, top_n = 20) {
  
  # rank by mean |SHAP|
  top_features <- shap_long %>%
    group_by(variable) %>%
    summarise(mean_abs = mean(abs(value)), .groups = "drop") %>%
    slice_max(mean_abs, n = top_n) %>%
    pull(variable)
  
  df <- shap_long %>%
    filter(variable %in% top_features) %>%
    mutate(variable = fct_reorder(variable, abs(value), .fun = mean))
  
  ggplot(df, aes(x = value, y = variable, colour = rfvalue)) +
    geom_jitter(
      height   = 0.35,
      size     = 0.8,
      alpha    = 0.55
    ) +
    geom_vline(xintercept = 0, linewidth = 0.5, colour = "grey40", linetype = "dashed") +
    shap_palette +
    labs(
      title    = paste0("SHAP Summary — ", toupper(cell)),
      subtitle = paste0("Top ", top_n, " features by mean |SHAP|"),
      x        = "SHAP value  (impact on model output)",
      y        = NULL
    ) +
    theme_shap
}

###############################################################
# 3.  Plot B — Mean |SHAP| importance bar chart
###############################################################
plot_shap_importance <- function(shap_long, cell, top_n = 20) {
  
  df <- shap_long %>%
    group_by(variable) %>%
    summarise(mean_abs_shap = mean(abs(value)), .groups = "drop") %>%
    slice_max(mean_abs_shap, n = top_n) %>%
    mutate(variable = fct_reorder(variable, mean_abs_shap))
  
  ggplot(df, aes(x = mean_abs_shap, y = variable)) +
    geom_col(alpha = 0.85, width = 0.7) +
    geom_text(
      aes(label = round(mean_abs_shap, 3)),
      hjust = -0.1, size = 3.2, colour = "black"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      x        = "Mean  SHAP value",
      y        = NULL
    ) +
    theme_shap +
    theme(legend.position = "none")
}

###############################################################
# 4.  Plot C — SHAP dependence plot for a single feature
###############################################################
plot_shap_dependence <- function(shap_long, cell, feature,
                                 colour_by = NULL) {
  #feature<-"gene_resp_rank"
  df <- shap_long %>% filter(variable == feature)
  
  colour_feature <- feature
  
  df_colour <- shap_long %>%
    filter(variable == colour_feature) %>%
    select(ID = ID, colour_val = rfvalue)       # rfvalue = scaled raw feature
  
  df <- df %>%
    rename(ID = ID) %>%
    left_join(df_colour, by = "ID")
  
  ggplot(df, aes(x = rfvalue, y = value, colour = colour_val)) +
    geom_point(size = 1.2, alpha = 0.6) +
    geom_smooth(method = "loess", colour = "black", linewidth = 0.8, se = TRUE) +
    scale_colour_gradient2(
      low = "#3B82F6", mid = "#F3F4F6", high = "#EF4444",
      midpoint = 0.5, name = colour_feature
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    labs(
      title    = paste0("SHAP Dependence: ", feature, " — ", toupper(cell)),
      x        = paste("Feature value:", feature),
      y        = paste("SHAP value:", feature)
    ) +
    theme_shap
}

###############################################################
# 5.  Multi-cell stacked importance comparison
###############################################################
plot_importance_faceted <- function(all_shap, top_n = 15) {
  
  df <- all_shap %>%
    group_by(cell_type, variable) %>%
    summarise(mean_abs_shap = mean(abs(value)), .groups = "drop") %>%
    group_by(cell_type) %>%
    slice_max(mean_abs_shap, n = top_n) %>%
    ungroup() %>%
    mutate(variable = fct_reorder(variable, mean_abs_shap))
  
  ggplot(df, aes(x = mean_abs_shap, y = variable, fill = cell_type)) +
    geom_col(show.legend = FALSE, width = 0.7, alpha = 0.85) +
    facet_wrap(~cell_type, scales = "free_y") +
    scale_fill_brewer(palette = "Set2") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(
      title    = "Feature Importance across Cell Types",
      subtitle = paste0("Top ", top_n, " features per cell type (mean |SHAP|)"),
      x        = "Mean |SHAP value|",
      y        = NULL
    ) +
    theme_shap +
    theme(strip.text = element_text(face = "bold"))
}

###############################################################
# 6.  Run everything and save
###############################################################
all_shap_list <- list()
shap_palette <- scale_colour_gradient2(
  low      = "#3B82F6",
  mid      = "#F3F4F6",
  high     = "#EF4444",
  midpoint = 0,
  name     = "Feature\nvalue"
)

theme_shap <- theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", size = 15, hjust = 0),
    plot.subtitle    = element_text(colour = "grey50", size = 11, hjust = 0),
    panel.grid.major = element_line(colour = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title.x     = element_text(margin = margin(t = 8)),
    legend.position  = "right"
  )


for (cell in target_cells) {
  #cell<-target_cells[1]
  message("Processing cell type: ", cell)
  
  shap_long <- build_shap_long(cell, locto_data)
  all_shap_list[[cell]] <- shap_long
  
  # ── Plot A: beeswarm
  p_bee <- plot_shap_beeswarm(shap_long, cell, top_n = 20)
  ggsave(plot=p_bee,
         file.path(output_dir, paste0(cell, "_shap_beeswarm.png")),
         width = 8, height = 7
         
  )

  # ── Plot B: importance bar
  p_imp <- plot_shap_importance(shap_long, cell, top_n = 20)
  ggsave(plot=p_imp,
         file.path(output_dir, paste0(cell, "_shap_importance.png")),
         width = 5, height = 4
  )
  
  # ── Plot C: dependence for top feature
  top_feat <- shap_long %>%
    group_by(variable) %>%
    summarise(m = mean(abs(value)), .groups = "drop") %>%
    slice_max(m, n = 1) %>%
    pull(variable)
  
  p_dep <- plot_shap_dependence(shap_long, cell, feature = top_feat)
  ggsave(plot=p_dep,
         file.path(output_dir, paste0(cell, "_shap_dependence_", top_feat, ".png")),
         width = 8, height = 7
         
  )
  
  
  message("  Saved plots for: ", cell)
}

# ── Plot D: multi-cell faceted comparison
all_shap <- bind_rows(all_shap_list)

p_facet <- plot_importance_faceted(all_shap, top_n = 15)
ggsave(
  file.path(output_dir, "all_cells_shap_importance_faceted.pdf"),
  p_facet,
  width  = 4 * ceiling(sqrt(length(target_cells))),
  height = 4 * ceiling(length(target_cells) / ceiling(sqrt(length(target_cells)))),
  device = cairo_pdf
)

message("\nAll SHAP ggplots saved to: ", output_dir)

###saving shap values for all cell types####
all_shap_list_big <- list()

for (cell in target_cells) {
  #cell<-target_cells[1]
  print(cell)
  message("Processing cell type: ", cell)
  
  shap_long_once <- build_shap_long(cell, locto_data)
  all_shap_list_big[[cell]] <- shap_long_once
}

all_shap_big <- bind_rows(all_shap_list_big)

saveRDS(all_shap_big,"/labs/smontgom/grps_smontgom/hkrupkin/GxE/plots/shap_ggplots/all_shap_big.RDS")



df_variable_importance_by_cecll_type <- all_shap_big %>%
  group_by(variable,cell_type) %>%
  summarise(mean_abs_shap = mean(abs(value)), .groups = "drop") %>%
  mutate(variable = fct_reorder(variable, mean_abs_shap))

head(df_variable_importance_by_cecll_type)
medians <- df_variable_importance_by_cecll_type %>%
  group_by(variable) %>%
  summarise(med = median(mean_abs_shap), .groups = "drop")

medians <- df_variable_importance_by_cecll_type %>%
  group_by(variable) %>%
  summarise(med = median(mean_abs_shap), .groups = "drop")

p_imp_by_cell_type <- df_variable_importance_by_cecll_type %>%
  ggplot(aes(x = variable, y = mean_abs_shap)) +
  geom_boxplot(alpha = 0.85, width = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.2, alpha = 0.5, colour = "#3B82F6") +
  geom_text(
    data     = medians,
    aes(x    = variable, y = med, label = round(med, 3)),
    hjust    = -0.2,       # nudge label just right of the median line (after coord_flip)
    size     = 3.2,
    colour   = "black",
    fontface = "bold"
  ) +
  stat_compare_means(
    label.y  = max(df_variable_importance_by_cecll_type$mean_abs_shap) * 0.5,
    size     = 3.5
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +  # room for p-value label
  labs(
    y = "Mean -SHAP value",
    x = NULL
  ) +
  theme_shap +
  theme(legend.position = "none") +
  coord_flip()

p_imp_by_cell_type

ggsave(
  plot = p_imp_by_cell_type,
  file.path(output_dir, paste0( "shap_importance_all_cell_types.png")),
  width = 5, height = 3.5
)


