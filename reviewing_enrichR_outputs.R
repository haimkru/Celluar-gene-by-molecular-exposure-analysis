###reviewing_enrichR_outputs####
#/oak/stanford/groups/smontgom/hkrupkin/GxE/Code/reviewing_enrichR_outputs
#auhtor:Haim Krupkin revieiwing maries work
#date:05/12/2026

#description:
#here we analyse the output fromm mamrie's enrichR anallysis
#the purpose of the original analysis was to provide annotations to whole WGCNA MODULES


###packages###
.libPaths(c("/scg/apps/software/r/4.3.3/lib", .libPaths()))

library(data.table)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(stringr)

#install.packages("readxl",lib="/oak/stanford/groups/smontgom/hkrupkin/Code")
#install.packages("tibble",lib="/oak/stanford/groups/smontgom/hkrupkin/Code")
#install.packages("dplyr",lib="/oak/stanford/groups/smontgom/hkrupkin/Code")
#install.packages("utf8", lib="/oak/stanford/groups/smontgom/hkrupkin/Code")
#install.packages("poolr", lib="/oak/stanford/groups/smontgom/hkrupkin/Code")

library(tibble,lib.loc="/oak/stanford/groups/smontgom/hkrupkin/Code")
library(dplyr,lib.loc="/oak/stanford/groups/smontgom/hkrupkin/Code")
#library(utf8,lib.loc="/oak/stanford/groups/smontgom/hkrupkin/Code")
library(readxl,lib.loc="/oak/stanford/groups/smontgom/hkrupkin/Code")
library(poolr,lib.loc="/oak/stanford/groups/smontgom/hkrupkin/Code")
library(stringr)
###analysis####


list.files("/oak/stanford/groups/smontgom/mahuynh/metabolites/MESA/data/associated_features/Exam 1 - proteins/enrichr",
           full.names=T,
           pattern="*.xlsx")

file <- "/oak/stanford/groups/smontgom/mahuynh/metabolites/MESA/data/associated_features/Exam 1 - proteins/enrichr/1yellow_members_enrichr.xlsx"

file <- "/oak/stanford/groups/smontgom/mahuynh/metabolites/MESA/data/associated_features/Exam 1 - proteins/enrichr/1yellow_members_enrichr.xlsx"

# extract color from filename
color <- gsub(".*1(\\w+)_members_enrichr\\.xlsx", "\\1", basename(file))

enrichr_results <- do.call(rbind,
                           lapply(excel_sheets(file), function(sheet) {
                             df <- read_xlsx(file, sheet = sheet)
                             df$sheet_name    <- sheet
                             df$file_name     <- basename(file)
                             df$module_color  <- color
                             df
                           })
)


###big rbinding####
enrichr_dir <- "/oak/stanford/groups/smontgom/mahuynh/metabolites/MESA/data/associated_features/Exam 1 - proteins/enrichr"

files <- list.files(enrichr_dir, full.names = TRUE, pattern = "*.xlsx")

enrichr_all <- do.call(rbind, lapply(files, function(file) {
  color <- gsub(".*1(\\w+)_members_enrichr\\.xlsx", "\\1", basename(file))
  do.call(rbind, lapply(excel_sheets(file), function(sheet) {
    df <- read_xlsx(file, sheet = sheet)
    df$sheet_name   <- sheet
    df$file_name    <- basename(file)
    df$module_color <- color
    df
  }))
}))

write.csv(enrichr_all, 
          "/oak/stanford/groups/smontgom/hkrupkin/GxE/Data/supplementary_tables/all_enrichR_outputs.csv",
          row.names = FALSE)

cat("Done! Dimensions:", dim(enrichr_all), "\n")

###lets do just an example plot for a singular term, maybe yellow cause its associated with CRP and IL6###

plot_enrichR_example_yellow <- enrichr_all %>%
  as.data.frame() %>%
  dplyr::filter(module_color == "yellow") %>%
  mutate(Genes_sorted = sapply(strsplit(Genes, ";"), function(x) paste(sort(trimws(x)), collapse = ";"))) %>%
  group_by(sheet_name, Genes_sorted) %>%
  slice_min(order_by = Adjusted.P.value, n = 1) %>%
  ungroup() %>%
  select(-Genes_sorted) %>%
  group_by(sheet_name) %>%
  slice_max(order_by = -log10(Adjusted.P.value), n = 5) %>%
  ungroup() %>%
  mutate(Term = str_remove(Term, "\\s*\\(GO:\\d+\\)")) %>%
  mutate(Term = str_remove(Term, "\\s*R-HSA-\\d+")) %>%
  mutate(Term = str_wrap(Term, width = 30)) %>%
  mutate(sheet_name = str_wrap(str_replace_all(sheet_name, "_", " "), width = 20)) %>%
  ggplot(aes(x = reorder(Term, -log10(Adjusted.P.value)), 
             y = -log10(Adjusted.P.value), 
             size = Odds.Ratio)) +
  facet_wrap(~sheet_name, scales = "free", nrow = 1) +
  geom_segment(aes(x = reorder(Term, -log10(Adjusted.P.value)), 
                   xend = reorder(Term, -log10(Adjusted.P.value)), 
                   y = 0, 
                   yend = -log10(Adjusted.P.value)),
               size = 0.8, color = "grey50") +
  geom_point(aes(color = -log10(Adjusted.P.value), size = Odds.Ratio)) +
  scale_color_gradient(low = "blue", high = "red") +
  coord_flip() +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 6, lineheight = 0.8),
    strip.text = element_text(size = 7, face = "bold"),
    panel.grid.major.y = element_blank(),
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8)
  ) +
  labs(x = "", y = "-log10(Adjusted P-value)", 
       color = "-log10(Adj. P)", size = "Odds Ratio")

ggsave(plot = plot_enrichR_example_yellow,
       filename = "/oak/stanford/groups/smontgom/hkrupkin/GxE/plots/plots_for_figures/plot_enrichR_example_yellow.png",
       width = 14,
       height = 4,
       dpi=900)


