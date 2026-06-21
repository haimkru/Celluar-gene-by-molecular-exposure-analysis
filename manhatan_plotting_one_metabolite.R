###manhatan_plotting_one_metabolite####
#/oak/stanford/groups/smontgom/hkrupkin/GxE/Code/manhatan_plotting_one_metabolite
#author: haim krupkin
#date: 03/20/2026
#descriptioN: make a manhatan plot for one tested thingy for GxE

##old code did it worng:
#zcat /labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/for_developments/output_bigs/mono_expression_MEgreenyellow_trans_qtl.tsv.gz \
# #| awk 'BEGIN{FS=OFS="\t"} NR==1 || $9 > 0.01' \
#> /labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/for_developments/output_bigs/mono_expression_MEgreenyellow_trans_qtl_filtered_pvalME_gt_0.01.tsv
#> 
#> 
#!/usr/bin/env Rscript
#right code did it right:
#zcat mono_expression_MEgreenyellow_trans_qtl.tsv.gz | awk 'NR==1 || $12 < 0.01' > mono_expression_MEgreenyellow_trans_qtl_filtered_12th_colun_pval_above_0.01.tsv
# Manhattan plot — MEgreenyellow trans-QTL (interaction p-value)
install.packages("qqman",lib = "/oak/stanford/groups/smontgom/hkrupkin/GxE/Code")


suppressPackageStartupMessages({
  library(data.table)
  library(qqman,lib.loc="/oak/stanford/groups/smontgom/hkrupkin/GxE/Code")
})

# ── CONFIG ────────────────────────────────────────────────────────────────────
PVAL_COL     <- "pval"
OUTPUT_FILE  <- "/labs/smontgom/grps_smontgom/hkrupkin/GxE/plots/plots_for_figures/manhatan_single_example_mono_expression_greenyellow_metabolite.png"
# ─────────────────────────────────────────────────────────────────────────────


# ── 1. LOAD & PARSE ───────────────────────────────────────────────────────────
cat("Loading data...\n")
#before 04/27/2026
#dt<-fread("/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/for_developments/output_bigs/mono_expression_MEgreenyellow_trans_qtl_filtered_pvalME_gt_0.01.tsv")
#after 04/27/2026

dt<-fread("/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/for_developments/output_bigs/mono_expression_MEgreenyellow_trans_qtl_filtered_12th_colun_pval_above_0.01.tsv")
#dt <- Expresion_hits_specific
cat(sprintf("  %s rows\n", format(nrow(dt), big.mark = ",")))

dt[, CHR := sub("^chr", "", sub("^(chr[^_]+)_.*", "\\1", variant_id))]
dt[, BP  := as.integer(sub("^chr[^_]+_(\\d+)_.*", "\\1", variant_id))]
dt[, P   := get(PVAL_COL)]
dt[, SNP := variant_id]

dt <- dt[CHR %in% c(as.character(1:22), "X") & !is.na(BP) & !is.na(P) & P > 0]
dt[CHR == "X", CHR := "23"]
dt[, CHR := as.integer(CHR)]
setorder(dt, CHR, BP)

gwas_df <- as.data.frame(dt[, .(SNP, CHR, BP, P)])
cat(sprintf("  After filtering: %s rows\n", format(nrow(gwas_df), big.mark = ",")))


# ── 2. PLOT ───────────────────────────────────────────────────────────────────
cat("Plotting...\n")
png(OUTPUT_FILE, width = 2000, height = 600, res = 200)
par(mar = c(4, 5, 4, 2))

manhattan(
  gwas_df,
  chr            = "CHR",
  bp             = "BP",
  p              = "P",
  snp            = "SNP",
  col            = c("grey30", "grey60"),
  cex            = 0.5,
  genomewideline = -log10(1e-5),
  suggestiveline = FALSE
)

dev.off()
cat(sprintf("Saved → %s\n", OUTPUT_FILE))