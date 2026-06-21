###script_to_plot_specific_colocs_GxE#####
#/oak/stanford/groups/smontgom/hkrupkin/GxE/Code/script_to_plot_specific_colocs_GxE
#author: Haim krupkin
#date: 04/27/2026 (at 11:55pm lols)


#description: this script loads the data of the big run of colocs from coloc_and_viewV3
#it then helps me chose one to manually plot, and then we plot it :)

#notes: none, this is part of the final analyses, So I am pretty excited to do it :)

###packages#####
.libPaths(c("/scg/apps/software/r/4.3.3/lib", .libPaths()))

library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
library(ggpubr)
library(ggplot2)
library(cowplot)
library(stringr)
####analysis
#first we must load the colocs
files <- list.files(
  #  "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/coloc_interaction_with_GWASes",
  "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/coloc_interaction_with_GWASesV2",
  #/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/coloc_interaction_with_GWASesV2
  pattern = "*all_coloc_seed_outputs_V2*",
  full.names = TRUE
)
files<-files[grepl(".rds",files)]
files
lslist_data <- lapply(files, function(f) {
  df <- readRDS(f)
  df$source_file <- f
  df
})

coloc_data <- rbindlist(lslist_data, use.names = TRUE, fill = TRUE)


candidates <- coloc_data[PP.H4 > 0.8 ]

# Sort by the current PP.H4 to see which one was 'closest' to working
candidates <- candidates[order(-PP.H4)]
#we hav efinished loading the intresting colocs
#ya habibi
###
idb_option<-candidates%>%
  dplyr::filter(gwas_name=="imputed_IBD.EUR.Crohns_Disease._active_loci")%>%
  pull(gxe_name)
idb_option<-idb_option[1]
ldl_options<-candidates%>%
  dplyr::filter(gwas_name=="imputed_MAGNETIC_LDL.C._active_loci")%>%
  pull(gxe_name)
ldl_options<-ldl_options[1] 
###
#important to be able to get the correct gene annotation
gtf_path <- "/oak/stanford/groups/smontgom/dnachun/projects/topmed/gencode.v39.annotation.gtf.gz"

# 2. Use a system command to extract the mapping (FAST)
# We only care about lines where the 3rd column is 'gene'
# We then extract gene_id and gene_name from the attributes
cmd <- paste0("zcat ", gtf_path, " | awk '$3 == \"gene\" {print $0}' | sed -E 's/.*gene_id \"([^ \"]+)\".*gene_name \"([^ \"]+)\".*/\\1\\t\\2/'")

# 3. Read into R
gene_map <- fread(cmd, col.names = c("gene_id", "gene_name"), header = FALSE)

# 4. Remove duplicates (GTF often has many entries per ID)
gene_map <- unique(gene_map)

# 5. Fast Lookup function
get_symbol <- function(ensembl_id, map_dt) {
  #gene_r
  #map_dt<-gene_map
  #ensembl_id<-gene_r
  # This works even if your ID has a version (e.g., .1) and the map doesn't, or vice-versa
  # Strip version for matching if necessary
  symbol <- map_dt[gene_id == ensembl_id, gene_name]
  symbol
  if (length(symbol) == 0) {
    # Try matching without version number
    clean_id <- gsub("\\..*", "", ensembl_id)
    symbol <- map_dt[gsub("\\..*", "", gene_id) == clean_id, gene_name]
  }
  
  return(if(length(symbol) > 0) symbol[1] else NA)
}






###now lets define a function to do the ploting for us :3


# Ensure these global paths match your environment
gwas_dir   <- "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/processed_rds_gwas"
eqtl_dir   <- "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/saving_for_coloc_susieV2"
WINDOW_SIZE <- 200000 

get_coloc_plot <- function(hit_row) {
  #hit_row<-hdl_row
  hit_row
  # 1. Parse Metadata from the hit row
  gene_r    <- as.character(hit_row$Gene)
  gwas_name <- as.character(hit_row$gwas_name)
  gxe_path  <- as.character(hit_row$gxe_name) # Uses the path stored during coloc
  lead_pos  <- as.numeric(hit_row$Lead_pos)
  
  # Parse region bounds
  parts   <- strsplit(as.character(hit_row$Region), "_")[[1]]
  chr_r   <- parts[1]
  start_r <- as.numeric(parts[2])
  end_r<-as.numeric(parts[3])
  # 2. Load GWAS Data
  gwas_path <- file.path(gwas_dir, paste0(gwas_name, ".rds"))
  if(!file.exists(gwas_path)) stop("GWAS file not found: ", gwas_path)
  gwas_full <- as.data.table(readRDS(gwas_path))
  
  gwas_sub  <- gwas_full[chromosome == chr_r & position %between% c(start_r, end_r)]
  gwas_sub[pvalue <= 0 | is.na(pvalue), pvalue := 1e-300]
  gwas_sub[, neglog10p := -log10(pvalue)]
  gwas_sub
  # 3. Load GxE Interaction Data
  gxe_path_full<-paste0("/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/saving_for_coloc_susieV2/",gxe_path,".rds")
  gxe_path_full
  if(!file.exists(gxe_path_full)) stop("GxE file not found: ", gxe_path_full)
  gxe_dt   <- as.data.table(readRDS(gxe_path_full))
  dim(gxe_dt)
  gxe_gene <- gxe_dt[phenotype_id == gene_r]
  gxe_gene
  dim(gxe_gene)
  # Merge with GWAS positions for plotting coordinates
  gxe_sub <- merge(gxe_gene, gwas_full[, .(panel_variant_id, chromosome, position)],
                   by.x = "variant_id", by.y = "panel_variant_id")
  gxe_sub <- gxe_sub[chromosome == chr_r & position %between% c(start_r, end_r)]
  gxe_sub[, neglog10p := -log10(pval)]
  gxe_sub[, neglog10p_main := -log10(pval_g          )]
  gxe_sub
  # 4. Load Main Effect (eQTL) Data
  # We extract cell type and modality from the GxE filename to find the ME file
  gxe_bn    <- basename(gxe_path)
  ct_mod    <- str_match(gxe_bn, "^([^_]+)_([^_]+)_")
  cell_type <- ct_mod[1, 2]
  modality  <- ct_mod[1, 3]
  
  # Use your existing logic to find the main_effect file
  
  
  
  # 5. Coordinate Alignment & Scaling
  gxe_sub
  end_r   <- max(gxe_sub$position,gwas_sub$position)
  end_r
  xmin_mb <- start_r / 1e6
  xmax_mb <- end_r / 1e6
  ymax    <- max(c(gwas_sub$neglog10p, gxe_sub$neglog10p), na.rm = TRUE) * 1.1
  
  # Identify Lead SNP row in each dataset for highlighting
  get_lead <- function(df) df[which.min(abs(position - lead_pos))]
  ls_gwas  <- get_lead(gwas_sub)
  ls_gxe   <- get_lead(gxe_sub)
  
  # 6. Build the 3 Panels
  th <- theme_bw() + theme(panel.grid.minor = element_blank(), plot.title = element_text(size = 10, face = "bold"))
  
  
  my_symbol <- get_symbol(gene_r, gene_map)
  my_symbol
  
  gwas_name<-str_replace_all(str_replace_all(gwas_name,"imputed_",""),"._active_loci","")
  module_name<-str_replace_all(str_replace_all(hdl_row$Seed_term,"ME",""),"_Eigen","")
  module_name
  
  p1 <- ggplot(gwas_sub, aes(position / 1e6, neglog10p)) +
    geom_point(size = 0.8, alpha = 0.5, color = "grey40") +
    geom_point(data = ls_gwas, aes(position / 1e6, neglog10p), size = 3, color = "red", shape = 18) +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "red", linewidth = 0.4) +
    coord_cartesian(xlim = c(xmin_mb, xmax_mb), ylim = c(0, ymax)) +
    labs(title = paste0("GWAS: ", gwas_name), y = expression(-log[10](p)), x = NULL) + th
  p1
  maximing_hight_ylim=ifelse(ymax>max(gxe_sub$neglog10p)+30,
                             max(gxe_sub$neglog10p),ymax)
  p2 <- ggplot(gxe_sub, aes(position / 1e6, neglog10p)) +
    geom_point(size = 0.8, alpha = 0.5, color = "steelblue") +
    geom_point(data = ls_gxe, aes(position / 1e6, neglog10p), size = 3, color = "red", shape = 18) +
    geom_hline(yintercept = -log10(1e-5), linetype = "dashed", color = "blue", linewidth = 0.4) +
    coord_cartesian(xlim = c(xmin_mb, xmax_mb), ylim = c(0, maximing_hight_ylim)) +
    labs(title = paste0("GxE Interaction: ", cell_type, " ", modality, " - ", my_symbol," - ",module_name), 
         y = expression(-log[10](p[int])), x = NULL) + th
  p2
  maximing_hight_ylim_main=ifelse(ymax>max(gxe_sub$neglog10p_main)+30,
                                  max(gxe_sub$neglog10p_main),ymax)
  p3 <- ggplot(gxe_sub, aes(position / 1e6, neglog10p_main)) +
    geom_point(size = 0.8, alpha = 0.5, color = "darkorange") +
    geom_point(data = ls_gxe, aes(position / 1e6, neglog10p_main), size = 3, color = "red", shape = 18) +
    geom_hline(yintercept = -log10(1e-5), linetype = "dashed", color = "blue", linewidth = 0.4) +
    coord_cartesian(xlim = c(xmin_mb, xmax_mb), ylim = c(0, maximing_hight_ylim_main)) +
    labs(title = paste0("G Main Effect:", cell_type, " ", modality, " - ", my_symbol), y = expression(-log[10](p[eQTL])), 
         x = paste0("Position on ", chr_r, " (Mb)")) +th
  p3
  # 7. Stack and Return
  combined  <- cowplot::plot_grid(p1, p2, p3, ncol = 1, align = "v", rel_heights = c(1, 1, 1.1))
  
  final_plot <- cowplot::ggdraw() +
    cowplot::draw_plot(combined, y = 0, height = 0.96)
  final_plot
  return(final_plot)
}

# 1. Prepare your candidates
# IBD Selection
ibd_row <- candidates[gwas_name == "imputed_IBD.EUR.Crohns_Disease._active_loci"][1]
# 2. Generate Plots
plot_ibd <- get_coloc_plot(ibd_row)
plot_ibd
ggsave(plot=plot_ibd,
       filename="/labs/smontgom/grps_smontgom/hkrupkin/GxE/plots/plots_for_figures/coloc_examples/test_ibd.png",
       width=10,
       height=10,
       dpi=300)



# LDL Selection
ldl_row <- candidates[gwas_name == "imputed_GLGC_Mc_LDL._active_loci"][1]
plot_ldl <- get_coloc_plot(ldl_row)
plot_ldl
hdl_row<- candidates[gwas_name == "imputed_GLGC_Mc_HDL._active_loci"][5]
plot_hdl <- get_coloc_plot(hdl_row)
plot_hdl

ggsave(plot=plot_hdl,
       filename="/labs/smontgom/grps_smontgom/hkrupkin/GxE/plots/plots_for_figures/coloc_examples/test_plot_hdl.png",
       width=7,
       height=5,
       dpi=300)
###loopdiloop####
output_dir <- "/labs/smontgom/grps_smontgom/hkrupkin/GxE/plots/plots_for_figures/coloc_examples/all_colocs"

# Create directory if it doesn't exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Loop through all hits in candidates
good_candidates<-candidates%>%
  dplyr::filter(Lead_SNP_matches_GWAS_lead==T & Lead_P_int<1^-8)

dim(good_candidates)

for (i in 1:nrow(good_candidates)) {
  
  # Extract relevant info from current hit
  hit <- good_candidates[i, ]
  
  # Generate plot (modify this based on your plotting function)
  p <- get_coloc_plot(hit)  # Replace with your actual plotting function
  
  # Save plot
  filename <- paste0(output_dir, "/", hit$Gene, "_", hit$gwas_name, "_coloc.png")  # Adjust naming as needed
  print("file name is :")
  print(filename)
  ggsave(filename, plot = p, width = 8, height = 6, dpi = 300)
  
  message(paste("Saved plot", i, "of", nrow(good_candidates)))
}




# 3. View or Save
print(plot_ibd)
#ggsave("IBD_coloc_plot.pdf", plot_ibd, width = 8, height = 9)
### identifying what are the componenets of the module####
# ── 0. Configuration ──────────────────────────────────────────────────────────
enriched_annoations<-fread("/labs/smontgom/grps_smontgom/hkrupkin/GxE/Code/METABOLITE_PATHWAY_ERNICHEMTN/outputV2/enrichment_ALL_MODULES.csv")
enriched_annoations
ANNOTATION_FILE<-fread("/labs/smontgom/grps_smontgom/mahuynh/metabolites/MESA/data/annotation_metabolites_met_ID/annotation_table2.csv")
MODULES_DIR     <- "/oak/stanford/groups/smontgom/mahuynh/metabolites/MESA/data/associated_features/Exam 1 - metabolites/Modules"
module_files <- list.files(MODULES_DIR, pattern = "_members\\.txt$", full.names = TRUE)
message(sprintf("  Found %d module files", length(module_files)))
list_of_metabolties<-c()
running_index<-0
for (module_file in module_files){
  running_index<-running_index+1
  print(running_index)
  print(module_file)
  module<-fread(module_file,header=F,sep="\t")
  module_name<-str_replace_all(basename(module_file),"_members.txt","")
  module_name<-substring(module_name, 2)
  module_name
  module$name<-module_name
  colnames(module)[1]<-"Metabolite_ID"
  list_of_metabolties[[module_name]]<-module
  
  
}
all_metabolties_in_all_annotations <- do.call(rbind, list_of_metabolties)
#View(all_metabolties_in_all_annotations)
all_metabolties_in_all_annotations_with_annotation<-merge(all_metabolties_in_all_annotations,ANNOTATION_FILE,
                                                          by.x="Metabolite_ID",
                                                          by.y="name")
