#Annotating_Methylation_DataV2.R
#/oak/stanford/groups/smontgom/hkrupkin/Code/Annotating_Methylation_DataV2.R
#author: Haim Krupkin
#date: 01/30/2026

#description: we improve uon the previous version as the previous version doesnt annotate methylation correctly
###packages####
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(ChIPseeker)
library(EnsDb.Hsapiens.v86)
library(tidyverse)

############

###anottating methylations#####
### input ####
all_hits <- fread("/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/all_cell_types_sig_hits.tsv")

Expression_Data <- all_hits %>% dplyr::filter(assay == "Expression")
Methylation_Data <- all_hits %>% dplyr::filter(assay == "Methylation")

locations_methylations <- fread("/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/methylation/epic_manifest.tsv")

Methylation_Data_with_positions <- merge(
  as.data.table(Methylation_Data),
  as.data.table(locations_methylations),
  by.x = "phenotype_id",
  by.y = "Probe_ID"
)


####
peaks <- GRanges(
  seqnames = Methylation_Data_with_positions$seqnames,
  ranges = IRanges(start = Methylation_Data_with_positions$start, 
                   end = Methylation_Data_with_positions$end)
)
seqlevelsStyle(peaks) <- "Ensembl"
peaks
peakAnno <- annotatePeak(
  peaks,
  tssRegion = c(-1000, 100),
  TxDb = EnsDb.Hsapiens.v86,
  level = "gene"  # Annotate at gene level
)

result <- as.data.frame(peakAnno)
table(result$annotation)
#result%>%View()
result <- result %>%
  mutate(
    feature_type = case_when(
      str_detect(annotation, "^Exon") ~ "exonic",
      str_detect(annotation, "^Intron") ~ "intronic",
      str_detect(annotation, "^Promoter") ~ "promoter",
      TRUE ~ annotation
    )
  )
table(result$feature_type)
result$geneId[1:5]
gene_features <- c("exonic","intronic","promoter","3' UTR","5' UTR","Downstream (<=300bp)")

result2 <- result %>%
  mutate(
    geneId = if_else(feature_type %in% gene_features & !is.na(geneId) & geneId != "", geneId, NA_character_)
  )
table(is.na(result2$geneId))

result
head(Expression_Data)
genes_and_positon_meth<-result2%>%
  dplyr::select(seqnames,start,end,geneId)
genes_and_positon_meth
genes_and_positon_meth<-genes_and_positon_meth%>%
  mutate(variant_pos=paste0(seqnames,":",start,end))
Methylation_Data_with_positions<-Methylation_Data_with_positions%>%
  mutate(chrom_position=str_replace_all(seqnames,"chr",""))%>%
  mutate(variant_pos=paste0(chrom_position,":",start,end))
Methylation_Data_with_positions_with_gene<-Methylation_Data_with_positions%>%
  merge(genes_and_positon_meth%>%distinct(),
        by.x="variant_pos",
        by.y="variant_pos",
        all=T)
#
all_hits_methylation_and_expression_with_methylation_genes <- bind_rows(Expression_Data,
                                                                        Methylation_Data_with_positions_with_gene)

readr::write_tsv(all_hits_methylation_and_expression_with_methylation_genes,
                 "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/all_hits_methylation_and_expression_with_methylation_genesV2.tsv")























