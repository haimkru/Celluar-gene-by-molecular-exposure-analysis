###Examining_GxE_Variants#######
#author:haim krupkin
#date: 01/06/2025
#description:this script loads everyone of the sites that are significnat

#pseudo code



#notes

###packages####
library(ggplot2)
library(ggpubr)
library(dplyr)
library(data.table)
library(readr)
library(stringr)
library(dplyr)
library(ggplot2)
library(stringr)
# Define a color palette based on the unique terms
color_palette <- c(
  "bisque4"        = "bisque4",
  "black"          = "black",
  "blue"           = "blue",
  "brown4"         = "brown4",
  "brown"          = "brown",
  "cyan"           = "cyan",
  "darkgreen"      = "darkgreen",
  "darkgrey"       = "darkgrey",
  "darkmagenta"    = "darkmagenta",
  "darkolivegreen" = "darkolivegreen",
  "darkorange2"    = "darkorange2",
  "darkorange"     = "darkorange",
  "darkred"        = "darkred",
  "darkslateblue"  = "darkslateblue",
  "darkturquoise"  = "darkturquoise",
  "floralwhite"    = "floralwhite",
  "green"          = "green",
  "greenyellow"    = "greenyellow",
  "grey60"         = "grey60",
  "ivory"          = "ivory",
  "lightcyan1"     = "lightcyan1",
  "lightcyan"      = "lightcyan",
  "lightgreen"     = "lightgreen",
  "lightsteelblue1" = "lightsteelblue1",
  "lightyellow"    = "lightyellow",
  "magenta"        = "magenta",
  "mediumpurple3"  = "mediumpurple3",
  "midnightblue"   = "midnightblue",
  "orange"         = "orange",
  "orangered4"     = "orangered4",
  "paleturquoise"  = "paleturquoise",
  "pink"           = "pink",
  "plum1"          = "plum1",
  "plum2"          = "plum2",
  "purple"         = "purple",
  "red"            = "red",
  "royalblue"      = "royalblue",
  "saddlebrown"    = "saddlebrown",
  "salmon4"        = "salmon4",
  "salmon"         = "salmon",
  "sienna3"        = "sienna3",
  "skyblue3"       = "skyblue3",
  "skyblue"        = "skyblue",
  "steelblue"      = "steelblue",
  "tan"            = "tan",
  "thistle1"       = "thistle1",
  "thistle2"       = "thistle2",
  "turquoise"      = "turquoise",
  "violet"         = "violet",
  "white"          = "white",
  "yellow"         = "yellow",
  "yellowgreen"    = "yellowgreen"
)

options(scipen = 999)


####
#loading file paths###
file_paths<-list.files("/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_outputs",
           full.names = T,
           pattern="sig.tsv")
###loading each file using a loop
all_hits<-data.frame()
file_number=0
for (file_path in file_paths){
  #file_path<-file_paths[1]
  print("file number : ")
  
  print(file_number)
  file_number=file_number+1
  print("out of total : ")
  print(length(file_paths))
  print(file_path)
  sig_hits<-fread(file_path)
  sig_hits$path<-file_path
  all_hits<-rbind(all_hits,sig_hits)
}

#saving data
all_hits <- all_hits %>%
  mutate(
    filename = basename(path),
    Omics = case_when(
      str_detect(filename, "^metabolome") ~ "Metabolome",
      str_detect(filename, "^proteome")   ~ "Proteome",
      TRUE ~ NA_character_
    ),
    assay = case_when(
      str_detect(filename, "_expression_")   ~ "Expression",
      str_detect(filename, "_methylation_")  ~ "Methylation",
      TRUE ~ NA_character_
    )
  )
all_hits <- all_hits %>%
  mutate(
    cell_type = case_when(
      str_detect(filename, "_bcell_") ~ "B Cells",
      str_detect(filename, "_bulk_")  ~ "Bulk",
      str_detect(filename, "_cd4t_")  ~ "CD4⁺ T",
      str_detect(filename, "_cd8t_")  ~ "CD8⁺ T",
      str_detect(filename, "_mono_")  ~ "Mono",
      str_detect(filename, "_neu_")   ~ "Neutro",
      str_detect(filename, "_nk_")    ~ "NK",
      TRUE ~ NA_character_
    )
  )
b<-all_hits
#all_hits%>%
#  dplyr::group_by(variant_id,phenotype_id)%>%
#  dplyr::mutate(pvalues_distinct=n_distinct(term))%>%
#  View()
example_df<-all_hits%>%
  head(n=10000)%>%
  mutate(effect_type = case_when(
    term == "pval_g" ~ "main_effect",
    grepl("^pval_g-ME", term) ~ "interaction_effect",
    grepl("^pval_ME", term) ~ "environment_effect",
    TRUE ~ NA_character_
  ))

all_hits<-all_hits%>%
  mutate(effect_type = case_when(
    term == "pval_g" ~ "main_effect",
    grepl("^pval_g-ME", term) ~ "interaction_effect",
    grepl("^pval_ME", term) ~ "environment_effect",
    TRUE ~ NA_character_
  ))



all_hits$term<-str_replace(all_hits$term,"pval_g-","")
all_hits$term<-str_replace(all_hits$term,"pval_","")
all_hits$term<-str_replace(all_hits$term,"_Protein","")
all_hits$term<-str_replace(all_hits$term,"_Metabolite","")

table(all_hits$Omics)
table(all_hits$assay)
table(all_hits$cell_type)


readr::write_tsv(all_hits,
                 "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/all_cell_types_sig_hits.tsv")

small_all_hits<-fread("/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/all_cell_types_sig_hits.tsv")

###
head(all_hits)
#gettting the ctual file type and such#
basename(unique(all_hits$path))

###
table_of_hits_per_cell_type_per_omics_per_assay<-all_hits%>%
  group_by(cell_type,assay,Omics)%>%
  dplyr::summarise(number_of_hits_in_cell_type_omic_assay=n())
  
plot_n_GxE_by_cell_type_assay_omics<-table_of_hits_per_cell_type_per_omics_per_assay%>%
  as.data.frame()%>%
  ggplot(aes(x=cell_type,y=number_of_hits_in_cell_type_omic_assay,fill=Omics))+
  geom_bar(stat="identity",position="dodge")+
  labs(x="Cell Type",y="Number of GxE Hits")+
  facet_wrap(~assay,scale="free_y",ncol=1)+
  theme_bw()+
  theme(
    axis.text.x  = element_text(color = "black", size = 12),
    axis.text.y  = element_text(color = "black", size = 12),
    axis.title.x = element_text(color = "black", size = 14),
    axis.title.y = element_text(color = "black", size = 14),
    strip.text   = element_text(color = "black", size = 13),
    
    axis.line  = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    
    legend.position = "bottom",
    legend.justification = "center"
  )
path_of_plots<-"/labs/smontgom/grps_smontgom/hkrupkin/GxE/plots/"
plot_n_GxE_by_cell_type_assay_omics

ggsave(plot=plot_n_GxE_by_cell_type_assay_omics,
       filename = paste0(path_of_plots,"plot_n_GxE_by_cell_type_assay_omics.png"),
       height=10,
       width=6,
       dpi=900)
###plotting per cell type per omics per assay per eigen####

table_of_hits_per_cell_type_per_omics_per_assay_per_eigen<-all_hits%>%
  group_by(cell_type,assay,Omics,term)%>%
  dplyr::summarise(number_of_hits_in_cell_type_omic_assay=n())

# Load necessary libraries

# Prepare the data

##plotting bulk###
data_filtered <- table_of_hits_per_cell_type_per_omics_per_assay_per_eigen %>%
  mutate(term = str_replace(term, "_Eigen", ""),
         term = str_replace(term, "ME", "")) %>%
  filter(cell_type == "Bulk")

# Extract unique terms for coloring
unique_terms <- unique(data_filtered$term)


# Create the ggplot
ggplot(data_filtered, aes(x = term, y = number_of_hits_in_cell_type_omic_assay, fill = term)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Eigen", y = "Number of GxE Hits For Bulk") +
  scale_fill_manual(values = color_palette) +  # Apply the custom color palette
  facet_wrap(~ assay + Omics, scales = "free") +  # Correct facet_wrap syntax
  theme_bw() +
  theme(
    axis.text.x  = element_text(color = "black", size = 12, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y  = element_text(color = "black", size = 12),
    axis.title.x = element_text(color = "black", size = 14),
    axis.title.y = element_text(color = "black", size = 14),
    strip.text   = element_text(color = "black", size = 13),
    axis.line    = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    legend.position = "none"
  )
###plotting monocytes###
data_filtered <- table_of_hits_per_cell_type_per_omics_per_assay_per_eigen %>%
  mutate(term = str_replace(term, "_Eigen", ""),
         term = str_replace(term, "ME", "")) %>%
  filter(cell_type == "Mono")

# Extract unique terms for coloring
unique_terms <- unique(data_filtered$term)
table_of_hits_per_cell_type_per_omics_per_assay_per_eigen %>%
  ggplot(aes(x = number_of_hits_in_cell_type_omic_assay)) +
  geom_histogram(bins = 30) +  # Specify the number of bins
  scale_x_log10(breaks = c(1, 10, 100, 1000, 10000, 100000)) +  # Custom breaks in a log scale
  labs(
    x = "Number of GxE per Cell Type Per Assay",  # Correct x-axis label
    y = "Frequency",  # Adding y-axis label for clarity
    title = "Histogram of GxE Hits per Cell Type Per Assay"  # Adding title for the plot
  ) +
  theme_bw()+
  facet_wrap(cell_type~assay,ncol=2)

# Create the ggplot
ggplot(data_filtered, aes(x = term, y = number_of_hits_in_cell_type_omic_assay, fill = term)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Eigen", y = "Number of GxE Hits For Mono") +
  scale_fill_manual(values = color_palette) +  # Apply the custom color palette
  facet_wrap(~ assay + Omics, scales = "free") +  # Correct facet_wrap syntax
  theme_bw() +
  theme(
    axis.text.x  = element_text(color = "black", size = 12, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y  = element_text(color = "black", size = 12),
    axis.title.x = element_text(color = "black", size = 14),
    axis.title.y = element_text(color = "black", size = 14),
    strip.text   = element_text(color = "black", size = 13),
    axis.line    = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    legend.position = "none"
  )
###

# Assume 'table_of_hits_per_cell_type_per_omics_per_assay_per_eigen' is your original data frame
# Let's calculate the ratios and filter for Bulk
data_filtered <- table_of_hits_per_cell_type_per_omics_per_assay_per_eigen %>%
  filter(cell_type == "Bulk") %>%
  mutate(
    # Assuming 'number_of_hits_in_cell_type_omic_assay' is the count for each category
    ratio = case_when(
      Omics == "Proteome" & assay == "Expression" ~ number_of_hits_in_cell_type_omic_assay / 
        sum(number_of_hits_in_cell_type_omic_assay[Omics == "Proteome" & assay == "Expression"]),
      Omics == "Metabolome" & assay == "Expression" ~ number_of_hits_in_cell_type_omic_assay /
        sum(number_of_hits_in_cell_type_omic_assay[Omics == "Metabolome" & assay == "Expression"]),
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(ratio)) %>%
  arrange(ratio) # Order by ratio

# Create the bar plot
ggplot(data_filtered,
       aes(x = term, y = number_of_hits_in_cell_type_omic_assay,fill=Omics)) +
  geom_bar(stat="identity")
  
data_filtered <- table_of_hits_per_cell_type_per_omics_per_assay_per_eigen %>%
  filter(cell_type == "Bulk") %>%
  group_by(term, Omics) %>%
  summarise(number_of_hits = sum(number_of_hits_in_cell_type_omic_assay, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(number_of_hits = ifelse(Omics == "Proteome", -number_of_hits, number_of_hits)) # Make Proteome negative, keep Metabolome positive
unique(data_filtered$term)



head(data_filtered)

# Create the basic pyramid plot
hits_pyramid <- ggplot(data_filtered, aes(x = term, fill = Omics, y = number_of_hits)) +
  geom_bar(stat = "identity") + 
  scale_y_continuous(labels = abs, limits = c(-max(abs(data_filtered$number_of_hits)), max(data_filtered$number_of_hits))) +  # Set y limits
  coord_flip() +  # Flip the coordinates for a horizontal bar chart
  theme_minimal() + 
  labs(
    x = "Eigen",
    y = "Number of GxE Hits For Bulk",
    fill = "Omics",
    title = "Population Pyramid of GxE Hits"
  )+
  theme(
    axis.text.x  = element_text(color = "black", size = 5),
    axis.text.y  = element_text(color = "black", size = 5),
    axis.title.x = element_text(color = "black", size = 14),
    axis.title.y = element_text(color = "black", size = 14),
    strip.text   = element_text(color = "black", size = 13),
    
    axis.line  = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    
    legend.position = "bottom",
    legend.justification = "center"
  )
hits_pyramid

###plotting the number of genes per eigenFeature*CellType*Ome
table_of_hits_per_cell_type_per_omics_per_assay_per_eigen<-all_hits%>%
  dplyr::select(cell_type,assay,Omics,phenotype_id)%>%
  distinct()%>%
  dplyr::filter(assay=="Expression")%>%
  group_by(cell_type,assay,Omics)%>%
  View()
  dplyr::summarise(number_of_hits_in_cell_type_omic_assay=n())
table_of_hits_per_cell_type_per_omics_per_assay_per_eigen%>%
  ggplot(aes(x=cell_type,y=number_of_hits_in_cell_type_omic_assay,fill=Omics))+
  geom_bar(stat="identity",position="dodge")

####

table_of_hits_per_cell_type_per_omics_per_assay_per_eigen<-small_all_hits%>%
  dplyr::filter(effect_type=="interaction_effect")%>%
  group_by(cell_type,assay,Omics,term)%>%
  dplyr::summarise(number_of_hits_in_cell_type_omic_assay=n_distinct(phenotype_id))

data_filtered <- table_of_hits_per_cell_type_per_omics_per_assay_per_eigen %>%
  mutate(term = str_replace(term, "_Eigen", ""),
         term = str_replace(term, "ME", "")) %>%
  dplyr::filter(cell_type == "Mono")

# Extract unique terms for coloring
unique_terms <- unique(data_filtered$term)
table_of_hits_per_cell_type_per_omics_per_assay_per_eigen %>%
  ggplot(aes(x = number_of_hits_in_cell_type_omic_assay)) +
  geom_histogram(bins = 30) +  # Specify the number of bins
  scale_x_log10(breaks = c(1, 10, 100, 1000, 10000, 100000)) +  # Custom breaks in a log scale
  labs(
    x = "Number of GxE per Cell Type Per Assay",  # Correct x-axis label
    y = "Frequency",  # Adding y-axis label for clarity
    title = "Histogram of GxE Hits per Cell Type Per Assay"  # Adding title for the plot
  ) +
  theme_bw()+
  facet_wrap(cell_type~assay,ncol=2)

# Create the ggplot
plo_v2_number_of_hits<-ggplot(data_filtered, aes(x = term, y = number_of_hits_in_cell_type_omic_assay, fill = term)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Module", y = "Number of GxE Hits For Mono") +
  scale_fill_manual(values = color_palette) +  # Apply the custom color palette
  facet_wrap(~ assay + Omics, scales = "free") +  # Correct facet_wrap syntax
  theme_bw() +
  theme(
    axis.text.x  = element_text(color = "black", size = 12, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y  = element_text(color = "black", size = 12),
    axis.title.x = element_text(color = "black", size = 14),
    axis.title.y = element_text(color = "black", size = 14),
    strip.text   = element_text(color = "black", size = 13),
    axis.line    = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    legend.position = "none"
  )
plo_v2_number_of_hits
ggsave(plot=plo_v2_number_of_hits,
       file="/labs/smontgom/grps_smontgom/hkrupkin/GxE/plots/plots_for_figures/plot_of_GxE_Hits_per_cell_type_per_omics_per_assay_per_eigenV2.png",
       width=15,
       height=7,
       dpi=900)  

###plotting number of hits sum across all cell types#######
###plotting the number of genes per eigenFeature*CellType*Ome
all_hits<-fread("/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/all_cell_types_sig_hits.tsv")

###getting number of hits per term###
count_per_MEM_discovery_rate<-all_hits%>%
  dplyr::select(cell_type,assay,Omics,phenotype_id,term)%>%
  group_by(cell_type,assay,Omics,term)%>%
  dplyr::summarise(number_of_hits_in_cell_type_omic_assay=n_distinct(phenotype_id))

count_per_MEM_discovery_rate<-count_per_MEM_discovery_rate%>%
  mutate(n_tests=ifelse(assay=="Expression",16750,159161))%>%
  mutate(percentage_per_mem=number_of_hits_in_cell_type_omic_assay/n_tests)
#View(count_per_MEM_discovery_rate)
library(ggbeeswarm)

aggregate_count_per_MEM_discovery_rate<-count_per_MEM_discovery_rate %>%
  mutate(percentage_per_mem=percentage_per_mem*100)%>%
  mutate(Assay=assay)%>%
  dplyr::select(-assay)%>%
  mutate(Omics=Omics)%>%
  ungroup()%>%
  dplyr::select(Assay,Omics,percentage_per_mem)%>%
  group_by(Assay, Omics) %>%
  dplyr::summarise(
    Median_Percent_of_Hits_per_Ome_Cell_Type = median(percentage_per_mem, na.rm = TRUE),
    Top75_Percent_of_Hits_per_Ome_Cell_Type = quantile(percentage_per_mem, 0.75, na.rm = TRUE),
    Top25_Percent_of_Hits_per_Ome_Cell_Type = quantile(percentage_per_mem, 0.25, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(aggregate_count_per_MEM_discovery_rate,
          "/oak/stanford/groups/smontgom/hkrupkin/GxE/Data/aggregate_count_per_MEM_discovery_rate.csv")

count_per_MEM_discovery_rate%>%
  dplyr::filter(term!="g")%>%
  mutate(facet_group = interaction(assay, Omics)) %>%
  group_by(facet_group) %>%
  arrange(number_of_hits_in_cell_type_omic_assay) %>%
  mutate(order = row_number()) %>%
  ungroup() %>%
  # Create a unique term label per facet (to handle same term in multiple facets)
  mutate(term_facet = paste(term, facet_group, sep = "___")) %>%
  mutate(term_facet = reorder(term_facet, order))%>%
  ggplot(aes(x =  Omics, y = percentage_per_mem)) +
  geom_boxplot()+
  geom_quasirandom(width = 0.2, alpha = 0.6)+
  labs(x = "Module", y = "Unique GxE Hits per Cell Types-Term") +
  scale_x_discrete(labels = function(x) gsub("___.*$", "", x)) +  # Strip the suffix for clean labels
  scale_fill_manual(values = color_palette) +
  scale_y_continuous(labels = scales::label_percent())+
  facet_wrap(~assay , scales = "free") +
  theme_bw() +
  theme(
    axis.text.x  = element_text(color = "black", size = 8),
    axis.text.y  = element_text(color = "black", size = 8),
    axis.title.x = element_text(color = "black", size = 8),
    axis.title.y = element_text(color = "black", size = 8),
    strip.text   = element_text(color = "black", size = 13),
    axis.line    = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    legend.position = "none"
  )

table_of_hits_per_cell_type_per_omics_per_assay_per_eigen<-all_hits%>%
  dplyr::select(cell_type,assay,Omics,phenotype_id)%>%
  distinct()%>%
  dplyr::filter(assay=="Expression")%>%
  group_by(cell_type,assay,Omics)%>%
  dplyr::summarise(number_of_hits_in_cell_type_omic_assay=n())
table_of_hits_per_cell_type_per_omics_per_assay_per_eigen%>%
  ggplot(aes(x=cell_type,y=number_of_hits_in_cell_type_omic_assay,fill=Omics))+
  geom_bar(stat="identity",position="dodge")

####

table_of_hits_per_per_omics_per_assay_per_eigen<-all_hits%>%
  dplyr::filter(effect_type=="interaction_effect")%>%
  dplyr::filter(cell_type!="Bulk")%>%
  group_by(assay,Omics,term)%>%
  dplyr::summarise(number_of_hits_in_cell_type_omic_assay=n_distinct(phenotype_id))

data_filtered <- table_of_hits_per_per_omics_per_assay_per_eigen %>%
  mutate(term = str_replace(term, "_Eigen", ""),
         term = str_replace(term, "ME", "")) 
# Extract unique terms for coloring
unique_terms <- unique(data_filtered$term)


# Create the ggplot
data_for_plot <- data_filtered %>%
  mutate(facet_group = interaction(assay, Omics)) %>%
  group_by(facet_group) %>%
  arrange(number_of_hits_in_cell_type_omic_assay) %>%
  mutate(order = row_number()) %>%
  ungroup() %>%
  # Create a unique term label per facet (to handle same term in multiple facets)
  mutate(term_facet = paste(term, facet_group, sep = "___")) %>%
  mutate(term_facet = reorder(term_facet, order))

plo_v2_number_of_hits <- data_for_plot %>%
  ggplot(aes(x = term_facet, y = number_of_hits_in_cell_type_omic_assay, fill = term)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Module", y = "Number of Unique GxE Hits Across Cell Types") +
  scale_x_discrete(labels = function(x) gsub("___.*$", "", x)) +  # Strip the suffix for clean labels
  scale_fill_manual(values = color_palette) +
  facet_wrap(~ assay + Omics, scales = "free") +
  theme_bw() +
  theme(
    axis.text.x  = element_text(color = "black", size = 12, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y  = element_text(color = "black", size = 12),
    axis.title.x = element_text(color = "black", size = 14),
    axis.title.y = element_text(color = "black", size = 14),
    strip.text   = element_text(color = "black", size = 13),
    axis.line    = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    legend.position = "none"
  )

plo_v2_number_of_hits

ggsave(plot = plo_v2_number_of_hits,
       file = "/labs/smontgom/grps_smontgom/hkrupkin/GxE/plots/plots_for_figures/plot_of_GxE_Hits_per_omics_per_assay_per_eigenV2.png",
       width = 15,
       height = 7,
       dpi = 900)

