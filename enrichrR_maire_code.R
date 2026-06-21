library(enrichR)
library(gtable, grid)
# Load the necessary libraries
library(ggplot2)
library(gridExtra)

#By default human genes are selected otherwise select your organism of choice.
websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}

#Then find the list of all available databases from Enrichr.
if (websiteLive) dbs <- listEnrichrDbs()

#Then, we select our databases of interest.
dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023", "Reactome_2022", "KEGG_2021_Human")

######################################################################
# We run enrichR on each of the modules for timestamp 1 and 5 to get the results
######################################################################

#### First, we load the data
datatime <- "1" #1 or 5
datatype <- "proteins"

#### We define the folder path where are located all the modules 
folder_path <- paste0("/oak/stanford/groups/smontgom/mahuynh/metabolites/MESA/data/associated_features/Exam ", datatime," - proteins/Modules")
# Use list.files to get a list of all files in the folder
file_list <- list.files(path = folder_path, full.names = TRUE)

cat(file_list[1])
#### We define a function to run enrichR on the top N genes of a module
run_enrichr <- function(file_path, top_n) {
  ############################## First, we read the file path
  # Initialize an empty list to store the genes
  gene_list <- character(0)
  # Check if the file exists
  if (!file.exists(file_path)) {
    cat("File not found.")
    return(NULL)
  }
  # Read the file line by line and add each gene to the list
  con <- file(file_path, open = "r")
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    gene_list <- c(gene_list, line)
  }
  close(con)
  ############################## Now, we run enrichr on that list
  enriched <- enrichr(gene_list, dbs)
  return(list(enriched = enriched, gene_list = gene_list))
}

# Define the vertical and horizontal spacing between genes
vertical_spacing <- 0.07
horizontal_spacing <- 0.09  # Adjust this value to increase horizontal space


for (file_path in file_list) {
  result <- run_enrichr(file_path, 100)
  enriched <- result$enriched
  gene_list <- result$gene_list
  save_path <- paste0("/oak/stanford/groups/smontgom/mahuynh/metabolites/MESA/data/associated_features/Exam ", datatime, " - proteins/enrichr")
  
  # We set the save name module_enrichr
  file_name <- basename(file_path)
  
  # Construct the save_name
  save_name <- gsub("\\..*", "_enrichr.RDS", file_name)
  
  # # Create a new PNG file for each set of 5 plots
  # save_file_path <- file.path(save_path, paste0(gsub("\\..*", "_enrichr", file_name), "_plots", ".pdf"))
  # pdf(file = save_file_path, width = 12, height = 7)  # Adjust dimensions and resolution as needed
  # 
  # genes_per_row = 10
  # # Add the gene list to the first page in rows
  # grid::grid.text(paste0("Gene List:", file_name), x = 0.2, y = 0.9)
  # y <- 0.85
  # x <- 0.1  # Initialize the x position
  # for (i in 1:length(gene_list)) {
  #   grid::grid.text(gene_list[i], x = x, y = y)
  #   x <- x + horizontal_spacing  # Increase horizontal spacing
  #   if (i %% genes_per_row == 0) {
  #     y <- y - vertical_spacing
  #     x <- 0.1  # Reset the x position for the next row
  #   }
  # }
  # 
  # p1 <- plotEnrich(enriched[[1]], showTerms = 20, numChar = 90, y = "Count", orderBy = "P.value", title = "GO_Molecular_Function_2023")
  # print(p1)
  # # Create a table from the dataframe and add it to the PDF page
  # # grid.table(enriched[["GO_Molecular_Function_2023"]])
  # 
  # p2 <- plotEnrich(enriched[[2]], showTerms = 20, numChar = 90, y = "Count", orderBy = "Adjusted.P.value", title = "GO_Cellular_Component_2023")
  # print(p2)
  # p3 <- plotEnrich(enriched[[3]], showTerms = 20, numChar = 90, y = "Count", orderBy = "Adjusted.P.value", title = "GO_Biological_Process_2023")
  # print(p3)
  # p4 <- plotEnrich(enriched[[4]], showTerms = 20, numChar = 90, y = "Count", orderBy = "Adjusted.P.value", title = "Reactome_2022")
  # print(p4)
  # p5 <- plotEnrich(enriched[[5]], showTerms = 20, numChar = 90, y = "Count", orderBy = "Adjusted.P.value", title = "KEGG_2021_Human")
  # print(p5)
  # 
  # dev.off()  
  
  # Save the enrichment result as an RDS file
  saveRDS(enriched, file.path(save_path, save_name))
  
  #Save each dataframe as a pdf
  printEnrich(enriched, prefix= paste0("/oak/stanford/groups/smontgom/mahuynh/metabolites/MESA/data/associated_features/Exam ", datatime, " - proteins/enrichr/", gsub("\\..*", "_enrichr", file_name)), outFile= "excel")
}