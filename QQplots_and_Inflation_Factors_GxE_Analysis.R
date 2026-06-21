#!/usr/bin/env Rscript

# ===========================================
# LOAD LIBRARIES
# ===========================================
suppressPackageStartupMessages({
  library(data.table)  # fast fread, rbindlist
  library(dplyr)       # mutate, case_when
  library(stringr)     # str_detect, str_extract
  library(tidyr)       # pivot_longer
})

data.table::setDTthreads(1)  # single-threaded to reduce memory spikes

# ===========================================
# HELPER FUNCTIONS
# ===========================================
ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

log_step <- function(msg) {
  message(sprintf("[%s] [STEP] %s", ts(), msg))
}

# memory report: R used + system free RAM (Linux)
mem_report <- function(stage) {
  g <- gc()
  used_gb <- sum(g[, "used"]) * 8 / 1024^2  # approx GB
  free_gb <- NA
  if (file.exists("/proc/meminfo")) {
    meminfo <- readLines("/proc/meminfo")
    mem_free <- as.numeric(str_extract(meminfo[grep("MemAvailable", meminfo)], "\\d+"))
    free_gb <- mem_free / 1024 / 1024
  }
  message(sprintf("[%s] [MEM] %-50s | used=%.3f GB | free_sys=%.3f GB",
                  ts(), stage, used_gb, free_gb))
}

# weighted median for inflation factor calculation
weighted_median <- function(x, w) {
  stopifnot(length(x) == length(w))
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  x[which(cw >= 0.5)[1]]
}

# ===========================================
# DISCOVER CSV FILES
# ===========================================
log_step("Discovering CSV files")
mem_report("before list.files")

all_files <- list.files(
  "/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_outputs",
  recursive = TRUE,
  full.names = TRUE,
  pattern = "\\.csv$"
)

log_step(sprintf("Found %d files", length(all_files)))
mem_report("after list.files")

df_paths <- data.table(file_paths = all_files)
df_paths[, basename := basename(file_paths)]

# ===========================================
# ANNOTATE METADATA
# ===========================================
log_step("Annotating metadata")
mem_report("before annotation")

df_paths[, Omics := fifelse(str_detect(file_paths, "metabolome"), "Metabolome",
                            fifelse(str_detect(file_paths, "proteome"), "Proteome", NA_character_))]
df_paths[, assay := fifelse(str_detect(file_paths, "expression"), "Expression",
                            fifelse(str_detect(file_paths, "methylation"), "Methylation", NA_character_))]
df_paths[, chunk := fifelse(str_detect(file_paths, "chunk_1"), "chunk1",
                            fifelse(str_detect(file_paths, "chunk_2"), "chunk2",
                                    fifelse(str_detect(file_paths, "chunk_3"), "chunk3", NA_character_)))]
df_paths[, raw_cell_type := str_extract(basename, "^[^_]+")]
df_paths[, cell_type := case_when(
  str_detect(raw_cell_type, "bcell") ~ "B Cells",
  str_detect(raw_cell_type, "bulk")  ~ "Bulk",
  str_detect(raw_cell_type, "cd4t")  ~ "CD4⁺ T",
  str_detect(raw_cell_type, "cd8t")  ~ "CD8⁺ T",
  str_detect(raw_cell_type, "mono")  ~ "Mono",
  str_detect(raw_cell_type, "neu")   ~ "Neutro",
  str_detect(raw_cell_type, "nk")    ~ "NK",
  TRUE ~ NA_character_
)]

mem_report("after annotation")

# ===========================================
# PARAMETER GRID
# ===========================================
param_grid <- df_paths %>%
  distinct(Omics, assay, cell_type, chunk)

log_step(sprintf("Parameter rows: %d", nrow(param_grid)))
mem_report("after param grid")

# ===========================================
# MAIN LOOP: per parameter set
# ===========================================
results <- list()
res_i <- 1

for (i in seq_len(nrow(param_grid))) {
  params <- param_grid[i, ]
  log_step(sprintf("PARAM %d/%d :: %s | %s | %s | %s", i, nrow(param_grid),
                   params$Omics, params$assay, params$cell_type, params$chunk))
  mem_report("start param")
  
  # select files for this param set
  files <- df_paths[Omics == params$Omics & assay == params$assay &
                      cell_type == params$cell_type & chunk == params$chunk, file_paths]
  
  if (length(files) == 0) {
    log_step("No files — skipping")
    next
  }
  
  log_step(sprintf("Using %d chromosome files:", length(files)))
  print(files)  # print all file paths
  
  # initialize lists to store per pval column medians and counts across chromosomes
  col_medians <- list()
  col_counts  <- list()
  
  # ===========================================
  # PROCESS EACH CHROMOSOME FILE
  # ===========================================
  for (f_i in seq_along(files)) {
    f <- files[f_i]
    log_step(sprintf("Reading file %d/%d: %s", f_i, length(files), basename(f)))
    mem_report("before fread")
    print(f)
    dt <- fread(f, showProgress = TRUE)  # read entire file
    mem_report(sprintf("after fread: %d rows x %d cols", nrow(dt), ncol(dt)))
    
    # detect pval columns
    pval_cols <- grep("^pval", names(dt), value = TRUE)
    log_step(sprintf("Detected %d pval columns: %s", length(pval_cols), paste(head(pval_cols,10), collapse=", ")))
    
    # convert to long format for group_by-like computation
    dt_long <- melt(dt, measure.vars = pval_cols,
                    variable.name = "pval_column", value.name = "pval")
    mem_report(sprintf("After pivot_longer: %d rows x %d cols", nrow(dt_long), ncol(dt_long)))
    
    # compute median chi-squared and row counts per pval column
    file_stats <- dt_long[, .(
      median_chi = median(qchisq(1 - pval, df = 1)),
      n_rows     = .N
    ), by = .(pval_column)]
    
    # accumulate into lists
    for (row in 1:nrow(file_stats)) {
      pcol <- file_stats$pval_column[row]
      col_medians[[pcol]] <- c(col_medians[[pcol]], file_stats$median_chi[row])
      col_counts[[pcol]]  <- c(col_counts[[pcol]], file_stats$n_rows[row])
    }
    
    # clean up
    rm(dt, dt_long, file_stats); gc()
    mem_report("after processing file")
  }
  
  # ===========================================
  # COMPUTE WEIGHTED MEDIAN LAMBDA GC per pval column
  # ===========================================
  for (pcol in names(col_medians)) {
    lambda_gc <- weighted_median(col_medians[[pcol]], col_counts[[pcol]]) / 0.4549364
    log_step(sprintf("Final lambda (%s) = %.6f | total rows = %d",
                     pcol, lambda_gc, sum(col_counts[[pcol]])))
    
    results[[res_i]] <- data.table(
      pval_column = pcol,
      lambda_gc   = lambda_gc,
      Omics       = params$Omics,
      assay       = params$assay,
      cell_type   = params$cell_type,
      chunk       = params$chunk,
      n_total     = sum(col_counts[[pcol]])
    )
    res_i <- res_i + 1
  }
  
  rm(col_medians, col_counts); gc()
  mem_report("end param")
}

# ===========================================
# COMBINE AND WRITE RESULTS
# ===========================================
log_step("Binding results")
mem_report("before rbindlist")

final <- rbindlist(results)
mem_report("after rbindlist")

fwrite(final,
       "/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/Outputs/all_inflation_factors.tsv",
       sep = "\t")

log_step(sprintf("DONE — rows written: %d", nrow(final)))
mem_report("final")
