# 1. Helper functions

suppressPackageStartupMessages({
  library(tidyverse)
})
# Output files -------------------------------------------------------------------------------------------------------
dir_results <- "result"
dir_figures <- file.path(dir_results, "figures")
dir_tables  <- file.path(dir_results, "tables")
dir_rds     <- file.path(dir_results, "rds")
dir_logs    <- file.path(dir_results, "logs")

for (d in c(dir_results, dir_figures, dir_tables, dir_rds, dir_logs)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Resume the pipeline/avoid recomputing completed steps----------------------------------------------------------------
# A cached result is reused ONLY if it was produced with the same parameters and the same upstream inputs. 
fingerprint_path <- function(path) paste0(path, ".hash")

read_fingerprint <- function(path) {
  sidecar <- fingerprint_path(path)
  if (file.exists(sidecar)) readLines(sidecar, warn = FALSE)[1] else NA_character_
}

write_fingerprint <- function(path, fingerprint) {
  writeLines(fingerprint, fingerprint_path(path))
}

fingerprint_matches <- function(path, fingerprint) {
  file.exists(path) && identical(read_fingerprint(path), fingerprint)
}

make_fingerprint <- function(params = list(), depends_on = character(0)) {
  upstream <- vapply(depends_on, function(p) {
    fp <- read_fingerprint(p)                      # upstream made by load_or_compute()
    if (is.na(fp) && file.exists(p)) fp <- rlang::hash_file(p)   # plain saveRDS() file: hash its contents
    fp
  }, character(1))
  rlang::hash(list(params = params, upstream = unname(upstream)))
}

load_or_compute <- function(path, compute_fun, params = list(),
                            depends_on = character(0), force = FORCE_RERUN) {
  fingerprint <- make_fingerprint(params, depends_on)

  if (!force && fingerprint_matches(path, fingerprint)) {
    message("  [cached] loading ", path)
    return(readRDS(path))
  }

  if (!force && file.exists(path)) {
    message("  [parameters or inputs changed] recomputing ", path)
  } else {
    message("  [computing] ", path)
  }
  result <- compute_fun()
  saveRDS(result, path)
  write_fingerprint(path, fingerprint)
  result
}

# Set plot format --------------------------------------------------------------
save_plot <- function(plot, filename, width = 8, height = 6) {
  ggplot2::ggsave(file.path(dir_figures, paste0(filename, ".pdf")),
                   plot = plot, width = width, height = height)
  ggplot2::ggsave(file.path(dir_figures, paste0(filename, ".png")),
                   plot = plot, width = width, height = height, dpi = 300)
}

# Extract samples IDs from FASTQ -----------------------------------------------
sample_id_from_filename <- function(filenames) {
  make.unique(sapply(strsplit(basename(filenames), "_"), `[`, 1))
}


# Metadata loading -------------------------------------------------------------
if (!file.exists(metadata_file)) {
  stop("metadata.csv not found: ", metadata_file)
}

metadata <- read.csv(
  metadata_file,
  colClasses = "character",
  fileEncoding = "UTF-8-BOM",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!group_column %in% colnames(metadata)) {
  stop("Column '", group_column, "' was not found in metadata.csv.")
}

if (!reference_group %in% metadata[[group_column]]) {
  stop("Reference group '", reference_group,
       "' was not found in column '", group_column, "'.")
}
