# 7. CHIMERA REMOVAL + ASV LENGTH FILTER

suppressPackageStartupMessages(library(dada2))

seqtab <- readRDS(file.path(dir_rds, "seqtab.rds"))

# Remove chimeras --------------------------------------------------------------

seqtab_nochim <- load_or_compute(
  file.path(dir_rds, "seqtab_nochim.rds"),
  function() {
    removeBimeraDenovo(
      seqtab,
      method = chimera_method,
      multithread = N_THREADS,
      verbose = TRUE
    )
  },
  params = list(
    method = chimera_method
  ),
  depends_on = file.path(dir_rds, "seqtab.rds")
)

message(
  "Fraction of reads kept after chimera removal: ",
  round(sum(seqtab_nochim) / sum(seqtab), 4)
)

# Filter ASVs by length --------------------------------------------------------

asv_lengths <- nchar(getSequences(seqtab_nochim))

write.csv(
  table(asv_lengths),
  file.path(dir_tables, "asv_length_distribution.csv"),
  row.names = FALSE
)

seqtab_final <- load_or_compute(
  file.path(dir_rds, "seqtab_final.rds"),
  function() {
    
    keep_by_length <- asv_lengths >= asv_min_len &
      asv_lengths <= asv_max_len
    
    message(
      sum(keep_by_length), " of ",
      length(keep_by_length),
      " ASVs kept within the ",
      asv_min_len, "-", asv_max_len, " bp window"
    )
    
    seqtab_nochim[, keep_by_length, drop = FALSE]
    
  },
  params = list(
    min_len = asv_min_len,
    max_len = asv_max_len
  ),
  depends_on = file.path(dir_rds, "seqtab_nochim.rds")
)