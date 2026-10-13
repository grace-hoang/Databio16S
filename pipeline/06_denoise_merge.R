# 6. DENOISE + MERGE PAIRS

suppressPackageStartupMessages(library(dada2))

fastq_info <- readRDS(file.path(dir_rds, "fastq_info_filtered.rds"))
err_fwd    <- readRDS(file.path(dir_rds, "err_fwd.rds"))
err_rev    <- readRDS(file.path(dir_rds, "err_rev.rds"))

# Dereplicate reads ------------------------------------------------------------

derep_fwd <- derepFastq(fastq_info$fnF.filt, verbose = FALSE)
derep_rev <- derepFastq(fastq_info$fnR.filt, verbose = FALSE)

names(derep_fwd) <- fastq_info$SampleID
names(derep_rev) <- fastq_info$SampleID

# Denoise ----------------------------------------------------------------------

dada_fwd <- load_or_compute(
  file.path(dir_rds, "dada_fwd.rds"),
  function() {
    dada(
      derep_fwd,
      err = err_fwd,
      multithread = N_THREADS
    )
  },
  depends_on = c(
    file.path(dir_rds, "err_fwd.rds"),
    file.path(dir_rds, "fastq_info_filtered.rds")
  )
)

dada_rev <- load_or_compute(
  file.path(dir_rds, "dada_rev.rds"),
  function() {
    dada(
      derep_rev,
      err = err_rev,
      multithread = N_THREADS
    )
  },
  depends_on = c(
    file.path(dir_rds, "err_rev.rds"),
    file.path(dir_rds, "fastq_info_filtered.rds")
  )
)

# Merge paired-end reads -------------------------------------------------------

merge_pairs_by_strategy <- function(strategy) {
  
  if (strategy == "merge") {
    return(
      mergePairs(
        dada_fwd,
        derep_fwd,
        dada_rev,
        derep_rev,
        minOverlap = mergepairs_min_overlap,
        maxMismatch = mergepairs_max_mismatch,
        verbose = TRUE
      )
    )
  }
  
  if (strategy == "concatenate") {
    return(
      mergePairs(
        dada_fwd,
        derep_fwd,
        dada_rev,
        derep_rev,
        justConcatenate = TRUE,
        verbose = TRUE
      )
    )
  }
  
  stop(
    "Unknown mergepairs_strategy: '",
    strategy,
    "'. Use 'merge' or 'concatenate'."
  )
}

mergers <- load_or_compute(
  file.path(dir_rds, "mergers.rds"),
  function() merge_pairs_by_strategy(mergepairs_strategy),
  params = list(
    strategy     = mergepairs_strategy,
    min_overlap  = mergepairs_min_overlap,
    max_mismatch = mergepairs_max_mismatch
  ),
  depends_on = c(
    file.path(dir_rds, "dada_fwd.rds"),
    file.path(dir_rds, "dada_rev.rds")
  )
)

# Build sequence table ---------------------------------------------------------

seqtab <- load_or_compute(
  file.path(dir_rds, "seqtab.rds"),
  function() makeSequenceTable(mergers),
  depends_on = file.path(dir_rds, "mergers.rds")
)

message(
  "Sequence table dimensions: ",
  paste(dim(seqtab), collapse = " x ")
)