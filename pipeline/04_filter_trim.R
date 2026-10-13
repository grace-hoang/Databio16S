# 4. FILTER AND TRIM

suppressPackageStartupMessages({
  library(dada2)
  library(ShortRead)
})

fastq_info <- readRDS(file.path(dir_rds, "fastq_info_primer_trimmed.rds"))

# Determine truncation length --------------------------------------------------

determine_trunc_length <- function(
    fastq_files,
    trunc_qmin,
    min_length_fraction = 0.90,
    safety_trim_bp = 0
) {
  
  quality_list <- lapply(fastq_files, function(f) {
    reads <- readFastq(f)
    as(quality(reads), "matrix")
  })
  
  max_cycle <- max(vapply(quality_list, ncol, integer(1)))
  
  quality_list <- lapply(quality_list, function(m) {
    if (ncol(m) < max_cycle) {
      m <- cbind(
        m,
        matrix(NA_real_,
               nrow = nrow(m),
               ncol = max_cycle - ncol(m))
      )
    }
    m
  })
  
  quality_matrix <- do.call(rbind, quality_list)
  
  read_fraction <- colMeans(!is.na(quality_matrix))
  
  supported_cycles <- which(read_fraction >= min_length_fraction)
  
  if (!length(supported_cycles)) {
    stop(
      "No sequencing cycle is reached by ",
      min_length_fraction * 100,
      "% of reads."
    )
  }
  
  max_supported_cycle <- max(supported_cycles)
  
  median_quality <- apply(
    quality_matrix[, seq_len(max_supported_cycle), drop = FALSE],
    2,
    median,
    na.rm = TRUE
  )
  
  first_drop <- which(median_quality < trunc_qmin)[1]
  
  if (is.na(first_drop)) {
    trunc_length <- max_supported_cycle
    if (safety_trim_bp > 0)
      trunc_length <- max(1, trunc_length - safety_trim_bp)
  } else {
    trunc_length <- max(1, first_drop - 1)
  }
  
  trunc_length
}

# Automatically choose truncation lengths -------------------------------------

trunc_len_fwd <- determine_trunc_length(
  fastq_info$fnF.trim,
  trunc_qmin,
  safety_trim_bp = 10
)

trunc_len_rev <- determine_trunc_length(
  fastq_info$fnR.trim,
  trunc_qmin
)

message(
  "Automatically chosen truncLen: forward = ",
  trunc_len_fwd,
  ", reverse = ",
  trunc_len_rev
)

writeLines(
  c(
    paste("trunc_len_fwd:", trunc_len_fwd),
    paste("trunc_len_rev:", trunc_len_rev)
  ),
  file.path(dir_logs, "trunc_len_chosen.txt")
)

# Filter and trim --------------------------------------------------------------

dir_filtered <- file.path(dir_results, "filtered")
dir.create(dir_filtered, showWarnings = FALSE)

fastq_info <- fastq_info %>%
  mutate(
    fnF.filt = file.path(dir_filtered, paste0(SampleID, "_F_filt.fastq.gz")),
    fnR.filt = file.path(dir_filtered, paste0(SampleID, "_R_filt.fastq.gz"))
  )

saveRDS(
  fastq_info,
  file.path(dir_rds, "fastq_info_filtered.rds")
)

filter_summary <- load_or_compute(
  file.path(dir_rds, "filter_summary.rds"),
  function() {
    
    filterAndTrim(
      fwd         = fastq_info$fnF.trim,
      filt        = fastq_info$fnF.filt,
      rev         = fastq_info$fnR.trim,
      filt.rev    = fastq_info$fnR.filt,
      truncLen    = c(trunc_len_fwd, trunc_len_rev),
      maxN        = max_n,
      maxEE       = max_ee,
      truncQ      = trunc_q,
      rm.phix     = remove_phix,
      compress    = TRUE,
      multithread = N_THREADS
    )
    
  },
  params = list(
    trunc_len_fwd = trunc_len_fwd,
    trunc_len_rev = trunc_len_rev,
    max_n         = max_n,
    max_ee        = max_ee,
    trunc_q       = trunc_q,
    remove_phix   = remove_phix
  ),
  depends_on = file.path(dir_rds, "fastq_info_primer_trimmed.rds")
)

rownames(filter_summary) <- fastq_info$SampleID

write.csv(
  filter_summary,
  file.path(dir_tables, "filter_summary.csv")
)