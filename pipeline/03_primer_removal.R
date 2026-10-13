# 3. PRIMER REMOVAL

suppressPackageStartupMessages({
  library(tidyverse)
  library(Biostrings)
})

fastq_info <- readRDS(file.path(dir_rds, "fastq_info.rds"))

# Check cutadapt -------------------------------------------------------------

cutadapt_command <- if (exists("CUTADAPT_COMMAND")) CUTADAPT_COMMAND else "cutadapt"
if (!file.exists(cutadapt_command) && Sys.which(cutadapt_command) == "") {
  stop("cutadapt is not installed or not available on PATH.")
}

# Output folder --------------------------------------------------------------

dir_primer_trimmed <- file.path(dir_results, "primer_trimmed")
dir.create(dir_primer_trimmed, showWarnings = FALSE, recursive = TRUE)

primer_removal_log <- file.path(dir_logs, "cutadapt_primer_removal.log")
rds_out <- file.path(dir_rds, "fastq_info_primer_trimmed.rds")

fastq_info <- fastq_info %>%
  mutate(
    fnF.trim = file.path(
      dir_primer_trimmed,
      paste0(SampleID, "_R1_trimmed.fastq.gz")
    ),
    fnR.trim = file.path(
      dir_primer_trimmed,
      paste0(SampleID, "_R2_trimmed.fastq.gz")
    )
  )

# Reverse-complement primers -------------------------------------------------

primer_rc <- function(seq) {
  as.character(reverseComplement(DNAString(seq)))
}

primer_fwd_rc <- primer_rc(primer_fwd)
primer_rev_rc <- primer_rc(primer_rev)

# Check whether primer trimming needs to be rerun -----------------------------

primer_fingerprint <- make_fingerprint(
  params = list(
    primer_fwd = primer_fwd,
    primer_rev = primer_rev,
    input_paths = c(fastq_info$fnF, fastq_info$fnR),
    input_sizes = file.size(c(fastq_info$fnF, fastq_info$fnR))
  ),
  depends_on = file.path(dir_rds, "fastq_info.rds")
)

trimmed_reads_current <-
  all(file.exists(fastq_info$fnF.trim, fastq_info$fnR.trim)) &&
  fingerprint_matches(primer_removal_log, primer_fingerprint)

# Run cutadapt ---------------------------------------------------------------

if (!trimmed_reads_current || FORCE_RERUN) {
  
  if (file.exists(primer_removal_log))
    file.remove(primer_removal_log)
  # A failed forced rerun must not retain the previous success fingerprint.
  if (file.exists(fingerprint_path(primer_removal_log)))
    file.remove(fingerprint_path(primer_removal_log))
  
  message("Running cutadapt on ", nrow(fastq_info), " samples...")
  
  threads_arg <- c(
    "-j",
    if (isTRUE(N_THREADS)) "0" else as.character(N_THREADS)
  )
  
  for (i in seq_len(nrow(fastq_info))) {
    
    args <- c(
      threads_arg,
      "-g", primer_fwd,
      "-G", primer_rev,
      "-a", primer_rev_rc,
      "-A", primer_fwd_rc,
      "--discard-untrimmed",
      "-n", "2",
      "-o", fastq_info$fnF.trim[i],
      "-p", fastq_info$fnR.trim[i],
      fastq_info$fnF[i],
      fastq_info$fnR[i]
    )
    
    result <- system2(
      cutadapt_command,
      vapply(args, shQuote, character(1)),
      stdout = TRUE,
      stderr = TRUE
    )
    
    cat(
      c(
        paste("### Sample:", fastq_info$SampleID[i]),
        result,
        ""
      ),
      file = primer_removal_log,
      sep = "\n",
      append = TRUE
    )
    status <- attr(result, "status")
    if (!is.null(status) && status != 0L) {
      stop("Cutadapt failed for sample ", fastq_info$SampleID[i],
           ". See ", primer_removal_log)
    }
    if (!all(file.exists(c(fastq_info$fnF.trim[i], fastq_info$fnR.trim[i])))) {
      stop("Cutadapt did not create both output files for ", fastq_info$SampleID[i])
    }
  }
  
  write_fingerprint(primer_removal_log, primer_fingerprint)
  
  message("Primer removal completed.")
  
} else {
  
  message("  [cached] primer-trimmed reads are up to date.")
  
}

# Save updated FASTQ information ---------------------------------------------

saveRDS(fastq_info, rds_out)
write_fingerprint(rds_out, primer_fingerprint)
