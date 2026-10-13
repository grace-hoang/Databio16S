# 2. QUALITY PROFILES ----------------------------------------------------------

suppressPackageStartupMessages({
  library(dada2)
  library(tidyverse)
})

fastq_info <- load_or_compute(
  file.path(dir_rds, "fastq_info.rds"),
  function() {
    
    # Find paired-end FASTQ files -----------------------------------------------
    
    pattern_F <- "(_R1|_1)(_001)?\\.(fastq|fq)(\\.gz)?$"
    pattern_R <- "(_R2|_2)(_001)?\\.(fastq|fq)(\\.gz)?$"
    
    # Python supplies exact selected files; never scan pipeline output.
    input_files <- if (exists("pipeline_input_files", inherits = TRUE)) {
      pipeline_input_files
    } else {
      list.files(project_dir, full.names = TRUE, recursive = TRUE)
    }
    fnF <- sort(input_files[grepl(pattern_F, basename(input_files))])
    fnR <- sort(input_files[grepl(pattern_R, basename(input_files))])
    
    sample_names <- sub(pattern_F, "", basename(fnF))
    
    if (!length(fnF) ||
        !identical(sample_names, sub(pattern_R, "", basename(fnR)))) {
      stop("FASTQ files not found or forward/reverse pairs do not match.")
    }
    
    # Load metadata -------------------------------------------------------------
    
    meta <- metadata
    colnames(meta)[1] <- "SampleID"
    
    if (anyDuplicated(meta$SampleID))
      stop("Duplicate SampleIDs found in metadata.csv.")
    
    missing <- setdiff(sample_names, meta$SampleID)
    
    if (length(missing))
      stop(
        "The following FASTQ samples are missing from metadata.csv:\n",
        paste(missing, collapse = ", ")
      )
    
    # Combine FASTQ paths with metadata -----------------------------------------
    
    fastq_info <- data.frame(
      SampleID = sample_names,
      fnF = fnF,
      fnR = fnR,
      stringsAsFactors = FALSE
    ) %>%
      inner_join(meta, by = "SampleID")
    
    message(
      "Successfully matched ",
      nrow(fastq_info),
      " samples with metadata.csv."
    )
    
    # Generate quality profile plots (first sample only) ------------------------
    
    save_plot(
      plotQualityProfile(fastq_info$fnF[1]) +
        ggtitle(paste("Forward:", fastq_info$SampleID[1])),
      "quality_profile_forward"
    )
    save_plot(
      plotQualityProfile(fastq_info$fnR[1]) +
        ggtitle(paste("Reverse:", fastq_info$SampleID[1])),
      "quality_profile_reverse"
    )
    fastq_info
  },
  params = list(
    metadata = metadata,
    project_dir = project_dir
  ),
  depends_on = c(
    metadata_file,
    if (exists("pipeline_input_files", inherits = TRUE)) pipeline_input_files else
      list.files(project_dir, pattern = "\\.(fastq|fq)(\\.gz)?$",
                 full.names = TRUE, recursive = TRUE)
  )
)
