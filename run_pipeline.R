# Called once from Python; all steps share the same R environment.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Usage: Rscript run_pipeline.R config.json [--check]")
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Missing jsonlite. Run install_packages.R in the selected R installation.")
}
cfg <- jsonlite::fromJSON(args[1], simplifyVector = TRUE)
step_files <- sort(list.files(cfg$pipeline_dir, pattern = "^[0-9]{2}_.*\\.R$",
                              full.names = TRUE))
if (length(step_files) != 17L ||
    !identical(substr(basename(step_files), 1, 2), sprintf("%02d", 0:16))) {
  stop("The pipeline directory must contain exactly one script per step, 00 to 16.")
}
source(step_files[1], local = .GlobalEnv, encoding = "UTF-8")
project_dir <- cfg$project_dir
metadata_file <- cfg$metadata_file
silva_ref_path <- cfg$silva_ref_path
group_column <- cfg$group_column
reference_group <- cfg$reference_group
primer_fwd <- cfg$primer_fwd
primer_rev <- cfg$primer_rev
asv_min_len <- cfg$asv_min_len
asv_max_len <- cfg$asv_max_len
trunc_qmin <- cfg$trunc_qmin
N_THREADS <- as.integer(cfg$threads)
FORCE_RERUN <- FALSE
rarefaction_mode <- cfg$rarefaction_mode
custom_rarefaction_depth <- cfg$custom_rarefaction_depth
CUTADAPT_COMMAND <- cfg$cutadapt
pipeline_input_files <- unlist(cfg$fastq_files, use.names = FALSE)
set.seed(SEED)

required <- c("tidyverse", "dada2", "ShortRead", "Biostrings", "phyloseq",
              "DECIPHER", "phangorn", "vegan", "Polychrome", "rlang")
if (isTRUE(cfg$run_ancombc)) required <- c(required, "ANCOMBC")
if (isTRUE(use_decontam)) required <- c(required, "decontam")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
message("R: ", R.version.string)
message("R libraries: ", paste(.libPaths(), collapse = "; "))
if (length(missing)) {
  stop("Missing R packages: ", paste(missing, collapse = ", "),
       ". Run install_packages.R in this R installation.")
}
if (!all(c("hash", "hash_file") %in% getNamespaceExports("rlang"))) {
  stop("Update rlang: this pipeline needs rlang::hash and rlang::hash_file.")
}
if (!file.exists(CUTADAPT_COMMAND) && Sys.which(CUTADAPT_COMMAND) == "") {
  stop("Cutadapt not found. Select its executable in the Python window.")
}
version <- system2(CUTADAPT_COMMAND, "--version", stdout = TRUE, stderr = TRUE)
if (!is.null(attr(version, "status")) && attr(version, "status") != 0L) {
  stop("Cutadapt version check failed: ", paste(version, collapse = "\n"))
}
message("Cutadapt: ", paste(version, collapse = " "))
if (!file.exists(metadata_file) || !file.exists(silva_ref_path) ||
    !length(pipeline_input_files) || !all(file.exists(pipeline_input_files))) {
  stop("A selected input file is missing.")
}
meta <- read.csv(metadata_file, colClasses = "character", check.names = FALSE,
                 stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
if (!group_column %in% names(meta)) stop("Group column missing in metadata.")
if (anyNA(meta[[1]]) || any(!nzchar(meta[[1]])) || anyDuplicated(meta[[1]])) {
  stop("The first metadata column must contain unique, nonempty SampleIDs.")
}
fwd_pattern <- "(_R1|_1)(_001)?\\.(fastq|fq)(\\.gz)?$"
rev_pattern <- "(_R2|_2)(_001)?\\.(fastq|fq)(\\.gz)?$"
fwd <- sort(pipeline_input_files[grepl(fwd_pattern, basename(pipeline_input_files))])
rev <- sort(pipeline_input_files[grepl(rev_pattern, basename(pipeline_input_files))])
sample_ids <- sub(fwd_pattern, "", basename(fwd))
if (!length(fwd) || anyDuplicated(sample_ids) ||
    !identical(sample_ids, sub(rev_pattern, "", basename(rev)))) {
  stop("FASTQ pairs must match, with unique sample names; use the same input folder.")
}
idx <- match(sample_ids, meta[[1]])
if (anyNA(idx)) stop("FASTQ SampleIDs missing from the first metadata column.")
groups <- meta[[group_column]][idx]
if (anyNA(groups) || any(!nzchar(trimws(groups))) ||
    length(unique(groups)) < 2L || any(table(groups) < 2L)) {
  stop("Full diversity comparisons need at least two groups and two samples per group.")
}
if (!reference_group %in% groups) stop("Reference group is absent from selected samples.")
message("Inputs validated: ", length(sample_ids), " paired-end samples.")
if ("--check" %in% args) {
  message("CHECK_OK")
  quit(save = "no", status = 0)
}
setwd(cfg$run_dir)
started <- Sys.time()
for (i in 2:length(step_files)) {
  if (i == 17L && !isTRUE(cfg$run_ancombc)) next
  message("\n[STEP ", basename(step_files[i]), "] ", format(Sys.time()))
  source(step_files[i], local = .GlobalEnv, encoding = "UTF-8")
  if (i == 11L && (phyloseq::ntaxa(physeq) < 3L ||
                    any(phyloseq::sample_sums(physeq) == 0))) {
    stop("Too few ASVs for a tree or a sample has no reads after taxonomic filtering. ",
         "Review primers, filtering and selected reference database.")
  }
}
capture.output(sessionInfo(), file = file.path(dir_logs, "sessionInfo.txt"))
message("\nPIPELINE_OK: ", normalizePath(dir_results, winslash = "/"))
message("Elapsed: ", round(as.numeric(difftime(Sys.time(), started, units = "mins")), 2),
        " minutes")
