# 0. Parameters
# Edit only this file before running the pipeline.

# Configuration ----------------------------------------------------------------

set.seed(1234)
SEED         <- 1234
N_THREADS    <- TRUE
FORCE_RERUN  <- FALSE

# Input ------------------------------------------------------------------------

# Change this to the path of your project directory containing the FASTQ files, metadata.csv, and Silva/ reference database.
project_dir  <- "C:/Users/Grace_Hoang/OneDrive/Desktop/bioinfor/Samples"

# Metadata path
metadata_file <- file.path(project_dir, "metadata.csv")

# Metadata column containing the groups to compare
# (e.g., "Group", "Treatment", or "Location")
group_column <- "Group"

# Reference level for differential abundance
reference_group <- "Lagoon_Amuricata"


# SILVA database ---------------------------------------------------------------

silva_ref_path <- file.path(
  project_dir,
  "Silva",
  "silva_nr99_v138.2_toGenus_trainset.fa.gz"
)

# Primers ----------------------------------------------------------------------

primer_fwd <- "GGATTAGATACCCTGGTA"
primer_rev <- "CRRCACGAGCTGACGAC"

# DADA2 ------------------------------------------------------------------------

trunc_qmin <- 30

mergepairs_strategy     <- "merge"
mergepairs_min_overlap  <- 12
mergepairs_max_mismatch <- 0

chimera_method <- "consensus"
max_ee         <- c(2, 2)
max_n          <- 0
trunc_q        <- 2
remove_phix    <- TRUE

# ASV length -------------------------------------------------------------------
# Calculate as: amplicon length - forward primer - reverse primer

asv_min_len <- 230
asv_max_len <- 280

# Taxonomic filtering ----------------------------------------------------------

target_kingdoms <- "Bacteria"
remove_unassigned <- TRUE

# Taxa to remove
# Always keep "Mitochondria" and "Chloroplast" in this list.
# If negative controls are available (use_decontam = TRUE),
# remove the laboratory contaminant genera below and let decontam identify them.
contaminant_genera <- c(
  "Mitochondria",
  "Chloroplast",
  "Brevibacterium",
  "Brachybacterium",
  "Dietzia",
  "Thermicanus",
  "Pseudomonas"
)

# Automated contaminant detection (requires negative controls)
use_decontam   <- FALSE
control_column <- "Sample_Type"
control_label  <- "Control"

# Rarefaction ------------------------------------------------------------------
# Choose rarefaction mode:
# 1. "min": Rarefies all samples to the minimum sequencing depth.
# 2. "custom": Rarefies all samples to a user-defined sequencing depth.
# 3. "auto_rarefy": Automatically selects an appropriate sequencing depth based on the data.

rarefaction_mode <- "auto_rarefy"

# If rarefaction_mode = "custom":
custom_rarefaction_depth <- 10000

# If rarefaction_mode = "auto_rarefy":
# Most users should leave these at their default values.

## Maximum sample loss (0.10 = 10%)
max_sample_loss <- 0.10

## Plateau gain threshold 
plateau_gain    <- 0.01

## Minimum samples reaching plateau
plateau_samples <- 0.90

## Plateau evaluation step size (reads)
read_step        <- 1000

# Diversity --------------------------------------------------------------------

alpha_measures <- c("Observed", "Chao1", "Shannon", "Simpson")
beta_distances <- c("bray", "wunifrac", "unifrac")
permanova_permutations <- 999

# Relative abundance -----------------------------------------------------------
# Taxa below these thresholds are grouped as "Other"

min_phylum_abundance <- 0.01
min_genus_abundance  <- 0.02

