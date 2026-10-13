# 10. BUILD THE PHYLOSEQ OBJECT
# Combines ASV table, taxonomy table, and metadata

suppressPackageStartupMessages({
  library(phyloseq)
  library(Biostrings)
})

seqtab_final <- readRDS(file.path(dir_rds, "seqtab_final.rds"))
taxonomy     <- readRDS(file.path(dir_rds, "taxonomy.rds"))
fastq_info   <- readRDS(file.path(dir_rds, "fastq_info.rds"))

physeq <- load_or_compute(
  file.path(dir_rds, "physeq.rds"),
  function() {
    
# Build phyloseq object -----------------------------------------------
    
    sample_data_df <- fastq_info
    rownames(sample_data_df) <- sample_data_df$SampleID
    
    physeq_raw <- phyloseq(
      otu_table(seqtab_final, taxa_are_rows = FALSE),
      tax_table(taxonomy),
      sample_data(sample_data_df)
    )
    
    ## Preserve original ASV sequences
    asv_sequences <- DNAStringSet(taxa_names(physeq_raw))
    names(asv_sequences) <- paste0("ASV", seq_along(asv_sequences))
    
    taxa_names(physeq_raw) <- names(asv_sequences)
    physeq_raw <- merge_phyloseq(physeq_raw, asv_sequences)
    
    writeXStringSet(
      asv_sequences,
      file.path(dir_tables, "asv_sequences.fasta")
    )
    
# Taxonomic filtering -------------------------------------------------
    
    tax_df <- as.data.frame(tax_table(physeq_raw))
    
    keep <- rep(TRUE, nrow(tax_df))
    
    if (!is.null(target_kingdoms)) {
      keep <- keep & tax_df$Kingdom %in% target_kingdoms
    }
    
    if (remove_unassigned) {
      keep <- keep & !is.na(tax_df$Kingdom)
    }
    
    keep <- keep &
      !(tax_df$Genus %in% contaminant_genera)
    
    physeq <- prune_taxa(keep, physeq_raw)
    
    message(
      sprintf(
        "Taxonomic filtering retained %d of %d taxa.",
        ntaxa(physeq),
        ntaxa(physeq_raw)
      )
    )
    
# Optional decontam --------------------------------------------------------
    
    if (use_decontam &&
        control_column %in% colnames(sample_data_df)) {
      
      message(">> Running automated decontamination...")
      
      library(decontam)
      
      is_control <- sample_data_df[[control_column]] == control_label
      
      if (sum(is_control) > 0) {
        
        contam_df <- isContaminant(
          seqtab_final,
          neg = is_control,
          method = "prevalence",
          threshold = 0.5
        )
        
        contam_asvs <- rownames(contam_df)[contam_df$isContaminant]
        
        message(
          "Removed ",
          length(contam_asvs),
          " contaminant ASVs."
        )
        write.csv(
          contam_df[contam_df$isContaminant, ],
          file.path(dir_tables, "decontam_removed_asvs.csv"),
          row.names = TRUE
        )
        physeq <- prune_taxa(
          !taxa_names(physeq) %in% contam_asvs,
          physeq
        )
        physeq <- subset_samples(
          physeq,
          !is_control
        )
      } else {
        warning(
          "No samples matched control_label = '",
          control_label,
          "'."
        )
      }
    }
    physeq
  },
  params = list(
    remove_unassigned = remove_unassigned,
    target_kingdoms   = target_kingdoms,
    contaminant_genera = contaminant_genera,
    use_decontam = use_decontam,
    control_column = control_column,
    control_label = control_label
  ),
  depends_on = c(
    file.path(dir_rds, "seqtab_final.rds"),
    file.path(dir_rds, "taxonomy.rds"),
    file.path(dir_rds, "fastq_info.rds")
  )
)