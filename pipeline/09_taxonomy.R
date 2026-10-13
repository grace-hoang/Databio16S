# 9. TAXONOMY ASSIGNMENT

suppressPackageStartupMessages(library(dada2))

if (!file.exists(silva_ref_path)) {
  stop(
    "SILVA reference not found at: ",
    silva_ref_path,
    ". Check `silva_ref_path` in 00_parameters.R."
  )
}

seqtab_final <- readRDS(file.path(dir_rds, "seqtab_final.rds"))

taxonomy_ranks <- c(
  "Kingdom",
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus"
)

taxonomy <- load_or_compute(
  file.path(dir_rds, "taxonomy.rds"),
  function() {
    
    assignTaxonomy(
      seqtab_final,
      refFasta    = silva_ref_path,
      taxLevels   = taxonomy_ranks,
      tryRC       = TRUE,
      multithread = N_THREADS
    )
    
  },
  params = list(
    ranks    = taxonomy_ranks,
    ref_path = normalizePath(silva_ref_path),
    ref_size = file.size(silva_ref_path),
    tryRC    = TRUE
  ),
  depends_on = file.path(dir_rds, "seqtab_final.rds")
)

taxonomy_out <- taxonomy
rownames(taxonomy_out) <- NULL

write.csv(
  taxonomy_out,
  file.path(dir_tables, "taxonomy.csv"),
  row.names = FALSE
)