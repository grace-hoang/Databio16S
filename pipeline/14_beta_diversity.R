# 14. BETA DIVERSITY

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(tidyverse)
})

physeq_rarefied <- readRDS(file.path(dir_rds, "physeq_rarefied.rds"))

# Extract metadata from phyloseq -----------------------------------------------
metadata_df <- as(sample_data(physeq_rarefied), "data.frame")
metadata_df$SampleID <- rownames(metadata_df)

# Remove samples with missing group information
metadata_df <- metadata_df %>%
  filter(!is.na(.data[[group_column]]))

# Keep phyloseq in the same sample order
physeq_rarefied <- prune_samples(metadata_df$SampleID, physeq_rarefied)
metadata_df <- metadata_df[match(sample_names(physeq_rarefied), metadata_df$SampleID), ]

metadata_df[[group_column]] <- droplevels(factor(metadata_df[[group_column]]))

distance_labels <- c(
  bray      = "Bray-Curtis",
  wunifrac  = "Weighted UniFrac",
  unifrac   = "Unweighted UniFrac"
)

permanova_results  <- list()
betadisper_results <- list()

for (dist_method in beta_distances) {
  
  dist_matrix <- load_or_compute(
    file.path(dir_rds, paste0("dist_", dist_method, ".rds")),
    function() {
      phyloseq::distance(physeq_rarefied, method = dist_method)
    },
    params = list(method = dist_method),
    depends_on = file.path(dir_rds, "physeq_rarefied.rds")
  )
  
  # PCoA -----------------------------------------------------------------------
  
  ordination <- ordinate(
    physeq_rarefied,
    method = "PCoA",
    distance = dist_matrix
  )
  
  pcoa_plot <- plot_ordination(
    physeq_rarefied,
    ordination,
    color = group_column
  ) +
    stat_ellipse(type = "norm", linetype = 2) +
    labs(title = paste("PCoA -", distance_labels[[dist_method]])) +
    theme_bw(base_size = 13)
  
  save_plot(pcoa_plot, paste0("pcoa_", dist_method))
  
  # PERMANOVA ------------------------------------------------------------------
  
  permanova_fit <- adonis2(
    as.formula(paste("dist_matrix ~", group_column)),
    data = metadata_df,
    permutations = permanova_permutations
  )
  
  permanova_results[[dist_method]] <-
    as.data.frame(permanova_fit) %>%
    rownames_to_column("term") %>%
    mutate(distance = distance_labels[[dist_method]], .before = 1)
  
  # Test homogeneity of dispersion ---------------------------------------------
  
  disper_fit <- betadisper(
    dist_matrix,
    metadata_df[[group_column]]
  )
  
  disper_test <- permutest(
    disper_fit,
    permutations = permanova_permutations
  )
  
  betadisper_results[[dist_method]] <-
    as.data.frame(disper_test$tab) %>%
    rownames_to_column("term") %>%
    mutate(distance = distance_labels[[dist_method]], .before = 1)
}

write.csv(
  bind_rows(permanova_results),
  file.path(dir_tables, "permanova_results.csv"),
  row.names = FALSE
)

write.csv(
  bind_rows(betadisper_results),
  file.path(dir_tables, "betadisper_results.csv"),
  row.names = FALSE
)