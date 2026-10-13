# 16. ANCOM-BC

suppressPackageStartupMessages({
  library(phyloseq)
  library(ANCOMBC)
  library(tidyverse)
})

physeq <- readRDS(file.path(dir_rds, "physeq.rds"))

# Check metadata ---------------------------------------------------------------
meta <- as(sample_data(physeq), "data.frame")

if (!(group_column %in% colnames(meta))) {
  stop("Column '", group_column, "' not found in metadata.csv.")
}

if (!(reference_group %in% meta[[group_column]])) {
  stop("Reference group '", reference_group,
       "' not found in metadata column '", group_column, "'.")
}

sample_data(physeq)[[group_column]] <-
  relevel(factor(sample_data(physeq)[[group_column]]),
          ref = reference_group)

# Run ANCOM-BC -----------------------------------------------------------------
ancombc_settings <- list(
  tax_level    = "Genus",
  p_adj_method = "holm",
  struc_zero   = TRUE,
  neg_lb       = TRUE,
  alpha        = 0.05,
  global       = FALSE
)

ancombc_fit <- load_or_compute(
  file.path(dir_rds, "ancombc_fit.rds"),
  function() {
    ancombc2(
      data = physeq,
      tax_level = ancombc_settings$tax_level,
      fix_formula = group_column,
      rand_formula = NULL,
      p_adj_method = ancombc_settings$p_adj_method,
      group = group_column,
      struc_zero = ancombc_settings$struc_zero,
      neg_lb = ancombc_settings$neg_lb,
      alpha = ancombc_settings$alpha,
      global = ancombc_settings$global
    )
  },
  params = c(
    list(
      group_column = group_column,
      reference_group = reference_group
    ),
    ancombc_settings
  ),
  depends_on = file.path(dir_rds, "physeq.rds")
)

write.csv(
  ancombc_fit$res,
  file.path(dir_tables, "ancombc_results.csv"),
  row.names = FALSE
)

# Volcano plot -----------------------------------------------------------------
lfc_col <- grep("^lfc_", names(ancombc_fit$res), value = TRUE)[1]
q_col   <- grep("^q_", names(ancombc_fit$res), value = TRUE)[1]

if (!is.na(lfc_col) && !is.na(q_col)) {
  
  volcano_data <- ancombc_fit$res %>%
    transmute(
      taxon,
      lfc = .data[[lfc_col]],
      q = .data[[q_col]],
      significant = q < 0.05
    )
  
  ancombc_plot <- ggplot(volcano_data,
                         aes(x = lfc,
                             y = -log10(q),
                             color = significant)) +
    geom_point(alpha = 0.7) +
    labs(
      x = paste0("Log fold change (", group_column,
                 ", ref = ", reference_group, ")"),
      y = expression(-log[10](q)),
      color = "q < 0.05"
    ) +
    theme_bw(base_size = 13)
  
  save_plot(ancombc_plot, "ancombc_volcano")
}