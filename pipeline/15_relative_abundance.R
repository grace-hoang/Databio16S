# 15. RELATIVE ABUNDANCE
# Computed on the original (non-rarefied) phyloseq object

suppressPackageStartupMessages({
  library(phyloseq)
  library(tidyverse)
  library(Polychrome)
})

physeq <- readRDS(file.path(dir_rds, "physeq.rds"))

# Check metadata ---------------------------------------------------------------
metadata_df <- as(sample_data(physeq), "data.frame")

if (!(group_column %in% colnames(metadata_df))) {
  stop(
    "Column '", group_column,
    "' was not found in metadata.csv.\nAvailable columns: ",
    paste(colnames(metadata_df), collapse = ", ")
  )
}

# Helper -----------------------------------------------------------------------
format_other_label <- function(cutoff) {
  pct <- if (cutoff <= 1) cutoff * 100 else cutoff
  paste0("Other (<", pct, "%)")
}

other_phylum_label <- format_other_label(min_phylum_abundance)
other_genus_label  <- format_other_label(min_genus_abundance)

# Phylum -----------------------------------------------------------------------
phylum_abundance <-
  physeq %>%
  tax_glom(taxrank = "Phylum", NArm = TRUE) %>%
  transform_sample_counts(function(x) x / sum(x)) %>%
  psmelt() %>%
  mutate(
    Phylum = ifelse(
      Abundance < min_phylum_abundance,
      other_phylum_label,
      as.character(Phylum)
    )
  ) %>%
  group_by(Sample, .data[[group_column]], Phylum) %>%
  summarise(
    Abundance = sum(Abundance),
    .groups = "drop"
  )

# Genus ------------------------------------------------------------------------
genus_abundance <-
  physeq %>%
  tax_glom(taxrank = "Genus", NArm = TRUE) %>%
  transform_sample_counts(function(x) x / sum(x)) %>%
  psmelt() %>%
  mutate(
    Genus = ifelse(
      Abundance < min_genus_abundance,
      other_genus_label,
      as.character(Genus)
    )
  ) %>%
  group_by(Sample, .data[[group_column]], Genus) %>%
  summarise(
    Abundance = sum(Abundance),
    .groups = "drop"
  )

# Save tables ------------------------------------------------------------------
write.csv(
  phylum_abundance,
  file.path(dir_tables, "relative_abundance_phylum.csv"),
  row.names = FALSE
)

write.csv(
  genus_abundance,
  file.path(dir_tables, "relative_abundance_genus.csv"),
  row.names = FALSE
)

# Color palettes ---------------------------------------------------------------
build_taxon_palette <- function(taxon_names, other_label, other_color = "#D3D3D3") {
  
  unique_taxa <- unique(as.character(taxon_names))
  
  if (other_label %in% unique_taxa) {
    
    regular_taxa <- sort(setdiff(unique_taxa, other_label))
    ordered_taxa <- c(regular_taxa, other_label)
    
    set.seed(SEED)
    
    colors <- Polychrome::createPalette(
      max(3, length(regular_taxa)),
      seedcolors = c("#EE3333", "#3377EE", "#33AA33"),
      range = c(20, 60)
    )[1:length(regular_taxa)]
    
    names(colors) <- regular_taxa
    
    palette <- c(colors, setNames(other_color, other_label))
    
  } else {
    
    ordered_taxa <- sort(unique_taxa)
    
    set.seed(SEED)
    
    palette <- Polychrome::createPalette(
      max(3, length(ordered_taxa)),
      seedcolors = c("#EE3333", "#3377EE", "#33AA33"),
      range = c(15,85)
    )[1:length(ordered_taxa)]
    
    names(palette) <- ordered_taxa
  }
  palette
}

phylum_palette <- build_taxon_palette(phylum_abundance$Phylum, other_phylum_label)
genus_palette  <- build_taxon_palette(genus_abundance$Genus, other_genus_label)

phylum_abundance$Phylum <- factor(phylum_abundance$Phylum,
                                  levels = names(phylum_palette))

genus_abundance$Genus <- factor(genus_abundance$Genus,
                                levels = names(genus_palette))

# Plot function ----------------------------------------------------------------
abundance_barplot <- function(data, taxrank, palette) {
  
  legend_ncol <- ceiling(length(palette) / 25)
  
  ggplot(
    data,
    aes(
      x = Sample,
      y = Abundance,
      fill = .data[[taxrank]]
    )
  ) +
    geom_bar(stat = "identity", width = 0.9) +
    facet_wrap(
      as.formula(paste0("~", group_column)),
      nrow = 1,
      scales = "free_x"
    ) +
    scale_fill_manual(
      values = palette,
      breaks = names(palette),
      name = taxrank
    ) +
    scale_y_continuous(
      labels = scales::percent,
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      x = "Sample",
      y = "Relative abundance"
    ) +
    guides(fill = guide_legend(ncol = legend_ncol, byrow = FALSE)) +
    theme_bw(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.4, "cm"),
      plot.margin = margin(10, 20, 10, 10)
    )
}

n_samples <- n_distinct(phylum_abundance$Sample)

# Save plots -------------------------------------------------------------------
phylum_width <- max(
  10,
  n_samples * 0.5 + ceiling(length(phylum_palette) / 25) * 2.5
)

save_plot(
  abundance_barplot(phylum_abundance, "Phylum", phylum_palette),
  "relative_abundance_phylum",
  width = phylum_width,
  height = 6.5
)

genus_width <- max(
  11,
  n_samples * 0.5 + ceiling(length(genus_palette) / 25) * 2.8
)

genus_height <- max(
  6.5,
  length(genus_palette) * 0.02
)

save_plot(
  abundance_barplot(genus_abundance, "Genus", genus_palette),
  "relative_abundance_genus",
  width = genus_width,
  height = genus_height
)