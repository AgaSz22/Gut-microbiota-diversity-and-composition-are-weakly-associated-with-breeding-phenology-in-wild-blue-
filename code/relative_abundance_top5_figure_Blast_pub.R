# =============================================================================
# Script:  phyla_relative_abundance.R
# Date:    28/07/2025
# Author:  Agnes Sz.
#
# Purpose:
#   Visualise the relative abundance of the top 5 bacterial phyla across
#   samples, grouping all remaining phyla as "Other".
#
# Notes:
#   - Chloroplasts (order) and Mitochondria (class) are excluded.
#   - Dataset has been pre-filtered to remove contaminants identified via
#     the decontam package (prevalence threshold = 0.57).
#   - Taxonomy source: BLAST-assigned taxonomy.
# =============================================================================


# ── 1. Libraries ──────────────────────────────────────────────────────────────

library(phyloseq)    # Microbiome data structures
library(ggplot2)     # Data visualisation
library(tidyverse)   # Core data manipulation (includes dplyr, tidyr, forcats)
library(taylor)      # Custom colour palettes


# ── 2. Configuration & Constants ─────────────────────────────────────────────

# File paths
path_asv      <- "../data/ASV_with_taxonomy_p57_Blast.csv"
path_metadata <- "../data/Microbial_metadata_cleaned.csv"

# Number of top phyla to display individually (remainder → "Other")
n_top_phyla <- 5

# Columns to remove after filtering (retain only phylum + sample counts)
cols_to_drop <- c(331, 333, 334, 335, 336, 337)

# Custom colour palette (one colour per top phylum + "Other")
palette_phyla <- c(
  "Bacillota"       = "#1AB788",
  "Bacteroidota"    = "#ffdc5d",
  "Actinomycetota"  = "#074151",
  "Pseudomonadota"  = "#F56A6A",
  "Other"           = "#A3334E",
  "Chlamydiota"     = "#00ffff"
)

# Ordered factor levels for stacked bar consistency
phylum_levels <- c(
  "Chlamydiota", "Other", "Bacteroidota",
  "Actinomycetota", "Bacillota", "Pseudomonadota"
)

# ── 3. Load Data ──────────────────────────────────────────────────────────────

asv_raw      <- read.csv(path_asv)
metadata_raw <- read.csv(path_metadata)


# ── 4. Clean Metadata ─────────────────────────────────────────────────────────

# Replace NAs in Stage and Capture_Date with sensible defaults for
# negative control samples
metadata <- metadata_raw %>%
  mutate(
    Stage        = replace_na(Stage, "Negative"),
    Capture_Date = replace_na(Capture_Date, 70)
  )


# ── 5. Filter ASV Table ───────────────────────────────────────────────────────

# Remove plastid-associated reads (Chloroplasts and Mitochondria) which are
# not of bacterial origin and would inflate non-bacterial relative abundances
asv_filtered <- asv_raw %>%
  filter(
    order != "Chloroplast",
    class != "Mitochondria"
  ) %>%
  select(-all_of(cols_to_drop))  # Drop non-phylum taxonomy columns


# ── 6. Reshape to Long Format ─────────────────────────────────────────────────

# Pivot so each row represents one phylum–sample combination
asv_long <- asv_filtered %>%
  pivot_longer(
    cols      = -phylum,
    names_to  = "Sample_ID",
    values_to = "Count"
  )


# ── 7. Merge with Metadata ────────────────────────────────────────────────────

asv_with_meta <- merge(asv_long, metadata, by = "Sample_ID")


# ── 8. Compute Relative Abundance ────────────────────────────────────────────

# 8a. Aggregate counts per phylum × sample (handles duplicate phylum rows)
asv_aggregated <- asv_with_meta %>%
  group_by(across(-Count)) %>%
  summarise(total_count = sum(Count), .groups = "drop")

# 8b. Calculate relative abundance within each sample.
#     Samples with zero total reads are assigned 0 to avoid division by zero.
asv_relabu <- asv_aggregated %>%
  group_by(Sample_ID) %>%
  mutate(
    sample_total  = sum(total_count, na.rm = TRUE),
    rel_abundance = if_else(sample_total == 0, 0, total_count / sample_total)
  ) %>%
  select(-sample_total) %>%
  ungroup()


# ── 9. Identify Top N Phyla ───────────────────────────────────────────────────

# Rank phyla by summed absolute abundance across all samples
top_phyla <- asv_relabu %>%
  group_by(phylum) %>%
  summarise(total_abundance = sum(total_count), .groups = "drop") %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = n_top_phyla) %>%
  pull(phylum)


# ── 10. Group Remaining Phyla as "Other" ──────────────────────────────────────

asv_relabu <- asv_relabu %>%
  mutate(
    phylum_grouped = if_else(phylum %in% top_phyla, phylum, "Other")
  )


# ── 11. Aggregate by Grouped Phylum ──────────────────────────────────────────

# Sum relative abundances within each Sample_ID × phylum_grouped combination
plot_data <- asv_relabu %>%
  group_by(Sample_ID, phylum_grouped) %>%
  summarise(rel_abundance = sum(rel_abundance), .groups = "drop") %>%
  merge(metadata, by = "Sample_ID")  # Re-attach metadata for faceting


# ── 12. Visualisation ─────────────────────────────────────────────────────────

# ── 12a. Ordered by Pseudomonadota abundance ─────────────────────────────────

# Determine sample order within each Stage based on Pseudomonadota abundance
pseudomonadota_order <- plot_data %>%
  filter(phylum_grouped == "Pseudomonadota") %>%
  group_by(Stage, Sample_ID) %>%
  summarise(pseudo_ra = sum(rel_abundance), .groups = "drop") %>%
  arrange(Stage, pseudo_ra) %>%
  group_by(Stage) %>%
  mutate(sample_id_ordered = factor(Sample_ID, levels = Sample_ID)) %>%
  ungroup()

# Join ordering back to the main plot dataset
plot_data_ordered <- plot_data %>%
  left_join(
    pseudomonadota_order %>% select(Sample_ID, Stage, sample_id_ordered),
    by = c("Sample_ID", "Stage")
  )

viz_ordered <- plot_data_ordered %>%
  ggplot(aes(x = sample_id_ordered, y = rel_abundance, fill = phylum_grouped)) +
  geom_bar(
    stat     = "identity",
    position = "fill",
    width    = 1,
    na.rm    = TRUE   # NAs expected for negative controls (no Pseudomonadota row)
  ) +
  facet_grid(~ Stage, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = palette_phyla) +
  theme_minimal() +
  labs(
    x     = "Sample ID",
    y     = "Relative abundance",
    fill  = "Phylum"
  ) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  )


# ── 12b. Aggregated barplot by Stage (excluding negative controls) ────────────

plot_data_stages <- plot_data %>%
  filter(Stage != "Negative") %>%
  group_by(Stage, phylum_grouped) %>%
  summarise(rel_abundance = mean(rel_abundance), .groups = "drop") %>%  # or sum()
  mutate(phylum_grouped = factor(phylum_grouped, levels = phylum_levels))

viz_by_stage <- ggplot(
  data = plot_data_stages,
  aes(x = as.factor(Stage), y = rel_abundance, fill = phylum_grouped)
) +
  geom_bar(stat = "identity", position = "fill",  alpha = 1) +
  scale_fill_manual(values = palette_phyla) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal(base_size = 15) +
  labs(
    x     = "Stage",
    y     = "Relative Abundance",
    fill  = "Phylum"
  ) +
  theme(legend.position = "none")

# Save the combined figure as an SVG for publication-quality vector output
ggsave(
  filename = "../figures/bp_phyla_blast.svg",
  plot     = viz_by_stage,
  device   = "svg",
  width    = 5,
  height   = 4,
  units    = "in"
)
