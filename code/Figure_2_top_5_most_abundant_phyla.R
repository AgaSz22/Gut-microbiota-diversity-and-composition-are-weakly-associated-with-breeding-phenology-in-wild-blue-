# =============================================================================
# Date:    09/07/2025
# Author:  Agnes Sz.
# Purpose: Visualise the relative abundance of the top 5 phyla across samples,
#          grouping all remaining phyla under "Other".
#          - All bacteria included; chloroplasts and mitochondria excluded.
#          - Dataset has been pre-filtered using decontam (prevalence threshold: 0.57).
# =============================================================================

#Libraries: 

library(tidyverse)


# ── Load Data ─────────────────────────────────────────────────────────────────

asv_raw <- read.csv("../data/ASV_with_taxonomy_p57.csv")
metadata_raw <- read.csv("../data/Microbial_metadata_cleaned.csv")


# ── Clean Metadata ────────────────────────────────────────────────────────────

metadata <- metadata_raw %>%
  mutate(
    Stage = replace_na(Stage, "Negative"),
    Capture_Date = replace_na(Capture_Date, 70) #dummy data for visualisation purposes
  )


# ── Filter ASV Table ──────────────────────────────────────────────────────────

# Retain only the Phylum column and sample count columns.

asv_filtered <- asv_raw %>%
  select(-c(331, 333, 334, 335, 336, 337)) #removing columns with other taxonomic groups


# ── Reshape to Long Format ────────────────────────────────────────────────────

asv_long <- asv_filtered %>%
  pivot_longer(
    cols = -Phylum,
    names_to = "Sample_ID",
    values_to = "Count"
  )


# ── Merge with Metadata ───────────────────────────────────────────────────────

asv_with_metadata <- merge(asv_long, metadata, by = "Sample_ID")


# ── Collapse Duplicate Phylum–Sample Rows ─────────────────────────────────────

asv_summarised <- asv_with_metadata %>%
  group_by(across(-Count)) %>%
  summarise(Total_count = sum(Count), .groups = "drop")


# ── Identify Top 5 Phyla ──────────────────────────────────────────────────────

# Rank phyla by total read count across all samples
top_5_phyla <- asv_summarised %>%
  group_by(Phylum) %>%
  summarise(Total_abundance = sum(Total_count), .groups = "drop") %>%
  arrange(desc(Total_abundance)) %>%
  slice_head(n = 5) %>%
  pull(Phylum)


# ── Group Remaining Phyla as "Other" ─────────────────────────────────────────

asv_summarised <- asv_summarised %>%
  mutate(Phylum_grouped = ifelse(Phylum %in% top_5_phyla, Phylum, "Other"))

# ── Aggregate by Sample and Grouped Phylum, Then Compute Relative Abundance ───

# Sum raw counts AFTER grouping into "Other", then compute RA from scratch.
# This avoids floating-point drift that causes per-sample totals to be
# slightly less than 1?

plot_phyla <- asv_summarised %>%
  group_by(Sample_ID, Phylum_grouped) %>%
  summarise(Total_count = sum(Total_count), .groups = "drop") %>%
  group_by(Sample_ID) %>%
  mutate(
    sample_total = sum(Total_count, na.rm = TRUE),
    Rel_abundance = ifelse(sample_total == 0, 0, Total_count / sample_total)
  ) %>%
  select(-sample_total) %>%
  ungroup()


# ── Handle Zero-Count Samples (Empty Negative Controls) ──────────────────────

# Separate zero-count and non-zero samples
non_empty <- plot_phyla %>%
  group_by(Sample_ID) %>%
  filter(sum(Rel_abundance) > 0) %>%
  ungroup()

# For zero-count samples, replace all rows with a single placeholder row
empty <- plot_phyla %>%
  group_by(Sample_ID) %>%
  filter(sum(Rel_abundance) == 0) %>%
  ungroup() %>%
  distinct(Sample_ID) %>%
  mutate(
    Phylum_grouped = "Empty",
    Total_count = 0,
    Rel_abundance = 1
  )

# Recombine
plot_phyla <- bind_rows(non_empty, empty)

# Re-attach metadata for faceting and ordering
plot_phyla <- merge(plot_phyla, metadata, by = "Sample_ID")


# ── Colour Palette ────────────────────────────────────────────────────────────

my_palette <- c(
  "Pseudomonadota" = "#F56A6A",
  "Bacillota" = "#1AB788",
  "Bacteroidota" = "#ffdc5d",
  "Actinomycetota" = "#074151",
  "Chlamydiota" = "#00ffff",
  "Other" = "#A3334E",
  "Empty" = "#D3D3D3"
)


# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 1 — Stacked Bar Chart (Ordered by Pseudomonadota Relative Abundance)
# ═══════════════════════════════════════════════════════════════════════════════

# Step 1: Determine sample order within each Stage based on Pseudomonadota RA,
# and append empty samples with Pseudo_RA = 0 so they get a valid
# factor level and are not rendered as NA by ggplot
pseudomonadota_order <- plot_phyla %>%
  filter(Phylum_grouped == "Pseudomonadota") %>%
  group_by(Stage, Sample_ID) %>%
  summarise(Pseudo_RA = sum(Rel_abundance), .groups = "drop") %>%
  
  # Append empty samples that have no Pseudomonadota rows
  bind_rows(
    plot_phyla %>%
      filter(Phylum_grouped == "Empty") %>%
      distinct(Sample_ID, Stage) %>%
      mutate(Pseudo_RA = 0)
  ) %>%
  
  arrange(Stage, Pseudo_RA) %>%
  group_by(Stage) %>%
  mutate(Sample_ID_ordered = factor(Sample_ID, levels = Sample_ID)) %>%
  ungroup()

# Step 2: Join ordering back to the main plotting data
pseudomona_abundance <- plot_phyla %>%
  left_join(
    pseudomonadota_order %>% select(Sample_ID, Stage, Sample_ID_ordered),
    by = c("Sample_ID", "Stage")
  )

# Step 3: Plot
pseudomona_abundance_viz <- pseudomona_abundance %>%
  ggplot(aes(x = Sample_ID_ordered, y = Rel_abundance, fill = Phylum_grouped)) +
  
  # Seamless bars: no border colour, no border width
  geom_bar(
    stat = "identity",
    position = "fill",
    width = 1,
    color = NA,
    linewidth = 0
  ) +
  
  facet_grid(~ Stage, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = my_palette) +
  
  # Remove internal padding so bars fill panels edge-to-edge
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  
  theme_minimal(base_size = 18) +
  
  labs(
    x = "Sample ID",
    y = "Relative Abundance",
    fill = "Phylum"
  ) +
  
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    # Remove grid lines that can bleed through bars
    panel.grid = element_blank(),
    panel.background = element_blank()
    
    # panel.spacing intentionally left at default to preserve
    # the white space separating Stage facet groups
  )

pseudomona_abundance_viz


# ── Export ────────────────────────────────────────────────────────────────────
#saved as png for the in-text figure
ggsave(
  filename = "../figures/top_5_bacterial_phyla_removed_contamination.svg",
  plot = pseudomona_abundance_viz,
  width = 14,
  height = 6,
  units = "in"
)
