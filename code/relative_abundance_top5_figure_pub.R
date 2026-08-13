# =============================================================================
# Date:    09/07/2025
# Author:  Agnes Sz.
# Purpose: Visualise the relative abundance of the top 5 phyla across samples,
#          grouping all remaining phyla under "Other".
#          - All bacteria included; chloroplasts and mitochondria excluded.
#          - Dataset has been pre-filtered using decontam (prevalence threshold: 0.57).
# =============================================================================

# ── Libraries ─────────────────────────────────────────────────────────────────

library(tidyverse)

# ── Load Data ─────────────────────────────────────────────────────────────────

asv_raw        <- read.csv("../data/ASV_with_taxonomy_p57.csv")
metadata_raw   <- read.csv("../data/Microbial_metadata_cleaned.csv")


# ── Clean Metadata ────────────────────────────────────────────────────────────

# Replace NA values with sensible defaults for downstream grouping
metadata <- metadata_raw %>%
  mutate(
    Stage        = replace_na(Stage, "Negative"),
    Capture_Date = replace_na(Capture_Date, 70)
  )


# Retain only the Phylum column and sample count columns.
# NOTE: Columns 331–337 are non-sample metadata columns; adjust indices if the
#       input file structure changes.
asv_filtered <- asv_raw %>%
  select(-c(331, 333, 334, 335, 336, 337))


# ── Reshape to Long Format ────────────────────────────────────────────────────

# Pivot so each row represents one phylum–sample observation
asv_long <- asv_filtered %>%
  pivot_longer(
    cols      = -Phylum,
    names_to  = "Sample_ID",
    values_to = "Count"
  )


# ── Merge with Metadata ───────────────────────────────────────────────────────

asv_with_metadata <- merge(asv_long, metadata, by = "Sample_ID")


# ── Compute Relative Abundance ────────────────────────────────────────────────

# Step 1: Collapse duplicate phylum–sample rows by summing counts
asv_summarised <- asv_with_metadata %>%
  group_by(across(-Count)) %>%
  summarise(Total_count = sum(Count), .groups = "drop")

# Step 2: Calculate relative abundance per sample.
#         Samples with a zero total are assigned 0 to avoid division by zero.
asv_relabu <- asv_summarised %>%
  group_by(Sample_ID) %>%
  mutate(
    sample_total = sum(Total_count, na.rm = TRUE),
    Rel_abundance = ifelse(sample_total == 0, 0, Total_count / sample_total)
  ) %>%
  select(-sample_total) %>%
  ungroup()


# ── Identify Top 5 Phyla ──────────────────────────────────────────────────────

# Rank phyla by total read count across all samples
top_5_phyla <- asv_relabu %>%
  group_by(Phylum) %>%
  summarise(Total_abundance = sum(Total_count), .groups = "drop") %>%
  arrange(desc(Total_abundance)) %>%
  slice_head(n = 5) %>%
  pull(Phylum)


# ── Group Remaining Phyla as "Other" ─────────────────────────────────────────

asv_relabu <- asv_relabu %>%
  mutate(Phylum_grouped = ifelse(Phylum %in% top_5_phyla, Phylum, "Other"))


# ── Aggregate Relative Abundance by Sample and Grouped Phylum ─────────────────

plot_phyla <- asv_relabu %>%
  group_by(Sample_ID, Phylum_grouped) %>%
  summarise(Rel_abundance = sum(Rel_abundance), .groups = "drop")

# Re-attach metadata for faceting and ordering in plots
plot_phyla <- merge(plot_phyla, metadata, by = "Sample_ID")


# ── Colour Palettes ───────────────────────────────────────────────────────────

# Named palette ensures consistent phylum–colour mapping across all plots
my_palette <- c(
  "Pseudomonadota"  = "#F56A6A",
  "Bacillota"       = "#1AB788",
  "Bacteroidota"    = "#ffdc5d",
  "Actinomycetota"  = "#074151",
  "Chlamydiota"     = "#00ffff",
  "Other"           = "#A3334E"
)


# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 1 — Stacked Bar Chart (Ordered by Pseudomonadota Relative Abundance)
# ═══════════════════════════════════════════════════════════════════════════════

# Step 1: Determine sample order within each Stage based on Pseudomonadota RA
pseudomonadota_order <- plot_phyla %>%
  filter(Phylum_grouped == "Pseudomonadota") %>%
  group_by(Stage, Sample_ID) %>%
  summarise(Pseudo_RA = sum(Rel_abundance), .groups = "drop") %>%
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

pseudomona_abundance_viz <- pseudomona_abundance %>%
  ggplot(aes(x = Sample_ID_ordered, y = Rel_abundance, fill = Phylum_grouped)) +
  
  # KEY FIX 1: Remove bar border seams within each group
  geom_bar(
    stat      = "identity",
    position  = "fill",
    width     = 1,
    color     = NA,        # no border colour
    linewidth = 0          # zero border width
  ) +
  
  facet_grid(~ Stage, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = my_palette) +
  
  # KEY FIX 2: Remove internal padding so bars fill the panel edge-to-edge
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  
  theme_minimal(base_size = 18) +
  
  labs(
    x    = "Sample ID",
    y    = "Relative Abundance",
    fill = "Phylum"
  ) +
  
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    
    # KEY FIX 3: Remove internal grid lines that bleed through bars
    panel.grid   = element_blank(),
    
    # ✅ panel.spacing is intentionally left at its default
    #    so the white space between Stage groups is preserved
  )

pseudomona_abundance_viz

# Save the combined figure as an SVG for publication-quality vector output
ggsave(
  filename = "../figures/top_5_bacterial_phyla_removed_contamination.svg",
  plot     = pseudomona_abundance_viz,
  device   = "svg",
  width    = 15,
  height   = 7,
  units    = "in"
)


# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 2 — Aggregated Bar Chart by Stage (Negative Controls Excluded)
# ═══════════════════════════════════════════════════════════════════════════════

# Exclude negative controls and aggregate relative abundance by Stage,
# so that each stage is represented by a single mean value per phylum
# rather than one bar segment per sample (which causes the striped effect)

plot_phyla_barplot <- plot_phyla %>%
  filter(Stage != "Negative") %>%
  group_by(Stage, Phylum_grouped) %>%
  summarise(Rel_abundance = mean(Rel_abundance), .groups = "drop") %>%
  mutate(
    Phylum_grouped = factor(
      Phylum_grouped,
      levels = c("Chlamydiota", "Other", "Bacteroidota",
                 "Actinomycetota", "Bacillota", "Pseudomonadota")
    )
  )

# Plot
bp_phyla_silva <- ggplot(
  data = plot_phyla_barplot,
  aes(x = Stage, y = Rel_abundance, fill = Phylum_grouped)
) +
  geom_bar(stat = "identity", position = "fill", alpha = 1) +
  scale_fill_manual(values = my_palette) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal(base_size = 15) + 
  labs(
    x    = "Stage",
    y    = "Relative abundance",
    fill = "Phylum",
  ) +
  theme(legend.position = "none")


bp_phyla_silva

# Save the combined figure as an SVG for publication-quality vector output
ggsave(
  filename = "../figures/bp_phyla_silvaaa.svg",
  plot     = bp_phyla_silva,
  device   = "svg",
  width    = 5,
  height   = 4,
  units    = "in"
)

# ── Summary Statistics by Stage and Phylum ────────────────────────────────────

plot_phyla_barplot_stats <- plot_phyla_barplot %>%
  group_by(Stage, Phylum_grouped) %>%
  summarise(Rel_abundance = sum(Rel_abundance), .groups = "drop")
