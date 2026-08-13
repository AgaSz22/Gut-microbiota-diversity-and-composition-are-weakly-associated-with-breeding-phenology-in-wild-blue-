# =============================================================================
# Author:  Agnes Szwarczynska
# Date:    29/09/2025
# Purpose: Visualise the relative abundance of the top 5 families across
#          sample groups, grouping all remaining families under "Other".
#          Data sourced from Bracken family-level metagenomic classifications.
# =============================================================================


# ── Libraries ─────────────────────────────────────────────────────────────────

library(dplyr)      # Data manipulation
library(tidyr)      # Reshaping (pivot_wider, pivot_longer)
library(ggplot2)    # Visualisation
library(reshape2)   # melt() for wide-to-long conversion
library(forcats)    # Factor reordering
library(phyloseq)   # Microbiome data structures
#library(DT)         # Interactive data tables


# ── Define File Paths and Group Labels ────────────────────────────────────────

# Base directory containing Bracken family-level output files
bracken_dir <- "../data"

# Named vector: group label → filename
bracken_files <- c(
  "Early_Chicks"  = "1-Early_Chicks_TCAATTCCGA-TGTAGATCCG.bracken",
  "Late_Chicks"   = "2-Late_Chicks_CGGAACACAA-AACTCTACCG.bracken",
  "Early_Females" = "3-Early_Females_CTATAACGGA-CGTGACATTC.bracken",
  "Late_Females"  = "4-Late_Females_AGATGTAAGG-CTGTACCAGG.bracken"
)

# ── Load and Merge Bracken Data ───────────────────────────────────────────────

# Read each Bracken file, attach its group label, and combine into one table.
# Column 1 = taxonomy name; column 6 = estimated read count (new_est_reads).
all_groups_bracken <- lapply(names(bracken_files), function(group) {
  
  filepath <- file.path(bracken_dir, bracken_files[group])
  df       <- read.table(filepath, header = TRUE, sep = "\t")
  df$Group <- group
  return(df)
  
}) %>%
  bind_rows() %>%
  select(1, 6, Group) %>%
  setNames(c("Phylum", "Count", "Group")) %>%  # ← Fixed here
  filter(Phylum != "Chordata")

head(all_groups_bracken)

# ── Compute Relative Abundance ────────────────────────────────────────────────

# Step 2: Calculate relative abundance per group.
#         Groups with zero total reads are assigned 0 to avoid division by zero.
bracken_relabu <- all_groups_bracken %>%
  group_by(Phylum, Group) %>%
  summarise(Total_Count = sum(Count), .groups = "drop") %>%
  mutate(
    group_total   = sum(Total_Count, na.rm = TRUE),
    Rel_Abundance = ifelse(group_total == 0, 0, Total_Count / group_total)
  ) %>%
  select(-group_total) %>%
  ungroup()


# ── Identify Top 5 Families ───────────────────────────────────────────────────

# Rank families by total read count across all groups
top_5_phyla <- bracken_relabu %>%
  group_by(Phylum) %>%
  summarise(Total_Abundance = sum(Total_Count), .groups = "drop") %>%
  arrange(desc(Total_Abundance)) %>%
  slice_head(n = 5) %>%
  pull(Phylum)

print(top_5_phyla)


# ── Group Remaining Families as "Other" ───────────────────────────────────────

bracken_relabu <- bracken_relabu %>%
  mutate(
    Phylum_grouped = ifelse(Phylum %in% top_5_phyla, Phylum, "Other")
  )

# ── Aggregate Relative Abundance by Group and Grouped Family ─────────────────

plot_phyla <- bracken_relabu %>%
  group_by(Group, Phylum_grouped) %>%
  summarise(Rel_Abundance = sum(Rel_Abundance), .groups = "drop")

# Add a Stage column derived from Group for faceting (mirrors phylum script)
# Early/Late stage is extracted from the group name

plot_phyla <- plot_phyla %>%
  mutate(
    Stage = case_when(
      grepl("Chicks", Group) ~ "Chicks",
      grepl("Females",  Group) ~ "Females",
      TRUE                  ~ "Unknown"
    )
  )

head(plot_phyla)

# ── Colour Palette ────────────────────────────────────────────────────────────

# Placeholder palette — replace hex codes with your preferred colours once
# top_5_families has been determined from the data.
# "Other" is always included as the catch-all category.

my_palette <- c(
  "Pseudomonadota"  = "#F56A6A",
  "Bacillota"       = "#1AB788",
  "Bacteroidota"    = "#ffdc5d",
  "Actinomycetota"  = "#074151",
  "Mycoplasmatota"  = "#a1b054",
  "Other"           = "#A3334E"
)

# ═════════════════════════════════════════════════════════════════════════════
# PLOT 1 — Aggregated Bar Chart by Stage
# ═════════════════════════════════════════════════════════════════════════════

# Aggregate relative abundance across all groups within each Stage

plot_phyla_by_stage <- plot_phyla %>%
  group_by(Stage, Phylum_grouped) %>%
  summarise(Rel_Abundance = sum(Rel_Abundance), .groups = "drop") %>%
  mutate(
    Phylum_grouped = factor(
      Phylum_grouped,
      levels = c("Other","Bacteroidota","Actinomycetota", "Mycoplasmatota", "Bacillota", "Pseudomonadota")   # "Other" always at the bottom
    )
  )

#Plot

plot_phyla_stage <- ggplot(
  data = plot_phyla_by_stage,
  aes(x = Stage, y = Rel_Abundance, fill = Phylum_grouped)
) +
  geom_bar(
    stat      = "identity",
    position  = "fill",
    width     = 0.9,
    color     = NA,
    linewidth = 0
  ) +
  scale_fill_manual(values = my_palette) +
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = c(0, 0),
    limits = c(0, 1)
  ) +
  scale_x_discrete(expand = c(0.3, 0.3)) +
  theme_minimal(base_size = 15) +
  labs(
    x    = NULL,
    y    = "Relative abundance",
    fill = "Phylum"
  ) +
  theme(
    axis.text.x        = element_text(size = 13),
    axis.text.y        = element_text(size = 11),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "none"
  )

plot_phyla_stage

# Save as SVG
ggsave(
  filename = "../figures/top_5_families_by_stage_metagenomics.svg",
  plot     = plot_phyla_stage,
  device   = "svg",
  width    = 5,
  height   = 4,
  units    = "in"
)
