library(tidyverse)
library(scales)

# ─────────────────────────────────────────────────────────────────────────────
# LOAD BOTH DATASETS
# ─────────────────────────────────────────────────────────────────────────────

blast_raw    <- read.csv("../data/ASV_with_taxonomy_p57_Blast_with_plastids.csv")
silva_raw    <- read.csv("../data/ASV_with_taxonomy_p57_SILVA_with_plastids.csv")
metadata_raw <- read.csv("../data/Microbial_metadata_cleaned.csv")

colnames(blast_raw)[colnames(blast_raw) == "phylum"] <- "Phylum"

blast_filtered <- blast_raw %>%
  select(-c(331, 333, 334, 335, 336, 337))

silva_filtered <- silva_raw %>%
  select(-c(331, 333, 334, 335, 336, 337))

# ── Clean Metadata ────────────────────────────────────────────────────────────

metadata <- metadata_raw %>%
  mutate(
    Stage        = replace_na(Stage, "Negative"),
    Capture_Date = replace_na(Capture_Date, 70)
  )

# ── Reshape to Long Format ────────────────────────────────────────────────────

blast_long <- blast_filtered %>%
  pivot_longer(
    cols      = -Phylum,
    names_to  = "Sample_ID",
    values_to = "Count"
  )

silva_long <- silva_filtered %>%
  pivot_longer(
    cols      = -Phylum,
    names_to  = "Sample_ID",
    values_to = "Count"
  )

# ── Join Stage info into both datasets ───────────────────────────────────────

blast_long <- blast_long %>%
  left_join(metadata %>% select(Sample_ID, Stage), by = "Sample_ID")

silva_long <- silva_long %>%
  left_join(metadata %>% select(Sample_ID, Stage), by = "Sample_ID")

# ─────────────────────────────────────────────────────────────────────────────
# PREPROCESSING FUNCTION
# ─────────────────────────────────────────────────────────────────────────────

preprocess_phyla <- function(long_df, dataset_label) {
  long_df %>%
    # Remove X prefix from Sample_ID if added by read.csv()
    mutate(Sample_ID = str_remove(Sample_ID, "^X")) %>%
    # Calculate relative abundance per sample
    group_by(Sample_ID) %>%
    mutate(Rel_abundance = Count / sum(Count, na.rm = TRUE)) %>%
    ungroup() %>%
    # Group phyla
    mutate(
      Dataset = dataset_label,
      Phylum_grouped = case_when(
        Phylum == "Pseudomonadota"  ~ "Pseudomonadota",
        Phylum == "Bacillota"       ~ "Bacillota",
        Phylum == "Actinomycetota"  ~ "Actinomycetota",
        Phylum == "Bacteroidota"    ~ "Bacteroidota",
        Phylum == "Chlamydiota"     ~ "Chlamydiota",
        TRUE                        ~ "Other"
      )
    ) %>%
    # Aggregate grouped phyla per sample
    group_by(Sample_ID, Stage, Dataset, Phylum_grouped) %>%
    summarise(Rel_abundance = sum(Rel_abundance, na.rm = TRUE), .groups = "drop")
}

plot_phyla_blast <- preprocess_phyla(blast_long, "BLAST")
plot_phyla_silva <- preprocess_phyla(silva_long, "SILVA")

# ─────────────────────────────────────────────────────────────────────────────
# COLOUR PALETTE
# ─────────────────────────────────────────────────────────────────────────────

my_palette <- c(
  "Pseudomonadota"  = "#F56A6A",
  "Bacillota"       = "#1AB788",
  "Bacteroidota"    = "#ffdc5d",
  "Actinomycetota"  = "#074151",
  "Mycoplasmatota"  = "#a1b054",
  "Other"           = "#A3334E",
  "Chlamydiota"     = "#00ffff"    
)


# ─────────────────────────────────────────────────────────────────────────────
# BAR CHART FUNCTION
# ─────────────────────────────────────────────────────────────────────────────

make_barplot <- function(plot_phyla, dataset_label) {
  
  # Aggregate: mean relative abundance per Stage x Phylum (exclude negatives)
  plot_data <- plot_phyla %>%
    filter(Stage != "Negative") %>%          # exclude negatives
    group_by(Stage, Phylum_grouped) %>%      # keep Stage separation
    summarise(Rel_abundance = mean(Rel_abundance, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      Phylum_grouped = factor(
        Phylum_grouped,
        levels = c(
           "Chlamydiota", "Other","Bacteroidota",
          "Actinomycetota", "Bacillota", "Pseudomonadota"
        )
      )
    )
  
  # Build plot
  ggplot(plot_data, aes(x = Stage, y = Rel_abundance, fill = Phylum_grouped)) +
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
      x     = NULL,
      y     = "Relative abundance",
      fill  = "Phylum"
    ) +
    theme(
      axis.text.x        = element_text(size = 13),
      axis.text.y        = element_text(size = 11),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.position    = "none"
    )
}
# ── Generate both plots ───────────────────────────────────────────────────────
bp_blast <- make_barplot(plot_phyla_blast, "BLAST")
bp_silva  <- make_barplot(plot_phyla_silva, "SILVA")

bp_blast
bp_silva
# ─────────────────────────────────────────────────────────────────────────────
# SAVE
# ─────────────────────────────────────────────────────────────────────────────

ggsave(
  filename = "../figures/bp_phyla_blast.svg",
  plot     = bp_blast,
  device   = "svg",
  width    = 5,
  height   = 4,
  units    = "in"
)

ggsave(
  filename = "../figures/bp_phyla_silva.svg",
  plot     = bp_silva,
  device   = "svg",
  width    = 5,
  height   = 4,
  units    = "in"
)
