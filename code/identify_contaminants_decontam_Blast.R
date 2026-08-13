# =============================================================================
# Script:  identify_contaminants_decontam_Blast.R
# Purpose: Identify and remove contaminant ASVs from the BLAST-assigned ASV
#          table using the decontam prevalence method, and export clean ASV
#          tables for downstream analysis.
# Author:  Agnes Sz.
# Date:    2025-07-28
# =============================================================================


# ── 1. Libraries ──────────────────────────────────────────────────────────────

library(tidyverse)   # Includes ggplot2, dplyr, readr
library(decontam)    # Contaminant identification
library(phyloseq)    # Phylogenetic sequencing data structures
library(Biostrings)  # DNA sequence handling


# ── 2. Load Input Data ────────────────────────────────────────────────────────

# ASV abundance table — BLAST taxonomy version (ASVs as rows, samples as columns)
asv_table <- read.csv("../data/ASV_table_Blast.csv", row.names = 1)

# Sample metadata
microbial_metadata <- read.csv("../data/Microbial_metadata_cleaned.csv", row.names = 1)

# BLAST reference taxonomy (converted to matrix as required by phyloseq)
tax_total <- read.csv("../data/tax.total_Blast.csv", row.names = 1)
tax_total <- as.matrix(tax_total)


# ── 3. Prepare Metadata ───────────────────────────────────────────────────────
# Negative control samples may have NA in Stage and Capture_Date.
# These are filled with placeholder values so phyloseq handles them correctly.

microbial_metadata$Stage[is.na(microbial_metadata$Stage)]               <- "Negative"
microbial_metadata$Capture_Date[is.na(microbial_metadata$Capture_Date)] <- 70  # Dummy ordinal day


# ── 4. Build phyloseq Object (with Contaminants) ──────────────────────────────
# NOTE: Unlike the SILVA pipeline, the BLAST taxonomy table does not use raw
#       nucleotide sequences as row names — ASV IDs are already meaningful.
#       The DNAStringSet step is therefore not needed here, and taxa_names()
#       already holds stable ASV identifiers that can be used directly for
#       contaminant removal without any renaming step.

phylo_with_cont <- phyloseq(
  otu_table(asv_table, taxa_are_rows = TRUE),
  tax_table(tax_total),
  sample_data(microbial_metadata)
)


# ── 5. Identify Contaminants via Prevalence Method ────────────────────────────
# Threshold of 0.57 is a relatively strict setting (default is 0.1 for the
# prevalence method). ASVs more prevalent in negative controls than in true
# samples are flagged as contaminants.
# See: Davis et al. (2018) Microbiome — doi:10.1186/s40168-018-0605-2

sample_data(phylo_with_cont)$is.neg <- sample_data(phylo_with_cont)$Stage == "Negative"

contamdf <- isContaminant(
  phylo_with_cont,
  method    = "prevalence",
  neg       = "is.neg",
  threshold = 0.57
)

# Summary of contaminant calls
table(contamdf$contaminant)

# Retrieve contaminant ASV names directly — no renaming was applied so these
# are stable identifiers that match row names in the ASV table
contaminant_indices <- which(contamdf$contaminant)
contaminant_names   <- taxa_names(phylo_with_cont)[contaminant_indices]


# ── 6. Export Contaminant List ────────────────────────────────────────────────

# Extract taxonomy of flagged contaminants
list_of_contaminants <- as.data.frame(
  as(tax_table(phylo_with_cont), "matrix")[contaminant_indices, ]
)

# Add ASV ID as an explicit column
list_of_contaminants$ASV_ID <- rownames(list_of_contaminants)

write.csv(
  list_of_contaminants,
  "../results/list_of_contaminants_Blast.csv",
  row.names = FALSE  # ASV_ID is now an explicit column, so row names are redundant
)

# Save contaminant indices for reproducibility
saveRDS(contaminant_indices, file = "../results/contaminants_indices_Blast.RData")


# ── 7. Build Decontaminated ASV Table ─────────────────────────────────────────
# Merge ASV table with taxonomy by sequence (row names), then remove
# contaminants by ASV name to avoid any index misalignment after merge().

# Add sequence as an explicit column for merging
taxa_df          <- as.data.frame(tax_total)
taxa_df$Sequence <- rownames(taxa_df)

asv_table$Sequence <- rownames(asv_table)

# Merge ASV counts with taxonomy on the shared Sequence column
asv_with_taxonomy <- merge(asv_table, taxa_df, by = "Sequence")

# Sanity check before removal
cat("Contaminants identified:  ", length(contaminant_names), "\n")
cat("ASVs before removal:      ", nrow(asv_with_taxonomy), "\n")

# Remove contaminants by matching ASV names against the Sequence column
# This is safe regardless of row order after merge()
asv_with_taxonomy <- asv_with_taxonomy %>%
  filter(!Sequence %in% contaminant_names)

cat("ASVs after removal:       ", nrow(asv_with_taxonomy), "\n")

# Confirm the difference between the initial number of ASVs and the number of ASVs after assigning the taxonomy (bigest loss) and removing contaminants
cat("ASVs removed:             ", 
    nrow(asv_table) - (nrow(asv_with_taxonomy) - length(contaminant_indices)), "\n")

# Drop the Sequence column before export
asv_with_taxonomy <- asv_with_taxonomy %>%
  select(-Sequence)

write.csv(
  asv_with_taxonomy,
  "../data/ASV_with_taxonomy_p57_Blast_with_plastids.csv",
  row.names = FALSE,
  quote     = FALSE
)


##### Removing plasmids

asv_with_taxonomy_all_with_phyla <- asv_with_taxonomy %>%
  filter(!is.na(phylum))
diff_phyla <- nrow(asv_with_taxonomy) - nrow(asv_with_taxonomy_all_with_phyla)
cat("The number of ASVs not identifiabe at the phylum level (removed):  ", diff_phyla, "\n")

plastids <- asv_with_taxonomy_all_with_phyla %>% subset(family == "Mitochondria" | order == "Chloroplast")

asv_with_taxonomy_no_plastids <- asv_with_taxonomy_all_with_phyla %>%
  filter(family != "Mitochondria" & order != "Chloroplast")
cat("The number of plastid ASVs removed:  ", nrow(plastids), "\n")


# ── 8. Export Decontaminated ASV Tables ───────────────────────────────────────
# Two versions exported:
#   - row.names = TRUE  → for phyloseq (requires ASV IDs as row names)
#   - row.names = FALSE → for relative abundance visualisation and general use

write.csv(
  asv_with_taxonomy_no_plastids,
  "../data/ASV_with_taxonomy_p57_row_names_Blast.csv",
  row.names = TRUE,
  quote     = FALSE
)

write.csv(
  asv_with_taxonomy_no_plastids,
  "../data/ASV_with_taxonomy_p57_Blast.csv",
  row.names = FALSE,
  quote     = FALSE
)


# ── 9. Prevalence Plot: Contaminants vs True Samples ─────────────────────────
# Visualises ASV prevalence in negative controls vs. true samples.
# Contaminants are expected to cluster toward the negative control axis.

ps_pa     <- transform_sample_counts(phylo_with_cont, function(abund) 1 * (abund > 0))
ps_pa_neg <- prune_samples(sample_data(ps_pa)$Stage == "Negative", ps_pa)
ps_pa_pos <- prune_samples(
  sample_data(ps_pa)$Stage == "Female" | sample_data(ps_pa)$Stage == "Chick",
  ps_pa
)

# ── Plot ───────────────────────────────────────────────────────────

prevalence_df <- data.frame(
  pa.pos      = taxa_sums(ps_pa_pos),
  pa.neg      = taxa_sums(ps_pa_neg),
  contaminant = contamdf$contaminant
)

plot_contaminants <- ggplot(prevalence_df, aes(x = pa.neg, y = pa.pos, color = contaminant)) +
  geom_point() +
  xlab("Prevalence (Negative Controls)") +
  ylab("Prevalence (True Samples)") +
  scale_color_manual(
    values = c("FALSE" = "grey60", "TRUE" = "#A3334E"),
    labels = c("FALSE" = "Non-contaminant", "TRUE" = "Contaminant")
  ) +
  theme_bw(base_size = 14)

print(plot_contaminants)

ggsave(
  "../results/contaminant_prevalence_plot_Blast.svg",
  plot   = plot_contaminants,
  width  = 210,
  height = 99,
  units  = "mm"
)
