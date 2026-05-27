# =============================================================================
# Script:  identify_contaminants_decontam_SILVA.R
# Purpose: Identify and remove contaminant ASVs from the ASV table using the
#          decontam prevalence method, and export clean ASV tables for
#          downstream analysis.
# Author:  Agnes Sz.
# Date:    2025-07-09
# =============================================================================


# ── 1. Libraries ──────────────────────────────────────────────────────────────

library(tidyverse)   # Includes ggplot2, dplyr, readr
library(decontam)    # Contaminant identification
library(phyloseq)    # Phylogenetic sequencing data structures
library(Biostrings)  # DNA sequence handling (used for ASV sequences)


# ── 2. Load Input Data ────────────────────────────────────────────────────────

# ASV abundance table (ASVs as rows, samples as columns)
asv_table <- read.csv("../data/ASV_table_cleaned.csv", row.names = 1)

# Sample metadata
microbial_metadata <- read.csv("../data/Microbial_metadata_cleaned.csv", row.names = 1)

# Reference taxonomy (converted to matrix as required by phyloseq)
tax_total <- read.csv("../data/tax.total.csv", row.names = 1)
tax_total <- as.matrix(tax_total)


# ── 3. Prepare Metadata ───────────────────────────────────────────────────────
# Negative control samples may have NA in Stage and Capture_Date.
# These are filled with placeholder values so phyloseq handles them correctly.

microbial_metadata$Stage[is.na(microbial_metadata$Stage)]               <- "Negative"
microbial_metadata$Capture_Date[is.na(microbial_metadata$Capture_Date)] <- 70  # Dummy ordinal day


# ── 4. Build phyloseq Object (with Contaminants) ──────────────────────────────

phylo_with_cont <- phyloseq(
  otu_table(asv_table, taxa_are_rows = TRUE),
  tax_table(tax_total),
  sample_data(microbial_metadata)
)

# Attach ASV sequences as a DNAStringSet
dna <- Biostrings::DNAStringSet(taxa_names(phylo_with_cont))
names(dna)      <- taxa_names(phylo_with_cont)
phylo_with_cont <- merge_phyloseq(phylo_with_cont, dna)

# ── ⚠️ Save the sequence → ASV ID mapping BEFORE renaming ────────────────────
# taxa_names() currently holds raw nucleotide sequences. Once we rename to
# ASV1, ASV2... the original sequences are lost from the names — so we
# capture the mapping here for use in contaminant removal later.

asv_id_map <- tibble(
  Sequence = taxa_names(phylo_with_cont),                          # Raw nucleotide sequence
  ASV_ID   = paste0("ASV", seq(ntaxa(phylo_with_cont)))            # Corresponding ASV ID
)

# Now safe to rename
taxa_names(phylo_with_cont) <- asv_id_map$ASV_ID


# ── 5. Identify Contaminants via Prevalence Method ────────────────────────────

sample_data(phylo_with_cont)$is.neg <- sample_data(phylo_with_cont)$Stage == "Negative"

contamdf <- isContaminant(
  phylo_with_cont,
  method    = "prevalence",
  neg       = "is.neg",
  threshold = 0.57
)

table(contamdf$contaminant)

contaminant_indices <- which(contamdf$contaminant)

# Retrieve the ASV IDs of contaminants, then map back to raw sequences
contaminant_asv_ids   <- taxa_names(phylo_with_cont)[contaminant_indices]
contaminant_sequences <- asv_id_map %>%
  filter(ASV_ID %in% contaminant_asv_ids) %>%
  pull(Sequence)


# ── 7. Build Decontaminated ASV Table ─────────────────────────────────────────

# Add sequence as an explicit column for merging
taxa_df          <- as.data.frame(tax_total)
taxa_df$Sequence <- rownames(taxa_df)       # Row names are raw nucleotide sequences

asv_table$Sequence <- rownames(asv_table)   # Row names are raw nucleotide sequences

# Merge ASV counts with taxonomy on the shared Sequence column
asv_with_taxonomy <- merge(asv_table, taxa_df, by = "Sequence")

# Sanity check before removal
cat("Contaminants identified:  ", length(contaminant_sequences), "\n")
cat("ASVs before removal:      ", nrow(asv_with_taxonomy), "\n")

# Remove contaminants by matching raw nucleotide sequences
asv_with_taxonomy <- asv_with_taxonomy %>%
  filter(!Sequence %in% contaminant_sequences)

cat("ASVs after removal:       ", nrow(asv_with_taxonomy), "\n")

# Confirm the difference between the initial number of ASVs and the number of ASVs after assigning the taxonomy (bigest loss) and removing contaminants
cat("ASVs removed:             ", 
    length(contaminant_sequences) - 
      (nrow(asv_with_taxonomy) - nrow(asv_table) + length(contaminant_sequences)), "\n")

# Drop the Sequence column before export
asv_with_taxonomy <- asv_with_taxonomy %>%
  select(-Sequence)
# ── 8. Export Decontaminated ASV Tables ───────────────────────────────────────
# Two versions are exported:
#   - row_names = TRUE  → for phyloseq (requires ASV IDs as row names)
#   - row_names = FALSE → for relative abundance visualisation and general use

write.csv(
  asv_with_taxonomy,
  "../data/ASV_with_taxonomy_p57_row_names.csv",
  row.names = TRUE,
  quote     = FALSE
)

write.csv(
  asv_with_taxonomy,
  "../data/ASV_with_taxonomy_p57.csv",
  row.names = FALSE,
  quote     = FALSE
)

# ── 8b Checking percentage of ASVs that were assigned at the level of... ───────────────────────────────────────

kk <- read.csv("../data/ASV_with_taxonomy_p57.csv") # 2 742

kk2 <- kk %>% filter(!is.na(Phylum)) # 2 739, 99,9%

kk3 <- kk %>% filter(!is.na(Class)) # 2 725, 99,4%

kk4 <- kk %>% filter(!is.na(Order)) # 2 682, 97,8%

kk5 <- kk %>% filter(!is.na(Family)) # 2 568, 93,7%

kk6 <- kk %>% filter(!is.na(Genus)) # 2 329, 84,9%

kk7 <- kk %>% filter(!is.na(Species)) # 1 242, 45,3%

# ── 9. Prevalence Plot: Contaminants vs True Samples ─────────────────────────
# Visualises ASV prevalence in negative controls vs. true samples.
# Contaminants (flagged by decontam) are expected to cluster toward the
# negative control axis. A clean dataset should show most ASVs along the
# true sample axis.

ps_pa     <- transform_sample_counts(phylo_with_cont, function(abund) 1 * (abund > 0))
ps_pa_neg <- prune_samples(sample_data(ps_pa)$Stage == "Negative", ps_pa)
ps_pa_pos <- prune_samples(
  sample_data(ps_pa)$Stage == "Female" | sample_data(ps_pa)$Stage == "Chick",
  ps_pa
)

# ── Plot ────────────────────────────────────────

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
  "../figures/contaminant_prevalence_plot.svg",
  plot   = plot_contaminants,
  width  = 210,
  height = 99,
  units  = "mm"
)
