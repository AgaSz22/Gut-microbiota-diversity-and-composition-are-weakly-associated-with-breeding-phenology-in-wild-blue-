# =============================================================================
# Script:  build_reference_taxonomy.R
# Purpose: Combine two SILVA-assigned taxonomy tables, remove duplicates,
#          and export a single clean reference taxonomy file.
# Author:  Agnes Szwarczynska
# Date:    2026-05-27
# =============================================================================


# ── 1. Libraries ──────────────────────────────────────────────────────────────

library(tidyverse)  # Includes readr and dplyr for import and wrangling


# ── 2. Load Taxonomy Tables ───────────────────────────────────────────────────
# Both files are SILVA nr99 v138.2 taxonomy assignments to species level,
# generated from two separate training set FASTA files (fa1 and fa2).

tax1 <- read.table(
  "../data/ASV_tables/08_assigned_taxonomy_silva_nr99_v138.2_toSpecies_trainset.fa1.csv",
  row.names = 1
)

tax2 <- read.table(
  "../data/ASV_tables/08_assigned_taxonomy_silva_nr99_v138.2_toSpecies_trainset.fa2.csv",
  row.names = 1
)


# ── 3. Combine & Deduplicate ──────────────────────────────────────────────────
# Row-bind both tables then remove any exactly duplicated rows.
# NOTE: If ASV IDs (row names) conflict between tax1 and tax2, consider
#       investigating with: intersect(rownames(tax1), rownames(tax2))

tax_combined <- rbind(tax1, tax2)
tax_total    <- unique(tax_combined)


# ── 4. Export Reference Taxonomy ──────────────────────────────────────────────
# row.names = TRUE preserves ASV IDs as the first column in the output file,
# which is required for downstream matching against ASV abundance tables.

write.csv(
  tax_total,
  "../data/tax_total.csv",
  row.names = TRUE
)