# =============================================================================
# Date:    26/05/2026
# Authors: Katy Maher & Agnes Sz.
# Purpose: Combine two corrected ASV count tables from separate Illumina runs
#          into a single merged output file.
# =============================================================================

rm(list = ls())
library(dplyr)


# =============================================================================
# SECTION 1: Load Corrected ASV Count Tables
# =============================================================================

counts1 <- read.table("../data/06_ASV_counts1corrected.tsv", row.names = NULL)  # Run 1 corrected
counts2 <- read.table("../data/06_ASV_counts2corrected.tsv", row.names = NULL)  # Run 2 corrected


# =============================================================================
# SECTION 2: Harmonise Column Names Across Runs
# =============================================================================

# Run 1 columns follow the pattern "X1_<id>.<suffix>"; reformat to "X<id>.<suffix>"
# Run 2 columns follow the pattern "X<id>_<suffix>";   reformat to "X<id>.<suffix>"
colnames(counts1) <- gsub("^X1_(\\d+\\..+)", "X\\1",  colnames(counts1))
colnames(counts2) <- gsub("^X(\\d+)_",       "X\\1.", colnames(counts2))


# =============================================================================
# SECTION 3: Merge Runs and Sum Counts
# =============================================================================

# Merge both runs by ASV ID (row.names), summing read counts across runs
counts_total <- aggregate(. ~ row.names, rbind(counts1, counts2), sum)


# =============================================================================
# SECTION 4: Standardise Column Names to Match Metadata
# =============================================================================

# Work on a copy to preserve the raw merged table
ASV_table_cleaned <- counts_total

# --- Helper function: rename columns to match metadata sample IDs ---
# Rules:
#   1. PCR control columns are kept as-is
#   2. Two specific columns with duplicate sample IDs are renamed manually
#   3. All other columns have their "X<plate>." prefix stripped
#   4. R-safe names are enforced, with leading "X" replaced by "S"
rename_asv_columns <- function(asv_table) {
  cols <- colnames(asv_table)
  
  new_names <- ifelse(
    grepl("PCR", cols),       cols,     # Rule 1: keep PCR controls unchanged
    ifelse(cols == "X304.70", "S70.2",  # Rule 2a: manual rename for duplicate sample
           ifelse(cols == "X312.71", "S71.2",  # Rule 2b: manual rename for duplicate sample
                  gsub("^X[^.]*\\.", "", cols) # Rule 3: strip plate prefix
           ))
  )
  
  new_names <- trimws(new_names)                     # Remove accidental whitespace
  new_names <- make.names(new_names, unique = TRUE)  # Enforce R-valid, unique names
  new_names <- gsub("^X", "S", new_names)            # Replace leading "X" with "S"
  
  return(new_names)
}

# Apply renaming
colnames(ASV_table_cleaned) <- rename_asv_columns(counts_total)


# =============================================================================
# SECTION 5: Save Output
# =============================================================================

# NOTE: row.names = FALSE is required for correct phyloseq import
write.csv(ASV_table_cleaned, "../data/06_ASV_counts_both_runs_corrected.csv", row.names = FALSE, quote = FALSE)