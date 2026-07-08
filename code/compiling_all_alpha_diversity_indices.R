
# ==============================================================================
# Date:    08/06/2026
# Author:  Agnes Szwarczynska
# Purpose: Combine all alpha diversity indices and metadata
# ==============================================================================

# --- Libraries ----------------------------------------------------------------

# library(decontam)  # Uncomment if decontamination step is needed
library(dplyr)

# --- Load Data ----------------------------------------------------------------

alpha_df_mean <- readRDS("../rds_objects/alpha_df_mean.rds")
faith_pd_df   <- readRDS("../rds_objects/faith_pd_df.rds")

# --- Preprocessing ------------------------------------------------------------

# Standardise Sample_ID column name in Faith's PD data frame
colnames(faith_pd_df)[1] <- "Sample_ID"

# Promote row names to a Sample_ID column in alpha diversity data frame
alpha_df_mean$Sample_ID <- rownames(alpha_df_mean)

# Select relevant diversity index columns alongside Sample_ID
diversity_indices <- alpha_df_mean[, c(1, 3, 5, 23)]

# Define samples to exclude
samples_to_exclude <- c("Extraction.Negative.8")

# Remove negative control sample from both data frames
diversity_indices <- diversity_indices %>%
  filter(!Sample_ID %in% samples_to_exclude)

faith_pd_df <- faith_pd_df %>%
  filter(!Sample_ID %in% samples_to_exclude)

# --- Load and Clean Metadata --------------------------------------------------

microbial_metadata <- read.csv("../data/Microbial_metadata_cleaned.csv", row.names = 1)

# Impute missing values with sensible defaults
microbial_metadata$Stage[is.na(microbial_metadata$Stage)]               <- "Negative"
microbial_metadata$Capture_Date[is.na(microbial_metadata$Capture_Date)] <- 70

# Promote row names to a Sample_ID column
microbial_metadata$Sample_ID <- rownames(microbial_metadata)

# Initialise Relative_to_frass column
microbial_metadata$Relative_to_frass <- NA_real_

# --- Merge Data Frames --------------------------------------------------------

all_diversity_indices <- diversity_indices %>%
  inner_join(faith_pd_df,       by = "Sample_ID") %>%
  inner_join(microbial_metadata, by = "Sample_ID")

# --- Feature Engineering ------------------------------------------------------

# Calculate capture date relative to frass (offset = 47 days)
# Note: vectorised operation used in place of a for-loop for efficiency
all_diversity_indices$Relative_to_frass <- all_diversity_indices$Capture_Date - 4 - 43

# --- Export -------------------------------------------------------------------

write.csv(all_diversity_indices, "../data/all_diversity_indices.csv", row.names = FALSE)