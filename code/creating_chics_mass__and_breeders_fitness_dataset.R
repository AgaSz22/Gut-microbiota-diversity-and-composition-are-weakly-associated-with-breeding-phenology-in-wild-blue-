# =============================================================================
# Author:  Agnes Szwarczynska
# Purpose: Create a combined dataset (2013 & 2024) to examine:
#            1. Trends in egg count & fledgling number vs. standardised
#               laying date
#            2. Chick mass vs. standardised hatching date
# Inputs:  - data/Chicks_2013.csv
#          - data/Chicks_2024.csv
#          - data/breeding_2013.csv
#          - data/breeding_2024.xlsx
# Outputs: - data/chicks_2013_mass_all.csv
#          - data/chicks_2024_mass_all.csv
#          - data/combined_chicks_dataset.csv
# =============================================================================


# =============================================================================
# 1. PACKAGES
# =============================================================================

# NOTE: tidyverse already includes dplyr, ggplot2, and readr,
#       so those three packages do not need to be loaded separately.
library(tidyverse)  # Includes: dplyr, ggplot2, readr, tidyr, etc.
library(readxl)     # For reading Excel files (.xlsx)


# =============================================================================
# 2. LOAD DATA
# =============================================================================

chicks_2013   <- read_csv("../data/Chicks_2013.csv")
chicks_2024   <- read_csv("../data/Chicks_2024.csv")
breeding_2013 <- read_csv("../data/breeding_2013.csv")
breeding_2024 <- read_excel("../data/breeding_2024.xlsx")


# =============================================================================
# 3. PROCESS 2013 DATA
# =============================================================================

# --- 3a. Filter chick data ---------------------------------------------------
# Filters applied:
#   - Blue Tit species only
#   - Valid weight: 0 < Weight < 13  
#   - Exclude manipulated nests

chicks_2013_mass <- chicks_2013[c(2, 3, 4, 5, 6, 10, 11, 12)] %>%
  dplyr::filter(
    Species == "BLUTI",
    !is.na(Weight) & Weight > 0 & Weight < 13,
    is.na(Manipulation) 
  )

# --- 3b. Filter breeding data ------------------------------------------------
# Filters applied:
#   - Used nestboxes only
#   - Blue Tit species only
#   - Remove clutch size outliers (ClutchSize < 20)
#   - Exclude manipulated nests

breeding_2013_mass <- breeding_2013[-c(1, 12)] %>%
  dplyr::filter(
    Nestbox_Fate == "used",
    Species      == "Blue Tit",
    ClutchSize   < 20,
    is.na(ManipulationID) 
  )

# --- 3c. Merge chick and breeding data, then drop redundant columns ----------

chicks_2013_mass_all <- merge(chicks_2013_mass, breeding_2013_mass, by = "BroodID") %>%
  select(-c(2, 3, 4))

# Rename columns for clarity and cross-year consistency
colnames(chicks_2013_mass_all)[6]  <- "Year"
colnames(chicks_2013_mass_all)[9]  <- "Species"
colnames(chicks_2013_mass_all)[12] <- "Hatching_Date"

# --- 3d. Save 2013 processed data --------------------------------------------

write.csv(chicks_2013_mass_all, "../data/chicks_2013_mass_all.csv", row.names = FALSE)


# =============================================================================
# 4. PROCESS 2024 DATA
# =============================================================================

# --- 4a. Filter chick data ---------------------------------------------------
# Note: The date column originally records the 14-day-old measurement date,
#       so we subtract 14 to recover the true hatching date.

chicks_2024_mass <- chicks_2024[c(2, 3, 4, 7, 11, 17, 21)]

# Rename columns for clarity before filtering
colnames(chicks_2024_mass)[5] <- "Hatching_Date"
colnames(chicks_2024_mass)[6] <- "Weight"
colnames(chicks_2024_mass)[7] <- "Manipulation"

chicks_2024_mass <- chicks_2024_mass %>%
  mutate(
    # Correct hatching date: the recorded date was at 14 days old
    Hatching_Date = Hatching_Date - 14
  ) %>%
  dplyr::filter(
    is.na(Manipulation) | Manipulation == "none",
    !is.na(Weight) & Weight > 0 & Weight < 13           
  )

# --- 4b. Filter breeding data ------------------------------------------------
# Filters applied:
#   - Used nestboxes only
#   - Blue Tit species only
#   - Unmanipulated nests only
#   - Remove fledgling number outliers (Fledged_ChickNumber < 18)
#   - Remove clutch size outliers (ClutchSize < 20)

breeding_2024_mass <- breeding_2024[c(2, 3, 4, 5, 6, 7, 10, 23, 25)] %>%
  mutate(
    # Convert to numeric BEFORE filtering, otherwise comparisons are string-based
    # and will produce incorrect results (e.g. "9" < "20" is FALSE as a string)
    ClutchSize          = as.numeric(ClutchSize),
    Fledged_ChickNumber = as.numeric(Fledged_ChickNumber)
  ) %>%
  dplyr::filter(
    Nestbox_Fate          == "used",
    Species               == "Blue Tit",
    is.na(ManipulationID) | ManipulationID == "none",
    Fledged_ChickNumber   < 19,
    ClutchSize            < 20
  )

# --- 4c. Merge chick and breeding data, then drop redundant columns ----------

chicks_2024_mass_all <- merge(chicks_2024_mass, breeding_2024_mass, by = "BroodID") %>%
  select(-c(2, 3))

# Rename columns for clarity and cross-year consistency
colnames(chicks_2024_mass_all)[6] <- "Year"
colnames(chicks_2024_mass_all)[7] <- "NestBoxID"

# --- 4d. Save 2024 processed data --------------------------------------------

write.csv(chicks_2024_mass_all, "../data/chicks_2024_mass_all.csv", row.names = FALSE)


# =============================================================================
# 5. COMBINE 2013 AND 2024 DATASETS
# =============================================================================

# Align the NestBoxID column name in 2013 data to match 2024 before binding
colnames(chicks_2013_mass_all)[7] <- "NestBoxID"

# Ensure LayingDate is numeric in the 2024 dataset for compatibility with 2013
chicks_2024_mass_all <- chicks_2024_mass_all %>%
  mutate(LayingDate = as.numeric(LayingDate))

# Row-bind both years into one combined dataset
combined_chick_mass <- bind_rows(chicks_2024_mass_all, chicks_2013_mass_all)

# --- Exclude second broods ---------------------------------------------------
# Keep only the first occurrence of each BroodID within a given year.
# Since each chick has its own row, we first identify the first BroodID per
# NestBoxID per year (i.e. the earliest LayingDate), then filter to keep only
# rows belonging to those first broods.

combined_chick_mass <- combined_chick_mass %>%
  filter(!is.na(LayingDate)) %>%                               # Exclude rows with no laying date
  group_by(Year, NestBoxID) %>%
  mutate(First_BroodID = BroodID[which.min(LayingDate)]) %>%  # Identify the first brood per nestbox per year
  ungroup() %>%
  filter(BroodID == First_BroodID) %>%                        # Keep only rows belonging to that first brood
  select(-First_BroodID)                                      # Drop the helper column                                     # Drop the helper column


# =============================================================================
# 6. STANDARDISE LAYING DATE, HATCHING DATE, AND WEIGHT
# =============================================================================
# All standardisations are performed within each year group.
# NA values in the date columns will produce NA in the standardised columns
# automatically — no need for a manual for-loop.

combined_chick_mass <- combined_chick_mass %>%
  group_by(Year) %>%
  mutate(
    # --- Standardised Laying Date ---
    Mean_LD         = mean(LayingDate, na.rm = TRUE),
    SD_LD           = sd(LayingDate,   na.rm = TRUE),
    Standardised_LD = (LayingDate - Mean_LD) / SD_LD,
    
    # --- Standardised Hatching Date ---
    Mean_HD         = mean(Hatching_Date, na.rm = TRUE),
    SD_HD           = sd(Hatching_Date,   na.rm = TRUE),
    Standardised_HD = (Hatching_Date - Mean_HD) / SD_HD,
    
    # --- Standardised Weight ---
    # Keeping your original scaling: (Weight - Mean_W) / Mean_W
    Mean_W          = mean(Weight, na.rm = TRUE),
    Standardised_W  = (Weight - Mean_W) / Mean_W
  ) %>%
  ungroup()


# =============================================================================
# 7. SAVE FINAL COMBINED DATASET
# =============================================================================

write.csv(combined_chick_mass, "../data/combined_chicks_dataset.csv", row.names = FALSE)
