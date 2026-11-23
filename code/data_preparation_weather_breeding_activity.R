# =============================================================================
# Script:  data_preparation_weather_breeding_activity.R
# Purpose: Prepare and export breeding distribution and weather data for the
#          2024 breeding season analysis.
# Author:  Agnes Szwarczynska
# Date:    2026-05-27
# =============================================================================


# ── 1. Libraries ──────────────────────────────────────────────────────────────

library(tidyverse)  # Includes dplyr, ggplot2, readr, tidyr, etc.
library(readxl)     # Reading Excel files


# ── 2. Helper Function: Compute Ordinal Day ───────────────────────────────────
# Converts Month + Day into an ordinal day relative to April 1st (Day 1),
# covering the core breeding season months (March–June).
#
# @param month  Integer month (3 = March, 4 = April, 5 = May, 6 = June)
# @param day    Integer day of the month
# @return       Integer ordinal day relative to April 1st, or NA

compute_ordinal_ld <- function(month, day) {
  case_when(
    month == 3 ~ day - 31,  # March:  Day -31 to -1 relative to April 1st
    month == 4 ~ day,       # April:  Day   1 to 30
    month == 5 ~ day + 30,  # May:    Day  31 to 61
    month == 6 ~ day + 61,  # June:   Day  62 to 91
    TRUE       ~ NA_real_   # Outside expected range → NA
  )
}


# ── 3. Load & Prepare Breeding Data ───────────────────────────────────────────

breeding_2024 <- read_excel("../data/breeding_2024.xlsx")

breeding_2024_distribution <- breeding_2024 %>%
  # Select only the relevant columns by name (more robust than column indices)
  select(2:7, 18, 25) %>%
  # Retain only used Blue Tit nestboxes from 2024
  filter(
    Nestbox_Fate == "used",
    Species      == "Blue Tit",
    Year         == 2024
  ) %>%
  drop_na() %>%
  # Count distinct broods per laying date
  group_by(LayingDate) %>%
  mutate(No_Nests = n_distinct(BroodID)) %>%
  ungroup() %>%
  # Count distinct broods reaching day 14 per date
  group_by(Day14_Date) %>%
  mutate(No_14Hatchings = n_distinct(BroodID)) %>%
  ungroup()

write.csv(
  breeding_2024_distribution,
  "../data/breeding_2024_distribution.csv",
  row.names = FALSE
)


# ── 4. Load & Prepare Weather Data ────────────────────────────────────────────

weather_2024 <- read_csv("../data/SilwoodWeatherDaily.csv") %>%
  # Strip the time component from the TIMESTAMP, keeping date only (DD/MM/YYYY)
  mutate(TIMESTAMP = sub(" .*", "", TIMESTAMP)) %>%
  rename(Date = TIMESTAMP) %>%
  # Split date string into numeric Day, Month, Year columns
  separate_wider_delim(
    Date,
    delim = "/",
    names = c("Day", "Month", "Year")
  ) %>%
  mutate(
    across(c(Day, Month, Year), as.numeric)
  ) %>%
  # Retain 2024 records only
  filter(Year == 2024) %>%
  # Compute ordinal day relative to April 1st
  mutate(Ordinal_LD = compute_ordinal_ld(Month, Day))

write.csv(
  weather_2024,
  "../data/weather_chapter_2.csv",
  row.names = FALSE
)
