# =============================================================================
# Script:  frass_breeding_weather_analysis.R
# Purpose: Visualise frass weight, nest activity, and weather data across
#          the 2024 breeding season using a shared ordinal day axis.
# Author:  Agnes Szwarczynska
# Date:    2026-05-27
# =============================================================================


# ── 1. Libraries ──────────────────────────────────────────────────────────────
# NOTE: ggplot2 and ggthemes are already included via tidyverse/explicit loads.
#       Removed duplicate library() calls.

library(tidyverse)   # Includes dplyr, ggplot2, readr, tidyr, etc.
library(readxl)      # Reading Excel files
library(ggthemes)    # Tufte theme and related extensions
library(cowplot)     # Combining multiple ggplots into one figure

# ── 2. Helper Function: Compute Ordinal Day ───────────────────────────────────
# Converts Month + Day into an ordinal day relative to April 1st (Day 0),
# covering the core breeding season months (March–June).
#
# Returns NA for any month outside the expected range.
#
# @param month  Integer month (3 = March, 4 = April, 5 = May, 6 = June)
# @param day    Integer day of the month
# @return       Integer ordinal day relative to April 1st, or NA

compute_ordinal_ld <- function(month, day) {
  case_when(
    month == 3 ~ day - 31,   # March:  Day -31 to -1 relative to April 1st
    month == 4 ~ day,        # April:  Day  1  to 30
    month == 5 ~ day + 30,   # May:    Day  31 to 61
    month == 6 ~ day + 61,   # June:   Day  62 to 91
    TRUE       ~ NA_real_    # Outside expected range → NA
  )
}


# ── 3. Load & Prepare Frass Data ──────────────────────────────────────────────

frass_2024 <- read_excel("../data/Frass_traps.xlsx")

frass_2024 <- frass_2024 %>%
  # Split ISO-format date string (YYYY-MM-DD) into separate numeric columns
  separate_wider_delim(
    Collection_date,
    delim  = "-",
    names  = c("Year", "Month", "Day")
  ) %>%
  mutate(
    across(c(Year, Month, Day), as.numeric),
    Ordinal_LD = compute_ordinal_ld(Month, Day)
  )


# ── 4. Summarise Frass Data ───────────────────────────────────────────────────

frass_mean <- frass_2024 %>%
  group_by(Ordinal_LD) %>%
  summarise(Mean_Frass_Weight = mean(Frass_Weight, na.rm = TRUE), .groups = "drop")

# Identify the ordinal day on which mean frass weight peaks
frass_peak_day <- frass_mean %>%
  filter(Mean_Frass_Weight == max(Mean_Frass_Weight, na.rm = TRUE)) %>%
  pull(Ordinal_LD)


# ── 5. Plot A — Mean Frass Weight ─────────────────────────────────────────────

plot_frass <- ggplot(frass_mean, aes(x = Ordinal_LD, y = Mean_Frass_Weight)) +
  geom_line(color = "#1AB788", linewidth = 1) +
  geom_vline(
    xintercept = frass_peak_day,
    linetype   = "dashed",
    color      = "grey50",
    linewidth  = 1,
    alpha      = 0.7
  ) +
  scale_x_continuous(
    name   = "Day of the breeding season",
    limits = c(-3, 65)
  ) +
  scale_y_continuous(
    name   = "Frass weight (g)",
    limits = c(0, 3)
  ) +
  theme_tufte(base_size = 14, ticks = TRUE)


# ── 6. Load & Prepare Weather Data ───────────────────────────────────────────
# Column names contain spaces and special characters; backtick-quote when used.

weather_2024 <- read_csv("../data/weather_chapter_2.csv") %>%
  filter(Ordinal_LD > -4 & Ordinal_LD < 66)

# Scale factor maps precipitation onto the temperature axis for dual-axis plot
temp_precip_scale <- max(weather_2024$`Air_Temp (Deg C) [Smp]`,    na.rm = TRUE) /
  max(weather_2024$`Rain_mm_Tot  (mm) [Tot]`, na.rm = TRUE)


# ── 7. Plot B — Temperature & Precipitation ───────────────────────────────────

plot_weather <- ggplot(weather_2024, aes(x = Ordinal_LD)) +
  
  # Air temperature (primary / left axis)
  geom_line(
    aes(y = `Air_Temp (Deg C) [Smp]`),
    color     = "#ffdc5d",
    linewidth = 1
  ) +
  
  # Precipitation scaled to temperature range (secondary / right axis)
  geom_line(
    aes(y = `Rain_mm_Tot  (mm) [Tot]` * temp_precip_scale),
    color     = "#00ffff",
    linewidth = 1
  ) +
  
  scale_y_continuous(
    name     = "Air temperature (°C)",
    sec.axis = sec_axis(
      transform = ~ . / temp_precip_scale,
      name      = "Precipitation (mm)"
    )
  ) +
  
  scale_x_continuous(name = "Day of the breeding season") +
  coord_cartesian(xlim = c(-3, 65), ylim = c(0, 25)) +
  
  theme_tufte(base_size = 14, ticks = TRUE) +
  theme(
    axis.title.y.left  = element_text(color = "black"),
    axis.title.y.right = element_text(color = "black"),
    axis.line          = element_blank()
  )


# ── 8. Load & Prepare Breeding Distribution Data ──────────────────────────────

breeding_2024 <- read_csv("../data/breeding_2024_distribution.csv")

# Count unique broods per laying date
laying_counts <- breeding_2024 %>%
  group_by(LayingDate) %>%
  summarise(No_Nests_Laid = n_distinct(BroodID), .groups = "drop") %>%
  rename(Date = LayingDate)

# Count unique broods reaching day-14 (chick stage) per date
chick_counts <- breeding_2024 %>%
  group_by(Day14_Date) %>%
  summarise(No_Nests_Chicks = n_distinct(BroodID), .groups = "drop") %>%
  rename(Date = Day14_Date)

# Merge both date series, keeping all dates from either source
nest_activity <- full_join(laying_counts, chick_counts, by = "Date") %>%
  mutate(Date = as.numeric(Date)) %>%
  filter(!is.na(Date)) %>%
  arrange(Date)


# ── 9. Plot C — Nest Activity ─────────────────────────────────────────────────

plot_nests <- ggplot(nest_activity, aes(x = Date)) +
  geom_line(
    aes(y = No_Nests_Laid),
    color     = "#074151",
    linewidth = 1
  ) +
  geom_line(
    aes(y = No_Nests_Chicks),
    color     = "#A3334E",
    linewidth = 1
  ) +
  scale_x_continuous(
    name   = "Day of the breeding season",
    limits = c(-3, 65)
  ) +
  ylab("Number of nests") +
  theme_tufte(base_size = 14, ticks = TRUE)


# ── 10. Combine Plots ─────────────────────────────────────────────────────────

combined_plot <- plot_grid(
  plot_frass,    # Panel A — Frass weight over season
  plot_nests,    # Panel B — Nest activity (laying vs. chick stage)
  plot_weather,  # Panel C — Temperature and precipitation
  labels  = c("A", "B", "C"),
  ncol    = 1,    # Stack panels vertically
  align   = "v",  # Align vertically
  axis    = "lr"  # Match left and right axis edges
)

# ── 11. Export Combined Plot as SVG ──────────────────────────────────────────
# Target size: 1/3 of A4 page height, full A4 width
# A4 dimensions: 210mm wide × 297mm tall  →  1/3 height ≈ 99mm

ggsave(
  filename = "../figures/frass_breeding_weather_combined.svg",
  plot     = combined_plot,
  device   = "svg",
  width    = 210,   # Full A4 width in mm
  height   = 230,    # 1/3 A4 height in mm
  units    = "mm"
)

print(combined_plot)


dev.off()

