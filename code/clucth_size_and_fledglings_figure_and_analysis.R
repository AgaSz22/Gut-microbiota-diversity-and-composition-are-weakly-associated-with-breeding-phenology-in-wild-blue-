# =============================================================================
# Script:  Clutch Size & Fledgling Analysis
# Author:  Agnes Szwarczynska
# Date:    11/05/2025
# Aim:     Model the relationship between standardised laying date and (1) clutch
#          size and (2) number of fledglings, using linear mixed-effects models
#          with NestBoxID and Year as random effects.
# =============================================================================


# =============================================================================
# 1. LIBRARIES
# =============================================================================

library(readr)       # CSV import
library(dplyr)       # Data manipulation
library(tidyverse)   # Core data science toolkit
library(lme4)        # Linear mixed-effects models
library(lmerTest)    # p-values for lme4 models
library(ggplot2)     # Plotting
library(cowplot)     # Multi-panel plot layout
library(ggthemes)    # ggplot2 themes (Tufte)
library(sjPlot)      # Model output tables


# =============================================================================
# 2. LOAD DATA
# =============================================================================

combined_chick_mass <- read.csv("../data/combined_chicks_dataset.csv")

# =============================================================================
# 3. EXPLORATORY CHECKS
# =============================================================================
# Quick visual check that key variables are roughly normally distributed

hist(combined_chick_mass$ClutchSize,          breaks = 30, main = "Clutch size distribution")
hist(combined_chick_mass$Fledged_ChickNumber,  breaks = 30, main = "Fledged chick number distribution")

# =============================================================================
# 4. CLUTCH SIZE ANALYSIS
# =============================================================================

# -----------------------------------------------------------------------------
# 4a. Aggregate to brood level
#     Each row becomes one brood; continuous variables are averaged across
#     chick-level rows, and categorical/ID variables take the first value.
# -----------------------------------------------------------------------------

brood_level_data_CS <- combined_chick_mass %>%
  group_by(BroodID) %>%
  summarise(
    NestBoxID      = first(NestBoxID),
    ClutchSize     = mean(ClutchSize,      na.rm = TRUE),
    Standardised_LD = mean(Standardised_LD, na.rm = TRUE),
    Year           = first(Year)
  ) %>%
  filter(
    !is.na(ClutchSize),
    !is.na(Standardised_LD),
    !is.na(Year)
  )

# -----------------------------------------------------------------------------
# 4b. Fit linear mixed-effects model
#     Fixed effect : Standardised laying date
#     Random effects: Year (inter-annual variation), NestBoxID (site variation)
# -----------------------------------------------------------------------------

clutch_model <- lmer(
  ClutchSize ~ Standardised_LD + (1 | Year) + (1 | NestBoxID),
  data = brood_level_data_CS
)

summary(clutch_model)
saveRDS(clutch_model, "../rds_objects/clutch_model.rds")

# Residual check — should be approximately normally distributed
hist(resid(clutch_model), main = "Clutch Model Residuals", xlab = "Residuals")

# -----------------------------------------------------------------------------
# 4c. Generate model-predicted (fixed-effects only) fitted values
# -----------------------------------------------------------------------------

brood_level_data_CS$Fitted <- predict(clutch_model, re.form = NA)

# -----------------------------------------------------------------------------
# 4d. Plot: Clutch Size ~ Standardised Laying Date
# -----------------------------------------------------------------------------

p1 <- ggplot(brood_level_data_CS, aes(x = Standardised_LD, y = ClutchSize)) +
  geom_point(color = "#AA336A", alpha = 0.4) +
  geom_line(
    aes(y = Fitted, group = 1),
    color    = "black",
    linewidth = 1.3,          # NOTE: 'size' deprecated in ggplot2 >= 3.4.0; use 'linewidth'
    linetype = "dashed"
  ) +
  xlim(-6, 6) +
  ylim(0, 20) +
  theme_tufte(base_size = 18) +
  labs(
    x = "Standardised egg-laying day (SD)",
    y = "Clutch size"
  )

# -----------------------------------------------------------------------------
# 4e. (Optional) Model output table — uncomment to generate HTML table
# -----------------------------------------------------------------------------

# tab1 <- tab_model(
#   clutch_model,
#   show.ci   = 0.95,
#   show.se   = TRUE,
#   show.stat = TRUE,
#   show.df   = TRUE,
#   show.re.var = TRUE
# )


# =============================================================================
# 5. FLEDGLING NUMBER ANALYSIS
# =============================================================================

# -----------------------------------------------------------------------------
# 5a. Compute standardised laying date within each year
#     (mean-centred and scaled by SD, calculated per year)
#     NOTE: This block standardises LayingDate for the combined_chick_mass
#     dataset, which does not yet contain Standardised_LD for fledgling rows.
# -----------------------------------------------------------------------------

combined_chick_mass <- combined_chick_mass %>%
  group_by(Year) %>%
  mutate(
    Mean_LD = mean(LayingDate, na.rm = TRUE),
    SD_LD   = sd(LayingDate,   na.rm = TRUE)
  ) %>%
  ungroup()

# Vectorised standardisation — replaces the original for-loop for clarity and speed.
# Rows with a missing LayingDate are assigned NA (handled downstream by filter()).
combined_chick_mass <- combined_chick_mass %>%
  mutate(
    Standardised_LD = if_else(
      is.na(LayingDate),
      NA_real_,
      (LayingDate - Mean_LD) / SD_LD
    )
  )

# -----------------------------------------------------------------------------
# 5b. Aggregate to brood level
# -----------------------------------------------------------------------------

brood_level_data_fledglings <- combined_chick_mass %>%
  group_by(BroodID) %>%
  summarise(
    NestBoxID          = first(NestBoxID),
    ClutchSize         = mean(ClutchSize,          na.rm = TRUE),
    Fledged_ChickNumber = mean(Fledged_ChickNumber, na.rm = TRUE),
    Standardised_LD    = mean(Standardised_LD,     na.rm = TRUE),
    Year               = first(Year)
  ) %>%
  filter(
    !is.na(ClutchSize),
    !is.na(Fledged_ChickNumber),
    !is.na(Standardised_LD),
    !is.na(Year)
  )

# -----------------------------------------------------------------------------
# 5c. Fit linear mixed-effects model
#     Fixed effect : Standardised laying date
#     Random effects: Year, NestBoxID
# -----------------------------------------------------------------------------

fledglings_model <- lmer(
  Fledged_ChickNumber ~ Standardised_LD + (1 | Year) + (1 | NestBoxID),
  data = brood_level_data_fledglings
)

summary(fledglings_model)
saveRDS(fledglings_model, "../rds_objects/fledglings_model.rds")

# Residual check
hist(resid(fledglings_model), main = "Fledglings model residuals", xlab = "Residuals")

# -----------------------------------------------------------------------------
# 5d. Generate fixed-effects fitted values
# -----------------------------------------------------------------------------

brood_level_data_fledglings$Fitted <- predict(fledglings_model, re.form = NA)

# -----------------------------------------------------------------------------
# 5e. Plot: Fledgling Number ~ Standardised Laying Date
# -----------------------------------------------------------------------------

p2 <- ggplot(brood_level_data_fledglings, aes(x = Standardised_LD, y = Fledged_ChickNumber)) +
  geom_point(color = "#F56A6A", alpha = 0.4) +
  geom_line(
    aes(y = Fitted, group = 1),
    color     = "black",
    linewidth  = 1.3,
    linetype  = "dashed"
  ) +
  xlim(-6, 6) +
  ylim(0, 18) +
  theme_tufte(base_size = 18) +
  labs(
    x = "Standardised egg-laying day (SD)",
    y = "Fledglings"
  )

# -----------------------------------------------------------------------------
# (Optional) Model output table — uncomment to generate HTML table
# -----------------------------------------------------------------------------

# tab2 <- tab_model(
#   fledglings_model,
#   show.ci     = 0.95,
#   show.se     = TRUE,
#   show.stat   = TRUE,
#   show.df     = TRUE,
#   show.re.var = TRUE
# )


# =============================================================================
# 6. COMBINED FIGURE
# =============================================================================

# Arrange p2 (Clutch Size) and p3 (Fledglings) in a single column with labels
combined_plot <- cowplot::plot_grid(
  p1, p2,
  labels = c("A", "B"),
  ncol   = 1
)

combined_plot

# Save the combined figure as an SVG for publication-quality vector output
ggsave(
  filename = "../figures/clutch_fledgling_combined.svg",
  plot     = combined_plot,
  device   = "svg",
  width    = 8,
  height   = 7,
  units    = "in"
)
