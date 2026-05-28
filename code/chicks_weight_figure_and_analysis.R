# =============================================================================
# Script:  Chick Weight ~ Hatching Date Quadratic Model
# Author:  Agnes
# Date:    11/05/2025
# Aim:     Fit a quadratic linear mixed-effects model of standardised chick
#          weight against standardised hatching date, and visualise the
#          relationship with a fitted curve.
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
# 3. DATA CLEANING
# =============================================================================

# Remove rows with missing values in variables required for modelling
combined_chick_mass <- combined_chick_mass %>%
  filter(
    !is.na(Weight),
    !is.na(Standardised_HD),
    !is.na(Year),
    !is.na(NestBoxID),
    !is.na(BroodID)
  )

# =============================================================================
# 4. EXPLORATORY CHECKS
# =============================================================================
# Visual checks for approximate normality and the raw relationship of interest

hist(
  combined_chick_mass$Standardised_W,
  breaks = 30,
  main   = "Standardised weight distribution",
  xlab   = "Standardised weight (SD)"
)

plot(
  combined_chick_mass$Standardised_HD,
  combined_chick_mass$Standardised_W,
  main = "Standardised weight vs standardised hatching date",
  xlab = "Standardised hatching date (SD)",
  ylab = "Standardised weight (SD)"
)

# =============================================================================
# 5. QUADRATIC MIXED-EFFECTS MODEL
# =============================================================================

# -----------------------------------------------------------------------------
# 5a. Fit the model
#     Fixed effects : Standardised_HD (linear) + Standardised_HD^2 (quadratic)
#     Random effects: Year (inter-annual variation),
#                     NestBoxID (nest site variation),
#                     BroodID (brood-level variation)
# -----------------------------------------------------------------------------

weight_model_sq <- lmer(
  Standardised_W ~ Standardised_HD + I(Standardised_HD^2) +
    (1 | Year) + (1 | NestBoxID) + (1 | BroodID),
  data = combined_chick_mass
)

summary(weight_model_sq)

# Residual check — should be approximately normally distributed
hist(
  resid(weight_model_sq),
  main = "Quadratic weight model residuals",
  xlab = "residuals"
)

# =============================================================================
# 6. VISUALISATION
# =============================================================================

# -----------------------------------------------------------------------------
# 6a. Extract fixed-effect coefficients from the model
#     These are used to define the fitted curve plotted over the raw data.
#     Using fixef() ensures the coefficients always stay in sync with the
#     model, rather than relying on hard-coded values.
# -----------------------------------------------------------------------------

model_coefs <- fixef(weight_model_sq)

b0 <- model_coefs["(Intercept)"]          # Intercept
b1 <- model_coefs["Standardised_HD"]      # Linear term
b2 <- model_coefs["I(Standardised_HD^2)"] # Quadratic term

# -----------------------------------------------------------------------------
# 6b. Define the fitted curve function (fixed effects only)
# -----------------------------------------------------------------------------

fitted_curve <- function(x) {
  b0 + b1 * x + b2 * x^2
}

# -----------------------------------------------------------------------------
# 6c. Build the plot
# -----------------------------------------------------------------------------

p1 <- ggplot(combined_chick_mass, aes(x = Standardised_HD, y = Standardised_W)) +
  geom_point(alpha = 0.4, color = "#174c56") +
  stat_function(
    fun      = fitted_curve,
    color    = "black",
    linewidth = 1.2,
    linetype = "dashed"
  ) +
  xlim(-6, 6) +
  ylim(-1, 1) +
  theme_tufte(base_size = 18) +
  labs(
    x = "Standardised hatching date (SD)",
    y = "Standardised chick weight (SD)"
  )

p1

# -----------------------------------------------------------------------------
# 6d. Save figure as SVG for publication-quality vector output
# -----------------------------------------------------------------------------

ggsave(
  filename = "../figures/weight_vs_hatching_date_quadratic.svg",
  plot     = p1,
  device   = "svg",
  width    = 8,
  height   = 6,
  units    = "in"
)


# =============================================================================
# 7. (OPTIONAL) MODEL OUTPUT TABLE
# =============================================================================
# Uncomment to generate an HTML summary table of the model output

# tab1 <- tab_model(
#   weight_model_sq,
#   show.ci     = 0.95,
#   show.se     = TRUE,
#   show.stat   = TRUE,
#   show.df     = TRUE,
#   show.re.var = TRUE,
#   digits      = 4
# )