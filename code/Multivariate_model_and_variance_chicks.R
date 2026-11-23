# =============================================================================
# Date:    09/06/2026
# Author:  Agnes Szwarczynska
# Purpose: Run a multivariate brms model and perform variance partitioning
# =============================================================================


# =============================================================================
# 1. LIBRARIES
# =============================================================================

library(ggplot2)
library(dplyr)
library(MCMCglmm)
library(brms)
library(cowplot)
library(sjPlot)
library(bayesplot)
library(rstan)
library(tidybayes)


# =============================================================================
# 2. DATA LOADING & PREPARATION
# =============================================================================

all_diversity_indices <- read.csv("../data/all_diversity_indices.csv")

# Define colour palette for nest identification
color_palette <- c(
  "#1AB788",
  "#ffdc5d",
  "#074151",
  "#F56A6A",
  "#A3334E",
  "#00ffff",
  "#B39EB5",
  "#FFC5D3",
  "#330C15"
)

# Assign colour indices to nests per sampling date
assign_colors_per_date <- function(data, date_col, nest_col) {
  data %>%
    group_by(across(all_of(date_col))) %>%
    mutate(
      NestID  = as.integer(factor(!!sym(nest_col))),
      ColorID = (NestID - 1) %% length(color_palette) + 1,
      ColorHex = color_palette[ColorID]
    ) %>%
    ungroup()
}

# Filter to chick stage and assign colours
alpha_diversity_chicks <- all_diversity_indices %>%
  filter(Stage == "Chick") %>%
  assign_colors_per_date(date_col = "Capture_Date", nest_col = "Nest")


set.seed(98)
# =============================================================================
# 3. PRIOR SPECIFICATION
# =============================================================================

priors_chicks <- c(
  
  # --- Intercepts ---
  prior(normal(84.5,   3.25), class = "Intercept", resp = "observed",                lb = 0),
  prior(normal(10,   0.30), class = "Intercept", resp = "FaithPD",                 lb = 0),
  prior(normal(12.5,   1.30), class = "Intercept", resp = "diversityinversesimpson", lb = 0),
  prior(normal(3,    0.15), class = "Intercept", resp = "diversityshannon",        lb = 0),
  
  # --- Slopes: observed ---
  prior(normal( 0.05, 0.50), class = "b", coef = "Relative_to_frass",    resp = "observed"),
  prior(normal( 0.02, 0.10), class = "b", coef = "IRelative_to_frassE2", resp = "observed"),
  
  # --- Slopes: inverse Simpson ---
  prior(normal(-0.55, 0.25), class = "b", coef = "Relative_to_frass",    resp = "diversityinversesimpson"),
  prior(normal( 0.04, 0.05), class = "b", coef = "IRelative_to_frassE2", resp = "diversityinversesimpson"),
  
  # --- Slopes: Shannon ---
  prior(normal(-0.04, 0.05),  class = "b", coef = "Relative_to_frass",    resp = "diversityshannon"),
  prior(normal( 0.00, 0.01), class = "b", coef = "IRelative_to_frassE2", resp = "diversityshannon"),
  
  # --- Slopes: Faith PD ---
  prior(normal(0.01, 0.05), class = "b", coef = "Relative_to_frass",    resp = "FaithPD"),
  prior(normal(0.00, 0.05), class = "b", coef = "IRelative_to_frassE2", resp = "FaithPD"),
  
  # --- Sigmas ---
  prior(normal(21.7,  1.25), class = "sigma", resp = "observed",                lb = 0),
  prior(normal(6.30,  0.50), class = "sigma", resp = "diversityinversesimpson",  lb = 0),
  prior(normal(0.5,   0.05), class = "sigma", resp = "diversityshannon",         lb = 0),
  prior(normal(2.00,  0.15), class = "sigma", resp = "FaithPD",                  lb = 0)
)


# =============================================================================
# 4. MODEL STRUCTURE
# =============================================================================

# Define multivariate model formula with correlated residuals
model_chicks <-
  bf(observed                  ~ I(Relative_to_frass^2) + Relative_to_frass + (1 | p | Nest),
     family = gaussian()) +
  bf(diversity_inverse_simpson ~ I(Relative_to_frass^2) + Relative_to_frass + (1 | p | Nest),
     family = gaussian()) +
  bf(diversity_shannon         ~ I(Relative_to_frass^2) + Relative_to_frass + (1 | p | Nest),
     family = gaussian()) +
  bf(Faith_PD                  ~ I(Relative_to_frass^2) + Relative_to_frass + (1 | p | Nest),
     family = gaussian()) +
  set_rescor(TRUE)  # Model residuals as correlated


# =============================================================================
# 5. MODEL FITTING
# =============================================================================

multivariate_model_chicks <- brm(
  model_chicks,
  data    = alpha_diversity_chicks,
  iter    = 20000,   # Total iterations per chain
  warmup  = 5000,   # Burn-in iterations
  chains  = 4,
  cores   = 4,
  prior   = priors_chicks,
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)


# =============================================================================
# 6. MODEL SUMMARY & SAVING
# =============================================================================

summary(multivariate_model_chicks)

saveRDS(multivariate_model_chicks, file = "../rds_objects/multivariate_model_chicks.rds")


# =============================================================================
# 7. MODEL DIAGNOSTICS
# =============================================================================

# Load saved model
multivariate_model_chicks <- readRDS("../rds_objects/multivariate_model_chicks.rds")

# Effective sample size plot
mcmc_plot(multivariate_model_chicks, type = "neff")

# Trace plots for key slope parameters
mcmc_trace(
  multivariate_model_chicks,
  pars = c(
    "b_observed_Relative_to_frass",
    "b_diversityinversesimpson_Relative_to_frass",
    "b_diversityshannon_Relative_to_frass",
    "b_FaithPD_Relative_to_frass"
  )
)

# HMC diagnostics
check_hmc_diagnostics(multivariate_model_chicks$fit)

# Posterior predictive checks per response variable
pp_check(multivariate_model_chicks, resp = "observed",                 ndraws = 200)
pp_check(multivariate_model_chicks, resp = "diversityinversesimpson",  ndraws = 200)
pp_check(multivariate_model_chicks, resp = "diversityshannon",         ndraws = 200)
pp_check(multivariate_model_chicks, resp = "FaithPD",                  ndraws = 200)

# Optional: check mean and SD recovery
# pp_check(multivariate_model_chicks, resp = "observed",
#          type = "stat_2d", stat = c("mean", "sd"))


# =============================================================================
# 8. VARIANCE DECOMPOSITION
# =============================================================================

vc <- VarCorr(multivariate_model_chicks, summary = TRUE) # Extract variance-covariance components

# Random-effect standard deviations (Nest)
random_sd <- vc$Nest$sd[, "Estimate"]
random_response <- sub("_Intercept$", "", rownames(vc$Nest$sd))

# Residual standard deviations
residual_sd <- vc$residual$sd[, "Estimate"]
residual_response <- sub("_Intercept$", "", rownames(vc$residual$sd))

var_random <- random_sd^2 # Convert standard deviations to variances
var_residual <- residual_sd^2

var_table <- data.frame(
  Response = random_response,
  Var_Random = unname(var_random)
)

var_table$Var_Residual <- var_residual[
  match(var_table$Response, residual_response) # Match residual variances to the corresponding responses
]

# Calculate total variance and variance proportions
var_table$Var_Total <- with(
  var_table,
  Var_Random + Var_Residual
)

var_table$Prop_Random <- with(
  var_table,
  Var_Random / Var_Total
)

var_table$Prop_Residual <- with(
  var_table,
  Var_Residual / Var_Total
)

# Display the results
print(var_table, digits = 3)

# =============================================================================
# 9. VISUALISATION
# =============================================================================

# Shared axis and theme settings
shared_y_expansion <- expansion(mult = c(0.05, 0.15))  # 5% padding below, 15% above
shared_point_size  <- 7
shared_alpha       <- 0.6
shared_base_size   <- 24

# --- Observed richness ---
plot_chick_observed <- ggplot(alpha_diversity_chicks, aes(x = Relative_to_frass, y = observed)) +
  geom_point(aes(color = ColorHex), alpha = shared_alpha, size = shared_point_size, show.legend = FALSE) +
  scale_color_identity() +
  labs(y = "Observed richness") +
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(
    limits = c(0, NA),
    n.breaks = 8,
    expand = shared_y_expansion
  ) +
  theme_minimal(base_size = shared_base_size) +
  theme(
    panel.grid   = element_blank(),
    axis.title.x = element_blank()
  )

# --- Faith's PD ---
plot_chick_faithPD <- ggplot(alpha_diversity_chicks, aes(x = Relative_to_frass, y = Faith_PD)) +
  geom_point(aes(color = ColorHex), alpha = shared_alpha, size = shared_point_size, show.legend = FALSE) +
  scale_color_identity() +
  labs(y = "Faith's PD") +
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(
    limits = c(0, NA),
    n.breaks = 8,
    expand = shared_y_expansion
  ) +
  theme_minimal(base_size = shared_base_size) +
  theme(
    panel.grid   = element_blank(),
    axis.title.x = element_blank()
  )

# --- Inverse Simpson ---
plot_chick_simpson <- ggplot(alpha_diversity_chicks, aes(x = Relative_to_frass, y = diversity_inverse_simpson)) +
  geom_point(aes(color = ColorHex), alpha = shared_alpha, size = shared_point_size, show.legend = FALSE) +
  scale_color_identity() +
  labs(y = "Inverse Simpson") +
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(
    limits = c(0, NA),
    n.breaks = 6,
    expand = shared_y_expansion
  ) +
  geom_abline(intercept = 12.08, slope = -0.56, color = "black", linetype = "dashed") +
  theme_minimal(base_size = shared_base_size) +
  theme(
    panel.grid   = element_blank(),
    axis.title.x = element_blank()
  )

# --- Shannon diversity ---
plot_chick_shannon <- ggplot(alpha_diversity_chicks, aes(x = Relative_to_frass, y = diversity_shannon)) +
  geom_point(aes(color = ColorHex), alpha = shared_alpha, size = shared_point_size, show.legend = FALSE) +
  scale_color_identity() +
  labs(x = "Days to frass peak", y = "Shannon") +
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(
    limits = c(0, NA),
    n.breaks = 9,
    expand = shared_y_expansion
  ) +
  geom_abline(intercept = 2.94, slope = -0.04, color = "black", linetype = "dashed") +
  theme_minimal(base_size = shared_base_size) +
  theme(panel.grid = element_blank())


# =============================================================================
# 10. PANEL ASSEMBLY & EXPORT
# =============================================================================

chicks_panel <- cowplot::plot_grid(
  plot_chick_observed,
  plot_chick_faithPD,
  plot_chick_simpson,
  plot_chick_shannon,
  ncol = 1
)

ggsave(
  filename = "../results/multivariate_chicks_panel.svg",
  plot     = chicks_panel,
  width    = 297,
  height   = 400,
  units    = "mm"
)
