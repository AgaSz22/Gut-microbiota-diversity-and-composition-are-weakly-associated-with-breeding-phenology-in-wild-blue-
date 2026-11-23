# =============================================================================
# Date:    09/06/2026
# Author:  Agnes Szwarczynska
# Purpose: Run a multivariate brms model on females data and perform variance partitioning
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


# =============================================================================
# 2. DATA LOADING & PREPARATION
# =============================================================================

all_diversity_indices <- read.csv("../data/all_diversity_indices.csv")

# Filter to female stage; remove 2 duplicates with lower number of reads
alpha_diversity_females <- all_diversity_indices %>%
  filter(Stage == "Female" & !Sample_ID %in% c("S71", "S70.2"))


set.seed(98)
# =============================================================================
# 3. PRIOR SPECIFICATION
# =============================================================================

priors_females <- c(
  
  # --- Intercepts ---
  prior(normal(67, 5), class = "Intercept", resp = "observed",                lb = 0),
  prior(normal(6.2,  1.2), class = "Intercept", resp = "diversityinversesimpson", lb = 0),
  prior(normal(2.3,  0.2), class = "Intercept", resp = "diversityshannon",        lb = 0),
  prior(normal(8.3,  0.5), class = "Intercept", resp = "FaithPD",                 lb = 0),
  
  # --- Slopes: observed ---
  prior(normal( 0.4, 0.7), class = "b", coef = "Laying_day",    resp = "observed"),
  prior(normal(-0.01, 0.05), class = "b", coef = "ILaying_dayE2", resp = "observed"),
  
  # --- Slopes: inverse Simpson ---
  prior(normal( 0.55, 0.2), class = "b", coef = "Laying_day",    resp = "diversityinversesimpson"),
  prior(normal(-0.02, 0.05), class = "b", coef = "ILaying_dayE2", resp = "diversityinversesimpson"),
  
  # --- Slopes: Shannon ---
  prior(normal( 0.02, 0.05), class = "b", coef = "Laying_day",    resp = "diversityshannon"),
  prior(normal(0.00, 0.01), class = "b", coef = "ILaying_dayE2", resp = "diversityshannon"),
  
  # --- Slopes: Faith PD ---
  prior(normal(-0.05, 0.1), class = "b", coef = "Laying_day",    resp = "FaithPD"),
  prior(normal( 0.00, 0.01), class = "b", coef = "ILaying_dayE2", resp = "FaithPD"),
  
  # --- Sigmas  ---
  prior(normal(36, 2), class = "sigma", resp = "observed",                lb = 0),
  prior(normal(7.8,  0.5), class = "sigma", resp = "diversityinversesimpson", lb = 0),
  prior(normal(1,  0.1), class = "sigma", resp = "diversityshannon",        lb = 0),
  prior(normal(3.2,  0.2), class = "sigma", resp = "FaithPD",                 lb = 0)
)

# =============================================================================
# 4. MODEL STRUCTURE
# =============================================================================

# Note: NAs present for two females missing egg-laying dates
# Define multivariate model formula with correlated residuals
model_females <-
  bf(observed                  ~ I(Laying_day^2) + Laying_day + (1 | p | Capture_Date),
     family = gaussian()) +
  bf(diversity_inverse_simpson ~ I(Laying_day^2) + Laying_day + (1 | p | Capture_Date),
     family = gaussian()) +
  bf(diversity_shannon         ~ I(Laying_day^2) + Laying_day + (1 | p | Capture_Date),
     family = gaussian()) +
  bf(Faith_PD                  ~ I(Laying_day^2) + Laying_day + (1 | p | Capture_Date),
     family = gaussian()) +
  set_rescor(TRUE)  # Model residuals as correlated


# =============================================================================
# 5. MODEL FITTING
# =============================================================================

multivariate_model_females <- brm(
  model_females,
  data    = alpha_diversity_females,
  iter    = 20000,   # Total iterations per chain
  warmup  = 5000,   # Burn-in iterations
  chains  = 4,
  cores   = 4,
  prior   = priors_females,
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)


# =============================================================================
# 6. MODEL SUMMARY & SAVING
# =============================================================================

summary(multivariate_model_females)

saveRDS(multivariate_model_females, file = "../rds_objects/multivariate_model_females.rds")


# =============================================================================
# 7. MODEL DIAGNOSTICS
# =============================================================================

# Load saved model
multivariate_model_females <- readRDS("../rds_objects/multivariate_model_females.rds")

# Effective sample size plot
mcmc_plot(multivariate_model_females, type = "neff")

# Trace plots for key slope parameters
mcmc_trace(
  multivariate_model_females,
  pars = c(
    "b_observed_Laying_day",
    "b_diversityinversesimpson_Laying_day",
    "b_diversityshannon_Laying_day",
    "b_FaithPD_Laying_day"
  )
)

# HMC diagnostics
check_hmc_diagnostics(multivariate_model_females$fit)

# Posterior predictive checks per response variable
pp_check(multivariate_model_females, resp = "observed",                ndraws = 200)
pp_check(multivariate_model_females, resp = "diversityinversesimpson", ndraws = 200)
pp_check(multivariate_model_females, resp = "diversityshannon",        ndraws = 200)
pp_check(multivariate_model_females, resp = "FaithPD",                 ndraws = 200)

# Optional: check mean and SD recovery
# pp_check(multivariate_model_females, resp = "observed",
#          type = "stat_2d", stat = c("mean", "sd"))


# =============================================================================
# 8. VARIANCE DECOMPOSITION
# =============================================================================

vc <- VarCorr(multivariate_model_females, summary = TRUE) # Extract variance-covariance components

# Random-effect standard deviations (Nest)
random_sd <- vc$Capture_Date$sd[, "Estimate"]
random_response <- sub("_Intercept$", "", rownames(vc$Capture_Date$sd))

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
plot_females_observed <- ggplot(alpha_diversity_females, aes(x = Laying_day, y = observed)) +
  geom_point(alpha = shared_alpha, size = shared_point_size, colour = "#F56A6A") +
  labs(x = "Egg-laying date", y = "Observed richness") +
  scale_x_continuous(n.breaks = 12) +
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
plot_females_faith <- ggplot(alpha_diversity_females, aes(x = Laying_day, y = Faith_PD)) +
  geom_point(alpha = shared_alpha, size = shared_point_size, colour = "#A3334E") +
  labs(x = "Egg-laying date", y = "Faith's PD") +
  scale_x_continuous(n.breaks = 12) +
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
plot_females_simpson <- ggplot(alpha_diversity_females, aes(x = Laying_day, y = diversity_inverse_simpson)) +
  geom_point(alpha = shared_alpha, size = shared_point_size, colour = "#B39EB5") +
  labs(x = "Egg-laying date", y = "Inverse Simpson") +
  geom_abline(intercept = 5.65, slope = 0.54, color = "black", linetype = "dashed") +
  scale_x_continuous(n.breaks = 12) +
  scale_y_continuous(
    limits = c(0, NA),
    n.breaks = 6,
    expand = shared_y_expansion
  ) +
  theme_minimal(base_size = shared_base_size) +
  theme(
    panel.grid   = element_blank(),
    axis.title.x = element_blank()
  )

# --- Shannon diversity ---
plot_females_shannon <- ggplot(alpha_diversity_females, aes(x = Laying_day, y = diversity_shannon)) +
  geom_point(alpha = shared_alpha, size = shared_point_size, colour = "#330C15") +
  labs(x = "Egg-laying date", y = "Shannon") +
  scale_x_continuous(n.breaks = 12) +
  scale_y_continuous(
    limits = c(0, NA),
    n.breaks = 9,
    expand = shared_y_expansion
  ) +
  theme_minimal(base_size = shared_base_size) +
  theme(panel.grid = element_blank())


# =============================================================================
# 10. PANEL ASSEMBLY & EXPORT
# =============================================================================

females_panel <- cowplot::plot_grid(
  plot_females_observed,
  plot_females_faith,
  plot_females_simpson,
  plot_females_shannon,
  ncol = 1  # Arrange plots vertically
)

ggsave(
  filename = "../results/multivariate_females_panel.svg",
  plot     = females_panel,
  width    = 297,
  height   = 400,
  units    = "mm"
)
