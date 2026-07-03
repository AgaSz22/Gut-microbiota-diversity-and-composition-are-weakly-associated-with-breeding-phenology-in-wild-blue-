# =============================================================================
# Purpose: Provide descriptive statistics about reads per sample and conduct iterative rarefaction 
# Author:  Agnes Sz.
# Date:    2025-07-09
# =============================================================================

# ── 1. Libraries ──────────────────────────────────────────────────────────────

# --- Phylogenetic & Microbiome Analysis ---
library(phyloseq)
library(microbiome)

# --- Data Wrangling ---
library(dplyr)

# --- Visualisation ---
library(ggplot2)
library(patchwork)
library(cowplot)

# --- Statistical Modelling ---
library(lme4)

# ── 2. Load Input Data ────────────────────────────────────────────────────────

phylo_clean    <- readRDS("../rds_objects/phyloseq_cleaned.rds")
phylo_with_cont <- readRDS("../rds_objects/phylo_with_cont.rds") 

# ── 3. Summarise read counts ────────────────────────────────────────────────────────

ps_no_negatives <- subset_samples(phylo_clean, Stage != "Negative")
ps_negatives <- subset_samples(phylo_clean, Stage == "Negative")
ps_cont_negatives <- subset_samples(phylo_with_cont, Stage == "Negative")

#' Summarise read counts for a phyloseq object
#' @param ps A phyloseq object
#' @param label A descriptive label for printing
summarise_readcounts <- function(ps, label) {
  rc <- microbiome::readcount(ps)
  cat(sprintf(
    "\n[%s]\n  Min:  %d\n  Mean: %.2f\n  Max:  %d\n",
    label, min(rc), mean(rc), max(rc)
  ))
  invisible(rc)
}

summarise_readcounts(ps_no_negatives,   "Non-negative samples (without contaminants)")
summarise_readcounts(ps_negatives,      "Negative controls (without contaminants)")
summarise_readcounts(ps_cont_negatives, "Negative controls (with contaminants)")

# ── 4. Pick rarefaction depth ────────────────────────────────────────────────────────

read_counts <- microbiome::readcount(ps_no_negatives) |>
  as.data.frame() |>
  setNames("Read_counts") |>
  tibble::rownames_to_column("Sample")


RAREFACTION_THRESHOLD <- 9000  

set.seed(46)
read_counts |>
  ggplot(aes(x = 1, y = Read_counts)) +
  geom_jitter(width = 0.2, alpha = 0.6, colour = "steelblue") +
  geom_hline(yintercept = RAREFACTION_THRESHOLD,
             linetype = "dashed", colour = "firebrick", linewidth = 0.8) +
  scale_y_log10(labels = scales::comma) +
  labs(
    title    = "Read count distribution (non-negative samples)",
    subtitle = paste("Dashed line = rarefaction threshold:", 
                     scales::comma(RAREFACTION_THRESHOLD)),
    x        = NULL,
    y        = "Read count (log10 scale)"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


# ── 5. Rarefy and save outputs────────────────────────────────────────────────────────

RAREFACTION_SIZE  <- 9000
RAREFACTION_ITERS <- 500
MASTER_SEED       <- 2605
OUTPUT_PATH       <- "../results/alpha_df_mean.rds"

set.seed(MASTER_SEED)
rngseed_vec <- sample(x = 1:100000, size = RAREFACTION_ITERS, replace = FALSE)

# --- Iterative Rarefaction ---
message(sprintf("Starting iterative rarefaction: %d iterations at depth %d",
                RAREFACTION_ITERS, RAREFACTION_SIZE))

# Initialise accumulator on iteration 1
alpha_df_sum <- phyloseq::rarefy_even_depth(
  phylo_clean,
  sample.size = RAREFACTION_SIZE,
  rngseed     = rngseed_vec[1],
  verbose     = FALSE
) |>
  microbiome::alpha(index = "all")

# Accumulate iterations 2:N
for (i in seq(2, RAREFACTION_ITERS)) {
  
  alpha_iter <- phyloseq::rarefy_even_depth(
    phylo_clean,
    sample.size = RAREFACTION_SIZE,
    rngseed     = rngseed_vec[i],
    verbose     = FALSE
  ) |>
    microbiome::alpha(index = "all")
  
  alpha_df_sum <- alpha_df_sum + alpha_iter
  
  if (i %% 50 == 0) message(sprintf("  Completed iteration %d / %d", i, RAREFACTION_ITERS))
}

# --- Average across iterations ---
alpha_df_mean <- alpha_df_sum / RAREFACTION_ITERS

# --- Save  ---
saveRDS(alpha_df_mean, file = OUTPUT_PATH)
