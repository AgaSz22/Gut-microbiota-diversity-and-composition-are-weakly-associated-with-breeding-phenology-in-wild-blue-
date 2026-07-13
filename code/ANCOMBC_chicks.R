# =============================================================================
# ANCOMBC2 Differential Abundance Analysis
# =============================================================================

library(ANCOMBC)
library(tidyverse)
library(DT)
library(phyloseq)

# =============================================================================
# 1. Load and Prepare Data
# =============================================================================

phylo_clean <- readRDS("../rds_objects/phyloseq_cleaned.rds")

# Extract metadata and compute Days_to_frass vectorised (no loop needed)
meta_data <- microbiome::meta(phylo_clean)
meta_data$Days_to_frass <- meta_data$Capture_Date - 4 - 43

# Assign updated metadata back to phyloseq object
sample_data(phylo_clean) <- sample_data(meta_data)

# Subset by life stage
# phylo_filtered_females <- subset_samples(phylo_clean, Stage == "Female")
phylo_filtered_chicks  <- subset_samples(phylo_clean, Stage == "Chick")
# phylo_filtered_both    <- subset_samples(phylo_clean, Stage %in% c("Chick", "Female"))

# =============================================================================
# 2. ANCOMBC2 Helper Function
# =============================================================================

#' Run ANCOMBC2 model, extract results, and save output
#'
#' @param data      A phyloseq object.
#' @param tax_level Taxonomic level string (e.g. "Species", "Genus").
#' @param lib_cut   Minimum library size cut-off (default 1000).
#' @param save_path File path for saving the RDS output.
#'
#' @return The full ancombc2 output object.
run_ancombc2 <- function(data, tax_level, lib_cut = 1000, save_path) {
  set.seed(123)
  
  output <- ancombc2(
    data          = data,
    tax_level     = tax_level,
    fix_formula   = "Days_to_frass",
    rand_formula  = "(1 | Nest)",
    p_adj_method  = "holm",
    pseudo_sens   = TRUE,
    prv_cut       = 0.05,
    lib_cut       = lib_cut,
    s0_perc       = 0.05,
    struc_zero    = FALSE,
    neg_lb        = TRUE,
    alpha         = 0.05,
    n_cl          = 4,
    verbose       = TRUE,
    global        = FALSE,
    pairwise      = FALSE,
    dunnet        = FALSE,
    trend         = FALSE,
    iter_control  = list(tol = 1e-2, max_iter = 50, verbose = TRUE),
    em_control    = list(tol = 1e-5, max_iter = 100),
    lme_control   = lme4::lmerControl(),
    mdfdr_control = list(fwer_ctrl_method = "holm", B = 100)
  )
  
  saveRDS(output, save_path)
  return(output)
}

# =============================================================================
# 3. Run ANCOMBC2 Models (Chicks) Across Taxonomic Levels
# =============================================================================

# Species
output_species_chicks  <- run_ancombc2(
  data      = phylo_filtered_chicks,
  tax_level = "Species",
  save_path = "../rds_objects/ANCOMBC_output_species_chicks.rds"
)
results_chicks_species <- output_species_chicks$res

# Genus
output_genus_chicks  <- run_ancombc2(
  data      = phylo_filtered_chicks,
  tax_level = "Genus",
  save_path = "../rds_objects/ANCOMBC_output_genus_chicks.rds"
)
results_chicks_genus <- output_genus_chicks$res

# Family
output_family_chicks  <- run_ancombc2(
  data      = phylo_filtered_chicks,
  tax_level = "Family",
  save_path = "../rds_objects/ANCOMBC_output_family_chicks.rds"
)
results_chicks_family <- output_family_chicks$res

# Class
output_class_chicks  <- run_ancombc2(
  data      = phylo_filtered_chicks,
  tax_level = "Class",
  save_path = "../rds_objects/ANCOMBC_output_class_chicks.rds"
)
results_chicks_class <- output_class_chicks$res

# Order — note: lib_cut differs from other levels (9000 vs 1000 and 5% vs 10% of taxa)
output_order_chicks  <- run_ancombc2(
  data      = phylo_filtered_chicks,
  tax_level = "Order",
  lib_cut   = 9000,
  save_path = "../rds_objects/ANCOMBC_output_order_chicks.rds"
)
results_chicks_order <- output_order_chicks$res

# Phylum
output_phylum_chicks  <- run_ancombc2(
  data      = phylo_filtered_chicks,
  tax_level = "Phylum",
  save_path = "../rds_objects/ANCOMBC_output_phylum_chicks.rds"
)
results_chicks_phylum <- output_phylum_chicks$res

# =============================================================================
# 4. Plot: Log Fold Changes at Order Level (Chicks ~ Days to Frass)
# =============================================================================

# --- Prepare data ---
df_DTF <- results_chicks_order %>%
  dplyr::select(taxon, ends_with("Days_to_frass"))

df_fig_DTF <- df_DTF %>%
  dplyr::filter(diff_Days_to_frass == 1) %>%
  dplyr::arrange(desc(lfc_Days_to_frass)) %>%
  dplyr::mutate(
    direct = ifelse(lfc_Days_to_frass > 0, "Positive LFC", "Negative LFC"),
    # Colour axis labels red if taxon did NOT pass sensitivity analysis
    color  = ifelse(passed_ss_Days_to_frass == 1, "black", "red")
  )

# Fix factor order so bars plot in descending LFC order
df_fig_DTF$taxon  <- factor(df_fig_DTF$taxon,  levels = df_fig_DTF$taxon)
df_fig_DTF$direct <- factor(df_fig_DTF$direct, levels = c("Positive LFC", "Negative LFC"))

# --- Build plot ---
fig_DTF <- df_fig_DTF %>%
  ggplot(aes(x = taxon, y = lfc_Days_to_frass, fill = direct)) +
  geom_bar(
    stat  = "identity",
    width = 0.7,
    color = "black"
    # position_dodge removed: each taxon has exactly one direction,
    # so dodging is unnecessary and misleading
  ) +
  geom_errorbar(
    aes(
      ymin = lfc_Days_to_frass - se_Days_to_frass,
      ymax = lfc_Days_to_frass + se_Days_to_frass
    ),
    width = 0.2,
    color = "black"
    # position_dodge removed for the same reason as above
  ) +
  labs(
    x     = NULL,
    y     = "Log fold change",
    title = "Log fold changes per unit increase in Days to Frass (Order level, Chicks)"
  ) +
  scale_fill_discrete(name = NULL) +
  theme_bw() +
  theme(
    plot.title       = element_text(hjust = 0.5),
    panel.grid.minor.y = element_blank(),
    # Red labels = failed sensitivity analysis; black = passed
    axis.text.x      = element_text(
      angle = 60,
      hjust = 1,
      color = df_fig_DTF$color
    )
  )

fig_DTF
