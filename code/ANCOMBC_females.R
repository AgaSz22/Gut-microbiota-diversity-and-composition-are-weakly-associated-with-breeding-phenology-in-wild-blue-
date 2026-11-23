# =============================================================================
# ANCOMBC2 Differential Abundance Analysis — Females ~ Laying Day
# =============================================================================

library(ANCOMBC)
library(tidyverse)
library(DT)
library(phyloseq)

# =============================================================================
# 1. Load and Prepare Data
# =============================================================================

phylo_clean <- readRDS("../rds_objects/phyloseq_cleaned_no_plasmids.rds")

# Add Sample_ID as explicit column (used for filtering)
sample_data(phylo_clean)$Sample_ID <- rownames(sample_data(phylo_clean))

# Subset females, removing problematic samples and rows with missing Laying_day
phylo_filtered_females <- subset_samples(
  phylo_clean,
  Stage == "Female" &
    !Sample_ID %in% c("S71", "S70.2") &
    !is.na(Laying_day)
)

# =============================================================================
# 2. ANCOMBC2 Helper Function
# =============================================================================

#' Run ANCOMBC2 model, extract results, and save output
#'
#' @param data      A phyloseq object.
#' @param tax_level Taxonomic level string (e.g. "Species", "Genus").
#' @param prv_cut   A proportion of samples containing given taxon.
#' @param lib_cut   Minimum library size cut-off (default 1000).
#' @param save_path File path for saving the RDS output.
#'
#' @return The full ancombc2 output object.
run_ancombc2 <- function(data, tax_level, prv_cut = 0.1, lib_cut = 4000, save_path) {
  set.seed(123)
  
  output <- ancombc2(
    data          = data,
    tax_level     = tax_level,
    fix_formula   = "Laying_day",
    rand_formula  = "(1 | Capture_Date)",
    p_adj_method  = "holm",
    pseudo_sens   = TRUE,
    prv_cut       = prv_cut,
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
# 3. Run ANCOMBC2 Models (Females) Across Taxonomic Levels
# =============================================================================

# Phylum
output_phylum_females  <- run_ancombc2(
  data      = phylo_filtered_females,
  tax_level = "Phylum",
  save_path = "../rds_objects/ANCOMBC_output_phylum_females.rds"
)
results_phylum_females <- output_phylum_females$res

# Class
output_class_females  <- run_ancombc2(
  data      = phylo_filtered_females,
  tax_level = "Class",
  save_path = "../rds_objects/ANCOMBC_output_class_females.rds"
)
results_class_females <- output_class_females$res

# Order
output_order_females  <- run_ancombc2(
  data      = phylo_filtered_females,
  tax_level = "Order",
  save_path = "../rds_objects/ANCOMBC_output_order_females.rds"
)
results_order_females <- output_order_females$res

# Family
output_family_females  <- run_ancombc2(
  data      = phylo_filtered_females,
  tax_level = "Family",
  save_path = "../rds_objects/ANCOMBC_output_family_females.rds"
)
results_family_females <- output_family_females$res

# Genus
output_genus_females  <- run_ancombc2(
  data      = phylo_filtered_females,
  tax_level = "Genus",
  save_path = "../rds_objects/ANCOMBC_output_genus_females.rds"
)
results_genus_females <- output_genus_females$res

# Species
output_species_females  <- run_ancombc2(
  data      = phylo_filtered_females,
  tax_level = "Species",
  save_path = "../rds_objects/ANCOMBC_output_species_females.rds"
)
results_species_females <- output_species_females$res

# =============================================================================
# 4. Plot: Log Fold Changes at Phylum Level (Females ~ Laying Day)
# =============================================================================

# --- Prepare data ---
df_LD <- results_phylum_females %>%
  dplyr::select(taxon, ends_with("Laying_day"))

df_fig_LD <- df_LD %>%
  dplyr::filter(diff_Laying_day == 1) %>%
  dplyr::arrange(desc(lfc_Laying_day)) %>%
  dplyr::mutate(
    direct = ifelse(lfc_Laying_day > 0, "Positive LFC", "Negative LFC"),
    # Colour axis labels red if taxon did NOT pass sensitivity analysis
    color  = ifelse(passed_ss_Laying_day == 1, "black", "red")
  )

# Fix factor order so bars plot in descending LFC order
df_fig_LD$taxon  <- factor(df_fig_LD$taxon,  levels = df_fig_LD$taxon)
df_fig_LD$direct <- factor(df_fig_LD$direct, levels = c("Positive LFC", "Negative LFC"))

# --- Build plot ---
fig_LD <- df_fig_LD %>%
  ggplot(aes(x = taxon, y = lfc_Laying_day, fill = direct)) +
  geom_bar(
    stat  = "identity",
    width = 0.7,
    color = "black"
    # position_dodge removed: each taxon has exactly one direction,
    # so dodging is unnecessary and misleading
  ) +
  geom_errorbar(
    aes(
      ymin = lfc_Laying_day - se_Laying_day,
      ymax = lfc_Laying_day + se_Laying_day
    ),
    width = 0.2,
    color = "black"
    # position_dodge removed for the same reason as above
  ) +
  labs(
    x     = NULL,
    y     = "Log fold change",
    title = "Log fold changes per unit increase in Laying Day (Phylum level, Females)"
  ) +
  scale_fill_discrete(name = NULL) +
  theme_bw() +
  theme(
    plot.title         = element_text(hjust = 0.5),
    panel.grid.minor.y = element_blank(),
    # Red labels = failed sensitivity analysis; black = passed
    axis.text.x        = element_text(
      angle = 60,
      hjust = 1,
      color = df_fig_LD$color
    )
  )

fig_LD