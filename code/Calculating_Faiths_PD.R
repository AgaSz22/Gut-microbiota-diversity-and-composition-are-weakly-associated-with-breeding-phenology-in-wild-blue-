# ==============================================================================
# Faith's Phylogenetic Diversity (PD) Calculation — Microbial Amplicon Data
# ==============================================================================
# Dependencies: BiocManager, ape, phangorn, phyloseq, picante, tidyverse, DECIPHER
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Install & Load Packages
# ------------------------------------------------------------------------------

# -- 1a. Installation (commented out — uncomment to run on a fresh environment) -

# # Install Bioconductor manager if not present
# if (!require("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
#
# # Install Bioconductor packages (skipped if already installed)
# BiocManager::install(
#   c("DECIPHER", "Biostrings", "phyloseq"),
#   update = FALSE,
#   ask    = FALSE
# )
#
# # Install CRAN packages if missing
# cran_packages <- c("ape", "phangorn", "picante", "lme4", "tidyverse")
# installed     <- rownames(installed.packages())
# to_install    <- cran_packages[!cran_packages %in% installed]
# if (length(to_install) > 0) install.packages(to_install)

# -- 1b. Load Libraries --------------------------------------------------------

library(ape)
library(phangorn)
library(phyloseq)
library(picante)
library(lme4)
library(tidyverse)
library(DECIPHER)
library(Biostrings)
library(parallel)


# ------------------------------------------------------------------------------
# 2. Load Data
# ------------------------------------------------------------------------------

# Load the cleaned phyloseq object produced by the upstream preprocessing script
phylo_clean <- readRDS("../rds_objects/phyloseq_cleaned.rds")


# ------------------------------------------------------------------------------
# 3. Build Phylogenetic Tree from Reference Sequences
# ------------------------------------------------------------------------------

# -- 3a. Extract & Align Reference Sequences -----------------------------------

# Pull reference sequences from the phyloseq object and label with taxon names
seqs        <- Biostrings::DNAStringSet(refseq(phylo_clean))
names(seqs) <- taxa_names(phylo_clean)

# Perform multiple sequence alignment using DECIPHER
# processors = NULL uses all available cores automatically
alignment <- DECIPHER::AlignSeqs(seqs, processors = NULL)

# -- 3b. Compute Distance Matrix -----------------------------------------------

# Convert the alignment to phangorn's phyDat format, then compute ML distances
phang_align <- phangorn::phyDat(as(alignment, "matrix"), type = "DNA")
dist_mat    <- phangorn::dist.ml(phang_align)

# -- 3c. Build Initial Neighbour-Joining Tree & Fit PML Model ------------------

# Construct a starting NJ tree from the distance matrix
tree_nj <- phangorn::NJ(dist_mat)

# Fit an initial PML model using the NJ tree as the starting topology
fit_pml <- phangorn::pml(tree_nj, data = phang_align)

# -- 3d. Optimise Tree Under GTR Model -----------------------------------------

# NOTE: GTR optimisation can be slow for large datasets.
#       To skip and use the unoptimised PML tree instead, comment out this block
#       and replace `fit_gtr` with `fit_pml` in all downstream steps.
fit_gtr <- phangorn::optim.pml(
  fit_pml,                    # Use the PML fit as the starting point
  model         = "GTR",      # General Time Reversible substitution model
  optInv        = TRUE,       # Optimise proportion of invariant sites
  optGamma      = TRUE,       # Optimise gamma rate heterogeneity
  rearrangement = "stochastic" # Stochastic nearest-neighbour interchange
)

# -- 3e. Root the Tree & Validate ----------------------------------------------

# Root the optimised tree using the first taxon as the outgroup
phy_tree_rooted <- ape::root(
  fit_gtr$tree,
  outgroup     = taxa_names(phylo_clean)[1],
  resolve.root = TRUE
)

# Confirm the tree is rooted before proceeding — Faith's PD requires a rooted tree
stopifnot(
  "Tree must be rooted before Faith's PD calculation" = is.rooted(phy_tree_rooted)
)

# -- 3f. Save Tree Objects -----------------------------------------------------

# Save GTR-optimised tree in both .Rdata and .rds formats
save(fit_gtr,         file = "../data/PD_phylotree.Rdata")
save(phy_tree_rooted, file = "../data/phy_tree_rooted.Rdata")

saveRDS(fit_gtr,         file = "../data/PD_phylotree.rds")
saveRDS(phy_tree_rooted, file = "../data/phy_tree_rooted.rds")

# Save the unoptimised PML tree for reference / comparison
save(fit_pml, file = "../data/PD_phylotree_pml.Rdata")


# ------------------------------------------------------------------------------
# 4. Rarefaction Setup
# ------------------------------------------------------------------------------

# Define rarefaction depth and number of iterations
rarefaction_size  <- 9000   # Minimum non-negative library size across samples
rarefaction_iters <- 500    # Number of rarefaction replicates

# Generate a reproducible vector of unique random seeds — one per iteration
set.seed(2605)
rngseed_vec <- sample(x = 1:100000, size = rarefaction_iters, replace = FALSE)


# ------------------------------------------------------------------------------
# 5. Rarefaction Loop — Compute Faith's PD per Iteration
# ------------------------------------------------------------------------------

# -- 5a. Parallel Setup --------------------------------------------------------

# Detect available cores; leave 2 free for the system and cap at 50
n_cores <- min(parallel::detectCores() - 2, 50)

# -- 5b. Run Rarefaction in Parallel -------------------------------------------

# For each seed: rarefy the phyloseq object, extract the OTU matrix,
# compute Faith's PD via picante, and return a named numeric vector of PD values
faith_pd_list <- mclapply(rngseed_vec, function(seed) {
  
  # Rarefy to the defined depth using the current iteration's seed
  pseq_rare <- phyloseq::rarefy_even_depth(
    phylo_clean,
    sample.size = rarefaction_size,
    rngseed     = seed,
    verbose     = FALSE
  )
  
  # Extract OTU table as a plain matrix; ensure samples are rows
  otu_mat <- as(otu_table(pseq_rare), "matrix")
  if (taxa_are_rows(pseq_rare)) {
    otu_mat <- t(otu_mat)
  }
  
  # Compute Faith's PD for all samples in this rarefied replicate
  pd_vals <- picante::pd(otu_mat, phy_tree_rooted, include.root = TRUE)$PD
  
  # Attach sample names so results can be aligned across iterations
  names(pd_vals) <- sample_names(pseq_rare)
  
  return(pd_vals)
  
}, mc.cores = n_cores)

# Save the raw per-iteration PD list before aggregation
save(faith_pd_list, file = "../data/faith_pd_list_initial.Rdata")

# -- 5c. Commented-Out Block: Reconstruct a Representative Rarefied Object -----

# # Reconstruct a representative rarefied phyloseq object using the first seed.
# # The tree and metadata are consistent across iterations; only OTU counts differ.
# pseq_rare <- phyloseq::rarefy_even_depth(
#   phylo_clean,
#   sample.size = rarefaction_size,
#   rngseed     = rngseed_vec[1],   # Use the first seed for reproducibility
#   verbose     = FALSE
# )
#
# phy_tree(pseq_rare) <- phy_tree_rooted


# ------------------------------------------------------------------------------
# 6. Aggregate Faith's PD Across Iterations
# ------------------------------------------------------------------------------

# -- 6a. Identify Stable Samples -----------------------------------------------

# Retain only samples that were present in every rarefaction iteration
# (some samples may be dropped if their library size falls below rarefaction_size)
rarefied_samples <- Reduce(intersect, lapply(faith_pd_list, names))

# -- 6b. Build Matrix & Compute Per-Sample Means -------------------------------

# Stack per-iteration PD vectors into a matrix (iterations × samples),
# then compute the column-wise mean to get a stable PD estimate per sample
faith_pd_matrix <- do.call(rbind, lapply(faith_pd_list, function(x) x[rarefied_samples]))
mean_faith_pd   <- colMeans(faith_pd_matrix)

# -- 6c. Assemble Final Data Frame ---------------------------------------------

# Combine sample IDs and their mean Faith's PD into a tidy data frame
faith_pd_df <- data.frame(
  SampleID = rarefied_samples,
  Faith_PD = mean_faith_pd[rarefied_samples]
)

# Save the final aggregated Faith's PD data frame for downstream analysis
saveRDS(faith_pd_df, "../rds_objects/faith_pd_df.rds")