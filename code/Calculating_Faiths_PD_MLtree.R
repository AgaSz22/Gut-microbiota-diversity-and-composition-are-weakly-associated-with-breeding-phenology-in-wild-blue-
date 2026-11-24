# Load all libraries
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

phylo_clean <- readRDS("../rds_objects/phyloseq_cleaned_no_plasmids.rds")

# ------------------------------------------------------------------------------
# 3. Build Phylogenetic Tree from Reference Sequences (ML Method via phangorn)
# ------------------------------------------------------------------------------

# Extract reference sequences from phyloseq object
# NOTE: Biostrings is required here — refseq() returns a DNAStringSet,
# and there is no phangorn equivalent for accessing phyloseq reference sequences
seqs        <- Biostrings::DNAStringSet(refseq(phylo_clean))
names(seqs) <- taxa_names(phylo_clean)

# Align sequences using DECIPHER
# NOTE: phangorn has no built-in alignment function, so DECIPHER is required here
alignment <- DECIPHER::AlignSeqs(seqs, processors = NULL)

# Convert aligned DNAStringSet -> DNAbin (ape) -> phyDat (phangorn)
# phangorn::as.phyDat() can coerce directly from DNAbin,
# avoiding the intermediate matrix conversion
phang_align <- phangorn::as.phyDat(
  ape::as.DNAbin(alignment),
  type = "DNA"
)

# --- Model Testing ---
# Use modelTest() to compare evolutionary models and select the best-fitting one
modeltest_results <- phangorn::modelTest(phang_align, model = "all")

# --- ML Tree Inference ---
# pml_bb() infers the ML tree using the best model identified by modelTest()
fit_ml <- phangorn::pml_bb(modeltest_results)

# --- Bootstrap Analysis ---
# Perform bootstrap analysis to assess node support
bstrees <- phangorn::bootstrap.pml(fit_ml)

# Separate the tree topology from the pml output object
ml_tree <- fit_ml$tree

# --- Root the ML Tree ---
# phangorn::midpoint() roots at the midpoint of the longest tip-to-tip path
# No outgroup required — fully native to phangorn
phy_tree_rooted <- phangorn::midpoint(ml_tree)

# Validate that the tree is rooted before proceeding
stopifnot("Tree must be rooted before Faith's PD calculation" = ape::is.rooted(phy_tree_rooted))

# --- Visualisation ---
# Plot the rooted ML tree with bootstrap support values overlaid
phangorn::plotBS(phy_tree_rooted, bstrees, main = "ML Tree with Bootstrap Support")

# Save outputs
save(fit_ml,          file = "../data/PD_phylotree_ML.Rdata")
save(phy_tree_rooted, file = "../data/phy_tree_rooted_ML.Rdata")

saveRDS(fit_ml,          file = "../data/PD_phylotree_ML.rds")
saveRDS(phy_tree_rooted, file = "../data/phy_tree_rooted_ML.rds")

# ------------------------------------------------------------------------------
# 4. Rarefaction Setup
# ------------------------------------------------------------------------------

rarefaction_size  <- 7000
rarefaction_iters <- 500

set.seed(2605)
rngseed_vec <- sample(x = 1:100000, size = rarefaction_iters, replace = FALSE)

# ------------------------------------------------------------------------------
# 5. Rarefaction Loop — Compute Faith's PD per Iteration
# ------------------------------------------------------------------------------

n_cores <- min(parallel::detectCores() - 2, 50)

faith_pd_list <- mclapply(rngseed_vec, function(seed) {
  
  pseq_rare <- phyloseq::rarefy_even_depth(
    phylo_clean,
    sample.size = rarefaction_size,
    rngseed     = seed,
    verbose     = FALSE
  )
  
  otu_mat <- as(otu_table(pseq_rare), "matrix")
  if (taxa_are_rows(pseq_rare)) {
    otu_mat <- t(otu_mat)
  }
  
  pd_vals <- picante::pd(otu_mat, phy_tree_rooted, include.root = TRUE)$PD
  
  names(pd_vals) <- sample_names(pseq_rare)
  
  return(pd_vals)
  
}, mc.cores = n_cores)

save(faith_pd_list, file = "../data/faith_pd_list_MLtree.Rdata")

# Identify samples that survived rarefaction in every iteration
rarefied_samples <- Reduce(intersect, lapply(faith_pd_list, names))

# Build matrix and compute means
faith_pd_matrix <- do.call(rbind, lapply(faith_pd_list, function(x) x[rarefied_samples]))
mean_faith_pd   <- colMeans(faith_pd_matrix)

# Build the dataframe
faith_pd_df <- data.frame(
  SampleID = rarefied_samples,
  Faith_PD = mean_faith_pd[rarefied_samples]
)

saveRDS(faith_pd_df, "../rds_objects/faith_pd_df_MLtree.rds")

pdf("../figures/ml_tree_bootstrap.pdf")
phangorn::plotBS(phy_tree_rooted, bstrees, main = "ML tree with bootstrap support")
dev.off()