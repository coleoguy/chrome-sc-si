library(ape)
library(diversitree)

## Read phenotype data
dat <- read.csv("../data/full_phenodat.csv", stringsAsFactors = FALSE)
state_col <- "mating_system"

tree_files <- c(
  aster = "../data/aster.tree",
  brass = "../data/brass.tree",
  fab   = "../data/fab.tree",
  sol   = "../data/sol.tree")

stats <- list()

for (i in 1:length(tree_files)) {
  clade <- names(tree_files)[i]
  cat("Processing", clade, "...\n")
  
  tr <- read.tree(tree_files[i])
  
  ## Before pruning
  match_before <- intersect(dat$genus_species, tr$tip.label)
  n_match_before <- length(match_before)
  dat_before <- dat[dat$genus_species %in% match_before, ]
  
  # intuitive counts: just count how many "si" and "sc"
  n_SI_before <- sum(dat_before[[state_col]] == "si")
  n_SC_before <- sum(dat_before[[state_col]] == "sc")
  
  ## Inline prune of tips descending from polytomies
  internal_nodes <- (Ntip(tr) + 1):(Ntip(tr) + tr$Nnode) # list of all internal nodes
  deg <- tabulate(tr$edge[, 1], nbins = max(tr$edge[, 1])) #count how many edges descend from each node
  poly_nodes <- internal_nodes[deg[internal_nodes] > 2] # find the polytomies
  
  poly_edges <- tr$edge[tr$edge[, 1] %in% poly_nodes, , drop = FALSE]
  tip_idx_to_drop <- poly_edges[poly_edges[, 2] <= Ntip(tr), 2]
  tips_to_drop <- tr$tip.label[tip_idx_to_drop]
  tr <- drop.tip(tr, tips_to_drop)
  
  ## After pruning
  match_after <- intersect(dat$genus_species, tr$tip.label)
  n_match_after <- length(match_after)
  dat_after <- dat[dat$genus_species %in% match_after, , drop = FALSE]
  
  n_SI_after <- sum(dat_after[[state_col]] == "si", na.rm = TRUE)
  n_SC_after <- sum(dat_after[[state_col]] == "sc", na.rm = TRUE)
  
  ## Save stats for this clade
  stats[[clade]] <- data.frame(
    clade           = clade,
    n_match_before  = n_match_before,
    n_match_after   = n_match_after,
    n_SI_before     = n_SI_before,
    n_SI_after      = n_SI_after,
    n_SC_before     = n_SC_before,
    n_SC_after      = n_SC_after,
    stringsAsFactors = FALSE
  )
}

summary_table <- do.call(rbind, stats)
print(summary_table)
