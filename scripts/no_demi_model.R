library(ape)
library(coda)
library(diversitree)
library(chromePlus)
library(phangorn)
library(parallel)

# --- Inputs -------------------------------------------------------------------
tree <- read.tree("../data/aster.tree") # change as needed 
tree <- nnls.tree(cophenetic(tree), tree, rooted = TRUE)
dat  <- read.csv("../data/full_phenodat.csv")  # SC/SI + haploid chrom data

# Filter dataset
dat <- dat[dat$Family == "Asteraceae", ] # change as needed
dat$Haploid.Chrom.. <- as.numeric(dat$Haploid.Chrom..)
dat <- dat[!is.na(dat$Haploid.Chrom..), ]

# Setup
n_reps <- 100  # number of tree replicates
n_cores <- 48 # use available cores minus one

run_rep <- function(i) {
  cat("Running replicate", i, "\n")
  
  # --- Randomly resolve polytomies and scale tree -----------------------------
  tree_i <- ape::multi2di(tree, random = TRUE)
  root_depth <- max(node.depth.edgelength(tree_i))
  tree_i$edge.length <- tree_i$edge.length / root_depth
  
  # --- Build states -----------------------------------------------------------
  range_chr <- c(pmax(1, min(dat$Haploid.Chrom..) - 1), max(dat$Haploid.Chrom..) + 1)
  
  state_df <- data.frame(
    Species = dat$genus_species,
    Haploid = dat$Haploid.Chrom..,
    p_SC    = ifelse(dat$mating_system == "sc", 1, 0),
    stringsAsFactors = FALSE)
  
  state_mat <- chromePlus::datatoMatrix(
    state_df,
    range = range_chr,
    hyper = TRUE)
  
  state_vec <- apply(state_mat, 1, which.max)
  names(state_vec) <- rownames(state_mat)
  state_vec <- state_vec[tree_i$tip.label]
  
  # --- Likelihood and MCMC ----------------------------------------------------
  lik0 <- diversitree::make.mkn(
    tree_i,
    states = state_vec,
    k      = ncol(state_mat),
    strict = FALSE,
    control = list(method = "ode", compiled = TRUE))
  
  lik <- chromePlus::constrainMkn(
    data       = state_mat,
    lik        = lik0,
    hyper      = TRUE,
    polyploidy = FALSE,
    constrain  = list(drop.demi = TRUE, drop.poly = FALSE))
  
  pilot <- diversitree::mcmc(
    lik,
    x.init      = runif(length(argnames(lik))),
    prior       = make.prior.exponential(2),
    nsteps      = 100,
    w           = 1,
    print.every = 0)
  
  w <- diff(apply(
    pilot[11:100, 2:(length(argnames(lik)) + 1)],
    2, quantile, c(0.05, 0.95)))
  
  chain_unit <- diversitree::mcmc(
    lik,
    x.init      = runif(length(argnames(lik))),
    prior       = make.prior.exponential(2),
    nsteps      = 1000,
    w           = w,
    print.every = 0)
  
  # Convert to per-Myr
  rate_cols <- !(colnames(chain_unit) %in% c("i", "p"))
  chain_unit[, rate_cols] <- chain_unit[, rate_cols] / root_depth
  
  # Save
  outfile <- paste0("../results/aster_nd_rep", i, ".csv") # change as needed
  write.csv(chain_unit, outfile, row.names = FALSE)
  
  # --- Cleanup to reduce memory load ------------------------------------------
  rm(tree_i, state_df, state_mat, state_vec, lik0, lik, pilot, chain_unit, w, range_chr, rate_cols, root_depth)
  gc()
  
  return(NULL)
}

mclapply(1:n_reps, run_rep, mc.cores = n_cores)
