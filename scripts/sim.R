# ==============================
# Neutral Δr Simulation - Asteraceae
# ==============================

library(ape)
library(coda)
library(diversitree)
library(chromePlus)
library(phangorn)
library(parallel)

# read in files
tree <- read.tree("../data/aster.tree") ## REPLACE
tree <- nnls.tree(cophenetic(tree), tree, rooted = TRUE)
dat <- read.csv("../data/full_phenodat.csv")

# filter to selected family
dat <- dat[dat$Family == "Asteraceae", ] ## REPLACE
matches <- intersect(dat$genus_species, tree$tip.label)
dat <- dat[dat$genus_species %in% matches, ]
dat$Haploid.Chrom.. <- as.numeric(dat$Haploid.Chrom..)
dat <- dat[!is.na(dat$Haploid.Chrom..), ]

# get frequencies of each mating system type
ms_emp <- table(dat$mating_system) # make table of mating system counts
freqs_emp <- ms_emp / sum(ms_emp) # make table of mating system frequency
emp_min_freq <- min(freqs_emp)   # get rarer state frequency

# Mean empirical transition rate 
q_emp <- 0.1607298   # REPLACE 

# function to simulate neutral binary trait 
simulate_neutral_binary <- function(tree, q = q_emp, emp_min_freq, tol = 0.10) {
  repeat {
    sim <- rTraitDisc(tree, model = matrix(c(-q, q, q, -q), 2, 2), states = c("0","1"))
    
    tab <- table(sim) # get table with counts for each tip state
    freqs <- tab / sum(tab) # get frequency of those states
    sim_min_freq <- min(freqs) # get the frequency of the rarer state
    
    if (sim_min_freq >= (emp_min_freq * (1 - tol)) && # check that it falls within our tol
        sim_min_freq <= (emp_min_freq * (1 + tol))) {
      return(as.numeric(sim))  # numeric vector (0/1)
    }
  }
}


# function to run simulations in parallel 
run_rep <- function(i) {
  cat("Running replicate", i, "\n")

  # resolve polytomies & scale tree
  tree_i <- ape::multi2di(tree, random = TRUE)
  root_depth <- max(node.depth.edgelength(tree_i))
  tree_i$edge.length <- tree_i$edge.length / root_depth
  
  # get neutral trait
  neutral_trait <- simulate_neutral_binary(tree_i, q = q_emp, 
                                           emp_min_freq = emp_min_freq, tol = 0.10)
  
  # build states
  range_chr <- c(min(dat$Haploid.Chrom..) - 1, max(dat$Haploid.Chrom..) + 1)
  
  state_df <- data.frame(
    Species = dat$genus_species,
    Haploid = dat$Haploid.Chrom..,
    Neutral = neutral_trait)
  
  state_mat <- chromePlus::datatoMatrix(
    state_df,
    range = range_chr,
    hyper = TRUE)
  
  state_vec <- apply(state_mat, 1, which.max)
  names(state_vec) <- rownames(state_mat)
  state_vec <- state_vec[tree_i$tip.label]
  
  # Likelihood & MCMC
  lik0 <- diversitree::make.mkn(
    tree_i, states = state_vec, 
    k = ncol(state_mat),
    strict = FALSE, 
    control = list(method = "ode", compiled = TRUE))
  
  lik <- chromePlus::constrainMkn(
    data = state_mat,
    lik = lik0,
    hyper = TRUE,
    polyploidy = FALSE,
    constrain = list(drop.demi = FALSE, drop.poly = FALSE))
  
  # pilot run
  pilot <- diversitree::mcmc(
    lik,
    x.init = runif(length(argnames(lik))),
    prior = make.prior.exponential(2),
    nsteps = 100,
    w = 1,
    print.every = 0)
  
  w <- diff(apply(
    pilot[11:100, 2:(length(argnames(lik)) + 1)],
    2, quantile, c(0.05, 0.95)))
  
  # main chain
  chain_unit <- diversitree::mcmc(
    lik,
    x.init      = runif(length(argnames(lik))),
    prior       = make.prior.exponential(2),
    nsteps      = 200,   
    w           = w,
    print.every = 0)
  
  # Convert to per-Myr
  rate_cols <- !(colnames(chain_unit) %in% c("i", "p"))
  chain_unit[, rate_cols] <- chain_unit[, rate_cols] / root_depth
  
  # save outputs
  write.csv(chain_unit, paste0("../results/aster_neutral_chain_rep", i, ".csv"), row.names = FALSE) ##REPLACE

  # cleanup
  rm(tree_i, root_depth, neutral_trait, range_chr,
     state_df, state_mat, state_vec, lik0, lik,
     pilot, w, chain_unit)
  gc()
  
  return(NULL)
}

# --- Run in parallel ----------------------------------------------------------
n_reps  <- 100   # adjust as needed
n_cores <- 1    # adjust based on your machine

mclapply(1:n_reps, run_rep, mc.cores = n_cores)
