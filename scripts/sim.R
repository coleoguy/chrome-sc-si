# ==============================
# Neutral Δr Simulation - Asteraceae
# ==============================

library(ape)
library(coda)
library(diversitree)
library(chromePlus)
library(phangorn)
library(parallel)

# --- Inputs -------------------------------------------------------------------
tree <- read.tree("../data/aster.tree") ## REPLACE
tree <- nnls.tree(cophenetic(tree), tree, rooted = TRUE)
dat  <- read.csv("../data/full_phenodat.csv")

# Filter to Solanaceae
dat <- dat[dat$Family == "Asteraceae", ] ## REPLACE
dat$Haploid.Chrom.. <- as.numeric(dat$Haploid.Chrom..)
dat <- dat[!is.na(dat$Haploid.Chrom..), ]

# --- Empirical frequencies of SC vs SI ----------------------------------------
tab_emp <- table(dat$mating_system)
freqs_emp <- tab_emp / sum(tab_emp)
emp_min_freq <- min(freqs_emp)   # rarer state frequency

# --- Mean empirical transition rate -------------------------------------------
q_emp <- 0.1607298   # REPLACE 

# --- Helper: simulate neutral binary trait ------------------------------------
simulate_neutral_binary <- function(tree, q = q_emp, emp_min_freq, tol = 0.10) {
  repeat {
    sim <- rTraitDisc(
      tree,
      model = matrix(c(-q, q, q, -q), 2, 2),
      states = c("0","1")
    )
    tab <- table(sim)
    freqs <- tab / sum(tab)
    
    if (length(freqs) == 2) {
      sim_min_freq <- min(freqs)
      # accept only if rarer state falls within ±10% of empirical value
      if (sim_min_freq >= (emp_min_freq * (1 - tol)) &&
          sim_min_freq <= (emp_min_freq * (1 + tol))) {
        return(as.numeric(sim))  # numeric vector (0/1)
      }
    }
  }
}

# --- Parallel replicate function ----------------------------------------------
run_rep <- function(i) {
  cat("Running replicate", i, "\n")

  # Resolve polytomies & scale tree
  tree_i <- ape::multi2di(tree, random = TRUE)
  root_depth <- max(node.depth.edgelength(tree_i))
  tree_i$edge.length <- tree_i$edge.length / root_depth
  
  # Neutral trait
  neutral_trait <- simulate_neutral_binary(
    tree_i,
    q = q_emp,
    emp_min_freq = emp_min_freq,
    tol = 0.10
  )
  
  # Build states
  range_chr <- c(pmax(1, min(dat$Haploid.Chrom..) - 1),
                 max(dat$Haploid.Chrom..) + 1)
  
  state_df <- data.frame(
    Species = dat$genus_species,
    Haploid = dat$Haploid.Chrom..,
    Neutral = neutral_trait,
    stringsAsFactors = FALSE
  )
  
  state_mat <- chromePlus::datatoMatrix(
    state_df,
    range = range_chr,
    hyper = TRUE
  )
  
  state_vec <- apply(state_mat, 1, which.max)
  names(state_vec) <- rownames(state_mat)
  state_vec <- state_vec[tree_i$tip.label]
  
  # Likelihood & MCMC
  lik0 <- diversitree::make.mkn(
    tree_i,
    states = state_vec,
    k      = ncol(state_mat),
    strict = FALSE,
    control = list(method = "ode", compiled = TRUE)
  )
  
  lik <- chromePlus::constrainMkn(
    data       = state_mat,
    lik        = lik0,
    hyper      = TRUE,
    polyploidy = FALSE,
    constrain  = list(drop.demi = FALSE, drop.poly = FALSE)
  )
  
  # pilot run
  pilot <- diversitree::mcmc(
    lik,
    x.init      = runif(length(argnames(lik))),
    prior       = make.prior.exponential(2),
    nsteps      = 100,
    w           = 1,
    print.every = 0
  )
  
  w <- diff(apply(
    pilot[11:100, 2:(length(argnames(lik)) + 1)],
    2, quantile, c(0.05, 0.95))
  )
  
  # main chain
  chain_unit <- diversitree::mcmc(
    lik,
    x.init      = runif(length(argnames(lik))),
    prior       = make.prior.exponential(2),
    nsteps      = 200,   
    w           = w,
    print.every = 0
  )
  
  # Convert to per-Myr
  rate_cols <- !(colnames(chain_unit) %in% c("i", "p"))
  chain_unit[, rate_cols] <- chain_unit[, rate_cols] / root_depth
  
  # Δr calculation
  delta_df <- data.frame(
    fusion  = chain_unit$desc1 - chain_unit$desc2,
    fission = chain_unit$asc1  - chain_unit$asc2,
    wgd     = chain_unit$pol1  - chain_unit$pol2,
    demi    = chain_unit$dem1  - chain_unit$dem2
  )
  
  delta_folded <- as.data.frame(lapply(delta_df, function(x) abs(x)))
  
  # Save outputs
  write.csv(chain_unit, paste0("../results/aster_neutral_chain_rep", i, ".csv"), row.names = FALSE) ##REPLACE
  write.csv(delta_df,  paste0("../results/aster_neutral_delta_rep", i, ".csv"), row.names = FALSE) ##REPLACE
  
  # cleanup
  rm(tree_i, root_depth, neutral_trait, range_chr,
     state_df, state_mat, state_vec, lik0, lik,
     pilot, w, chain_unit, delta_df, delta_folded)
  gc()
  
  return(NULL)
}

# --- Run in parallel ----------------------------------------------------------
n_reps  <- 100   # adjust as needed
n_cores <- 48    # adjust based on your machine

mclapply(1:n_reps, run_rep, mc.cores = n_cores)
