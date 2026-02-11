### Megan Copeland
### Plot Simulation Results with folded null

# load library
library(dplyr)

# read in mean delta r values from empirical data (see plot_full.r)
emp_means <- read.csv("emp_means_pruned.csv")

# get families and paths to thei results
families <- list(
  Asteraceae    = "../pruned.tr_sim/aster_sims/aster_neutral_chain_rep",
  Brassicaceae  = "../pruned.tr_sim/brass_sims/brass_neutral_chain_rep",
  Solanaceae    = "../pruned.tr_sim/sol_sims/sol_neutral_chain_rep",
  Fabaceae      = "../pruned.tr_sim/fab_sims/fab_neutral_chain_rep")

# loop over families
for (i in 1:length(families)) {
  
  # get current family name
  curr_fam <- families[[i]]  
  rep_means <- list()
  
  # loop over sim replicates
  for (j in 1:100) {
    # read file and toss out half of rows as burn in
    df <- read.csv(paste0(curr_fam, j, ".csv"))
    n <- nrow(df)
    df_postburn <- df[(n/2 + 1):n, ]
    
    # delta r calculation
    delta_df <- data.frame(
      fusion  = df_postburn$desc1 - df_postburn$desc2,
      fission = df_postburn$asc1  - df_postburn$asc2,
      wgd     = df_postburn$pol1  - df_postburn$pol2,
      demi    = df_postburn$dem1  - df_postburn$dem2)
    
    # take mean of delta r
    mns <- colMeans(delta_df)
    rep_means[[j]] <- mns
  }
  
  # rbind mean delta r from replicates, take absolute value of those (fold null)
  mean_df <- as.data.frame(do.call(rbind, rep_means))
  mean_df_abs <- as.data.frame(lapply(mean_df, abs))
  
  # get empirical mean for the current family
  emp_sub <- emp_means[emp_means$Family == names(families)[i],]
  
  par(mfrow=c(2,2), mar=c(4,4,3,1))
  
  # loop over transition parameters for the current family
  for (k in 1:length(mean_df_abs)) {
    # get the current parameter, grab values for that it, get density 
    curr_param <- names(mean_df_abs)[k]
    x <- mean_df_abs[[curr_param]]
    dens <- density(x)
    
    # take absolute value of the empirical mean delta r for that transition
    mean_val <- abs(emp_sub$MeanDeltaR[emp_sub$Transition == curr_param])
    
    # calculate p-value
    r <- sum(x >= mean_val)
    n_sim <- length(x)
    pval <- (r + 1) / (n_sim + 1)
    
    # get min and max values for plotting
    x_min <- min(dens$x, mean_val)
    x_max <- max(dens$x, mean_val)
    y_max <- max(dens$y)
    
    # plot
    plot(dens,
         main = bquote(.(names(families)[i]) ~ "|" * Delta * "r|" ~ .(curr_param)),
         xlab = expression("|" * Delta * "r|"),
         ylab = "Density",
         type = "n",
         xlim = c(x_min, x_max),
         ylim = c(0, y_max * 1.2))
    
    polygon(dens, col = "lightblue", border = "black", lwd = 1.2)
    abline(v = mean_val, col = "red", lwd = 2)
    
    # p-value annotation (always top-left)
    legend("topleft",
           legend = paste0("p = ", signif(pval, 3)),
           bty = "n", text.col = "red", cex = 0.9)
  }
}
