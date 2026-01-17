# =============================
# Neutral vs Empirical Δr (PDF-safe Δr + p-values)
# =============================

library(dplyr)

emp_means <- read.csv("emp_means_pruned.csv")

fam_map <- list(
  Asteraceae    = "../pruned_sim/aster_sims/aster_neutral_chain_rep",
  Brassicaceae  = "../pruned_sim/brass_sims/brass_neutral_chain_rep",
  Solanaceae    = "../pruned_sim/sol_sims/sol_neutral_chain_rep",
  Fabaceae      = "../pruned_sim/fab_sims/fab_neutral_chain_rep"
)

# Loop over families
for (fam in names(fam_map)) {
  
  cat("Processing:", fam, "\n")
  rep_means <- list()
  
  for (i in 1:100) {
    file <- paste0(fam_map[[fam]], i, ".csv")
    df <- read.csv(file, check.names = FALSE)
    n <- nrow(df)
    df_post <- df[(n/2 + 1):n, ]
    
    # Δr calculation
    delta_df <- data.frame(
      fusion  = df_post$desc1 - df_post$desc2,
      fission = df_post$asc1  - df_post$asc2,
      wgd     = df_post$pol1  - df_post$pol2,
      demi    = df_post$dem1  - df_post$dem2
    )
    mns <- colMeans(delta_df, na.rm = TRUE)
    rep_means[[i]] <- mns
  }
  
  # Combine replicates
  mean_df <- as.data.frame(do.call(rbind, rep_means))
  mean_df_abs <- as.data.frame(lapply(mean_df, abs))
  emp_sub <- subset(emp_means, Family == fam)
  
  par(mfrow=c(2,2), mar=c(4,4,3,1))
  
  for (col in names(mean_df_abs)) {
    x <- mean_df_abs[[col]]
    d <- density(x, na.rm=TRUE)
    mean_val <- abs(emp_sub$MeanDeltaR[emp_sub$Transition == col])
    
    # Monte Carlo p-value
    r <- sum(x >= mean_val)
    n_sim <- length(x)
    pval <- (r + 1) / (n_sim + 1)
    
    x_min <- min(d$x, mean_val, na.rm=TRUE)
    x_max <- max(d$x, mean_val, na.rm=TRUE)
    y_max <- max(d$y, na.rm=TRUE)
    
    # Plot (PDF-safe Δr)
    plot(d,
         main = bquote(.(fam) ~ "|" * Delta * "r|" ~ .(col)),
         xlab = expression("|" * Delta * "r|"),
         ylab = "Density",
         type = "n",
         xlim = c(x_min, x_max),
         ylim = c(0, y_max * 1.2))
    
    polygon(d, col = "lightblue", border = "black", lwd = 1.2)
    abline(v = mean_val, col = "red", lwd = 2)
    
    # p-value annotation (always top-left)
    legend("topleft",
           legend = paste0("p = ", signif(pval, 3)),
           bty = "n", text.col = "red", cex = 0.9)
  }
}
