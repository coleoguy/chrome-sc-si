# *Author*    Megan Copeland, Heath Blackmon
# *Purpose*   Aggregate MCMC runs per family, sample post-burn, then plot each family manually
# *Inputs*    ../results/{aster_nd, brass_nd, sol_nd, fab_nd}/*_rep#.csv

library(coda)

# =========================
# Config
# =========================
base_dir   <- "../results"
dirs       <- c("aster_nd","brass_nd","sol_nd","fab_nd")
fam_labels <- c("Asteraceae","Brassicaceae","Solanaceae","Fabaceae")

samples_per_run <- 10   # rows to randomly sample after burn-in

param_pairs <- list(
  fusion  = c("desc1","desc2"),
  fission = c("asc1","asc2"),
  wgd     = c("pol1","pol2")
)

cols <- c(
  fusion  = rgb(228,26,28,125, maxColorValue=255),
  fission = rgb(55,126,184,125, maxColorValue=255),
  wgd     = rgb(77,175,74,125, maxColorValue=255)
)

# Store results per family
res_list <- list()
ci_list  <- list()
props_all <- data.frame()

# =========================
# Loop over families: collect sampled Δ
# =========================
for(i in 1:length(dirs)) {
  family <- fam_labels[i]
  res_dir <- file.path(base_dir, dirs[i])
  files <- list.files(res_dir, pattern="_rep\\d+\\.csv$", full.names=TRUE)
  if(length(files) == 0) {
    cat("No files found for", family, "\n")
    next
  }
  
  res <- data.frame()
  
  for(f in files) {
    chain <- read.csv(f, check.names=FALSE)
    n <- nrow(chain)
    if(n < 2) next
    
    # Drop burn-in
    post_idx <- (floor(n/2) + 1):n
    if(length(post_idx) < samples_per_run) next
    idx <- sample(post_idx, samples_per_run)
    
    deltas <- data.frame(matrix(nrow=samples_per_run, ncol=0))
    for(pn in names(param_pairs)) {
      pp <- param_pairs[[pn]]
      if(all(pp %in% names(chain))) {
        v1 <- as.numeric(chain[[pp[1]]])[idx]
        v2 <- as.numeric(chain[[pp[2]]])[idx]
        deltas[[pn]] <- v1 - v2
      }
    }
    res <- rbind(res, deltas)
  }
  
  res <- res[, colSums(!is.na(res)) > 0, drop=FALSE]
  if(ncol(res) == 0) next
  
  ci <- HPDinterval(as.mcmc(res))
  props <- colMeans(res > 0, na.rm=TRUE)
  
  res_list[[family]] <- res
  ci_list[[family]]  <- ci
  
  props_all <- rbind(props_all,
                     data.frame(Family=family,
                                param=names(props),
                                p_delta_gt0=as.numeric(props)))
  
  # Console summary
  cat("\nFamily:", family, "| Files processed:", length(files), "\n")
  for(nm in colnames(res)) {
    cat(sprintf("%-8s  n=%d,  P(Δ>0)=%.3f,  HPD=[%.5f, %.5f]\n",
                nm, sum(!is.na(res[[nm]])), mean(res[[nm]] > 0, na.rm=TRUE),
                ci[nm,1], ci[nm,2]))
  }
}

# =========================
# Plot families individually (manual tweaking here)
# =========================

# Asteraceae
res <- res_list[["Asteraceae"]]
ci  <- ci_list[["Asteraceae"]]

dens_list <- apply(res, 2, density, na.rm=TRUE)
y_max <- max(sapply(dens_list, function(d) max(d$y)))

plot(dens_list[[1]], xlim=c(-0.08,0.8),
     ylim=c(-0.3*y_max, 1.05*y_max), main="Asteraceae",
     xlab=expression(Delta[r]~"rate statistic"))
abline(v=0, lty=2, lwd=2)

for(i in 1:ncol(res)){
  polygon(density(res[,i]), col=cols[i])
}
ys <- seq(from=-0.1*y_max, by=-0.1*y_max, length.out=ncol(res))
for(i in 1:ncol(res)){
  lines(ci[i,], rep(ys[i],2), lwd=5, col=cols[i])
}
legend("topright", legend=colnames(res),
       fill=cols[colnames(res)], border=NA, bty="n")


# Fabaceae
res <- res_list[["Fabaceae"]]
ci  <- ci_list[["Fabaceae"]]

dens_list <- apply(res, 2, density, na.rm=TRUE)
y_max <- max(sapply(dens_list, function(d) max(d$y)))

plot(dens_list[[1]], xlim=c(-0.02,0.03),
     ylim=c(-0.3*y_max, 1.05*y_max), main="Fabaceae",
     xlab=expression(Delta[r]~"rate statistic"))
abline(v=0, lty=2, lwd=2)

for(i in 1:ncol(res)){
  polygon(density(res[,i]), col=cols[i])
}
ys <- seq(from=-0.1*y_max, by=-0.1*y_max, length.out=ncol(res))
for(i in 1:ncol(res)){
  lines(ci[i,], rep(ys[i],2), lwd=5, col=cols[i])
}
legend("topright", legend=colnames(res),
       fill=cols[colnames(res)], border=NA, bty="n")


# Brassicaceae
res <- res_list[["Brassicaceae"]]
ci  <- ci_list[["Brassicaceae"]]

dens_list <- apply(res, 2, density, na.rm=TRUE)
y_max <- max(sapply(dens_list, function(d) max(d$y)))

plot(dens_list[[1]], xlim=c(-0.03,0.25),
     ylim=c(-0.3*y_max, 1.05*y_max), main="Brassicaceae",
     xlab=expression(Delta[r]~"rate statistic"))
abline(v=0, lty=2, lwd=2)

for(i in 1:ncol(res)){
  polygon(density(res[,i]), col=cols[i])
}
ys <- seq(from=-0.1*y_max, by=-0.1*y_max, length.out=ncol(res))
for(i in 1:ncol(res)){
  lines(ci[i,], rep(ys[i],2), lwd=5, col=cols[i])
}
legend("topright", legend=colnames(res),
       fill=cols[colnames(res)], border=NA, bty="n")


# Solanaceae
res <- res_list[["Solanaceae"]]
ci  <- ci_list[["Solanaceae"]]

dens_list <- apply(res, 2, density, na.rm=TRUE)
y_max <- max(sapply(dens_list, function(d) max(d$y)))

plot(dens_list[[1]], xlim=c(-0.01,0.08),
     ylim=c(-0.3*y_max, 1.05*y_max), main="Solanaceae",
     xlab=expression(Delta[r]~"rate statistic"))
abline(v=0, lty=2, lwd=2)

for(i in 1:ncol(res)){
  polygon(density(res[,i]), col=cols[i])
}
ys <- seq(from=-0.1*y_max, by=-0.1*y_max, length.out=ncol(res))
for(i in 1:ncol(res)){
  lines(ci[i,], rep(ys[i],2), lwd=5, col=cols[i])
}
legend("topright", legend=colnames(res),
       fill=cols[colnames(res)], border=NA, bty="n")
