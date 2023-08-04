# pred-svcTPGOcc.R: this script contains code to predict occurrence
#                   across the eastern US for a single-species using 
#                   a multi-season spatially-varying coefficients occupancy
#                   model. 
# Author: Jeffrey W. Doser
rm(list = ls())
library(spOccupancy)

# Get info from command line ----------------------------------------------
# This code is to extract the current species name from the command line
# to easily run the script for different species
args <- commandArgs(trailingOnly = TRUE)
# Current species
curr.sp <- args[1]
# Alternatively, if not running the script from the command line:
# curr.sp <- 'BWWA'
if(length(args) == 0) base::stop('Need to tell spOccupancy the current species')
print(curr.sp)

# Read in full BBS data ---------------------------------------------------
load("data/case-study-1/final-spOccupancy-data.rda")
# Loads the data.list object, which is formatted for spOccupancy.
# Species info
sp.names <- dimnames(data.list$y)[[1]]
sp.indx <- which(sp.names == curr.sp)
data.list$y <- data.list$y[sp.indx, , , ]
# Only use locations within the breeding range of the species.
data.list$y <- data.list$y[data.list$range.ind[sp.indx, ] == 1, , ]
data.list$occ.covs$years <- data.list$occ.covs$years[data.list$range.ind[sp.indx, ] == 1, ]
data.list$occ.covs$BCR <- data.list$occ.covs$BCR[data.list$range.ind[sp.indx, ] == 1]
data.list$occ.covs$tmax <- data.list$occ.covs$tmax[data.list$range.ind[sp.indx, ] == 1]
data.list$det.covs$day <- data.list$det.covs$day[data.list$range.ind[sp.indx, ] == 1, ]
data.list$det.covs$obs <- data.list$det.covs$obs[data.list$range.ind[sp.indx, ] == 1, ]
data.list$det.covs$rep.val <- data.list$det.covs$rep.val[data.list$range.ind[sp.indx, ] == 1, , ]
data.list$det.covs$obs.first.year <- data.list$det.covs$obs.first.year[data.list$range.ind[sp.indx, ] == 1, ]
data.list$det.covs$year.det <- data.list$det.covs$year.det[data.list$range.ind[sp.indx, ] == 1, ]
data.list$coords <- data.list$coords[data.list$range.ind[sp.indx, ] == 1, ]
data.list$range.ind <- NULL

# Load in model results
load(paste("results/case-study-1/svcTPGOcc-", curr.sp, '.rda', sep = ''))

# Load prediction data ----------------------------------------------------
load('data/case-study-1/pred-coordinates.rda')
J.0 <- nrow(coords.0)
n.time.max <- ncol(data.list$y)
X.0 <- array(NA, dim = c(J.0, n.time.max, dim(out$X)[3]))
X.0[, , 1] <- 1
unique.years <- unique(c(data.list$occ.covs$years))
for (t in 1:n.time.max) {
  X.0[, t, 2] <- unique.years[t]
}
# Make sure to scale everything by values used to fit the model. 
X.0[, , 2] <- (X.0[, , 2] - mean(c(data.list$occ.covs$years))) / sd(c(data.list$occ.covs$years))
t.cols <- 1:n.time.max
X.0[, , 3] <- (tmax.0 - mean(c(data.list$occ.covs$tmax))) / sd(c(data.list$occ.covs$tmax))

# Predict at a subet of values at a time to speed things up ---------------
vals <- split(1:J.0, ceiling(seq_along(1:J.0) / 100))
psi.samples <- array(NA, dim = c(out$n.post * out$n.chains, J.0, n.time.max))
svc.samples <- array(NA, dim = c(out$n.post * out$n.chains, J.0))
for (j in 1:length(vals)) {
  print(paste("Currently on set ", j, " out of ", length(vals), sep = ''))
  curr.indx <- vals[[j]]
  # Can change the number of threads as desired.
  out.pred <- predict(out, X.0[curr.indx, , , drop = FALSE], coords.0[curr.indx, ],
        	      t.cols = t.cols, n.omp.threads = 5, verbose = FALSE)
  psi.samples[, curr.indx, ] <- out.pred$psi.0.samples
  tmp <- getSVCSamples(out, out.pred)
  svc.samples[, curr.indx] <- tmp[[2]]
  rm(tmp, out.pred)
  gc()
}

# Save quantiles and other summary information to reduce size of resulting objects
psi.quants <- apply(psi.samples, c(2, 3), quantile, c(0.025, 0.5, 0.975))
trend.prob.pos <- apply(svc.samples, 2, function(a) mean(a > 0))
trend.quants <- apply(svc.samples, 2, quantile, c(0.025, 0.5, 0.975))

# Save to hard drive ------------------------------------------------------
save(psi.quants, trend.quants, coords.0, trend.prob.pos,
     file = paste('results/case-study-1/predict-svcTPGOcc-', curr.sp, '.rda', sep = ''))
