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

# Read in BBS data from 
load("data/spOcc-bbs-data.rda")
# Loads the data.list object, which is formatted for spOccupancy.
# Species info
sp.names <- dimnames(data.list$y)[[1]]
sp.indx <- which(sp.names == curr.sp)
data.list$y <- data.list$y[sp.indx, , , ]

# Load in model results
load(paste("results/svcTPGOcc-", curr.sp, '.R', sep = ''))

# Load prediction data ----------------------------------------------------
load('data/pred-coordinates.rda')
J.0 <- nrow(coords.0)
n.time.max <- ncol(data.list$y)
X.0 <- array(NA, dim = c(J.0, n.time.max, dim(out$X)[3]))
X.0[, , 1] <- 1
unique.years <- unique(c(data.list$occ.covs$years))
for (t in 1:n.time.max) {
  X.0[, t, 2] <- unique.years[t]
}
X.0[, , 2] <- (X.0[, , 2] - mean(c(data.list$occ.covs$years))) / sd(c(data.list$occ.covs$years))
t.cols <- 1:n.time.max

vals <- split(1:J.0, ceiling(seq_along(1:J.0) / 100))
psi.samples <- array(NA, dim = c(out$n.post * out$n.chains, J.0, n.time.max))
w.samples <- array(NA, dim = c(out$n.post * out$n.chains, J.0, length(out$svc.cols)))
for (j in 1:length(vals)) {
  print(paste("Currently on set ", j, " out of ", length(vals), sep = ''))
  curr.indx <- vals[[j]]
  # Can change the number of threads as desired.
  out.pred <- predict(out, X.0[curr.indx, , , drop = FALSE], coords.0[curr.indx, ],
        	      t.cols = t.cols, n.omp.threads = 5, verbose = FALSE)
  psi.samples[, curr.indx, ] <- out.pred$psi.0.samples
  w.samples[, curr.indx, ] <- aperm(out.pred$w.0.samples, c(1, 3, 2))
}

psi.quants <- apply(psi.samples, c(2, 3), quantile, c(0.025, 0.5, 0.975))
w.trend.quants <- apply(w.samples[, , 2], 2, quantile, c(0.025, 0.5, 0.975))

save(psi.quants, w.trend.quants, coords.0, 
     file = paste('results/predict-svcTPGOcc-', curr.sp, sep = ''))
