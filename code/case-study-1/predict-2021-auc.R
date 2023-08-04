# predict-2021-auc.R: script to predict occupancy probability for each of the 
#                     51 eastern forest bird species in 2021 as an assessment
#                     of out-of-sample predictive performance in the future
#                     (i.e., a short term forecast).
# Author: Jeffrey W. Doser
rm(list = ls())
library(spOccupancy)
# Package used to calculate WAIC
library(pROC)

# Get info from command line ----------------------------------------------
# This code is to extract the current species name from the command line
# to easily run the script for different species
args <- commandArgs(trailingOnly = TRUE)
# Current species
curr.sp <- args[1]
# Alternatively, if not running the script from the command line:
# curr.sp <- 'YBCU'
if(length(args) == 0) base::stop('Need to tell spOccupancy the current species')
print(curr.sp)

# Read in full BBS data ---------------------------------------------------
load("data/case-study-1/final-spOccupancy-data.rda")
# Loads the data.list object, which is formatted for a multi-species model, 
# so we need to do a bit of manipulation below.
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

# Read in 2021 data -------------------------------------------------------
load("data/case-study-1/spOccupancy-2021-data.rda")
# Loads the data.list object, which is formatted for spOccupancy.
data.list.2021$y <- data.list.2021$y[sp.indx, , ]
# Only use locations within the breeding range of the species.
data.list.2021$y <- data.list.2021$y[data.list.2021$range.ind[sp.indx, ] == 1, ]
data.list.2021$occ.covs <- data.list.2021$occ.covs[data.list.2021$range.ind[sp.indx, ] == 1, ]
data.list.2021$det.covs$day <- data.list.2021$det.covs$day[data.list.2021$range.ind[sp.indx, ] == 1]
data.list.2021$det.covs$rep.val <- data.list.2021$det.covs$rep.val[data.list.2021$range.ind[sp.indx, ] == 1, ]
data.list.2021$coords <- data.list.2021$coords[data.list.2021$range.ind[sp.indx, ] == 1, ]
data.list.2021$range.ind <- NULL

# Read in estimated occupancy values for 2021 -----------------------------
# These are the actual occupancy values predicted from a model that only uses
# data from 2021.
load(paste0("results/case-study-1/spPGOcc-2021-", curr.sp, ".R"))
# Calculate median values as the "True" value.
z.2021.median <- apply(z.samples, 2, median)

# Constant Model ----------------------------------------------------------
load(paste0("results/case-study-1/stPGOcc-", curr.sp, ".rda"))
coords.0 <- data.list.2021$coords
J.0 <- nrow(coords.0)
# Predicting for 1 year
n.time.max <- 1
X.0 <- array(NA, dim = c(J.0, n.time.max, dim(out$X)[3]))
X.0[, , 1] <- 1
X.0[, , 2] <- 2021
# Scaling covariate values by those values used to fit the model.
X.0[, , 2] <- (X.0[, , 2] - mean(c(data.list$occ.covs$years))) / 
	       sd(c(data.list$occ.covs$years))
tmax.0 <- data.list.2021$occ.covs$tmax
X.0[, , 3] <- (tmax.0 - mean(c(data.list$occ.covs$tmax))) / sd(c(data.list$occ.covs$tmax))

out.pred <- predict(out, X.0, coords.0, forecast = TRUE)

# Calculate AUC -----------------------------------------------------------
z.0.samples <- out.pred$z.0.samples
# Will calculate an AUC value for each species and each iteration of the 
# MCMC in order to get an AUC estimate with uncertainty.
n.samples <- out$n.post * out$n.chains
auc.constant.vals <- rep(NA, n.samples)
for (j in 1:n.samples) {
  auc.constant.vals[j] <- auc(response = z.2021.median, predictor = z.0.samples[j, , ])
} # j (iteration)

rm(out, out.pred)
gc()

# Max temp interaction model ----------------------------------------------	
load(paste0("results/case-study-1/stPGOcc-tmax-", curr.sp, ".rda"))
X.0 <- array(NA, dim = c(J.0, n.time.max, dim(out$X)[3]))
X.0[, , 1] <- 1
X.0[, , 2] <- 2021
# Scaling covariate values by those values used to fit the model.
X.0[, , 2] <- (X.0[, , 2] - mean(c(data.list$occ.covs$years))) /
	       sd(c(data.list$occ.covs$years))
tmax.0 <- data.list.2021$occ.covs$tmax
X.0[, , 3] <- (tmax.0 - mean(c(data.list$occ.covs$tmax))) / sd(c(data.list$occ.covs$tmax))
X.0[, , 4] <- X.0[, , 2] * X.0[, , 3]

out.pred <- predict(out, X.0, coords.0, forecast = TRUE)

# Calculate AUC -----------------------------------------------------------
z.0.samples <- out.pred$z.0.samples
# Will calculate an AUC value for each species and each iteration of the 
# MCMC in order to get an AUC estimate with uncertainty.
n.samples <- out$n.post * out$n.chains
auc.tmax.vals <- rep(NA, n.samples)
for (j in 1:n.samples) {
  auc.tmax.vals[j] <- auc(response = z.2021.median, predictor = z.0.samples[j, , ])
} # j (iteration)

rm(out, out.pred)
gc()

# BCR model ---------------------------------------------------------------
load(paste0("results/case-study-1/stPGOcc-bcr-", curr.sp, ".rda"))
X.0 <- array(NA, dim = c(J.0, n.time.max, dim(out$X)[3]))
X.0[, , 1] <- 1
year.scaled <- (2021 - mean(c(data.list$occ.covs$years))) /
	       sd(c(data.list$occ.covs$years))
tmax.0 <- data.list.2021$occ.covs$tmax
# Scaling covariate values by those values used to fit the model.
X.0[, , 2] <- (tmax.0 - mean(c(data.list$occ.covs$tmax))) / sd(c(data.list$occ.covs$tmax))
unique.bcrs <- sort(unique(data.list$occ.covs$BCR))
X.bcr <- matrix(0, J.0, length(unique.bcrs))
for (i in 1:ncol(X.bcr)) {
  X.bcr[which(data.list.2021$occ.covs$BCR == unique.bcrs[i]), i] <- 1 
}
X.0[, , -c(1:2)] <- X.bcr * year.scaled

out.pred <- predict(out, X.0, coords.0, forecast = TRUE)

# Calculate AUC -----------------------------------------------------------
z.0.samples <- out.pred$z.0.samples
# Will calculate an AUC value for each species and each iteration of the 
# MCMC in order to get an AUC estimate with uncertainty.
n.samples <- out$n.post * out$n.chains
auc.bcr.vals <- rep(NA, n.samples)
for (j in 1:n.samples) {
  auc.bcr.vals[j] <- auc(response = z.2021.median, predictor = z.0.samples[j, , ])
} # j (iteration)

rm(out, out.pred)
gc()

# SVC model ---------------------------------------------------------------
load(paste0("results/case-study-1/svcTPGOcc-", curr.sp, ".rda"))
X.0 <- array(NA, dim = c(J.0, n.time.max, dim(out$X)[3]))
X.0[, , 1] <- 1
X.0[, , 2] <- 2021
# Scaling covariate values by those values used to fit the model.
X.0[, , 2] <- (X.0[, , 2] - mean(c(data.list$occ.covs$years))) /
	       sd(c(data.list$occ.covs$years))
tmax.0 <- data.list.2021$occ.covs$tmax
X.0[, , 3] <- (tmax.0 - mean(c(data.list$occ.covs$tmax))) / sd(c(data.list$occ.covs$tmax))

out.pred <- predict(out, X.0, coords.0, forecast = TRUE)

# Calculate AUC -----------------------------------------------------------
z.0.samples <- out.pred$z.0.samples
# Will calculate an AUC value for each species and each iteration of the
# MCMC in order to get an AUC estimate with uncertainty.
n.samples <- out$n.post * out$n.chains
auc.svc.vals <- rep(NA, n.samples)
for (j in 1:n.samples) {
  auc.svc.vals[j] <- auc(response = z.2021.median, predictor = z.0.samples[j, , ])
} # j (iteration)

auc.df <- data.frame(constant = auc.constant.vals,
		     tmax = auc.tmax.vals, 
		     bcr = auc.bcr.vals,
		     svc = auc.svc.vals)
apply(auc.df, 2, mean)

# Save AUC values to hard drive -------------------------------------------
save(auc.df, file = paste0('results/case-study-1/auc-', curr.sp, '.rda'))
