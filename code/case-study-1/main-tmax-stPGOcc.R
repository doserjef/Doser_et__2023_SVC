# main-tmax-stPGOcc.R: this script contains code to run a spatial multi-season
#                      occupancy model for each of the 51 forest bird species
#                      in which we fit a model with a trend that varies 
#                      with 30-year average maximum temperature
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
# curr.sp <- 'EABL'
if(length(args) == 0) base::stop('Need to tell spOccupancy the current species')
print(curr.sp)

# Format data for spOccupancy ---------------------------------------------
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

# Run the model -----------------------------------------------------------
# Can change as desired.
# Really big
n.batch <- 4000
n.burn <- 50000
n.thin <- 50
n.chains <- 3

# Specify priors and initial values. 
tuning.list <- list(phi = 0.5, rho = 0.5, sigma.sq = 0.3)
dist.bbs <- dist(data.list$coords)
min.dist <- min(dist.bbs)
low.dist <- quantile(dist.bbs, 0.15)
max.dist <- max(dist.bbs)
mean.dist <- mean(dist.bbs)
inits <- list(phi = 3 / mean.dist, sigma.sq = 1, sigma.sq.p = 1, 
              rho = 0, sigma.sq.t = 1)

priors <- list(phi.unif = c(3 / max.dist, 3 / low.dist))

out <- stPGOcc(occ.formula = ~ scale(years) * scale(tmax),
                 det.formula = ~ scale(day) + I(scale(day)^2) + 
			         scale(rep.val) + I(scale(rep.val)^2) +
                                 scale(year.det) + I(scale(year.det)^2),
                 data = data.list,
	         inits = inits,
		 priors = priors,
                 n.batch = n.batch,
                 batch.length = 25,
                 accept.rate = 0.43,
                 cov.model = "exponential",
                 tuning = tuning.list,
	         ar1 = TRUE,
                 n.omp.threads = 1, # Change as desired
                 verbose = TRUE,
                 NNGP = TRUE,
                 n.neighbors = 5,
                 n.report = 1,
                 n.burn = n.burn,
                 n.thin = n.thin,
                 n.chains = n.chains) 

# Save result -------------------------------------------------------------
save(out, file = paste("results/case-study-1/stPGOcc-tmax-", curr.sp, ".rda", sep = ''))

