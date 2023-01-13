# main-stPGOcc.R: this script contains code to run a spatial multi-season
#                 occupancy model for each of the 51 forest bird species
#                 in which we assume a constant trend across the eastern
#                 U.S. (the "Constant" model).
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
load("data/spOcc-bbs-data.rda")
# Loads the data.list object, which is formatted for spOccupancy.
# Species info
sp.names <- dimnames(data.list$y)[[1]]
sp.indx <- which(sp.names == curr.sp)
data.list$y <- data.list$y[sp.indx, , , ]
# Only use locations within the breeding range of the species.
# Load in species range information based on bird life data. 
load("data/bird-life-processed.rda")
data.list$y <- data.list$y[range.ind[sp.indx, ] == 1, , ]
data.list$occ.covs$years <- data.list$occ.covs$years[range.ind[sp.indx, ] == 1, ]
data.list$occ.covs$BCR <- data.list$occ.covs$BCR[range.ind[sp.indx, ] == 1]
data.list$det.covs$day <- data.list$det.covs$day[range.ind[sp.indx, ] == 1, ]
data.list$det.covs$obs <- data.list$det.covs$obs[range.ind[sp.indx, ] == 1, ]
data.list$det.covs$obs.first.year <- data.list$det.covs$obs.first.year[range.ind[sp.indx, ] == 1, ]
data.list$det.covs$year.det <- data.list$det.covs$year.det[range.ind[sp.indx, ] == 1, ]
data.list$coords <- data.list$coords[range.ind[sp.indx, ] == 1, ]

# Add in replicate as a covariate effect on detection
rep.val <- array(NA, dim = dim(data.list$y))
for (i in 1:dim(rep.val)[3]) {
  rep.val[, , i] <- i
}
data.list$det.covs$rep.val <- rep.val

# Run the model -----------------------------------------------------------
# Can change as desired.
# Really big
n.batch <- 2000
n.burn <- 30000
n.thin <- 20
# Big
# n.batch <- 500
# n.burn <- 7500
# n.thin <- 5
# Moderate
# n.batch <- 300
# n.burn <- 5000 
# n.thin <- 5
# Small
# n.batch <- 100
# n.burn <- 2000
# n.thin <- 1
# Smaller
# n.batch <- 50
# n.burn <- 1000
# n.thin <- 1
# Real Small
# n.batch <- 1
# n.burn <-0 
# n.thin <- 1

# Specify priors and initial values. 
tuning.list <- list(phi = 0.5, rho = 0.5, sigma.sq = 0.3)
dist.bbs <- dist(data.list$coords)
min.dist <- min(dist.bbs)
low.dist <- quantile(dist.bbs, 0.25)
max.dist <- max(dist.bbs)
mean.dist <- mean(dist.bbs)
inits <- list(phi = 3 / mean.dist, sigma.sq = 1, sigma.sq.p = 1, 
              rho = 0, sigma.sq.t = 1)

priors <- list(sigma.sq.unif = c(0, 10), 
               phi.unif = c(3 / max.dist, 3 / low.dist))

out <- stPGOcc(occ.formula = ~ scale(years),
                 det.formula = ~ scale(day) + scale(I(day^2)) + 
			         scale(rep.val) + scale(I(rep.val^2)) +
                                 factor(year.det) + (1 | obs),
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
                 n.chains = 1) 

# Save result -------------------------------------------------------------
save(out, file = paste("results/stPGOcc-", curr.sp, ".R", sep = ''))

