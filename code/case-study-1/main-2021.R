# main-2021.R: this script contains code to run a spatial occupancy model
#              for each of the 51 forest bird species in 2021, which we use 
#              to generate estimates of true occurrence in 2021 that are subsequently
#              used for an assessment of model predictive performance for the 
#              different SVC trend models with AUC.
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
# curr.sp <- 'RHWO'
if(length(args) == 0) base::stop('Need to tell spOccupancy the current species')
print(curr.sp)

# Format data for spOccupancy ---------------------------------------------
load("data/case-study-1/spOccupancy-2021-data.rda")
# Loads the data.list object, which is formatted for spOccupancy.
# Species info
sp.names <- dimnames(data.list.2021$y)[[1]]
sp.indx <- which(sp.names == curr.sp)
data.list.2021$y <- data.list.2021$y[sp.indx, , ]
# Only use locations within the breeding range of the species.
data.list.2021$y <- data.list.2021$y[data.list.2021$range.ind[sp.indx, ] == 1, ]
data.list.2021$occ.covs <- data.list.2021$occ.covs[data.list.2021$range.ind[sp.indx, ] == 1, ]
data.list.2021$det.covs$day <- data.list.2021$det.covs$day[data.list.2021$range.ind[sp.indx, ] == 1]
data.list.2021$det.covs$rep.val <- data.list.2021$det.covs$rep.val[data.list.2021$range.ind[sp.indx, ] == 1, ]
data.list.2021$coords <- data.list.2021$coords[data.list.2021$range.ind[sp.indx, ] == 1, ]
data.list.2021$range.ind <- NULL

# Run the model -----------------------------------------------------------
n.batch <- 4000
n.burn <- 50000
n.thin <- 50
n.chains <- 1

# Specify priors and initial values.
tuning.list <- list(phi = 0.5)
dist.bbs <- dist(data.list.2021$coords)
min.dist <- min(dist.bbs)
low.dist <- quantile(dist.bbs, 0.15)
max.dist <- max(dist.bbs)
mean.dist <- mean(dist.bbs)
inits <- list(phi = 3 / mean.dist, sigma.sq = 1, sigma.sq.p = 1)

priors <- list(phi.unif = c(3 / max.dist, 3 / low.dist))

out <- spPGOcc(occ.formula = ~ scale(tmax),
                 det.formula = ~ scale(day) + I(scale(day)^2) +
			         scale(rep.val) + I(scale(rep.val)^2),
                 data = data.list.2021,
	         inits = inits,
		 priors = priors,
                 n.batch = n.batch,
                 batch.length = 25,
                 accept.rate = 0.43,
                 cov.model = "exponential",
                 tuning = tuning.list,
                 n.omp.threads = 1, # Change as desired
                 verbose = TRUE,
                 NNGP = TRUE,
                 n.neighbors = 5,
                 n.report = 100,
                 n.burn = n.burn,
                 n.thin = n.thin,
                 n.chains = n.chains)

# Save z.samples only, which is what is used for a predictive performance assessment
z.samples <- out$z.samples

# Save result -------------------------------------------------------------
save(z.samples, file = paste("results/case-study-1/spPGOcc-2021-", curr.sp, ".R", sep = ''))

