# main-constant-GRSP.R: this script fits the constant model for the 
#                       Grasshopper Sparrow case study
# Author: Jeffrey W. Doser
rm(list = ls())
library(spOccupancy)

# Load data ---------------------------------------------------------------
load("data/case-study-2/GRSP-spOcc-data.rda")

# Set values for running the model ----------------------------------------
n.batch <- 4000
n.burn <- 50000
n.thin <- 50
n.chains <- 3

# Restricting prior bound on phi
dist.bbs <- dist(data.GRSP$coords)
low.dist <- quantile(dist.bbs, .1)
high.dist <- quantile(dist.bbs, 0.99)
priors <- list(phi.unif = c(a = 3 / high.dist, b = 3 / low.dist))

tuning.list <- list(phi = c(0.5), rho = 0.5, sigma.sq = 1)
inits.list <- list(beta = 0, sigma.sq = 1)

out <- stPGOcc(occ.formula = ~ scale(grass.dev) + scale(crop.dev) + 
                                 scale(tmax) + I(scale(tmax)^2) + 
				 scale(grass.mean) + scale(crop.mean),
                det.formula = ~ scale(day) + I(scale(day)^2) + 
                                factor(replicate) + (1 | year.det),
		ar1 = TRUE,
                data = data.GRSP,
                n.batch = n.batch,
                batch.length = 25,
		priors = priors,
		inits = inits.list, 
                accept.rate = 0.43,
                cov.model = "exponential",
                tuning = tuning.list,
                n.omp.threads = 3,
                verbose = TRUE,
                NNGP = TRUE,
                n.neighbors = 5,
                n.report = 1,
                n.burn = n.burn,
                n.thin = n.thin,
                n.chains = n.chains) 

save(out, file = 'results/case-study-2/constant-GRSP-model-results.rda')
