# main-svc-GRSP.R: script to run the SVC model for GRSP case study
# Author: Jeffrey W. Doser
rm(list = ls())
library(spOccupancy)

# Load data ---------------------------------------------------------------
load("data/case-study-2/GRSP-spOcc-data.rda")

# Set values for running the model ----------------------------------------
# n.batch <- 20000
# n.burn <- 200000
# n.thin <- 300
# n.chains <- 3
# Bigger
# n.batch <- 4000
# n.burn <- 70000
# n.thin <- 5
# Big
n.batch <- 4000
n.burn <- 50000
n.thin <- 50
n.chains <- 3

# Restricting prior bound on phi
dist.bbs <- dist(data.GRSP$coords)
low.dist <- quantile(dist.bbs, .1)
high.dist <- quantile(dist.bbs, 0.99)
priors <- list(phi.unif = list(a = 3 / high.dist, b = 3 / low.dist)) 

tuning.list <- list(phi = c(0.5), rho = 0.5, sigma.sq = 1)
inits.list <- list(beta = 0, sigma.sq = 1)
data.GRSP$occ.covs$grass.dev <- data.GRSP$occ.covs$grass.dev + abs(min(data.GRSP$occ.covs$grass.dev))
data.GRSP$occ.covs$crop.dev <- data.GRSP$occ.covs$crop.dev + abs(min(data.GRSP$occ.covs$crop.dev))

out <- svcTPGOcc(occ.formula = ~ grass.dev + crop.dev + 
                                 scale(tmax) + I(scale(tmax)^2) + 
				 scale(grass.mean) + scale(crop.mean),
                det.formula = ~ scale(day) + I(scale(day)^2) + 
                                factor(replicate) + (1 | year.det),
                svc.cols = c(1, 2, 3),
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

save(out, file = 'results/case-study-2/svc-GRSP-model-results.rda')
