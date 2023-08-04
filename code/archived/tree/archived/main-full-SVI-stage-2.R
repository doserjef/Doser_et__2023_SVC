# main-uni-gaussian: fits a Gaussian univariate SVC model analogous to spSVC 
rm(list = ls())
library(spAbundance)

# Load data ---------------------------------------------------------------
load('data/tree/stage-2-data.rda')

# Prep the model ----------------------------------------------------------
# Priors
curr.indx <- which(data.list.2$y != 0)
dist.fia <- dist(data.list.2$coords[curr.indx, ])
mean.dist <- mean(dist.fia)
min.dist <- 5
max.dist <- max(dist.fia)
prior.list <- list(beta.normal = list(mean = 0, var = 100),
		   phi.unif = list(a = 3 / max.dist, b = 3 / min.dist), 
                   sigma.sq = list(a = 2, b = 2), 
                   tau.sq = c(2, 1))
inits.list <- list(beta = 0, sigma.sq = 5, tau.sq = 0.2, phi = 3 / mean.dist)
tuning.list <- list(phi = c(0.3))

# Biggest
n.batch <- 2000
batch.length <- 25
n.burn <- 30000
n.thin <- 20
n.chains <- 1

# Big
# n.batch <- 500
# batch.length <- 25
# n.burn <- 7500
# n.thin <- 5

# n.batch <- 50
# batch.length <- 25
# n.burn <- 500
# n.thin <- 1

out <- svcAbund(formula = ~ scale(vpd) + I(scale(vpd)^2) + 
		            scale(ppt) + I(scale(ppt)^2) + 
			    scale(elev) + I(scale(elev)^2), 
		  data = data.list.2, priors = prior.list, 
		  inits = inits.list, tuning = tuning.list, svc.cols = c(1),
	          n.neighbors = 5, cov.model = 'exponential', NNGP = TRUE,
	          n.batch = n.batch, batch.length = batch.length, family = 'Gaussian-hurdle',
	          n.burn = n.burn, accept.rate = 0.43, n.thin = n.thin, 
	          n.chains = n.chains, n.report = 1, n.omp.threads = 5)

save(out, file = 'results/tree/stage-2-sugar-maple-full-SVI-samples.rda')

