# main-uni-gaussian: fits a Gaussian univariate SVC model analogous to spSVC 
rm(list = ls())
library(tidyverse)
library(sf)
library(spAbundance)

# Load data ---------------------------------------------------------------
load('data/tree/spAbundance-data.rda')

# Prep the model ----------------------------------------------------------
# Priors
dist.fia <- dist(data.list$coords)
mean.dist <- mean(dist.fia)
min.dist <- min(dist.fia)
max.dist <- max(dist.fia)
prior.list <- list(beta.normal = list(mean = 0, var = 100),
		   phi.unif = list(a = 3 / max.dist, b = 3 / min.dist), 
                   sigma.sq = list(a = 2, b = 2), 
                   tau.sq = c(2, 1))
# Based on an initial run with 12,500 iterations.
inits.list <- list(beta = c(0.4, 0.04), 
		   sigma.sq = c(.5, .009), 
		   tau.sq = 0.6, 
		   phi = c(.94, .02))
tuning.list <- list(phi = c(1))

# Biggest
n.batch <- 4000
batch.length <- 25
n.burn <- 50000
n.thin <- 25 
n.chains <- 1

# Big
# n.batch <- 500
# batch.length <- 25
# n.burn <- 7500
# n.thin <- 5
# n.chains <- 1

# data.list$covs$tcc <- data.list$covs$tcc / 100
data.list$y <- sqrt(data.list$y)

out <- svcAbund(formula = ~ tcc,
		  data = data.list, priors = prior.list, 
		  inits = inits.list, tuning = tuning.list, svc.cols = c(1, 2),
	          n.neighbors = 5, cov.model = 'exponential', NNGP = TRUE,
	          n.batch = n.batch, batch.length = batch.length, family = 'Gaussian',
	          n.burn = n.burn, accept.rate = 0.43, n.thin = n.thin, 
	          n.chains = n.chains, n.report = 1, n.omp.threads = 1) # TODO: 

summary(out)
y.rep.means <- apply(out$y.rep.samples, 2, mean)
plot(data.list$y, y.rep.means, pch = 19)
abline(0, 1)

save(out, file = 'results/tree/vermont-SVC-results.rda')

