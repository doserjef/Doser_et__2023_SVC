# main-uni-gaussian: fits a Gaussian univariate SVC model analogous to spSVC 
rm(list = ls())
library(tidyverse)
library(sf)
library(spAbundance)

# Load data ---------------------------------------------------------------
load('data/tree/spAbundance-data.rda')

# Generate a hold-out set and remove 25% of locations ---------------------
set.seed(19191)
J <- nrow(data.list$coords)
load('data/tree/pred-indx.rda')
tcc.pred <- data.list$covs[pred.indx, ]
coords.0 <- data.list$coords[pred.indx, ]
y.0 <- data.list$y[pred.indx]
data.list$y <- data.list$y[-pred.indx]
data.list$coords <- data.list$coords[-pred.indx, ]
data.list$covs <- data.list$covs[-pred.indx, , drop = FALSE]

# Prep the model ----------------------------------------------------------
# Priors
min.dist <- 236.6716
max.dist <- 2132.5627
prior.list <- list(beta.normal = list(mean = 0, var = 100),
		   phi.unif = c(a = 3 / max.dist, b = 3 / min.dist), 
                   sigma.sq = c(a = 2, b = 2), 
                   tau.sq = c(2, 1))
# Based on an initial run with 12,500 iterations.
# inits.list <- list(beta = c(6, 0.04), 
# 		   sigma.sq = c(.5), 
# 		   tau.sq = 0.6, 
# 		   phi = c(.5))
inits.list <- list(sigma.sq = 1)
tuning.list <- list(phi = c(0.5))

# Biggest
n.batch <- 5000
batch.length <- 25
n.burn <- 50000
n.thin <- 75
n.chains <- 1

# Big
# n.batch <- 500
# batch.length <- 25
# n.burn <- 7500
# n.thin <- 5
# n.chains <- 1

data.list$y <- sqrt(data.list$y)
data.list$covs$tcc <- data.list$covs$tcc / 100

out <- abund(formula = ~ tcc,
	     data = data.list, priors = prior.list, 
	     inits = inits.list, tuning = tuning.list,
	     n.batch = n.batch, batch.length = batch.length, family = 'Gaussian',
	     n.burn = n.burn, accept.rate = 0.43, n.thin = n.thin, 
	     n.chains = n.chains, n.report = 1)

summary(out)
# y.rep.means <- apply(out$y.rep.samples, 2, mean)
# plot(data.list$y, y.rep.means, pch = 19)
# abline(0, 1)

# Predict at hold-out locations -------------------------------------------
J.0 <- nrow(coords.0)
X.0 <- matrix(1, J.0, ncol(out$X))
colnames(X.0) <- c("(Intercept)", "tcc")
X.0[, 2] <- tcc.pred
out.pred <- predict(out, X.0)

# y.0.means <- apply(out.pred$y.0.samples^2, 2, mean)
# plot(y.0, y.0.means, pch = 19)
# abline(0, 1)

save(out, out.pred, file = 'results/tree/hold-out-east-us-non-spatial-results.rda')

