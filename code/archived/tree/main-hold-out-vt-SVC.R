# main-uni-gaussian: fits a Gaussian univariate SVC model analogous to spSVC 
rm(list = ls())
library(tidyverse)
library(sf)
library(spAbundance)

# Load data ---------------------------------------------------------------
load('data/tree/spAbundance-data.rda')

# Generate a hold-out set and remove 25% of locations ---------------------
load('data/tree/pred-indx.rda')
tcc.pred <- data.list$covs[pred.indx, ]
coords.0 <- data.list$coords[pred.indx, ]
y.0 <- data.list$y[pred.indx]
data.list$y <- data.list$y[-pred.indx]
data.list$coords <- data.list$coords[-pred.indx, ]
data.list$covs <- data.list$covs[-pred.indx, , drop = FALSE]


library(MBA)
library(fields)
library(geoR)
par(mfrow=c(1,2))
out.lm <- lm(sqrt(data.list$y) ~ data.list$covs$tcc)
curr.res <- out.lm$residuals
vario.1.raw <- variog(coords = data.list$coords, data = sqrt(data.list$y))
par(mfrow = c(1, 2))
plot(vario.1.raw, pch=16)
vario.1.res <- variog(coords = data.list$coords, data = curr.res)
plot(vario.1.res, pch = 16)
par(mfrow = c(1, 1))


# Prep the model ----------------------------------------------------------
dist.fia <- dist(data.list$coords)
# Priors
# bound.dists <- quantile(dist.fia, c(0.05, 0.95))
# 0.05: 236.6716
# 0.95: 2132.5627
# min.dist <- bound.dists[1]
# max.dist <- bound.dists[2]
min.dist <- min(dist.fia)
max.dist <- max(dist.fia)
prior.list <- list(beta.normal = list(mean = 0, var = 100),
		   phi.unif = list(a = 3 / max.dist, b = 3 / min.dist), 
                   sigma.sq = list(a = 2, b = 2), 
                   tau.sq = c(2, 1))
tuning.list <- list(phi = c(0.3))

# Biggest
n.batch <- 4000
batch.length <- 25
n.burn <- 50000
n.thin <- 50 
n.chains <- 1

# Big
n.batch <- 500
batch.length <- 25
n.burn <- 7500
n.thin <- 5
n.chains <- 1

# n.batch <- 100
# batch.length <- 25
# n.burn <- 500
# n.thin <- 2
# n.chains <- 1

data.list$y <- sqrt(data.list$y)
data.list$covs$tcc <- data.list$covs$tcc / 100

out <- svcAbund(formula = ~ tcc,
		data = data.list, priors = prior.list, 
	        tuning = tuning.list, svc.cols = c(1, 2),
	        n.neighbors = 5, cov.model = 'exponential', NNGP = TRUE,
	        n.batch = n.batch, batch.length = batch.length, family = 'Gaussian',
	        n.burn = n.burn, accept.rate = 0.43, n.thin = n.thin, 
	        n.chains = n.chains, n.report = 1, n.omp.threads = 1)

summary(out)
y.rep.means <- apply(out$y.rep.samples, 2, mean)
plot(data.list$y, y.rep.means, pch = 19)
abline(0, 1)

# Predict at hold-out locations -------------------------------------------
J.0 <- nrow(coords.0)
X.0 <- matrix(1, J.0, ncol(out$X))
X.0[, 2] <- tcc.pred / 100
# X.0[, 2] <- (tcc.pred / 100 - mean(data.list$covs$tcc)) / sd(data.list$covs$tcc)
colnames(X.0) <- c("(Intercept)", "tcc")
out.pred <- predict(out, X.0, coords.0, n.omp.threads = 1)

# png("tmp.png")
y.0.means <- apply(out.pred$y.0.samples^2, 2, mean)
plot(y.0, y.0.means, pch = 19)
abline(0, 1)
# dev.off()

save(out, out.pred, file = 'results/tree/hold-out-vermont-SVC-results.rda')

