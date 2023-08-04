rm(list = ls())
library(MASS)
library(spBayes)
library(tidyverse)
library(viridis)
library(ggpubr)
library(devtools)
load_all("~/Dropbox/DFKZ21/spOccupancy")
# Simulate the occupancy process ------------------------------------------
# set.seed(500)
# Intercept + linear + quadratic + stratum + interaction + missing interaction
# Subroutines -----------------------------------------------------------
logit <- function(theta, a = 0, b = 1){log((theta-a)/(b-theta))}
logit.inv <- function(z, a = 0, b = 1){b-(b-a)/(1+exp(z))}
# Simulate a data set based on the actual covariate data that you have.

# Main covariate
beta.0 <- 0
beta.linear <- 0.5
beta.quadratic <- -0.5
beta.int <- -0.5
# sigma.sq <- 3
# phi <- 3 / .8
# Sigma <- mkSpCov(coords, as.matrix(sigma.sq), as.matrix(0), as.matrix(phi), 'exponential')
# x <- mvrnorm(1, rep(0, J), Sigma)
load("data/rel-import/EATO-spOccupancy-data.rda")
coords <- data.EATO$coords
# coord.x <- runif(nrow(data.EATO$coords))
# coord.y <- runif(nrow(data.EATO$coords))
# coords <- cbind(coord.x, coord.y)
# coords <- coords
J <- nrow(coords)
curr.points <- sample(1:J, 500, replace = FALSE)
coords <- coords[curr.points, ]
forest <- data.EATO$occ.covs$forest[curr.points]
tmax <- data.EATO$occ.covs$tmax[curr.points]
J <- length(curr.points)
# forest <- rnorm(J)
# tmax <- rnorm(J)

# Form detection covariate (if any) -------------------------------------
n.rep <- rep(4, J)
n.rep.max <- max(n.rep)
alpha <- c(0.5, -0.3)
n.alpha <- length(alpha)
X.p <- array(NA, dim = c(J, n.rep.max, n.alpha))
# Get index of surveyed replicates for each site. 
rep.indx <- list()
for (j in 1:J) {
  rep.indx[[j]] <- sample(1:n.rep.max, n.rep[j], replace = FALSE)
}
X.p[, , 1] <- 1
if (n.alpha > 1) {
  for (i in 2:n.alpha) {
    for (j in 1:J) {
      X.p[j, rep.indx[[j]], i] <- rnorm(n.rep[j])
    } # j
  } # i
}

# Simulate spatial random effect ----------------------------------------
svc.cols <- c(1)
p.svc <- length(svc.cols)
cov.model <- 'exponential'
w.mat <- matrix(NA, J, p.svc)
# phi <- c(3 / mean(dist(coords)))
phi <- c(3 / quantile(dist(coords), 0.7))
sigma.sq <- 1
theta <- as.matrix(phi)
for (i in 1:p.svc) {
  Sigma <- mkSpCov(coords, as.matrix(sigma.sq[i]), as.matrix(0), theta[i, ], cov.model)
  # Random spatial process
  w.mat[, i] <- mvrnorm(1, rep(0, J), Sigma)
}
w <- c(t(w.mat))

# Latent Occupancy Process ----------------------------------------------
# TODO: change for different forms
psi <- logit.inv(beta.0 + w + beta.linear * scale(tmax) + beta.quadratic * I(scale(tmax)^2) + 
                 beta.linear * scale(forest) + beta.quadratic * I(scale(forest)^2)) 
# psi <- logit.inv(beta.0 + w + beta.linear * scale(forest) + beta.quadratic * scale(forest)^2)
weights <- rep(1, J)
z <- rbinom(J, weights, psi)

# Data Formation --------------------------------------------------------
p <- matrix(NA, nrow = J, ncol = n.rep.max)
y <- matrix(NA, nrow = J, ncol = n.rep.max)
for (j in 1:J) {
  p[j, rep.indx[[j]]] <- logit.inv(X.p[j, rep.indx[[j]], ] %*% as.matrix(alpha))
  y[j, rep.indx[[j]]] <- rbinom(n.rep[j], 1, p[j, rep.indx[[j]]] * z[j])
} # j


data.list <- list(y = y,
		  occ.covs = data.frame(forest = forest, 
				    tmax = tmax),
		  det.covs = list(det.cov.1 = X.p[, , 2]),
                  coords = coords, 
                  weights = weights)
dist.coords <- dist(coords)
priors <- list(phi.unif = list(a = c(3 / max(dist.coords), 3 / quantile(dist.coords, 0.5), 
				     3 / quantile(dist.coords, 0.5)),
			       b = c(3 / quantile(dist.coords, 0.6), 3 / quantile(dist.coords, 0.5), 
				     3 / quantile(dist.coords, 0.5))))
priors <- list(phi.unif = list(a = 3 / quantile(dist.coords, 0.95), b = 3 / quantile(dist.coords, 0.1)))
Sigma <- matrix(c(1, 0.8, 0.1, 0.8, 1, 0.1, 0.1, 0.1, 1), 3, 3)
# lambda <- matrix(c(1, 0.8, 0, 1), 2, 2)
lambda <- t(chol(Sigma))
inits <- list(lambda = lambda)

n.batch <- 800
batch.length <- 25
n.chains <- 1
n.burn <- 10000
n.thin <- 10

n.batch <- 400
batch.length <- 25
n.chains <- 1
n.burn <- 5000
n.thin <-5 
# This works, can easily recover all parameters here, so there
# clearly isn't a fundamental problem with the data locations.
# out.2 <- spPGOcc(occ.formula = ~ scale(tmax) + I(scale(tmax)^2),
# 		det.formula = ~ det.cov.1,
# 		data = data.list,
# # 		priors = priors,
# 		tuning = list(phi = 0.5),
# 		n.batch = n.batch,
# 		batch.length = batch.length,
# 		n.chains = n.chains,
# 		n.neighbors = 5,
# 		n.burn = n.burn,
# 		n.thin = n.thin,
#                 n.report = 10)

# psi.means.2 <- apply(out.2$psi.samples, 2, mean)
# plot(psi, psi.means.2, pch = 19)
# abline(0, 1)

# n.batch <- 1600
# batch.length <- 25
# n.chains <- 1
# n.burn <- 20000
# n.thin <- 10

out <- svcPGOcc(occ.formula = ~ scale(tmax) + scale(forest),
		  det.formula = ~ det.cov.1,
		data = data.list,
  		priors = priors,
		inits = inits, 
		svc.cols = c(1, 2, 3),
		tuning = list(phi = 0.5, sigma.sq = 0.1, lambda = 0.5),
		n.batch = n.batch,
		batch.length = batch.length,
		n.chains = n.chains,
		n.neighbors = 5,
		n.burn = n.burn,
		n.thin = n.thin,
                n.report = 10)
svc.samples <- getSVCSamples(out)
int.means <- apply(svc.samples[[1]], 2, mean)
tmax.means <- apply(svc.samples[[2]], 2, mean)
forest.means <- apply(svc.samples[[3]], 2, mean)
# svc.means <- svc.samples$x[10, ]

plot(w, int.means, pch = 19)
abline(0, 1)

# Tmax
svc.true <- beta.linear + 2 * beta.quadratic * scale(tmax)

plot(svc.true, tmax.means, pch = 19)
abline(0, 1)

# Forest
forest.svc.true <- beta.linear + 2 * beta.quadratic * scale(forest)

plot(forest.svc.true, forest.means, pch = 19)
abline(0, 1)

psi.means <- apply(out$psi.samples, 2, mean)
plot(psi, psi.means, pch = 19)
abline(0, 1)

plot.df <- data.frame(x = coords[, 1], 
		      y = coords[, 2], 
		      true = svc.true, 
		      int = int.means,
		      psi = psi,
		      int.true = w,
		      svc.est = tmax.means)
svc.est.map <- ggplot(data = plot.df, aes(x = x, y = y, color = svc.est)) +
  geom_point(size = 5) +
  scale_color_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC',
  	               na.value = NA) +
  theme_bw(base_size = 10) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = 'Easting', y = 'Northing', title = 'SVC Estimate', color = '') +
  guides(color = 'none') +
  theme(text = element_text(family="LM Roman 10"),
        axis.ticks.x = element_blank(),
	plot.background = element_rect(color = 'transparent', fill = NA),
	axis.text.x = element_blank(),
	axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

svc.true.map <- ggplot(data = plot.df, aes(x = x, y = y, color = svc.true)) +
  geom_point(size = 5) +
  scale_color_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC',
  	               na.value = NA) +
  theme_bw(base_size = 10) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = 'Easting', y = 'Northing', title = 'SVC True', color = '') +
  guides(color = 'none') +
  theme(text = element_text(family="LM Roman 10"),
        axis.ticks.x = element_blank(),
	plot.background = element_rect(color = 'transparent', fill = NA),
	axis.text.x = element_blank(),
	axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggarrange(svc.est.map, svc.true.map)

int.est.map <- ggplot(data = plot.df, aes(x = x, y = y, color = int)) +
  geom_point(size = 5) +
  scale_color_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC',
  	               na.value = NA) +
  theme_bw(base_size = 10) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = 'Easting', y = 'Northing', title = 'Intercept Estimate', color = '') +
  guides(color = 'none') +
  theme(text = element_text(family="LM Roman 10"),
        axis.ticks.x = element_blank(),
	plot.background = element_rect(color = 'transparent', fill = NA),
	axis.text.x = element_blank(),
	axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

int.true.map <- ggplot(data = plot.df, aes(x = x, y = y, color = int.true)) +
  geom_point(size = 5) +
  scale_color_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC',
  	               na.value = NA) +
  theme_bw(base_size = 10) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = 'Easting', y = 'Northing', title = 'Intercept Truth', color = '') +
  guides(color = 'none') +
  theme(text = element_text(family="LM Roman 10"),
        axis.ticks.x = element_blank(),
	plot.background = element_rect(color = 'transparent', fill = NA),
	axis.text.x = element_blank(),
	axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggarrange(int.est.map, int.true.map)
