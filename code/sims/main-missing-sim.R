# main-missing-sim.R: script to perform the simulation study comparing the
#                     five candidate models in their ability to predict
#                     a species-environment relationship from the "Missing"
#                     scenario.
# Author: Jeffrey W. Doser
rm(list = ls())
library(spOccupancy)

# Functions ---------------------------------------------------------------
logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

# Parameters for simulations ----------------------------------------------
# Number of data sets simulated for each scenario
n.sims <- 50
# Set seed to generate same values
set.seed(1111)
# Random seeds for each data set
my.seeds <- sample(1:100000, n.sims, replace = FALSE)
# Simulation scenarios
model.types <- c('linear', 'quadratic', 'stratum', 'interaction', 'svc')
# Total number of candidate models
n.models <- length(model.types)
# Number of spatial locations
J.x <- 20
J.y <- 20
J <- J.x * J.y

# Simulation setup --------------------------------------------------------
waic.vals <- array(NA, dim = c(n.sims, n.models))
beta.effect.true.vals <- array(NA, dim = c(n.sims, J))
beta.effect.mean.vals <- array(NA, dim = c(n.sims, J, n.models))

# MCMC Info ---------------------------
n.samples <- 15000
batch.length <- 25
n.batch <- n.samples / batch.length
n.burn <- 10000
n.thin <- 5
n.chains <- 1
accept.rate <- 0.43

# Set parameters and values -----------------------------------------------
# Matrix of spatial locations
s.x <- seq(0, 1, length.out = J.x)
s.y <- seq(0, 1, length.out = J.y)
coords <- as.matrix(expand.grid(s.x, s.y))
# Get strata for each cell
strata <- ifelse(coords[, 1] <= .33 & coords[, 2] <= .33, 1,
		 ifelse(coords[, 1] <= .33 & coords[, 2] <= .67, 2,
			ifelse(coords[, 1] <= .33 & coords[, 2] <= 1, 3,
			       ifelse(coords[, 1] <= .67 & coords[, 2] <= .33, 4,
				      ifelse(coords[, 1] <= .67 & coords[, 2] <= .67, 5,
					     ifelse(coords[, 1] <= .67 & coords[, 2] <= 1, 6,
						    ifelse(coords[, 1] <= 1 & coords[, 2] <= .33, 7,
							   ifelse(coords[, 1] <= 1 & coords[, 2] <= .67, 8,
								  9))))))))
# Main covariate
x.1 <- seq(from = -5, to = 5, length.out = J)
x.interaction <- c(t(matrix(x.1, J.x, J.y)))
beta.0 <- 0
beta.linear <- 0
beta.quadratic <- -0.2
beta.strata <- runif(length(unique(strata)), -1, 1)
beta.interaction <- 0.4
beta.miss.int <- 0.5
sigma.sq <- 2
phi <- 3 / .8
Sigma <- spBayes::mkSpCov(coords, as.matrix(sigma.sq), as.matrix(0), as.matrix(phi), 'exponential')
x.miss.int <- MASS::mvrnorm(1, rep(0, J), Sigma)

# Do the simulations ------------------------------------------------------
for (j in 1:n.sims) {
  print(paste("Currently on simulation set ", j, " out of ", n.sims, sep = ''))
  set.seed(my.seeds[j])
  # Simulate the data set -------------------------------------------------
  logit <- function(theta, a = 0, b = 1){log((theta-a)/(b-theta))}
  logit.inv <- function(z, a = 0, b = 1){b-(b-a)/(1+exp(z))}
  
  # Form detection covariate (if any) -------------------------------------
  # Assume constant detection
  n.rep <- rep(3, J)
  n.rep.max <- max(n.rep)
  alpha <- c(0.5)
  n.alpha <- length(alpha)
  X.p <- array(NA, dim = c(J, n.rep.max, n.alpha))
  # Get index of surveyed replicates for each site. 
  rep.indx <- list()
  for (k in 1:J) {
    rep.indx[[k]] <- sample(1:n.rep.max, n.rep[k], replace = FALSE)
  }
  X.p[, , 1] <- 1
  if (n.alpha > 1) {
    for (i in 2:n.alpha) {
      for (k in 1:J) {
        X.p[k, rep.indx[[k]], i] <- rnorm(n.rep[k])
      } # j
    } # i
  }
  
  # Latent Occupancy Process ----------------------------------------------
  psi <- logit.inv(beta.0 + beta.linear * x.1 + beta.miss.int * x.miss.int * x.1)
  z <- rbinom(J, 1, psi)
  
  # Data Formation --------------------------------------------------------
  p <- matrix(NA, nrow = J, ncol = n.rep.max)
  y <- matrix(NA, nrow = J, ncol = n.rep.max)
  for (k in 1:J) {
    p[k, rep.indx[[k]]] <- logit.inv(X.p[k, rep.indx[[k]], ] %*% as.matrix(alpha))
    y[k, rep.indx[[k]]] <- rbinom(n.rep[k], 1, p[k, rep.indx[[k]]] * z[k])
  } # j
  data.list <- list(y = y, 
		    occ.covs = data.frame(x = x.1, 
					  x.interaction = x.interaction, 
					  strata = strata), 
		    coords = coords)
  beta.effect.true.vals[j, ] <- beta.linear + beta.miss.int * x.miss.int

  for (i in 1:n.models) {
    print(paste("Currently on model ", i, " out of ", n.models, sep = ''))
    curr.model <- model.types[i]
    phi.tune <- 0.5
    # Priors
    prior.list <- list(beta.normal = list(mean = 0, var = 2.72), 
    		       alpha.normal = list(mean = 0, var = 2.72), 
    		       sigma.sq.ig = list(a = 2, b = 1), 
                       phi.unif = list(a = 3 / 1, b = 3 / 0.1)) 
    # Starting values
    z.init <- apply(y, 1, function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
    inits.list <- list(beta = 0, alpha = 0, sigma.sq = 1, phi = 3 / 0.5,
    		       z = z.init)
    # Tuning
    tuning.list <- list(phi = phi.tune) 
    # Determine occupancy formula
    if (curr.model == 'linear') {
      occ.formula <- ~ x
    } else if (curr.model == 'quadratic') {
      occ.formula <- ~ x + I(x^2)
    } else if (curr.model == 'stratum') {
      occ.formula <- ~ x + x:factor(strata)
    } else if (curr.model == 'interaction') {
      occ.formula <- ~ x + x:x.interaction
    } else if (curr.model == 'svc') {
      occ.formula <- ~ x
    }

    if (curr.model != 'svc') {
      out <- PGOcc(occ.formula = occ.formula,
		   det.formula = ~ 1,
		   data = data.list,
		   n.samples = n.batch * batch.length,
		   inits = inits.list,
		   priors = prior.list,
		   verbose = TRUE,
		   n.report = 1000,
		   n.burn = n.burn, 
		   n.thin = n.thin,
		   n.chains = 1)

    } else {
      n.samples <- 30000
      batch.length <- 25
      n.batch <- n.samples / batch.length
      n.burn <- 20000
      n.thin <- 10
      out <- svcPGOcc(occ.formula = occ.formula,
		      det.formula = ~ 1,
		      svc.cols = 2,
		      data = data.list,
		      n.batch = n.batch,
		      batch.length = batch.length,
		      inits = inits.list,
		      priors = prior.list,
		      accept.rate = 0.43, 
		      cov.model = 'exponential', 
		      tuning = tuning.list,
		      n.omp.threads = 2, # Change as desired. 
		      verbose = TRUE,
		      NNGP = TRUE,
		      n.neighbors = 5, 
		      n.report = 25,
		      n.burn = n.burn, 
		      n.thin = n.thin,
		      n.chains = 1)

    }
    waic.vals[j, i] <- waicOcc(out)[3]
    if (curr.model == 'linear') {
      beta.effect.mean.vals[j, , i] <- mean(out$beta.samples[, 2])
    }
    if (curr.model == 'quadratic') {
      beta.effect.mean.vals[j, , i] <- mean(out$beta.samples[, 2]) + mean(out$beta.samples[, 3]) * x.1
    }
    if (curr.model == 'stratum') {
      beta.mean.vals <- apply(out$beta.samples[, -1], 2, mean) 
      beta.mean.vals[-1] <- beta.mean.vals[-1] + beta.mean.vals[1]
      beta.effect.mean.vals[j, , i] <- beta.mean.vals[strata]
    }
    if (curr.model == 'interaction') {
      beta.effect.mean.vals[j, , i] <- mean(out$beta.samples[, 2]) + mean(out$beta.samples[, 3]) * x.interaction
    }
    if (curr.model == 'svc') {
      beta.effect.mean.vals[j, , i] <- apply(getSVCSamples(out)[[1]], 2, mean)
    }
  } # i (n.scenarios)
} # j (n.sims)

# Save results ------------------------------------------------------------
save(waic.vals, beta.effect.true.vals, beta.effect.mean.vals, 
     file = 'results/sim-true-missing.rda')
