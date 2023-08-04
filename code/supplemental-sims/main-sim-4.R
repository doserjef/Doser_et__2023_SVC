# main-sim-4.R: code to perform the fourth simulation study to assess the 
#               effects of varying detection probability and number of 
#               replicate surveys in a spatially-varying coefficient occupancy
#               model. 
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
set.seed(10101)
# Spatial locations
J.x <- 20
J.y <- 20
J <- J.x * J.y
# Different number of replicates
n.rep.vals <- c(2, 3, 4)
# Occurrence coefficients -------------
beta <- c(0, 0)
# Detection coefficients --------------
alpha.vals <- c(logit(0.1), logit(0.3), logit(0.5), logit(0.7))
# Spatial parameters ------------------
sp <- TRUE
svc.cols <- c(1, 2)
p.svc <- length(svc.cols)
cov.model <- 'exponential'
# Different spatial variance values
sigma.sq <- c(0.5, 0.3)
phi <- c(3 / 0.5, 3 / .7)
# Total number of simulation scenarios
n.scenarios <- length(n.rep.vals) * length(alpha.vals)
scenario.vals <- expand.grid(n.rep = n.rep.vals, alpha = alpha.vals)
# Random seeds for each data set and simulation scenario
my.seeds <- sample(1:100000, n.sims * n.scenarios, replace = FALSE)

# Simulation setup --------------------------------------------------------
psi.true <- array(NA, dim = c(J, n.sims, n.scenarios))
w.true <- array(NA, dim = c(J, p.svc, n.sims, n.scenarios))
psi.mean.samples <- array(NA, dim = c(J, n.sims, n.scenarios))
psi.low.samples <- array(NA, dim = c(J, n.sims, n.scenarios))
psi.high.samples <- array(NA, dim = c(J, n.sims, n.scenarios))
w.mean.samples <- array(NA, dim = c(J, p.svc, n.sims, n.scenarios))
w.low.samples <- array(NA, dim = c(J, p.svc, n.sims, n.scenarios))
w.high.samples <- array(NA, dim = c(J, p.svc, n.sims, n.scenarios))

# MCMC Info ---------------------------
n.samples <- 15000
batch.length <- 25
n.batch <- n.samples / batch.length
n.burn <- 10000
n.thin <- 5
n.chains <- 3
accept.rate <- 0.43

# Simulate Data -----------------------------------------------------------
seed.indx <- 0
for (j in 1:n.sims) { 
  print(paste("Currently on simulation set ", j, " out of ", n.sims, sep = ''))
  for (i in 1:n.scenarios) {
    seed.indx <- seed.indx + 1
    set.seed(my.seeds[seed.indx])
    print(paste("Currently on scenario ", i, " out of ", n.scenarios, sep = ''))
    dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = rep(scenario.vals$n.rep[i], J), 
		  beta = beta, alpha = scenario.vals$alpha[i], 
		  sp = sp, svc.cols = svc.cols, cov.model = cov.model, 
		  sigma.sq = sigma.sq, phi = phi)
    psi.true[, j, i] <- dat$psi 
    w.true[1:J, 1:p.svc, j, i] <- dat$w
    # Prep the data for spOccupancy -------------------------------------------
    # Site x Replicate
    y <- dat$y
    # Occurrence Covariates
    X <- dat$X
    # Detection Covariates
    X.p <- dat$X.p
    # Coordinates
    coords <- dat$coords
    # Package all data into a list
    occ.covs <- X
    colnames(occ.covs) <- c('int', 'occ.cov.1')
    data.list <- list(y = y, occ.covs = occ.covs, coords = coords)
    # Priors
    prior.list <- list(beta.normal = list(mean = 0, var = 2.72), 
    		       alpha.normal = list(mean = 0, var = 2.72), 
    		       sigma.sq.ig = list(a = 2, b = 1), 
                       phi.unif = list(a = 3 / 1, b = 3 / 0.05)) 
    # Starting values
    z.init <- apply(y, 1, function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
    inits.list <- list(beta = 0, alpha = 0, sigma.sq = 1, phi = 3 / 0.5,
    		       z = z.init)
    # Tuning
    tuning.list <- list(phi = 1) 
    # SVC model
    out <- svcPGOcc(occ.formula = ~ occ.cov.1, 
		    det.formula = ~ 1,
		    svc.cols = svc.cols,
		    data = data.list,
		    n.batch = n.batch,
		    batch.length = batch.length,
		    inits = inits.list,
		    priors = prior.list,
		    accept.rate = 0.43, 
		    cov.model = 'exponential', 
		    tuning = tuning.list,
		    n.omp.threads = 4, # Change as desired
		    verbose = TRUE,
		    NNGP = TRUE,
		    n.neighbors = 5, 
		    n.report = 25,
		    n.burn = n.burn, 
		    n.thin = n.thin,
		    n.chains = 1) 
    psi.mean.samples[, j, i] <- apply(out$psi.samples, 2, mean)
    psi.low.samples[, j, i] <- apply(out$psi.samples, 2, quantile, 0.025)
    psi.high.samples[, j, i] <- apply(out$psi.samples, 2, quantile, 0.975)
    w.mean.samples[1:J, 1:p.svc, j, i] <- t(apply(out$w.samples, c(2, 3), mean))
    w.low.samples[1:J, 1:p.svc, j, i] <- t(apply(out$w.samples, c(2, 3), 
						   quantile, 0.025))
    w.high.samples[1:J, 1:p.svc, j, i] <- t(apply(out$w.samples, c(2, 3), 
						    quantile, 0.975))
  } # i (n.scenarios)
} # j (n.sims)

save(psi.mean.samples, psi.low.samples, psi.high.samples, scenario.vals, psi.true, 
     w.mean.samples, w.low.samples, w.high.samples, w.true,
     file = 'results/supplemental-sims/sim-4-results.rda')
