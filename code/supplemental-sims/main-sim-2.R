# main-sim-2.R: code to perform the second supplemental simulation study to assess how 
#               sample size and the number of time periods influences 
#               the bias and precision of estimates from a multiseason
#               SVC occupancy model.
# Author: Jeffrey W. Doser
rm(list = ls())
library(devtools)
library(spOccupancy)

# Functions ---------------------------------------------------------------
logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

# Parameters for simulations ----------------------------------------------
# Number of data sets simulated for each scenario
n.sims <- 50
# Set seed to generate same values
set.seed(99382)
# Random seeds for each data set
my.seeds <- sample(1:100000, n.sims, replace = FALSE)
# Simulation scenarios
J.vals <- c(100, 200, 400, 800, 1600)
n.time.vals <- c(5, 10, 15)
p.svc.vals <- c(2, 3, 4)
# Total number of simulation scenarios
n.scenarios <- length(J.vals) * length(n.time.vals) * length(p.svc.vals)
scenario.vals <- expand.grid(J = J.vals, n.time = n.time.vals, p.svc = p.svc.vals)
# Constant parameters across all scenarios
alpha <- c(logit(0.4), -0.5)
p.RE <- list()
psi.RE <- list()
# Spatial parameters ------------------
sp <- TRUE
cov.model <- 'exponential'

# Simulation setup --------------------------------------------------------
psi.true <- array(NA, dim = c(max(J.vals), max(n.time.vals), n.sims, n.scenarios))
w.true <- array(NA, dim = c(max(J.vals), max(p.svc.vals), n.sims, n.scenarios))
psi.mean.samples <- array(NA, dim = c(max(J.vals), max(n.time.vals), n.sims, n.scenarios))
psi.low.samples <- array(NA, dim = c(max(J.vals), max(n.time.vals), n.sims, n.scenarios))
psi.high.samples <- array(NA, dim = c(max(J.vals), max(n.time.vals), n.sims, n.scenarios))
w.mean.samples <- array(NA, dim = c(max(J.vals), max(p.svc.vals), n.sims, n.scenarios))
w.low.samples <- array(NA, dim = c(max(J.vals), max(p.svc.vals), n.sims, n.scenarios))
w.high.samples <- array(NA, dim = c(max(J.vals), max(p.svc.vals), n.sims, n.scenarios))

# MCMC Info ---------------------------
n.samples <- 15000
batch.length <- 25
n.batch <- n.samples / batch.length
n.burn <- 10000
n.thin <- 5
n.chains <- 3
accept.rate <- 0.43

# Simulate Data -----------------------------------------------------------
for (j in 1:n.sims) {
  print(paste("Currently on simulation set ", j, " out of ", n.sims, sep = ''))
  set.seed(my.seeds[j])
  for (i in 1:n.scenarios) {
    print(paste("Currently on scenario ", i, " out of ", n.scenarios, sep = ''))
    J.curr <- scenario.vals$J[i]
    n.time.curr <- scenario.vals$n.time[i]
    p.svc.curr <- scenario.vals$p.svc[i]
    # Note: this is hardcoded.
    if (J.curr == 100) {
      J.x <- 10
      J.y <- 10 
      phi.tune <- 1.25
    } else if (J.curr == 200) {
      J.x <- 20
      J.y <- 10 
      phi.tune <- 1
    } else if (J.curr == 400) {
      J.x <- 20
      J.y <- 20
      phi.tune <- 0.75
    } else if (J.curr == 800) {
      J.x <- 40
      J.y <- 20
      phi.tune <- 0.5
    } else if (J.curr == 1600) {
      J.x <- 40
      J.y <- 40
      phi.tune <- 0.35
    }
    n.rep <- matrix(NA, J.curr, n.time.curr)
    for (k in 1:J.curr) {
      n.rep[k, 1:n.time.curr] <- rep(4, n.time.curr)
    }
    beta <- c(logit(0.4), 0, 0, 0, 0)
    # Note: this is hardcoded
    if (p.svc.curr == 1) {
      svc.cols <- 1
    } else if (p.svc.curr == 2) {
      svc.cols <- 1:2
    } else if (p.svc.curr == 3) {
      svc.cols <- 1:3
    } else if (p.svc.curr == 4) {
      svc.cols <- 1:4
    } else if (p.svc.curr == 5) {
      svc.cols <- 1:5
    }
    sigma.sq <- rep(0.75, p.svc.curr) 
    # sigma.sq <- runif(p.svc.curr, 0.2, 1.5)
    # phi <- runif(p.svc.curr, 3 / 1, 3 / 0.1)
    phi <- rep(3 / 0.5, p.svc.curr)
    trend <- FALSE 
    sp.only <- 0
    rho <- 0.7
    sigma.sq.t <- 1.3
    ar1 <- TRUE
    dat <- simTOcc(J.x = J.x, J.y = J.y, n.time = rep(n.time.curr, J.curr), beta = beta, alpha = alpha, 
		   n.rep = n.rep, sp = sp, svc.cols = svc.cols, cov.model = cov.model, 
		   trend = trend, sp.only = sp.only, sigma.sq = sigma.sq, phi = phi, 
		   ar1 = ar1, rho = rho, sigma.sq.t = sigma.sq.t)
    psi.true[1:J.curr, 1:n.time.curr, j, i] <- dat$psi
    w.true[1:J.curr, 1:p.svc.curr, j, i] <- dat$w
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
    occ.covs <- list(int = X[, , 1], 
                     occ.cov.1 = X[, , 2], 
                     occ.cov.2 = X[, , 3], 
                     occ.cov.3 = X[, , 4],
                     occ.cov.4 = X[, , 5]) 
    det.covs <- list(det.cov.1 = X.p[, , , 2])
    data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs, coords = coords)
    # Priors
    prior.list <- list(beta.normal = list(mean = 0, var = 2.72), 
    		       alpha.normal = list(mean = 0, var = 2.72), 
    		       sigma.sq.ig = list(a = 2, b = 1.5), 
		       sigma.sq.t.ig = c(2, 1),
                       phi.unif = list(a = 3 / 1, b = 3 / 0.1)) 
    # Starting values
    z.init <- apply(y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
    inits.list <- list(beta = beta, alpha = alpha, sigma.sq = sigma.sq, phi = phi,
    		       z = z.init, rho = rho, sigma.sq.t = sigma.sq.t)
    # Tuning
    tuning.list <- list(phi = phi.tune, rho = 0.5) 
    # SVC model
    out <- svcTPGOcc(occ.formula = ~ occ.cov.1 + occ.cov.2 + occ.cov.3 + occ.cov.4, 
		     det.formula = ~ det.cov.1,
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
		     ar1 = TRUE,
		     NNGP = TRUE,
		     n.neighbors = 5, 
		     n.report = 25,
		     n.burn = n.burn, 
		     n.thin = n.thin,
		     n.chains = 1)
    psi.mean.samples[1:J.curr, 1:n.time.curr, j, i] <- apply(out$psi.samples, c(2, 3), mean)
    psi.low.samples[1:J.curr, 1:n.time.curr, j, i] <- apply(out$psi.samples, c(2, 3), quantile, 0.025)
    psi.high.samples[1:J.curr, 1:n.time.curr, j, i] <- apply(out$psi.samples, c(2, 3), quantile, 0.975)
    w.mean.samples[1:J.curr, 1:p.svc.curr, j, i] <- t(apply(out$w.samples, c(2, 3), mean))
    w.low.samples[1:J.curr, 1:p.svc.curr, j, i] <- t(apply(out$w.samples, c(2, 3), 
						   quantile, 0.025))
    w.high.samples[1:J.curr, 1:p.svc.curr, j, i] <- t(apply(out$w.samples, c(2, 3), 
						    quantile, 0.975))
  } # i (n.scenarios)
} # j (n.sims)

save(psi.mean.samples, psi.low.samples, psi.high.samples, 
     w.mean.samples, w.low.samples, w.high.samples, 
     psi.true, w.true, scenario.vals, file = 'results/supplemental-sims/sim-2-results.rda')
