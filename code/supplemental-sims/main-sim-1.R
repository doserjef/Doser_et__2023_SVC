# main-sim-1-occ.R: code to perform the first supplemental simulation study to assess how 
#                   sample size influences the bias and precision of 
#                   estimates from SVC models.
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
set.seed(1996)
# Random seeds for each data set
my.seeds <- sample(1:100000, n.sims, replace = FALSE)
# Simulation scenarios
J.vals <- c(200, 400, 800, 1200, 1600, 2000, 4000, 6000)
p.svc.vals <- c(2, 3, 4)
# Total number of simulation scenarios
n.scenarios <- length(J.vals) * length(p.svc.vals)
scenario.vals <- expand.grid(J = J.vals, p.svc = p.svc.vals)
# Constant parameters across all scenarios
alpha <- c(logit(0.4), -0.5)
p.RE <- list()
psi.RE <- list()
# Spatial parameters ------------------
sp <- TRUE
cov.model <- 'exponential'

# Simulation setup --------------------------------------------------------
psi.true <- array(NA, dim = c(max(J.vals), n.sims, n.scenarios))
w.true <- array(NA, dim = c(max(J.vals), max(p.svc.vals), n.sims, n.scenarios))
psi.mean.samples <- array(NA, dim = c(max(J.vals), n.sims, n.scenarios))
psi.low.samples <- array(NA, dim = c(max(J.vals), n.sims, n.scenarios))
psi.high.samples <- array(NA, dim = c(max(J.vals), n.sims, n.scenarios))
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
    p.svc.curr <- scenario.vals$p.svc[i]
    # Note this is hardcoded.
    if (J.curr == 200) {
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
    } else if (J.curr == 1200) {
      J.x <- 40
      J.y <- 30
      phi.tune <- 0.35
    } else if (J.curr == 1600) {
      J.x <- 40
      J.y <- 40
      phi.tune <- 0.35
    } else if (J.curr == 2000) {
      J.x <- 40
      J.y <- 50
      phi.tune <- 0.2
    } else if (J.curr == 4000) {
      J.x <- 80
      J.y <- 50
      phi.tune <- 0.1
    } else if (J.curr == 6000) {
      J.x <- 100
      J.y <- 60
      phi.tune <- 0.1
    }
    n.rep <- rep(5, J.curr)
    beta <- c(logit(0.4), 0, 0, 0, 0)
    # Note this is hardcoded
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
    dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha, 
		  sp = sp, svc.cols = svc.cols, cov.model = cov.model, 
		  sigma.sq = sigma.sq, phi = phi)
    psi.true[1:J.curr, j, i] <- c(dat$psi)
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
    occ.covs <- X
    colnames(occ.covs) <- c('int', 'occ.cov.1', 'occ.cov.2', 'occ.cov.3', 'occ.cov.4')
    det.covs <- list(det.cov.1 = X.p[, , 2])
    data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs, coords = coords)
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
    tuning.list <- list(phi = phi.tune) 
    # SVC model
    out <- svcPGOcc(occ.formula = ~ occ.cov.1 + occ.cov.2 + occ.cov.3 + occ.cov.4, 
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
		    n.omp.threads = 5, # Change as desired. 
		    verbose = TRUE,
		    NNGP = TRUE,
		    n.neighbors = 5, 
		    n.report = 25,
		    n.burn = n.burn, 
		    n.thin = n.thin,
		    n.chains = 1)
    psi.mean.samples[1:J.curr, j, i] <- apply(out$psi.samples, 2, mean)
    psi.low.samples[1:J.curr, j, i] <- apply(out$psi.samples, 2, quantile, 0.025)
    psi.high.samples[1:J.curr, j, i] <- apply(out$psi.samples, 2, quantile, 0.975)
    w.mean.samples[1:J.curr, 1:p.svc.curr, j, i] <- t(apply(out$w.samples, c(2, 3), mean))
    w.low.samples[1:J.curr, 1:p.svc.curr, j, i] <- t(apply(out$w.samples, c(2, 3), 
						   quantile, 0.025))
    w.high.samples[1:J.curr, 1:p.svc.curr, j, i] <- t(apply(out$w.samples, c(2, 3), 
						    quantile, 0.975))
  } # i (n.scenarios)
} # j (n.sims)

save(psi.mean.samples, psi.low.samples, psi.high.samples, 
     w.mean.samples, w.low.samples, w.high.samples, 
     psi.true, w.true, scenario.vals, file = 'results/supplemental-sims/sim-1-results.rda')
