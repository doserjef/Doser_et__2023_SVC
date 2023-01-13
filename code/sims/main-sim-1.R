# main-sim-1.R: code to perform the first simulation study to assess the 
#               consequences of failing to account for a spatially-varying 
#               coefficient in an occupancy model for assessing model fit, 
#               out-of-sample performance, and coverage rates of the latent
#               occurrence values.
# Author: Jeffrey W. Doser
rm(list = ls())
library(spOccupancy)

# Parameters for simulations ----------------------------------------------
# Number of data sets simulated for each scenario
n.sims <- 50
# Set seed to generate same values
set.seed(10101)
# Spatial locations
J.x <- 20
J.y <- 20
J <- J.x * J.y
# Replicates
# n.rep <- sample(5, J, replace = TRUE)
n.rep <- rep(5, J)
# Occurrence coefficients -------------
beta <- c(0, 0)
# Detection coefficients --------------
alpha <- c(-0.2, 0.4)
# Spatial parameters ------------------
sp <- TRUE
svc.cols <- c(1, 2)
cov.model <- 'exponential'
# Different spatial variance values
sigma.sq.vals <- c(0.1, 0.5, 1, 2)
phi.vals <- c(3 / 0.1, 3 / 0.5, 3 / 0.8)
# Values for spatially-varying intercept
sigma.sq.int <- 1
phi.int <- 3 / 0.4
# Total number of simulation scenarios
n.scenarios <- length(sigma.sq.vals) * length(phi.vals)
scenario.vals <- expand.grid(sigma.sq = sigma.sq.vals, phi = phi.vals)
# Fitting an SVC, SVI, and non-spatial occupancy model. 
n.models <- 3
# Random seeds for each data set and simulation scenario
my.seeds <- sample(1:100000, n.sims * n.scenarios, replace = FALSE)

# Simulation setup --------------------------------------------------------
psi.true <- array(NA, dim = c(J, n.sims, n.scenarios))
psi.mean.samples <- array(NA, dim = c(J, n.sims, n.scenarios, n.models))
psi.low.samples <- array(NA, dim = c(J, n.sims, n.scenarios, n.models))
psi.high.samples <- array(NA, dim = c(J, n.sims, n.scenarios, n.models))
waic.vals <- array(NA, dim = c(n.sims, n.scenarios, n.models))
deviance.vals <- array(NA, dim = c(n.sims, n.scenarios, n.models))

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
    dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha, 
		  sp = sp, svc.cols = svc.cols, cov.model = cov.model, 
		  sigma.sq = c(sigma.sq.int, scenario.vals$sigma.sq[i]), 
		  phi = c(phi.int, scenario.vals$phi[i]))
    psi.true[, j, i] <- dat$psi 
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
    tuning.list <- list(phi = 1) 
    # SVC model
    out <- svcPGOcc(occ.formula = ~ occ.cov.1, 
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
		    NNGP = TRUE,
		    n.neighbors = 5, 
		    n.report = 25,
		    n.burn = n.burn, 
		    n.thin = n.thin,
		    n.chains = 1, 
                    k.fold = 4, 
                    k.fold.threads = 4)
    psi.mean.samples[, j, i, 1] <- apply(out$psi.samples, 2, mean)
    psi.low.samples[, j, i, 1] <- apply(out$psi.samples, 2, quantile, 0.025)
    psi.high.samples[, j, i, 1] <- apply(out$psi.samples, 2, quantile, 0.975)
    waic.vals[j, i, 1] <- waicOcc(out)[3]
    deviance.vals[j, i, 1] <- out$k.fold.deviance

    # SVI model 
    prior.list.2 <- prior.list
    prior.list.2$phi.unif <- c(3 / 1, 3 / 0.05)
    prior.list.2$sigma.sq.ig <- c(2, 1)
    out.2 <- spPGOcc(occ.formula = ~ occ.cov.1, 
		     det.formula = ~ det.cov.1,
		     data = data.list,
		     n.batch = n.batch,
		     batch.length = batch.length,
		     inits = inits.list,
		     priors = prior.list.2,
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
		     n.chains = 1, 
                     k.fold = 4, 
                     k.fold.threads = 4)
    psi.mean.samples[, j, i, 2] <- apply(out.2$psi.samples, 2, mean)
    psi.low.samples[, j, i, 2] <- apply(out.2$psi.samples, 2, quantile, 0.025)
    psi.high.samples[, j, i, 2] <- apply(out.2$psi.samples, 2, quantile, 0.975)
    waic.vals[j, i, 2] <- waicOcc(out.2)[3]
    deviance.vals[j, i, 2] <- out.2$k.fold.deviance

    # Non-spatial model 
    out.3 <- PGOcc(occ.formula = ~ occ.cov.1, 
		   det.formula = ~ det.cov.1,
		   data = data.list,
		   n.samples = n.batch * batch.length,
		   inits = inits.list,
		   priors = prior.list.2,
		   n.omp.threads = 1, # Change as desired
		   verbose = TRUE,
		   n.report = 1000,
		   n.burn = n.burn, 
		   n.thin = n.thin,
		   n.chains = 1, 
                   k.fold = 4, 
                   k.fold.threads = 4)
    psi.mean.samples[, j, i, 3] <- apply(out.3$psi.samples, 2, mean)
    psi.low.samples[, j, i, 3] <- apply(out.3$psi.samples, 2, quantile, 0.025)
    psi.high.samples[, j, i, 3] <- apply(out.3$psi.samples, 2, quantile, 0.975)
    waic.vals[j, i, 3] <- waicOcc(out.3)[3]
    deviance.vals[j, i, 3] <- out.3$k.fold.deviance
  } # i (n.scenarios)
} # j (n.sims)

save(psi.mean.samples, psi.low.samples, psi.high.samples, waic.vals,
     deviance.vals, scenario.vals, psi.true, 
     file = 'results/sim-1-results.rda')
