# extract-estimates.R: this code extracts the complete spatially-varying 
#                      effect of changes in grassland cover and cropland 
#                      cover on Grasshopper Sparrow occurrence for all of the 
#                      different candidate models. Note the resulting results
#                      files are too large, and so this script will not run
#                      with what is available on GitHub.
# Author: Jeffrey W. Doser
rm(list = ls())
library(spOccupancy)

# Load the data used to fit the model -------------------------------------
load('data/case-study-2/GRSP-spOcc-data.rda')

# Constant model ----------------------------------------------------------
load('results/case-study-2/constant-GRSP-model-results.rda')
beta.names <- colnames(out$beta.samples)
beta.star.grass.dev.samples <- out$beta.samples[, which(beta.names == 'scale(grass.dev)')]
beta.star.crop.dev.samples <- out$beta.samples[, which(beta.names == 'scale(crop.dev)')]
# Save effects to hard drive ----------------------------------------------
save(beta.star.grass.dev.samples, beta.star.crop.dev.samples, 
     file = 'results/case-study-2/constant-GRSP-svc-samples.rda')

# Interaction land cover model --------------------------------------------
load('results/case-study-2/int-lulc-GRSP-model-results.rda')
beta.names <- colnames(out$beta.samples)
beta.grass.dev <- out$beta.samples[, which(beta.names == 'scale(grass.dev)')]
beta.grass.dev.land <- out$beta.samples[, which(beta.names == 'scale(grass.dev):scale(grass.mean)')]
beta.crop.dev <- out$beta.samples[, which(beta.names == 'scale(crop.dev)')]
beta.crop.dev.land <- out$beta.samples[, which(beta.names == 'scale(crop.dev):scale(crop.mean)')]
# Get the interaction covariates to form full covariate effect
load('data/case-study-2/GRSP-pred-data.rda')
# Get average grassland and cropland values scaled by values used to fit model
# Get means and sds of values used to fit model
grass.mean <- mean(data.GRSP$occ.covs$grass.mean)
grass.sd <- sd(data.GRSP$occ.covs$grass.mean)
crop.mean <- mean(data.GRSP$occ.covs$crop.mean)
crop.sd <- sd(data.GRSP$occ.covs$crop.mean)
grass.mean.scaled <- (occ.pred.covs$grass.mean - grass.mean) / grass.sd
crop.mean.scaled <- (occ.pred.covs$crop.mean - crop.mean) / crop.sd
# Form full covariate effect
n.samples <- out$n.post * out$n.chains
J <- nrow(coords.0)
beta.star.grass.dev.samples <- matrix(NA, n.samples, J)
beta.star.crop.dev.samples <- matrix(NA, n.samples, J)
for (i in 1:n.samples) {
  beta.star.grass.dev.samples[i, ] <- beta.grass.dev[i] + beta.grass.dev.land[i] * grass.mean.scaled
  beta.star.crop.dev.samples[i, ] <- beta.crop.dev[i] + beta.crop.dev.land[i] * crop.mean.scaled
}
# Save effects to hard drive ----------------------------------------------
save(beta.star.grass.dev.samples, beta.star.crop.dev.samples, 
     file = 'results/case-study-2/int-lulc-GRSP-svc-samples.rda')

# Interaction tmax model --------------------------------------------
load('results/case-study-2/int-tmax-GRSP-model-results.rda')
beta.names <- colnames(out$beta.samples)
beta.grass.dev <- out$beta.samples[, which(beta.names == 'scale(grass.dev)')]
beta.crop.dev <- out$beta.samples[, which(beta.names == 'scale(crop.dev)')]
beta.grass.dev.tmax <- out$beta.samples[, which(beta.names == 'scale(grass.dev):scale(tmax)')]
beta.crop.dev.tmax <- out$beta.samples[, which(beta.names == 'scale(crop.dev):scale(tmax)')]
# Get the interaction covariates to form full covariate effect
load('data/case-study-2/GRSP-pred-data.rda')
# Get means and sds of values used to fit model
tmax.mean <- mean(data.GRSP$occ.covs$tmax)
tmax.sd <- sd(data.GRSP$occ.covs$tmax)
tmax.scaled <- (occ.pred.covs$tmax - tmax.mean) / tmax.sd
# Form full covariate effect
n.samples <- out$n.post * out$n.chains
J <- nrow(coords.0)
beta.star.grass.dev.samples <- matrix(NA, n.samples, J)
beta.star.crop.dev.samples <- matrix(NA, n.samples, J)
for (i in 1:n.samples) {
  beta.star.grass.dev.samples[i, ] <- beta.grass.dev[i] + beta.grass.dev.tmax[i] * tmax.scaled
  beta.star.crop.dev.samples[i, ] <- beta.crop.dev[i] + beta.crop.dev.tmax[i] * tmax.scaled
}
# Save effects to hard drive ----------------------------------------------
save(beta.star.grass.dev.samples, beta.star.crop.dev.samples, 
     file = 'results/case-study-2/int-tmax-GRSP-svc-samples.rda')


# Full model --------------------------------------------------------------
load('results/case-study-2/full-GRSP-model-results.rda')

# Load the prediction data ------------------------------------------------
load('data/case-study-2/GRSP-pred-data.rda')
# Standardize/scale variables by those used to fit the model
J.0 <- nrow(coords.0)
n.time <- ncol(out$X)
X.0 <- array(1, dim = c(J.0, ncol(out$X), dim(out$X)[3]))
dimnames(X.0)[[3]] <- c('(Intercept)', 'grass.dev', 'crop.dev',
			'scale(tmax)', 'I(scale(tmax)^2)', 'scale(grass.mean)',
			'scale(crop.mean)', 'grass.dev:scale(grass.mean)',
			'crop.dev:scale(crop.mean)', 'grass.dev:scale(tmax)',
			'crop.dev:scale(tmax)')
# Get means and sds of values used to fit model
grass.dev.val <- abs(min(data.GRSP$occ.covs$grass.dev))
crop.dev.val <- abs(min(data.GRSP$occ.covs$crop.dev))
tmax.mean <- mean(data.GRSP$occ.covs$tmax)
tmax.sd <- sd(data.GRSP$occ.covs$tmax)
grass.mean <- mean(data.GRSP$occ.covs$grass.mean)
grass.sd <- sd(data.GRSP$occ.covs$grass.mean)
crop.mean <- mean(data.GRSP$occ.covs$crop.mean)
crop.sd <- sd(data.GRSP$occ.covs$crop.mean)
for (t in 1:n.time) {
  X.0[, t, 'grass.dev'] <- occ.pred.covs$grass.dev[, t] + grass.dev.val
  X.0[, t, 'crop.dev'] <- occ.pred.covs$crop.dev[, t] + crop.dev.val
  X.0[, t, 'scale(tmax)'] <- (occ.pred.covs$tmax - tmax.mean) / tmax.sd
  X.0[, t, 'I(scale(tmax)^2)'] <- X.0[, t, 'scale(tmax)']^2
  X.0[, t, 'scale(grass.mean)'] <- (occ.pred.covs$grass.mean - grass.mean) / grass.sd
  X.0[, t, 'scale(crop.mean)'] <- (occ.pred.covs$crop.mean - crop.mean) / crop.sd
  X.0[, t, 'grass.dev:scale(grass.mean)'] <- X.0[, t, 'grass.dev'] * X.0[, t, 'scale(grass.mean)']
  X.0[, t, 'crop.dev:scale(crop.mean)'] <- X.0[, t, 'crop.dev'] * X.0[, t, 'scale(crop.mean)']
  X.0[, t, 'grass.dev:scale(tmax)'] <- X.0[, t, 'grass.dev'] * X.0[, t, 'scale(tmax)']
  X.0[, t, 'crop.dev:scale(tmax)'] <- X.0[, t, 'crop.dev'] * X.0[, t, 'scale(tmax)']
}

# Load the prediction results ---------------------------------------------
load('results/case-study-2/full-pred-results.rda')

# Get spatially-varying effect of the two variables at each site ----------
beta.names <- colnames(out$beta.samples)
beta.grass.dev <- out$beta.samples[, which(beta.names == 'grass.dev')]
beta.grass.dev.land <- out$beta.samples[, which(beta.names == 'grass.dev:scale(grass.mean)')]
beta.grass.dev.tmax <- out$beta.samples[, which(beta.names == 'grass.dev:scale(tmax)')]
beta.crop.dev <- out$beta.samples[, which(beta.names == 'crop.dev')]
beta.crop.dev.land <- out$beta.samples[, which(beta.names == 'crop.dev:scale(crop.mean)')]
beta.crop.dev.tmax <- out$beta.samples[, which(beta.names == 'crop.dev:scale(tmax)')]
n.samples <- out$n.post * out$n.chains
J <- nrow(coords.0)
beta.star.grass.dev.samples <- matrix(NA, n.samples, J)
beta.star.crop.dev.samples <- matrix(NA, n.samples, J)
# Effect of grassland: grass.dev + grass.dev:scale(grass) * X[, scale(grass.mean)] + 
#                      grass.dev:scale(tmax) * X[, scale(tmax)] + w
# Effect of cropland: crop.dev + crop.dev:scale(crop) * X[, scale(crop.mean)] + 
#                      crop.dev:scale(tmax) * X[, scale(tmax)] + w
for (i in 1:n.samples) {
  print(paste0("Currently on ", i, " out of ", n.samples))
  beta.star.grass.dev.samples[i, ] <- beta.grass.dev[i] +
                                      beta.grass.dev.land[i] * X.0[, 1, 'scale(grass.mean)'] + 
				      beta.grass.dev.tmax[i] * X.0[, 1, 'scale(tmax)'] + 
				      w.samples[i, which(dimnames(out$X.w)[[3]] == 'grass.dev'), ]
  beta.star.crop.dev.samples[i, ] <- beta.crop.dev[i] +
                                      beta.crop.dev.land[i] * X.0[, 1, 'scale(crop.mean)'] + 
				      beta.crop.dev.tmax[i] * X.0[, 1, 'scale(tmax)'] + 
				      w.samples[i, which(dimnames(out$X.w)[[3]] == 'crop.dev'), ]

}
# Save effects to hard drive ----------------------------------------------
save(beta.star.grass.dev.samples, beta.star.crop.dev.samples, 
     file = 'results/case-study-2/full-GRSP-svc-samples.rda')

# SVC model ---------------------------------------------------------------
# Load the model fit object -----------------------------------------------
load('results/case-study-2/svc-GRSP-model-results.rda')

# Load the prediction data ------------------------------------------------
load('data/case-study-2/GRSP-pred-data.rda')
J.0 <- nrow(coords.0)
n.time <- ncol(out$X)

# Load the prediction results ---------------------------------------------
load('results/case-study-2/svc-pred-results.rda')

# Get spatially-varying effect of the two variables at each site ----------
beta.names <- colnames(out$beta.samples)
beta.grass.dev <- out$beta.samples[, which(beta.names == 'grass.dev')]
beta.crop.dev <- out$beta.samples[, which(beta.names == 'crop.dev')]
n.samples <- out$n.post * out$n.chains
J <- nrow(coords.0)
beta.star.grass.dev.samples <- matrix(NA, n.samples, J)
beta.star.crop.dev.samples <- matrix(NA, n.samples, J)
for (i in 1:n.samples) {
  print(paste0("Currently on ", i, " out of ", n.samples))
  beta.star.grass.dev.samples[i, ] <- beta.grass.dev[i] +
				      w.samples[i, which(dimnames(out$X.w)[[3]] == 'grass.dev'), ]
  beta.star.crop.dev.samples[i, ] <- beta.crop.dev[i] +
				     w.samples[i, which(dimnames(out$X.w)[[3]] == 'crop.dev'), ]
}
# Save effects to hard drive ----------------------------------------------
save(beta.star.grass.dev.samples, beta.star.crop.dev.samples, 
     file = 'results/case-study-2/svc-GRSP-svc-samples.rda')
