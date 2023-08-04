# pred-svc.R: script to generate SVC predictions for the SVC model for the 
#             Grasshopper Sparrow case study.
# Author: Jeffrey W. Doser
rm(list = ls())
library(spOccupancy)

# Load the data used to fit the model -------------------------------------
load('data/case-study-2/GRSP-spOcc-data.rda')

# Load the model fit object -----------------------------------------------
load('results/case-study-2/svc-GRSP-model-results.rda')

# Load the prediction data ------------------------------------------------
load('data/case-study-2/GRSP-pred-data.rda')
# Standardize/scale variables by those used to fit the model
J.0 <- nrow(coords.0)
n.time <- ncol(out$X)
X.0 <- array(1, dim = c(J.0, ncol(out$X), dim(out$X)[3]))
dimnames(X.0)[[3]] <- c('(Intercept)', 'grass.dev', 'crop.dev', 
			'scale(tmax)', 'I(scale(tmax)^2)', 'scale(grass.mean)', 
			'scale(crop.mean)')
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
}
vals <- split(1:J.0, ceiling(seq_along(1:J.0) / 500))
n.samples <- out$n.post * out$n.chains
p.svc <- length(out$svc.cols)
w.samples <- array(NA, dim = c(n.samples, p.svc, J.0))
for (j in 1:length(vals)) {
  print(paste("Currently on set ", j, " out of ", length(vals), sep = ''))
  curr.indx <- vals[[j]]
  out.pred <- predict(out, X.0[curr.indx, 1:2, ], coords.0[curr.indx, ],
        	      t.cols = 1:2, n.omp.threads = 10, verbose = FALSE)
  w.samples[, , curr.indx] <- out.pred$w.0.samples
}
save(w.samples, file = 'results/case-study-2/svc-pred-results.rda')
