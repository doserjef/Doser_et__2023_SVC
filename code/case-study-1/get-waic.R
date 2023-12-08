# get-waic.R: script to extract the WAIC from each of the four candidate
#             models for each species. This avoids having to read in the 
#             spOccupancy model object every time to calculate WAIC.
# Author: Jeffrey W. Doser
rm(list = ls())
library(spOccupancy)

# Read in the data --------------------------------------------------------
load("data/case-study-1/final-spOccupancy-data.rda")
# Loads the data.list object, which is formatted for spOccupancy.
# Species info
sp.names <- dimnames(data.list$y)[[1]]
N <- length(sp.names)

# Load full model results -------------------------------------------------
# Data frame to hold WAIC results
waic.df <- data.frame(constant = rep(NA, N),
		      tmax = rep(NA, N),
		      bcr = rep(NA, N),
		      svc = rep(NA, N))
for (i in 1:length(sp.names)) {
  print(paste0("Currently on species ", i, " out of ", length(sp.names)))
  # Constant
  load(paste0("results/case-study-1/stPGOcc-", sp.names[i], ".rda"))
  waic.df$constant[i] <- waicOcc(out)[3]
  # Max temperature
  load(paste0("results/case-study-1/stPGOcc-tmax-", sp.names[i], ".rda"))
  waic.df$tmax[i] <- waicOcc(out)[3]
  # BCR
  load(paste0("results/case-study-1/stPGOcc-bcr-", sp.names[i], ".rda"))
  waic.df$bcr[i] <- waicOcc(out)[3]
  # SVC
  load(paste0("results/case-study-1/svcTPGOcc-", sp.names[i], ".rda"))
  waic.df$svc[i] <- waicOcc(out)[3]
  rm(out)
  gc()
}
save(waic.df, file = 'results/case-study-1/waic-results.rda')
