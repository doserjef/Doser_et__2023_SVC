# get-auc.R: script that extracts the mean AUC value for each species from 
#            all four candidate models in the eastern forest bird case study.
# Author: Jeffrey W. Doser
rm(list = ls())

# Read in the data in spOccupancy format
load("data/case-study-1/final-spOccupancy-data.rda")
# Species names
sp.names <- dimnames(data.list$y)[[1]]
# Number of species
N <- length(sp.names)

# Data frame to hold mean AUC value for each species.
auc.mean.df <- data.frame(sp = sp.names, 
			  constant = rep(NA, N),
                          tmax = rep(NA, N),
                          bcr = rep(NA, N),
                          svc = rep(NA, N))
# Extract the AUC value for each species from the posterior samples of the 
# AUC values.
for (i in 1:length(sp.names)) {
  print(paste0("Currently on species ", i, " out of ", length(sp.names)))
  load(paste0("results/case-study-1/auc-", sp.names[i], ".rda"))
  auc.mean.df[i, -1] <- apply(auc.df, 2, mean)
}
save(auc.mean.df, file = 'results/case-study-1/auc-results.rda')
