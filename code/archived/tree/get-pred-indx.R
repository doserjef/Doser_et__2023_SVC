rm(list = ls())

# Load data ---------------------------------------------------------------
load('data/tree/spAbundance-data.rda')

# Generate a hold-out set and remove 25% of locations ---------------------
set.seed(19191)
J <- nrow(data.list$coords)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
save(pred.indx, file = 'data/tree/pred-indx.rda')

