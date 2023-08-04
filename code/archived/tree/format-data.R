rm(list = ls())

load("data/tree/vermont-bio-coords.rda")
load('data/tree/tree-canopy-cover.rda')

data.list <- list(y = y[!is.na(tcc)], 
		       covs = data.frame(tcc = tcc[!is.na(tcc)]),
		       coords = coords[!is.na(tcc), ])

save(data.list, file = 'data/tree/spAbundance-data.rda')
