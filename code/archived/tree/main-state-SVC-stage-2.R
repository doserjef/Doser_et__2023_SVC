# main-uni-gaussian: fits a Gaussian univariate SVC model analogous to spSVC 
rm(list = ls())
library(tidyverse)
library(sf)
library(spAbundance)

# Get info from command line ----------------------------------------------
# This code is to extract the current species name from the command line
# to easily run the script for different species
args <- commandArgs(trailingOnly = TRUE)
# Current species
curr.state <- args[1]
# Alternatively, if not running the script from the command line:
curr.state <- 'new york'
if(length(args) == 0) base::stop('Need to tell spAbundance the current state')
print(curr.state)

# Load data ---------------------------------------------------------------
load('data/tree/stage-2-data.rda')

# Filter to one state of interest -----------------------------------------
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
usa <- usa %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs")
state <- usa %>% filter(ID == 'new york')
coords.sf <- st_as_sf(as.data.frame(data.list.2$coords),
			   coords = c('x', 'y'),
			   crs = st_crs(usa))
site.indx <- which(sapply(st_within(coords.sf, st_make_valid(state)), length) == 1)
data.list.2$y <- data.list.2$y[site.indx]
data.list.2$coords <- data.list.2$coords[site.indx, ]
data.list.2$covs <- data.list.2$covs[site.indx, ]
data.list.2$z <- data.list.2$z[site.indx]

# Remove 25% of locations -------------------------------------------------
load(paste0("data/tree/pred-indx-", str_replace_all(curr.state, ' ', '-'), '.rda'))
data.list.2$y <- data.list.2$y[-pred.indx]
data.list.2$coords <- data.list.2$coords[-pred.indx, ]
data.list.2$covs <- data.list.2$covs[-pred.indx, ]
data.list.2$z <- data.list.2$z[-pred.indx]

# Prep the model ----------------------------------------------------------
# Priors
curr.indx <- which(data.list.2$y != 0)
dist.fia <- dist(data.list.2$coords[curr.indx, ])
mean.dist <- mean(dist.fia)
min.dist <- 5
max.dist <- max(dist.fia)
prior.list <- list(beta.normal = list(mean = 0, var = 100),
		   phi.unif = list(a = 3 / max.dist, b = 3 / min.dist), 
                   sigma.sq.ig = list(a = 2, b = 2), 
                   tau.sq = c(2, 1))
inits.list <- list(beta = 0, sigma.sq = 5, tau.sq = 0.2, phi = 3 / mean.dist)
tuning.list <- list(phi = c(0.3))

# Biggest
n.batch <- 2000
batch.length <- 25
n.burn <- 30000
n.thin <- 20
n.chains <- 1

# Big
# n.batch <- 500
# batch.length <- 25
# n.burn <- 7500
# n.thin <- 5

# n.batch <- 50
# batch.length <- 25
# n.burn <- 500
# n.thin <- 1

data.list.2$covs$vpd <- data.list.2$covs$vpd / sd(data.list.2$covs$vpd)
data.list.2$covs$ppt <- data.list.2$covs$ppt / sd(data.list.2$covs$ppt)

out <- svcAbund(formula = ~ vpd + ppt + scale(elev) + I(scale(elev)^2), 
		  data = data.list.2, priors = prior.list, 
		  inits = inits.list, tuning = tuning.list, svc.cols = c(1, 2, 3),
	          n.neighbors = 5, cov.model = 'exponential', NNGP = TRUE,
	          n.batch = n.batch, batch.length = batch.length, family = 'Gaussian-hurdle',
	          n.burn = n.burn, accept.rate = 0.43, n.thin = n.thin, 
	          n.chains = n.chains, n.report = 1, n.omp.threads = 5)

save(out, file = paste('results/tree/stage-2-sugar-maple-', str_replace(curr.state, ' ', '-'), 
		       '-SVC-samples.rda', sep = ''))

