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
load('data/tree/all-trees-data.rda')

# Filter to one state of interest -----------------------------------------
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
usa <- usa %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs")
state <- usa %>% filter(ID == 'new york')
coords.sf <- st_as_sf(as.data.frame(data.list$coords),
			   coords = c('x', 'y'),
			   crs = st_crs(usa))
site.indx <- which(sapply(st_within(coords.sf, st_make_valid(state)), length) == 1)
data.list$y <- data.list$y[site.indx]
data.list$coords <- data.list$coords[site.indx, ]
data.list$covs <- data.list$covs[site.indx, ]

# Remove 25% of locations -------------------------------------------------
load(paste0("data/tree/pred-indx-", str_replace_all(curr.state, ' ', '-'), '.rda'))
data.list$y <- data.list$y[-pred.indx]
data.list$coords <- data.list$coords[-pred.indx, ]
data.list$covs <- data.list$covs[-pred.indx, ]
data.list$z <- data.list$z[-pred.indx]

# TODO: transform the y's
data.list$y <- log(data.list$y)

load("results/tree/stage-2-biomasss-new-york-SVI-samples.rda")


y.rep.means <- apply(out$y.rep.samples, 2, mean)


curr.indx <- which(data.list$y != 0)
dist.fia <- dist(data.list$coords[curr.indx, ])
mean.dist <- mean(dist.fia)
min.dist <- 5
max.dist <- max(dist.fia)
prior.list <- list(beta.normal = list(mean = 0, var = 100),
		   phi.unif = list(a = 3 / max.dist, b = 3 / min.dist),
                   sigma.sq = list(a = 2, b = 2),
                   tau.sq = c(2, 1))
inits.list <- list(beta = 0, sigma.sq = 5, tau.sq = 0.2, phi = 3 / mean.dist)
tuning.list <- list(phi = c(0.3))

# Compare with spNNGP
# library(spNNGP)
# starting <- list(phi = 3 / mean.dist, sigma.sq = 5, tau.sq = 0.2)
# tuning <- list(phi = 0.3, sigma.sq = 0.5)
# priors <- list(phi.Unif = c(3 / max.dist, 3 / min.dist), sigma.sq.IG = c(2, 2), 
# 	       tau.sq.IG = c(2, 1))
# out.spNNGP <- spNNGP(data.list$y ~ 1, coords = data.list$coords, starting = starting, 
# 		     method = 'latent', n.neighbors = 5, tuning = tuning, priors = priors, 
# 		     cov.model = 'exponential', n.samples = 12500)
# y.quants <- fitted(out.spNNGP)$y.rep.quants
