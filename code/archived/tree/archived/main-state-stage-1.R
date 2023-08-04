# main-full-stage-1.R: this script fits a spatially-explicit GLM for sugar
#                      maple across the US.
rm(list = ls())
library(tidyverse)
library(sf)
library(spOccupancy)

# Get info from command line ----------------------------------------------
# This code is to extract the current species name from the command line
# to easily run the script for different species
args <- commandArgs(trailingOnly = TRUE)
# Current species
curr.state <- args[1]
# Alternatively, if not running the script from the command line:
# TODO: 
curr.state <- 'new york'
if(length(args) == 0) base::stop('Need to tell spOccupancy the current state')
print(curr.state)

# Load data ---------------------------------------------------------------
load("data/tree/stage-1-data.rda")

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
data.list$weights <- data.list$weights[site.indx]

# Remove 25% of locations -------------------------------------------------
set.seed(4040)
J <- nrow(data.list$coords)
pred.indx <- sample(1:J, round(J * .25), replace = FALSE)
save(pred.indx, file = paste0('data/tree/pred-indx-', str_replace_all(curr.state, ' ', '-'), 
			      '.rda'))
data.list$y <- data.list$y[-pred.indx]
data.list$coords <- data.list$coords[-pred.indx, ]
data.list$covs <- data.list$covs[-pred.indx, ]
data.list$weights <- data.list$weights[-pred.indx]


# Prep the model ----------------------------------------------------------
# Priors
dist.fia <- dist(data.list$coords)
mean.dist <- mean(dist.fia) # 1783.395
# min.dist <- min(dist.fia) # 0.007832157
min.dist <- 5
max.dist <- max(dist.fia) # 4642.168
# Might need further restrictions on phi later on. 
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   phi.unif = list(a = 3 / max.dist, b = 3 / min.dist), 
                   sigma.sq.unif = list(0, 10))
inits.list <- list(beta = 0, phi = 3 / mean.dist, sigma.sq = 1)
tuning.list <- list(phi = 0.05, sigma.sq = 0.1)
n.neighbors <- 5

cov.model <- "exponential"
n.chains <- 1

# Biggest
n.batch <- 2000
batch.length <- 25
n.burn <- 30000
n.thin <- 20

# Medium
# n.batch <- 200
# batch.length <- 25
# n.burn <- 3000
# n.thin <- 5
# 
# n.batch <- 1
# batch.length <- 25
# n.burn <- 0
# n.thin <- 1

out <- svcPGBinom(formula = ~ scale(tmin) + I(scale(tmin)^2) + 
                              scale(ppt) + I(scale(ppt)^2) + 
			      scale(elev) + I(scale(elev)^2), 
		  data = data.list, priors = prior.list, 
		  inits = inits.list, tuning = tuning.list, svc.cols = c(1),
	          n.neighbors = n.neighbors, cov.model = cov.model, NNGP = TRUE,
	          n.batch = n.batch, batch.length = batch.length, 
	          n.burn = n.burn, accept.rate = 0.43, n.thin = n.thin, 
	          n.chains = n.chains, n.report = 1, n.omp.threads = 5)

save(out, file = paste0('results/tree/stage-1-sugar-maple-', 
		   str_replace_all(curr.state, ' ', '-'), '-samples.rda'))

