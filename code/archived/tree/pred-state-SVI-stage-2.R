rm(list = ls())
library(tidyverse)
library(sf)
library(spOccupancy)
library(spAbundance)

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
load("data/tree/stage-2-data.rda")

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
covs.0 <- data.list.2$covs[pred.indx, ]
coords.0 <- data.list.2$coords[pred.indx, ]
y.0 <- data.list.2$y[pred.indx]
z.0 <- data.list.2$z[pred.indx]
data.list.2$y <- data.list.2$y[-pred.indx]
data.list.2$coords <- data.list.2$coords[-pred.indx, ]
data.list.2$covs <- data.list.2$covs[-pred.indx, ]
data.list.2$z <- data.list.2$z[-pred.indx]

# Get prediction data -----------------------------------------------------
load(paste0("results/tree/stage-2-sugar-maple-", str_replace_all(curr.state, ' ', '-'), 
	    '-SVI-samples.rda'))
X.0 <- matrix(1, length(pred.indx), ncol(out$X))
colnames(X.0) <- c('intercept', 'vpd', 'vpd.2', 'ppt', 'ppt.2', 'elev', 'elev.2')
X.0[, 'vpd'] <- (covs.0[, 'vpd'] - mean(data.list.2$covs[, 'vpd'])) / sd(data.list.2$covs[, 'vpd'])
X.0[, 'ppt'] <- (covs.0[, 'ppt'] - mean(data.list.2$covs[, 'ppt'])) / sd(data.list.2$covs[, 'ppt'])
X.0[, 'elev'] <- (covs.0[, 'elev'] - mean(data.list.2$covs[, 'elev'])) / sd(data.list.2$covs[, 'elev'])
X.0[, 'vpd.2'] <- X.0[, 'vpd']^2
X.0[, 'ppt.2'] <- X.0[, 'ppt']^2
X.0[, 'elev.2'] <- X.0[, 'elev']^2

# Load predicted first stage samples
# TODO: testing
# load(paste0('results/tree/pred-stage-1-sugar-maple-', 
# 		                str_replace_all(curr.state, ' ', '-'), 
# 		                '-samples.rda'))
z.0.samples <- matrix(z.0, nrow(out$beta.samples), length(y.0), byrow = TRUE)
out.pred <- predict(out, X.0, coords.0, z.0.samples = z.0.samples)
save(out.pred, file = paste0('results/tree/pred-stage-2-sugar-maple-', 
		                str_replace_all(curr.state, ' ', '-'), 
		                '-SVI-samples.rda'))

