# format-spOccupancy-data.R: this script formats all the component parts
#                            of an occupancy model in the list format 
#                            required for spOccupancy. 
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(stars)
library(AOI)
library(climateR)
library(sf)
library(raster)
library(rasterVis)

# Get coordinates and y and det.covs --------------------------------------
load("data/case-study-2/bbs-data-y-det-covs.rda")
# Load in r
routes <- st_read(dsn = "data/BBS/route-shapefile/")
coords.sf <- st_as_sf(as.data.frame(coords),
		      coords = c('X', 'Y'), 
		      crs = st_crs(routes))

coords.sf.albers <- coords.sf %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs")
# Get coordinates in Albers Equal Area (in kilometers)
coords.albers <- st_coordinates(coords.sf.albers)
coords.sf.lat.long <- coords.sf %>%
  st_transform(crs = 4326)

# Get 30 year TMAX normal -------------------------------------------------
period <- "19812010"
# Get max temperature ---------------
tmp <- getTerraClimNormals(coords.sf.lat.long, 'tmax', period, 1:12)[[1]]
tmp.2 <- extract(tmp, coords.sf.lat.long)
# Get max temp over years
tmax <- apply(tmp.2[, 2:ncol(tmp.2)], 1, max)
# Impute one value with it's neighbors that's not picked up by TerraClim
# NOTE: Hardcoded
tmax[2521] <- 30.9

# Get grassland and cropland deviations -----------------------------------
J <- nrow(coords)
years <- 1970:2019
n.years <- length(years)
grass <- matrix(0, J, n.years)
cropland <- matrix(0, J, n.years)
for (t in 1:n.years) {
  load(paste("data/case-study-2/yearly-lulc/lulc-vals-", years[t], sep = ''))
  grass[, t] <- grass.cov
  cropland[, t] <- crop.cov
}
grassland.mean <- apply(grass, 1, mean)
cropland.mean <- apply(cropland, 1, mean)
grass.dev <- apply(grass, 2, function(a) a - grassland.mean)
cropland.dev <- apply(cropland, 2, function(a) a - cropland.mean)

# Format occupancy covariates for spOccupancy -----------------------------
occ.covs <- list(grass.dev = grass.dev, 
		 crop.dev = cropland.dev, 
		 crop.mean = cropland.mean,
		 grass.mean = grassland.mean, 
		 tmax = tmax,
		 grass = grass,
		 crop = cropland)
# Only use sites that are surveyed at least 5 years
keep.indx <- which(apply(y[1, , , 1], 1, function(a) sum(!is.na(a)) >=5))

# Get individual data set for each species and filter to species range ----
# Note these ranges come from BirdLife, which require a data agreement, and 
# so it is not available on GitHub.
load("data/case-study-2/species-BirdLife-ranges.rda")
# GRSP --------------------------------
indx.tmp <- unlist(c(st_contains(GRSP.range, coords.sf.albers)))
indx <- indx.tmp[which(indx.tmp %in% keep.indx)]
data.GRSP <- list(y = y, 
		  occ.covs = occ.covs, 
		  det.covs = det.covs, 
		  coords = coords.albers)
sp.names <- dimnames(y)[[1]]
data.GRSP$y <- y[which(sp.names == 'GRSP'), indx, , ]
data.GRSP$occ.covs$grass.dev <- data.GRSP$occ.covs$grass.dev[indx, ]
data.GRSP$occ.covs$crop.dev <- data.GRSP$occ.covs$crop.dev[indx, ]
data.GRSP$occ.covs$grass.mean <- data.GRSP$occ.covs$grass.mean[indx]
data.GRSP$occ.covs$crop.mean <- data.GRSP$occ.covs$crop.mean[indx]
data.GRSP$occ.covs$tmax <- data.GRSP$occ.covs$tmax[indx]
data.GRSP$occ.covs$grass <- data.GRSP$occ.covs$grass[indx, ]
data.GRSP$occ.covs$crop <- data.GRSP$occ.covs$crop[indx, ]
data.GRSP$det.covs$day <- data.GRSP$det.covs$day[indx, ]
data.GRSP$det.covs$obs <- data.GRSP$det.covs$obs[indx, ]
data.GRSP$det.covs$obs.first.year <- data.GRSP$det.covs$obs.first.year[indx, ]
data.GRSP$det.covs$year.det <- data.GRSP$det.covs$year.det[indx, ]
data.GRSP$coords <- data.GRSP$coords[indx, ]
data.GRSP$det.covs$replicate <- array(NA, dim = dim(data.GRSP$y))
data.GRSP$det.covs$replicate[, , 1] <- 1
data.GRSP$det.covs$replicate[, , 2] <- 2
data.GRSP$det.covs$replicate[, , 3] <- 3
data.GRSP$det.covs$replicate[, , 4] <- 4
data.GRSP$det.covs$replicate[, , 5] <- 5

# Save to hard drive ------------------------------------------------------
save(data.GRSP, file = 'data/case-study-2/GRSP-spOcc-data.rda')

