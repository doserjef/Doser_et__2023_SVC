# format-prediction-data.R: this script formats the prediction data for 
#                           use in the grasshopper sparrow case study.
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(stars)
library(AOI)
library(climateR)
library(sf)
library(raster)
library(rasterVis)

# Load the prediction grid ------------------------------------------------
load("data/case-study-2/pred-coordinates.rda")

# Convert coordinates to sf
coords.sf <- st_as_sf(data.frame(coords.0),
		      coords = c("X", "Y"),
		      crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs")

# Convert coordinates to lat long for extracting TMAX
coords.lat.long <- coords.sf %>%
  st_transform(crs = 4326)
# Get 30 year TMAX normal -------------------------------------------------
period <- "19812010"
# Get max temperature ---------------
tmp <- getTerraClimNormals(coords.lat.long, 'tmax', period, 1:12)[[1]]
tmp.2 <- extract(tmp, coords.lat.long)
# Get max temp over years
tmax <- apply(tmp.2[, 2:ncol(tmp.2)], 1, max)

# Get grassland and cropland deviations -----------------------------------
J <- nrow(coords.0)
years <- 1970:2019
n.years <- length(years)
grass <- matrix(0, J, n.years)
cropland <- matrix(0, J, n.years)
for (t in 1:n.years) {
  load(paste("data/case-study-2/yearly-lulc/lulc-pred-vals-", years[t], sep = ''))
  grass[, t] <- grass.cov
  cropland[, t] <- crop.cov
}
grassland.mean <- apply(grass, 1, mean)
cropland.mean <- apply(cropland, 1, mean)
grass.dev <- apply(grass, 2, function(a) a - grassland.mean)
cropland.dev <- apply(cropland, 2, function(a) a - cropland.mean)

# Save covariates ---------------------------------------------------------
occ.pred.covs <- list(grass.dev = grass.dev,
		      crop.dev = cropland.dev,
		      crop.mean = cropland.mean,
		      grass.mean = grassland.mean,
		      tmax = tmax,
		      grass = grass,
		      crop = cropland)

save(occ.pred.covs, coords.0, file = 'data/case-study-2/GRSP-pred-data.rda')
