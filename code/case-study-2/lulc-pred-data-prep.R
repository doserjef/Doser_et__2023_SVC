# lulc-pred-data-prep.R: this script extracts the grassland and cropland 
#                        cover covariates across the continental US from 
#                        USGS EROS for use in predicting effects of land cover
#                        change on Grasshopper Sparrow across the US.
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(stars)
library(sf)

# Get filename to process -------------------------------------------------
# Grab the name from the commandline when running the script. Alternatively, 
# could comment this out and just write the file name manually. 
# Name should be full path name relative to the working directory
file.name <- commandArgs(trailingOnly = TRUE)
# TODO: 
# file.name <- 'CONUS_Backcasting_y1970.tif'
if(length(file.name) == 0) base::stop('Need to give the file name to process')

# Load formatted BBS data -------------------------------------------------
# Reads in object called coords.0
load("data/case-study-2/pred-coordinates.rda")

coords.sf <- st_as_sf(data.frame(coords.0),
		      coords = c("X", "Y"),
		      crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs")

# Calculate covariates within a buffer of the route. 
buffer.radius <- 500

# Land Cover Classes ------------------
# 1: Water
# 2: Developed
# 3: Mechanically Disturbed National Forests
# 4: Mechanically Disturbed Other Public Lands
# 5: Mechanically Disturbed Private
# 6: Mining
# 7: Barren
# 8: Deciduous Forest
# 9: Evergreen Forest
# 10: Mixed Forest
# 11: Grassland
# 12: Shrubland
# 13: Cropland
# 14: Hay/Pasture Land
# 15: Herbaceous Wetland
# 16: Woody Wetland
# 17: Perennial Ice/Snow

# Load rasters ------------------------------------------------------------
lulc.dir <- "~/eros-lulc/"
lulc.curr <- read_stars(paste0(lulc.dir, file.name))
# Convert coordinates to the coordinates of the raster. 
coords.ACEA <- coords.sf %>%
  st_transform(crs = st_crs(lulc.curr))
# Buffer the coordinates
coords.ACEA.buffered <- st_buffer(coords.ACEA, dist = buffer.radius)
sum.cover <- function(a, val) {
  return(sum(a == val))
}
prop.grass <- function(a) {
  mean(a %in% 11, na.rm = TRUE)
}
prop.crop <- function(a) {
  mean(a %in% c(13, 14), na.rm = TRUE)
}
# Create matrix to store all values for 3 forest classes and 1 developed class
J <- nrow(coords.sf)
# Loop through all the sites
vals <- split(1:J, ceiling(seq_along(1:J)/1000))
grass.cov <- rep(0, J)
crop.cov <- rep(0, J)
for (j in 1:length(vals)) {
  print(j)
  coords.curr <- coords.ACEA.buffered[vals[[j]], ]
  tmp <- aggregate(lulc.curr, by = coords.curr, FUN = prop.grass)
  grass.cov[vals[[j]]] <- tmp[[1]]
  tmp <- aggregate(lulc.curr, by = coords.curr, FUN = prop.crop)
  crop.cov[vals[[j]]] <- tmp[[1]]
}

# Save results ------------------------------------------------------------
# Extract year from the file name
curr.year <- parse_number(file.name)
save(grass.cov, crop.cov, file = paste("data/case-study-2/yearly-lulc-500/lulc-pred-vals-", 
				       curr.year, sep = ''))
