# tmax-data-prep.R: this script extracts 30 year climate normals for 
#                   maximum temperature at the BBS locations used in the model. 
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(sf)
library(sp)
library(stars)
library(AOI)
library(climateR)
library(raster)
library(rasterVis)

# Load data and coordinates -----------------------------------------------
load("data/case-study-1/spOcc-bbs-data.rda")

# Get prediction coordinates ----------------------------------------------
my.proj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"
coords.sf <- st_as_sf(data.frame(data.list$coords),
		      coords = c("X", "Y"),
		      crs = my.proj)

coords.lat.long <- coords.sf %>%
  st_transform(crs = '+proj=longlat +datum=WGS84')

# Get climate normals from TerraClim --------------------------------------
# Time period to extract
period <- "19812010"
J <- nrow(coords.sf)
# Variables to extract
tmax <- rep(NA, J)
tmp <- getTerraClimNormals(coords.lat.long, 'tmax', period, 1:12)[[1]]
plt_dat <- extract(tmp, coords.lat.long)
tmax <- apply(plt_dat[,2:ncol(plt_dat)], 1, max)

save(tmax, file = 'data/case-study-1/tmax-data.rda')
