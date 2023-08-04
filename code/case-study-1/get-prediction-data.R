# get-prediction-data.R: this script extracts the grid for prediction of
#                        bird occurrence trends across the eastern US.
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(sf)
library(sp)
library(FedData)
library(stars)
library(AOI)
library(climateR)
library(raster)
library(rasterVis)

# Get prediction coordinates ----------------------------------------------
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
# Restrict to east of the 100th meridian
usa.bbox <- st_bbox(usa)
usa.bbox[1] <- -100
usa.bbox <- as.vector(usa.bbox)
sf_use_s2(FALSE)
east.us <- st_crop(st_make_valid(usa), xmin = usa.bbox[1], ymin = usa.bbox[2], 
                   xmax = usa.bbox[3], ymax = usa.bbox[4])
east.us <- east.us %>%
  st_transform(st_crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"))

# Grid the area for predictions. 
# The dx x dy indicates the resolution in terms of km
grid.pred <- st_as_stars(st_bbox(east.us), dx = 20, dy = 20)
# Convert to data frame
coords.pred <- as.data.frame(grid.pred, center = TRUE)
# Convert coordinates to an sf object
coords.pred.sf <- st_as_sf(coords.pred, 
			   coords = c('x', 'y'), 
			   crs = st_crs(east.us))

# Intersect with region of interest
coords.pred.sf <- st_intersection(coords.pred.sf, st_make_valid(east.us))
coords.0 <- as.data.frame(st_coordinates(coords.pred.sf))

bcrs <- st_read(dsn = "data/case-study-1/BCR_Terrestrial/", layer = "BCR_Terrestrial_master")

bcrs.albers <- bcrs %>%
  st_transform(crs = st_crs(coords.pred.sf))

# Join points toBCRs
tmp <- st_join(coords.pred.sf, bcrs.albers, join = st_intersects)

bcr.pred.factor <- tmp$BCR

# Get climate normals from TerraClim --------------------------------------
# Time period to extract
period <- "19812010"
coords.pred.lat.long <- coords.pred.sf %>%
  st_transform(4326)
J <- nrow(coords.pred.lat.long)
# Variables to extract
tmax <- rep(NA, J)
tmp <- getTerraClimNormals(coords.pred.lat.long, 'tmax', period, 1:12)[[1]]
plt_dat <- extract(tmp, coords.pred.lat.long)
tmax.0 <- apply(plt_dat[,2:ncol(plt_dat)], 1, max)

# Remove sites with NAs 
bad.indx <- which(is.na(tmax.0) | is.na(bcr.pred.factor))
coords.0 <- coords.0[-bad.indx, ]
tmax.0 <- tmax.0[-bad.indx]
bcr.pred.factor <- bcr.pred.factor[-bad.indx]

# Save prediction coordinates
save(coords.0, bcr.pred.factor, tmax.0, file = "data/case-study-1/pred-coordinates.rda")
