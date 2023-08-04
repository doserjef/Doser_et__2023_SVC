# get-prediction-coords.R: this script extracts the grid for prediction across
#                        the US
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(sf)
library(sp)
library(stars)

# Get prediction coordinates ----------------------------------------------
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
# Restrict to east of the 100th meridian
usa.bbox <- st_bbox(usa)
usa.bbox[1] <- -100
usa.bbox <- as.vector(usa.bbox)
sf_use_s2(FALSE)
usa <- usa %>%
  st_transform(st_crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"))
# Predict effects across GRSP range
load("data/case-study-2/species-BirdLife-ranges.rda")
GRSP.bbox <- st_bbox(GRSP.range)


# Grid the area for predictions. 
# The dx x dy indicates the resolution in terms of km
grid.pred <- st_as_stars(GRSP.bbox, dx = 10, dy = 10)
# Convert to data frame
coords.pred <- as.data.frame(grid.pred, center = TRUE)
# Convert coordinates to an sf object
coords.pred.sf <- st_as_sf(coords.pred, 
			   coords = c('x', 'y'), 
			   crs = st_crs(GRSP.bbox))

# Intersect with region of interest
coords.pred.sf <- st_intersection(coords.pred.sf, st_make_valid(GRSP.range))
coords.0 <- as.data.frame(st_coordinates(coords.pred.sf))

# Save prediction coordinates
save(coords.0, file = "data/case-study-2/pred-coordinates.rda")


