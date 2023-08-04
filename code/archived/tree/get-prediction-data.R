# get-prediction-data.R: this script extracts the grid for prediction across
#                        Vermont.
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(sf)
library(sp)
library(FedData)
library(stars)

# Get prediction coordinates ----------------------------------------------
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
usa <- usa %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs")
vermont <- usa %>%
  dplyr::filter(ID == 'michigan')

# Grid the area for predictions. 
# The dx x dy indicates the resolution in terms of km
grid.pred <- st_as_stars(st_bbox(vermont), dx = 1, dy = 1)
# Convert to data frame
coords.pred <- as.data.frame(grid.pred, center = TRUE)
# Convert coordinates to an sf object
coords.pred.sf <- st_as_sf(coords.pred, 
			   coords = c('x', 'y'), 
			   crs = st_crs(vermont))

# Intersect with region of interest
coords.pred.sf <- st_intersection(coords.pred.sf, st_make_valid(vermont))
coords.0 <- as.data.frame(st_coordinates(coords.pred.sf))

# Get percent forest data -------------------------------------------------
# TODO: there's a problem here that's causing an error with the calculation of forest
#       cover in certain locations... Need to figure out what the deal is with that.
# Create matrix to store all values for grassland and cropland variables
J <- nrow(coords.pred.sf)
# Loop through all the sites
vals <- split(1:J, ceiling(seq_along(1:J)/100))
# Function to get proportion forest cover
props <- function(a, na.rm = TRUE) {
  my.sum <- sum(!is.na(a))
  prop.for <- sum(a %in% c(41, 42, 43), na.rm = na.rm) / my.sum
  return(prop.for)
}
forest.0 <- rep(NA, J)
for (j in 1:length(vals)) {
  print(paste0("Currently on ", j, " out of ", length(vals)))
  coords.curr <- coords.pred.sf[vals[[j]], ]
  nlcd.dat <- tryCatch({get_nlcd(template = coords.curr, label = paste0('biomass-', j), year = 2019)}, error = function(e) NULL)
  # Calculating forest cover in 1km.
  if (!is.null(nlcd.dat)) {
    forest.0[vals[[j]]] <- raster::extract(nlcd.dat, coords.curr, buffer = 1000, fun = props)
  }
}

# Save prediction coordinates
save(coords.0, file = "data/pred-coordinates.rda")
