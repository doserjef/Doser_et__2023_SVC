rm(list = ls())
library(tidyverse)
library(sf)
library(MODIStsp)

# Get boundary of region of interest (aka Vermont) ------------------------
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
usa <- usa %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs")
vermont <- usa %>%
  dplyr::filter(ID == 'vermont')

