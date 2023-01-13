# birdlife-species-ranges.R: this script calculates the species ranges
#                            (including a 50km buffer) for use when 
#                            defining the area of modeling for each of the 51
#                            forest bird species in the case study. 
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(coda)
library(sf)
library(lwgeom)

# Read in species name data for eastern forest bird community -------------
load("data/spOcc-bbs-data.rda")
sp.names <- dimnames(data.list$y)[[1]]
comm.group.dat <- read.csv("data/bird-species-table-bateman.csv")
scientific.names <- comm.group.dat %>%
  filter(Code %in% sp.names) %>%
  select(Code, name = Scientific.Name) %>%
  mutate(name = tolower(name)) %>%
  unique()

# Three species with differing classifications that you need to adjust
# Nashville Warbler
scientific.names$name[which(scientific.names$Code == 'NAWA')] <- "leiothlypis ruficapilla"
# Pileaated Woodpecker
scientific.names$name[which(scientific.names$Code == 'PIWO')] <- "hylatomus pileatus"
# Red-Cockaded Woodpecker 
scientific.names$name[which(scientific.names$Code == 'RCWO')] <- "leuconotopicus borealis"

# Read in extracted bird-life data ----------------------------------------
# Object read in is called bird.life
load("data/bird-life-data.rda")

# Only extract data for the eastern US
coords.sf <- st_as_sf(as.data.frame(data.list$coords),
		      coords = c("X", "Y"),
		      crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs")

# Map of US states
usa <- st_as_sf(maps::map("usa", fill = TRUE, plot = FALSE))
# Restrict to east of the 100th meridian
usa.bbox <- st_bbox(usa)
usa.bbox[1] <- -100
usa.bbox <- as.vector(usa.bbox)
sf_use_s2(FALSE)
east.us <- st_crop(st_make_valid(usa), xmin = usa.bbox[1], ymin = usa.bbox[2],
                   xmax = usa.bbox[3], ymax = usa.bbox[4])
east.us <- east.us %>%
  st_transform(st_crs(coords.sf))

# Extract necessary information -------------------------------------------
usa <- st_as_sf(maps::map("usa", fill = TRUE, plot = FALSE))
# Restrict to east of the 100th meridian
usa.bbox <- st_bbox(usa)
usa.bbox[1] <- -100
usa.bbox <- as.vector(usa.bbox)
sf_use_s2(FALSE)
east.us <- st_crop(st_make_valid(usa), xmin = usa.bbox[1], ymin = usa.bbox[2],
                   xmax = usa.bbox[3], ymax = usa.bbox[4])
east.us <- east.us %>%
  st_transform(st_crs(coords.sf))
# sf object for the area within the species range.
within.sf.list <- list()
# sf object for the area outside the species range. 
outside.sf.list <- list()
bird.life.sp <- unique(bird.life$sci_name)
# Percentage of breeding range within the eastern US study region. 
percent.area <- rep(NA, length(bird.life.sp))
for (i in 1:length(bird.life.sp)) {
  print(i)
  tmp.indx <- which(scientific.names$name == tolower(bird.life.sp[i]))
  indx <- which(sp.names == scientific.names$Code[tmp.indx])
  tmp <- bird.life %>%
    filter(sci_name == bird.life.sp[i])
  # Need a few special filters for some species where the main filter
  # doesn't work
  if (sp.names[indx] == 'RCWO') {
    curr.full <- tmp %>%
      filter(presence == 1)
  } else if (sp.names[indx] == 'BLJA') {
    curr.full <- tmp %>%
      filter(presence == 1)
  } else if (nrow(tmp) > 1) {
    curr.full <- tmp %>%
      filter(seasonal %in% c(1, 2))
  } else {
    curr.full <- tmp
  }
  curr.sp.range <- curr.full %>%
    st_transform(st_crs(coords.sf))
  curr.sp.range <- st_union(st_make_valid(curr.sp.range))
  filtered.range <- st_intersection(curr.sp.range, east.us[1, 1])
  percent.area[indx] <- as.numeric(st_area(filtered.range) / st_area(curr.sp.range))
  filtered.range <- st_buffer(filtered.range, 50)
  filtered.range <- st_intersection(filtered.range, east.us[1, 1])
  diff.range <- st_difference(east.us[1, 1], filtered.range)
  within.sf.list[[indx]] <- filtered.range
  outside.sf.list[[indx]] <- diff.range
}

# To check it all makes sense
# for (i in 1:nrow(scientific.names)) {
#   print(i)
#   plot(within.sf.list[[i]], main = sp.names[i])
#   Sys.sleep(2)
# }

N <- nrow(data.list$y)
J <- ncol(data.list$y)
range.ind <- matrix(0, N, J)
for (i in 1:N) {
  indx <- unlist(c(st_contains(within.sf.list[[i]], coords.sf)))
  range.ind[i, indx] <- 1
}

# Save out processed bird life data
save(percent.area, within.sf.list, outside.sf.list, 
     range.ind, file = "data/bird-life-processed.rda")

