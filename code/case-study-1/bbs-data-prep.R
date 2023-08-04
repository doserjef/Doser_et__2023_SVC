# bbs-data-prep.R: this script takes the raw BBS data and preps it for analysis.
# Author: Jeffrey W. Doser

rm(list = ls())
library(tidyverse)
library(lubridate)
library(sp)
library(raster)
library(sf)
library(stars)
# For linking bird codes with species info.
library(wildlifeR)

# Read in BBS Data --------------------------------------------------------
# This code extracts the BBS data directly from BBS downloaded data. It is 
# memory intensive, and so it is commented out and I read in the saved object
# bbs.dat below that contains the data for the species of interest across the 
# entire continental US from 1972-present. 
# These are the 10-stop summary data
# bbs.dat <- list.files(path = "data/BBS/States/", full.names = TRUE) %>%
#   lapply(read.csv) %>%
#   bind_rows()
# # Get associated route data
# route.dat <- read.csv("data/BBS/routes.csv")
# # Note that route id is nested within state. 
# # Join BBS data with route data
# bbs.dat <- left_join(bbs.dat, route.dat, by = c('Route', 'CountryNum', 'StateNum'))
# # Grab BBS data from 1970 and later, as well as routes with their starting points 
# # east of the 100th meridian. 
# bbs.dat <- bbs.dat %>%
#   filter(Year >= 1970, RPID == 101, Longitude > -100)
# # Select columns of interest
# bbs.dat <- bbs.dat %>%
#   dplyr::select(RouteDataID, RouteName, Latitude, Longitude, Year, AOU, starts_with("Count"))
# # Fill in implicit zeros. For this situation, gives a value for each species in 
# # each existing combination of RouteDataID, Latitude, Longitude, and Year
# # RouteDataID is a unique identifier for each combo of CountryNum, StateNum, Route, RPID, and Year
# bbs.dat <- bbs.dat %>%
#   complete(AOU, nesting(RouteDataID, RouteName, Latitude, Longitude, Year))
# # Replace NAs with 0s for all columns at once.
# bbs.dat <- bbs.dat %>%
#   replace(is.na(.), 0)
# # Extract community of interior forest obliage birds
# aou.info <- AOU_species_codes %>%
#   mutate(AOU = spp.num)
# bbs.dat <- left_join(bbs.dat, aou.info, by = c("AOU"))
# # Work with community of eastern forest birds
# comm.group.dat <- read.csv("data/bird-species-table-bateman.csv")
# my.sp.code <- comm.group.dat %>%
#   filter(Group %in% c("Eastern.Forests"))
# bbs.dat <- bbs.dat %>%
#   mutate(alpha.code = as.character(alpha.code)) %>%	
#   filter(alpha.code %in% my.sp.code$Code)
# # 66 eastern forest birds
# # Save the BBS data to avoid having to run the above memory intensive code a bunch.
# save(bbs.dat, file = "data/case-study-1/bbs-raw-data.rda")

# Load in the above BBS data ----------------------------------------------
load("data/case-study-1/bbs-raw-data.rda")
# Get associated weather data
weather.dat <- read.csv("data/BBS/weather.csv")

# Join data with Weather data
# Get date in proper format
weather.dat <- weather.dat %>% 
  unite('date', sep = '-', Year, Month, Day, remove = FALSE)
weather.dat$date <- as.Date(weather.dat$date, tz = "America/New_York")
# Get julian date of each survey
weather.dat$julian <- as.numeric(format(weather.dat$date, '%j'))
weather.covs <- weather.dat %>%
  filter(RouteDataID %in% unique(bbs.dat$RouteDataID))

# Get data in format for spOccupancy --------------------------------------
bbs.dat <- left_join(bbs.dat, weather.covs, by = c('RouteDataID', 'Year'))
# Sort data by species, Year, Longitude, Latitude
bbs.dat <- bbs.dat %>%
  arrange(AOU, Year, Longitude, Latitude)
# bbs.dat %>%
#   group_by(RouteDataID) %>%
#   summarize(n.coords = n_distinct(Longitude)) %>%
#   arrange(desc(n.coords))
coords <- unique(bbs.dat[, c('Longitude', 'Latitude')])
# Number of routes
J <- nrow(coords)
# Create a site index for easy linking. 
coords$site.index <- 1:J
# Join the site index with the full data
bbs.dat <- bbs.dat %>%
  left_join(coords, by = c('Longitude', 'Latitude'))
# Just grab from bbs.dat for a single species. I know, not very elegant...
dat.one.sp <- bbs.dat %>%
  filter(AOU == 2280)
route.names <- dat.one.sp$RouteName
site.curr <- dat.one.sp$site.index
n.tmp <- length(route.names)

# Get route locations -----------------------------------------------------
routes <- st_read(dsn = "data/BBS/route-shapefile/")
# Convert to sf
coords.sf <- st_as_sf(coords,
		      coords = c("Longitude", "Latitude"),
		      crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
coords.sf <- coords.sf %>%
  st_transform(st_crs(routes))

ggplot() +
  geom_sf(data = routes) +
  geom_sf(data = coords.sf)

# Find the route nearest to starting location provided in the main BBS data set.
nearest.route <- st_nearest_feature(x = coords.sf, y = routes)

# Extract the routes that are someones nearest route.
routes.sub <- routes[nearest.route, ]

# Get indices for the nearest route to each data point, only keeping those
# that have the same route name as the main data set, as this is an indicator
# for determining if the nearest route is actually the route, or is just the nearest
# route with the full shapefile that is not actually the same route.
good.indices <- rep(FALSE, J)
for (j in 1:J) {
  print(j)
  tmp <- dat.one.sp %>%
    filter(site.index == j)
  curr.names <- unique(tmp$RouteName)
  for (t in 1:length(curr.names)) {
    if (curr.names[t] == routes.sub[j, ]$RTENAME) {
      good.indices[j] <- TRUE
    }
  }
}

# Extract the final routes
routes.final <- routes.sub[good.indices, ]
# Get the centroid of each route for use in spatial models.
route.center <- st_centroid(routes.final)

route.vals <- which(good.indices)
coords <- st_coordinates(route.center)
# Only use data from the sites where you have the full route centroid.
bbs.dat <- bbs.dat %>%
  filter(site.index %in% route.vals)
J <- nrow(coords)

# Detection Covariates ----------------
dat.one.sp <- bbs.dat %>%
  filter(AOU == 2280)
# Get the first year of observation for each observer
obs.first.year <- dat.one.sp %>%
  group_by(ObsN) %>%
  summarize(first.year = first(Year))
# Number of observers
n.obs <- nrow(obs.first.year)
dat.one.sp$obs.first.year <- 0
# Takes a few seconds
for (i in 1:n.obs) {
  tmp <- unlist(c((obs.first.year[i, ])))
  indx <- which((dat.one.sp$ObsN == tmp[1]) & (dat.one.sp$Year == tmp[2]))
  dat.one.sp$obs.first.year[indx] <- 1
}
# Number of species
N <- n_distinct(bbs.dat$AOU)
sp.codes <- unique(bbs.dat$AOU)
sp.names <- AOU_species_codes$alpha.code[which(AOU_species_codes$spp.num %in% sp.codes)]
# Number of years
n.years <- n_distinct(bbs.dat$Year)
unique.years <- unique(bbs.dat$Year)
# Starting coordinates of each route
# Number of stop replicates (here reduced to 5, not using the full 50 stop data)
K <- 5
# Four dimensional array of detection-nondetection data.
y <- array(NA, dim = c(N, J, n.years, K))
# Detection covariates (two-dimensional arrays)
day <- matrix(NA, J, n.years)
tod <- matrix(NA, J, n.years)
obs <- matrix(NA, J, n.years)
obs.first.year <- matrix(NA, J, n.years)
year.det <- matrix(NA, J, n.years)
# Add the first year of observer indicator to the full N species data set. 
bbs.dat$obs.first.year <- rep(dat.one.sp$obs.first.year, times = N)
# Unique site indices for linking
site.ordered.indx <- sort(unique(bbs.dat$site.index))
# Calculate
for (j in 1:J) {
  print(j)
  tmp <- bbs.dat %>%
    filter(site.index == site.ordered.indx[j]) %>%
    arrange(Year)
  tmp.2 <- tmp %>%
    group_by(Year) %>%
    summarize(julian = first(julian), 
              ObsN = first(ObsN), 
              obs.first.year = first(obs.first.year))
  tmp.years <- which(unique.years %in% tmp.2$Year)
  day[j, tmp.years] <- tmp.2$julian 
  obs[j, tmp.years] <- tmp.2$ObsN
  obs.first.year[j, tmp.years] <- tmp.2$obs.first.year
  year.det[j, tmp.years] <- tmp.2$Year
}
# Takes 1 hour or so. 
for (j in 1:J) {
  print(paste("Currently on site ", j, " out of ", J, sep = ''))
  # Get detection-nondetection data for each species. 
  # group_by is necessary to only grab data from one observer for the rare 
  # occassions where there are multiple observers for a given route
  # in a given year. 
  for (i in 1:N) {
    tmp <- bbs.dat %>%
      filter(AOU == sp.codes[i], site.index == site.ordered.indx[j]) %>%
      arrange(Year) %>%
      group_by(Year) %>%
      dplyr::select(Year, starts_with('Count')) %>%
      summarize(across(.cols = everything(), .fns = first))    
    tmp.years <- which(unique.years %in% tmp$Year)
    tmp.2 <- as.matrix(tmp[, c('Count10', 'Count20', 'Count30', 'Count40',
  				'Count50')])
    y[i, j, tmp.years, ] <- tmp.2 
    y[i, j, tmp.years, ] <- ifelse(y[i, j, tmp.years, ] > 0, 1, 0)
  }
}

# Filter to only grab years of interest.
# Last two decades (2000 - 2020)
my.years <- 31:50
y <- y[, , my.years, ]
site.indx <- which(apply(y, 2, sum, na.rm = TRUE) != 0)
y <- y[, site.indx, , ]
det.covs <- list(day = day[site.indx, my.years], 
		 obs = obs[site.indx, my.years], 
		 obs.first.year = obs.first.year[site.indx, my.years], 
		 year.det = year.det[site.indx, my.years])
coords.sf <- st_as_sf(data.frame(coords[site.indx, ]),
		      coords = c("X", "Y"),
		      crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
# Albers equal area across contiguous US.
coords.sf.albers <- coords.sf %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs")
# Get coordinates in Albers Equal Area (in kilometers)
coords.albers <- st_coordinates(coords.sf.albers)
J <- nrow(coords.albers)

# Format data for spOccupancy ---------------------------------------------
years <- t(matrix(rep(2000:2019, J), length(my.years), J))
occ.covs <- list(years = years)
sp.names <- as.character(sp.names)
dimnames(y)[[1]] <- sp.names
data.list <- list(y = y, 
		  occ.covs = occ.covs, 
		  det.covs = det.covs,
		  coords = coords.albers)

# Add in BCRs to the occurrence covariates --------------------------------
bcrs <- st_read(dsn = "data/case-study-1/BCR_Terrestrial/", layer = "BCR_Terrestrial_master")
my.proj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"
coords.sf <- st_as_sf(as.data.frame(data.list$coords),
		      coords = c("X", "Y"),
		      crs = my.proj)

bcrs.albers <- bcrs %>%
  st_transform(crs = my.proj)

# Join points toBCRs
tmp <- st_join(coords.sf, bcrs.albers, join = st_intersects)

bcr.factor <- tmp$BCR
# Fill in values that are messed up
bcr.factor[559] <- 29 
bcr.factor[836] <- 30
bcr.factor[929] <- 30 
bcr.factor[931] <- 14
bcr.factor[1332] <- 26
bcr.factor[1829] <- 28

data.list$occ.covs$BCR <- bcr.factor

# Save the data -----------------------------------------------------------
save(data.list, file = "data/case-study-1/spOcc-bbs-data.rda")

