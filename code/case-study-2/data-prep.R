# data-prep.R: this script takes the raw BBS data and preps it for analysis. Note
#              that this script does not generate the covariate data. This script
#              will not run, as the raw BBS files are not included here.
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
# This code extracts the BBS data directly from BBS downloaded data.
bbs.dat <- list.files(path = "data/BBS/States/", full.names = TRUE) %>%
  lapply(read.csv) %>%
  bind_rows()
# Get associated route data
route.dat <- read.csv("data/BBS/routes.csv")
# Note that route id is nested within state.
# Join BBS data with route data
bbs.dat <- left_join(bbs.dat, route.dat, by = c('Route', 'CountryNum', 'StateNum'))
# Grab BBS data from 1990 and later
bbs.dat <- bbs.dat %>%
  filter(Year >= 1970, RPID == 101)
# Select columns of interest
bbs.new <- bbs.dat %>%
  dplyr::select(RouteDataID, RouteName, Latitude, Longitude, Year, AOU, starts_with("Count"))
rm(bbs.dat)
gc()
# Fill in Missing NA values
bbs.new <- bbs.new %>%
  complete(AOU, nesting(RouteDataID, RouteName, Latitude, Longitude, Year))
# Filter to only get species of interest
# Grasshopper sparrow: GRSP, 5460
# Vesper sparrow: VESP, 5400
bbs.dat <- bbs.new %>%
  filter(AOU %in% c(5460, 5400))
rm(bbs.new)
gc()
# Replace NAs with 0s for all columns at once.
bbs.dat <- bbs.dat %>%
  replace(is.na(.), 0)
# Add in four letter alpha-numeric code
bbs.dat$alpha.code <- ifelse(bbs.dat$AOU == 5460, 'GRSP', 'VESP')

# Load in the above BBS data ----------------------------------------------
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
  filter(AOU == 5460)
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
  filter(AOU == 5460)
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
# Takes 5 or 10 min or so.
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

det.covs <- list(day = day, 
		 obs = obs, 
		 obs.first.year = obs.first.year, 
		 year.det = year.det)

# Save the data -----------------------------------------------------------
dimnames(y)[[1]] <- as.character(sp.names)
save(coords, y, det.covs, routes.final, file = "data/case-study-2/bbs-data-y-det-covs.rda")
