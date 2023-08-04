# get-bbs-2021-data.R: script to get the 2021 BBS data for use in 
#                      the predictive performance assessment with AUC.
#                      Note the data objects required for this are too large
#                      for GitHub, so this script will not run.
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
# For bird life extraction
library(lwgeom)
# For getting TMAX from TerraClimate
library(AOI)
library(climateR)
library(raster)
library(rasterVis)

# Read in BBS Data --------------------------------------------------------
# This code extracts the BBS data directly from BBS downloaded data. It is
bbs.dat.full <- list.files(path = "data/case-study-1/bbs-2022/States/", full.names = TRUE) %>%
  lapply(read.csv) %>%
  bind_rows()
# Get associated route data
route.dat <- read.csv("data/case-study-1/bbs-2022/routes.csv")
# Note that route id is nested within state.
# Join BBS data with route data
bbs.dat.full <- left_join(bbs.dat.full, route.dat, by = c('Route', 'CountryNum', 'StateNum'))
# Grab all BBS data in 2020 east of the 100th meridian
bbs.dat.full <- bbs.dat.full %>%
  filter(Year == 2021, RPID == 101, Longitude > -100)
# Select columns of interest
bbs.dat.full <- bbs.dat.full %>%
  dplyr::select(RouteDataID, RouteName, Latitude, Longitude, Year, AOU, starts_with("Count"), 
                -CountryNum)
# Fill in implicit zeros. For this situation, gives a value for each species in
# each existing combination of RouteDataID, Latitude, Longitude, and Year
# RouteDataID is a unique identifier for each combo of CountryNum, StateNum, Route, RPID, and Year
bbs.dat.full <- bbs.dat.full %>%
  complete(AOU, nesting(RouteDataID, RouteName, Latitude, Longitude, Year))
# Replace NAs with 0s for all columns at once.
bbs.dat.full <- bbs.dat.full %>%
  replace(is.na(.), 0)
# Extract community of interior forest obliage birds
aou.info <- AOU_species_codes %>%
  mutate(AOU = spp.num)
bbs.dat.full <- left_join(bbs.dat.full, aou.info, by = c("AOU"))
# Work with community of eastern forest birds
comm.group.dat <- read.csv("data/bird-species-table-bateman.csv")
my.sp.code <- comm.group.dat %>%
  filter(Group %in% c("Eastern.Forests"))
bbs.dat <- bbs.dat.full %>%
  mutate(alpha.code = as.character(alpha.code)) %>%
  filter(alpha.code %in% my.sp.code$Code)

# Clear out big files and free memory
rm(bbs.dat.full, route.dat)
gc()

# Get associated weather data
weather.dat <- read.csv("data/case-study-1/bbs-2022/weather.csv")

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
# Three dimensional array of detection-nondetection data.
y <- array(NA, dim = c(N, J, K))
# Detection covariates 
day <- rep(NA, J)
# Unique site indices for linking
site.ordered.indx <- sort(unique(bbs.dat$site.index))
# Calculate
for (j in 1:J) {
  tmp <- bbs.dat %>%
    filter(site.index == site.ordered.indx[j])
  day[j] <- first(tmp$julian)
}
for (j in 1:J) {
  print(paste("Currently on site ", j, " out of ", J, sep = ''))
  # Get detection-nondetection data for each species. 
  # group_by is necessary to only grab data from one observer for the rare 
  # occassions where there are multiple observers for a given route
  # in a given year. 
  y[, j, ] <- bbs.dat %>% 
    filter(site.index == site.ordered.indx[j]) %>% 
    arrange(AOU) %>% 
    dplyr::select(starts_with('Count'), -starts_with('CountryNum')) %>%
    as.matrix()
}

coords.sf <- st_as_sf(data.frame(coords),
		      coords = c("X", "Y"),
		      crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
# Albers equal area across contiguous US.
coords.sf.albers <- coords.sf %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs")
# Get coordinates in Albers Equal Area (in kilometers)
coords.albers <- st_coordinates(coords.sf.albers)
J <- nrow(coords.albers)

# Only grab species used to fit the full trend models ---------------------
load("data/case-study-1/final-spOccupancy-data.rda")
sp.names.fit <- dimnames(data.list$y)[[1]]
tmp <- which(sp.names %in% sp.names.fit)
y <- y[tmp, , ]
y <- ifelse(y > 0, 1, 0)
dimnames(y)[[1]] <- sp.names.fit

# Format data for spOccupancy ---------------------------------------------
det.covs <- list(day = day)
data.list.2021 <- list(y = y,
		       occ.covs = list(),
		       det.covs = det.covs,
		       coords = coords.albers)

# Add in BCRs to the occurrence covariates --------------------------------
bcrs <- st_read(dsn = "data/case-study-1/BCR_Terrestrial/", layer = "BCR_Terrestrial_master")
my.proj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"
coords.sf <- st_as_sf(as.data.frame(data.list.2021$coords),
		      coords = c("X", "Y"),
		      crs = my.proj)

bcrs.albers <- bcrs %>%
  st_transform(crs = my.proj)

# Join points toBCRs
tmp <- st_join(coords.sf, bcrs.albers, join = st_intersects)

bcr.factor <- tmp$BCR
# Fill in values that are messed up
bcr.factor[834] <- 31
bcr.factor[1099] <- 13
bcr.factor[1166] <- 30

data.list.2021$occ.covs$BCR <- bcr.factor

# Add in replicate as a covariate effect on detection
rep.val <- array(NA, dim = dim(data.list.2021$y[1, , ]))
for (i in 1:dim(rep.val)[2]) {
  rep.val[, i] <- i
}
data.list.2021$det.covs$rep.val <- rep.val

# Determine which points fall within the species range --------------------
load('data/case-study-1/spOcc-bbs-data.rda')
load("data/case-study-1/bird-life-processed.rda")
sp.names.big <- dimnames(data.list$y)[[1]]
sp.names <- dimnames(data.list.2021$y)[[1]]
names(within.sf.list) <- sp.names.big
range.ind <- matrix(0, nrow(data.list.2021$y), ncol(data.list.2021$y))
for (i in 1:nrow(data.list.2021$y)) {
  curr.indx <- which(names(within.sf.list) == sp.names[i])
  indx <- unlist(c(st_contains(within.sf.list[[curr.indx]], coords.sf)))
  range.ind[i, indx] <- 1
}
data.list.2021$range.ind <- range.ind

# Get maximum temperature -------------------------------------------------
coords.lat.long <- coords.sf %>%
  st_transform(crs = 4326)
# Time period to extract
period <- "19812010"
J <- nrow(coords.sf)
# Variables to extract
tmax <- rep(NA, J)
tmp <- getTerraClimNormals(coords.lat.long, 'tmax', period, 1:12)[[1]]
plt_dat <- extract(tmp, coords.lat.long)
tmax <- apply(plt_dat[,2:ncol(plt_dat)], 1, max)
data.list.2021$occ.covs$tmax <- tmax
data.list.2021$occ.covs <- as.data.frame(data.list.2021$occ.covs)

# Save results for single-season model in spOccupancy ---------------------
save(data.list.2021, file = 'data/case-study-1/spOccupancy-2021-data.rda')
