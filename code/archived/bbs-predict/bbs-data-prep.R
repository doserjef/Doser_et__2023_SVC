# bbs-data-prep.R: this script takes the raw BBS data and preps it for analysis. Note
#                  that this script does not generate the covariate data.
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
#   filter(Year == 2019, RPID == 101)
# # Select columns of interest
# bbs.dat <- bbs.dat %>%
#   dplyr::select(RouteDataID, RouteName, Latitude, Longitude, AOU, starts_with("Count")) %>%
#   # Data set only include US records, so can remove CountryNum
#   dplyr::select(-CountryNum)
# # Fill in implicit zeros. For this situation, gives a value for each species in 
# # each existing combination of RouteDataID, Latitude, Longitude, and Year
# # RouteDataID is a unique identifier for each combo of CountryNum, StateNum, Route, RPID, and Year
# bbs.dat <- bbs.dat %>%
#   complete(AOU, nesting(RouteDataID, RouteName, Latitude, Longitude))
# # Replace NAs with 0s for all columns at once.
# bbs.dat <- bbs.dat %>%
#   replace(is.na(.), 0)
# # Save the BBS data to avoid having to run the above memory intensive code a bunch.
# save(bbs.dat, file = "data/bbs-predict/bbs-raw-data.rda")

# Load in the above BBS data ----------------------------------------------
load("data/bbs-predict/bbs-raw-data.rda")
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

# Get data in format for spAbundance ---------------------------------------
bbs.dat <- left_join(bbs.dat, weather.covs, by = c('RouteDataID'))
# Sort data by species, Year, Longitude, Latitude
bbs.dat <- bbs.dat %>%
  arrange(AOU, Longitude, Latitude)
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

# Detection Covariates ----------------
J <- nrow(coords)
day <- rep(NA, J)
obs <- rep(NA, J)
# Unique site indices for linking
site.ordered.indx <- sort(unique(bbs.dat$site.index))
# Calculate
for (j in 1:J) {
  tmp <- bbs.dat %>%
    filter(site.index == site.ordered.indx[j])
  day[j] <- first(tmp$julian)
  obs[j] <- first(tmp$ObsN)
}

# Get observed species richness -------
y <- bbs.dat %>%
  mutate(totalCount = ifelse(Count10 + Count20 + Count30 + Count40 + Count50 > 0, 1, 0)) %>%
  group_by(site.index) %>%
  summarize(obs.rich = sum(totalCount))

# Save the data -----------------------------------------------------------
save(coords, y, obs, day, file = "data/bbs-predict/bbs-rich-data.rda")
