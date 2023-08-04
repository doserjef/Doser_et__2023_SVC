rm(list = ls())
library(spAbundance)
library(tidyverse)
library(sf)

# Get info from command line ----------------------------------------------
# This code is to extract the current species name from the command line
# to easily run the script for different species
args <- commandArgs(trailingOnly = TRUE)
# Current species
curr.state <- args[1]
# Alternatively, if not running the script from the command line:
curr.state <- 'new york'
if(length(args) == 0) base::stop('Need to tell spAbundance the current state')
print(curr.state)

# Load data ---------------------------------------------------------------
load('data/tree/stage-2-data.rda')

# Filter to one state of interest -----------------------------------------
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
usa <- usa %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs")
state <- usa %>% filter(ID == 'new york')
coords.sf <- st_as_sf(as.data.frame(data.list.2$coords),
			   coords = c('x', 'y'),
			   crs = st_crs(usa))
site.indx <- which(sapply(st_within(coords.sf, st_make_valid(state)), length) == 1)
data.list.2$y <- data.list.2$y[site.indx]
data.list.2$coords <- data.list.2$coords[site.indx, ]
data.list.2$covs <- data.list.2$covs[site.indx, ]
data.list.2$z <- data.list.2$z[site.indx]

# Remove 25% of locations -------------------------------------------------
load(paste0("data/tree/pred-indx-", str_replace_all(curr.state, ' ', '-'), '.rda'))
data.list.2$y <- data.list.2$y[-pred.indx]
data.list.2$coords <- data.list.2$coords[-pred.indx, ]
data.list.2$covs <- data.list.2$covs[-pred.indx, ]
data.list.2$z <- data.list.2$z[-pred.indx]

# Load results ------------------------------------------------------------
# SVC model
load(paste0("results/tree/stage-2-sugar-maple-", 
	   str_replace_all(curr.state, ' ', '-'), "-SVC-samples.rda"))
out.svc <- out
summary(out.svc)
# SVI model
load(paste0("results/tree/stage-2-sugar-maple-", 
	   str_replace_all(curr.state, ' ', '-'), "-SVI-samples.rda"))
out.svi <- out
summary(out.svi)

y.rep.means.svc <- apply(out.svc$y.rep.samples, 2, mean)
y.rep.means.svi <- apply(out.svi$y.rep.samples, 2, mean)
y.true <- data.list.2$y

par(mfrow = c(1, 2))
plot(y.true, y.rep.means.svc, pch = 19, main = 'SVC')
abline(0, 1)
plot(y.true, y.rep.means.svi, pch = 19, main = 'SVI')
abline(0, 1)
par(mfrow = c(1, 1))



