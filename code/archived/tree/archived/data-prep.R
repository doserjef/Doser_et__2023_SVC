rm(list = ls())
library(tidyverse)
library(sf)

# Load data ---------------------------------------------------------------
# Don't move these somewhere else
load("~/../andy/data_actual/actual_forest_plot_spp_data.RData")
# Spatial coordinates
load("~/../andy/data_actual/actual_forest_plot_coords.RData")
# Covariates
load("~/../andy/data_actual/X_vars.rda")

coords <- actual_forest_plot_coords %>%
  select(x = LON, y = LAT) %>%
  as.matrix()
# Convert to matrix
# Number of sites
J <- nrow(coords)
# Number of species
N <- n_distinct(actual_forest_plot_spp_data$SPCD)
# Species ids
sp.id <- unique(actual_forest_plot_spp_data$SPCD)
# Get data in species x site matrix.
y <- actual_forest_plot_spp_data %>%
  pull(SP_PRESENT)
y <- matrix(y, N, J)
rownames(y) <- sp.id
y.bio <- actual_forest_plot_spp_data %>%
  pull(BIO_ACRE)
y.bio <- matrix(y.bio, N, J)
rownames(y.bio) <- sp.id

# Format data for spOccupancy ---------------------------------------------
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
usa <- usa %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs")
coords.sf <- st_as_sf(data.frame(coords),
		      coords = c('x', 'y'),
		      crs = st_crs(usa))

# Reset all the variables
J <- nrow(coords)
# 318 is sugar maple. That is the species we're going to focus on. 
curr.sp.indx <- which(sp.id == 318) 
curr.indx <- 1:J
data.list <- list(y = y[curr.sp.indx, ],
		    coords = coords,
		    covs = X)
data.list$weights <- rep(1, J)
save(data.list, file = 'data/tree/stage-1-data.rda')

# Format Stage 2 data for spAbundance -------------------------------------
# Get ecoregion for splitting up data into smaller sets.
ecr <- st_read(dsn = "~/DFISD22/data/ecoregions/", layer = "us_eco_l3")
coords.sf <- st_as_sf(as.data.frame(data.list$coords),
		      coords = c("x", "y"),
		      crs = st_crs(usa))

# Get the state that each
ecrs.albers <- ecr %>%
  st_transform(crs = st_crs(usa))

# Join points toBCRs
tmp <- st_join(coords.sf, ecrs.albers, join = st_intersects)
# Sites where species was observed
# site.indx <- which(data.list$y == 1)
data.list.2 <- list(y = ifelse(y.bio[curr.sp.indx, ] == 0, 0, log(y.bio[curr.sp.indx, ])),
		    coords = coords,
		    covs = as.data.frame(X), 
                    z = data.list$y)
# data.list.2$covs$ecoregionL3 <- as.numeric(tmp$US_L3CODE)[site.indx]
# Fill in missing indices for ecoregion
# miss.indx <- which(is.na(data.list.2$covs$ecoregionL3))
# # NOTE: hardcoded
# data.list.2$covs$ecoregionL3[c(692, 715, 793, 947)] <- 50
# data.list.2$covs$ecoregionL3[c(6415, 6426)] <- 51

save(data.list.2, file = 'data/tree/stage-2-data.rda')
