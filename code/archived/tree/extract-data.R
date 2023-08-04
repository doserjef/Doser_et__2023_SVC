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

# Subset data to Vermont --------------------------------------------------
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
vermont <- usa
# vermont <- usa %>% filter(ID == 'maine')
# Restrict to east of the 100th meridian
# usa.bbox <- st_bbox(usa)
# usa.bbox[1] <- -100
# usa.bbox <- as.vector(usa.bbox)
# sf_use_s2(FALSE)
# east.us <- st_crop(st_make_valid(usa), xmin = usa.bbox[1], ymin = usa.bbox[2], 
#                    xmax = usa.bbox[3], ymax = usa.bbox[4])
vermont <- vermont %>%
  st_transform(st_crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"))
coords.sf <- st_as_sf(data.frame(coords),
		      coords = c('x', 'y'),
		      crs = st_crs(vermont))

# Indexing vector for the coordinates out vermont
vermont.indx <- st_within(coords.sf, vermont)
vermont.indx <- which(sapply(vermont.indx, function(a) ifelse(length(a) == 0, 0, 1)) == 1)
coords.vermont <- coords[vermont.indx, ]
# Reset all the variables
J <- nrow(coords.vermont)
y.bio.vermont <- y.bio[, vermont.indx]
y.vermont <- y[, vermont.indx]
X.vermont <- X[vermont.indx, ]
# Grab the 10 most common species in Vermont
sp.indx.order <- rev(order(apply(y.bio.vermont, 1, sum, na.rm = TRUE)))
curr.sp.indx <- sp.indx.order[1:10]
curr.indx <- 1:J

y <- apply(y.bio.vermont, 2, sum)
coords <- coords.vermont

# Randomly subset to a couple thousand data points
curr.indx <- sample(1:length(y), 2000, replace = FALSE)

y <- y[curr.indx]
coords <- coords[curr.indx, ]

save(y, coords, file = "data/tree/vermont-bio-coords.rda")

