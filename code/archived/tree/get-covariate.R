rm(list = ls())

library(tidyverse)
library(magrittr)
library(stars)
library(sf)
library(FedData)

# Load data ---------------------------------------------------------------
load('data/tree/vermont-bio-coords.rda')

# Change order to help speed things up
ord <- order(coords[, 1])
# coords <- coords[ord, ]

# Filter to NE, area of interest ------------------------------------------
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
usa <- usa %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs")
coords.sf <- st_as_sf(as.data.frame(coords),
			   coords = c('x', 'y'),
			   crs = st_crs(usa))

coords.sf.albers <- coords.sf

# Get percent forest cover from NLCD --------------------------------------
# Create matrix to store all values for grassland and cropland variables
J <- nrow(coords.sf.albers)
# Loop through all the sites
vals <- split(1:J, ceiling(seq_along(1:J)/100))
# # Function to get proportion forest cover
# props <- function(a, na.rm = TRUE) {
#   my.sum <- sum(!is.na(a))
#   prop.for <- sum(a %in% c(41, 42, 43), na.rm = na.rm) / my.sum
#   return(prop.for)
# }
my.mean <- function(a) {
  mean(a, na.rm = TRUE)
}
tcc <- rep(NA, J)
for (j in 1:length(vals)) {
  print(paste0("Currently on ", j, " out of ", length(vals)))
  coords.curr <- coords.sf.albers[vals[[j]], ]
  nlcd.dat <- get_nlcd(template = coords.curr, label = paste0('canopy-', j), year = 2016, 
                       dataset = 'canopy')
  # Calculating forest cover in 1km. 
  # forest[vals[[j]]] <- raster::extract(nlcd.dat, coords.curr, buffer = 50, fun = props)
  tcc[vals[[j]]] <- raster::extract(nlcd.dat, coords.curr, buffer = 50, fun = my.mean)
}

# Save results ------------------------------------------------------------
# Extract year from the file name
tcc <- tcc[order(ord)]
tcc <- ifelse(is.na(tcc), 0, tcc)
save(tcc, file = 'data/tree/tree-canopy-cover.rda')


# Exploratory variogram analysis
library(MBA)
library(fields)
library(geoR)
par(mfrow=c(1,2))
out.lm <- lm(sqrt(y) ~ tcc)
curr.res <- out.lm$residuals
vario.1.raw <- variog(coords = coords, data = sqrt(y))
par(mfrow = c(1, 2))
plot(vario.1.raw, pch=16)
vario.1.res <- variog(coords = coords[which(!is.na(tcc)), ], data = curr.res)
plot(vario.1.res, pch = 16)
par(mfrow = c(1, 1))



# An exploratory plot, motivated by May et al. 2022
library(MASS)
library(ggplot2)
library(viridis)
#> Loading required package: viridisLite
theme_set(theme_bw(base_size = 16))

# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
dat <- data.frame(tcc = tcc, y = y)
dat$density <- get_density(dat$tcc, dat$y, n = 100)
ggplot(dat) + geom_point(aes(tcc, sqrt(y), color = density)) + scale_color_viridis()
