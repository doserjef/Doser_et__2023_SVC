rm(list = ls())
library(spAbundance)
library(tidyverse)
library(sf)
library(scoringRules)

# Load data ---------------------------------------------------------------
load('data/tree/spAbundance-data.rda')

# Filter to one state of interest -----------------------------------------
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
# Restrict to east of the 100th meridian
usa.bbox <- st_bbox(usa)
usa.bbox[1] <- -100
usa.bbox <- as.vector(usa.bbox)
sf_use_s2(FALSE)
east.us <- st_crop(st_make_valid(usa), xmin = usa.bbox[1], ymin = usa.bbox[2],
                   xmax = usa.bbox[3], ymax = usa.bbox[4])
east.us <- east.us %>%
  st_transform(st_crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"))
coords.sf <- st_as_sf(data.frame(data.list$coords),
		      coords = c('x', 'y'),
		      crs = st_crs(east.us))

# Load results ------------------------------------------------------------
# SVC model
# load('results/tree/ne-SVC-results.rda')
# out.svc <- out
# summary(out.svc)
# # SVI model
# load('results/tree/ne-SVI-results.rda')
# out.svi <- out
# summary(out.svi)
# 
# y.rep.means.svc <- apply(out.svc$y.rep.samples^2, 2, mean)
# y.rep.means.svi <- apply(out.svi$y.rep.samples^2, 2, mean)
# y.true <- data.list$y
# 
# par(mfrow = c(1, 2))
# plot(y.true, y.rep.means.svc, pch = 19, main = 'SVC')
# abline(0, 1)
# plot(y.true, y.rep.means.svi, pch = 19, main = 'SVI')
# abline(0, 1)
# par(mfrow = c(1, 1))

# Compare prediction results ----------------------------------------------
load('results/tree/hold-out-vermont-SVC-results.rda')
out.pred.svc <- out.pred
load('results/tree/hold-out-vermont-SVI-results.rda')
out.pred.svi <- out.pred
load('data/tree/pred-indx.rda')
y.0 <- data.list$y[pred.indx]

# Calculate CRPS ----------------------
crps.svc <- crps_sample(sqrt(y.0), t(out.pred.svc$y.0.samples))
crps.svi <- crps_sample(sqrt(y.0), t(out.pred.svi$y.0.samples))

mean(crps.svc)
mean(crps.svi)

# Covergage ---------------------------
y.rep.svc.quants <- apply(out.pred.svc$y.0.samples, 2, quantile, c(0.025, 0.975))
y.rep.svi.quants <- apply(out.pred.svi$y.0.samples, 2, quantile, c(0.025, 0.975))
svc.coverage <- mean(sqrt(y.0) >= y.rep.svc.quants[1, ] & sqrt(y.0) <= y.rep.svc.quants[2, ])
svi.coverage <- mean(sqrt(y.0) >= y.rep.svi.quants[1, ] & sqrt(y.0) <= y.rep.svi.quants[2, ])
svc.coverage
svi.coverage

# Calculate RMSPE ---------------------
y.0.med.svc <- apply(out.pred.svc$y.0.samples, 2, median)
y.0.med.svi <- apply(out.pred.svi$y.0.samples, 2, median)
rmspe.svc <- sqrt(mean((y.0.med.svc - sqrt(y.0))^2))
rmspe.svi <- sqrt(mean((y.0.med.svi - sqrt(y.0))^2))
rmspe.svc
rmspe.svi

# Plot the estimates ------------------
par(mfrow = c(1, 2))
plot(sqrt(y.0), y.0.med.svc, pch = 19)
abline(0, 1)
plot(sqrt(y.0), y.0.med.svi, pch = 19)
abline(0, 1)
par(mfrow = c(1, 1))

# Credible interval widhts ------------
mu.0.quants.svc <- apply(out.pred.svc$mu.0.samples, 2, quantile, c(0.025, 0.5, 0.975))
mu.0.quants.svi <- apply(out.pred.svi$mu.0.samples, 2, quantile, c(0.025, 0.5, 0.975))
mu.0.svc.width <- mu.0.quants.svc[3, ] - mu.0.quants.svc[1, ]
mu.0.svi.width <- mu.0.quants.svi[3, ] - mu.0.quants.svi[1, ]
plot(mu.0.svc.width, mu.0.svi.width, pch = 19)
abline(0, 1)
