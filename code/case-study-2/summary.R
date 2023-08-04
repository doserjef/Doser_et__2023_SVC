# summary.R: script to summarize results from the Grasshopper Sparrow 
#            case study.
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(sf)
library(spOccupancy)
library(viridis)
library(pals)
library(ggpubr)
library(stars)

# Load in data sets -------------------------------------------------------
load("data/case-study-2/GRSP-spOcc-data.rda")
load("data/case-study-2/GRSP-pred-data.rda")
load("data/case-study-2/species-BirdLife-ranges.rda")

# Set up plotting stuff ---------------------------------------------------
my.crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"

# Map of US states
usa <- st_as_sf(maps::map("usa", fill = TRUE, plot = FALSE))
usa <- usa %>%
  st_transform(my.crs)

# Compare WAIC across models ----------------------------------------------
load('results/case-study-2/waic-results.rda')
waic.df <- waic.df %>%
  mutate(waic.diff = waic - min(waic))
waic.df

# Covariates plot ---------------------------------------------------------
plot.df <- data.frame(x = coords.0[, 1], 
		      y = coords.0[, 2],
		      grass.change = occ.pred.covs$grass[, 50] - occ.pred.covs$grass[, 1],
		      crop.change = occ.pred.covs$crop[, 50] - occ.pred.covs$crop[, 1],
                      grass.mean = occ.pred.covs$grass.mean,
		      tmax = occ.pred.covs$tmax,
                      crop.mean = occ.pred.covs$crop.mean)
plot.df <- st_as_stars(plot.df, dims = c('x', 'y'))
grass.change.plot <- ggplot() +
  geom_stars(data = plot.df, aes(x = x, y = y, fill = grass.change), interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0, col = 'grey') +
  geom_sf(data = GRSP.range, alpha = 0, col = 'black') +
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC',
  	               na.value = NA) +
  theme_bw(base_size = 15) +
  labs(x = "Longitude", y = "Latitude", fill = "",
       title = '(C) Change in grassland area 1970-2019') +
  theme(legend.position = c(0.93, 0.29),
        legend.background = element_rect(fill = NA),
	text = element_text(family = 'LM Roman 10'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 12),
	plot.title = element_text(size = 15),
        legend.text = element_text(size = 10))
crop.change.plot <- ggplot() +
  geom_stars(data = plot.df, aes(x = x, y = y, fill = crop.change), interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0, col = 'grey') +
  geom_sf(data = GRSP.range, alpha = 0, col = 'black') +
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC',
  	               na.value = NA) +
  theme_bw(base_size = 15) +
  labs(x = "Longitude", y = "Latitude", fill = "",
       title = '(D) Change in cropland area 1970-2019') +
  theme(legend.position = c(0.93, 0.29),
        legend.background = element_rect(fill = NA),
	text = element_text(family = 'LM Roman 10'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 12),
	plot.title = element_text(size = 15),
        legend.text = element_text(size = 10))
grass.mean.plot <- ggplot() +
  geom_stars(data = plot.df, aes(x = x, y = y, fill = grass.mean), interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0, col = 'grey') +
  geom_sf(data = GRSP.range, alpha = 0, col = 'black') +
  scale_fill_viridis_c(na.value = NA) +
  theme_bw(base_size = 15) +
  labs(x = "Longitude", y = "Latitude", fill = "",
       title = '(A) Average grassland area 1970-2019') +
  theme(legend.position = c(0.93, 0.29),
        legend.background = element_rect(fill = NA),
	text = element_text(family = 'LM Roman 10'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 12),
	plot.title = element_text(size = 15),
        legend.text = element_text(size = 10))
crop.mean.plot <- ggplot() +
  geom_stars(data = plot.df, aes(x = x, y = y, fill = crop.mean), interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0, col = 'grey') +
  geom_sf(data = GRSP.range, alpha = 0, col = 'black') +
  scale_fill_viridis_c(na.value = NA) +
  theme_bw(base_size = 15) +
  labs(x = "Longitude", y = "Latitude", fill = "",
       title = '(B) Average cropland area 1970-2019') +
  theme(legend.position = c(0.93, 0.29),
        legend.background = element_rect(fill = NA),
	text = element_text(family = 'LM Roman 10'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 12),
	plot.title = element_text(size = 15),
        legend.text = element_text(size = 10))
# Generate Supplemental Information S3 Figure 1
ggarrange(grass.mean.plot, crop.mean.plot, grass.change.plot, crop.change.plot, 
	  nrow = 2, ncol = 2)
ggsave(file = 'figures/case-study-2/Supp-Info-S3-Figure-1.png', device = 'png', 
       units = 'in', height = 8, width = 12, bg = 'white')

# Generate maps of the covariate effects ----------------------------------
# Constant model ----------------------
load("results/case-study-2/constant-GRSP-svc-samples.rda")
grass.dev.quantiles <- quantile(beta.star.grass.dev.samples, c(0.025, 0.5, 0.975)) 
crop.dev.quantiles <- quantile(beta.star.crop.dev.samples, c(0.025, 0.5, 0.975)) 
plot.df <- data.frame(x = coords.0[, 1], 
		      y = coords.0[, 2],
		      grass.med = grass.dev.quantiles[2], 
		      grass.low = grass.dev.quantiles[1],
		      grass.high = grass.dev.quantiles[3],
		      crop.med = crop.dev.quantiles[2], 
		      crop.low = crop.dev.quantiles[1],
		      crop.high = crop.dev.quantiles[3])
plot.df <- st_as_stars(plot.df, dims = c('x', 'y'))
grass.constant.plot <- ggplot() +
  geom_stars(data = plot.df, aes(x = x, y = y, fill = grass.med), interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0, col = 'grey') +
  geom_sf(data = GRSP.range, alpha = 0, col = 'black') +
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC',
  	               na.value = NA, limits = c(-0.5, 0.5)) +
  theme_bw(base_size = 15) +
  labs(x = "Longitude", y = "Latitude", fill = "",
       title = 'Constant Model Grassland Effect') +
  guides(fill = 'none') +
  theme(legend.position = c(0.87, 0.25),
        legend.background = element_rect(fill = NA),
	text = element_text(family = 'LM Roman 10'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 12),
	plot.title = element_text(size = 15),
        legend.text = element_text(size = 12))
crop.constant.plot <- ggplot() +
  geom_stars(data = plot.df, aes(x = x, y = y, fill = crop.med), interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0, col = 'grey') +
  geom_sf(data = GRSP.range, alpha = 0, col = 'black') +
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC',
  	               na.value = NA, limits = c(-0.5, 0.5)) +
  theme_bw(base_size = 15) +
  labs(x = "Longitude", y = "Latitude", fill = "",
       title = 'Constant Model Cropland Effect') +
  guides(fill = 'none') +
  theme(legend.position = c(0.87, 0.25),
        legend.background = element_rect(fill = NA),
	text = element_text(family = 'LM Roman 10'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 12),
	plot.title = element_text(size = 15),
        legend.text = element_text(size = 12))
rm(beta.star.grass.dev.samples, beta.star.crop.dev.samples)
gc()
# Landcover interaction model ---------
load("results/case-study-2/int-lulc-GRSP-svc-samples.rda")
grass.dev.quantiles <- apply(beta.star.grass.dev.samples, 2, quantile, c(0.025, 0.5, 0.975))
crop.dev.quantiles <- apply(beta.star.crop.dev.samples, 2, quantile, c(0.025, 0.5, 0.975))
grass.prob.pos <- apply(beta.star.grass.dev.samples, 2, function(a) mean(a > 0))
crop.prob.pos <- apply(beta.star.crop.dev.samples, 2, function(a) mean(a > 0))
plot.df <- data.frame(x = coords.0[, 1], 
		      y = coords.0[, 2],
		      grass.med = grass.dev.quantiles[2, ], 
		      grass.low = grass.dev.quantiles[1, ],
		      grass.high = grass.dev.quantiles[3, ],
		      crop.med = crop.dev.quantiles[2, ], 
		      crop.low = crop.dev.quantiles[1, ],
		      crop.high = crop.dev.quantiles[3, ],
		      grass.prob.pos = grass.prob.pos,
		      crop.prob.pos = crop.prob.pos) 
plot.df <- st_as_stars(plot.df, dims = c('x', 'y'))
grass.int.lulc.plot <- ggplot() +
  geom_stars(data = plot.df, aes(x = x, y = y, fill = grass.med), interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0, col = 'grey') +
  geom_sf(data = GRSP.range, alpha = 0, col = 'black') +
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC',
  	               na.value = NA) +
  theme_bw(base_size = 15) +
  labs(x = "Longitude", y = "Latitude", fill = "",
       title = '(A) Grassland Change x Grassland Amount') +
  guides(fill = 'none') +
  theme(legend.position = c(0.87, 0.25),
        legend.background = element_rect(fill = NA),
	text = element_text(family = 'LM Roman 10'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 12),
	plot.title = element_text(size = 15),
        legend.text = element_text(size = 12))
crop.int.lulc.plot <- ggplot() +
  geom_stars(data = plot.df, aes(x = x, y = y, fill = crop.med), interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0, col = 'grey') +
  geom_sf(data = GRSP.range, alpha = 0, col = 'black') +
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC',
  	               na.value = NA) +
  theme_bw(base_size = 15) +
  labs(x = "Longitude", y = "Latitude", fill = "",
       title = '(D) Cropland Change x Cropland Amount') +
  guides(fill = 'none') +
  theme(legend.position = c(0.87, 0.25),
        legend.background = element_rect(fill = NA),
	text = element_text(family = 'LM Roman 10'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 12),
	plot.title = element_text(size = 15),
        legend.text = element_text(size = 12))
rm(beta.star.grass.dev.samples, beta.star.crop.dev.samples)
gc()
# TMAX interaction model --------------
load("results/case-study-2/int-tmax-GRSP-svc-samples.rda")
grass.dev.quantiles <- apply(beta.star.grass.dev.samples, 2, quantile, 
			     c(0.025, 0.5, 0.975), na.rm = TRUE)
crop.dev.quantiles <- apply(beta.star.crop.dev.samples, 2, quantile, 
			    c(0.025, 0.5, 0.975), na.rm = TRUE)
grass.prob.pos <- apply(beta.star.grass.dev.samples, 2, function(a) mean(a > 0, na.rm = TRUE))
crop.prob.pos <- apply(beta.star.crop.dev.samples, 2, function(a) mean(a > 0, na.rm = TRUE))
plot.df <- data.frame(x = coords.0[, 1], 
		      y = coords.0[, 2],
		      grass.med = grass.dev.quantiles[2, ], 
		      grass.low = grass.dev.quantiles[1, ],
		      grass.high = grass.dev.quantiles[3, ],
		      crop.med = crop.dev.quantiles[2, ], 
		      crop.low = crop.dev.quantiles[1, ],
		      crop.high = crop.dev.quantiles[3, ],
		      grass.prob.pos = grass.prob.pos,
		      crop.prob.pos = crop.prob.pos) 
plot.df <- st_as_stars(plot.df, dims = c('x', 'y'))
grass.int.tmax.plot <- ggplot() +
  geom_stars(data = plot.df, aes(x = x, y = y, fill = grass.med), interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0, col = 'grey') +
  geom_sf(data = GRSP.range, alpha = 0, col = 'black') +
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC',
  	               na.value = NA) +
  theme_bw(base_size = 15) +
  labs(x = "Longitude", y = "Latitude", fill = "",
       title = '(B) Grassland Change x Max Temp') +
  guides(fill = 'none') +
  theme(legend.position = c(0.87, 0.25),
        legend.background = element_rect(fill = NA),
	text = element_text(family = 'LM Roman 10'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 12),
	plot.title = element_text(size = 15),
        legend.text = element_text(size = 12))
crop.int.tmax.plot <- ggplot() +
  geom_stars(data = plot.df, aes(x = x, y = y, fill = crop.med), interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0, col = 'grey') +
  geom_sf(data = GRSP.range, alpha = 0, col = 'black') +
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC',
  	               na.value = NA) +
  theme_bw(base_size = 15) +
  labs(x = "Longitude", y = "Latitude", fill = "",
       title = '(E) Cropland Change x Max Temp') +
  guides(fill = 'none') +
  theme(legend.position = c(0.87, 0.25),
        legend.background = element_rect(fill = NA),
	text = element_text(family = 'LM Roman 10'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 12),
	plot.title = element_text(size = 15),
        legend.text = element_text(size = 12))
rm(beta.star.grass.dev.samples, beta.star.crop.dev.samples)
gc()
# SVC model ---------------------------
load("results/case-study-2/svc-GRSP-svc-samples.rda")
grass.dev.quantiles <- apply(beta.star.grass.dev.samples, 2, quantile, c(0.025, 0.5, 0.975))
crop.dev.quantiles <- apply(beta.star.crop.dev.samples, 2, quantile, c(0.025, 0.5, 0.975))
grass.prob.pos <- apply(beta.star.grass.dev.samples, 2, function(a) mean(a > 0))
crop.prob.pos <- apply(beta.star.crop.dev.samples, 2, function(a) mean(a > 0))
plot.df <- data.frame(x = coords.0[, 1], 
		      y = coords.0[, 2],
		      grass.med = grass.dev.quantiles[2, ], 
		      grass.low = grass.dev.quantiles[1, ],
		      grass.high = grass.dev.quantiles[3, ],
		      crop.med = crop.dev.quantiles[2, ], 
		      crop.low = crop.dev.quantiles[1, ],
		      crop.high = crop.dev.quantiles[3, ],
		      grass.prob.pos = grass.prob.pos,
		      crop.prob.pos = crop.prob.pos, 
                      x.grass.dev.end = occ.pred.covs$grass.dev[, 50], 
                      x.grass.mean = occ.pred.covs$grass.mean,
                      x.crop.dev.end = occ.pred.covs$crop.dev[, 50], 
                      x.crop.mean = occ.pred.covs$crop.mean)
plot.df <- st_as_stars(plot.df, dims = c('x', 'y'))
grass.svc.plot <- ggplot() +
  geom_stars(data = plot.df, aes(x = x, y = y, fill = grass.med), interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0, col = 'grey') +
  geom_sf(data = GRSP.range, alpha = 0, col = 'black') +
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC',
  	               na.value = NA) +
  theme_bw(base_size = 15) +
  labs(x = "Longitude", y = "Latitude", fill = "",
       title = 'SVC Model Grassland Effect') +
  guides(fill = 'none') +
  theme(legend.position = c(0.87, 0.25),
        legend.background = element_rect(fill = NA),
	text = element_text(family = 'LM Roman 10'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 12),
	plot.title = element_text(size = 15),
        legend.text = element_text(size = 12))
crop.svc.plot <- ggplot() +
  geom_stars(data = plot.df, aes(x = x, y = y, fill = crop.med), interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0, col = 'grey') +
  geom_sf(data = GRSP.range, alpha = 0, col = 'black') +
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC',
  	               na.value = NA) +
  theme_bw(base_size = 15) +
  guides(fill = 'none') +
  labs(x = "Longitude", y = "Latitude", fill = "",
       title = 'SVC Model Cropland Effect') +
  theme(legend.position = c(0.87, 0.25),
        legend.background = element_rect(fill = NA),
	text = element_text(family = 'LM Roman 10'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 12),
	plot.title = element_text(size = 15),
        legend.text = element_text(size = 12))
rm(beta.star.grass.dev.samples, beta.star.crop.dev.samples)
gc()
# Full model --------------------------
load("results/case-study-2/full-GRSP-svc-samples.rda")
grass.dev.quantiles <- apply(beta.star.grass.dev.samples, 2, quantile, c(0.025, 0.5, 0.975), na.rm = TRUE)
crop.dev.quantiles <- apply(beta.star.crop.dev.samples, 2, quantile, c(0.025, 0.5, 0.975), na.rm = TRUE)
grass.prob.pos <- apply(beta.star.grass.dev.samples, 2, function(a) mean(a > 0, na.rm = TRUE))
crop.prob.pos <- apply(beta.star.crop.dev.samples, 2, function(a) mean(a > 0, na.rm = TRUE))
full.plot.df <- data.frame(x = coords.0[, 1], 
		      y = coords.0[, 2],
		      grass.med = grass.dev.quantiles[2, ], 
		      grass.low = grass.dev.quantiles[1, ],
		      grass.high = grass.dev.quantiles[3, ],
		      crop.med = crop.dev.quantiles[2, ], 
		      crop.low = crop.dev.quantiles[1, ],
		      crop.high = crop.dev.quantiles[3, ],
		      grass.prob.pos = grass.prob.pos,
		      crop.prob.pos = crop.prob.pos, 
                      x.grass.dev.end = occ.pred.covs$grass.dev[, 50], 
                      x.grass.mean = occ.pred.covs$grass.mean,
                      x.crop.dev.end = occ.pred.covs$crop.dev[, 50], 
                      x.crop.mean = occ.pred.covs$crop.mean)
full.plot.df <- st_as_stars(full.plot.df, dims = c('x', 'y'))
grass.full.plot <- ggplot() +
  geom_stars(data = full.plot.df, aes(x = x, y = y, fill = grass.med), interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0, col = 'grey') +
  geom_sf(data = GRSP.range, alpha = 0, col = 'black') +
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC',
  	               na.value = NA) +
  theme_bw(base_size = 15) +
  labs(x = "Longitude", y = "Latitude", fill = "",
       title = '(C) Grassland Change SVC + Interactions') +
  guides(fill = 'none') +
  theme(legend.position = c(0.87, 0.25),
        legend.background = element_rect(fill = NA),
	text = element_text(family = 'LM Roman 10'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 12),
	plot.title = element_text(size = 15),
        legend.text = element_text(size = 12))
crop.full.plot <- ggplot() +
  geom_stars(data = full.plot.df, aes(x = x, y = y, fill = crop.med), interpolate = TRUE) +
  geom_sf(data = usa, alpha = 0, col = 'grey') +
  geom_sf(data = GRSP.range, alpha = 0, col = 'black') +
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC',
  	               na.value = NA) +
  theme_bw(base_size = 15) +
  labs(x = "Longitude", y = "Latitude", fill = "",
       title = '(F) Cropland Change SVC + Interactions') +
  guides(fill = 'none') +
  theme(legend.position = c(0.87, 0.25),
        legend.background = element_rect(fill = NA),
	text = element_text(family = 'LM Roman 10'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 12),
	plot.title = element_text(size = 15),
        legend.text = element_text(size = 12))
rm(beta.star.grass.dev.samples, beta.star.crop.dev.samples)
gc()

# Generate Figure 4
ggarrange(grass.int.lulc.plot, grass.int.tmax.plot, grass.full.plot, 
	  crop.int.lulc.plot, crop.int.tmax.plot, crop.full.plot, nrow = 2, ncol = 3)
ggsave(file = 'figures/case-study-2/Figure-4.png', device = 'png', height = 7, 
       width = 15, units = 'in', bg = 'white')


