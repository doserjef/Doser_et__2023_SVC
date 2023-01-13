# summary.R: this script summarizes the case study results in which we 
#            assess spatially-varying trends in 51 forest birds across the 
#            eastern US. 
rm(list = ls())
library(sf)
library(spOccupancy)
library(tidyverse)
library(ggpubr)
library(ggthemes)
library(viridis)
library(coda)
library(stars)
library(pals)

# Load in data used to fit model ------------------------------------------
load("data/spOcc-bbs-data.rda")
load("data/bird-life-processed.rda")
sp.names <- dimnames(data.list$y)[[1]]
sp.names <- sp.names[percent.area > 0.5]
within.sf.list <- within.sf.list[percent.area > 0.5]
N <- length(sp.names)

# Average number of years each route is surveyed. 
n.surveys.per.route <- apply(data.list$y[1, , , 1], 1, function(a) sum(!is.na(a)))
mean(n.surveys.per.route)
sd(n.surveys.per.route)

# Generate table of WAIC values -------------------------------------------
# This code reads in the individual results file from each model fit and 
# calculates the WAIC using the waicOcc() function in spOccupancy. It takes 
# few minutes to read in all the files, so it is commented out and the resulting
# file is saved in "results/model-comparison-results.rda", which we load in 
# below. 
# n.models <- 3
# waic.vals <- matrix(NA, N, n.models)
# colnames(waic.vals) <- c('constant', 'bcr', 'svc')
# rownames(waic.vals) <- sp.names
# for (i in 1:N) {
#   print(i)
#   load(paste("results/stPGOcc-", sp.names[i], ".R", sep = ''))
#   waic.vals[i, 1] <- waicOcc(out)[3]
#   load(paste("results/stPGOcc-bcr-", sp.names[i], ".R", sep = ''))
#   waic.vals[i, 2] <- waicOcc(out)[3]
#   load(paste("results/svcTPGOcc-", sp.names[i], ".R", sep = ''))
#   waic.vals[i, 3] <- waicOcc(out)[3]
# }
#  
# waic.deviance.df <- data.frame(sp.code = sp.names, 
#  		               waic = waic.vals - waic.vals[, 3]) 
# save(waic.deviance.df, file = 'results/model-comparison-results.rda')
# Load in the resulting file. 
load('results/model-comparison-results.rda')
# Table S3 in Appendix S2
waic.deviance.df %>%
  select(sp.code, waic.constant, waic.bcr, waic.svc) %>%
  arrange(sp.code)

# Stuff for plotting ------------------------------------------------------
coords.sf <- st_as_sf(as.data.frame(data.list$coords),
		      coords = c("X", "Y"),
		      crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs")
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
usa.no.states <- st_as_sf(maps::map("usa", fill = TRUE, plot = FALSE))
# Restrict to east of the 100th meridian
usa.bbox <- st_bbox(usa)
usa.bbox[1] <- -100
usa.bbox <- as.vector(usa.bbox)
sf_use_s2(FALSE)
# Full data
east.us <- st_crop(st_make_valid(usa), xmin = usa.bbox[1], ymin = usa.bbox[2], 
                   xmax = usa.bbox[3], ymax = usa.bbox[4])
east.us <- east.us %>%
  st_transform(st_crs(coords.sf))
east.us.no.states <- st_crop(st_make_valid(usa.no.states), xmin = usa.bbox[1], ymin = usa.bbox[2], 
                   xmax = usa.bbox[3], ymax = usa.bbox[4])
east.us.no.states <- east.us.no.states %>%
  st_transform(st_crs(coords.sf))


# BCRS, in case you want to visualize the maps with these
bcrs <- st_read(dsn = "data/BCR_Terrestrial/", layer = "BCR_Terrestrial_master")
my.proj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"
coords.sf <- st_as_sf(as.data.frame(data.list$coords),
		      coords = c("X", "Y"),
		      crs = my.proj)

bcrs.albers <- bcrs %>%
  st_transform(crs = my.proj)

# Clip bcrs to region of ne.states
bcrs.ne <- st_intersection(st_make_valid(bcrs.albers), st_make_valid(east.us))

bcrs.ne.grouped <- bcrs.ne %>%
  group_by(BCR) %>%
  summarize(geometry = st_union(geometry))

# Prediction coordinates
load('data/pred-coordinates.rda')
coords.0.sf <- st_as_sf(as.data.frame(coords.0),
		      coords = c("X", "Y"),
		      crs = st_crs(coords.sf))

# Generate plots for Figure 4 in manuscript -------------------------------
# Example species 1: Tufted Titmouse
curr.sp <- "TUTI"
load(paste('results/predict-svcTPGOcc-', curr.sp, sep = ''))
load(paste("results/svcTPGOcc-", curr.sp, ".R", sep = ''))
out.TUTI <- out
trend.TUTI.quants <- w.trend.quants + quantile(out.TUTI$beta.samples[, 2], c(0.025, 0.5, 0.975))
# Only plot the predicted trends within the range of TUTI
indx <- unlist(c(st_contains(within.sf.list[[which(sp.names == 'TUTI')]], coords.0.sf)))
trend.TUTI.quants[, -indx] <- NA

# Example species 2: Eastern Wood-pewee
curr.sp <- "EAWP"
load(paste('results/predict-svcTPGOcc-', curr.sp, sep = ''))
load(paste("results/svcTPGOcc-", curr.sp, ".R", sep = ''))
out.EAWP <- out
trend.EAWP.quants <- w.trend.quants + quantile(out.EAWP$beta.samples[, 2], c(0.025, 0.5, 0.975))
# Only plot the predicted trends within the range of EAWP 
indx <- unlist(c(st_contains(within.sf.list[[which(sp.names == 'EAWP')]], coords.0.sf)))
trend.EAWP.quants[, -indx] <- NA

# Example species 3: Wood Thrush
curr.sp <- "WOTH"
load(paste('results/predict-svcTPGOcc-', curr.sp, sep = ''))
load(paste("results/svcTPGOcc-", curr.sp, ".R", sep = ''))
out.WOTH <- out
trend.WOTH.quants <- w.trend.quants + quantile(out.WOTH$beta.samples[, 2], c(0.025, 0.5, 0.975))
# Only plot the predicted trends within the range of WOTH
indx <- unlist(c(st_contains(within.sf.list[[which(sp.names == 'WOTH')]], coords.0.sf)))
trend.WOTH.quants[, -indx] <- NA

# Example species 4: Blue-winged Warbler
curr.sp <- "BWWA"
load(paste('results/predict-svcTPGOcc-', curr.sp, sep = ''))
load(paste("results/svcTPGOcc-", curr.sp, ".R", sep = ''))
out.BWWA <- out
trend.BWWA.quants <- w.trend.quants + quantile(out.BWWA$beta.samples[, 2], c(0.025, 0.5, 0.975))
# Only plot the predicted trends within the range of BWWA
indx <- unlist(c(st_contains(within.sf.list[[which(sp.names == 'BWWA')]], coords.0.sf)))
trend.BWWA.quants[, -indx] <- NA

# Put all values in a data frame for prediction. 
plot.df <- data.frame(trend.TUTI = trend.TUTI.quants[2, ],
		      ci.width.TUTI = trend.TUTI.quants[3, ] - trend.TUTI.quants[1, ],
		      trend.EAWP = trend.EAWP.quants[2, ],
		      ci.width.EAWP = trend.EAWP.quants[3, ] - trend.EAWP.quants[1, ],
		      trend.WOTH = trend.WOTH.quants[2, ],
		      ci.width.WOTH = trend.WOTH.quants[3, ] - trend.WOTH.quants[1, ],
		      trend.BWWA = trend.BWWA.quants[2, ],
		      ci.width.BWWA = trend.BWWA.quants[3, ] - trend.BWWA.quants[1, ],
		      x = coords.0[, 1],
		      y = coords.0[, 2])

# Generate maps
pred.stars <- st_as_stars(plot.df, dims = c('x', 'y'))
w.trend.TUTI.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = trend.TUTI), interpolate = TRUE) +
  geom_sf(data = east.us, alpha = 0, col = 'grey') +
  geom_sf(data = within.sf.list[[which(sp.names == 'TUTI')]], alpha = 0, col = 'black') + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA) + 
  theme_bw(base_size = 18) +
  labs(x = "Longitude", y = "Latitude", fill = "Trend\n(logit scale)", 
       title = '(A) Tufted Titmouse') +
  theme(legend.position = c(0.87, 0.22), 
        legend.background = element_rect(fill = NA), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_text(size = 14),
	plot.title = element_text(size = 18),
        legend.text = element_text(size = 12))

w.trend.EAWP.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = trend.EAWP), interpolate = TRUE) +
  geom_sf(data = east.us, alpha = 0, col = 'grey') +
  geom_sf(data = within.sf.list[[which(sp.names == 'EAWP')]], alpha = 0, col = 'black') + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA) + 
  theme_bw(base_size = 18) +
  labs(x = "Longitude", y = "Latitude", fill = "Trend\n(logit scale)", 
       title = '(B) Eastern Wood-pewee') +
  theme(legend.position = c(0.87, 0.22), 
        legend.background = element_rect(fill = NA), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_text(size = 14),
	plot.title = element_text(size = 18),
        legend.text = element_text(size = 12))

w.trend.WOTH.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = trend.WOTH), interpolate = TRUE) +
  geom_sf(data = east.us, alpha = 0, col = 'grey') +
  geom_sf(data = within.sf.list[[which(sp.names == 'WOTH')]], alpha = 0, col = 'black') + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA) + 
  theme_bw(base_size = 18) +
  labs(x = "Longitude", y = "Latitude", fill = "Trend\n(logit scale)", 
       title = '(C) Wood Thrush') +
  theme(legend.position = c(0.87, 0.22), 
        legend.background = element_rect(fill = NA), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_text(size = 14),
	plot.title = element_text(size = 18),
        legend.text = element_text(size = 12))

w.trend.BWWA.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = trend.BWWA), interpolate = TRUE) +
  geom_sf(data = east.us, alpha = 0, col = 'grey') +
  geom_sf(data = within.sf.list[[which(sp.names == 'BWWA')]], alpha = 0, col = 'black') + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA) + 
  theme_bw(base_size = 18) +
  labs(x = "Longitude", y = "Latitude", fill = "Trend\n(logit scale)", 
       title = '(D) Blue-winged Warbler') +
  theme(legend.position = c(0.87, 0.22), 
        legend.background = element_rect(fill = NA), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_text(size = 14),
	plot.title = element_text(size = 18),
        legend.text = element_text(size = 12))

ggarrange(w.trend.TUTI.plot, w.trend.EAWP.plot, w.trend.WOTH.plot, 
	  w.trend.BWWA.plot, nrow = 2, ncol = 2, common.legend = FALSE)
ggsave(device = 'pdf', filename = 'figures/Fig-4.pdf', height = 12, width = 11,
       units = 'in')

# Calculate species-specific figures for appendix -------------------------
# These are the species codes that we model in the data set. 
#  [1] "AMWO" "MIKI" "RSHA" "BWHA" "EASO" "YBCU" "BBCU" "RCWO" "RHWO" "RBWO"
# [11] "CWWI" "CHSW" "RTHU" "GCFL" "EAPH" "EAWP" "ACFL" "BLJA" "OROR" "BAOR"
# [21] "FISP" "BACS" "EATO" "NOCA" "RBGR" "INBU" "SCTA" "SUTA" "YTVI" "WEVI"
# [31] "PROW" "SWWA" "WEWA" "BWWA" "GWWA" "NOPA" "CERW" "YTWA" "PIWA" "PRAW"
# [41] "LOWA" "KEWA" "HOWA" "GRCA" "BRTH" "CARW" "BHNU" "TUTI" "CACH" "WOTH"
# [51] "EABL"

# Data frame to hold information on trends for each species for use later on 
trend.sp.summaries <- data.frame(sp = sp.names, 
				 prop.pos = rep(NA, N), 
				 prop.low.pos = rep(NA, N),
				 prop.high.pos = rep(NA, N))

# For loop that goes over each species to generate the figures shown in Appendix S2
for (i in 1:N) {
  print(i)
  curr.sp <- sp.names[i]
  load(paste("results/svcTPGOcc-", curr.sp, ".R", sep = ''))
  out.svc <- out
  svc.samples <- getSVCSamples(out.svc)
  # Determine the proportion of the spatially-varying trends that are positive
  # for each MCMC sample
  svc.props.samples <- apply(svc.samples[[2]], 1, function(a) mean(a > 0))
  # Calculate median and 95% credible interval for the proportion of positive
  # trends across the study region.
  svc.quants <- quantile(svc.props.samples, c(0.025, 0.5, 0.975))
  trend.sp.summaries$prop.pos[i] <- svc.quants[2] 
  trend.sp.summaries$prop.low.pos[i] <- svc.quants[1]
  trend.sp.summaries$prop.high.pos[i] <- svc.quants[3]
  load(paste('results/predict-svcTPGOcc-', curr.sp, sep = ''))
  indx <- unlist(c(st_contains(within.sf.list[[i]], coords.0.sf)))
  trend.quants <- w.trend.quants + quantile(out.svc$beta.samples[, 2], c(0.025, 0.5, 0.975))
  trend.quants[, -indx] <- NA
  
  plot.df <- data.frame(trend = trend.quants[2, ],
		        ci.width = trend.quants[3, ] - trend.quants[1, ],
		        x = coords.0[, 1],
		        y = coords.0[, 2])
  pred.stars <- st_as_stars(plot.df, dims = c('x', 'y'))

  w.trend.plot <- ggplot() +
    geom_stars(data = pred.stars, aes(x = x, y = y, fill = trend), interpolate = TRUE) +
    geom_sf(data = east.us, alpha = 0, col = 'grey') +
    geom_sf(data = within.sf.list[[i]], alpha = 0, col = 'black') + 
    # geom_sf(data = bcrs.ne.grouped, alpha = 0) +
    scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
    	               na.value = NA) + 
    theme_bw(base_size = 18) +
    labs(x = "Longitude", y = "Latitude", fill = "Trend\n(logit scale)", 
         title = paste("(A) ", curr.sp, " Mean", sep = '')) +
    theme(legend.position = c(0.87, 0.22), 
          legend.background = element_rect(fill = NA), 
          axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.title = element_text(size = 14),
          plot.title = element_text(size = 18),
          legend.text = element_text(size = 12))
  
  w.ci.plot <- ggplot() +
    geom_stars(data = pred.stars, aes(x = x, y = y, fill = ci.width), interpolate = TRUE) +
    geom_sf(data = east.us, alpha = 0, col = 'grey') +
    geom_sf(data = within.sf.list[[i]], alpha = 0, col = 'black') + 
    # geom_sf(data = bcrs.ne.grouped, alpha = 0) +
    scale_fill_gradientn(colors = ocean.tempo(1000), na.value = NA,
                           guide = guide_colourbar(title.position="top", reverse = FALSE)) +
    # scale_fill_viridis_c(na.value = NA) +
    theme_bw(base_size = 18) +
    labs(x = "Longitude", y = "Latitude", fill = "CI Width\n(logit scale)", 
         title = paste("(B) ", curr.sp, " 95% CI Width", sep = '')) +
    theme(legend.position = c(0.87, 0.22), 
          legend.background = element_rect(fill = NA), 
          axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.title = element_text(size = 14),
          plot.title = element_text(size = 18),
          legend.text = element_text(size = 12))
  
  plot(ggarrange(w.trend.plot, w.ci.plot, ncol = 2))
  ggsave(device = 'pdf', filename = paste('figures/species-figs/', curr.sp, '-trend.pdf', sep = ''),
	 height = 6, width = 11, units = 'in')
}

# Summarize species-specific trends ---------------------------------------
# Generate a good order to make the summary plot display in order of species
# with predominately positive trends to species with predominately negative trends.
plot.order <- sp.names[order(trend.sp.summaries$prop.pos)]

# Generate Figure 3 in manuscript. 
trend.sp.summaries %>%
  mutate(sp = factor(sp, levels = plot.order, ordered = TRUE)) %>%
    ggplot(aes(x = prop.pos, y = sp, fill = prop.pos)) +
    geom_vline(xintercept = 0.5, lty = 2) +
    geom_segment(aes(x = prop.low.pos, y = sp, xend = prop.high.pos, yend = sp), 
		 lineend = 'butt', linewidth = 1, col = 'lightgray') + 
    geom_point(size = 4, pch = 21) +
    scale_fill_gradient2(midpoint = 0.5, low = '#B2182B', mid = 'white', high = '#2166AC',
    	               na.value = NA) +
    scale_color_gradient2(midpoint = 0.5, low = '#B2182B', mid = 'white', high = '#2166AC',
    	               na.value = NA) +
    scale_x_continuous(limits = c(0, 1)) +
    theme_classic(base_size = 18) +
    guides(fill = 'none', color = 'none') +
    labs(x = 'Proportion of sites with positive trend',
	 y = 'Species', fill = '')
  ggsave(file = 'figures/Fig-3.pdf', width = 7, height = 10)

# Compare select species with constant and BCR models ---------------------
# Gray Catbird ------------------------
curr.sp <- "GRCA"
# SVC model ---------------------------------------------------------------
load(paste('results/predict-svcTPGOcc-', curr.sp, sep = ''))
load(paste("results/svcTPGOcc-", curr.sp, ".R", sep = ''))
out.GRCA <- out
trend.GRCA.quants <- w.trend.quants + quantile(out.GRCA$beta.samples[, 2], c(0.025, 0.5, 0.975))
# Only plot the predicted trends within the range of GRCA 
indx <- unlist(c(st_contains(within.sf.list[[which(sp.names == 'GRCA')]], coords.0.sf)))
trend.GRCA.quants[, -indx] <- NA

# Constant model
load(paste("results/stPGOcc-", curr.sp, ".R", sep = ''))
out.GRCA.constant <- out
trend.GRCA.constant <- rep(median(out$beta.samples[, 2]), nrow(coords.0))
trend.GRCA.constant[-indx] <- NA

# BCR model
load(paste('results/stPGOcc-bcr-', curr.sp, '.R', sep = ''))
out.GRCA.bcr <- out
trend.GRCA.bcr.ests <- apply(out$beta.samples[, -1], 2, median)
x.names <- dimnames(out.GRCA.bcr$X)[[3]]
bcr.est <- sort(parse_number(x.names[-1]))
trend.GRCA.bcr <- rep(NA, nrow(coords.0))
for (i in 1:nrow(coords.0)) {
  if (length(which(bcr.est == bcr.pred.factor[i])) == 1) {
    trend.GRCA.bcr[i] <- trend.GRCA.bcr.ests[which(bcr.est == bcr.pred.factor[i])]
  }
}
trend.GRCA.bcr[-indx] <- NA


plot.df <- data.frame(trend = trend.GRCA.quants[2, ],
      	        ci.width = trend.GRCA.quants[3, ] - trend.GRCA.quants[1, ],
		trend.constant = trend.GRCA.constant,
		trend.bcr = trend.GRCA.bcr,
      	        x = coords.0[, 1],
      	        y = coords.0[, 2])
pred.stars <- st_as_stars(plot.df, dims = c('x', 'y'))
min.trend <- min(plot.df$trend, plot.df$trend.constant, plot.df$trend.bcr, na.rm = TRUE)
max.trend <- max(plot.df$trend, plot.df$trend.constant, plot.df$trend.bcr, na.rm = TRUE)

svc.GRCA.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = trend), interpolate = TRUE) +
  geom_sf(data = east.us, alpha = 0, col = 'grey') +
  geom_sf(data = within.sf.list[[which(sp.names == 'GRCA')]], alpha = 0, col = 'black') + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA, limits = c(min.trend, max.trend)) + 
  theme_bw(base_size = 18) +
  labs(x = "Longitude", y = "Latitude", fill = "Trend\n(logit scale)", 
       title = '(C) Gray Catbird SVC') +
  theme(legend.position = c(0.87, 0.22), 
        legend.background = element_rect(fill = NA), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_text(size = 14),
	plot.title = element_text(size = 18),
        legend.text = element_text(size = 12))

constant.GRCA.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = trend.constant), interpolate = TRUE) +
  geom_sf(data = east.us, alpha = 0, col = 'grey') +
  geom_sf(data = within.sf.list[[which(sp.names == 'GRCA')]], alpha = 0, col = 'black') + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA, 
		       limits = c(min.trend, max.trend)) + 
  theme_bw(base_size = 18) +
  labs(x = "Longitude", y = "Latitude", fill = "Trend\n(logit scale)", 
       title = '(A) Gray Catbird Constant') +
  theme(legend.position = c(0.87, 0.22), 
        legend.background = element_rect(fill = NA), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_text(size = 14),
	plot.title = element_text(size = 18),
        legend.text = element_text(size = 12))
bcr.GRCA.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = trend.bcr), interpolate = TRUE) +
  geom_sf(data = bcrs.ne.grouped, alpha = 0, col = 'grey') +
  geom_sf(data = within.sf.list[[which(sp.names == 'GRCA')]], alpha = 0, col = 'black') + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA, 
		       limits = c(min.trend, max.trend)) + 
  theme_bw(base_size = 18) +
  labs(x = "Longitude", y = "Latitude", fill = "Trend\n(logit scale)", 
       title = '(B) Gray Catbird BCR') +
  theme(legend.position = c(0.87, 0.22), 
        legend.background = element_rect(fill = NA), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_text(size = 14),
	plot.title = element_text(size = 18),
        legend.text = element_text(size = 12))

# Eastern Phoebe
curr.sp <- "EAPH"
# SVC model ---------------------------------------------------------------
load(paste('results/predict-svcTPGOcc-', curr.sp, sep = ''))
load(paste("results/svcTPGOcc-", curr.sp, ".R", sep = ''))
out.EAPH <- out
trend.EAPH.quants <- w.trend.quants + quantile(out.EAPH$beta.samples[, 2], c(0.025, 0.5, 0.975))
# Only plot the predicted trends within the range of EAPH 
indx <- unlist(c(st_contains(within.sf.list[[which(sp.names == 'EAPH')]], coords.0.sf)))
trend.EAPH.quants[, -indx] <- NA

# Constant model
load(paste("results/stPGOcc-", curr.sp, ".R", sep = ''))
out.EAPH.constant <- out
trend.EAPH.constant <- rep(median(out$beta.samples[, 2]), nrow(coords.0))
trend.EAPH.constant[-indx] <- NA

# BCR model
load(paste('results/stPGOcc-bcr-', curr.sp, '.R', sep = ''))
out.EAPH.bcr <- out
trend.EAPH.bcr.ests <- apply(out$beta.samples[, -1], 2, median)
x.names <- dimnames(out.EAPH.bcr$X)[[3]]
bcr.est <- sort(parse_number(x.names[-1]))
trend.EAPH.bcr <- rep(NA, nrow(coords.0))
for (i in 1:nrow(coords.0)) {
  if (length(which(bcr.est == bcr.pred.factor[i])) == 1) {
    trend.EAPH.bcr[i] <- trend.EAPH.bcr.ests[which(bcr.est == bcr.pred.factor[i])]
  }
}
trend.EAPH.bcr[-indx] <- NA



plot.df <- data.frame(trend = trend.EAPH.quants[2, ],
      	        ci.width = trend.EAPH.quants[3, ] - trend.EAPH.quants[1, ],
		trend.constant = trend.EAPH.constant,
		trend.bcr = trend.EAPH.bcr,
      	        x = coords.0[, 1],
      	        y = coords.0[, 2])
pred.stars <- st_as_stars(plot.df, dims = c('x', 'y'))
min.trend <- min(plot.df$trend, plot.df$trend.constant, plot.df$trend.bcr, na.rm = TRUE)
max.trend <- max(plot.df$trend, plot.df$trend.constant, plot.df$trend.bcr, na.rm = TRUE)

svc.EAPH.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = trend), interpolate = TRUE) +
  geom_sf(data = east.us, alpha = 0, col = 'grey') +
  geom_sf(data = within.sf.list[[which(sp.names == 'EAPH')]], alpha = 0, col = 'black') + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA, limits = c(min.trend, max.trend)) + 
  theme_bw(base_size = 18) +
  labs(x = "Longitude", y = "Latitude", fill = "Trend\n(logit scale)", 
       title = '(F) Eastern Phoebe SVC') +
  theme(legend.position = c(0.87, 0.22), 
        legend.background = element_rect(fill = NA), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_text(size = 14),
	plot.title = element_text(size = 18),
        legend.text = element_text(size = 12))

constant.EAPH.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = trend.constant), interpolate = TRUE) +
  geom_sf(data = east.us, alpha = 0, col = 'grey') +
  geom_sf(data = within.sf.list[[which(sp.names == 'EAPH')]], alpha = 0, col = 'black') + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA, 
		       limits = c(min.trend, max.trend)) + 
  theme_bw(base_size = 18) +
  labs(x = "Longitude", y = "Latitude", fill = "Trend\n(logit scale)", 
       title = '(D) Eastern Phoebe Constant') +
  theme(legend.position = c(0.87, 0.22), 
        legend.background = element_rect(fill = NA), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_text(size = 14),
	plot.title = element_text(size = 18),
        legend.text = element_text(size = 12))
bcr.EAPH.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = trend.bcr), interpolate = TRUE) +
  geom_sf(data = bcrs.ne.grouped, alpha = 0, col = 'grey') +
  geom_sf(data = within.sf.list[[which(sp.names == 'EAPH')]], alpha = 0, col = 'black') + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA, 
		       limits = c(min.trend, max.trend)) + 
  theme_bw(base_size = 18) +
  labs(x = "Longitude", y = "Latitude", fill = "Trend\n(logit scale)", 
       title = '(E) Eastern Phoebe BCR') +
  theme(legend.position = c(0.87, 0.22), 
        legend.background = element_rect(fill = NA), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_text(size = 14),
	plot.title = element_text(size = 18),
        legend.text = element_text(size = 12))

# Scarlett Tanager
curr.sp <- "SCTA"
# SVC model ---------------------------------------------------------------
load(paste('results/predict-svcTPGOcc-', curr.sp, sep = ''))
load(paste("results/svcTPGOcc-", curr.sp, ".R", sep = ''))
out.SCTA <- out
trend.SCTA.quants <- w.trend.quants + quantile(out.SCTA$beta.samples[, 2], c(0.025, 0.5, 0.975))
# Only plot the predicted trends within the range of SCTA 
indx <- unlist(c(st_contains(within.sf.list[[which(sp.names == 'SCTA')]], coords.0.sf)))
trend.SCTA.quants[, -indx] <- NA

# Constant model
load(paste("results/stPGOcc-", curr.sp, ".R", sep = ''))
out.SCTA.constant <- out
trend.SCTA.constant <- rep(median(out$beta.samples[, 2]), nrow(coords.0))
trend.SCTA.constant[-indx] <- NA

# BCR model
load(paste('results/stPGOcc-bcr-', curr.sp, '.R', sep = ''))
out.SCTA.bcr <- out
trend.SCTA.bcr.ests <- apply(out$beta.samples[, -1], 2, median)
x.names <- dimnames(out.SCTA.bcr$X)[[3]]
bcr.est <- sort(parse_number(x.names[-1]))
trend.SCTA.bcr <- rep(NA, nrow(coords.0))
for (i in 1:nrow(coords.0)) {
  if (length(which(bcr.est == bcr.pred.factor[i])) == 1) {
    trend.SCTA.bcr[i] <- trend.SCTA.bcr.ests[which(bcr.est == bcr.pred.factor[i])]
  }
}
trend.SCTA.bcr[-indx] <- NA



plot.df <- data.frame(trend = trend.SCTA.quants[2, ],
      	        ci.width = trend.SCTA.quants[3, ] - trend.SCTA.quants[1, ],
		trend.constant = trend.SCTA.constant,
		trend.bcr = trend.SCTA.bcr,
      	        x = coords.0[, 1],
      	        y = coords.0[, 2])
pred.stars <- st_as_stars(plot.df, dims = c('x', 'y'))
min.trend <- min(plot.df$trend, plot.df$trend.constant, plot.df$trend.bcr, na.rm = TRUE)
max.trend <- max(plot.df$trend, plot.df$trend.constant, plot.df$trend.bcr, na.rm = TRUE)

svc.SCTA.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = trend), interpolate = TRUE) +
  geom_sf(data = east.us, alpha = 0, col = 'grey') +
  geom_sf(data = within.sf.list[[which(sp.names == 'SCTA')]], alpha = 0, col = 'black') + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA, limits = c(min.trend, max.trend)) + 
  theme_bw(base_size = 18) +
  labs(x = "Longitude", y = "Latitude", fill = "Trend\n(logit scale)", 
       title = '(I) Scarlet Tanager SVC') +
  theme(legend.position = c(0.87, 0.22), 
        legend.background = element_rect(fill = NA), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_text(size = 14),
	plot.title = element_text(size = 18),
        legend.text = element_text(size = 12))

constant.SCTA.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = trend.constant), interpolate = TRUE) +
  geom_sf(data = east.us, alpha = 0, col = 'grey') +
  geom_sf(data = within.sf.list[[which(sp.names == 'SCTA')]], alpha = 0, col = 'black') + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA, 
		       limits = c(min.trend, max.trend)) + 
  theme_bw(base_size = 18) +
  labs(x = "Longitude", y = "Latitude", fill = "Trend\n(logit scale)", 
       title = '(G) Scarlet Tanager Constant') +
  theme(legend.position = c(0.87, 0.22), 
        legend.background = element_rect(fill = NA), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_text(size = 14),
	plot.title = element_text(size = 18),
        legend.text = element_text(size = 12))
bcr.SCTA.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = trend.bcr), interpolate = TRUE) +
  geom_sf(data = bcrs.ne.grouped, alpha = 0, col = 'grey') +
  geom_sf(data = within.sf.list[[which(sp.names == 'SCTA')]], alpha = 0, col = 'black') + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA, 
		       limits = c(min.trend, max.trend)) + 
  theme_bw(base_size = 18) +
  labs(x = "Longitude", y = "Latitude", fill = "Trend\n(logit scale)", 
       title = '(H) Scarlet Tanager BCR') +
  theme(legend.position = c(0.87, 0.22), 
        legend.background = element_rect(fill = NA), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_text(size = 14),
	plot.title = element_text(size = 18),
        legend.text = element_text(size = 12))

ggarrange(constant.GRCA.plot, bcr.GRCA.plot, svc.GRCA.plot, 
          constant.EAPH.plot, bcr.EAPH.plot, svc.EAPH.plot,
	  constant.SCTA.plot, bcr.SCTA.plot, svc.SCTA.plot,
	  ncol = 3, nrow = 3, common.legend = TRUE, legend = 'bottom')
ggsave(device = 'pdf', filename = 'figures/Fig-5.pdf', 
       height = 15, width = 13, units = 'in')
