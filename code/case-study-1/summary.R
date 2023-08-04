# summary.R: this script summarizes Case Study 1 results in which we 
#            assess spatially-varying trends in 51 forest birds across the 
#            eastern US. 
# Author: Jeffrey W. Doser
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
load("data/case-study-1/spOcc-bbs-data.rda")
load("data/case-study-1/bird-life-processed.rda")
sp.names <- dimnames(data.list$y)[[1]]
sp.names <- sp.names[percent.area > 0.5]
within.sf.list <- within.sf.list[percent.area > 0.5]
N <- length(sp.names)

# Load data list used to fit model
load("data/case-study-1/final-spOccupancy-data.rda")

# Average number of years each route is surveyed. 
n.surveys.per.route <- apply(data.list$y[1, , , 1], 1, function(a) sum(!is.na(a)))
mean(n.surveys.per.route)
sd(n.surveys.per.route)

# Assess AUC predictive performance for 2021 -----------------------------
load("results/case-study-1/auc-results.rda")
# This is Table 3 in Supp Info S2.
auc.mean.df %>%
  arrange(sp)

# Assess model fit with WAIC ----------------------------------------------
load("results/case-study-1/waic-results.rda")

# This is Table 2 in Supp Info S2
waic.df$sp <- sp.names
waic.df$constant <- waic.df$constant - waic.df$svc
waic.df$tmax <- waic.df$tmax - waic.df$svc
waic.df$bcr <- waic.df$bcr - waic.df$svc
waic.df$svc <- waic.df$svc - waic.df$svc
waic.df %>% arrange(sp)

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
load('data/case-study-1/pred-coordinates.rda')
coords.0.sf <- st_as_sf(as.data.frame(coords.0),
		      coords = c("X", "Y"),
		      crs = st_crs(coords.sf))

# Generate summary figure of species trends -------------------------------
trend.summaries <- data.frame(sp = sp.names, 
			      prop.lowest = NA,
			      prop.low = NA, 
			      prop.med = NA,
			      prop.high = NA,
			      prop.highest = NA)
ordering.var <- rep(NA, N)
trend.meds <- matrix(NA, N, nrow(coords.0.sf))
# Matrix to hold species-specific trends, only at locations within that species range
for (i in 1:N) {
  curr.sp <- sp.names[i]
  load(paste('results/case-study-1/predict-svcTPGOcc-', curr.sp, '.rda', sep = ''))
  indx <- unlist(c(st_contains(within.sf.list[[i]], coords.0.sf)))
  trend.prob.pos[-indx] <- NA
  trend.quants[, -indx] <- NA
  trend.meds[i, ] <- trend.quants[2, ]
  ordering.var[i] <- mean(trend.prob.pos, na.rm = TRUE)
  trend.summaries$prop.lowest[i] <- mean(trend.prob.pos <= .2, na.rm = TRUE)
  trend.summaries$prop.low[i] <- mean(trend.prob.pos > .2 & trend.prob.pos <= .4, na.rm = TRUE)
  trend.summaries$prop.med[i] <- mean(trend.prob.pos > .4 & trend.prob.pos <= .6, na.rm = TRUE)
  trend.summaries$prop.high[i] <- mean(trend.prob.pos > .6 & trend.prob.pos <= .8, na.rm = TRUE)
  trend.summaries$prop.highest[i] <- mean(trend.prob.pos > .8, na.rm = TRUE)
}

plot.df <- data.frame(prop = c(trend.summaries$prop.lowest, trend.summaries$prop.low, 
			       trend.summaries$prop.med, trend.summaries$prop.high,
			       trend.summaries$prop.highest),
		      type = factor(rep(c('Strong Negative', 'Moderate Negative', 'No Support',
				   'Moderate Positive', 'Strong Positive'), each = N),
				    levels = c('Strong Negative', 'Moderate Negative',
					       'No Support', 'Moderate Positive',
					       'Strong Positive')),
		      species = rep(sp.names, times = 5),
                      mean.val = rep(ordering.var, times = 5))
plot.df <- plot.df %>%
  arrange(desc(mean.val))
plot.df$species <- factor(plot.df$species, levels = unique(plot.df$species), order = TRUE)

# Generate Figure 2
trend.plot <- ggplot(plot.df, aes(x = species, y = prop, fill = type)) +
  geom_bar(stat = 'identity', width = 1, color = 'grey') +
  theme_bw(base_size = 15) +
  # scale_fill_viridis_d() +
  scale_fill_brewer(palette = 'RdBu') +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks.x = element_blank(),
	legend.title = element_text(size = 10),
	plot.title = element_text(size = 15),
	legend.text = element_text(size = 10),
	text = element_text(family  = 'LM Roman 10'),
 	legend.position = 'top') +
  labs(x = 'Species', y = 'Proportion of Locations',
       fill = '')
# ggsave(plot = trend.plot, file = 'figures/case-study-1/Figure-2.png', 
#        units = 'in', device = 'png', width = 12, height = 5)

# Generate summary of the interaction -------------------------------------
interaction.df <- data.frame(sp = sp.names, 
			     med = NA, 
			     low = NA, 
			     high = NA, 
			     prob.neg = NA)
for (i in 1:N) {
  print(i)
  curr.sp <- sp.names[i]
  load(paste("results/case-study-1/stPGOcc-tmax-", curr.sp, ".rda", sep = ''))
  beta.tmax.int <- out$beta.samples[, ncol(out$beta.samples)]
  interaction.df$med[i] <- median(beta.tmax.int)
  interaction.df$low[i] <- quantile(beta.tmax.int, 0.025)
  interaction.df$high[i] <- quantile(beta.tmax.int, 0.975)
  interaction.df$prob.neg[i] <- mean(beta.tmax.int < 0)
}

plot.order <- sp.names[order(interaction.df$med)]

# Generate Supp Info S2 Figure 1.
interaction.plot <- interaction.df %>%
  mutate(sp = factor(sp, levels = plot.order, ordered = TRUE)) %>%
    ggplot(aes(x = med, y = sp, fill = prob.neg)) +
    geom_vline(xintercept = 0, lty = 2) +
    geom_segment(aes(x = low, y = sp, xend = high, yend = sp),
		 lineend = 'butt', linewidth = 1, col = 'lightgray') +
    geom_point(size = 4, pch = 21) +
    scale_fill_gradient2(midpoint = 0.5, high = '#B2182B', mid = 'white', low = '#2166AC',
    	               na.value = NA) +
    theme_classic(base_size = 17) +
    labs(x = 'Interaction Effect Size',
	 y = 'Species', fill = 'P(effect < 0)') +
    theme(text = element_text(family = 'LM Roman 10'), 
          legend.position = c(0.80, 0.22), 
          legend.background = element_rect(fill = NA))
  ggsave(plot = interaction.plot, 
	 file = 'figures/case-study-1/Supp-Info-S2-Figure-1.png', width = 7, height = 9)

# Calculate species-specific figures for appendix -------------------------
# This code generates Supp Info S2 Figures 2-52.
# Data frame to hold information on trends for each species for use later on 
trend.sp.summaries <- data.frame(sp = sp.names, 
				 prop.pos = rep(NA, N), 
				 prop.low.pos = rep(NA, N),
				 prop.high.pos = rep(NA, N))

# For loop that goes over each species to generate the figures shown in Appendix S2
for (i in 1:N) {
  print(i)
  curr.sp <- sp.names[i]
  load(paste("results/case-study-1/svcTPGOcc-", curr.sp, ".rda", sep = ''))
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
  load(paste('results/case-study-1/predict-svcTPGOcc-', curr.sp, '.rda', sep = ''))
  indx <- unlist(c(st_contains(within.sf.list[[i]], coords.0.sf)))
  trend.quants[, -indx] <- NA
  trend.prob.pos[-indx] <- NA
  psi.quants[, -indx, ] <- NA
  
  plot.df <- data.frame(trend = trend.quants[2, ],
		        ci.width = trend.quants[3, ] - trend.quants[1, ],
			psi.med = apply(psi.quants[2, , ], 1, median),
			prob.pos = trend.prob.pos,
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
    theme_bw(base_size = 15) +
    labs(x = "Longitude", y = "Latitude", fill = "Trend\n(logit scale)", 
         title = paste("(A) ", curr.sp, " Mean", sep = '')) +
    theme(legend.position = c(0.89, 0.22), 
          legend.background = element_rect(fill = NA), 
          axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.title = element_text(size = 10),
          text = element_text(family="LM Roman 10"),
          plot.title = element_text(size = 15),
          legend.text = element_text(size = 10))
  
  w.p.trend.plot <- ggplot() +
    geom_stars(data = pred.stars, aes(x = x, y = y, fill = prob.pos)) +
    geom_sf(data = east.us, alpha = 0, col = 'grey') +
    geom_sf(data = within.sf.list[[i]], alpha = 0, col = 'black') + 
    # geom_sf(data = bcrs.ne.grouped, alpha = 0) +
    scale_fill_steps2(midpoint = 0.5, low = '#B2182B', mid = 'white', high = '#2166AC',
    	               na.value = NA, limits = c(0, 1), n.breaks = 6) +
    # scale_fill_viridis_c(na.value = NA) +
    theme_bw(base_size = 15) +
    labs(x = "Longitude", y = "Latitude", fill = "P(trend > 0)", 
         title = paste("(B) ", curr.sp, " P(trend > 0)", sep = '')) +
    theme(legend.position = c(0.89, 0.22), 
          legend.background = element_rect(fill = NA), 
          axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.title = element_text(size = 10),
          text = element_text(family="LM Roman 10"),
          plot.title = element_text(size = 15),
          legend.text = element_text(size = 10))

  psi.plot <- ggplot() +
    geom_stars(data = pred.stars, aes(x = x, y = y, fill = psi.med), interpolate = TRUE) +	
    geom_sf(data = east.us, alpha = 0, col = 'grey') +
    geom_sf(data = within.sf.list[[i]], alpha = 0, col = 'black') + 
    scale_fill_gradientn("", colors = ocean.tempo(5), limits = c(0, 1),
                             guide = guide_colourbar(title.position="top", reverse = FALSE), 
        		   na.value = NA) +
    theme_bw(base_size = 15) +
    labs(x = "Longitude", y = "Latitude", fill = "", 
         title = paste("(C) ", curr.sp, " Occurrence Probability", sep = '')) +
    theme(legend.position = c(0.89, 0.22), 
          legend.background = element_rect(fill = NA), 
          text = element_text(family = 'LM Roman 10'),
          axis.text.x = element_text(angle = 45, hjust = 1), 
          legend.title = element_text(size = 12),
          plot.title = element_text(size = 15),
          legend.text = element_text(size = 10))
    
    trend.plot <- ggarrange(w.trend.plot, w.p.trend.plot, psi.plot, ncol = 3)
    ggsave(plot = trend.plot, device = 'png', 
  	 filename = paste('figures/case-study-1/species-figs/', curr.sp, '-trend.png', sep = ''),
  	 height = 6, width = 14, units = 'in', bg = 'white')
}

# Compare estimated trends across models ----------------------------------
# Gray Catbird ------------------------
curr.sp <- "GRCA"
# SVC model ---------------------------------------------------------------
load(paste('results/case-study-1/predict-svcTPGOcc-', curr.sp, '.rda', sep = ''))
trend.GRCA.quants <- trend.quants 
# Only plot the predicted trends within the range of GRCA 
indx <- unlist(c(st_contains(within.sf.list[[which(sp.names == 'GRCA')]], coords.0.sf)))
trend.GRCA.quants[, -indx] <- NA

# Constant model ----------------------------------------------------------
load(paste("results/case-study-1/stPGOcc-", curr.sp, ".rda", sep = ''))
out.GRCA.constant <- out
trend.GRCA.constant <- rep(median(out$beta.samples[, 2]), nrow(coords.0))
trend.GRCA.constant[-indx] <- NA

# BCR model ---------------------------------------------------------------
load(paste('results/case-study-1/stPGOcc-bcr-', curr.sp, '.rda', sep = ''))
out.GRCA.bcr <- out
# NOTE: hardcoded
bad.cols <- c(1, 2)
trend.GRCA.bcr.ests <- apply(out$beta.samples[, -bad.cols], 2, median)
x.names <- dimnames(out.GRCA.bcr$X)[[3]]
bcr.est <- sort(parse_number(x.names[-bad.cols]))
trend.GRCA.bcr <- rep(NA, nrow(coords.0))
for (i in 1:nrow(coords.0)) {
  if (length(which(bcr.est == bcr.pred.factor[i])) == 1) {
    trend.GRCA.bcr[i] <- trend.GRCA.bcr.ests[which(bcr.est == bcr.pred.factor[i])]
  }
}
trend.GRCA.bcr[-indx] <- NA

# TMAX model --------------------------------------------------------------
load(paste('results/case-study-1/stPGOcc-tmax-', curr.sp, '.rda', sep = ''))
out.GRCA.tmax <- out
fit.tmax <- data.list$occ.covs$tmax[data.list$range.ind[which(sp.names == 'GRCA'), ] == 1]
tmax.0.scaled <- (tmax.0 - mean(fit.tmax)) / sd(fit.tmax)
# NOTE: hardcoded
bad.indx <- c(1, 3)
beta.tmax.ests <- apply(out.GRCA.tmax$beta.samples[, -bad.indx], 2, median)
trend.GRCA.tmax <- beta.tmax.ests[1] + beta.tmax.ests[2] * tmax.0.scaled
trend.GRCA.tmax[-indx] <- NA


plot.df <- data.frame(trend = trend.GRCA.quants[2, ],
      	        ci.width = trend.GRCA.quants[3, ] - trend.GRCA.quants[1, ],
		trend.constant = trend.GRCA.constant,
		trend.bcr = trend.GRCA.bcr,
		trend.tmax = trend.GRCA.tmax, 
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
  theme_bw(base_size = 15) +
  labs(fill = "Trend\n(logit scale)", 
       title = '(D) GRCA SVC') +
  theme(legend.position = c(0.87, 0.22), 
        legend.background = element_rect(fill = NA), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(family="LM Roman 10"),
        legend.title = element_text(size = 14),
	plot.title = element_text(size = 15),
        legend.text = element_text(size = 12))

tmax.GRCA.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = trend.tmax), interpolate = TRUE) +
  geom_sf(data = east.us, alpha = 0, col = 'grey') +
  geom_sf(data = within.sf.list[[which(sp.names == 'GRCA')]], alpha = 0, col = 'black') + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA, limits = c(min.trend, max.trend)) + 
  theme_bw(base_size = 15) +
  labs(fill = "Trend\n(logit scale)", 
       title = '(C) GRCA TMAX') +
  theme(legend.position = c(0.87, 0.22), 
        legend.background = element_rect(fill = NA), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(family="LM Roman 10"),
        legend.title = element_text(size = 14),
	plot.title = element_text(size = 15),
        legend.text = element_text(size = 12))

constant.GRCA.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = trend.constant), interpolate = TRUE) +
  geom_sf(data = east.us, alpha = 0, col = 'grey') +
  geom_sf(data = within.sf.list[[which(sp.names == 'GRCA')]], alpha = 0, col = 'black') + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA, 
		       limits = c(min.trend, max.trend)) + 
  theme_bw(base_size = 15) +
  labs(fill = "Trend\n(logit scale)", 
       title = '(A) GRCA Linear') +
  theme(legend.position = c(0.87, 0.22), 
        legend.background = element_rect(fill = NA), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_text(size = 14),
        text = element_text(family="LM Roman 10"),
	plot.title = element_text(size = 15),
        legend.text = element_text(size = 12))
bcr.GRCA.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = trend.bcr), interpolate = TRUE) +
  geom_sf(data = bcrs.ne.grouped, alpha = 0, col = 'grey') +
  geom_sf(data = within.sf.list[[which(sp.names == 'GRCA')]], alpha = 0, col = 'black') + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA, 
		       limits = c(min.trend, max.trend)) + 
  theme_bw(base_size = 15) +
  labs(fill = "Trend\n(logit scale)", 
       title = '(B) GRCA Strata') +
  theme(legend.position = c(0.87, 0.22), 
        legend.background = element_rect(fill = NA), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(family="LM Roman 10"),
        legend.title = element_text(size = 14),
	plot.title = element_text(size = 15),
        legend.text = element_text(size = 12))

# Eastern Phoebe ----------------------
curr.sp <- "EAPH"
# SVC model ---------------------------------------------------------------
load(paste('results/case-study-1/predict-svcTPGOcc-', curr.sp, '.rda', sep = ''))
trend.EAPH.quants <- trend.quants 
# Only plot the predicted trends within the range of EAPH 
indx <- unlist(c(st_contains(within.sf.list[[which(sp.names == 'EAPH')]], coords.0.sf)))
trend.EAPH.quants[, -indx] <- NA

# Constant model ----------------------------------------------------------
load(paste("results/case-study-1/stPGOcc-", curr.sp, ".rda", sep = ''))
out.EAPH.constant <- out
trend.EAPH.constant <- rep(median(out$beta.samples[, 2]), nrow(coords.0))
trend.EAPH.constant[-indx] <- NA

# BCR model ---------------------------------------------------------------
load(paste('results/case-study-1/stPGOcc-bcr-', curr.sp, '.rda', sep = ''))
out.EAPH.bcr <- out
# NOTE: hardcoded
bad.cols <- c(1, 2)
trend.EAPH.bcr.ests <- apply(out$beta.samples[, -bad.cols], 2, median)
x.names <- dimnames(out.EAPH.bcr$X)[[3]]
bcr.est <- sort(parse_number(x.names[-bad.cols]))
trend.EAPH.bcr <- rep(NA, nrow(coords.0))
for (i in 1:nrow(coords.0)) {
  if (length(which(bcr.est == bcr.pred.factor[i])) == 1) {
    trend.EAPH.bcr[i] <- trend.EAPH.bcr.ests[which(bcr.est == bcr.pred.factor[i])]
  }
}
trend.EAPH.bcr[-indx] <- NA

# TMAX model --------------------------------------------------------------
load(paste('results/case-study-1/stPGOcc-tmax-', curr.sp, '.rda', sep = ''))
out.EAPH.tmax <- out
fit.tmax <- data.list$occ.covs$tmax[data.list$range.ind[which(sp.names == 'EAPH'), ] == 1]
tmax.0.scaled <- (tmax.0 - mean(fit.tmax)) / sd(fit.tmax)
# NOTE: hardcoded
bad.indx <- c(1, 3)
beta.tmax.ests <- apply(out.EAPH.tmax$beta.samples[, -bad.indx], 2, median)
trend.EAPH.tmax <- beta.tmax.ests[1] + beta.tmax.ests[2] * tmax.0.scaled
trend.EAPH.tmax[-indx] <- NA


plot.df <- data.frame(trend = trend.EAPH.quants[2, ],
      	        ci.width = trend.EAPH.quants[3, ] - trend.EAPH.quants[1, ],
		trend.constant = trend.EAPH.constant,
		trend.bcr = trend.EAPH.bcr,
		trend.tmax = trend.EAPH.tmax, 
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
  theme_bw(base_size = 15) +
  labs(fill = "Trend\n(logit scale)", 
       title = '(H) EAPH SVC') +
  theme(legend.position = c(0.87, 0.22), 
        legend.background = element_rect(fill = NA), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(family="LM Roman 10"),
        legend.title = element_text(size = 14),
	plot.title = element_text(size = 15),
        legend.text = element_text(size = 12))

tmax.EAPH.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = trend.tmax), interpolate = TRUE) +
  geom_sf(data = east.us, alpha = 0, col = 'grey') +
  geom_sf(data = within.sf.list[[which(sp.names == 'EAPH')]], alpha = 0, col = 'black') + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA, limits = c(min.trend, max.trend)) + 
  theme_bw(base_size = 15) +
  labs(fill = "Trend\n(logit scale)", 
       title = '(G) EAPH TMAX') +
  theme(legend.position = c(0.87, 0.22), 
        legend.background = element_rect(fill = NA), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(family="LM Roman 10"),
        legend.title = element_text(size = 14),
	plot.title = element_text(size = 15),
        legend.text = element_text(size = 12))

constant.EAPH.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = trend.constant), interpolate = TRUE) +
  geom_sf(data = east.us, alpha = 0, col = 'grey') +
  geom_sf(data = within.sf.list[[which(sp.names == 'EAPH')]], alpha = 0, col = 'black') + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA, 
		       limits = c(min.trend, max.trend)) + 
  theme_bw(base_size = 15) +
  labs(fill = "Trend\n(logit scale)", 
       title = '(E) EAPH Linear') +
  theme(legend.position = c(0.87, 0.22), 
        legend.background = element_rect(fill = NA), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_text(size = 14),
        text = element_text(family="LM Roman 10"),
	plot.title = element_text(size = 15),
        legend.text = element_text(size = 12))
bcr.EAPH.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = trend.bcr), interpolate = TRUE) +
  geom_sf(data = bcrs.ne.grouped, alpha = 0, col = 'grey') +
  geom_sf(data = within.sf.list[[which(sp.names == 'EAPH')]], alpha = 0, col = 'black') + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA, 
		       limits = c(min.trend, max.trend)) + 
  theme_bw(base_size = 15) +
  labs(fill = "Trend\n(logit scale)", 
       title = '(F) EAPH Strata') +
  theme(legend.position = c(0.87, 0.22), 
        legend.background = element_rect(fill = NA), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(family="LM Roman 10"),
        legend.title = element_text(size = 14),
	plot.title = element_text(size = 15),
        legend.text = element_text(size = 12))

# Wood Thrush
curr.sp <- "WOTH"
# SVC model ---------------------------------------------------------------
load(paste('results/case-study-1/predict-svcTPGOcc-', curr.sp, '.rda', sep = ''))
trend.WOTH.quants <- trend.quants 
# Only plot the predicted trends within the range of WOTH 
indx <- unlist(c(st_contains(within.sf.list[[which(sp.names == 'WOTH')]], coords.0.sf)))
trend.WOTH.quants[, -indx] <- NA

# Constant model ----------------------------------------------------------
load(paste("results/case-study-1/stPGOcc-", curr.sp, ".rda", sep = ''))
out.WOTH.constant <- out
trend.WOTH.constant <- rep(median(out$beta.samples[, 2]), nrow(coords.0))
trend.WOTH.constant[-indx] <- NA

# BCR model ---------------------------------------------------------------
load(paste('results/case-study-1/stPGOcc-bcr-', curr.sp, '.rda', sep = ''))
out.WOTH.bcr <- out
# NOTE: hardcoded
bad.cols <- c(1, 2)
trend.WOTH.bcr.ests <- apply(out$beta.samples[, -bad.cols], 2, median)
x.names <- dimnames(out.WOTH.bcr$X)[[3]]
bcr.est <- sort(parse_number(x.names[-bad.cols]))
trend.WOTH.bcr <- rep(NA, nrow(coords.0))
for (i in 1:nrow(coords.0)) {
  if (length(which(bcr.est == bcr.pred.factor[i])) == 1) {
    trend.WOTH.bcr[i] <- trend.WOTH.bcr.ests[which(bcr.est == bcr.pred.factor[i])]
  }
}
trend.WOTH.bcr[-indx] <- NA

# TMAX model --------------------------------------------------------------
load(paste('results/case-study-1/stPGOcc-tmax-', curr.sp, '.rda', sep = ''))
out.WOTH.tmax <- out
fit.tmax <- data.list$occ.covs$tmax[data.list$range.ind[which(sp.names == 'WOTH'), ] == 1]
tmax.0.scaled <- (tmax.0 - mean(fit.tmax)) / sd(fit.tmax)
# NOTE: hardcoded
bad.indx <- c(1, 3)
beta.tmax.ests <- apply(out.WOTH.tmax$beta.samples[, -bad.indx], 2, median)
trend.WOTH.tmax <- beta.tmax.ests[1] + beta.tmax.ests[2] * tmax.0.scaled
trend.WOTH.tmax[-indx] <- NA


plot.df <- data.frame(trend = trend.WOTH.quants[2, ],
      	        ci.width = trend.WOTH.quants[3, ] - trend.WOTH.quants[1, ],
		trend.constant = trend.WOTH.constant,
		trend.bcr = trend.WOTH.bcr,
		trend.tmax = trend.WOTH.tmax, 
      	        x = coords.0[, 1],
      	        y = coords.0[, 2])
pred.stars <- st_as_stars(plot.df, dims = c('x', 'y'))
min.trend <- min(plot.df$trend, plot.df$trend.constant, plot.df$trend.bcr, na.rm = TRUE)
max.trend <- max(plot.df$trend, plot.df$trend.constant, plot.df$trend.bcr, na.rm = TRUE)

svc.WOTH.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = trend), interpolate = TRUE) +
  geom_sf(data = east.us, alpha = 0, col = 'grey') +
  geom_sf(data = within.sf.list[[which(sp.names == 'WOTH')]], alpha = 0, col = 'black') + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA, limits = c(min.trend, max.trend)) + 
  theme_bw(base_size = 15) +
  labs(fill = "Trend\n(logit scale)", 
       title = '(L) WOTH SVC') +
  theme(legend.position = c(0.87, 0.22), 
        legend.background = element_rect(fill = NA), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(family="LM Roman 10"),
        legend.title = element_text(size = 14),
	plot.title = element_text(size = 15),
        legend.text = element_text(size = 12))

tmax.WOTH.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = trend.tmax), interpolate = TRUE) +
  geom_sf(data = east.us, alpha = 0, col = 'grey') +
  geom_sf(data = within.sf.list[[which(sp.names == 'WOTH')]], alpha = 0, col = 'black') + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA, limits = c(min.trend, max.trend)) + 
  theme_bw(base_size = 15) +
  labs(fill = "Trend\n(logit scale)", 
       title = '(K) WOTH TMAX') +
  theme(legend.position = c(0.87, 0.22), 
        legend.background = element_rect(fill = NA), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(family="LM Roman 10"),
        legend.title = element_text(size = 14),
	plot.title = element_text(size = 15),
        legend.text = element_text(size = 12))

constant.WOTH.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = trend.constant), interpolate = TRUE) +
  geom_sf(data = east.us, alpha = 0, col = 'grey') +
  geom_sf(data = within.sf.list[[which(sp.names == 'WOTH')]], alpha = 0, col = 'black') + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA, 
		       limits = c(min.trend, max.trend)) + 
  theme_bw(base_size = 15) +
  labs(fill = "Trend\n(logit scale)", 
       title = '(I) WOTH Linear') +
  theme(legend.position = c(0.87, 0.22), 
        legend.background = element_rect(fill = NA), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_text(size = 14),
        text = element_text(family="LM Roman 10"),
	plot.title = element_text(size = 15),
        legend.text = element_text(size = 12))
bcr.WOTH.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = trend.bcr), interpolate = TRUE) +
  geom_sf(data = bcrs.ne.grouped, alpha = 0, col = 'grey') +
  geom_sf(data = within.sf.list[[which(sp.names == 'WOTH')]], alpha = 0, col = 'black') + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA, 
		       limits = c(min.trend, max.trend)) + 
  theme_bw(base_size = 15) +
  labs(fill = "Trend\n(logit scale)", 
       title = '(J) WOTH Strata') +
  theme(legend.position = c(0.87, 0.22), 
        legend.background = element_rect(fill = NA), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(family="LM Roman 10"),
        legend.title = element_text(size = 14),
	plot.title = element_text(size = 15),
        legend.text = element_text(size = 12))

comparison.plot <- ggarrange(constant.GRCA.plot, bcr.GRCA.plot, tmax.GRCA.plot, svc.GRCA.plot, 
          constant.EAPH.plot, bcr.EAPH.plot, tmax.EAPH.plot, svc.EAPH.plot,
	  constant.WOTH.plot, bcr.WOTH.plot, tmax.WOTH.plot, svc.WOTH.plot,
	  ncol = 4, nrow = 3, common.legend = TRUE, legend = 'bottom')
# Generate Figure 3
ggsave(plot = comparison.plot, device = 'png', 
       filename = 'figures/case-study-1/Figure-3.png', 
       height = 11, width = 13, units = 'in', bg = 'white')

