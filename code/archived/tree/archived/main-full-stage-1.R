# main-full-stage-1.R: this script fits a spatially-explicit GLM for sugar
#                      maple across the US.
rm(list = ls())
library(tidyverse)
library(sf)
library(spOccupancy)

# Load data ---------------------------------------------------------------
load("data/tree/stage-1-data.rda")

# Prep the model ----------------------------------------------------------
# Priors
# dist.fia <- dist(data.list$coords)
# mean.dist <- mean(dist.fia) # 1783.395
# min.dist <- min(dist.fia) # 0.007832157
# max.dist <- max(dist.fia) # 4642.168
mean.dist <- 1783.395
min.dist <- 0.007832157
max.dist <- 4642.168
# Might need further restrictions on phi later on. 
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
		   phi.unif = list(a = 3 / max.dist, b = 3 / min.dist), 
                   sigma.sq.ig = list(2, 5))
inits.list <- list(beta = 0, phi = 3 / mean.dist, sigma.sq = 1)
tuning.list <- list(phi = 0.05)
n.neighbors <- 5
# Should probably explore this more thoroughly. 
n.factors <- 1
cov.model <- "exponential"
n.chains <- 1

# Biggest
n.batch <- 2000
batch.length <- 25
n.burn <- 30000
n.thin <- 20

# Medium
n.batch <- 200
batch.length <- 25
n.burn <- 3000
n.thin <- 5

n.batch <- 1
batch.length <- 25
n.burn <- 0
n.thin <- 1

out <- svcPGBinom(formula = ~ scale(tmin) + I(scale(tmin)^2) + 
                              scale(ppt) + I(scale(ppt)^2) + 
			      scale(elev) + I(scale(elev)^2), 
		  data = data.list, priors = prior.list, 
		  inits = inits.list, tuning = tuning.list, svc.cols = c(1),
	          n.neighbors = n.neighbors, cov.model = cov.model, NNGP = TRUE,
	          n.batch = n.batch, batch.length = batch.length, 
	          n.burn = n.burn, accept.rate = 0.43, n.thin = n.thin, 
	          n.chains = n.chains, n.report = 1, n.omp.threads = 10)

svc.samples <- getSVCSamples(out)
psi.samples <- out$psi.samples

save(svc.samples, psi.samples, 
     file = 'results/stage-1-sugar-maple-full-samples.rda') 

# w.means <- apply(out$w.samples[, 1, ], 2, mean)
# psi.ci.width <- psi.high - psi.low
# y.rep.means <- apply(out$y.rep.samples, 2, mean)
# 
# plot.df <- data.frame(w.mean = w.means, 
# 		      psi.mean = psi.means, 
# 		      psi.ci.width = psi.ci.width,
# 		      x = coords[curr.indx, 1], 
# 		      y = coords[curr.indx, 2])
# 
# plot.sf <- st_as_sf(plot.df,
# 		    coords = c("x", "y"),
# 		    crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs")
# 
# ggplot(plot.sf) +
#   geom_sf(aes(col = w.mean)) +
#   scale_color_viridis() +
#   theme_bw(base_size = 18)
# 
# psi.mean.plot <- ggplot(plot.sf) +
#   geom_sf(aes(col = psi.mean)) +
#   scale_color_viridis() +
#   theme_bw(base_size = 18)
# 
# psi.ci.width.plot <- ggplot(plot.sf) +
#   geom_sf(aes(col = psi.ci.width)) +
#   scale_color_viridis() +
#   theme_bw(base_size = 18)
# 
# ggarrange(psi.mean.plot, psi.ci.width.plot, labels = c('Mean', 'CI Width'))

