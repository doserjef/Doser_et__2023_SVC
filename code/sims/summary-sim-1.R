# summary-sim-1.R: this script summarizes results from the first simulation study
#                  in which we essentially perform a proof of concept for the SVC
#                  occupancy model across varying types of nonstationarity in the 
#                  covariate effect. Additionally, this script contains summary
#                  results included in the appendix on the same simulation study
#                  done without imperfect detection for binary and binomial data.
rm(list = ls())
library(tidyverse)
library(viridis)
library(devtools)
library(ggthemes)
library(spOccupancy)

# SVC occupancy model results ---------------------------------------------
# Load results ------------------------------------------------------------
load("results/sim-1-results.rda")
# Number of sites
J <- dim(psi.mean.samples)[1]
# Number of simulations for each scenario/model combo
n.sims <- dim(psi.mean.samples)[2]
# Number of scenarios
n.scenarios <- dim(psi.mean.samples)[3]
# Number of models
n.models <- dim(psi.mean.samples)[4]

# Simulated spatial variance parameters
sigma.sq.vals <- c(0.1, 0.5, 1, 2)
phi.vals <- c(3 / 0.1, 3 / 0.3, 3 / 0.8)
scenario.vals <- expand.grid(sigma.sq = sigma.sq.vals, phi = phi.vals)

# Assessment of model deviance --------------------------------------------
# Average deviance score across all simulations
deviance.avg <- apply(deviance.vals, c(2, 3), mean)
deviance.low <- apply(deviance.vals, c(2, 3), quantile, 0.025)
deviance.high <- apply(deviance.vals, c(2, 3), quantile, 0.975)
plot.deviance <- data.frame(phi = factor(rep(scenario.vals$phi, times = n.models)),
			    sigma.sq = factor(rep(scenario.vals$sigma.sq, times = n.models)),
			    deviance = c(deviance.avg) - deviance.avg[, 1],
			    type = factor(rep(c('SVC', 'SVI', 'OCC'), each = n.scenarios), 
			                  levels = c('OCC', 'SVI', 'SVC'), ordered = TRUE))
levels(plot.deviance$phi) <- c('High Range', 'Medium Range', 'Low Range')
levels(plot.deviance$sigma.sq) <- c('Low Variance', 'Medium Low Variance', 'Medium High Variance', 'High Variance')

# This is part of Table 1
plot.deviance %>%
  arrange(desc(phi), sigma.sq, type)

# Assessment of WAIC ------------------------------------------------------
# Average WAIC score across all simulations
waic.avg <- apply(waic.vals, c(2, 3), mean)
plot.waic <- data.frame(phi = factor(rep(scenario.vals$phi, times = n.models)),
			    sigma.sq = factor(rep(scenario.vals$sigma.sq, times = n.models)),
			    waic = c(waic.avg) - waic.avg[, 1],
			    type = factor(rep(c('SVC', 'SVI', 'OCC'), each = n.scenarios), 
			                  levels = c('OCC', 'SVI', 'SVC'), ordered = TRUE))
levels(plot.waic$phi) <- c('High Range', 'Medium Range', 'Low Range')
levels(plot.waic$sigma.sq) <- c('Low Variance', 'Medium Low Variance', 'Medium High Variance', 'High Variance')

# This is the other part of Table 1
plot.waic %>%
  arrange(desc(phi), sigma.sq, type)

# Absolute bias -----------------------------------------------------------
psi.abs.bias <- array(NA, dim = c(n.sims, n.scenarios))
for (i in 1:n.sims) {
  for (j in 1:n.scenarios) {
    # psi.abs.bias[i, j] <- mean(abs(psi.mean.samples[1:scenario.vals$J[j], i, j] - psi.true[1:scenario.vals$J[j], i, j]) / psi.true[1:scenario.vals$J[j], i, j])
    psi.abs.bias[i, j] <- mean(abs(psi.mean.samples[, i, j, 1] - psi.true[, i, j]))
  } # j
} # i

psi.abs.bias.avg <- apply(psi.abs.bias, 2, mean)
plot.psi.abs.bias <- data.frame(phi = scenario.vals$phi,
				sigma.sq = factor(scenario.vals$sigma.sq),
				psi.abs.bias = psi.abs.bias.avg)
plot.psi.abs.bias

# Calculate average CI width ----------------------------------------------
psi.ci.width <- array(NA, dim = c(n.sims, n.scenarios))
for (i in 1:n.sims) {
  for (j in 1:n.scenarios) {
    psi.ci.width[i, j] <- mean(psi.high.samples[, i, j, 1] - psi.low.samples[, i, j, 1])
  } # j
} # i

psi.ci.width.avg <- apply(psi.ci.width, 2, mean)
plot.psi.ci.width <- data.frame(phi = scenario.vals$phi,
				sigma.sq = factor(scenario.vals$sigma.sq),
				psi.ci.width = psi.ci.width.avg)

# SVC GLM results ---------------------------------------------------------
# Load results ------------------------------------------------------------
load("results/sim-1-glm-results.rda")

# Absolute bias -----------------------------------------------------------
psi.abs.glm.bias <- array(NA, dim = c(n.sims, n.scenarios))
for (i in 1:n.sims) {
  for (j in 1:n.scenarios) {
    # psi.abs.bias[i, j] <- mean(abs(psi.mean.samples[1:scenario.vals$J[j], i, j] - psi.true[1:scenario.vals$J[j], i, j]) / psi.true[1:scenario.vals$J[j], i, j])
    psi.abs.glm.bias[i, j] <- mean(abs(psi.mean.samples[, i, j, 1] - psi.true[, i, j]))
  } # j
} # i

psi.abs.glm.bias.avg <- apply(psi.abs.glm.bias, 2, mean)
plot.psi.abs.glm.bias <- data.frame(phi = scenario.vals$phi,
				    sigma.sq = factor(scenario.vals$sigma.sq),
				    psi.abs.bias = psi.abs.glm.bias.avg)
plot.psi.abs.glm.bias

# Calculate average CI width ----------------------------------------------
psi.ci.width.glm <- array(NA, dim = c(n.sims, n.scenarios))
for (i in 1:n.sims) {
  for (j in 1:n.scenarios) {
    psi.ci.width.glm[i, j] <- mean(psi.high.samples[, i, j, 1] - psi.low.samples[, i, j, 1])
  } # j
} # i

psi.ci.width.avg.glm <- apply(psi.ci.width.glm, 2, mean)
plot.psi.ci.width.glm <- data.frame(phi = scenario.vals$phi,
				    sigma.sq = factor(scenario.vals$sigma.sq),
				    psi.ci.width = psi.ci.width.avg.glm)

# Compare occupancy model with GLM ----------------------------------------
# Bias
(plot.psi.abs.bias$bias.diff <- plot.psi.abs.bias$psi.abs.bias - 
                                plot.psi.abs.glm.bias$psi.abs.bias)
# Percentage change in bias 
mean(plot.psi.abs.bias$bias.diff / plot.psi.abs.bias$psi.abs.bias * 100)
# CI Width
(plot.psi.ci.width$ci.diff <- plot.psi.ci.width$psi.ci.width - 
                              plot.psi.ci.width.glm$psi.ci.width)
mean(plot.psi.ci.width$ci.diff / plot.psi.ci.width$psi.ci.width * 100)

# Figure of spatial range and spatial variance ----------------------------
# This code creates Figure S1, which gives an example of how the spatial
# variance and spatial decay influence the resulting spatially-varying
# coefficient. 
# Spatial locations
J.x <- 50
J.y <- 50
J <- J.x * J.y
# Replicates
# n.rep <- sample(5, J, replace = TRUE)
n.rep <- rep(5, J)
# Occurrence coefficients -------------
beta <- c(0, 0)
# Detection coefficients --------------
alpha <- c(-0.2, 0.4)
# Spatial parameters ------------------
sp <- TRUE
svc.cols <- c(1, 2)
cov.model <- 'exponential'
# Values for spatially-varying intercept
sigma.sq.int <- 1
phi.int <- 3 / 0.4
plot.data <- list()
for (i in 1:nrow(scenario.vals)) {
  set.seed(150)
  print(i)
  dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
                sp = sp, svc.cols = svc.cols, cov.model = cov.model,
                sigma.sq = c(sigma.sq.int, scenario.vals$sigma.sq[i]),
                phi = c(phi.int, scenario.vals$phi[i]))
  
  plot.data[[i]] <- data.frame(x = dat$coords[, 1], 
  		               y = dat$coords[, 2], 
  		               w = dat$w[, 2])
}

plot.data.df <- do.call(rbind, plot.data)
plot.data.df <- plot.data.df %>%
  mutate(phi = factor(rep(scenario.vals$phi, each = J), 
		      levels = c(30, 10, 3.75), ordered = TRUE), 
	 sigma.sq = factor(rep(scenario.vals$sigma.sq, each = J)))
lower.val <- min(sapply(plot.data, function(a) min(a$w)))
upper.val <- max(sapply(plot.data, function(a) max(a$w)))
sigma.sq.labs <- c('$\\sigma^2 = 0.1$', '$\\sigma^2 = 0.5$', 
		   '$\\sigma^2 = 1.0$', '$\\sigma^2 = 2.0$')
names(sigma.sq.labs) <- c('0.1', '0.5', '1', '2')
phi.labs <- c('$\\phi = 3.75$', '$\\phi = 10$', '$\\phi = 30$')
names(phi.labs) <- c('3.75', '10', '30')


tikzDevice::tikz(file = 'figures/theta-examples.tex', width = 8, height = 5)
my.plot <- ggplot(data = plot.data.df, aes(x = x, y = y, fill = w)) + 
    scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
      	               na.value = NA, limits = c(lower.val, upper.val)) + 
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) + 
    geom_raster() + 
    facet_grid(phi ~ sigma.sq, 
	       labeller = labeller(phi = phi.labs, sigma.sq = sigma.sq.labs)) + 
    theme_light(base_size = 14) + 
    labs(x = 'Easting', y = 'Northing', fill = 'Effect') +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          strip.text.x = element_text(color = 'black'), 
          strip.text.y = element_text(color = 'black')) 
my.plot
dev.off()  
