rm(list = ls())
library(spOccupancy)
library(tidyverse)
library(ggpubr)
library(MASS)
library(spBayes)
library(viridis)

# Simulate the occupancy process ------------------------------------------
set.seed(500)
# Intercept + linear + quadratic + stratum + interaction + missing interaction
J.x <- 50
J.y <- 50
J <- J.x * J.y
# Subroutines -----------------------------------------------------------
logit <- function(theta, a = 0, b = 1){log((theta-a)/(b-theta))}
logit.inv <- function(z, a = 0, b = 1){b-(b-a)/(1+exp(z))}

# Matrix of spatial locations
s.x <- seq(0, 1, length.out = J.x)
s.y <- seq(0, 1, length.out = J.y)
coords <- as.matrix(expand.grid(s.x, s.y))

# Get strata for each cell
strata <- ifelse(coords[, 1] <= .33 & coords[, 2] <= .33, 1,
		 ifelse(coords[, 1] <= .33 & coords[, 2] <= .67, 2,
			ifelse(coords[, 1] <= .33 & coords[, 2] <= 1, 3,
			       ifelse(coords[, 1] <= .67 & coords[, 2] <= .33, 4,
				      ifelse(coords[, 1] <= .67 & coords[, 2] <= .67, 5,
					     ifelse(coords[, 1] <= .67 & coords[, 2] <= 1, 6,
						    ifelse(coords[, 1] <= 1 & coords[, 2] <= .33, 7,
							   ifelse(coords[, 1] <= 1 & coords[, 2] <= .67, 8,
								  9))))))))

# Main covariate
beta.0 <- 0
beta.linear <- 0.5
beta.quadratic <- -0.5
beta.strata <- rnorm(n_distinct(strata))
beta.interaction <- 0.5
beta.miss.int <- 0.5
# Main covariate
x.1 <- seq(from = -5, to = 5, length.out = J)
x.interaction <- c(t(matrix(x.1, J.x, J.y)))
sigma.sq <- 3
phi <- 3 / .8
Sigma <- mkSpCov(coords, as.matrix(sigma.sq), as.matrix(0), as.matrix(phi), 'exponential')
x.miss.int <- mvrnorm(1, rep(0, J), Sigma)

# Form detection covariate (if any) -------------------------------------
n.rep <- rep(4, J)
n.rep.max <- max(n.rep)
alpha <- c(0.5, -0.3)
n.alpha <- length(alpha)
X.p <- array(NA, dim = c(J, n.rep.max, n.alpha))
# Get index of surveyed replicates for each site. 
rep.indx <- list()
for (j in 1:J) {
  rep.indx[[j]] <- sample(1:n.rep.max, n.rep[j], replace = FALSE)
}
X.p[, , 1] <- 1
if (n.alpha > 1) {
  for (i in 2:n.alpha) {
    for (j in 1:J) {
      X.p[j, rep.indx[[j]], i] <- rnorm(n.rep[j])
    } # j
  } # i
}

# Simulate spatial random effect ----------------------------------------
svc.cols <- c(1)
p.svc <- length(svc.cols)
cov.model <- 'exponential'
w.mat <- matrix(NA, J, p.svc)
phi <- c(3 / 0.5)
sigma.sq <- 0.8
theta <- as.matrix(phi)
for (i in 1:p.svc) {
  Sigma <- mkSpCov(coords, as.matrix(sigma.sq[i]), as.matrix(0), theta[i, ], cov.model)
  # Random spatial process
  w.mat[, i] <- mvrnorm(1, rep(0, J), Sigma)
}
X.w <- matrix(1, J, 1) 
w <- c(t(w.mat))

# Latent Occupancy Process ----------------------------------------------
psi <- logit.inv(beta.0 + w + beta.linear * x.1 + beta.quadratic * x.1^2 + 
		 beta.strata[strata] * x.1 + beta.interaction * x.1 * x.interaction + 
		 beta.miss.int * x.1 * x.miss.int)
z <- rbinom(J, 1, psi)

# Data Formation --------------------------------------------------------
p <- matrix(NA, nrow = J, ncol = n.rep.max)
y <- matrix(NA, nrow = J, ncol = n.rep.max)
for (j in 1:J) {
  p[j, rep.indx[[j]]] <- logit.inv(X.p[j, rep.indx[[j]], ] %*% as.matrix(alpha))
  y[j, rep.indx[[j]]] <- rbinom(n.rep[j], 1, p[j, rep.indx[[j]]] * z[j])
} # j

# Plot of the full covariate effect --------------------------------------
full.cov.effect <- beta.linear + beta.quadratic * x.1 + beta.strata[strata] + 
                   beta.interaction * x.interaction + beta.miss.int * x.miss.int
plot.df <- data.frame(easting = coords[, 1],
		      northing = coords[, 2],
                      beta.linear = beta.linear,
		      x.1 = x.1, 
		      x.interaction = x.interaction, 
		      x.miss.int = x.miss.int,
		      beta.quadratic = beta.quadratic * x.1,
		      beta.strata = beta.strata[strata],
		      beta.interaction = beta.interaction * x.interaction,
		      beta.miss.int = beta.miss.int * x.miss.int,
		      beta.full = full.cov.effect)
linear.map <- ggplot(data = plot.df, aes(x = easting, y = northing, fill = beta.linear)) + 
  geom_raster() + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA, limits = c(-5, 5)) + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw(base_size = 10) + 
  labs(x = 'Easting', y = 'Northing', fill = '', title = 'Linear') +
  guides(fill = 'none') +
  theme(text = element_text(family="LM Roman 10"), 
        axis.ticks.x = element_blank(), 
	axis.text.x = element_blank(),
	#plot.background = element_rect(fill = 'cornsilk', color = 'black'),
	axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave(plot = linear.map, file = 'figures/main-fig-linear.png', device = 'png', 
       units = 'in', width = 2, height = 2)
quadratic.map <- ggplot(data = plot.df, aes(x = easting, y = northing, fill = beta.quadratic)) + 
  geom_raster() + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA) + 
  theme_bw(base_size = 10) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = 'Easting', y = 'Northing', title = 'Quadratic', fill = '') +
  guides(fill = 'none') +
  theme(text = element_text(family="LM Roman 10"), 
        axis.ticks.x = element_blank(), 
	axis.text.x = element_blank(),
	axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave(plot = quadratic.map, file = 'figures/main-fig-quadratic.png', device = 'png', 
       units = 'in', width = 2, height = 2)

strata.map <- ggplot(data = plot.df, aes(x = easting, y = northing, fill = beta.strata)) + 
  geom_raster() + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA) + 
  theme_bw(base_size = 10) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill = 'none') +
  labs(x = 'Easting', y = 'Northing', title = 'Stratum', fill = '') +
  theme(text = element_text(family="LM Roman 10"), 
        axis.ticks.x = element_blank(), 
	axis.text.x = element_blank(),
	axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave(plot = strata.map, file = 'figures/main-fig-strata.png', device = 'png', 
       units = 'in', width = 2, height = 2)
interaction.map <- ggplot(data = plot.df, aes(x = easting, y = northing, fill = beta.interaction)) + 
  geom_raster() + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA) + 
  theme_bw(base_size = 10) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill = 'none') +
  labs(x = 'Easting', y = 'Northing', fill = '', title = 'Interaction') +
  theme(text = element_text(family="LM Roman 10"), 
	axis.text.x = element_blank(),
	axis.text.y = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank())
ggsave(plot = interaction.map, file = 'figures/main-fig-interaction.png', device = 'png', 
       units = 'in', width = 2, height = 2)
miss.int.map <- ggplot(data = plot.df, aes(x = easting, y = northing, fill = beta.miss.int)) + 
  geom_raster() + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA) + 
  theme_bw(base_size = 10) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill = 'none') +
  labs(x = 'Easting', y = 'Northing', fill = '', title = 'Unknown interactions') +
  theme(text = element_text(family="LM Roman 10"), 
	axis.text.x = element_blank(),
	axis.text.y = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank())
ggsave(plot = miss.int.map, file = 'figures/main-fig-miss-int.png', device = 'png', 
       units = 'in', width = 2, height = 2)
full.effect.map <- ggplot(data = plot.df, aes(x = easting, y = northing, fill = beta.full)) + 
  geom_raster() + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA) + 
  theme_bw(base_size = 10) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = 'Easting', y = 'Northing', title = 'Full Effect', fill = '') +
  guides(fill = 'none') +
  theme(text = element_text(family="LM Roman 10"), 
        axis.ticks.x = element_blank(), 
	axis.text.x = element_blank(),
	axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave(plot = full.effect.map, file = 'figures/main-fig-full-effect.png', device = 'png', 
       units = 'in', width = 2, height = 2)

cov.1.map <- ggplot(data = plot.df, aes(x = easting, y = northing, fill = x.1)) + 
  geom_raster() + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw(base_size = 10) + 
  labs(x = 'Easting', y = 'Northing', fill = '', title = 'Covariate') +
  guides(fill = 'none') + 
  theme(text = element_text(family="LM Roman 10"), 
        axis.ticks.x = element_blank(), 
	# plot.background = element_rect(fill = 'gray', color = 'black'),
	axis.text.x = element_blank(),
	axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

cov.int.map <- ggplot(data = plot.df, aes(x = easting, y = northing, fill = x.interaction)) + 
  geom_raster() + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw(base_size = 10) + 
  labs(x = 'Easting', y = 'Northing', fill = '', title = 'Interacting Effect') +
  guides(fill = 'none') + 
  theme(text = element_text(family="LM Roman 10"), 
        axis.ticks.x = element_blank(), 
	# plot.background = element_rect(fill = 'gray', color = 'black'),
	axis.text.x = element_blank(),
	axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

cov.miss.int.map <- ggplot(data = plot.df, aes(x = easting, y = northing, fill = x.miss.int)) + 
  geom_raster() + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw(base_size = 10) + 
  labs(x = 'Easting', y = 'Northing', fill = '', title = 'Missing interactions') +
  guides(fill = 'none') + 
  theme(text = element_text(family="LM Roman 10"), 
        axis.ticks.x = element_blank(), 
	# plot.background = element_rect(fill = 'gray', color = 'black'),
	axis.text.x = element_blank(),
	axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

# Fit an SVC model --------------------------------------------------------
data.list <- list(y = y, 
		  occ.covs = data.frame(x.1 = x.1),
		  det.covs = list(det.cov.1 = X.p[, , 2]), 
                  coords = coords)
priors <- list(phi.unif = list(a = 3 / .9, b = 3 / .1))

n.batch <- 800
batch.length <- 25
n.chains <- 1
n.burn <- 10000
n.thin <- 10
out <- svcPGOcc(occ.formula = ~ x.1, 
		det.formula = ~ det.cov.1, 
		data = data.list,
		priors = priors, 
		svc.cols = c(1, 2),
		tuning = list(phi = 0.5), 
		n.batch = n.batch,
		batch.length = batch.length, 
		n.chains = n.chains, 
		n.neighbors = 10,
		n.burn = n.burn, 
		n.thin = n.thin, 
                n.report = 10)
svc.samples <- getSVCSamples(out)
svc.means <- apply(svc.samples$x.1, 2, mean)

plot.df$svc.est <- svc.means
svc.est.map <- ggplot(data = plot.df, aes(x = easting, y = northing, fill = svc.est)) + 
  geom_raster() + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA) + 
  theme_bw(base_size = 10) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = 'Easting', y = 'Northing', title = 'SVC Estimate', fill = '') +
  guides(fill = 'none') +
  theme(text = element_text(family="LM Roman 10"), 
        axis.ticks.x = element_blank(), 
	plot.background = element_rect(fill = 'transparent', color = NA),
	axis.text.x = element_blank(),
	axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave(plot = svc.est.map, file = 'figures/main-fig-svc-est.png', device = 'png', 
       units = 'in', width = 2, height = 2, bg = NULL)

# Put it all together -----------------------------------------------------
# ggarrange(linear.curve, quadratic.curve, strata.curve, interaction.curve, 
# 	  svc.curve, linear.map, quadratic.map, strata.map, interaction.map, 
# 	  svc.map, ncol = 5, nrow = 2)
# ggsave(file = 'figures/main-figure.png', device = 'png', units = 'in', width = 10, height = 5)
# ggarrange(linear.curve, linear.map, 
# 	  quadratic.curve, quadratic.map, 
# 	  strata.curve, strata.map, 
# 	  interaction.curve, interaction.map, 
# 	  svc.curve, svc.map, ncol = 2, nrow = 5)
# ggsave(file = 'figures/main-figure.png', device = 'png', units = 'in', height = 12, width = 5)
