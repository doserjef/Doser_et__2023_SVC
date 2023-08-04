# summary-sim.R: script to summarize the simulation results presented in the 
#                "How do spatially-varying coefficient models compare to 
#                simpler alternatives?" section. This script produces Figure 
#                1 in the manuscript.
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(ggthemes)

# Read in simulations -----------------------------------------------------
# The five different candidate models
model.types <- c('linear', 'quadratic', 'stratum', 'interaction', 'svc')
# Linear scenario
load("results/sims/sim-true-linear.rda")
beta.means.linear <- beta.effect.mean.vals
beta.true.linear <- beta.effect.true.vals
waic.linear <- waic.vals
# Quadratic scenario
load("results/sims/sim-true-quadratic.rda")
beta.means.quadratic <- beta.effect.mean.vals
beta.true.quadratic <- beta.effect.true.vals
waic.quadratic <- waic.vals
# Stratum scenario
load("results/sims/sim-true-strata.rda")
beta.means.strata <- beta.effect.mean.vals
beta.true.strata <- beta.effect.true.vals
waic.strata <- waic.vals
# Interaction scenario
load("results/sims/sim-true-interaction.rda")
beta.means.interaction <- beta.effect.mean.vals
beta.true.interaction <- beta.effect.true.vals
waic.interaction <- waic.vals
# Missing interaction scenario
load("results/sims/sim-true-missing.rda")
beta.means.missing <- beta.effect.mean.vals
beta.true.missing <- beta.effect.true.vals
waic.missing <- waic.vals
# Full scenario
load("results/sims/sim-true-full.rda")
beta.means.full <- beta.effect.mean.vals
beta.true.full <- beta.effect.true.vals
waic.full <- waic.vals

# Compare WAIC values across models ---------------------------------------
# These values are those shown in Supplemental Information S1 Table 1
apply(waic.linear, 2, mean)
apply(waic.quadratic, 2, mean)
apply(waic.strata, 2, mean)
apply(waic.interaction, 2, mean)
apply(waic.missing, 2, mean)
apply(waic.full, 2, mean)

# Generate Figure 1 as summary of simulation study ------------------------
# Get coordinates of the sites
J <- ncol(beta.means.linear)
J.x <- sqrt(J)
J.y <- sqrt(J)
# Matrix of spatial locations
s.x <- seq(0, 1, length.out = J.x)
s.y <- seq(0, 1, length.out = J.y)
coords <- as.matrix(expand.grid(s.x, s.y))

# Mean values across simulations and the true values
linear.means <- apply(beta.means.linear, c(2, 3), mean)
linear.true <- beta.true.linear[1, ]
quadratic.means <- apply(beta.means.quadratic, c(2, 3), mean)
quadratic.true <- beta.true.quadratic[1, ]
strata.means <- apply(beta.means.strata, c(2, 3), mean)
strata.true <- beta.true.strata[1, ]
interaction.means <- apply(beta.means.interaction, c(2, 3), mean)
interaction.true <- beta.true.interaction[1, ]
missing.means <- apply(beta.means.missing, c(2, 3), mean)
missing.true <- beta.true.missing[1, ]
full.means <- apply(beta.means.full, c(2, 3), mean)
full.true <- beta.true.full[1, ]

# NOTE: hardcoded
# Number of simulation scenarios
n.scenarios <- 6 # linear, quadratic, strata, interaction, missing, full
# Number of models (plus the truth)
n.models <- 6 # linear, quadratic, strata, interaction, SVC, TRUE

# Put everything in a data frame for plotting
plot.df <- data.frame(easting = rep(coords[, 1], times = n.scenarios * n.models),
		      northing = rep(coords[, 2], times = n.scenarios * n.models),
		      val = c(c(linear.means), c(linear.true), 
				c(quadratic.means), c(quadratic.true), 
				c(strata.means), c(strata.true),
				c(interaction.means), c(interaction.true),
				c(missing.means), c(missing.true),
				c(full.means), c(full.true)),
		      model = factor(rep(rep(c('Linear Estimate', 'Quadratic Estimate', 
					       'Stratum Estimate', 'Interaction Estimate', 
					       'SVC Estimate', 'True'), each = J), 
					 times = n.scenarios), 
				     levels = c('True', 'Linear Estimate', 
						'Quadratic Estimate', 'Stratum Estimate', 
						'Interaction Estimate', 'SVC Estimate'), 
				     ordered = TRUE), 
                      scenario = factor(rep(c('Linear', 'Quadratic', 'Stratum', 
				              'Interaction', 'Missing Interaction', 
					      'Full'), each = J * n.models), 
					levels = c('Linear', 'Quadratic', 'Stratum', 
						   'Interaction', 'Missing Interaction', 
						   'Full'), ordered = TRUE))
min.val <- min(plot.df$val)
max.val <- max(plot.df$val)

ggplot(data = plot.df, aes(x = easting, y = northing, fill = val)) + 
  geom_raster() + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA, limits = c(min.val, max.val)) + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_light(base_size = 14) + 
  guides(fill = 'none') +
  facet_grid(scenario ~ model) +
  labs(x = 'Model', y = 'True Relationship') +
  theme(text = element_text(family="LM Roman 10"), 
        axis.ticks.x = element_blank(), 
        strip.text.y = element_text(color = 'black'),
        strip.text.x = element_text(color = 'black'), 
	axis.text.x = element_blank(),
	axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
# Save to hard drive if desired
ggsave(file = 'figures/Figure-1.png', units = 'in', device = 'png', width = 11, height = 10)
