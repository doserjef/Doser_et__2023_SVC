# summary-sim-4.R: this script summarizes the fourth simulation study from
#                  the manuscript to assess how sample size and the number 
#                  of time periods influences the bias and precision of an SVC
#                  trend from a multiseason SVC occupancy model.
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(viridis)
library(devtools)
library(ggthemes)
library(pals)
library(ggpubr)
library(spOccupancy)

# Functions ---------------------------------------------------------------
logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

# Load results ------------------------------------------------------------
load("results/sim-4-results.rda")
# Simulation scenarios
n.sims <- ncol(psi.mean.samples)
J.vals <- c(100, 200, 400, 800, 1600)
# Setting these to characters here makes the plotting nicer
n.time.vals <- c('5 years', '10 years', '15 years')
# Total number of simulation scenarios
n.scenarios <- length(J.vals) * length(n.time.vals)
scenario.vals <- expand.grid(J = J.vals, n.time = n.time.vals)

# Calculate SVC bias ------------------------------------------------------
w.abs.bias <- array(NA, dim = c(1, n.sims, n.scenarios))
for (i in 1:n.sims) {
  for (j in 1:n.scenarios) {
    indx.1 <- 1:scenario.vals$J[j]
    w.abs.bias[1, i, j] <- mean(abs(abs(w.mean.samples[indx.1, 2, i, j] - w.true[indx.1, 2, i, j])))
  } # j
} # i

w.abs.bias.avg <- apply(w.abs.bias, 3, mean, na.rm = TRUE)
plot.w.abs.bias <- data.frame(J = factor(scenario.vals$J),
			      n.time = factor(scenario.vals$n.time),
                              w.abs.bias = w.abs.bias.avg)

plot.w.abs.bias.2 <- plot.w.abs.bias %>%
  group_by(J, n.time) %>%
  summarize(w.abs.bias = mean(w.abs.bias))

# This is part of Table S6
plot.w.abs.bias.2 %>%
  arrange(J, n.time) %>%
  print(n = nrow(.))

# Calculate SVC coverage --------------------------------------------------
w.coverage.rates <- array(NA, dim = c(1, n.sims, n.scenarios))
for (i in 1:n.sims) {
  for (j in 1:n.scenarios) {
    indx.1 <- 1:scenario.vals$J[j]
    w.coverage.rates[1, i, j] <- mean(ifelse((w.true[indx.1, 2, i, j] > w.low.samples[indx.1, 2, i, j]) &
					      (w.true[indx.1, 2, i, j] < w.high.samples[indx.1, 2, i, j]),
				      1, 0))
  } # j
} # i

w.coverage.avg <- apply(w.coverage.rates, 3, mean, na.rm = TRUE)
plot.w.coverage <- data.frame(J = factor(scenario.vals$J),
			      n.time = factor(scenario.vals$n.time),
                              w.coverage = w.coverage.avg)

plot.w.coverage.2 <- plot.w.coverage %>%
  group_by(J, n.time) %>%
  summarize(w.coverage = mean(w.coverage))

# This is part of Table S6
plot.w.coverage.2 %>%
  arrange(J, n.time) %>%
  print(n = nrow(.))

# Calculate SVC CI Width -------------------------------------------------
w.ci.widths <- array(NA, dim = c(1, n.sims, n.scenarios))
for (i in 1:n.sims) {
  for (j in 1:n.scenarios) {
    indx.1 <- 1:scenario.vals$J[j]
    w.ci.widths[1, i, j] <- mean(w.high.samples[indx.1, 2, i, j]  - w.low.samples[indx.1, 2, i, j])
  } # j
} # i

w.ci.width.avg <- apply(w.ci.widths, 3, mean, na.rm = TRUE)
plot.w.ci.width <- data.frame(J = factor(scenario.vals$J),
			      n.time = factor(scenario.vals$n.time),
                              w.ci.width = w.ci.width.avg)

plot.w.ci.width.2 <- plot.w.ci.width %>%
  group_by(J, n.time) %>%
  summarize(w.ci.width = mean(w.ci.width))

# This is part of Table S6
plot.w.ci.width.2 %>%
  arrange(J, n.time) %>%
  print(n = nrow(.))

# Correlation of SVC trend with the truth ---------------------------------
# NOTE: this is hardcoded to just look at the spatially-varying trend, not 
#       also the intercept.
w.cor <- array(NA, dim = c(1, n.sims, n.scenarios))
for (i in 1:n.sims) {
  for (j in 1:n.scenarios) {
    indx.1 <- 1:scenario.vals$J[j]
    p.svc.curr <- scenario.vals$p.svc[j]
    w.cor[1, i, j] <- cor(w.mean.samples[indx.1, 2, i, j], w.true[indx.1, 2, i, j])
  } # j
} # i

w.cor.avg <- apply(w.cor, 3, mean, na.rm = TRUE)
plot.w.cor <- data.frame(J = scenario.vals$J,
			 n.time = factor(scenario.vals$n.time),
                         w.cor = w.cor.avg)

occ.cor.data <- plot.w.cor %>%
  group_by(J, n.time) %>%
  summarize(avg.cor = mean(w.cor)) %>%
  mutate(type = 'Occupancy')

# Create Figure S2
cor.w.plot <- ggplot(occ.cor.data, aes(x = J, y = avg.cor, fill = n.time)) +
  geom_line(aes(group = n.time, col = n.time), size = 1) +
  geom_point(size = 3, pch = 21) +
  scale_color_colorblind() +
  scale_fill_colorblind() +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(size = 18)) +
  labs(x = 'Number of sites', y = 'Correlation coefficient', fill = 'Number of time periods',
       col = 'Number of time periods') + 
  theme(legend.position = 'bottom')
cor.w.plot
ggsave(file = 'figures/Fig-S2.pdf', units = 'in', device = 'pdf', width = 7, height = 5)

