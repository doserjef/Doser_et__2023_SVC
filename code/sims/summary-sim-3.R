# summary-sim-3.R: this script summarizes the third simulation study from
#                  the manuscript in which we assessed the reliability of 
#                  estimates from an SVC multi-season occupancy model under varying 
#                  numbers of spatial locations, varying amounts of 
#                  spatially-varying coefficients, and varying numbers of 
#                  primary time periods
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
load("results/sim-3-results.rda")
# Simulation scenarios
n.sims <- ncol(psi.mean.samples)
J.vals <- c(100, 200, 400, 800, 1600)
# Setting these to characters here makes the plotting nicer
n.time.vals <- c('5 years', '10 years', '15 years')
p.svc.vals <- c(1, 2, 3)
# Total number of simulation scenarios
n.scenarios <- length(J.vals) * length(n.time.vals) * length(p.svc.vals)
scenario.vals <- expand.grid(J = J.vals, n.time = n.time.vals, p.svc = p.svc.vals)

# Calculate occurrence probability bias -----------------------------------
psi.abs.bias <- array(NA, dim = c(n.sims, n.scenarios))
for (i in 1:n.sims) {
  for (j in 1:n.scenarios) {
    psi.abs.bias[i, j] <- mean(abs(psi.mean.samples[1:scenario.vals$J[j], 1:scenario.vals$n.time[j], i, j] - psi.true[1:scenario.vals$J[j], 1:scenario.vals$n.time[j], i, j]))
  } # j
} # i

psi.abs.bias.avg <- apply(psi.abs.bias, 2, mean)
plot.psi.abs.bias <- data.frame(J = scenario.vals$J,
				n.time = factor(scenario.vals$n.time),
				p.svc = factor(scenario.vals$p.svc),
				psi.abs.bias = psi.abs.bias.avg)
plot.psi.abs.bias.2 <- plot.psi.abs.bias %>%
  group_by(J, n.time, p.svc) %>%
  summarize(psi.abs.bias = mean(psi.abs.bias))

# This is part of Table S4
plot.psi.abs.bias.2 %>%
  arrange(p.svc, J, n.time) %>%
  print(n = nrow(.))

psi.bias.plot <- ggplot(plot.psi.abs.bias.2, aes(x = J, y = psi.abs.bias, fill = p.svc)) + 
  geom_line(aes(group = p.svc, col = p.svc), size = 1) +
  geom_point(size = 3, pch = 21) + 
  facet_wrap(vars(n.time)) + 
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme_bw(base_size = 18) + 
  theme(plot.title = element_text(size = 18)) +
  labs(x = 'Number of sites', y = 'Absolute bias', fill = 'Number of SVCs', 
       col = 'Number of SVCs', title = '(A) Occupancy probability bias')
psi.bias.plot

# Calculate average coverage rates ----------------------------------------
psi.coverage.rates <- array(NA, dim = c(n.sims, n.scenarios))
for (i in 1:n.sims) {
  for (j in 1:n.scenarios) {
    indx <- 1:scenario.vals$J[j]
    indx.2 <- 1:scenario.vals$n.time[j]
    psi.coverage.rates[i, j] <- mean(ifelse((psi.true[indx, indx.2, i, j] > psi.low.samples[indx, indx.2, i, j]) & 
						  (psi.true[indx, indx.2, i, j] < psi.high.samples[indx, indx.2, i, j]), 1, 0))
  } # j
} # i

psi.coverage.avg <- apply(psi.coverage.rates, 2, mean)
plot.psi.coverage <- data.frame(J = factor(scenario.vals$J),
				n.time = factor(scenario.vals$n.time),
				p.svc = factor(scenario.vals$p.svc),
				psi.coverage = psi.coverage.avg)
plot.psi.coverage.2 <- plot.psi.coverage %>%
  group_by(J, n.time, p.svc) %>%
  summarize(psi.coverage = mean(psi.coverage))

# This is part of Table S4
plot.psi.coverage.2 %>%
  arrange(p.svc, J, n.time) %>%
  print(n= nrow(.))

# Calculate average CI width ----------------------------------------------
psi.ci.width <- array(NA, dim = c(n.sims, n.scenarios))
for (i in 1:n.sims) {
  for (j in 1:n.scenarios) {
    indx <- 1:scenario.vals$J[j]
    indx.2 <- 1:scenario.vals$n.time[j]
    psi.ci.width[i, j] <- mean(psi.high.samples[indx, indx.2, i, j] - psi.low.samples[indx, indx.2, i, j])
  } # j
} # i

psi.ci.width.avg <- apply(psi.ci.width, 2, mean)
plot.psi.ci.width <- data.frame(J = scenario.vals$J,
				n.time = factor(scenario.vals$n.time),
				p.svc = factor(scenario.vals$p.svc),
				psi.ci.width = psi.ci.width.avg)
plot.psi.ci.width.2 <- plot.psi.ci.width %>%
  group_by(J, n.time, p.svc) %>%
  summarize(psi.ci.width = mean(psi.ci.width))

psi.ci.width.plot <- ggplot(plot.psi.ci.width.2, aes(x = J, y = psi.ci.width, fill = p.svc)) + 
  geom_line(aes(group = p.svc, col = p.svc), size = 1) +
  geom_point(size = 3, pch = 21) + 
  facet_wrap(vars(n.time)) + 
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme_bw(base_size = 18) + 
  theme(plot.title = element_text(size = 18)) +
  labs(x = 'Number of sites', y = '95% CI Width', fill = 'Number of SVCs', 
       col = 'Number of SVCs', title = '(B) Occupancy probability uncertainty')
psi.ci.width.plot

# This is part of Table S4
plot.psi.ci.width.2 %>%
  arrange(p.svc, J, n.time) %>%
  print(n= nrow(.))

# Calculate SVC bias ------------------------------------------------------
p.svc.max <- max(p.svc.vals)
w.abs.bias <- array(NA, dim = c(p.svc.max, n.sims, n.scenarios))
for (i in 1:n.sims) {
  for (j in 1:n.scenarios) {
    indx.1 <- 1:scenario.vals$J[j]
    p.svc.curr <- scenario.vals$p.svc[j]
    for (k in 1:p.svc.curr) {
    # w.abs.bias[k, i, j] <- mean(abs(abs(w.mean.samples[indx.1, k, i, j] - w.true[indx.1, k, i, j]) / w.true[indx.1, k, i, j]))
    w.abs.bias[k, i, j] <- mean(abs(abs(w.mean.samples[indx.1, k, i, j] - w.true[indx.1, k, i, j])))
    } #k
  } # j
} # i

w.abs.bias.avg <- apply(w.abs.bias, 3, mean, na.rm = TRUE)
plot.w.abs.bias <- data.frame(J = factor(scenario.vals$J),
			      n.time = factor(scenario.vals$n.time),
                              p.svc = factor(scenario.vals$p.svc),
                              w.abs.bias = w.abs.bias.avg)

plot.w.abs.bias.2 <- plot.w.abs.bias %>%
  group_by(J, p.svc, n.time) %>%
  summarize(w.abs.bias = mean(w.abs.bias))

# This is part of Table S5
plot.w.abs.bias.2 %>%
  arrange(p.svc, J, n.time) %>%
  print(n = nrow(.))

# Calculate SVC coverage --------------------------------------------------
p.svc.max <- max(p.svc.vals)
w.coverage.rates <- array(NA, dim = c(p.svc.max, n.sims, n.scenarios))
for (i in 1:n.sims) {
  for (j in 1:n.scenarios) {
    indx.1 <- 1:scenario.vals$J[j]
    p.svc.curr <- scenario.vals$p.svc[j]
    for (k in 1:p.svc.curr) {
      w.coverage.rates[k, i, j] <- mean(ifelse((w.true[indx.1, k, i, j] > w.low.samples[indx.1, k, i, j]) &
					      (w.true[indx.1, k, i, j] < w.high.samples[indx.1, k, i, j]),
				      1, 0))
    } #k
  } # j
} # i

w.coverage.avg <- apply(w.coverage.rates, 3, mean, na.rm = TRUE)
plot.w.coverage <- data.frame(J = factor(scenario.vals$J),
			      n.time = factor(scenario.vals$n.time),
                              p.svc = factor(scenario.vals$p.svc),
                              w.coverage = w.coverage.avg)

plot.w.coverage.2 <- plot.w.coverage %>%
  group_by(J, p.svc, n.time) %>%
  summarize(w.coverage = mean(w.coverage))

# This is part of Table S5
plot.w.coverage.2 %>%
  arrange(p.svc, J, n.time) %>%
  print(n = nrow(.))

# Calculate SVC CI Width -------------------------------------------------
p.svc.max <- max(p.svc.vals)
w.ci.widths <- array(NA, dim = c(p.svc.max, n.sims, n.scenarios))
for (i in 1:n.sims) {
  for (j in 1:n.scenarios) {
    indx.1 <- 1:scenario.vals$J[j]
    p.svc.curr <- scenario.vals$p.svc[j]
    for (k in 1:p.svc.curr) {
      w.ci.widths[k, i, j] <- mean(w.high.samples[indx.1, k, i, j]  - w.low.samples[indx.1, k, i, j])
    } #k
  } # j
} # i

w.ci.width.avg <- apply(w.ci.widths, 3, mean, na.rm = TRUE)
plot.w.ci.width <- data.frame(J = factor(scenario.vals$J),
			      n.time = factor(scenario.vals$n.time),
                              p.svc = factor(scenario.vals$p.svc),
                              w.ci.width = w.ci.width.avg)

plot.w.ci.width.2 <- plot.w.ci.width %>%
  group_by(J, p.svc, n.time) %>%
  summarize(w.ci.width = mean(w.ci.width))

# This is part of Table S5
plot.w.ci.width.2 %>%
  arrange(p.svc, J, n.time) %>%
  print(n = nrow(.))

# Correlation of SVC estimates with truth ---------------------------------
p.svc.max <- max(p.svc.vals)
w.cor <- array(NA, dim = c(p.svc.max, n.sims, n.scenarios))
for (i in 1:n.sims) {
  for (j in 1:n.scenarios) {
    indx.1 <- 1:scenario.vals$J[j]
    p.svc.curr <- scenario.vals$p.svc[j]
    for (k in 1:p.svc.curr) {
    w.cor[k, i, j] <- cor(w.mean.samples[indx.1, k, i, j], w.true[indx.1, k, i, j])
    } #k
  } # j
} # i

w.cor.avg <- apply(w.cor, 3, mean, na.rm = TRUE)
plot.w.cor <- data.frame(J = scenario.vals$J,
			 n.time = factor(scenario.vals$n.time),
                         p.svc = factor(scenario.vals$p.svc),
                         w.cor = w.cor.avg)

occ.cor.data <- plot.w.cor %>%
  group_by(p.svc, J, n.time) %>%
  summarize(avg.cor = mean(w.cor)) %>%
  mutate(type = 'Occupancy')

cor.w.plot <- ggplot(occ.cor.data, aes(x = J, y = avg.cor, fill = p.svc)) +
  geom_line(aes(group = p.svc, col = p.svc), size = 1) +
  geom_point(size = 3, pch = 21) +
  facet_wrap(vars(n.time)) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(size = 18)) +
  labs(x = 'Number of sites', y = 'Correlation coefficient', fill = 'Number of SVCs',
       col = 'Number of SVCs', title = '(C) SVC Accuracy')

# Create Figure 2.
ggarrange(psi.bias.plot, psi.ci.width.plot, cor.w.plot, nrow = 3, common.legend = TRUE,
          legend = 'bottom')
ggsave(file = 'figures/Fig-2.pdf', units = 'in', device = 'pdf', width = 13, height = 12)

