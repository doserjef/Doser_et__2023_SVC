# summary-sim-1.R: this script summarizes the first supplemental simulation study from
#                  the manuscript in which we assessed the reliability of 
#                  estimates from an SVC occupancy model under varying 
#                  numbers of spatial locations and varying amounts of 
#                  spatially-varying coefficients. 
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(viridis)
library(devtools)
library(ggthemes)
library(ggpubr)
library(pals)
library(spOccupancy)

# Functions ---------------------------------------------------------------
logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

# Occupancy simulations ---------------------------------------------------
# Load results ------------------------------------------------------------
load("results/supplemental-sims/sim-1-results.rda")
# Simulation scenarios
n.sims <- ncol(psi.mean.samples)
J.vals <- c(200, 400, 800, 1200, 1600, 2000, 4000, 6000)
p.svc.vals <- c(2, 3, 4)
# Total number of simulation scenarios
n.scenarios <- length(J.vals) * length(p.svc.vals)
scenario.vals <- expand.grid(J = J.vals, p.svc = p.svc.vals)
scenario.vals$svc.plot <- scenario.vals$p.svc - 1

# Calculate occurrence probability bias -----------------------------------
psi.abs.bias <- array(NA, dim = c(n.sims, n.scenarios))
for (i in 1:n.sims) {
  for (j in 1:n.scenarios) {
    # psi.abs.bias[i, j] <- mean(abs(psi.mean.samples[1:scenario.vals$J[j], i, j] - psi.true[1:scenario.vals$J[j], i, j]) / psi.true[1:scenario.vals$J[j], i, j])
    psi.abs.bias[i, j] <- mean(abs(psi.mean.samples[1:scenario.vals$J[j], i, j] - psi.true[1:scenario.vals$J[j], i, j]))
  } # j
} # i

psi.abs.bias.avg <- apply(psi.abs.bias, 2, mean)
plot.psi.abs.bias <- data.frame(J = scenario.vals$J,
				p.svc = factor(scenario.vals$svc.plot),
				psi.abs.bias = psi.abs.bias.avg)
plot.psi.abs.bias.2 <- plot.psi.abs.bias %>%
  group_by(J, p.svc) %>%
  summarize(psi.abs.bias = mean(psi.abs.bias))

psi.bias.plot <- ggplot(plot.psi.abs.bias.2, aes(x = J, y = psi.abs.bias, fill = p.svc)) + 
  geom_line(aes(group = p.svc, col = p.svc), size = 1) +
  geom_point(size = 3, pch = 21) + 
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme_bw(base_size = 18) + 
  theme(plot.title = element_text(size = 18), 
	text = element_text(family = 'LM Roman 10')) +
  labs(x = 'Number of sites', y = 'Absolute bias', fill = 'Number of SVCs', 
       col = 'Number of SVCs', title = '(A) Occupancy probability bias')
psi.bias.plot

# Averaged across number of sites.
plot.psi.abs.bias %>%
  group_by(J) %>%
  summarize(avg.bias = mean(psi.abs.bias))
# Averaged across number of SVCs
plot.psi.abs.bias %>%
  group_by(p.svc) %>%
  summarize(avg.bias = mean(psi.abs.bias))
# Grouped by SVC and Number of sites. This is part of Supp Info S4 Table 1.
plot.psi.abs.bias %>%
  group_by(p.svc, J) %>%
  summarize(avg.bias = mean(psi.abs.bias)) %>%
  print(n = nrow(.))

# Calculate average coverage rates ----------------------------------------
psi.coverage.rates <- array(NA, dim = c(n.sims, n.scenarios))
for (i in 1:n.sims) {
  for (j in 1:n.scenarios) {
    indx <- 1:scenario.vals$J[j]
    psi.coverage.rates[i, j] <- mean(ifelse((psi.true[indx, i, j] > psi.low.samples[indx, i, j]) & 
						  (psi.true[indx, i, j] < psi.high.samples[indx, i, j]), 1, 0))
  } # j
} # i

psi.coverage.avg <- apply(psi.coverage.rates, 2, mean)
plot.psi.coverage <- data.frame(J = factor(scenario.vals$J),
				p.svc = factor(scenario.vals$svc.plot),
				psi.coverage = psi.coverage.avg)

plot.psi.coverage.2 <- plot.psi.coverage %>%
  group_by(J, p.svc) %>%
  summarize(psi.coverage = mean(psi.coverage))

# Averaged across number of sites.
plot.psi.coverage %>%
  group_by(J) %>%
  summarize(avg.cov = mean(psi.coverage))
# Averaged across number of SVCs
plot.psi.coverage %>%
  group_by(p.svc) %>%
  summarize(avg.cov = mean(psi.coverage))
# Grouped by SVC and Number of sites. This is part of Supp Info S4 Table 1
plot.psi.coverage %>%
  group_by(p.svc, J) %>%
  summarize(avg.cov = mean(psi.coverage)) %>%
  print(n = nrow(.))

# Calculate average CI width ----------------------------------------------
psi.ci.width <- array(NA, dim = c(n.sims, n.scenarios))
for (i in 1:n.sims) {
  for (j in 1:n.scenarios) {
    indx <- 1:scenario.vals$J[j]
    psi.ci.width[i, j] <- mean(psi.high.samples[1:scenario.vals$J[j], i, j] - psi.low.samples[1:scenario.vals$J[j], i, j])
  } # j
} # i

psi.ci.width.avg <- apply(psi.ci.width, 2, mean)
plot.psi.ci.width <- data.frame(J = scenario.vals$J,
				p.svc = factor(scenario.vals$svc.plot),
				psi.ci.width = psi.ci.width.avg)

plot.psi.ci.width.2 <- plot.psi.ci.width %>%
  group_by(J, p.svc) %>%
  summarize(psi.ci.width = mean(psi.ci.width))

psi.ci.width.plot <- ggplot(plot.psi.ci.width.2, aes(x = J, y = psi.ci.width, fill = p.svc)) + 
  geom_line(aes(group = p.svc, col = p.svc), size = 1) +
  geom_point(size = 3, pch = 21) + 
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme_bw(base_size = 18) + 
  theme(plot.title = element_text(size = 18), 
	text = element_text(family = 'LM Roman 10')) +
  labs(x = 'Number of sites', y = '95% CI Width', fill = 'Number of SVCs', 
       col = 'Number of SVCs', title = '(B) Occupancy probability uncertainty')
psi.ci.width.plot

# Averaged across number of sites.
plot.psi.ci.width %>%
  group_by(J) %>%
  summarize(avg.width = mean(psi.ci.width))
# Averaged across number of SVCs
plot.psi.ci.width %>%
  group_by(p.svc) %>%
  summarize(avg.width = mean(psi.ci.width))
# Grouped by SVC and Number of sites. This is part of Supp Info S4 Table 1.
plot.psi.ci.width %>%
  group_by(p.svc, J) %>%
  summarize(avg.width = mean(psi.ci.width)) %>%
  print(n = nrow(.))

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
                              p.svc = factor(scenario.vals$svc.plot),
                              w.abs.bias = w.abs.bias.avg)

plot.w.abs.bias.2 <- plot.w.abs.bias %>%
  group_by(J, p.svc) %>%
  summarize(w.abs.bias = mean(w.abs.bias))

# Averaged across number of sites.
plot.w.abs.bias %>%
  group_by(J) %>%
  summarize(avg.bias = mean(w.abs.bias))
# Averaged across number of SVCs
plot.w.abs.bias %>%
  group_by(p.svc) %>%
  summarize(avg.bias = mean(w.abs.bias))
# Grouped by SVC and Number of sites. This is part of Supp Info S4 Table 2 
plot.w.abs.bias %>%
  group_by(p.svc, J) %>%
  summarize(avg.bias = mean(w.abs.bias)) %>%
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
                              p.svc = factor(scenario.vals$svc.plot),
                              w.coverage = w.coverage.avg)

plot.w.coverage.2 <- plot.w.coverage %>%
  group_by(J, p.svc) %>%
  summarize(w.coverage = mean(w.coverage))
# This is part of Supp Info S4 Table 2
plot.w.coverage.2 %>%
  arrange(p.svc, J) %>%
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
                              p.svc = factor(scenario.vals$svc.plot),
                              w.ci.width = w.ci.width.avg)

plot.w.ci.width.2 <- plot.w.ci.width %>%
  group_by(J, p.svc) %>%
  summarize(w.ci.width = mean(w.ci.width))
# This is part of Supp Info S4 Table 2
plot.w.ci.width.2 %>%
  arrange(p.svc, J) %>%
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
                              p.svc = factor(scenario.vals$svc.plot),
                              w.cor = w.cor.avg)

# Averaged across number of sites.
plot.w.cor %>%
  group_by(J) %>%
  summarize(avg.cor = mean(w.cor))
# Averaged across number of SVCs
plot.w.cor %>%
  group_by(p.svc) %>%
  summarize(avg.cor = mean(w.cor))
# Grouped by SVC and Number of sites
occ.cor.data <- plot.w.cor %>%
  group_by(p.svc, J) %>%
  summarize(avg.cor = mean(w.cor)) %>%
  mutate(type = 'Occupancy')

cor.plot.data <- plot.w.cor %>%
  group_by(J) %>%
  summarize(avg.cor = mean(w.cor))


cor.w.plot <- ggplot(occ.cor.data, aes(x = J, y = avg.cor, fill = p.svc)) + 
  geom_line(aes(group = p.svc, col = p.svc), size = 1) +
  geom_point(size = 3, pch = 21) + 
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme_bw(base_size = 18) + 
  theme(plot.title = element_text(size = 18), 
	text = element_text(family = 'LM Roman 10')) +
  labs(x = 'Number of sites', y = 'Correlation coefficient', fill = 'Number of SVCs', 
       col = 'Number of SVCs', title = '(C) SVC Accuracy')

# Supplemental Information S4 Figure 1
ggarrange(psi.bias.plot, psi.ci.width.plot, cor.w.plot, nrow = 1, common.legend = TRUE, 
          legend = 'bottom')
ggsave(file = 'figures/Supp-Info-S4-Figure-1.png', units = 'in', device = 'png', width = 15, height = 6, 
       bg = 'white')
