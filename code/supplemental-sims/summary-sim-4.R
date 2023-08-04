# summary-sim-4.R: this script summarizes results from the fourth supplemental 
#                  simulation study that shows how estimates from SVC 
#                  occupancy models varies across differing levels of replication
#                  and different detection probabilities. 
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(viridis)
library(devtools)
library(ggthemes)
library(spOccupancy)

logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

# SVC occupancy model results ---------------------------------------------
# Load results ------------------------------------------------------------
load("results/supplemental-sims/sim-4-results.rda")
# Number of sites
J <- dim(psi.mean.samples)[1]
# Number of simulations for each scenario
n.sims <- dim(psi.mean.samples)[2]
# Number of scenarios
n.scenarios <- dim(psi.mean.samples)[3]

# Simulated values that vary from simulation to simulation. 
n.rep.vals <- c(2, 3, 4)
alpha.vals <- c(logit(0.1), logit(0.3), logit(0.5), logit(0.7))
scenario.vals <- expand.grid(n.rep = n.rep.vals, alpha = alpha.vals)

# Absolute bias -----------------------------------------------------------
# Occurrence --------------------------
psi.abs.bias <- array(NA, dim = c(n.sims, n.scenarios))
for (i in 1:n.sims) {
  for (j in 1:n.scenarios) {
    psi.abs.bias[i, j] <- mean(abs(psi.mean.samples[, i, j] - psi.true[, i, j]))
  } # j
} # i

psi.abs.bias.avg <- apply(psi.abs.bias, 2, mean)
plot.psi.abs.bias <- data.frame(n.rep = scenario.vals$n.rep,
				alpha = factor(logit.inv(scenario.vals$alpha)),
				psi.abs.bias = psi.abs.bias.avg)
# This is part of Supp Info S4 Table 6
plot.psi.abs.bias %>%
  arrange(n.rep, alpha)
# SVC ---------------------------------
w.abs.bias <- array(NA, dim = c(n.sims, n.scenarios))
for (i in 1:n.sims) {
  for (j in 1:n.scenarios) {
    w.abs.bias[i, j] <- mean(abs(w.mean.samples[, 2, i, j] - w.true[, 2, i, j]))
  } # j
} # i

w.abs.bias.avg <- apply(w.abs.bias, 2, mean)
plot.w.abs.bias <- data.frame(n.rep = scenario.vals$n.rep,
				alpha = factor(logit.inv(scenario.vals$alpha)),
				w.abs.bias = w.abs.bias.avg)
# This is part of Supp Info S4 Table 7
plot.w.abs.bias %>%
  arrange(n.rep, alpha)

# Calculate average CI width ----------------------------------------------
# Occurrence --------------------------
psi.ci.width <- array(NA, dim = c(n.sims, n.scenarios))
for (i in 1:n.sims) {
  for (j in 1:n.scenarios) {
    psi.ci.width[i, j] <- mean(psi.high.samples[, i, j] - psi.low.samples[, i, j])
  } # j
} # i

psi.ci.width.avg <- apply(psi.ci.width, 2, mean)
plot.psi.ci.width <- data.frame(n.rep = scenario.vals$n.rep,
				alpha = factor(logit.inv(scenario.vals$alpha)),
				psi.ci.width = psi.ci.width.avg)
# This is part of Supp Info S4 Table 6
plot.psi.ci.width %>%
  arrange(n.rep, alpha)
# SVC ---------------------------------
w.ci.width <- array(NA, dim = c(n.sims, n.scenarios))
for (i in 1:n.sims) {
  for (j in 1:n.scenarios) {
    w.ci.width[i, j] <- mean(w.high.samples[, 2, i, j] - w.low.samples[, 2, i, j])
  } # j
} # i

w.ci.width.avg <- apply(w.ci.width, 2, mean)
plot.w.ci.width <- data.frame(n.rep = scenario.vals$n.rep,
				alpha = factor(logit.inv(scenario.vals$alpha)),
				w.ci.width = w.ci.width.avg)
# This is part of Supp Info S4 Table 7
plot.w.ci.width %>%
  arrange(n.rep, alpha)
