rm(list = ls())
library(spBayes)

data(PM10.dat)
PM10.mod <- PM10.dat[!is.na(PM10.dat$pm10.obs), ]
PM10.pred <- PM10.dat[is.na(PM10.dat$pm10.obs), ]

d.max <- max(iDist(PM10.mod[, c('x.coord', 'y.coord')]))

r <- 2

priors <- list('phi.Unif' = list(rep(3 / (.75 * d.max), r), rep(3 / (.001 * d.max), r)), 
	       'sigma.sq.IG' = list(rep(2, r), rep(1, r)), 
	       'tau.sq.IG' = c(2, 1))

starting <- list('phi' = rep(3 / (.1 * d.max), r), 'sigma.sq' = rep(1, r), 'tau.sq' = 1)
tuning <- list(phi = rep(0.1, r), sigma.sq = rep(0.05, r), tau.sq = 0.1)
n.samples <- 40000
out <- spSVC(pm10.obs ~ pm10.ctm, coords = c('x.coord', 'y.coord'), 
	     data = PM10.mod, starting = starting, svc.cols = c(1, 2),
	     tuning = tuning, priors = priors, cov.model = 'exponential', 
	     n.samples = n.samples, n.report = 1000, n.omp.threads = 4)

out.recover <- spRecover(out, start = floor(0.75 * n.samples), thin = 2, 
			 n.omp.threads = 4)

summary(out.recover$p.beta.recover.samples)
summary(out.recover$p.theta.recover.samples)

y.rep.means <- apply(out.recover$p.y.samples, 1, mean)

plot(PM10.mod$pm10.obs, y.rep.means, pch = 19)
abline(0, 1)
