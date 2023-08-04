rm(list = ls())
library(spAbundance)

load('data/bbs-predict/spAbundance-data.rda')

out <- svcAbund(formula = ~ scale(tmax) + scale(ppt), 
		data = data.list,
		n.batch = 400,
		batch.length = 25,
		NNGP = TRUE, 
		n.neighbors = 5,
		svc.cols = c(1, 2, 3),
		n.report = 1)
	

# Exploratory variogram analysis
library(MBA)
library(fields)
library(geoR)
par(mfrow=c(1,2))
out.lm <- lm(data.list$y ~ scale(data.list$covs$tmax) + scale(data.list$covs$ppt))
curr.res <- out.lm$residuals
vario.1.raw <- variog(coords = data.list$coords, data = data.list$y)
par(mfrow = c(1, 2))
plot(vario.1.raw, pch=16)
vario.1.res <- variog(coords = data.list$coords, data = curr.res)
plot(vario.1.res, pch = 16)
par(mfrow = c(1, 1))

