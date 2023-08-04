library(MBA)
library(fields)
library(geoR)
par(mfrow=c(1,2))

load("data/tree/vermont-bio-coords.rda")
vario.1.raw <- variog(coords = coords, data = sqrt(y))
# par(mfrow = c(1, 2))
plot(vario.1.raw, pch=16)

