# get-WAIC.R: this script extracts the WAIC for each of the five
#             candidate models and saves the results. 
# Author: Jeffrey W. Doser
rm(list = ls())
library(spOccupancy)

n.cand <- 5
waic.df <- data.frame(model = c('Constant', 'Int-land', 'Int-tmax', 
				'SVC', 'Full'), 
		      waic = NA)
load('results/case-study-2/constant-GRSP-model-results.rda')
waic.df$waic[1] <- waicOcc(out)[3]
load('results/case-study-2/int-lulc-GRSP-model-results.rda')
waic.df$waic[2] <- waicOcc(out)[3]
load('results/case-study-2/int-tmax-GRSP-model-results.rda')
waic.df$waic[3] <- waicOcc(out)[3]
load('results/case-study-2/svc-GRSP-model-results.rda')
waic.df$waic[4] <- waicOcc(out)[3]
load('results/case-study-2/full-GRSP-model-results.rda')
waic.df$waic[5] <- waicOcc(out)[3]

save(waic.df, file = 'results/case-study-2/waic-results.rda')
