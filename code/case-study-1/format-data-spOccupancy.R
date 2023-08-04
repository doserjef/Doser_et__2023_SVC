# format-data-spOccupancy.R: this script does some final formatting of the 
#                            data for use in spOccupancy. Note that the data 
#                            formatted for use with a multi-species, multi-season
#                            SVC occupancy model, and then the object is subset 
#                            individually for each species. Note this script will
#                            not run successfully with the files available on GitHub, 
#                            as the BirdLife data is proprietary and requires a 
#                            data sharing agreement to access
# Author: Jeffrey W. Doser
rm(list = ls())
library(spOccupancy)

# BBS data and detection covariates
load("data/case-study-1/spOcc-bbs-data.rda")
# Processed bird-life data
load("data/case-study-1/bird-life-processed.rda")
# Max temperature data from PRISM.
load("data/case-study-1/tmax-data.rda")

# Only keep species with >50% overlap with the eastern US study region. 
keep.sp <- which(percent.area > 0.5)

data.list$y <- data.list$y[keep.sp, , , ]
data.list$range.ind <- range.ind[keep.sp, ]

# Add in replicate as a covariate effect on detection
rep.val <- array(NA, dim = dim(data.list$y[1, , , ]))
for (i in 1:dim(rep.val)[3]) {
  rep.val[, , i] <- i
}
data.list$det.covs$rep.val <- rep.val

data.list$occ.covs$tmax <- tmax

# Save to final object ready for use in spOccupancy
save(data.list, file = 'data/case-study-1/final-spOccupancy-data.rda')
