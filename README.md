# Guidelines for the use of spatially-varying coefficients in species distribution models

### In Review

### Jeffrey W. Doser, Marc K&eacute;ry, Andrew O. Finley, Sarah P. Saunders, Aaron S. Weed, Elise F. Zipkin

### spOccupancy Package [Website](https://www.jeffdoser.com/files/spoccupancy-web/) and [Repository](https://github.com/doserjef/spOccupancy/)

### Please contact the first author for questions about the code or data used in the manuscript: Jeffrey W. Doser (doserjef@msu.edu)

---------------------------------

## Abstract

###  Aim 

Species distribution models (SDMs) are increasingly applied across macroscales using widely available detection-nondetection data sources. However, assumptions of stationarity in species-environment relationships or population trends inherent to most SDM techniques are frequently violated at broad spatial scales. Bayesian spatially-varying coefficient (SVC) models can readily account for nonstationarity, yet their use is relatively scarce, due, in part, to a gap in understanding both the data requirements needed to fit SVC SDMs, as well as the inferential benefits of applying a more complex modeling framework.

### Innovation

Using simulations, we present guidelines and recommendations for fitting single-season and multi-season SVC SDMs. We display the inferential benefits of SVC SDMs using an empirical case study assessing spatially-varying trends of 51 forest birds in the eastern US from 2000-2019. We provide user-friendly and computationally efficient software to fit SVC SDMs in the spOccupancy R package. 

### Main conclusions

While all datasets are unique, we recommend a minimum sample size of ${\sim}500$ spatial locations when fitting single-season SVC SDMs, while for multi-season SVC SDMs, ${\sim}100$ sites is adequate for even moderate amounts of temporal replication (e.g., 5 years). Within our case study, we found 88\% (45 of 51) of species had strong support for spatially-varying occurrence trends. Further, SVC SDMs revealed spatial patterns in occurrence trends that were not evident in simpler models that assumed a constant trend or separate trends across distinct ecoregions. We suggest five guidelines: (1) only fit single-season SVC SDMs with more than ${\sim}500$ sites; (2) consider using informative priors on spatial parameters to improve spatial process estimates; (3) use data from multiple seasons if available; (4) use model selection to compare SVC SDMs with simpler alternatives; and (5) develop simulations to assess the reliability of inferences. These guidelines provide a comprehensive foundation for using SVC SDMs to evaluate the presence and impact of nonstationary environmental factors that drive species distributions at macroscales.

## Repository Directory

All code and resulting model objects were created and saved using spOccupancy v0.5.2. See the [spOccupancy Website](https://www.jeffdoser.com/files/spoccupancy-web/) and [Repository](https://github.com/doserjef/spOccupancy) for multiple vignettes for how to get started with `spOccupancy`.   

### [code/case-study](./code/case-study)

Contains all code used in the BBS case study

+ `bbs-data-prep.R`: preps raw BBS data for analysis.
+ `bird-life-extract.R`: extracts data from BirdLife International for determination of species ranges.
+ `birdlife-species-ranges.R`: calculates the species ranges based on BirdLife International
+ `get-prediction-data.R`: generates the spatial locations for prediction across the eastern US.
+ `main-bcr-stPGOcc.R`: runs the "BCR" model, which is a multi-season spatial occupancy model where the trend is allowed to vary across BCRs. 
+ `main-stPGOcc.R`: runs the "Constant" model, which is a multi-season spatial occupancy model where the trend is assumed constant across the eastern US.
+ `main-svcTPGOcc.R`: runs the "SVC" model, which is a multi-season spatially-varying coefficients occupancy model where the trend is allowed to vary continuously across the eastern US.
+ `pred-svcTPGOcc.R`: predicts species-specific occurrence and occurrence trends across the eastern US for the SVC model.
+ `summary.R`: summarizes results from all model fits and generates all figures in the manuscript.

### [code/sims](./code/sims)

Contains all code used in the simulation study. 

+ `main-sim-1.R`: runs a simulation study that serves as a proof of concept for SVC occupancy models.
+ `main-sim-1-glm.R`: runs a simulation study that serves as a proof of concept for SVC GLMs without imperfect detection.
+ `main-sim-2.R`: runs a simulation study to assess sample size requirements for single-season SVC models.
+ `main-sim-3.R`: runs a simulation study to assess sample size requirements for multi-season SVC models.
+ `main-sim-4.R`: runs a simulation study to assess sample size requirements for assessing spatially-varying population trends.
+ `main-sim-5.R`: runs a simulation study to assess the effect of detection probability and the amount of replication on estimates from an SVC occupancy model.
+ `summary-sim-1.R`: script that summarizes results from the first case study.
+ `summary-sim-2.R`: script that summarizes results from the second case study.
+ `summary-sim-3.R`: script that summarizes results from the third case study.
+ `summary-sim-4.R`: script that summarizes results from the fourth case study.
+ `summary-sim-5.R`: script that summarizes results from the fifth case study.

### [data](./data/)

Contains data used for BBS case study.

+ `bird-life-processed.rda`: processed BirdLife International data for use in determining the area over which to model each species.
+ `bird-species-table-bateman.csv`: CSV file from Bateman et al. (2020) containing information on species classifications that was used to select the bird community.
+ `pred-coordinates.rda`: coordinates that form a grid across the eastern US that we used for predicting spatially-varying occurrence trends.
+ `spOcc-bbs-data.rda`: BBS data in the required format for fitting multi-season occupancy models in `spOccupancy`. 

