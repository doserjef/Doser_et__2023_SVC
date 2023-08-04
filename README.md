# Guidelines for the use of spatially-varying coefficients in species distribution models

### In Review

### Jeffrey W. Doser, Marc K&eacute;ry, Sarah P. Saunders, Andrew O. Finley, Brooke L. Bateman, Joanna Grand, Shannon Reault, Aaron S. Weed, Elise F. Zipkin

### spOccupancy Package [Website](https://www.jeffdoser.com/files/spoccupancy-web/) and [Repository](https://github.com/doserjef/spOccupancy/)

### Please contact the first author for questions about the code or data used in the manuscript: Jeffrey W. Doser (doserjef@msu.edu)

---------------------------------

## Abstract

###  Aim 

Species distribution models (SDMs) are increasingly applied across macroscales using detection-nondetection data. Such models typically assume that a single set of regression coefficients can adequately describe species-environment relationships and/or population trends. However, such relationships often show nonlinear and/or spatially-varying patterns that arise from complex interactions with abiotic and biotic processes that operate at different scales. Spatially-varying coefficient (SVC) models can readily account for variability in the effects of environmental covariates. Yet, their use in ecology is relatively scarce due to gaps in understanding the inferential benefits that SVC models can provide compared to simpler frameworks.

### Innovation

Here we demonstrate the inferential benefits of SVC SDMs, with a particular focus on how this approach can be used to generate and test ecological hypotheses regarding the drivers of spatial variability in population trends and species-environment relationships. We illustrate the inferential benefits of SVC SDMs with simulations and two case studies: one that assesses spatially-varying trends of 51 forest bird species in the eastern US over two decades and a second that evaluates spatial variability in the effects of five decades of land cover change on Grasshopper Sparrow (*Ammodramus savannarum*) occurrence across the continental US.

### Main conclusions

We found strong support for SVC SDMs compared to simpler alternatives in both empirical case studies. Factors operating at fine spatial scales, accounted for by the SVCs, were the primary divers of spatial variability in forest bird occurrence trends. Additionally, SVCs revealed complex species-habitat relationships with grassland and cropland area for Grasshopper Sparrow, providing nuanced insights into how future land use change may shape its distribution. These applications display the utility of SVC SDMs to help reveal the environmental factors that drive species distributions across both local and broad scales. We conclude by discussing the potential applications of SVC SDMs in ecology and conservation.

## Repository Directory

All code and resulting model objects were created and saved using spOccupancy v0.7.0. See the [spOccupancy Website](https://www.jeffdoser.com/files/spoccupancy-web/) and [Repository](https://github.com/doserjef/spOccupancy) for multiple vignettes for how to get started with `spOccupancy`.   

### [code/case-study-1](./code/case-study-1)

Contains all code used in the first case study assessing trends of 51 forest bird species across the eastern US

+ `bbs-data-prep.R`: preps raw BBS data for analysis.
+ `bird-life-extract.R`: extracts data from BirdLife International for determination of species ranges.
+ `birdlife-species-ranges.R`: calculates the species ranges based on BirdLife International
+ `format-data-spOccupancy.R`: formats the prepared data into `spOccupancy` format.
+ `get-auc.R`: extracts AUC for each of the candidate models for each species in the assessment of model predictive performance.
+ `get-bbs-2022-data.R`: script to extract BBS data from 2022 for the assessment of out-of-sample predictive performance.
+ `get-prediction-data.R`: generates the spatial locations for prediction across the eastern US.
+ `get-waic.R`: script to extract WAIC from all candidate models for all species.
+ `main-2021.R`: script to run a single-season spatial occupancy model with BBS data from 2021 for use in an assessment of model predictive performance.
+ `main-bcr-stPGOcc.R`: script to run the Strata candidate model, where the strata are Bird Conservation Regions. 
+ `main-stPGOcc.R`: script to run the constant linear candidate model, where the temporal trend is assumed constant over the entire modeled region. 
+ `main-svcTPGOcc.R`: script to run the SVC model, where the temporal trend is allowed to vary spatially. 
+ `main-tmax-stPGOcc.R`: script to run the interaction model with maximum temperature interacting with the temporal trend.
+ `predict-2021-auc.R`: script to predict occupancy probability for each species in 2021 as an assessment of out-of-sample predictive performance.
+ `pred-svcTPGOcc.R`: predicts species-specific occurrence and occurrence trends across the eastern US for the SVC model.
+ `summary.R`: summarizes results from all model fits and generates all figures for this case study.
+ `tmax-data-prep.R`: script to extract the maximum temperature data for use as an interacting variable.

### [code/case-study-2](./code/case-study-2)

Contains all code used in the second case study assessing the impacts of land use change on Grasshopper Sparrow from 1970-2019.

+ `data-prep.R`: prepares the raw BBS data for analysis.
+ `extract-estimates.R`: script to extract the complete coefficient effects across all prediction locations for each of the candidate models.
+ `format-prediction-data.R`: script to format all prediction data together for use with predict functions in `spOccupancy`.
+ `format-spOccupancy-data.R`: script to format all data into the list format required for fitting models in `spOccupancy`.
+ `get-prediction-coords.R`: extracts the prediction coordinates across the US.
+ `get-WAIC.R`: script to calculate the WAIC values for all candidate models.
+ `lulc-data-prep.R`: script to prepare and extract the land cover data from USGS EROS at the BBS routes.
+ `lulc-pred-data-prep.R`: script to prepare and extract the land cover data from USGS EROS across all prediction locations.
+ `main-constant-GRSP.R`: script to fit the constant linear model that assumes constant effects of grassland and cropland cover change on Grasshopper Sparrow occurrence.
+ `main-full-GRSP.R`: script to fit the Full model that includes SVCs and interactive effects of grassland and cropland cover change on Grasshopper Sparrow occurrence.
+ `main-int-lulc-GRSP.R`: script that fits the interaction model that allows the effects of cropland and grassland cover to interact with the average amount of habitat over the 50 year period.
+ `main-int-tmax-GRSP.R`: script that fits the interaction model that allows the effects of cropland and grassland cover to interact with maximum temperature over the 50 year period.  
+ `main-svc-GRSP.R`: script that fits the SVC model for the Grasshopper Sparrow case study.
+ `pred-full.R`: script to predict occurrence probability and the spatially-varying effects of land cover change on Grasshopper sparrow across the continental US using the full model. 
+ `pred-svc.R`: script to predict occurrence probability and the spatially-varying effects of land cover change on Grasshopper sparrow across the continental US using the SVC model. 
+ `summary.R`: summarizes results from all model fits and generates all figures for this case study.

### [code/sims](./code/sims)

Contains all code used in the simulation study in the main text of the manuscript.

+ `main-full-sim.R`: script to run simulations for all candidate models when the species-environment relationship was generated from the "Full" scenario.
+ `main-interaction-sim.R`: script to run simulations for all candidate models when the species-environmetn relationship was assumed to interact with a known variable.
+ `main-linear-sim.R`: script to run simulations for all candidate models when the species-environment was generated as a linear relationship.
+ `main-missing-sim.R`: script to run simulations for all candidate models when the species-environment relationships was generated as an interaction with a missing variable.
+ `main-quadratic-sim.R`: script to run simulations for all candidate models when the species-environment relationship was generated as a a quadratic relationship.
+ `main-strata-sim.R`: script to run simulations for all candidate models when the species-environment relationship was generated as varying across nine pre-specified strata.
+ `summary-sim.R`: script to summarize the results from all simulation scenarios and to generate Figure 1in the manuscript.

### [code/supplemental-sims](./code/supplemental-sims) 

Contains all code used to run the simulation studies in Supplemenetal Information S4.

+ `main-sim-1.R`: runs a simulation study to assess sample size impacts on single-season SVC models.
+ `main-sim-2.R`: runs a simulation study to assess sample size impacts on multi-season SVC models.
+ `main-sim-3.R`: runs a simulation study to assess sample size impacts for assessing spatially-varying population trends.
+ `main-sim-4.R`: runs a simulation study to assess the effect of detection probability and the amount of replication on estimates from an SVC occupancy model.
+ `summary-sim-1.R`: script that summarizes results from the first supplemental simulation study.
+ `summary-sim-2.R`: script that summarizes results from the second supplemental simulation study.
+ `summary-sim-3.R`: script that summarizes results from the third supplemental simulation study.
+ `summary-sim-4.R`: script that summarizes results from the fourth supplemental simulation study.

### [data](./data/)

Contains processed data for the two empirical case studies.

+ `bird-species-table-bateman.csv`: CSV file from Bateman et al. (2020) containing information on species classifications that was used to select the bird community.
+ `case-study-1/`: contains processed data files for Case Study 1 including:
     + `pred-coordinates.rda`: coordinates that form a grid across the eastern US that we used for predicting spatially-varying occurrence trends.
     + `spOcc-bbs-data.rda`: BBS data in the required format for fitting multi-season occupancy models in `spOccupancy`. 
     + `final-spOccupancy-data.rda`: full data list required for fitting the models with `spOccupancy`
     + `spOccupancy-2021-data.rda`: data list used to fit a single-season spatial occupancy model with BBS data from 2021 for use as a predictive performance assessment.
     + `tmax-data.rda`: maximum temperature data calculated from TerraClimate.
+ `case-study-2/`: contains processed data files for Case Study 2 including: 
     + `bbs-data-y-det-covs.rda`: data object consisting of the detection covariates and detection-nondetection data for use in a `spOccupancy` model.
     + `GRSP-pred-data.rda`: data used for prediction across the continental US.
     + `GRSP-spOcc-data.rda`: data list formatted for fitting occupancy models in `spOccupancy`.
     + `pred-coordinates.rda`: prediction coordinates for predicting ocupancy probability across the continental US.


