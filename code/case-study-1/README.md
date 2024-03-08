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
