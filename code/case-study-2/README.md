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
