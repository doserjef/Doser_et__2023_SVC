library(AOI)
library(climateR)
library(sf)
library(raster)
library(rasterVis)
library(elevatr)
library(ggplot2)
library(cowplot)

load("data/bbs-predict-bbs-rich-data.rda")
coords.sf <- st_as_sf(coords, 
		      coords = c('Longitude', 'Latitude'),
		      crs = 4326)

period <- "19812010"
cvars <- c("tmax","tmin","ppt")
nc.vars <- length(cvars)
J <- nrow(coords.sf)
X <- matrix(NA, J, nc.vars + 1)

for (i in 1:nc.vars){
    cdat = getTerraClimNormals(coords.sf, cvars[i], period, 1:12)[[1]]
    plt_dat = extract(cdat, coords.sf)
    if (cvars[i] == "tmax"){
        val = apply(plt_dat[,2:ncol(plt_dat)], 1, max)
    }else if (cvars[i] == "tmin"){
        val = apply(plt_dat[,2:ncol(plt_dat)], 1, min)
    }else if (cvars[i] == "ppt"){
        val = apply(plt_dat[,2:ncol(plt_dat)], 1, sum)
    }else if (cvars[i] %in% c("pet","aet","def")){
        val = apply(plt_dat[,5:9], 1, sum)
    }else {
        val = apply(plt_dat[,5:9], 1, mean)
    }
    X[,i] = val
}
#### Download elevation data
## Source: Amazon Web Services (AWS) Terrain Tiles (https://registry.opendata.aws/terrain-tiles/)
## Citation: Terrain Tiles was accessed on INSERT DATE from https://registry.opendata.aws/terrain-tiles.
## Data accessed on: August 17, 2022

elev = get_elev_point(coords.sf, src = "aws")
X[,nc.vars + 1] = elev$elevation

# Format data for spAbundance ---------------------------------------------
my.proj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"
coords.aea <- coords.sf %>%
  st_transform(my.proj)
coords <- st_coordinates(coords.aea)
data.list <- list(y = y$obs.rich, 
		  covs = data.frame(tmax = X[, 1], 
				    tmin = X[, 2], 
				    ppt = X[, 3],
				    elev = X[, 4], 
				    day = day, 
				    obs = obs),
		  coords = coords)
save(data.list, file = 'data/bbs-predict/spAbundance-data.rda')
