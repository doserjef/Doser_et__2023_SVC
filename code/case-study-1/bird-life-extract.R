# bird-life-extract.R: this script extracts data from Bird Life International
#                      for determination of species ranges. Note the BirdLife
#                      data must be obtained via a data sharing agreement, and 
#                      so this will not run using the available files on GitHub.
#                      This code extracts data for a community of eastern forest bird species.
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(coda)
library(spOccupancy)
library(sf)
library(lwgeom)

# Read in species name data for eastern forest bird community -------------
load("data/case-study-1/spOcc-bbs-data.rda")
sp.names <- dimnames(data.list$y)[[1]]
comm.group.dat <- read.csv("data/bird-species-table-bateman.csv")
scientific.names <- comm.group.dat %>%
  filter(Code %in% sp.names) %>%
  select(Code, name = Scientific.Name) %>%
  mutate(name = tolower(name)) %>%
  unique()

# Three species with differing classifications that you need to adjust
# Nashville Warbler
scientific.names$name[which(scientific.names$Code == 'NAWA')] <- "leiothlypis ruficapilla"
# Pileaated Woodpecker
scientific.names$name[which(scientific.names$Code == 'PIWO')] <- "hylatomus pileatus"
# Red-Cockaded Woodpecker 
scientific.names$name[which(scientific.names$Code == 'RCWO')] <- "leuconotopicus borealis"

# Read in bird-life data --------------------------------------------------
# What layers are in the GDB? 
st_layers(dsn = 'data/case-study-1/BOTW/BOTW.gdb')
# Accipiter brevipes
# Takes a few minutes
# bird.life <- st_read(dsn = 'data/BOTW/BOTW.gdb', query = "SELECT * FROM All_Species limit 2000" )
bird.life <- st_read(dsn = 'data/BOTW/BOTW.gdb', 
		     query = "SELECT * 
		              FROM All_Species
		              WHERE sci_name IN ('Antrostomus carolinensis', 'Archilochus colubris', 
						 'Baeolophus bicolor', 'Buteo lineatus', 
						 'Buteo platypterus', 'Cardinalis cardinalis', 
						 'Catharus fuscescens', 'Chaetura pelagica', 
						 'Coccyzus americanus', 'Coccyzus erythropthalmus', 
						 'Contopus virens', 'Cyanocitta cristata', 
						 'Dumetella carolinensis', 'Empidonax virescens', 
						 'Falco columbarius', 'Geothlypis formosa',
						 'Helmitheros vermivorum', 'Hylatomus pileatus', 
						 'Hylocichla mustelina', 'Icteria virens', 
						 'Icterus galbula', 'Icterus spurius', 
						 'Ictinia mississippiensis', 'Leiothlypis ruficapilla', 
						 'Leuconotopicus borealis', 'Limnothlypis swainsonii', 
						 'Megascops asio', 'Melanerpes carolinus', 
						 'Melanerpes erythrocephalus', 'Mniotilta varia', 
						 'Myiarchus crinitus', 'Parkesia motacilla', 
						 'Passerina cyanea', 'Peucaea aestivalis', 
						 'Pheucticus ludovicianus', 'Pipilo erythrophthalmus', 
						 'Piranga olivacea', 'Piranga rubra', 
						 'Poecile carolinensis', 'Protonotaria citrea', 
						 'Quiscalus quiscula', 'Sayornis phoebe', 
						 'Scolopax minor', 'Seiurus aurocapilla', 
						 'Setophaga americana', 'Setophaga caerulescens', 
						 'Setophaga cerulea', 'Setophaga citrina', 
						 'Setophaga discolor', 'Setophaga dominica', 
						 'Setophaga pensylvanica', 'Setophaga pinus', 
						 'Setophaga virens', 'Sialia sialis', 
						 'Sitta carolinensis', 'Sitta pusilla', 
						 'Sphyrapicus varius', 'Spizella pusilla', 
						 'Thryothorus ludovicianus', 'Toxostoma rufum', 
						 'Vermivora chrysoptera', 'Vermivora cyanoptera', 
						 'Vireo flavifrons', 'Vireo griseus', 
						 'Vireo olivaceus', 'Vireo solitarius')")
save(bird.life, file = 'data/case-study-1/bird-life-data.rda')

