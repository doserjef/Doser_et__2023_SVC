rm(list = ls())
library(spOccupancy)
library(sf)
library(viridis)
library(pals)
library(ggpubr)
library(tidyverse)

# Read in results file ----------------------------------------------------
load('results/tree/stage-1-sugar-maple-full-samples.rda')

# Read in data ------------------------------------------------------------
load('data/tree/stage-1-data.rda')

svi.means <- apply(svc.samples[[1]], 2, mean)
psi.meds <- apply(psi.samples, 2, median)
psi.ci.width <- apply(psi.samples, 2, quantile, 0.975) - apply(psi.samples, 2, quantile, 0.025)

plot.df <- data.frame(psi.med = psi.meds,
		      psi.ci.width = psi.ci.width,
		      x = data.list$coords[, 1],
		      y = data.list$coords[, 2])

plot.sf <- st_as_sf(plot.df,
		    coords = c("x", "y"),
		    crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs")
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
usa <- usa %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs")


# Make a summary plot -----------------------------------------------------
psi.plot <- ggplot(plot.sf) +
  geom_sf(aes(col = psi.med), size = 0.1) +
  scale_color_gradientn("", colors = ocean.tempo(1000), limits = c(0, 1),
                         guide = guide_colourbar(title.position="top", reverse = FALSE),
			 na.value = NA) +
  geom_sf(data = usa, fill = NA, color=alpha("grey", 0.75)) +
  theme_bw(base_size = 14) +
  theme(text = element_text(family="LM Roman 10", size=10),
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
	# axis.text.y = element_text(size = 18),
        legend.key.size = unit(0.5, 'cm'),
        legend.direction="horizontal",
        legend.box="vertical",
        legend.position=c(0.175,0.12),
        legend.background = element_rect(fill = NA)) +
    labs(title = 'Sugar Maple Occurrence Probability')
psi.ci.plot <- ggplot(plot.sf) +
  geom_sf(aes(col = psi.ci.width), size = 0.1) +
  scale_color_gradientn("", colors = ocean.tempo(1000), limits = c(0, 1),
                         guide = guide_colourbar(title.position="top", reverse = FALSE),
			 na.value = NA) +
  geom_sf(data = usa, fill = NA, color=alpha("grey", 0.75)) +
  theme_bw(base_size = 14) +
  theme(text = element_text(family="LM Roman 10", size=10),
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
	# axis.text.y = element_text(size = 18),
        legend.key.size = unit(0.5, 'cm'),
        legend.direction="horizontal",
        legend.box="vertical",
        legend.position=c(0.175,0.12),
        legend.background = element_rect(fill = NA)) +
    labs(title = 'Sugar Maple 95% Credible Interval Width')
plot(ggarrange(psi.plot, psi.ci.plot))

