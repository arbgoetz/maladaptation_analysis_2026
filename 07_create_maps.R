# Set working directory
dir<-"C:/Users/alexg/Downloads/31937499/"

# Source custom helper functions
source(paste0(dir, "00_convenience_functions.R"))

# load libraries ----------------------------------------------------------

library(mgcv)
library(terra)
library(sf)
library(stars)
library(tidyterra)
library(scales)
library(gstat)
library(concaveman)
library(tidyverse)
library(viridis)
library(cowplot)
library(ragg)
library(prism)
library(ggspatial)
library(cowplot)
library(tidyterra)
conflicts_prefer(terra::mask)

theme_set(theme_bw())


# load source data --------------------------------------------------------


calif<-vect(paste0(dir, "03_spatial_prediction_figures/ca_state/CA_state.shp")) %>%
  project("epsg:4326")
gardens<-read.csv(paste0(dir, "03_spatial_prediction_figures/gardencoords.csv"))
dat_clim_pc<-read.csv(paste0(dir, "03_spatial_prediction_figures/maternal_climate.csv"))
qrange_r<-vect(paste0(dir, "03_spatial_prediction_figures/qlobata_range.shp")) %>% ## calculated in 01b_spatial_calculations.R
  project("epsg:4326")
qrange2<-st_read(paste0(dir, "03_spatial_prediction_figures/qlobata_range2.shp")) ## calculated in 01b_spatial_calculations.R

tmax_sum_current<-rast(paste0(dir, "03_spatial_prediction_figures/tmax_sum.tif")) %>% crop(qrange_r, mask = T)

# Fig 1: map of origin sites ----------------------------------------------

ggplot() + 
  geom_spatraster(data = tmax_sum_current) + 
  #scale_fill_gradient(na.value = "transparent", low = "white", high = muted("red"), name = "Annual maximum\ntemperature\n(30-year normal)") +
  #scale_fill_distiller(palette = "Reds", na.value = "transparent", trans = "reverse") +
  scale_fill_viridis(option = "turbo", na.value = "transparent", name = "Summer maximum\ntemperature\n(°C, 1991-2020\naverage)", begin = 0, end = .75) +
  geom_spatvector(data = qrange_r, fill = "transparent", linewidth = .25, color = "black") +
  geom_sf(data = calif, fill = "transparent", linewidth = .5, color = "black") +
  geom_point(data = dat_clim_pc, aes(x = Longitude, y = Latitude), size = 1) + 
  geom_point(data = gardens, aes(x = Longitude, y = Latitude, color = Site), shape = 15, size = 5, show.legend = FALSE) +
  geom_point(data = gardens, aes(x = Longitude, y = Latitude), color = "black", shape = 22, size = 5) +
  ggplot2::annotate("label", y = c(38.739841, 39.70868), x = c(-120.742112 + .25, -121.782664 + .25), 
                    label = c("Placerville garden", "Chico garden") , color="black",
                    size=4 , hjust = 0) +
  ggplot2::annotate("text", y = -Inf, x = -Inf, label = "Temperature data from PRISM (prism.oregonstate.edu)", hjust = -.075, vjust = -1, size = 3, fontface = "italic") +
  theme_bw() +
  #ggspatial::annotation_scale(location = "tr", style = "ticks") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())
#ggsave("~/git/demographyMS_ecolett/figures/2026/fig1.png", height = 8, width = 8, scale = .75, units = "in", type = "cairo")


# Fig 4: comparison of management scenarios -------------------------------

## load in rasters
mrf_pred_plot<-rast(paste0(dir, "03_spatial_prediction_figures/predictions/mrf_pred_plot.tif"))


## bounding box for map
bbox<-tibble(
  "lat" = c(33.5, 41),
  "lon" = c(-124.482016861124, -117.5)
)

## create map
ggplot() +
  geom_spatvector(data = calif) +
  geom_spatraster_contour_filled(data = mrf_pred_plot, na.rm = F
                                 #, breaks = c(-50, -30, -10, seq(0, 70, by = 10))
                                 , breaks = seq(-60, 60, by = 20)
  ) +
  #geom_spatraster(data = mrf_pred_plot) +
  geom_spatvector(data = qrange_r, fill = "transparent", linewidth = .2, color = "grey1") +
  # geom_spatvector(data = mat_coords) +
  facet_wrap(~lyr) +
  #scale_fill_discrete(palette = colorRampPalette(pal_brewer( palette = "YlGn", direction = 1)(9))(10),  na.value = "transparent") +
  #scale_fill_discrete(palette = colorRampPalette(pal_brewer( palette = "YlGn", direction = 1)(9))(6),  na.value = "transparent") +
  scale_fill_discrete(palette = colorRampPalette(c("#FFFFE5", "#F7FCBD", "#DFF2A7", "#8FD082", "#1B7D40", "#004529")),  na.value = "transparent") +
  geom_spatraster_contour(data = mrf_pred_plot, breaks = 0, color = "grey") +
  #scale_fill_distiller(palette = "YlGn", breaks = seq(-30, 70, by = 10), na.value = "transparent") +
  coord_sf(xlim = bbox$lon, ylim = bbox$lat, expand = FALSE) +
  labs(
    lyr = c(
      "Local acorns", "Future climate match", "Performance-based future climate match"
    ),
    fill = "Relative fitness\nchange (%) by\n2100, RCP 4.5"
  ) + theme(strip.background =element_rect(fill="transparent"))
#ggsave("~/git/demographyMS_ecolett/figures/2026/fig4_v3.png", width = 9.5, height = 4)



# Fig S6: MFF spatial interpolation ---------------------------------------

interp_ifg<-rast(paste0(dir, "03_spatial_prediction_figures/predictions/interp_ifg.tif"))
interp_cso<-rast(paste0(dir, "03_spatial_prediction_figures/predictions/interp_cso.tif"))

cowplot::plot_grid(
  ggplot() +
    geom_sf(data = calif, fill = "white", linewidth = .5, show.legend = FALSE) +
    geom_spatraster(data = mask(interp_ifg, qrange2), aes(fill = var1.pred), show.legend = FALSE) +
    scale_fill_viridis(name = "Relative fitness (2024)", na.value = "transparent", 
                       limits = c(.05, .6)
                       #limits = c(0, 1)
                       ) +
    geom_spatvector(data = qrange_r, fill = "transparent", linewidth = .25, color = "black") +

    
    theme_bw() +
    theme(legend.position = "none",
          legend.margin = margin(0, 0, 0, 0),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) +
    labs(
         title = "Placerville garden"),
  
  ggplot() +
    geom_sf(data = calif, fill = "white", linewidth = .5) +
    geom_spatraster(data = mask(interp_cso, qrange2), aes(fill = var1.pred)) +
    scale_fill_viridis(name = "Performance (2024)", na.value = "transparent", 
                       limits = c(.05, .6)
                       #limits = c(0, 1)
    ) +
    geom_spatvector(data = qrange_r, fill = "transparent", linewidth = .25, color = "black") +

   theme_bw() +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank()) +
    labs(tdiffcat = "",
         title = "Chico garden"),
  rel_heights = c(1, 1), align = "h", axis = "tb",
  ncol = 2
)
#ggsave("~/git/demographyMS_ecolett/figures/figs6.png", height = 5, width = 12, units = "in", type = "cairo")

