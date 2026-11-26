# Set working directory
dir<-"C:/Users/alexg/Downloads/31937499/"

# Source custom helper functions
source(paste0(dir, "01_convenience_functions.R"))

# load libraries ----------------------------------------------------------

library(mgcv)        # Generalized Additive Models
library(terra)       # Spatial raster data manipulation
library(sf)          # Simple features for vector spatial data
library(stars)       # Spatiotemporal arrays for raster data
library(scales)      # Scale functions for visualization
library(gstat)       # Geostatistical modeling
library(concaveman)  # Concave hull polygons
library(tidyverse)   # Data manipulation suite
library(viridis)     # Color palettes
library(ragg)        # Graphics device
library(ggspatial)   # Spatial visualization with ggplot2
library(tidyterra)   # Tidy methods for terra objects
library(lme4)        # Mixed-effects models
library(gratia)      # GAM visualization
library(automap)     # Automated spatial interpolation
conflicted::conflicts_prefer(terra::resample)
conflicted::conflicts_prefer(dplyr::filter)

# Load required packages for cross-validation
library(modelr)
library(mgcViz)
library(DHARMa)

theme_set(theme_bw())

# load data ---------------------------------------------------------------

## Hypothetical range polygon for Quercus lobata (valley oak)
qrange_r <- vect(paste0(dir, "03_spatial_prediction_figures/qlobata_range.shp")) %>%
  project("epsg:4326")  # WGS84 coordinate system

## Relative fitness empirical data from common garden experiment
##2024
relfit24 <- read.csv(paste0(dir, "mff24_input.csv")) %>%
  mutate(
    Accession = as.factor(Accession),
    Locality = as.factor(Locality),
    Year = as.factor(Year),
    Site = as.factor(Site)
  )
##2017
relfit17<-read.csv(paste0(dir, "03_spatial_prediction_figures/mff17_input.csv")) %>%
  mutate(
    Accession = as.factor(Accession),
    Locality = as.factor(Locality),
    Year = as.factor(Year),
    Site = as.factor(Site)
  )

## Climate data from maternal (origin) sites
dat_clim_pc <- read.csv(paste0(dir, "03_spatial_prediction_figures/maternal_climate.csv")) %>%
  mutate(Accession = as.factor(Accession))

# Create spatial points object from maternal coordinates
mat_coords <- st_as_sf(dplyr::select(dat_clim_pc, Longitude, Latitude, Accession), 
                       coords = c("Longitude", "Latitude"), 
                       crs = 4326)

## Climate rasters - historical (1961-1990), current, and future projections
# All cropped to Q. lobata range
tmax_sum_61 <- rast(paste0(dir, "03_spatial_prediction_figures/tmax_sum_61.tif")) %>% 
  crop(qrange_r, mask = T)  # Historical summer max temperature
tmax_sum_current <- rast(paste0(dir, "03_spatial_prediction_figures/tmax_sum.tif")) %>% 
  crop(qrange_r, mask = T)  # Current summer max temperature
tmaxpred_sum_45 <- rast(paste0(dir, "03_spatial_prediction_figures/tmaxpred_sum_45.tif")) %>% 
  crop(qrange_r, mask = T)  # Future (2045) summer max temperature
bio11_61 <- rast(paste0(dir, "03_spatial_prediction_figures/bio11_61.tif")) %>% 
  crop(qrange_r, mask = T)  # Historical winter temperature
bio11_current <- rast(paste0(dir, "03_spatial_prediction_figures/bio11_current.tif")) %>% 
  crop(qrange_r, mask = T)  # Current winter temperature
bio18_61 <- rast(paste0(dir, "03_spatial_prediction_figures/bio18_61.tif")) %>% 
  crop(qrange_r, mask = T)  # Historical warm quarter precipitation
bio18_current <- rast(paste0(dir, "03_spatial_prediction_figures/bio18_current.tif")) %>% 
  crop(qrange_r, mask = T)  # Current warm quarter precipitation
ca_dem <- rast(paste0(dir, "03_spatial_prediction_figures/ca_dem.tif")) %>% 
  crop(qrange_r, mask = T)  # Digital elevation model
Latitude <- rast(paste0(dir, "03_spatial_prediction_figures/Latitude.tif")) %>% 
  crop(qrange_r, mask = T)  # Latitude raster
Longitude <- rast(paste0(dir, "03_spatial_prediction_figures/Longitude.tif")) %>% 
  crop(qrange_r, mask = T)  # Longitude raster
tdiff_rast_ca <- rast(paste0(dir, "03_spatial_prediction_figures/tdiff_rast_ca.tif")) %>% 
  crop(qrange_r, mask = T)  # Temperature difference: future - historical

# fit phenotype-environment association model on point data -----------------------------

# Extract historical climate values at maternal origin sites
relfit24_61 <- bind_cols(terra::extract(c(tmax_sum_61, bio11_61, bio18_61), mat_coords)) %>%
  rename("Accession" = ID,
         "tmax_sum_61" = mean,
         "bio11_61" = `bigfile[, varname]...3`,
         "bio18_61" = `bigfile[, varname]...4`) %>%
  mutate(Accession = as.factor(Accession)) %>%
  right_join(relfit24, by = "Accession")

# Phenotype-Environment Association (PEA) model
# Models relative fitness as a function of climate and geography at origin sites
# This establishes the "baseline" fitness distribution across the landscape
gam_dist <- bam(
  mrf ~ 
    s(tmax_sum_61, bs = "cr") +           # Summer max temperature
    s(bio11_61, bs = "cr") +              # Winter temperature
    s(bio18_61, bs = "cr") +              # Warm quarter precipitation
    s(Elevation..m., bs = "cr") +         # Elevation
    s(Latitude, bs = "cr") +              # Latitude
    s(Longitude, bs = "cr") +             # Longitude
    ti(Latitude, Longitude, bs = "cr"),   # Spatial interaction
  data = relfit24_61,
  method = "fREML",                       # Fast REML
  discrete = FALSE
)

write_rds(gam_dist, paste0(dir, "model_outputs/models/gam_dist.rds"))

# create averaged transfer distance model for prediction ---------------------------

# Average fitness by locality to reduce noise for spatial predictions
relfit24_avg <- relfit24 %>% 
  summarize(tdiff2 = mean(tdiff2, na.rm = T), 
            mrf = mean(mrf, na.rm = T), 
            mrf_2014 = mean(mrf_2014, na.rm = T), 
            .by = c(Locality))

# Simplified model for predicting fitness under climate change
gam_rf24_avg <- bam(mrf ~ s(tdiff2, bs = "cr") + s(mrf_2014, bs = "cr"), 
                    family = "tw", 
                    data = relfit24_avg)
write_rds(gam_rf24_avg, paste0(dir, "model_outputs/models/gam_rf24_avg.rds"))

# climate raster stack to predict initial rel fitness dist from PEA ----------------

# Stack all environmental predictors needed for PEA model
clim_stack <- c(
  tmax_sum_61, 
  bio11_61, 
  bio18_61, 
  ca_dem,
  Latitude,
  Longitude
)
# Rename layers to match variable names in gam_dist model
names(clim_stack) <- c("tmax_sum_61", "bio11_61", "bio18_61", 
                       "Elevation..m.", "Latitude", "Longitude")

# predict initial rel fitness dist across range ---------------------------
## Used as the baseline for all three scenarios

# Predict current fitness across entire range using PEA model
# This represents the "pre-adaptation" fitness landscape
mrf_dist <- terra::predict(object = clim_stack, 
                           model = gam_dist, 
                           fun = predict.gam, 
                           type = "response", 
                           na.rm = TRUE)

# raster stack to predict scenario 1 --------------------------------------

clim_stack1 <- c(
  tdiff_rast_ca,  # Climate mismatch by 2100
  mrf_dist       # Initial (current) fitness distribution
  #mrf_dist_rand
)

names(clim_stack1) <- c("tdiff2", "mrf_2014")

# predict scenario 1 ------------------------------------------------------
# Fitness in 2100 if populations remain in place

mrf_pred_dist_null <- predict(clim_stack1, 
                              model = gam_rf24_avg, 
                              fun = predict.gam, 
                              type = "response") %>% 
  crop(qrange_r, mask = T)


# generate starting distributions for future climate match scenarios ------

set.seed(3333)
relfit24_vect<- relfit24 %>% ## Convert common garden data to a spatial vector of each maternal tree's fitness relative to its origin location. This allows us to select potential seed sources based on their performance in the common gardens
  right_join(select(relfit17, Accession, "mrf_2017" = mrf), by = "Accession") %>% 
  summarize(.by = Accession,  ## One value per maternal tree (average of values from each common garden)
            mrf = mean(mrf, na.rm = T),
            mrf_2017 = mean(mrf_2017, na.rm = T),
            Latitude = mean(Latitude, na.rm = T),
            Longitude = mean(Longitude, na.rm = T),
            tdiff = mean(tdiff, na.rm = T),
            tdiff2 = mean(tdiff2, na.rm = T),
            mrf_2014 = mean(mrf_2014, na.rm = T)
  ) %>% 
  slice_sample(n = 330) %>% vect( geom = c("Longitude", "Latitude"), crs = "epsg:4326")


mrf_dist2 <- interpNear(mrf_dist, ## interpolate 2024 fitness in gardens across range
                            relfit24_vect, 
                            radius = 1, 
                            field = relfit24_vect$mrf) %>% 
  mask(qrange_r) %>% 
  resample(tmax_sum_61)

mrf_dist2_17<-interpNear(mrf_dist, ## interpolate 2017 fitness in gardens across range
                         relfit24_vect, 
                         radius = 1, 
                         field = relfit24_vect$mrf_2017) %>% 
  mask(qrange_r) %>% 
  resample(tmax_sum_61)

t_dist <- c(mrf_dist2, mrf_dist2_17, tmax_sum_61) %>% ## Create distribution of fitness and temperature across range as points, to allow for radiused spatial join
  as.points() %>% 
  st_as_sf() %>% 
  mutate(id = row_number())

# For each pixel, find nearby pixels (within 25km) that are hotter than the planting site
t_dist2 <- st_join(t_dist, t_dist, 
                   join = st_is_within_distance, 
                   dist = 25000) %>% 
  filter(mean.x < mean.y) %>%  # Keep only hotter source populations
  group_by(id.x)


# Select a random source population from a warmer site within 25km for each destination
set.seed(333)
t_dist_temp <- t_dist2 %>%
  slice_sample(n = 1)

# Convert back to raster - this represents fitness of climate-informed "pre-adapted" sources
t_rast_temp <- t_dist_temp %>% 
  st_rasterize() %>%
  rast()

# Select the top performer, according to 2017 common garden data (lyr1.1.y)
t_dist_top<-t_dist2 %>%
  slice_max(n = 1, order_by = lyr1.1.y)

# Convert top performer point grid back into raster - this represents fitness of climate and fitness-informed "pre-adapted" sources
t_rast_top<-t_dist_top %>% 
  st_rasterize() %>%
  rast()

# Raster stack for predicting future climate match approach ---------------
clim_stack2 <- c(
  tdiff_rast_ca,
  resample(t_rast_temp$lyr1.y_lyr1.y, tdiff_rast_ca)  # Fitness from warmer sources
 )

names(clim_stack2) <- c("tdiff2", "mrf_2014")

# predict scenario 2 ------------------------------------------------------
# Fitness in 2100 if seeds are selected based on match to future climate

mrf_pred_dist_tdiff <- predict(clim_stack2, 
                               gam_rf24_avg, 
                               fun = predict.gam, 
                               type = "response") %>% 
  crop(qrange_r, mask = T) 



# Raster stack for predicting performance-based future climate match ------

clim_stack3 <- c(
  tdiff_rast_ca,
  resample(t_rast_top$lyr1.y_lyr1.y, tdiff_rast_ca)
)
names(clim_stack3) <- c("tdiff2", "mrf_2014")

# predict scenario 3 ------------------------------------------------------
# Fitness in 2100 if seeds are selected based on match to future climate AND performance in common gardens

mrf_pred_dist_top <- predict(clim_stack3, 
                             gam_rf24_avg, 
                             fun = predict.gam, 
                             type = "response") %>% 
  crop(qrange_r, mask = T)



# Assemble scenarios for map ----------------------------------------------

mrf_dist_crop<-mrf_dist_crop<-crop(mrf_dist, qrange_r, mask = T) ## starting distribution, cropped to species range limits
mrf_pred_plot <- c(
  ((mrf_pred_dist_null - mrf_dist_crop) / mrf_dist_crop * 100),   # Scenario 1: No action
  ((mrf_pred_dist_tdiff - mrf_dist_crop) / mrf_dist_crop * 100),  # Scenario 2: Climate match
  ((mrf_pred_dist_top - mrf_dist_crop) / mrf_dist_crop * 100)     # Scenario 3: Top performers
)

names(mrf_pred_plot) <- c("Local acorns", 
                          "Future climate match", 
                          "Performance-based future climate match")
writeRaster(mrf_pred_plot, paste0(dir, "model_outputs/predictions/mrf_pred_plot.tif"), overwrite = T)


# model cross-validation --------------------------------------------------

set.seed(999)
# phenotype env association -----------------------------------------------
cv_dat_gam_dist<-crossv_kfold(relfit24_61, k = 10) %>%
  mutate(
    train = map(train, as_tibble),
    test = map(test, as_tibble))

cv_gam_dist<-cv_dat_gam_dist %>%
  mutate(
    mod = map(train, \(df) bam(
      mrf ~ 
        s(tmax_sum_61, bs = "cr") +           # Summer max temperature
        s(bio11_61, bs = "cr") +              # Winter temperature
        s(bio18_61, bs = "cr") +              # Warm quarter precipitation
        s(Elevation..m., bs = "cr") +         # Elevation
        s(Latitude, bs = "cr") +              # Latitude
        s(Longitude, bs = "cr") +             # Longitude
        ti(Latitude, Longitude, bs = "cr"),   # Spatial interaction
      data = as_tibble(df),
      method = "fREML",                       # Fast REML
      discrete = FALSE))) %>%
  mutate(rmse = map2_dbl(mod, test, \(mod, df) rmse(model = mod, data = df)))

cv_gam_dist %>%
  select(rmse) %>%
  summarize(mean_rmse = mean(rmse), sd_rmse = sd(rmse))

# 0.158  0.0119


# rf summary --------------------------------------------------------------

relfit24_avg <- relfit24 %>%
  summarize(tdiff2 = mean(tdiff2, na.rm = T),
            mrf = mean(mrf, na.rm = T),
            mrf_2014 = mean(mrf_2014, na.rm = T),
            .by = c(Locality))

cv_dat_gam_rf24_avg<-crossv_kfold(relfit24_avg, k = 10) %>%
  mutate(
    train = map(train, as_tibble),
    test = map(test, as_tibble))

cv_gam_rf24_avg<-cv_dat_gam_rf24_avg %>%
  mutate(
    mod = map(train, \(df) bam(mrf ~ s(tdiff2, bs = "cr") + s(mrf_2014, bs = "cr"),
                               family = "tw",
                               data = as_tibble(df),
                               method = "fREML",                       # Fast REML
                               discrete = FALSE))) %>%
  mutate(rmse = map2_dbl(mod, test, \(mod, df) rmse(model = mod, data = df)))

cv_gam_rf24_avg %>%
  select(rmse) %>%
  summarize(mean_rmse = mean(rmse), sd_rmse = sd(rmse))

# 1.41  0.0531


