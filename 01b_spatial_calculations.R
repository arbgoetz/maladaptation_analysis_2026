# Set working directory
dir<-"C:/Users/alexg/Downloads/31937499/"

# load libraries ----------------------------------------------------------
library(terra)
library(sf)
library(tidyterra)
library(scales)
library(gstat)
library(concaveman)
library(tidyverse)
library(ragg)
library(ggspatial)

conflicted::conflicts_prefer(terra::mask)
conflicted::conflicts_prefer(dplyr::summarize)

## all climate rasters were obtained from: https://adaptwest.databasin.org/pages/adaptwest-climatena/ except "bio11_lgm", which was obtained from: https://www.chelsa-climate.org/datasets/chelsa-trace21k-centennial-bioclim
## Citation: AdaptWest Project. 2022. Gridded current and projected climate data for North America at 1km resolution, generated using the ClimateNA v7.30 software (T. Wang et al., 2022). Available at adaptwest.databasin.org.
## refer to full raster filenames (below) for time period, GCM, SSP scenario, and specific variable - all filenames are unchanged from the defaults upon downloading

## california boundary from 
calif<-vect(paste0(dir, "03_spatial_prediction_figures/ca_state/CA_state.shp")) %>%
  project("epsg:4326")


## maternal origin climate data California Open Data (https://data.ca.gov/@california-department-of-technology/ca-geographic-boundaries)
dat_clim_pc<-read.csv(paste0(dir, "03_spatial_prediction_figures/maternal_climate.csv")) %>%
  mutate(Accession = as.factor(Accession))


sf_use_s2(TRUE)
mat_coords<-st_as_sf(dplyr::select(dat_clim_pc, Longitude, Latitude, Accession), coords = c("Longitude", "Latitude"), crs = 4326)
gardens<-read.csv(paste0(dir, "03_spatial_prediction_figures/gardencoords.csv"))

# generate hypothesized sp range ------------------------------------------
qrange<- mat_coords %>%
  concaveman() %>%
  st_transform("+proj=eqc") %>%
  st_buffer(25000) %>%
  st_intersection(st_transform(st_as_sf(calif), "+proj=eqc")) %>%
  st_transform("epsg:4326")

st_write(qrange, paste0(dir, "03_spatial_prediction_figures/qlobata_range.shp"))
qrange<-st_read(paste0(dir, "03_spatial_prediction_figures/qlobata_range.shp"))

qrange_r<-vect("~/git/demographyMS/intermediate_data/qlobata_range.shp") %>% ## spatvector version of qrange for use with terra spatrasters
  project("epsg:4326")

##buffered around maternal tree sites

qrange2<- mat_coords %>%
  st_union() %>%
  st_transform("+proj=eqc") %>%
  st_buffer(25000) %>%
  st_as_sf() %>%
  st_intersection(st_transform(st_as_sf(calif), "+proj=eqc")) %>%
  st_transform("epsg:4326")
st_write(qrange2, paste0(dir, "03_spatial_prediction_figures/qlobata_range2.shp"))

# rasters

## digital elevation model from USGS 3D Elevation Program (https://portal.opentopography.org/datasetMetadata?otCollectionID=OT.012021.4269.3)
ca_dem<-rast("C:/Users/alexg/Box/UCLA postdoc - Copy/geospatial/elevation.tif") %>%
  project("epsg:4326") %>%
  crop(calif, mask = T)
writeRaster(ca_dem, paste0(dir, "03_spatial_prediction_figures/ca_dem.tif"))

# read and project bioclim rasters ----------------------------------------

## bioclim_11 (mean winter temp) 1961-1990
bio11_61<-rast("C:/Users/alexg/Box/UCLA postdoc - Copy/geospatial/Normal_1961_1990_bioclim/Normal_1961_1990/Normal_1961_1990_bioclim/Normal_1961_1990_Tave_wt.tif") %>%
  project("epsg:4326") %>%
  crop(calif, mask= T)
writeRaster(bio11_61, paste0(dir, "03_spatial_prediction_figures/bio11_61.tif"))

## bioclim_11 (mean winter temp) 1991-2020
bio11_current <- rast("C:/Users/alexg/Box/UCLA postdoc - Copy/geospatial/Normal_1991_2020_bioclim/Normal_1991_2020_bioclim/Normal_1991_2020_Tave_wt.tif") %>%
  project("epsg:4326") %>%
  crop(calif, mask = T)
writeRaster(bio11_current, paste0(dir, "03_spatial_prediction_figures/bio11_current.tif"))


## bioclim_11 (mean winter temp) 2070-2100, RCP 4.5
bio11_45 <- rast("C:/Users/alexg/Box/UCLA postdoc - Copy/geospatial/ensemble_8GCMs_ssp245_2071_2100/ensemble_8GCMs_ssp245_2071_2100_bioclim/ensemble_8GCMs_ssp245_2071_2100_Tave_wt.tif") %>%
  project("epsg:4326") %>%
  crop(calif, mask= T)
writeRaster(bio11_45, paste0(dir, "03_spatial_prediction_figures/bio11_45.tif"))

##  summer precip 1961-1990
bio18_61<-rast("C:/Users/alexg/Box/UCLA postdoc - Copy/geospatial/Normal_1961_1990_bioclim/Normal_1961_1990/Normal_1961_1990_bioclim/Normal_1961_1990_PPT_sm.tif") %>%
  project("epsg:4326") %>%
  crop(calif, mask= T)
writeRaster(bio18_61, paste0(dir, "03_spatial_prediction_figures/bio18_61.tif"))

bio18_current<-rast("C:/Users/alexg/Box/UCLA postdoc - Copy/geospatial/Normal_1991_2020_bioclim/Normal_1991_2020_bioclim/Normal_1991_2020_MSP.tif") %>%
  project("epsg:4326") %>%
  crop(calif, mask = T)
writeRaster(bio18_current, paste0(dir, "03_spatial_prediction_figures/bio18_current.tif", overwrite = T))

# calculate tmax summary rasters -----------------------------------

## tmax 1961-1990

tmax6_61<-terra::rast("C:/Users/alexg/Box/UCLA postdoc - Copy/geospatial/Normal_1961_1990_monthly/Normal_1961_1990/Normal_1961_1990_monthly/Normal_1961_1990_Tmax06.tif") %>%
  project("epsg:4326") %>%
  crop(calif, mask= T)
tmax7_61<-terra::rast("C:/Users/alexg/Box/UCLA postdoc - Copy/geospatial/Normal_1961_1990_monthly/Normal_1961_1990/Normal_1961_1990_monthly/Normal_1961_1990_Tmax07.tif") %>%
  project("epsg:4326") %>%
  crop(calif, mask= T)
tmax8_61<-terra::rast("C:/Users/alexg/Box/UCLA postdoc - Copy/geospatial/Normal_1961_1990_monthly/Normal_1961_1990/Normal_1961_1990_monthly/Normal_1961_1990_Tmax08.tif") %>%
  project("epsg:4326") %>%
  crop(calif, mask= T)

tmax_sum_61<-mean(c(tmax6_61, tmax7_61, tmax8_61))
writeRaster(tmax_sum_61, paste0(dir, "03_spatial_prediction_figures/tmax_sum_61.tif", overwrite = TRUE))


## tmax 1991-2020 (for Fig 1)
tmax6<-rast("C:/Users/alexg/Box/UCLA postdoc - Copy/geospatial/Normal_1991_2020_monthly/Normal_1991_2020/Normal_1991_2020_monthly/Normal_1991_2020_Tmax06.tif") %>%
  project("epsg:4326") %>%
  crop(calif, mask= T)
tmax7<-rast("C:/Users/alexg/Box/UCLA postdoc - Copy/geospatial/Normal_1991_2020_monthly/Normal_1991_2020/Normal_1991_2020_monthly/Normal_1991_2020_Tmax07.tif") %>%
  project("epsg:4326") %>%
  crop(calif, mask= T)
tmax8<-rast("C:/Users/alexg/Box/UCLA postdoc - Copy/geospatial/Normal_1991_2020_monthly/Normal_1991_2020/Normal_1991_2020_monthly/Normal_1991_2020_Tmax08.tif") %>%
  project("epsg:4326") %>%
  crop(calif, mask= T)

tmax_sum<-mean(c(tmax6, tmax7, tmax8))
writeRaster(tmax_sum, paste0(dir, "03_spatial_prediction_figures/tmax_sum.tif", overwrite = T))
tmax_sum<-rast(paste0(dir, "03_spatial_prediction_figures/tmax_sum.tif"))

## tmax 2070-2100 (RCP 4.5)
tmaxpred6_45<-rast("C:/Users/alexg/Box/UCLA postdoc - Copy/geospatial/ensemble_8GCMs_ssp245_2071_2100_monthly/ensemble_8GCMs_ssp245_2071_2100/ensemble_8GCMs_ssp245_2071_2100_monthly/ensemble_8GCMs_ssp245_2071_2100_Tmax06.tif") %>%
  project("epsg:4326") %>%
  crop(calif, mask = T)
tmaxpred7_45<-rast("C:/Users/alexg/Box/UCLA postdoc - Copy/geospatial/ensemble_8GCMs_ssp245_2071_2100_monthly/ensemble_8GCMs_ssp245_2071_2100/ensemble_8GCMs_ssp245_2071_2100_monthly/ensemble_8GCMs_ssp245_2071_2100_Tmax07.tif") %>%
  project("epsg:4326") %>%
  crop(calif, mask = T)
tmaxpred8_45<-rast("C:/Users/alexg/Box/UCLA postdoc - Copy/geospatial/ensemble_8GCMs_ssp245_2071_2100_monthly/ensemble_8GCMs_ssp245_2071_2100/ensemble_8GCMs_ssp245_2071_2100_monthly/ensemble_8GCMs_ssp245_2071_2100_Tmax08.tif") %>%
  project("epsg:4326") %>%
  crop(calif, mask = T)

stack<-c(tmaxpred6_45, tmaxpred7_45, tmaxpred8_45)

tmaxpred_sum_45<-mean(stack)
writeRaster(tmaxpred_sum_45, paste0(dir, "03_spatial_prediction_figures/tmaxpred_sum_45.tif", overwrite = TRUE))
tmaxpred_sum_45<-rast(paste0(dir, "03_spatial_prediction_figures/tmaxpred_sum_45.tif"))



## calculate lat and long from tmax raster
Longitude <- init(tmax_sum_61, 'x')
writeRaster(Longitude, paste0(dir, "03_spatial_prediction_figures/Longitude.tif", overwrite = T))
Latitude <- init(tmax_sum_61, 'y')
writeRaster(Latitude, paste0(dir, "03_spatial_prediction_figures/Latitude.tif", overwrite = T))

## calculate tdiff between 1961-1990 and 2070-2100
tdiff_rast_ca<-(tmaxpred_sum_45 - tmax_sum_61)
writeRaster(tdiff_rast_ca, paste0(dir, "03_spatial_prediction_figures/tdiff_rast_ca.tif", overwrite = T))
tdiff_rast_ca<-rast(paste0(dir, "03_spatial_prediction_figures/tdiff_rast_ca.tif"))

## calculate past, present, future temperature summaries for figures


### RCP 4.5 - winter average temp

bio11_45_gardenave <- bio11_45 %>%
  terra::extract(st_as_sf(gardens, coords = c("Longitude", "Latitude"), crs = 4326)) 
bio11_45_gardenave<- mean(c(9.6, 10.7)) 


### RCP 4.5 - summer max temp

tmaxpred_sum_gardenave<-tmaxpred_sum_45 %>%
  terra::extract(st_as_sf(gardens, coords = c("Longitude", "Latitude"), crs = 4326)) ## find predicted future temperatures at the common gardens in 2070-2100
tmaxpred_sum_gardenave<-mean(c(34.83333,37.4)) ## average of the two gardens - to annotate plot


###  Last Glacial Maximum (CHELSA) - winter average temp
bio11_lgm <- rast("C:/Users/alexg/Box/UCLA postdoc - Copy/geospatial/chelsa/CHELSA_PMIP_CCSM4_BIO_11.tif")
bio11_lgm_gardenave<-bio11_lgm %>%
  terra::extract(st_as_sf(gardens, coords = c("Longitude", "Latitude"), crs = 4326)) %>%
  mutate(CHELSA_PMIP_CCSM4_BIO_11 = (CHELSA_PMIP_CCSM4_BIO_11/10) - 273.15) %>%
  summarize(mean = mean(CHELSA_PMIP_CCSM4_BIO_11)) %>% pull() ##-0.7


tibble(
  "name" = c("bio11_45_gardenave", "tmaxpred_sum_gardenave", "bio11_lgm_gardenave"),
  "value" = c(bio11_45_gardenave, tmaxpred_sum_gardenave, bio11_lgm_gardenave)
) %>%
  write.csv(paste0(dir, "03_spatial_prediction_figures/garden_aves_rast.csv"))

## interpolate performance values for Fig S6

ifg_lm<-gstat(formula = mrf ~ 1, locations = ~Longitude+Latitude, data = filter(relfit_all, Site == "IFG", Year == 2024, !is.na(mrf)))
cso_lm<-gstat(formula = mrf ~ 1, locations = ~Longitude+Latitude, data = filter(relfit_all, Site == "Chico", Year == 2024, !is.na(mrf)))


interp_ifg<-terra::interpolate(tmax_sum_current, model=ifg_lm, xyNames = c("Longitude", "Latitude")) %>%
  project("epsg:4326") %>%
  crop(calif) %>%
  mask(calif)

interp_cso<-terra::interpolate(tmax_sum_current, model=cso_lm, xyNames = c("Longitude", "Latitude")) %>%
  project("epsg:4326") %>%
  crop(calif) %>%
  mask(calif)


## variogram - check interpolation
plot(variogram(mrf ~ 1, locations = ~Longitude+Latitude, data = select(filter(relfit_all, Site == "IFG", Year == 2024, !is.na(mrf)), mrf, Latitude, Longitude, PC1, PC2) ))
plot(variogram(mrf ~ 1, locations = ~Longitude+Latitude, data = select(filter(relfit_all, Site == "Chico", Year == 2024, !is.na(mrf)), mrf, Latitude, Longitude, PC1, PC2) ))

writeRaster(interp_ifg, paste0(dir, "03_spatial_prediction_figures/predictions/interp_ifg.tif"))
writeRaster(interp_cso, paste0(dir, "03_spatial_prediction_figures/predictions/interp_cso.tif"))
