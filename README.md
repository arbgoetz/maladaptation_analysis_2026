Scripts and underlying data to reproduce the main analyses in "Persistent climate maladaptation decreases performance of a keystone tree species (Quercus lobata)" (Goetz et al.).

Data can be downloaded from the following FigShare repository: 10.6084/m9.figshare.31937499

Prior to running analyses, save the main source datasets ("rgrwide_input.csv", "rgryr_input.csv", and "mff24_input.csv"), along with script "00_convenience_functions.R", in the working directory. Refer to "Data column descriptions.xlsx" for descriptions of columns in the main source datasets.

Files are organized into folders by analysis step. To replicate the complete data pipeline, run all scripts in the order shown below.

Contents of folder "00_process_raw_data":

To calculate the standardized height metric estimated from allometry, run script "01a_fit_models_allometry_annotated.R", which loads the following raw data files: "census_main_core.csv", "biomass_ifg.xlsx", and "biomass_cso.xlsx"

To conduct the remaining data cleaning and calculate key variables (relative growth rate, relative fitness, and temperature transfer distances), run script "01c_process_raw_data.R", which loads the following data files: "census_main_core.csv", "prism_annual_max.csv", "prism_monthly.csv", "prism_annual_mean.csv", "maternal_climate.csv", "allometry_pred_hts.csv", and "budburst_date1.csv".

Script "01b_spatial_calculations.R" is included to demonstrate how spatial files were summarized for analysis. Raw data were collected from external sources (AdaptWest and CHELSA-TRACE21k). Specific file names in the script are unchanged from the raw downloads, and contain details on the source and contents of each raster.

Contents of folder "01_heritability_rgr_analyses":

To reproduce calculation and analysis of heritability of growth, run "03b_calculate_heritability.R" in R, which loads "rgrwide_input.csv"

To reproduce analyses of climate effects on relative growth rates, including model cross-validation, run "03_fit_models_RGR_annotated.R" in R, which loads "rgrwide_input.csv", "rgryr_input.csv", "prismtmax.csv", "prismtmax2.csv", and "allometry_pred_hts.csv" (contained in folder "00_process_raw_data").

Contents of folder "02_relative_fitness_analysis":

To reproduce analyses of climate effects on relative fitness, including model cross-validation, run "04_fit_models_MFF_annotated.R" in R, which loads "mff24_input.csv"

Contents of folder "03_spatial_prediction_figures":

To reproduce spatial predictions of future relative fitness, including model cross-validation, run "05_predict_AGF_outcomes_annotated_revised.R" in R, which loads all raster (.tif) files, along with "qlobata_range.shp", "gardencoords.csv", "mff24_input.csv", "garden_aves_rast.csv", "site_clim_avgs.csv", and "maternal_climate.csv"

To reproduce figures, run "06_create_figures.R" and "07_create_maps.R" scripts in R, which load all data in "predictions" subfolder, as well as "site_clim_avgs.csv", "garden_aves_rast.csv", "rgryr_input.csv", "rgrwide_input.csv", "dat_fig3c.csv", "maternal_climate.csv", "CA_state.shp", "qlobata_range.shp", "qlobata_range2.shp", and "gardencoords.csv"
