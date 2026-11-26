# Set working directory
dir<-"C:/Users/alexg/Downloads/31937499/"

# Source custom helper functions
source(paste0(dir, "00_convenience_functions.R"))

# Install specific version of tidygam package (version 0.2.2 required for compatibility)
#devtools::install_version("tidygam", version = "0.2.2")

# load libraries ----------------------------------------------------------

# Load required packages for data manipulation, modeling, and analysis
library(mgcv)        # Generalized Additive Models (GAMs)
library(gratia)      # Helper functions for GAM visualization
library(tidygam)     # Tidy interface for GAM predictions (version 0.2.0 required)
library(tidyverse)   # Suite of data manipulation packages (dplyr, ggplot2, etc.)
library(stringr)
library(conflicted)  # Handle function name conflicts between packages
library(GauPro)      # Gaussian Process regression for smoothing predictions

# Load required packages for cross-validation
set.seed(999)
library(modelr)
library(mgcViz)
library(DHARMa)

# Resolve function conflicts by preferring dplyr versions
conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)

# load in data ------------------------------------------------------------

relfit24 <- read.csv(paste0(dir, "mff24_input.csv")) %>%
  mutate(Accession = as.factor(Accession),
         Locality = as.factor(Locality))

tdiffrange <- seq(-20, 10, length.out = 100)

# fit GAMs ----------------------------------------------------------------

## 2024 cumulative fitness analysis

# Fit GAM: fitness ~ temperature difference + initial fitness + random effects
# Uses Tweedie distribution for positive continuous data with zeros
gam_rf24 <- bam(mrf ~ s(tdiff, bs = "cr", k = 5) +  # Smooth effect of temp difference
                  s(mrf_2014, bs = "cr") +           # Smooth effect of initial fitness
                  s(Accession, bs = "re") +          # Accession random effect
                  s(Locality, bs = "re") +           # Locality random effect
                  Site,                               # Fixed effect of site
                family = "tw",                        # Tweedie family
                data = filter(relfit24, mrf > 0))    # Only positive fitness values

# Generate predictions for temperature difference effect
p_gamrf24 <- predict_gam(gam_rf24, 
                         series = "tdiff", 
                         exclude_terms = c("s(Accession)", "s(Locality)", "Site"), 
                         length_out = 50,
                         tran_fun = exp)  # Back-transform from log scale

# Smooth predictions using Gaussian Process
gpc24 <- predict(gpkm(X = p_gamrf24$tdiff, Z = p_gamrf24$mrf), tdiffrange, se.fit = TRUE)
gpc24$tdiff <- tdiffrange
gpc24$mrf <- gpc24$mean
gpc24$lower <- gpc24$mrf - 1.96 * gpc24$se  # 95% CI
gpc24$upper <- gpc24$mrf + 1.96 * gpc24$se

# Save model and predictions
write_rds(gam_rf24, paste0(dir, "model_outputs/models/cumulative_2014_2024/mff24_gam.rda"))
# gam_rf24 <- read_rds(paste0(dir, "model_outputs/models/cumulative_2014_2024/mff24_gam.rda"))
write.csv(p_gamrf24, paste0(dir, "model_outputs/predictions/cumulative_2014_2024/mff24_pred.csv"))
write.csv(gpc24, paste0(dir, "model_outputs/predictions/cumulative_2014_2024/mff24_pred_gau.csv"))

## Alternative model using average temperature metric
gam_rf24_2 <- bam(mrf ~ s(tdiff2, bs = "cr", k = 12) +
                    s(mrf_2014, bs = "cr") +
                    s(Accession, bs = "re") + 
                    s(Locality, bs = "re") +
                    Site,
                  family = "tw",
                  data = filter(relfit24, mrf > 0))

# Generate predictions
p_gamrf24_2 <- predict_gam(gam_rf24_2, 
                           series = "tdiff2", 
                           exclude_terms = c("s(Accession)", "s(Locality)", "Site"), 
                           length_out = 50,
                           tran_fun = exp)

# Save model and predictions
write_rds(gam_rf24_2, paste0(dir, "model_outputs/models/cumulative_2014_2024/mff24_tave_gam.rda"))
write.csv(p_gamrf24_2, paste0(dir, "model_outputs/predictions/cumulative_2014_2024/mff24_tave_pred.csv"))

# Smooth predictions with Gaussian Process
gpc24_2 <- predict(gpkm(X = p_gamrf24_2$tdiff2, Z = p_gamrf24_2$mrf), tdiffrange, se.fit = TRUE)
gpc24_2$tdiff2 <- tdiffrange
gpc24_2$mrf <- gpc24_2$mean
gpc24_2$lower <- gpc24_2$mrf - 1.96 * gpc24_2$se
gpc24_2$upper <- gpc24_2$mrf + 1.96 * gpc24_2$se

#write.csv(gpc24_2, paste0(dir, "model_outputs/predictions/cumulative_2014_2024/mff24_tave_pred_gau.csv"))
max_mrf <- max(gpc24_2$mrf)  # Maximum predicted fitness: 0.9640845

## Model using winter temperature difference
gam_rf24_win <- bam(mrf ~ s(tdiff_win, bs = "cr", k = 12) +
                      s(mrf_2014, bs = "cr") +
                      s(Accession, bs = "re") + 
                      s(Locality, bs = "re") +
                      Site,
                    family = "tw",
                    data = filter(relfit24, mrf > 0))

# Generate predictions for winter temperature effect
p_gamrf24_win <- predict_gam(gam_rf24_win, 
                             series = "tdiff_win", 
                             exclude_terms = c("s(Accession)", "s(Locality)", "Site"), 
                             length_out = 50,
                             tran_fun = exp)

# Save model and predictions
write_rds(gam_rf24_win, paste0(dir, "model_outputs/models/cumulative_2014_2024/mff24_win_gam.rda"))
write.csv(p_gamrf24_win, paste0(dir, "model_outputs/predictions/cumulative_2014_2024/mff24_win_pred.csv"))

# geographic analysis -----------------------------------------------------

# Linear model: fitness ~ climate PCs + geographic location (IFG site)
lm_clim_ifg <- lm(mrf ~ scale(PC1) + scale(PC2) + scale(Latitude) + scale(Longitude), 
                  data = filter(relfit24, Site == "IFG"))
write_rds(lm_clim_ifg, paste0(dir, "model_outputs/models/mffyr_geog_lm_ifg.rda"))

# Linear model: fitness ~ climate PCs + geographic location (Chico site)
lm_clim_cso <- lm(mrf ~ scale(PC1) + scale(PC2) + scale(Latitude) + scale(Longitude), 
                  data = filter(relfit24, Site == "Chico"))
write_rds(lm_clim_cso, paste0(dir, "model_outputs/models/mffyr_geog_lm_cso.rda"))


# model cross-validation --------------------------------------------------

# mff 2024 ----------------------------------------------------------------


cv_dat_mff24<-crossv_kfold(filter(relfit24, mrf > 0), k = 10) %>%
  mutate(
    train = map(train, as_tibble),
    test = map(test, as_tibble))
cv_mff24<-cv_dat_mff24 %>%
  mutate(
    mod = map(train, \(df) bam(mrf ~ s(tdiff, bs = "cr", k = 5) +
                                 s(mrf_2014, bs = "cr") +
                                 s(Accession, bs = "re") + 
                                 s(Locality, bs = "re") +
                                 Site,
                               family = "tw", data = as_tibble(df))
    )) %>%
  mutate(rmse = map2_dbl(mod, test, \(mod, df) rmse(model = mod, data = df)))


cv_mff24 %>%
  select(rmse) %>%
  summarize(mean_rmse = mean(rmse), sd_rmse = sd(rmse))

#  1.40  0.0174

# mff 2024_tdiffave -------------------------------------------------------

cv_dat_mff24<-crossv_kfold(filter(relfit24, mrf > 0), k = 10) %>%
  mutate(
    train = map(train, as_tibble),
    test = map(test, as_tibble))
cv_mff24_2<-cv_dat_mff24 %>%
  mutate(
    mod = map(train, \(df) bam(mrf ~ s(tdiff2, bs = "cr", k = 12) +
                                 s(mrf_2014, bs = "cr") +
                                 s(Accession, bs = "re") + 
                                 s(Locality, bs = "re") +
                                 Site,
                               family = "tw", data = as_tibble(df))
    )) %>%
  mutate(rmse = map2_dbl(mod, test, \(mod, df) rmse(model = mod, data = df)))


cv_mff24_2 %>%
  select(rmse) %>%
  summarize(mean_rmse = mean(rmse), sd_rmse = sd(rmse))

#  1.40  0.0163

# mff 2024 win ------------------------------------------------------------

cv_mff24_win<-cv_dat_mff24 %>%
  mutate(
    mod = map(train, \(df) bam(mrf ~ s(tdiff_win, bs = "cr", k = 12) +
                                 s(mrf_2014, bs = "cr") +
                                 s(Accession, bs = "re") + 
                                 s(Locality, bs = "re") +
                                 Site,
                               family = "tw", data = as_tibble(df))
    )) %>%
  mutate(rmse = map2_dbl(mod, test, \(mod, df) rmse(model = mod, data = df)))


cv_mff24_win %>%
  select(rmse) %>%
  summarize(mean_rmse = mean(rmse), sd_rmse = sd(rmse))

#  1.40  0.0182
