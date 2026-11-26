# load libraries ----------------------------------------------------------

library(tidyverse)
library(readxl)

conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)


# Set working directory
dir<-"C:/Users/alexg/Downloads/31937499/"

# Source custom helper functions
source(paste0(dir, "00_convenience_functions.R"))


# Process datasets for RGR analyses ---------------------------------------


# Raw census data with plant measurements
long_core<-read.csv(paste0(dir, "00_process_raw_data/census_main_core.csv"))
long_core <- long_core %>%
  rename("ht" = Height_cm) %>%
  mutate(Accession = as.numeric(Accession)) 

# Annual maximum temperature data from PRISM. Refer to 01b
prismtmax <- read.csv(paste0(dir, "00_process_raw_data/prism_annual_max.csv")) %>%
  group_by(Site) %>%
  mutate(Year = as.integer(Year),
         tmaxprev = lag(tmax_sum),        # Previous year's temperature
         pptprev = lag(ppt),              # Previous year's precipitation
         tmaxrel = tmax_sum/mean(tmax_sum), # Relative temperature
         tmaxrelprev = lag(tmaxrel))      # Previous year's relative temp

# Monthly climate data from PRISM
prismmonthly <- read.csv(paste0(dir, "00_process_raw_data/prism_monthly.csv")) %>%
  separate(Date, into = c("Year", "Month"), sep = "-") %>%
  mutate(Year = as.integer(Year),
         Month = as.integer(Month)) %>%
  rename("ppt_mm" = ppt..mm.,
         "tmax" = tmax..degrees.C.,
         "tmin" = tmin..degrees.C.,
         "tave" = tmean..degrees.C.)

# Calculate average temperature for the 3 coldest months per year (winter temperature)
prismtmax <- prismmonthly %>%
  mutate(
    # Assign winter months to the following year (Dec 2014 → winter 2015)
    Year_winter = case_when(
      Month == 12 ~ lead(Year),
      Month == 11 ~ lead(Year, n = 2L),
      .default = Year
    )
  ) %>% 
  group_by(Site, Year_winter) %>%
  slice_min(order_by = tave, n = 3) %>%  # Select 3 coldest months
  summarize(tave_win = mean(tave)) %>%   # Average of coldest 3 months
  ungroup() %>%
  right_join(prismtmax, by = c("Site", "Year_winter" = "Year")) %>%
  rename("Year" = Year_winter)

write.csv(prismtmax, paste0(dir, "01_heritability_rgr_analyses/prismtmax.csv"))

# Annual mean temperature data (alternative temperature metric)
prismtmax2 <- read.csv(paste0(dir, "00_process_raw_data/prism_annual_mean.csv")) %>%
  group_by(Site) %>%
  mutate(Year = as.factor(Year)) %>%
  rename("tmax_sum_2" = tmax_sum)

write.csv(prismtmax2, paste0(dir, "01_heritability_rgr_analyses/prismtmax2.csv"))

# Calculate average climate at garden sites (2014-2024)
prismtmax_numyr_2 <- prismtmax2 %>% mutate(Year = as.numeric(as.character(Year)))
ifg_tmax2 <- mean(filter(prismtmax_numyr_2, Site == "IFG", Year >= 2014, Year <= 2024)$tmax_sum_2)
cso_tmax2 <- mean(filter(prismtmax_numyr_2, Site == "Chico", Year >= 2014, Year <= 2024)$tmax_sum_2)

# Maternal climate data (climate at origin sites), from Basin Characterization Model (Flint et al. 2013). Refer to Browne et al. (2019) for details.
dat_clim_pc <- read.csv(paste0(dir, "03_spatial_prediction_figures/maternal_climate.csv")) %>%
  mutate(Accession = as.factor(Accession))

# Predicted heights from allometry models
## Calculated in 01a_fit_models_allometry_annotated.R
allom_dat_pred3 <- read.csv(paste0(dir, "00_process_raw_data/allometry_pred_hts.csv")) %>%
mutate(Year = as.integer(Year),
         Accession = as.factor(Accession))

# Budburst phenology data (date of first budburst)
dat_d1_bb <- read.csv(paste0(dir, "00_process_raw_data/budburst_date1.csv")) %>%
  separate(Budburst_d1_ID, c("Accession", "progeny", "Year"), sep = "-") %>%
  mutate(Year = as.integer(Year)) %>%
  filter(!is.na(Year)) %>%
  unite("Accession_progeny", Accession, progeny, sep = "-", remove = FALSE)

# Sequence for temperature difference predictions
tdiffrange <- seq(-20, 10, length.out = 100)

# RGR data transformations --------------------------------------------------------

# Calculate average maximum temperatures at garden sites
prismtmax_numyr <- prismtmax %>% mutate(Year = as.numeric(as.character(Year)))
ifg_tmax <- mean(filter(prismtmax_numyr, Site == "IFG", Year >= 2014, Year <= 2024)$tmax_sum)
cso_tmax <- mean(filter(prismtmax_numyr, Site == "Chico", Year >= 2014, Year <= 2024)$tmax_sum)

# Calculate average winter temperatures at garden sites
ifg_tavewin <- mean(filter(prismtmax_numyr, Site == "IFG", Year >= 2014, Year <= 2024)$tave_win)
cso_tavewin <- mean(filter(prismtmax_numyr, Site == "Chico", Year >= 2014, Year <= 2024)$tave_win)

# Prepare initial dataset for RGR calculations
dat_rgr_init <- long_core %>%
  filter(!is.na(Accession_progeny), !is.na(Accession), Year != 2013) %>%
  arrange(Accession_progeny) %>%
  group_by(Accession_progeny) %>%
  fill(Site, Section, Block, Row, Line, Border, .direction = "up") %>%  # Fill missing metadata
  mutate(Year = as.integer(Year)) %>%
  # Join predicted heights
  left_join(select(allom_dat_pred3, Accession_progeny, Year, ht_complete), 
            by = c("Accession_progeny", "Year")) %>%
  mutate(ht_complete = case_when(
    Year == 2013 | Year == 2014 ~ ht,  # Use measured heights for early years
    .default = ht_complete
  )) %>%
  # Filter out border plants, gopher damage, and problematic individuals
  filter(Border == 0, Gopher_damage != 1, Accession_progeny != "245-15") %>%
  unite("SiteBlock", c(Site, Block), sep = "", remove = FALSE)

# Convert to factors
dat_rgr_init$SiteBlock <- as.factor(dat_rgr_init$SiteBlock)
dat_rgr_init$Locality <- as.factor(dat_rgr_init$Locality)
dat_rgr_init$Accession <- as.factor(dat_rgr_init$Accession)

# Calculate cumulative RGR over different time periods (wide format)
dat_rgr_htc <- dat_rgr_init %>%
  pivot_wider(id_cols = c(Site, Block, Border, Accession_progeny, Accession, Locality), 
              names_from = Year, 
              values_from = c(ht_complete, Alive_new)) %>% 
  mutate(
    # Relative Growth Rate = (ln(ht2) - ln(ht1)) / time_interval
    RGRh1424 = (log(ht_complete_2024) - log(ht_complete_2014)) / 10,  # 2014-2024
    RGRh1417 = (log(ht_complete_2017) - log(ht_complete_2014)) / 3,   # 2014-2017
    RGRh1720 = (log(ht_complete_2020) - log(ht_complete_2017)) / 3,   # 2017-2020
    RGRh2023 = (log(ht_complete_2023) - log(ht_complete_2020)) / 3,   # 2020-2023
    Locality = as.factor(Locality),
    Block = as.factor(Block),
    Accession = as.factor(Accession)
  ) %>%
  left_join(dat_clim_pc, by = c("Accession", "Locality")) %>%
  mutate(
    # Calculate climate mismatch: garden temperature - maternal temperature
    tdiff = case_when(Site == "IFG" ~ ifg_tmax - tmax_sum,
                      Site == "Chico" ~ cso_tmax - tmax_sum),
    tdiff2 = case_when(Site == "IFG" ~ ifg_tmax2 - tmax_sum,
                       Site == "Chico" ~ cso_tmax2 - tmax_sum),
    tdiff_win = case_when(Site == "IFG" ~ ifg_tavewin - bioclim_11,
                          Site == "Chico" ~ cso_tavewin - bioclim_11)
  ) %>%
  unite("SiteBlock", c(Site, Block), sep = "", remove = FALSE)

allom_dat_pred3 <- allom_dat_pred3 %>% mutate(Year = as.integer(as.character(Year)))

# Calculate year-to-year RGR (long format) for annual analysis
dat_rgr7 <- long_core %>% 
  filter(Year != 2025) %>%
  arrange(Accession_progeny, Site) %>% 
  group_by(Accession_progeny, Site) %>% 
  filter(!any(Alive_new == 0)) %>%  # Keep only plants that survived entire study
  mutate(Year = as.integer(as.character(Year))) %>% 
  filter(!is.na(Accession_progeny), !is.na(Accession)) %>% 
  left_join(select(allom_dat_pred3, Accession_progeny, Year, ht_complete), 
            by = c("Accession_progeny", "Year")) %>% 
  fill(Site, Section, Block, Row, Line, Border, .direction = "up") %>% 
  arrange(Site, Accession_progeny, Year) %>% 
  group_by(Accession_progeny) %>% 
  mutate(
    # Calculate annual RGR
    rgrh = (log(ht_complete) - lag(log(ht_complete))) / (Year - lag(Year)),
    ht1 = lag(ht_complete),  # Previous year's height
    # Rate of change in RGR (acceleration)
    rgrria = (rgrh - lag(rgrh)) / lag(rgrh),
    Accession = as.factor(Accession),
    Accession_progeny = as.factor(Accession_progeny),
    Block = as.factor(Block)
  ) %>% 
  filter(!any(Gopher_damage == 1)) %>%  # Exclude plants with any gopher damage
  filter(Border != 1, Accession_progeny != "245-15") %>% 
  # Keep only positive growth (except early years where establishment is expected)
  filter((Year == 2013 | Year == 2014 | rgrh > 0), 
         (ht_complete > ht1 | Year == 2013 | Year == 2014)) %>% 
  left_join(dat_clim_pc, by = c("Accession", "Locality")) %>% 
  left_join(select(dat_d1_bb, Accession_progeny, Year, date1), 
            by = c("Year", "Accession_progeny")) %>%
  left_join(prismtmax_numyr, by = c("Year", "Site")) %>%
  mutate(
    # Yearly climate mismatch: garden temp that year - maternal temp
    tdiffyr = tmax_sum.y - tmax_sum.x,
    tdiffyr_win = tave_win - bioclim_11,
    rgrhst = rgrh / ht1  # RGR standardized by size
  ) %>%
  left_join(prismtmax_numyr_2, by = c("Year", "Site")) %>%
  mutate(tdiffyr2 = tmax_sum_2 - tmax_sum.x) %>% 
  # Join initial height from 2014
  left_join(select(filter(long_core, Year == 2014), Accession_progeny, ht), 
            by = "Accession_progeny") %>%
  rename("ht" = ht.x,
         "ht_2014" = ht.y) %>%
  ungroup() %>%
  filter(Year != 2013) %>% 
  unite("SiteBlock", c(Site, Block), sep = "", remove = FALSE)


# Save processed datasets
write.csv(dat_rgr7, paste0(dir, "rgryr_input.csv"))
write.csv(dat_rgr_htc, paste0(dir, "rgrwide_input.csv"))



# Process data for relative fitness analyses ------------------------------


# Calculate average climate conditions at each garden site (2014-2024)
prismtmax_numyr <- prismtmax %>% mutate(Year = as.numeric(as.character(Year)))
prismtmax_numyr_2 <- prismtmax2 %>% mutate(Year = as.numeric(as.character(Year)))

# IFG site climate averages
ifg_tmax <- mean(filter(prismtmax_numyr, Site == "IFG", Year >= 2014, Year <= 2024)$tmax_sum)
ifg_tmax_ave <- mean(filter(prismtmax_numyr_2, Site == "IFG", Year >= 2014, Year <= 2024)$tmax_sum_2)
ifg_ppt <- mean(filter(prismtmax_numyr, Site == "IFG", Year >= 2014, Year <= 2024)$ppt)
ifg_tavewin <- mean(filter(prismtmax_numyr, Site == "IFG", Year >= 2014, Year <= 2024)$tave_win)

# Chico site climate averages
cso_tmax <- mean(filter(prismtmax_numyr, Site == "Chico", Year >= 2014, Year <= 2024)$tmax_sum)
cso_tmax_ave <- mean(filter(prismtmax_numyr_2, Site == "Chico", Year >= 2014, Year <= 2024)$tmax_sum_2)
cso_ppt <- mean(filter(prismtmax_numyr, Site == "Chico", Year >= 2014, Year <= 2024)$ppt)
cso_tavewin <- mean(filter(prismtmax_numyr, Site == "Chico", Year >= 2014, Year <= 2024)$tave_win)

tibble(
  "name" = c("ifg_tmax", "ifg_tmax_ave", "ifg_ppt", "ifg_tavewin", "cso_tmax", "cso_tmax_ave", "cso_ppt", "cso_tavewin"),
  "var" = c(ifg_tmax, ifg_tmax_ave, ifg_ppt, ifg_tavewin, cso_tmax, cso_tmax_ave, cso_ppt, cso_tavewin)
) %>%
  write.csv(paste0(dir, "03_spatial_prediction_figures/site_clim_avgs.csv"))

allom_dat_pred3<-allom_dat_pred3 %>% mutate(Year = as.factor(Year))

# Load and clean census data for core experiment
long_core_relfit <- read.csv(paste0(dir, "00_process_raw_data/census_main_core.csv")) %>%
  rename("ht" = Height_cm) %>%
  mutate(Year = as.factor(Year), Accession = as.factor(Accession)) %>%
  arrange(Accession_progeny) %>%
  group_by(Accession_progeny) %>%
  fill(Site, Section, Block, Row, Line, Border, .direction = "up") %>%  # Fill missing metadata
  ungroup() %>%
  # Filter out border plants, missing IDs, and specific problematic individuals
  filter(!is.na(Accession_progeny), !is.na(Accession), !is.na(Block), Border == 0,  
         Accession_progeny != "245-15", !str_detect(Accession_progeny, "[:alpha:]")) %>%
  # Join removal status from 2021 census
  left_join(select(filter(mutate(read.csv(paste0(dir, "00_process_raw_data/census_main_core.csv")), 
                                 Year = as.factor(Year)), Year == 2021), 
                   Accession_progeny, Removed), by = "Accession_progeny") %>% 
  filter(Removed.y != 1, Year != 2013) %>%  # Exclude removed plants and 2013 data
  # Join predicted heights from allometry
  left_join(select(allom_dat_pred3, Accession_progeny, Year, ht_complete), 
            by = c("Accession_progeny", "Year")) %>%
  mutate(ht_complete = case_when(
    is.na(ht_complete) ~ 0,  # Dead/missing plants get height of 0
    .default = ht_complete
  ))

write.csv(long_core_relfit, paste0(dir, "00_process_raw_data/long_core_relfit.csv"))

# calculate multiplicative fitness function -------------------------------

# Count initial number of progeny per accession at IFG in 2014
accnum_ifg <- long_core_relfit %>%
  filter(Year == 2014, Site == "IFG") %>%
  arrange(Block, Accession) %>%
  group_by(Block, Accession) %>%
  summarize(n_initial = n())

# Count initial number of progeny per accession at Chico in 2014
accnum_cso <- long_core_relfit %>%
  filter(Year == 2014, Site == "Chico") %>%
  arrange(Block, Accession) %>%
  group_by(Block, Accession) %>%
  summarize(n_initial = n())

# Count survivors in 2024 and compare to initial counts
alive_24 <- long_core_relfit %>%
  filter(Year == 2024, Site == "IFG", Alive_new == 1) %>%
  summarize(n_24 = n(), .by = c(Site, Block, Accession)) %>%
  right_join(summarize(filter(long_core_relfit, Year == 2014), n_initial = n(), 
                       .by = c(Site, Block, Accession)), 
             by = c("Site", "Block", "Accession"))

## IFG site fitness calculations

# Calculate maximum fitness (wmax) per block per year at IFG
# This is the denominator for relative fitness
wmaxblock_ifg <- long_core_relfit %>%
  filter(Site == "IFG") %>% 
  filter(!(is.na(ht_complete) & Alive_new == 1)) %>%  # Remove alive plants missing height
  arrange(Block, Year) %>%
  group_by(Block, Year) %>%
  # Individual fitness (wij) = height if alive, 0 if dead
  mutate(wij = case_when(
    Alive_new == 1 ~ ht_complete,
    Alive_new == 0 ~ 0,
    .default = NA)) %>% 
  group_by(Year, Accession_progeny, Gopher_damage) %>% 
  filter(Gopher_damage != 1) %>%  # Exclude gopher-damaged plants
  ungroup() %>%
  # Sum fitness across all individuals of an accession in a block
  group_by(Year, Block, Accession) %>%
  summarize(wj = sum(wij),
            htj = sum(ht_complete)) %>%
  left_join(accnum_ifg, by = c("Block", "Accession")) %>%
  # Normalize by initial count (per capita fitness)
  mutate(wj = wj / n_initial,
         htj = htj / n_initial) %>%
  ungroup() %>%
  # Find maximum fitness in each block-year combination
  group_by(Year, Block) %>%
  filter(!is.na(wj), !is.na(htj)) %>% 
  summarize(wmax = max(wj),
            htmax = max(htj))

# Calculate relative fitness for each accession at IFG
# Relative fitness = accession fitness / maximum fitness in that block
relfit_ifg <- long_core_relfit %>%
  filter(Site == "IFG") %>%
  filter(!(is.na(ht_complete) & Alive_new == 1)) %>% 
  arrange(Block, Year) %>%
  group_by(Block, Year) %>%
  mutate(wij = case_when(
    Alive_new == 1 ~ ht_complete,
    Alive_new == 0 ~ 0,
    .default = NA)) %>% 
  group_by(Year, Accession_progeny, Gopher_damage) %>%
  filter(Gopher_damage != 1) %>%
  ungroup() %>%
  group_by(Year, Block, Accession) %>%
  summarize(wj = sum(wij),
            htj = sum(ht_complete)) %>% 
  left_join(accnum_ifg, by = c("Block", "Accession")) %>% 
  left_join(wmaxblock_ifg, by = c("Year", "Block")) %>%
  mutate(wj = wj / n_initial,
         htj = htj / n_initial) %>% 
  # Calculate relative fitness (rf) and relative height (rht)
  mutate(rf = wj / wmax,
         rht = htj / htmax) %>% 
  # Average across blocks for each accession-year
  group_by(Year, Accession) %>% 
  summarize(mrf = mean(rf),    # Mean relative fitness
            mrht = mean(rht)) %>%  # Mean relative height
  # Join climate data from maternal (origin) sites
  left_join(dat_clim_pc, by = "Accession") %>%
  ungroup() %>%
  arrange(Year, Accession) %>%
  group_by(Accession) %>%
  mutate(latlong = Latitude * Longitude,  # Interaction term
         lagmrf = lag(mrf))  # Previous year's fitness

## Chico site fitness calculations (same process as IFG)

# Calculate maximum fitness per block per year at Chico
wmaxblock_cso <- long_core_relfit %>%
  filter(Site == "Chico") %>%
  filter(!(is.na(ht_complete) & Alive_new == 1)) %>%
  arrange(Block, Year) %>%
  group_by(Block, Year) %>%
  mutate(wij = case_when(
    Alive_new == 1 ~ ht_complete,
    Alive_new == 0 ~ 0,
    .default = NA)) %>% 
  group_by(Year, Accession, Gopher_damage) %>% 
  filter(Gopher_damage != 1) %>% 
  group_by(Year, Block, Accession) %>%
  summarize(wj = sum(wij),
            htj = sum(ht_complete)) %>%
  left_join(accnum_cso, by = c("Block", "Accession")) %>%
  mutate(wj = wj / n_initial,
         htj = htj / n_initial) %>%
  ungroup() %>%
  group_by(Year, Block) %>%
  filter(!is.na(wj), !is.na(htj)) %>%
  summarize(wmax = max(wj),
            htmax = max(htj))

# Calculate relative fitness for each accession at Chico
relfit_cso <- long_core_relfit %>%
  filter(Site == "Chico") %>%
  filter(!(is.na(ht_complete) & Alive_new == 1)) %>%
  arrange(Block, Year) %>%
  group_by(Block, Year) %>%
  mutate(wij = case_when(
    Alive_new == 1 ~ ht_complete,
    Alive_new == 0 ~ 0,
    .default = NA)) %>% 
  group_by(Year, Accession, Gopher_damage) %>% 
  filter(Gopher_damage != 1) %>%
  ungroup() %>%
  group_by(Year, Block, Accession) %>%
  summarize(wj = sum(wij),
            htj = sum(ht_complete)) %>% 
  left_join(accnum_cso, by = c("Block", "Accession")) %>% 
  left_join(wmaxblock_cso, by = c("Year", "Block")) %>% 
  mutate(wj = wj / n_initial,
         htj = htj / n_initial) %>% 
  mutate(rf = wj / wmax,
         rht = htj / htmax) %>%
  group_by(Year, Accession) %>% 
  summarize(mrf = mean(rf),
            mrht = mean(rht)) %>%
  left_join(dat_clim_pc, by = "Accession") %>%
  ungroup() %>%
  arrange(Year, Accession) %>%
  group_by(Accession) %>%
  mutate(latlong = Latitude * Longitude,
         lagmrf = lag(mrf))

prismtmax$Year <- as.factor(prismtmax$Year)

## Combine both sites into single dataset
relfit_all <- bind_rows(relfit_cso, relfit_ifg, .id = "Site") %>% 
  mutate(Site = case_when(Site == 2 ~ "IFG",
                          Site == 1 ~ "Chico")) %>%
  left_join(prismtmax, by = c("Site", "Year")) %>%
  # Calculate climate mismatch variables
  mutate(tdiffyr = tmax_sum.y - tmax_sum.x,  # Garden temp - maternal temp (yearly)
         pdiffyr = ppt.y - ppt.x,              # Garden precip - maternal precip (yearly)
         # Latitude and longitude differences from garden sites
         latdiff = case_when(
           Site == "IFG" ~ 38.74 - Latitude,
           Site == "Chico" ~ 39.71 - Latitude,
           .default = NA
         ),
         longdiff = case_when(
           Site == "IFG" ~ -120.742112 - Longitude,
           Site == "Chico" ~ -121.782664 - Longitude,
           .default = NA
         ),
         latlongdiff = latdiff * longdiff,  # Geographic mismatch interaction
         Accession = as.factor(Accession),
         Locality = as.factor(Locality),
         Site = as.factor(Site)) %>%
  ungroup()


# Filter to 2017 data and calculate climate differences
relfit17 <- relfit_all %>%
  ungroup() %>%
  filter(Year == 2017) %>%
  mutate(tdiff2 = case_when(Site == "IFG" ~ ifg_tmax - tmax_sum.x,
                            Site == "Chico" ~ cso_tmax - tmax_sum.x),
         tdiff2 = case_when(Site == "IFG" ~ ifg_tmax_ave - tmax_sum.x,
                            Site == "Chico" ~ cso_tmax_ave - tmax_sum.x),
         pdiff = case_when(Site == "IFG" ~ ifg_ppt - ppt.x,
                           Site == "Chico" ~ cso_ppt - ppt.x)) %>%
  left_join(select(filter(relfit_all, Year == 2014), Accession, mrf, Site), by = c("Accession", "Site")) %>%
  rename("mrf" = mrf.x, "mrf_2014" = mrf.y)
write.csv(relfit17, paste0(dir, "mff17_input.csv"))

# Filter to 2024 data and calculate climate differences
relfit24 <- relfit_all %>%
  ungroup() %>%
  filter(Year == 2024) %>%
  mutate(
    # Temperature difference: garden average - maternal site (maximum metric)
    tdiff = case_when(Site == "IFG" ~ ifg_tmax - tmax_sum.x,
                      Site == "Chico" ~ cso_tmax - tmax_sum.x),
    # Temperature difference: garden average - maternal site (mean metric)
    tdiff2 = case_when(Site == "IFG" ~ ifg_tmax_ave - tmax_sum.x,
                       Site == "Chico" ~ cso_tmax_ave - tmax_sum.x),
    # Precipitation difference
    pdiff = case_when(Site == "IFG" ~ ifg_ppt - ppt.x,
                      Site == "Chico" ~ cso_ppt - ppt.x),
    # Winter temperature difference
    tdiff_win = case_when(Site == "IFG" ~ ifg_tavewin - bioclim_11,
                          Site == "Chico" ~ cso_tavewin - bioclim_11)
  ) %>%
  # Join initial fitness from 2014 as a covariate
  left_join(select(filter(relfit_all, Year == 2014), Accession, mrf, Site), 
            by = c("Accession", "Site")) %>%
  rename("mrf" = mrf.x, "mrf_2014" = mrf.y)

write.csv(relfit24, paste0(dir, "mff24_input.csv"))


## data transformations for Figure 3c
highlow<- ## top and bottom 10% of performers in 2017, the median year for tmax
  bind_rows(
    slice_max(filter(relfit_all, Year == 2017, Site == "Chico"), order_by = tdiffyr, prop = .1),
    slice_min(filter(relfit_all, Year == 2017, Site == "Chico"), order_by = tdiffyr, prop = .1),
    slice_max(filter(relfit_all, Year == 2017, Site == "IFG"), order_by = tdiffyr, prop = .1),
    slice_min(filter(relfit_all, Year == 2017, Site == "IFG"), order_by = tdiffyr, prop = .1)
  ) %>% mutate(group = case_when( 
    tdiffyr < 4 ~ "warm",
    tdiffyr > 4 ~ "cool"
  ),
  Accession = as.factor(Accession)) %>%
  select(Accession, group) %>%
  distinct()

dat_fig3c<-long_core_relfit %>%
  mutate(Accession = as.factor(Accession)) %>%
  filter(!is.na(Accession_progeny), !((Year == 2018 | Year == 2020) & Site == "Chico") ) %>%
  summarize(.by = c(Site, Year, Accession), ht_complete = mean(ht_complete)) %>%
  #filter(Accession %in% highlow$Accession) %>%
  inner_join(highlow, by = "Accession") %>%  ## compare the top and bottom 10% in all years
  mutate(Year = as.integer(as.character(Year)),
         Site = as.factor(Site)) %>%
  mutate(Site = factor(Site, levels = rev(levels(Site))))

write.csv(dat_fig3c, paste0(dir, "03_spatial_prediction_figures/dat_fig3c.csv"))
