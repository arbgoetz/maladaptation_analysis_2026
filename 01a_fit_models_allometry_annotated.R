# Set working directory
dir<-"C:/Users/alexg/Downloads/31937499/"

# Source custom helper functions
source(paste0(dir, "00_convenience_functions.R"))


# load in and clean data ------------------------------------------------------------

# Load main census data with height measurements
long_core <- read.csv(paste0(dir, "00_process_raw_data/census_main_core.csv")) %>%
  rename("ht" = Height_cm)  # Standardize height column name

## Placerville garden (IFG site) biomass data
biomass_ifg <- read_xlsx(paste0(dir, "00_process_raw_data/biomass_ifg.xlsx")) %>%
  separate(Coordinate, into = c("Row", "Line"), sep = " _ ") %>%  # Split coordinate into row/line
  separate(Tree_ID, into = c("Locality", "Accession"), sep = " _ ") %>%  # Split tree ID
  unite("Accession_progeny", c(Accession, Progeny), sep = "-", remove = FALSE) %>%  # Create unique ID
  # Rename biomass measurement columns
  rename("shoot_kg" = shoot_mass_kg,
         "root_kg" = root_mass_kg,
         "cookie_g" = cookie_mass_g) %>%
  mutate(
    diam_y_mm = as.numeric(diam_y_mm),  # Convert to numeric
    Site = "IFG",
    Section = "IFG",
    Accession_progeny = apleadzero(Accession_progeny),  # Add leading zeros (custom function)
    # Calculate maximum diameter (stems are measured in two perpendicular directions)
    diam_max = case_when(
      diam_x_mm > diam_y_mm ~ diam_x_mm,
      diam_y_mm > diam_x_mm ~ diam_y_mm,
      .default = NA),
    # Calculate cross-sectional area assuming ellipse
    area = (diam_x_mm / 2) * (diam_y_mm / 2) * pi) %>%
  select(-section) %>%
  # Remove records with missing IDs or problematic individual
  filter(!is.na(Accession), !is.na(Accession_progeny),
         Accession_progeny != "289-03")

## Chico garden (CSO site) biomass data
biomass_cso <- read_xlsx(paste0(dir, "00_process_raw_data/biomass_cso.xlsx")) %>%
  separate(Coordinate, into = c("Row", "Line"), sep = " _ ") %>%
  separate(`Tree ID`, into = c("Locality", "Accession"), sep = " _ ") %>%
  unite("Accession_progeny", c(Accession, Progeny), sep = "-", remove = FALSE) %>%
  # Rename columns to match IFG dataset
  rename("diam_x_mm" = cookie_mm_x,
         "diam_y_mm" = cookie_mm_y,
         "ht" = shoot_cm) %>%
  mutate(Accession_progeny = apleadzero(Accession_progeny),
         Site = "Chico",
         # Standardize section names
         Section = case_when(
           section == "n" ~ "North",
           section == "b" ~ "Block1",
           section == "x" ~ "North annex",
           .default = NA
         ),
         bark_mm = as.numeric(bark_mm),
         # Calculate maximum diameter
         diam_max = case_when(
           diam_x_mm > diam_y_mm ~ diam_x_mm,
           diam_y_mm > diam_x_mm ~ diam_y_mm,
           .default = NA
         ),
         # Calculate cross-sectional area
         area = (diam_x_mm / 2) * (diam_y_mm / 2) * pi) %>%
  # Remove unnecessary columns
  select(-section, -`...16`, -stem_pre, -stem_post, -root_pre, -root_post) %>%
  # Filter out missing IDs, zero heights, and outliers
  filter(!is.na(Accession), !is.na(Accession_progeny),
         ht > 0,
         !(ht > 200 & shoot_kg < .05),  # Remove tall plants with suspiciously low biomass
         !(log(ht) < 4.8 & log(shoot_kg > -1)),  # Remove biologically implausible combinations
         Accession_progeny != "373-11")  # Remove problematic individual

# Combine both sites' biomass data
biomass_all <- bind_rows(biomass_ifg, biomass_cso) %>%
  # Join multistem status from 2018 census data
  left_join(select(filter(long_core, Year == 2018), Accession_progeny, Multistem_yn), 
            by = "Accession_progeny") %>%
  mutate(Multistem_yn = as.factor(Multistem_yn))

# models ------------------------------------------------------------------

## Model 1: Predict basal diameter from DBH (diameter at breast height)
# This allows conversion between the two diameter measurement types
bdfromdbh <- lm(BDiam_mm ~ DBH_1, data = long_core)

## Model 2: Predict height from diameter using log-log relationship
# Log transformation linearizes the allometric power-law relationship
htfromdbh_lin <- lm(log(ht) ~ log(diam_max), 
                    data = filter(biomass_all, ht > 0, !is.na(diam_max)))

## Model 3: Predict aboveground biomass from basal diameter
# Log-log relationship is standard for biomass allometry
agbfrombd_all <- lm(log(shoot_kg) ~ log(diam_max), data = biomass_all)  # R² = 0.84
write_rds(agbfrombd_all, paste0(dir, "00_process_raw_data/agbfrombd_all.rds"))

# add modeled columns to data frame ---------------------------------------

# Apply allometric models to fill in missing data
allom_dat2 <- long_core %>%
  mutate(Year = as.factor(Year),
         # Calculate or estimate maximum diameter
         diam_max = case_when(
           # Manual correction for specific individual
           Accession_progeny == "028-10" & Year == 2018 ~ 55,
           # Convert DBH to basal diameter using fitted model coefficients
           # Formula from bdfromdbh: BDiam = (DBH * 1.42) + 27.94
           is.na(BDiam_mm) ~ (DBH_1 * 1.42) + 27.94,
           # Use measured basal diameter if available
           !is.na(BDiam_mm) ~ BDiam_mm,
           .default = NA
         ),
         # Clean height data (999 used as missing value code)
         ht = case_when(
           ht >= 999 ~ NA,
           .default = ht
         )) 

# Create complete height dataset by filling missing values with predictions
allom_dat_pred3 <- allom_dat2 %>%
  mutate(
    # Predict height from diameter using fitted model
    # exp() back-transforms from log scale
    ht_pred_dbh = exp(predict(htfromdbh_lin, allom_dat2)),
    # Use measured height if available, otherwise use predicted height
    ht_complete = coalesce(ht, ht_pred_dbh)
  ) %>% 
  as_tibble()

# Save complete dataset with predicted heights
write.csv(allom_dat_pred3, paste0(dir, "00_process_raw_data/allometry_pred_hts.csv"))
