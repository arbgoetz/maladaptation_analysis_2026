
# Convenience functions used in data cleaning/analysis --------------------




# read in libraries -------------------------------------------------------


library(tidyverse)
library(readxl)
library(janitor)
library(stringr)



# Functions for adjusting accession_progeny ID ----------------------------

apfix<-function(x) { ## remove any trailing dashes or locality IDs (for compatibility)
  case_when(
    str_sub(x, 1, 1)=="-" ~ sub('.', '', x),
    nchar(sub('-[^-]*$', '', x))>3~ sub('.*?-', '', x), 
    .default=as.character(x)
  ) 
  
}

apleadzero<-function(x) { ## add leading zeroes to accession_progeny to force format XXX-XX, return NA if not possible
  case_when(
    x == "Tree Removed" | x == "no tree" | x == "no tag" | x == "blank" ~ x,
    !is.na(as.numeric(sub('-.*', '', x))) & !is.na(as.numeric(sub('.*?-', '', x))) ~ str_c(
      sprintf("%03d", as.numeric(sub('-.*', '', x))), 
      "-",
      sprintf("%02d", as.numeric(sub('.*?-', '', x))) 
    ),
    !is.na(as.numeric(sub('-.*', '', x))) & grepl("a", sub('.*?-', '', x)) ~ str_c(
      sprintf("%03d", as.numeric(sub('-.*', '', x))), 
      "-",
      sprintf("%02d", as.numeric(sub('.*?-', '', sub('a.*?', '', x)))), "a"
    ),
    !is.na(as.numeric(sub('-.*', '', x))) & grepl("b", sub('.*?-', '', x)) ~ str_c(
      sprintf("%03d", as.numeric(sub('-.*', '', x))), 
      "-",
      sprintf("%02d", as.numeric(sub('.*?-', '', sub('b.*?', '', x)))), "b"
    ),
    .default=NA)
}



apaddmissing <- function(x, from) { ## "from" is the dataset to get missing a_ps from
  x %>% left_join(
    select(
      filter(from, !is.na(from$Accession_progeny)), Section, Row, Line, Block, Accession_progeny, Accession, Progeny), 
    by=c("Section", "Row", "Line")
  ) %>%
    mutate(
      Accession = case_when(
        (is.na(Accession.x) | is.na(Progeny.x)) & !is.na(Accession.y) ~ Accession.y,
        .default = Accession.x
      ),
      Progeny = case_when(
        (is.na(Accession.x) | is.na(Progeny.x)) & !is.na(Progeny.y) ~ Progeny.y,
        .default = Progeny.x
      ),
      Accession_progeny = case_when(
        (is.na(Accession.x) | is.na(Progeny.x)) & !is.na(Accession_progeny.y) ~ Accession_progeny.y,
        .default = Accession_progeny.x),
      .before = Height_cm
    ) %>%
    mutate(Block = case_when(
      is.na(Block.x) & !is.na(Block.y) ~ Block.y,
      .default = Block.x
    ), .before = Row) %>%
    select(-Accession_progeny.x, -Accession_progeny.y, -Accession.x, -Accession.y, -Progeny.x, -Progeny.y, -Block.x, -Block.y )
}

apaddmissing_lite <- function(x, from) { ## "from" is the dataset to get missing a_ps from
  x %>% left_join(
    select(
      filter(from, !is.na(from$Accession_progeny)), Section, Row, Line, Accession_progeny), 
    by=c("Section", "Row", "Line")
  ) %>%
    mutate(Accession_progeny = case_when(
        (is.na(Accession_progeny.x)) & !is.na(Accession_progeny.y) ~ Accession_progeny.y,
        .default = Accession_progeny.x)
    ) %>%
    select(-Accession_progeny.x, -Accession_progeny.y, )
}


# Function for adding missing Locality ------------------------------------

## Note: requires "dat_maternal" ("Maternal tree locations and sample ID", linked in data availability summary)
locaddmissing <- function(x) {
  x %>%
    left_join(select(dat_maternal, Accession, Locality), by="Accession", na_matches = "never") %>%
    mutate(Locality = case_when(
      !is.na(Locality.y) ~ Locality.y,
      .default = Locality.x),
      .before = "Accession") %>%
    select(-Locality.x, -Locality.y)
}



# Function for adding missing Block (CHICO ONLY) (improved) ----------------------------
## Doesn't require anything to be loaded in

blockaddmissing<-function(x) {
  x %>%
    row_fac() %>%
    mutate(Block = case_when(
      Section == "Block1" ~ "1",
    Section == "North annex" ~ "2",
    Section == "North" & Row >= "A" & Row <= "O" ~ "3", 
    Section == "North" & Row >= "P" & Row <= "Y" ~ "4",
    Section == "North" & Row >= "Z" & Row <= "AL" ~ "5",
    .default = Block
    ))
    
}

# Function for converting Row to factor -----------------------------------
row_fac<-function(x) {
  lts<-bind_rows(tibble(let = LETTERS[1:26]), tibble(let = paste0("A",LETTERS[1:26]))) %>% #creating a tibble with rows A-AZ
    mutate(let = factor(let))
  x %>%
    mutate(Row = factor(Row, levels = lts$let, ordered = TRUE))
}



# Function for standardizing border label (new) ---------------------------
## fixing some incorrect ones in 2017a and then applying that to all years
bdrfix<-function(x) {
  dat_17a <- row_fac(dat_17a)
  x <- row_fac(x)
  
  link_tbl <- select(dat_17a, Section, Row, Line, Border) %>% ## link table is based on map: https://docs.google.com/spreadsheets/d/1h3T4wFQ43G3kNoldgOI1DPbz0JKmKAHd/edit#gid=481360690
    mutate(Border = case_when(
      Section == "North" & Row == "A" ~ 1,
      Section == "North" & Row == "AL" ~ 1,
      Section == "North" & Line == 1 ~ 1,
      Section == "North" & Line == 69 ~ 1,
      Section == "North" & Row == "O" & Line >= 53 ~ 1,
      Section == "North" & Row == "AJ" & Line >= 44 ~ 1,
      Section == "North" & Row >= "O" & Line == 53 ~ 1,
      Section == "North" & Row == "AK" & Line == 44 ~ 1,
      
      Section == "North annex" & Row == "A" ~ 1,
      Section == "North annex" & Row == "Z" ~ 1,
      Section == "North annex" & Line == 1 ~ 1,
      Section == "North annex" & Line == 24 ~ 1,
      Section == "North annex" & Row == "N" & Line >= 20 ~ 1,
      Section == "North annex" & Row == "X" & Line >= 12 ~ 1,
      Section == "North annex" & Row == "V" & Line >= 15 ~ 1,
      Section == "North annex" & Row >= "N" & Line == 20 ~ 1,
      Section == "North annex" & Row == "W" & Line == 15 ~ 1,
      Section == "North annex" & Row == "Y" & Line == 12 ~ 1,
      
      Section == "Block1" & Row == "A" ~ 1,
      Section == "Block1" & Row == "Q" ~ 1,
      Section == "Block1" & Line == 1 ~ 1,
      Section == "Block1" & Line == 47 ~ 1,
      Section == "Block1" & Line == 1 ~ 1,
      Section == "Block1" & Row == "B" & Line >= 32 ~ 1,
      .default = 0
    ))
  
  x %>%
    left_join(link_tbl, by=c("Section", "Row", "Line")) %>%
    mutate(Border = case_when(
      !is.na(Border.y) ~ Border.y,
       .default = NA),
       .after = Line) %>%
    select(-Border.x, -Border.y)
}

# Function for survival back-dating (clumsy) -------------------------------


resurrect<-function(x) {
  x %>%
    filter(!is.na(Accession)) %>%
    mutate(Alive = replace_na(Alive, 0)) %>%
    select(Accession_progeny, Year, Alive) %>%
    pivot_wider(id_cols = c(Accession_progeny), values_from = Alive, names_from = Year, names_prefix = "Alive_", id_expand = FALSE) %>%
    mutate(Alive_new_2024 = case_when(
      is.na(Alive_2024) ~ 0,
      .default = Alive_2024
    )) %>%
    mutate(Alive_new_2023 = case_when(
      Alive_new_2024 == 1 ~ 1,
      Alive_new_2024 == 0 ~ Alive_2023,
      is.na(Alive_new_2024) ~ Alive_2023
    ))  %>%
    mutate(Alive_new_2022 = case_when(
      Alive_new_2023 == 1 ~ 1,
      Alive_new_2023 == 0 ~ Alive_2022,
      is.na(Alive_new_2023) ~ Alive_2022
    ))  %>%
    mutate(Alive_new_2021 = case_when(
      Accession_progeny == "362-04" ~ 0,
      Accession_progeny == "641-10" ~ 0,
      Alive_new_2022 == 1 ~ Alive_new_2022,
      Alive_new_2022 == 0 ~ Alive_2021,
      is.na(Alive_new_2022) ~ Alive_2021
    )) %>%
    mutate(Alive_new_2020 = case_when(
      Alive_new_2021 == 1 ~ Alive_new_2021,
      Alive_new_2021 == 0 ~ Alive_2020,
      is.na(Alive_new_2021) ~ Alive_2020
    )) %>%
    mutate(Alive_new_2019 = case_when( ##need to verify that everything works for 2020 and 2019
      Alive_new_2020 == 1 ~ Alive_new_2020,
      Alive_new_2020 == 0 ~ Alive_2019,
      is.na(Alive_new_2020) ~ Alive_2019
    )) %>%
    mutate(Alive_new_2018 = case_when(
      Alive_new_2019 == 1 ~ Alive_new_2019,
      Alive_new_2019 == 0 ~ Alive_2018,
      is.na(Alive_new_2019) ~ Alive_2018
    )) %>%
    mutate(Alive_new_2017 = case_when(
      Accession_progeny == "062-04" ~ 0,
      Accession_progeny == "022-08" ~ 0,
      Alive_new_2018 == 1 ~ Alive_new_2018,
      Alive_new_2018 == 0 ~ Alive_2017,
      is.na(Alive_new_2018) ~ Alive_2017
    )) %>%
    mutate(Alive_new_2016 = case_when(
      Accession_progeny == "381-04" ~ 0,
      Accession_progeny == "022-08" ~ 0,
      Accession_progeny == "474-13" ~ 1,
      Accession_progeny == "477-12" ~ 1,
      is.na(Alive_2016)  ~ case_when(
        Alive_new_2017 == 1 & Alive_2015 == 1 ~ 1,
        Alive_new_2017 == 0 & Alive_2015 == 0 ~ 0,
        .default = 1),
      Alive_new_2017 == 1 ~ Alive_new_2017,
      Alive_new_2017 == 0 | is.na(Alive_new_2017) ~ Alive_2016,
      .default = NA
    )
    ) %>%
    mutate(Alive_new_2015 = case_when(
      Accession_progeny == "381-04" ~ 0,
      Alive_new_2016 == 1 ~ Alive_new_2016,
      Alive_new_2016 == 0 ~ Alive_2015,
      is.na(Alive_new_2016) ~ Alive_2015
    )) %>%
    mutate(Alive_new_2014 = case_when(
      Alive_new_2015 == 1 ~ Alive_new_2015,
      Alive_new_2015 == 0 ~ Alive_2014,
      is.na(Alive_new_2015) ~ Alive_2014 
    )) %>%
    mutate(Alive_new_2013 = case_when(
      Alive_new_2014 == 1 ~ Alive_new_2014,
      Alive_new_2014 == 0 ~ Alive_2013,
      is.na(Alive_new_2014) ~ Alive_2013
    )) %>%
    pivot_longer(cols=c(Alive_new_2013, 
                        Alive_new_2014, 
                        Alive_new_2015, 
                        Alive_new_2016, 
                        Alive_new_2017, 
                        Alive_new_2018, 
                        Alive_new_2019, 
                        Alive_new_2020,
                        Alive_new_2021, 
                        Alive_new_2022,
                        Alive_new_2023
                        ,
                        Alive_new_2024
    ), 
    names_to = "Year", names_prefix = "Alive_new_",  values_to = "Alive_new") %>%
    
    select(Accession_progeny, Year, Alive_new) %>%
    #mutate(Year) = as.numeric(Year) %>%
    right_join(x, by = c("Accession_progeny", "Year")) %>%
    mutate(Alive_new = case_when(
      is.na(Accession) ~ 0, ##changed from ~ Alive - the only ones with NA accession were recorded as Alive=1 for some reason
      (Year == 2013 | Year == 2014) & (!is.na(Height_cm) & !is.na(BDiam_mm)) ~ 1,
      .default = Alive_new
    )) %>%
    relocate(Accession_progeny, .after = Accession) %>%
    relocate(Year, .after = Problem_fixed) %>%
    relocate(Alive_new, .after = Alive) %>%
    arrange(Year, Accession_progeny)
  
}


# Function for seeing what percentage of height values represent actual measurements  --------
ht_check <- function(yr, site) {
  (dim(filter(long_core, Year == yr, Site == site, Height_cm < 999))[1]) / dim(filter(long_core, Year == yr, Site == site, !is.na(Height_cm)))[1]
}

ht_checkv<-Vectorize(ht_check)

# ht_checks<-tibble(Year = c(2013, 2014, 2015, 2016, 2017, 2018, 2019, 2021, 2022, 2023)) %>%
#   mutate(htpctIFG = ht_checkv(Year, "IFG"),
#          htpctCSO = ht_checkv(Year, "Chico"))
# View(ht_checks)

# Year htpctIFG htpctCSO
# 1  2013  NaN      NaN    
# 2  2014    1        1    
# 3  2015    1        1    
# 4  2016    1        1    
# 5  2017    1        1    
# 6  2018    1        0.848
# 7  2019    1        0.827
# 8  2021    0.560    0.285
# 9  2022    0.268    0.124
# 10 2023    0.434    0.156

##create cohorts from the set tdiff breaks

set_cohorts <- function(x, site) {
  set.seed(1000)
  
  ranges <- case_when(
    site == "Chico" ~  c(0, 0, 2, 3, 5, 5, 7, 7),
    site == "IFG" ~ c(-2, -2, 0, 0, 2, 2, 5, 5)
  )
  yr<-case_when(
    site == "Chico" ~ 2023,
    site == "IFG" ~ 2019
  )
  
  x %>%
    ungroup() %>%
      filter(Year == yr, Site == site) %>%
      mutate(cohort = case_when(
        tdiffyr <= ranges[1] ~ "A",
        tdiffyr > ranges[2] & tdiffyr <= ranges[3] ~ "B",
        tdiffyr > ranges[4] & tdiffyr <= ranges[5] ~ "C",
        tdiffyr > ranges[6] & tdiffyr <= ranges[7] ~ "D",
        tdiffyr > ranges[8] ~ "E",
        .default = NA
      )) %>%
      select(Accession, cohort) %>% 
      group_by(cohort) %>%
      slice_sample(n = 50) %>%
      ungroup() %>%
      right_join(filter(x, Site == site), by = "Accession")  %>% 
      filter(!is.na(cohort))
  return(x)
}

rgr_mean <- function(yr_init, yr_end, data) {
  ifg_tmax<-mean(filter(prismtmax, Site == "IFG", Year >= yr_init, Year <= yr_end)$tmax_sum)
  cso_tmax<-mean(filter(prismtmax, Site == "Chico", Year >= yr_init, Year <= yr_end)$tmax_sum)
  data2<-data %>%
    filter(Year == yr_init | Year == yr_end, Alive_new == 1) %>%
    mutate(loght = log(ht_complete)) %>%
    arrange(Accession_progeny, Year) %>%
    group_by(Accession_progeny) %>%
    summarize(dif = diff(loght)) %>%
    ungroup() %>%
    mutate(rgr = dif / (yr_end - yr_init)) %>%
    left_join(filter(data, Year == yr_end), by = "Accession_progeny") %>%
    select(-dif) %>%
    filter(!is.na(rgr), rgr > 0) %>%
    left_join(select(filter(long_core, Year == 2014), Accession_progeny, ht), by = "Accession_progeny") %>%
    rename("ht" = ht.x,
           "ht_2014" = ht.y) %>%
    left_join(dat_clim_pc, by = c("Accession", "Locality")) %>%
    mutate(tdiff = case_when(Site == "IFG" ~ ifg_tmax - tmax_sum,
                             Site == "Chico" ~ cso_tmax - tmax_sum))
  data2$Locality<-as.factor(data2$Locality)
  mean <- mean(data2$rgr)
 return(data2)
}
