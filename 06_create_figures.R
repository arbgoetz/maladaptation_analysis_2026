# Set working directory
dir<-"C:/Users/alexg/Downloads/31937499/"

# Source custom helper functions
source(paste0(dir, "00_convenience_functions.R"))

# Fig 1: map of origin sites ----------------------------------------------

## see 07_create_maps.R

# Fig 2: cumulative growth metrics ----------------------------------------


# libraries ---------------------------------------------------------------

library(readr)
library(ggplot2)
library(mgcv)
library(gratia)
library(tidygam)
library(multcomp)
library(viridis)
library(conflicted)
library(forcats)
library(terra)
library(tidyterra)
library(sf)
library(cowplot)
library(ggeffects)
library(flextable)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("annotate", "ggplot2")
library(GauPro)
library(Hmisc)
conflicts_prefer(dplyr::lag, dplyr::summarize)

theme_set(theme_bw())

# data --------------------------------------------------------------------

## variables

site_clim_avgs<-read.csv(paste0(dir, "03_spatial_prediction_figures/site_clim_avgs.csv")) ## calculated in 01b_spatial_calculations.R
ifg_tmax<-site_clim_avgs[1,3]
ifg_tmax_ave<-site_clim_avgs[2,3]
ifg_tavewin<-site_clim_avgs[4,3]
cso_tmax<-site_clim_avgs[5,3]
cso_tmax_ave<-site_clim_avgs[6,3]
cso_tavewin<-site_clim_avgs[8,3]
site_clim_avgs_r<-read.csv(paste0(dir, "03_spatial_prediction_figures/garden_aves_rast.csv")) ## calculated in 01b_spatial_calculations.R
bio11_45_gardenave<-site_clim_avgs_r[1,3]
tmaxpred_sum_gardenave<-site_clim_avgs_r[2,3]
bio11_lgm_gardenave<-site_clim_avgs_r[3,3]

tdiffrange<-seq(-20, 10, length.out = 100)



## input data
dat_rgr7<-read.csv(paste0(dir, "rgryr_input.csv"))
dat_rgr_htc<-read.csv(paste0(dir, "rgrwide_input.csv"))

## model predictions
### cumulative
gpc24<-read.csv(paste0(dir, "03_spatial_prediction_figures/predictions/cumulative_2014_2024/mff24_pred_gau.csv"))
gpc24_2<-read.csv(paste0(dir, "03_spatial_prediction_figures/predictions/cumulative_2014_2024/mff24_tave_pred_gau.csv"))
gp1424<-read.csv(paste0(dir, "03_spatial_prediction_figures/predictions/cumulative_2014_2024/rgr24_pred_gau.csv"))
gp1424_2 <- read.csv(paste0(dir, "03_spatial_prediction_figures/predictions/cumulative_2014_2024/rgr24_tave_pred_gau.csv"))

p_gamrf24<-read.csv(paste0(dir, "03_spatial_prediction_figures/predictions/cumulative_2014_2024/mff24_pred.csv"))
p_gamrf24_2<-read.csv(paste0(dir, "03_spatial_prediction_figures/predictions/cumulative_2014_2024/mff24_tave_pred.csv"))
p_gamrf24_win<-read.csv(paste0(dir, "03_spatial_prediction_figures/predictions/cumulative_2014_2024/mff24_win_pred.csv"))
p_1424_gam<-read.csv(paste0(dir, "03_spatial_prediction_figures/predictions/cumulative_2014_2024/rgr24_pred.csv"))
p_1424_gam2<-read.csv(paste0(dir, "03_spatial_prediction_figures/predictions/cumulative_2014_2024/rgr24_tave_pred.csv"))
p_1424_gam2_win<-read.csv(paste0(dir, "03_spatial_prediction_figures/predictions/cumulative_2014_2024/rgr24_win_pred.csv"))

### between-year
gau_pred_r<-read.csv(paste0(dir, "03_spatial_prediction_figures/predictions/single_year/rgryr_pred_gau.csv")) %>%
  mutate(Year = as.factor(Year)) %>%
  mutate(Year_bytemp = fct_reorder(Year, mtemp))
gau_pred_r2<-read.csv(paste0(dir, "03_spatial_prediction_figures/predictions/single_year/rgryr_tave_pred_gau.csv"))

p_gam1_temp<-read.csv(paste0(dir, "03_spatial_prediction_figures/predictions/single_year/rgryr_pred.csv"))
p_gam1_temp2<-read.csv(paste0(dir, "03_spatial_prediction_figures/predictions/single_year/rgryr_tave_pred.csv")) %>%
  mutate(Year_bytemp = fct_reorder(as.factor(Year), mtemp))





# Fig 1 -------------------------------------------------------------------

## see 07_create_maps.R

# Fig 2 -------------------------------------------------------------------

rgr_pts<-dat_rgr_htc %>%
  filter(RGRh1424 > 0) %>%
  summarize(.by = c(Accession),
            tdiff = mean(tdiff, na.rm = T),
            tdiff2 = mean(tdiff2, na.rm = T),
            tdiff_win = mean(tdiff_win, na.rm = T),
            RGRh1424 = mean(RGRh1424, na.rm = T))

##A

fig_rgr24<-ggplot(gp1424, aes(x = tdiff, y = RGRh1424)) + 
  
  geom_point(data = filter(rgr_pts, RGRh1424 > 0), aes(x = tdiff, y = RGRh1424), color = "grey") +
  geom_line(linetype = "dashed", color = "#c5283d", linewidth = .8) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), color = NA, fill = "#c5283d", alpha = .1) +
  geom_line(data = p_1424_gam, color = "#c5283d", linewidth = .8) + 
  geom_ribbon(data = p_1424_gam, aes(ymin = lower_ci, ymax = upper_ci), color = NA, fill = "#c5283d", alpha = 0.2)  +
  #geom_point(data = rgr_pts, aes(x = tdiff, y = RGRh1424), color = "grey") +
  scale_color_viridis(begin = .25, end = .75, name = "Tmin diff. (°C)", option = "magma") +
  #scale_fill_viridis(begin = .25, end = .75, name = "Annual maximum\ntemperature,\naveraged across\ngardens (°C)", option = "magma") +
  scale_x_continuous(limits = c(-11, 9)) + 
  #scale_y_continuous(limits = c(0, 1)) +
  geom_vline(xintercept = -5.2, linewidth = .2) + ## tdiff from LGM (source: Browne et al. 2019 PNAS)
  geom_vline(xintercept = tmaxpred_sum_gardenave - mean(c(ifg_tmax, cso_tmax)), linewidth = .2) +
  geom_vline(xintercept = 0, linetype = "longdash") +
  theme_bw() +
  labs(x = "Maximum summer temperature, transfer distance (°C)",
       y = "Relative growth rate 2014-2024")


##B

fig_rgr24_win<-ggplot(p_1424_gam2_win, aes(x = tdiff_win, y = RGRh1424)) + 
  geom_point(data = rgr_pts, aes(x =  tdiff_win, y = RGRh1424), color = "grey") +
  geom_line(color = "#255f85", linewidth = .8) + geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "#255f85", alpha = .2) + 
  
  #geom_point(data = filter(dat_rgr_htc, RGRh1424 > 0), aes(x = tdiff_win, y = RGRh1424)) +
  geom_vline(xintercept = 0, linetype = "longdash") + 
  geom_vline(xintercept = bio11_45_gardenave - mean(c(ifg_tavewin, cso_tavewin)), linewidth = .2) +
  geom_vline(xintercept = bio11_lgm_gardenave - mean(c(ifg_tavewin, cso_tavewin)), linewidth = .2 ) +
  labs(x = "Mean temperature of coldest quarter, transfer distance (°C)",
       y = "Relative growth rate 2014-2024") +
  theme_bw()



##C

rf_pts<-relfit24 %>%
  filter(mrf > 0) %>%
  summarize(.by = Accession,
            mrf = mean(mrf, na.rm = T),
            tdiff = mean(tdiff, na.rm = T),
            tdiff2 = mean(tdiff2, na.rm = T),
            tdiff_win = mean(tdiff_win, na.rm = T))


fig_rf24<-ggplot(gpc24, aes(x = tdiff, y = mrf)) + 
  geom_point(data = rf_pts, aes(x = tdiff, y = mrf
                                #, color = tdiff_win
  ), color = "grey") +
  geom_line(linetype = "dashed", color = "#c5283d", linewidth = .8) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), color = NA, fill = "#c5283d", alpha = .1) +
  geom_line(data = p_gamrf24, color = "#c5283d", linewidth = .8) + 
  geom_ribbon(data = p_gamrf24, aes(ymin = lower_ci, ymax = upper_ci), color = NA, fill = "#c5283d", alpha = 0.2)  +
  
  #scale_color_viridis( name = "Tmin diff. (°C)", option = "magma") +
  #scale_fill_viridis(begin = .25, end = .75, name = "Annual maximum\ntemperature,\naveraged across\ngardens (°C)", option = "magma") +
  scale_x_continuous(limits = c(-11, 9)) + 
  #cale_y_continuous(limits = c(0, 1)) +
  geom_vline(xintercept = -5.2, linewidth = .2) +
  geom_vline(xintercept = tmaxpred_sum_gardenave - mean(c(ifg_tmax, cso_tmax)), linewidth = .2) +
  geom_vline(xintercept = 0, linetype = "longdash") +
  theme_bw() +
  labs(x = "Maximum summer temperature, transfer distance (°C)",
       y = "Relative fitness, 2024")


##D

fig_rf24_win<-ggplot(p_gamrf24_win, aes(x = tdiff_win, y = mrf)) + 
  geom_point(data = rf_pts, aes(x = tdiff_win, y = mrf
                                #, color = tdiff_win
  ), color = "grey") +
  geom_line(color = "#255f85", linewidth = .8) + 
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "#255f85", alpha = .2) +
  geom_vline(xintercept = 0, linetype = "longdash") + 
  geom_vline(xintercept = bio11_45_gardenave - mean(c(ifg_tavewin, cso_tavewin)), linewidth = .2) +
  geom_vline(xintercept = bio11_lgm_gardenave - mean(c(ifg_tavewin, cso_tavewin)), linewidth = .2 ) +
  labs(x = "Mean temperature of coldest quarter, transfer distance (°C)",
       y = "Relative fitness, 2024") +
  theme_bw()


## Create and save multi-panel figure

cowplot::plot_grid(
  fig_rgr24 + scale_y_continuous(limits = c(.17, .28)) + 
    theme(axis.title.x=element_blank(),
          #axis.text.x=element_blank(),
          #axis.ticks.x=element_blank()
    ),
  fig_rgr24_win + scale_y_continuous(limits = c(.17, .28)) +
    scale_x_continuous(limits = c(-11, 9)) +
    theme(axis.title.x=element_blank(),
          #axis.text.x=element_blank(),
          #axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
    ),
  fig_rf24 + scale_y_continuous(limits = c(.25, 1.05)), 
  fig_rf24_win +
    scale_x_continuous(limits = c(-11, 9)) + 
    scale_y_continuous(limits = c(.25, 1.05)) +
    theme(axis.title.y = element_blank())
  ,
  
  # theme(axis.title.y=element_blank(),
  
  align = "hv", labels = "AUTO", label_x = .15, label_y = .95
)
#ggsave("~/git/demographyMS/figures/2026/fig2.png", height = 7.5, width = 10.12, units = "in", type = "cairo")



# Fig 3: yearly growth metrics --------------------------------------------


# Fig 3 ---------------------------------------------------------------



ggplot(gau_pred_r, aes(x = tdiffyr, y = rgrh, group = Year, color = mtemp)) + 
  geom_line(linetype = "dashed", linewidth = .8) +
  geom_ribbon(aes(ymin = lower, ymax = upper, group = Year, fill = mtemp), color = NA, alpha = .1) +
  geom_line(data = filter(p_gam1_temp, Year != 2014, Year != 2018, tdiffyr > -2), aes(color = mtemp), linewidth = .8) + 
  geom_ribbon(data = filter(p_gam1_temp, Year != 2014, Year != 2018, tdiffyr > -2), aes(ymin = lower_ci, ymax = upper_ci, group = Year, fill = mtemp), color = NA, alpha = .2)  +
  # scale_linetype_manual(values = c("y" = "solid", "n" = "dotted")) +
  # scale_linewidth_manual(values = c("y" = .8, "n" = .6)) +
  #scale_alpha_manual(values = c("y" = .2, "n" = 0)) +
  scale_color_viridis(begin = .25, end = .75, name = "Annual maximum\ntemperature,\naveraged across\ngardens (°C)", option = "magma") +
  scale_fill_viridis(begin = .25, end = .75, name = "Annual maximum\ntemperature,\naveraged across\ngardens (°C)", option = "magma") +
  scale_x_continuous(limits = c(-7, 9)) + 
  coord_cartesian(ylim = c(.1, .65)) +
  facet_wrap(~Year_bytemp, nrow = 1) +
  geom_vline(xintercept = 0, linetype = "longdash") +
  geom_text(aes(label = Year, x = 8.5, y = .6, hjust = 1, vjust = 0), color = "black", check_overlap = TRUE) +
  theme_bw() +
  labs(x = "Yearly temperature transfer distance (°C)",
       y = "Relative growth rate\n(single-year)",
       facet = "Year") +
  guides(linetype = "none",
         linewidth = "none",
         alpha = "none") +
  theme(strip.text = element_blank(),
        legend.position = "bottom")
#ggsave("~/git/demographyMS/figures/2026/fig3a.png", width = 12, height = 2.75)

dat_rgr_avg_yr<-dat_rgr7 %>%
  filter(!is.na(Accession_progeny),  !((Year == 2018 | Year == 2020) & Site == "Chico"), !(Year %in% c(2014, 2015, 2016)) ) %>%
  mutate(Accession = as.factor(Accession)) %>%
  summarize(
    .by = c(Site, Year),
    tmax_sum.y = mean(tmax_sum.y),
    rgrh = mean(rgrh)
  )

ggplot(dat_rgr_avg_yr, aes(x = tmax_sum.y, y = rgrh, color = tmax_sum.y)) + 
  stat_smooth(method = "lm", color = "black") +
  geom_point(aes(shape = Site), size = 3) +
  scale_color_viridis(begin = .25, end = .75, option = "magma") +
  labs(
    x = "Annual maximum temperature at garden (°C)",
    y = "Mean annual relative growth rate"
  ) +
  guides(color = "none",
         shape = "none")
##ggsave("~/git/demographyMS/figures/2026/fig3b_2.png", width = 5, height = 4, scale = .9)


#read.csv(paste0(dir, "03_spatial_prediction_figures/dat_fig3c.csv")
  ggplot( dat_fig3c, aes(x = Year, y = ht_complete, color = group, group = group)) + 
  stat_summary(geom = "line") +  
  stat_summary(geom = "pointrange") +
  facet_wrap(~Site) +
  scale_color_manual(values = c( "#2166AC",  "#B2182B")) + 
  #scale_color_viridis_d(option = "magma", begin = .3, end = .8, direction = -1) +
  #geom_vline(data = mutate(prismtmax, Year = as.factor(Year)), aes(xintercept = Year, color = tmax_sum)) +
  scale_x_continuous(minor_breaks = c(2014:2024)) +
  #scale_color_viridis_d() +
  labs(y = "Height\'") +
  guides(color = "none") +
  theme(strip.text = element_blank())

#ggsave("~/git/demographyMS/figures/2026/fig3c.png", width = 7.5, height = 4, scale = .9)

  

# Table S11 ----------------------------------------------------------------

  
  dat_rgr_lag<- dat_rgr7 %>%
    filter(!is.na(Accession_progeny),  !((Year == 2018 | Year == 2020) & Site == "Chico"), !(Year %in% c(2015, 2016)) ) %>%
    group_by(Accession_progeny) %>%
    mutate(Accession = as.factor(Accession),
           lag_temp = lag(tmax_sum.y)) %>%
    ungroup()
  

  
  ## table s11a
  lm(rgrh ~ tdiffyr*tmax_sum.y + tdiffyr*lag_temp + as.factor(Block), data = filter(dat_rgr_lag, Site == "IFG")) %>% 
    anova() %>%
    as_flextable()
  
  ## table s11b
  lm(rgrh ~ tdiffyr*tmax_sum.y + tdiffyr*lag_temp + as.factor(Block), data = filter(dat_rgr_lag, Site == "Chico")) %>%
    anova() %>%
    as_flextable()  

# Fig S1: precipitation over time -----------------------------------------

mppt<-prismtmax %>% group_by(Year) %>% summarize(ppt = mean(ppt)) %>%
  mutate(Site = "Mean")

prismtmax %>%
  mutate(Site = case_when(
    Site == "Chico" ~ "Warmer",
    Site == "IFG" ~ "Cooler"
  )) %>%
  #bind_rows(mppt) %>%
  ggplot(aes(x = Year, y = ppt, color = Site, group = Site)) +
  geom_line() +
  geom_line(data = mppt, aes(x = Year, y = ppt, color = Site, group = Site),   linetype = "dashed")+
  scale_color_discrete(direction = -1) +
  theme_bw() +
  #scale_x_continuous(breaks = c(2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024)) +
  labs(x = "Year",
       y = "Total annual precipitation (mm)",
       color = "Garden")
#ggsave("~/git/demographyMS/figures/figS1.png", height = 3, width = 6, units = "in")

# Fig S2: RGR and MFF tave cumulative ---------------------------------------------

##A


fig_rgr24_tave<-ggplot(gp1424_2, aes(x = tdiff2, y = RGRh1424)) + 
  geom_line(linetype = "dashed", color = "#c5283d", linewidth = .8) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), color = NA, fill = "#c5283d", alpha = .1) +
  geom_line(data = p_1424_gam2, color = "#c5283d", linewidth = .8) + 
  geom_ribbon(data = p_1424_gam2, aes(ymin = lower_ci, ymax = upper_ci), color = NA, fill = "#c5283d", alpha = 0.2)  +
  scale_x_continuous(limits = c(-20, 9)) + 
  geom_vline(xintercept = -5.2, linewidth = .2) +
  geom_vline(xintercept = tmaxpred_sum_gardenave - mean(c(ifg_tmax_ave, cso_tmax_ave)), linewidth = .2) +
  geom_vline(xintercept = 0, linetype = "longdash") +
  theme_bw() +
  labs(x = "Maximum temperature, transfer distance (°C)",
       y = "Relative growth rate 2014-2024")


##B


fig_rf24_2<-ggplot(gpc24_2, aes(x = tdiff2, y = mrf)) + 
  geom_line(linetype = "dashed", color = "#c5283d") + 
  geom_ribbon(aes(ymin = lower, ymax = upper), color = NA, fill = "#c5283d", alpha = .1) +
  geom_line(data = p_gamrf24_2, color = "#c5283d") + 
  geom_ribbon(data = p_gamrf24_2, aes(ymin = lower_ci, ymax = upper_ci), color = NA, fill = "#c5283d", alpha = 0.2)  +
  #scale_color_viridis(begin = .25, end = .75, name = "Annual maximum\ntemperature,\naveraged across\ngardens (°C)", option = "magma") +
  #scale_fill_viridis(begin = .25, end = .75, name = "Annual maximum\ntemperature,\naveraged across\ngardens (°C)", option = "magma") +
  scale_x_continuous(limits = c(-20, 9)) + 
  #cale_y_continuous(limits = c(0, 1)) +
  geom_vline(xintercept = -5.2, linewidth = .2) +
  geom_vline(xintercept = tmaxpred_sum_gardenave - mean(c(ifg_tmax_ave, cso_tmax_ave)), linewidth = .2) +
  geom_vline(xintercept = 0, linetype = "longdash") +
  theme_bw() +
  labs(x = "Maximum temperature, transfer distance (°C)",
       y = "Relative fitness, 2024")

cowplot::plot_grid(fig_rgr24_tave, fig_rf24_2, align = "hv", labels = c("A", "B")) 
#ggsave("~/git/demographyMS/figures/figS2.png", height = 3, width = 6, units = "in", type = "cairo")

# Fig S3: Averaged tmax RGR fig single-yr -------------------------------------------

#fig_rgr_ave<-
  ggplot(filter(p_gam1_temp2, Year != 2014, Year != 2018), aes(x = tdiffyr2, y = rgrh)) + 
  geom_line(aes(color = mtemp), linewidth = .8) + 
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, group = Year, fill = mtemp), alpha = .2)  +
  geom_line(data = gau_pred_r2, aes(color = mtemp), linetype = "dashed", linewidth = .8) + 
  geom_ribbon(data = gau_pred_r2, aes(ymin = lower, ymax = upper, group = Year, fill = mtemp), color = NA, alpha = .1) +
  scale_color_viridis(begin = .25, end = .75, name = "Annual maximum\ntemperature,\naveraged across\ngardens (°C)", option = "magma") +
  scale_fill_viridis(begin = .25, end = .75, name = "Annual maximum\ntemperature,\naveraged across\ngardens (°C)", option = "magma") +
  scale_x_continuous(limits = c(-20, 9)) + 
  coord_cartesian(ylim = c(0, .75)) +
  facet_wrap(~Year_bytemp) +
  geom_vline(xintercept = 0, linetype = "longdash") +
  geom_vline(xintercept = -5.2, linewidth = .2) +
  geom_vline(xintercept = tmaxpred_sum_gardenave - mean(c(ifg_tmax_ave, cso_tmax_ave)), linewidth = .2) +
  theme_bw() +
  guides(linetype = "none", linewidth = "none", alpha = "none") +
  labs(x = "Yearly temperature transfer distance (°C)",
       y = "Relative growth rate (single-year)",
       facet = "Year")

#ggsave("~/git/demographyMS/figures/figS3.png", height = 6, width = 8, units = "in", type = "cairo")

# Fig S4: temp correlation ------------------------------------------------

garden_tmin<-read.csv(paste0(dir, "03_spatial_prediction_figures/garden_clim_monthly.csv")) %>%
  separate(Date, c("Year", "Month"), sep = "-") %>%
  filter(Month %in% c("12", "01", "02")) %>%
  mutate(gYear = case_when(
    Month == "12" ~ as.character(as.numeric(Year) + 1),
    .default = Year
  )) %>%
  group_by(Site, gYear) %>%
  summarize(across(c(ppt, tmax, tmin), list(mean = mean, max = max), .names = "{.col}_{.fn}_decjanfeb")) %>%
  filter(gYear != "2014") %>%
  mutate(gYear = as.integer(gYear)) %>%
  ungroup() %>%
  summarize(bioclim_11 = mean(tmin_mean_decjanfeb), .by = Site) %>%
  mutate(tmax_sum = c(cso_tmax_ave, ifg_tmax_ave),
         Elevation..m. = c(70, 822))

ggplot(dat_clim_pc, aes(y = bioclim_11, x = tmax_sum, color = Elevation..m., fill = Elevation..m.)) + 
  geom_point() + 
  geom_point(data = garden_tmin, shape = 24, size = 3, color = "black") +
  scale_color_viridis() + 
  scale_fill_viridis() +
  labs(y = "Winter mean temp (°C)",
       x = "Summer max. temp (°C)",
       color = "Elevation (m)",
       fill = "Elevation (m)") 
#ggsave("~/git/demographyMS/figures/figS4.png",height = 4, width = 5, units = "in")

