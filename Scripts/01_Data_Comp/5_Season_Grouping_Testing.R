## Testing different stratification cutoffs for season grouping 
## DWH - summer 2022

#load libraries 
library(tidyverse)
library(plotly)
library(rLakeAnalyzer)

### read in bathy data for schmidt stability calculations 
edi_bathy <- read_csv("./Data/EDI2023/Reservoir_bathy_EDI.csv")

bathy_schmidt <- edi_bathy %>% 
  filter(Reservoir == 'FCR') %>% 
  dplyr::select(Depth_m, SA_m2) %>% 
  dplyr::rename(depths = 1,
         areas = 2)


### read in Catwalk and hobo thermistors, and interpolated temp profiles for early 2018
hobo_tempprofs <- read_csv("./Data/Model_Input/2015_18/FCR_2015_18_TempProfiles_hobos.csv") 

interp2018_profiles <- read_csv("./Data/Model_Input/2015_18/FCR_2018_TempProfiles_hobo_glm_ctdysi_scc.csv") %>% 
  filter(dateTime > ymd_hms("2018-01-15 00:00:00"),
         dateTime < ymd_hms("2018-08-29 00:00:00"))

catwalk_tempprofs <- read_csv("./Data/Model_Input/2018_22/FCR_2018_22_TempProfiles.csv") 

tempprofs <- bind_rows(hobo_tempprofs, interp2018_profiles, catwalk_tempprofs) 


### Calculate density and temperature differences for profiles and assign startifciation definitions
strat_calcs <- tempprofs %>% 
  # mutate(Diff_0.1_9 = temp0.1 - temp9.0,
  #        Strat_1C_surfto9 = ifelse(Diff_0.1_9 < 1, "Not_Strat", "Strat")) %>% 
  # mutate(Diff_Density_surfto9 = water.density(temp9.0) - water.density(temp0.1),
  #        Strat_Dens0.3_surfto9 = ifelse(Diff_Density_surfto9 > 0.3, "Strat", "Not_Strat")) %>% 
  mutate(Diff_C_1_8 = temp1.0 - temp8.0,
         Diff_Dens_1to8 = water.density(temp8.0) - water.density(temp1.0),
         Strat_1C_1to8 = ifelse(Diff_C_1_8 < 1, "Not_Strat", "Strat")) %>%
  mutate(Strat_Dens0.05_1to8 = ifelse(Diff_Dens_1to8 > 0.05, "Strat", "Not_Strat")) %>% 
  mutate(Strat_Dens0.1_1to8 = ifelse(Diff_Dens_1to8 > 0.1, "Strat", "Not_Strat")) %>% 
  mutate(Strat_Dens0.5_1to8 = ifelse(Diff_Dens_1to8 > 0.5, "Strat", "Not_Strat"))

### Calculate schmidt stability 
schmidt_data <- strat_calcs %>% 
  rename(datetime = dateTime,
       wtr_00.1 = temp0.1,
       wtr_00.8 = temp0.8,
       wtr_01.0 = temp1.0,
       wtr_01.6 = temp1.6,
       wtr_02.0 = temp2.0,
       wtr_02.8 = temp2.8,
       wtr_03.0 = temp3.0,
       wtr_03.8 = temp3.8,
       wtr_04.0 = temp4.0,
       wtr_05.0 = temp5.0,
       wtr_06.0 = temp6.0,
       wtr_06.2 = temp6.2,
       wtr_07.0 = temp7.0,
       wtr_08.0 = temp8.0,
       wtr_09.0 = temp9.0,
       wtr_09.3 = temp9.3) %>% 
  dplyr::select(datetime, wtr_00.1, wtr_00.8, wtr_01.0, wtr_01.6, wtr_02.0, wtr_02.8, wtr_03.0, wtr_03.8, wtr_04.0,
         wtr_05.0, wtr_06.0, wtr_06.2, wtr_07.0, wtr_08.0, wtr_09.0, wtr_09.3)

schmidt <- ts.schmidt.stability(schmidt_data, bathy_schmidt, na.rm = T)
#plot(schmidt$datetime, schmidt$schmidt.stability)

# buoy <- ts.buoyancy.freq(schmidt_data, na.rm = F) #NAs giving shorter TS and na.rm = T gives error
# plot(buoy$datetime, buoy$n2)

### join Schmidt stability to previous stratification metrics 
strat_calcsA <- left_join(strat_calcs, schmidt, by = c("dateTime" = "datetime"))
#strat_calcsA <- left_join(exo_strat_calcsA, buoy, by = c("dateTime" = "datetime"))


### calculate daily values of schmidt and density and temp differences 
#then look at rate of change of schmidt, and see how all these variables compare 
daily_ROC <- strat_calcsA %>% 
  mutate(Date = as.Date(dateTime)) %>% 
  group_by(Date) %>% 
  summarise(schmidt_daily = median(schmidt.stability, na.rm = T),
            # buoy_daily = median(n2, na.rm = T),
            Diff_Dens_1to8 = median(Diff_Dens_1to8, na.rm = T),
            Diff_C_1_8 = median(Diff_C_1_8, na.rm = T)
  ) %>% 
  mutate(ROC_schmidt = ifelse(lag(is.na(schmidt_daily)), NA, schmidt_daily - lag(schmidt_daily))
         # ROC_buoy = ifelse(lag(is.na(buoy_daily)), NA, buoy_daily - lag(buoy_daily))
  ) %>% 
  mutate(Schmidt_percent_max = (schmidt_daily / max(schmidt_daily, na.rm = T) ) * 100
  ) %>% 
  mutate(ROC_schmidt = ifelse(ROC_schmidt > 40, NA, ROC_schmidt)) #remove ROC when 2018 gap happens 

#plot(daily_ROC$Date, daily_ROC$schmidt_daily)
#plot(daily_ROC$Date, daily_ROC$Schmidt_percent_max)
#plot(daily_ROC$Date, daily_ROC$ROC_schmidt)


daily_ROC_long <- daily_ROC %>% 
  pivot_longer(-1)

roc_plots <- daily_ROC_long %>% 
  mutate(new_name = ifelse(name == "Diff_Dens_1to8", "A - Density Diff 1 - 8m ", NA),
         new_name = ifelse(name == "Diff_C_1_8", "B - Temp Diff 1 - 8m ", new_name),
         new_name = ifelse(name == "schmidt_daily", "C - schmidt stability", new_name),
         new_name = ifelse(name == "Schmidt_percent_max", "D - Schmidt % max ", new_name),
         new_name = ifelse(name == "ROC_schmidt", "E - ROC schmidt stability", new_name),
  ) %>%
  ggplot(aes(x = Date, y = value))+
  geom_point()+
  facet_wrap(~new_name, scales = "free_y", nrow = 6)+
  theme_classic()

#roc_plots
#ggplotly(roc_plots)
# ggsave(filename = "./Figures/Schmidt_0.1dens_ROC_figure.png",
#        roc_plots, device = "png", width = 300, height = 300, units = "mm")


dens_diff <- daily_ROC_long %>% 
  filter(name == "Diff_Dens_1to8") %>% 
  ggplot(aes(x = Date, y = value))+
  geom_point()+
  labs(y = "Density Diff")+
  geom_hline(yintercept = 0.5)+
  geom_hline(yintercept = 0.3)+
  geom_hline(yintercept = 0.1)+
  geom_hline(yintercept = 0.05)+
  theme_classic()

#dens_diff
#ggplotly(dens_diff)
# ggsave(filename = "./Figures/DensityDiff_1to8m.png",
#        dens_diff, device = "png", width = 300, height = 300, units = "mm")

schmidt_per_plot <- daily_ROC_long %>% 
  filter(name == "Schmidt_percent_max") %>% 
  ggplot(aes(x = Date, y = value))+
  geom_point()+
  labs(y = "Schmidt % of max")+
  geom_hline(yintercept = 2)+
  geom_hline(yintercept = 1)+
  geom_hline(yintercept = 0.1)+
  theme_classic()

schmidt_per_plot
# ggplotly(schmidt_per_plot)


### Make data frame that can be joined to metab data to determine seasons 
head(daily_ROC)

strat_metrics <- daily_ROC %>% 
  select(Date, Diff_C_1_8, Diff_Dens_1to8, Schmidt_percent_max) %>% 
  dplyr::mutate(
    Strat_Dens_0.05 = ifelse(Diff_Dens_1to8 > 0.05, "Strat", "Not_Strat"),
    Strat_Dens_0.1 = ifelse(Diff_Dens_1to8 > 0.1, "Strat", "Not_Strat"),
    Strat_Dens_0.3 = ifelse(Diff_Dens_1to8 > 0.3, "Strat", "Not_Strat"),
    Strat_Dens_0.5 = ifelse(Diff_Dens_1to8 > 0.5, "Strat", "Not_Strat"),
    Strat_C_1 = ifelse(Diff_C_1_8 > 0.5, "Strat", "Not_Strat"),
    Strat_schmidtmax_2 = ifelse(Schmidt_percent_max > 2, "Strat", "Not_Strat"),
    Strat_schmidtmax_5 = ifelse(Schmidt_percent_max > 5, "Strat", "Not_Strat")
    )

write.csv(strat_metrics, "./Data/strat_metrics_9nov23.csv", row.names = F)


### Plot different strat metrics 
strat_metrics %>% 
  pivot_longer(-c(1:4)) %>% 
  dplyr::mutate(value_numeric = ifelse(value == "Strat", 1, 0)) %>% 
  filter(name %in% c("Strat_C_1", "Strat_Dens_0.1", "Strat_schmidtmax_2")) %>% 
  ggplot()+
  geom_point(aes(x = Date, y = value_numeric, color = name))



### Put schmidt stability toghether for envi drivers 
head(daily_ROC)

schmidt_ts <- daily_ROC %>% 
  select(Date, schmidt_daily)

write.csv(schmidt_ts, "./Data/schmidt_timeseries_9nov23.csv", row.names = F)












