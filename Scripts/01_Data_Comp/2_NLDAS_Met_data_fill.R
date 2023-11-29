#### Script for making met driver data from NDLAS data 
#DWH - September 2022

#load libraries 
library(tidyverse)
library(plotly)
library(LakeMetabolizer) #using to convert SW to PAR  
library(ggpmisc) #equations on ggplot
library(Metrics) #rmse

#### Read in and format data ####
#Read in data 
nldas <- read_csv("./Data/Carey2022_glm/FCR_GLM_NLDAS_010113_123119_GMTadjusted.csv")

nldas_filt <- nldas %>% 
  filter(time >= ymd_hms("2015-11-09 00:00:00"))

fcrmet <- read_csv("./Data/EDI2023/FCR_Met_final_2015_2022.csv")

#make met hourly
fcrmet_filt <- fcrmet %>% 
  select(DateTime, ShortwaveRadiationUp_Average_W_m2, WindSpeed_Average_m_s) %>% 
  filter(DateTime >= ymd_hms("2015-11-09 00:00:00"),
         DateTime < ymd_hms("2020-01-01 00:00:00")) %>% 
  mutate(Date = as.Date(DateTime),
         Hour = hour(DateTime)) %>% 
  group_by(Date, Hour) %>% 
  summarise(SW_hourly = mean(ShortwaveRadiationUp_Average_W_m2, na.rm = T),
            Wind_hourly = mean(WindSpeed_Average_m_s, na.rm = T),
            .groups = 'drop') %>% 
  mutate(Datetime = make_datetime(year = year(Date), month = month(Date), day = day(Date), 
                                  hour = Hour, tz ="UTC")) %>% 
  select(Datetime, SW_hourly, Wind_hourly)

joinedSW <- left_join(nldas_filt, fcrmet_filt, by = c("time" = "Datetime"))


#### Comparing NLDAS to FCR shortwave data ####

#compare hourly met to nldas 
SW_nldas_fcrmet<- joinedSW %>% 
  select(time, SW_hourly, ShortWave) %>% 
  rename(FCRmet = SW_hourly,
         NLDAS = ShortWave) %>% 
  ggplot(aes(x = NLDAS , y = FCRmet))+
  geom_point()+
  ggtitle("Shortwave FCR met ~ NLDAS")+
  labs(x = "NLDAS Shortwave", y = "FCR met station Shortwave")+
  stat_poly_line(method = "lm", linewidth = 2)+
  stat_poly_eq(formula=y~x, label.x = "left", label.y="top", parse=TRUE, inherit.aes = F,
               aes(x = NLDAS, y = FCRmet, label=paste(..adj.rr.label..,..p.value.label..,sep="~~~"),size=3))+
  stat_poly_eq(label.y = 0.85,
               aes(label = after_stat(eq.label))) +
  geom_abline(slope = 1, intercept = 0, linewidth=2, linetype =2, col = "red")+ #1:1 line
  theme_classic()

SW_nldas_fcrmet

dir.create("./Figures/")

ggsave(filename = "./Figures/Fig_S2_Shortwave_NLDAS_FCRmet.png",
       SW_nldas_fcrmet, device = "png", width = 180, height = 150, units = "mm")

#### Comparing NLDAS to FCR WindSpeed data ####

Wind_nldas_fcrmet<- joinedSW %>% 
  select(time, Wind_hourly, WindSpeed) %>% 
  rename(FCRmet = Wind_hourly,
         NLDAS = WindSpeed) %>% 
  ggplot(aes(x = NLDAS , y = FCRmet))+
  geom_point()+
  ggtitle("WindSpeed FCR met ~ NLDAS")+
  labs(x = "NLDAS WindSpeed", y = "FCR met station WindSpeed")+
  stat_poly_line(method = "lm", linewidth = 2)+
  stat_poly_eq(formula=y~x, label.x = "left", label.y="top", parse=TRUE, inherit.aes = F,
               aes(x = NLDAS, y = FCRmet, label=paste(..adj.rr.label..,..p.value.label..,sep="~~~"),size=3))+
  stat_poly_eq(label.y = 0.8,
               aes(label = after_stat(eq.label))) +
  geom_abline(slope = 1, intercept = 0, size=2, linetype =2, col = "red")+ #1:1 line
  theme_classic()

Wind_nldas_fcrmet

ggsave(filename = "./Figures/Fig_S3_WindSpeed_NLDAS_FCRmet.png",
       Wind_nldas_fcrmet, device = "png", width = 180, height = 150, units = "mm")



#### Make input files for SW and WindSpeed based on NLDAS ####

sw_lm <- summary(lm(joinedSW$SW_hourly~joinedSW$ShortWave))
wind_lm <- summary(lm(joinedSW$Wind_hourly~joinedSW$WindSpeed))


model_input <- joinedSW %>% 
  select(time, SW_hourly, Wind_hourly, ShortWave, WindSpeed) %>% 
  mutate(SW_regress = ( coefficients(sw_lm)[1] + (coefficients(sw_lm)[2]* ShortWave) ) ,
         WindSpeed_regress = ( coefficients(wind_lm)[1] + (coefficients(wind_lm)[2] * WindSpeed) )) %>% 
  mutate(PAR_from_SWregress = sw.to.par.base(SW_regress, coeff = 2.114)) %>% 
  rename(SW_hourly_metstation = SW_hourly,
         Wind_hourly_metstation = Wind_hourly,
         ShortWave_nldas = ShortWave,
         WindSpeed_nldas = WindSpeed)

#write csv
write.csv(model_input, "./Data/Generated_Data/Met_input_from_NLDAS.csv", row.names = F)


#check rmse and model fit for SW and Wind
model_input_tests <- model_input %>% 
  mutate(Diff_SW = SW_hourly_metstation - ShortWave_nldas,
         Diff_Wind = Wind_hourly_metstation - WindSpeed_nldas) %>% 
  filter(!is.na(SW_hourly_metstation),
         !is.na(Wind_hourly_metstation))
  
plot(model_input_tests$time, model_input_tests$Diff_SW)

round(rmse(model_input_tests$SW_hourly_metstation, model_input_tests$ShortWave_nldas), digits = 3)


plot(model_input_tests$time, model_input_tests$Diff_Wind)

round(rmse(model_input_tests$Wind_hourly_metstation, model_input_tests$WindSpeed_nldas), digits = 3)




