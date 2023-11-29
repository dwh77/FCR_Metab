### Compile driver data for Metab envi driver analysis 


### load packages needed ####
library(tidyverse) #for data filtering and compiling
library(LakeMetabolizer) #using to convert Shortwave to PAR for metabolism model input 
library(ggpmisc)
library(Metrics)


### Compile EXO and thermistor data ####

#read in data from EDI
catwalk_EDI <- read_csv("./Data/EDI2023/FCR_Catwalk_2018_2022.csv", col_types = cols(.default = "d", Reservoir = "c", DateTime = "T"))

#Exo daily 
exo_daily <- catwalk_EDI %>% 
  mutate(Date = as.Date(DateTime)) %>% 
  group_by(Date) %>%
  dplyr::summarize(daily_EXO_fdom = mean(EXOfDOM_QSU_1, na.rm = TRUE),
            daily_EXO_chla_ugL = mean(EXOChla_ugL_1, na.rm = TRUE)) %>% 
  filter(Date > ymd("2018-08-28"),
         Date < ymd("2022-03-02")) 


### compile met variables and summarize to 10 or 15 for exo or wvwa #### 

#read in data from EDI
met_EDI <- read_csv("./Data/EDI2023/FCR_Met_final_2015_2022.csv")

#daily met 
met_daily <- met_EDI %>% 
  mutate(Date = as.Date(DateTime)) %>% 
  group_by(Date) %>%
  dplyr::summarize(daily_airtemp = mean(AirTemp_C_Average, na.rm = T),
            daily_rain_mm = sum(Rain_Total_mm, na.rm = TRUE)) %>% 
  filter(Date > ymd("2015-11-10"),
         Date < ymd("2022-03-02")) 


#### MET NLDAS ####
nldas <- read_csv("./Data/Carey2022_glm/FCR_GLM_NLDAS_010113_123119_GMTadjusted.csv")

nldas_daily <- nldas %>% 
  mutate(Date = as.Date(time)) %>% 
  group_by(Date) %>%
  dplyr::summarize(daily_airtemp_nldas = mean(AirTemp, na.rm = T),
            daily_rain_mm_nldas = sum(Rain, na.rm = TRUE)) %>% #nldas precip is in units of kg/m2 but this equals mm
  filter(Date > ymd("2015-11-10"),
         Date < ymd("2022-03-02"))

nldas_met_daily_offset <- left_join(nldas_daily, met_daily, by = "Date")

airtemp <- nldas_met_daily_offset %>% 
  ggplot(aes(x = daily_airtemp_nldas , y = daily_airtemp))+
  geom_point()+
  ggtitle("AirTemp FCR met ~ NLDAS")+
  # ylim(0,1000)+
  # xlim(0,1000)+
  labs(x = "NLDAS Air Temp", y = "FCR met station AirTemp")+
  stat_poly_line(method = "lm", linewidth = 2)+
  # geom_smooth(aes(x = NLDAS, y = FCRmet), method = "lm", inherit.aes = F)+
  stat_poly_eq(formula=y~x, label.x = "left", label.y="top", parse=TRUE, inherit.aes = F,
               aes(x = daily_airtemp_nldas, y = daily_airtemp, label=paste(..adj.rr.label..,..p.value.label..,sep="~~~"),size=3))+
  stat_poly_eq(label.y = 0.8,
               aes(label = after_stat(eq.label))) +
  geom_abline(slope = 1, intercept = 0, linewidth=2, linetype =2, col = "red")+ #1:1 line
  theme_classic()

airtemp


rain <- nldas_met_daily_offset %>% 
  ggplot(aes(x = daily_rain_mm_nldas , y = daily_rain_mm))+
  geom_point()+
  ggtitle("Rain FCR met ~ NLDAS")+
  # ylim(0,1000)+
  # xlim(0,1000)+
  labs(x = "NLDAS Rain", y = "FCR met station Rain")+
  stat_poly_line(method = "lm", linewidth = 2)+
  # geom_smooth(aes(x = NLDAS, y = FCRmet), method = "lm", inherit.aes = F)+
  stat_poly_eq(formula=y~x, label.x = "left", label.y="top", parse=TRUE, inherit.aes = F,
               aes(x = daily_rain_mm_nldas, y = daily_rain_mm, label=paste(..adj.rr.label..,..p.value.label..,sep="~~~"),size=3))+
  stat_poly_eq(label.y = 0.8,
               aes(label = after_stat(eq.label))) +
  geom_abline(slope = 1, intercept = 0, linewidth=2, linetype =2, col = "red")+ #1:1 line
  theme_classic()

rain


nldasmet_regress <- nldas_met_daily_offset %>% 
  mutate(AirTemp_regress = (0.488 + 0.94 * daily_airtemp_nldas),
         Rain_regress = (0.412 + 33.5 * daily_rain_mm_nldas)
         ) %>% 
  mutate(AirTemp_fin = ifelse(is.na(daily_airtemp), AirTemp_regress, daily_airtemp),
         Rain_fin = ifelse(is.na(daily_rain_mm), Rain_regress, daily_rain_mm)) %>% 
  select(Date, AirTemp_fin, Rain_fin)

met_2020_2022 <- met_daily %>% 
  filter(Date >= ymd("2020-01-01")) %>% 
  rename(AirTemp_fin = daily_airtemp,
         Rain_fin = daily_rain_mm)

MET_FIN <- rbind(nldasmet_regress, met_2020_2022)

#### Hydrology ####
inflow <- read_csv("./Data/EDI2023/FCR_weir_2022.csv")

inf_daily <- inflow %>% 
  select(DateTime, WVWA_Flow_cms, VT_Flow_cms) %>% 
  filter(DateTime > ymd_hms("2015-11-10 00:00:00"),
         DateTime < ymd("2022-03-02")) %>%         
  mutate(Date = as.Date(DateTime)) %>% 
  group_by(Date) %>% 
  dplyr::summarise(WVWA_Flow_cms_daily_mean = mean(WVWA_Flow_cms, na.rm = TRUE),
            VT_Flow_cms_daily_mean = mean(VT_Flow_cms, na.rm = TRUE)) %>% 
  mutate(WRT_days_daily_wvwa = ( (3.1E5/WVWA_Flow_cms_daily_mean) * (1/60) * (1/60) * (1/24) ),
         WRT_days_daily_vt = ( (3.1E5/VT_Flow_cms_daily_mean) * (1/60) * (1/60) * (1/24) )
         ) 


inf_daily_fin <- inf_daily %>% 
  mutate(Fin_WRT_days = ifelse(is.na(WRT_days_daily_wvwa), WRT_days_daily_vt, WRT_days_daily_wvwa),
         Fin_flow_cms = ifelse(is.na(WVWA_Flow_cms_daily_mean), VT_Flow_cms_daily_mean , WVWA_Flow_cms_daily_mean)
         ) %>% 
  select(Date, Fin_flow_cms) 


#WRT stats 
inf_wrt <- inf_daily %>% 
  mutate(Fin_WRT_days = ifelse(is.na(WRT_days_daily_wvwa), WRT_days_daily_vt, WRT_days_daily_wvwa),
         Fin_flow_cms = ifelse(is.na(WVWA_Flow_cms_daily_mean), VT_Flow_cms_daily_mean , WVWA_Flow_cms_daily_mean)
  )%>% 
  #filter(Date > ymd("2019-01-01")) %>% 
  select(Date, Fin_flow_cms, Fin_WRT_days) 

plot(inf_wrt$Date, inf_wrt$Fin_WRT_days)
summary(inf_wrt$Fin_WRT_days)
sd(inf_wrt$Fin_WRT_days, na.rm = T)
hist(inf_wrt$Fin_WRT_days, na.rm = T)



#### Chemsitry #### 
chem <- read_csv("./Data/EDI2023/FCR_chem_2022.csv")

chemA <- chem %>% 
  filter(Reservoir == "FCR",
         Site == 50,
         Depth_m %in% c(0.8, 1.6)) %>% 
  mutate(Date = as.Date(DateTime)) %>% 
  filter(Date > ymd("2015-11-10")) %>% 
  select(Date, TN_ugL, TP_ugL, DOC_mgL, SRP_ugL, NO3NO2_ugL, NH4_ugL)


chemB <- chemA %>% 
  group_by(Date) %>% 
  summarise(across(.cols = where(is.numeric), .fns = ~mean(.x, na.rm = TRUE)), .groups = "drop") 


#### Secchi - Kd #### 
secchi <- read_csv("./Data/EDI2023/FCR_secchi_2022.csv")

secchiA <- secchi %>% 
  filter(Reservoir == "FCR",
         Site == 50) %>% 
  mutate(Date = as.Date(DateTime)) %>% 
  select(Date, Secchi_m) %>% 
  mutate(Kd = 1.7 / Secchi_m) %>% 
  filter(Date > ymd("2015-11-09"),
         Date < ymd("2022-03-01"))
  

#### Kd for PAR and CTD #### 
## CTD from Abby - https://github.com/CareyLabVT/Reservoirs/blob/master/Data/DataNotYetUploadedToEDI/Raw_CTD/CTD_code/Light%20attenuation.Rmd
ctd <- read.csv("./Data/EDI2023/CTD_final_2013_2022.csv") #load saved data

ctdKD <- ctd %>% 
  filter(Reservoir=="FCR")%>%
  filter(Site == 50)%>%
  mutate(Date = ymd_hms(DateTime)) %>% 
  filter(!is.na(PAR_umolm2s))%>%
  group_by(Date)%>%
  filter(sum(Depth_m<0)>0)%>% #Filter so there is at least one measurement out of the water
  mutate(I0 = mean(PAR_umolm2s[Depth_m<0], na.rm = T))%>% #I0 is the mean of surface measures
  filter(Depth_m>0)%>%
  dplyr::summarize(I0 = unique(I0),
            k = coef(lm(I(log(PAR_umolm2s)-log(I0))~ 0 + Depth_m)), #From light attenuation eq.
            r2 = summary(lm(I(log(PAR_umolm2s)-log(I0))~ 0 + Depth_m))$r.squared, #R2
            Zeu = min(Depth_m[PAR_umolm2s<I0/100]), #Euphotic zone = 1% of surface light
            Zeu_0.1 = min(Depth_m[PAR_umolm2s<I0/1000]))%>% #Or 0.1% of surface light
  filter(r2>0.9)


par <- read.csv("./Data/EDI2023/YSI_PAR_profiles_2013-2022.csv") %>% 
  filter(Reservoir=="FCR")%>%
  filter(Site == 50)%>%
  mutate(Date = ymd_hms(DateTime)) 

parKD <- par %>% 
  select(Date, Depth_m, PAR_umolm2s) %>% 
  filter(!is.na(PAR_umolm2s)) %>% 
  filter(PAR_umolm2s > 0) %>% #removing 0 PAR values since we can't compute w/ log functions below
  group_by(Date)%>%
  mutate(I0 = mean(PAR_umolm2s[Depth_m<0], na.rm = T))%>% #I0 is the mean of surface measures
  filter(Depth_m>0,
         !is.na(I0))%>% #removing dates that didn't have an above surface PAR since that breaks code below
  dplyr::summarize(I0 = unique(I0),
            k = coef(lm(I(log(PAR_umolm2s)-log(I0))~ 0 + Depth_m)), #From light attenuation eq.
            r2 = summary(lm(I(log(PAR_umolm2s)-log(I0))~ 0 + Depth_m))$r.squared, #R2
            Zeu = min(Depth_m[PAR_umolm2s<I0/100]), #Euphotic zone = 1% of surface light
            Zeu_0.1 = min(Depth_m[PAR_umolm2s<I0/1000]))%>% #Or 0.1% of surface light
  filter(r2>0.9)
  

#### Put KD data together ####
secchi_Kd_fin <- secchiA %>% rename(Kd_secchi = Kd)

ctdKD_fin <- ctdKD %>% 
  mutate(k = -k,
         Date = as.Date(Date)) %>%
  group_by(Date) %>% 
  summarise(Kd_ctd = mean(k))
  
parKD_fin <- parKD %>% 
  mutate(k = -k,
         Date = as.Date(Date)) %>%
  group_by(Date) %>% 
  summarise(Kd_par = mean(k))

kd_joinA <- full_join(secchi_Kd_fin, parKD_fin, by = "Date")
kd_joinB <- full_join(kd_joinA, ctdKD_fin, by = "Date")

kd_joinC <- kd_joinB %>% 
  mutate(Kd_PARorCTD = ifelse(!is.na(Kd_ctd), Kd_ctd, Kd_par))

kd_joinD <- kd_joinC %>% 
  mutate(Kd_fin = ifelse(is.na(Kd_PARorCTD), Kd_secchi, Kd_PARorCTD)) %>% 
  arrange(Date)

kd_FIN <- kd_joinD %>% 
  select(Date, Kd_fin) %>% 
  rename(Kd = Kd_fin) %>% 
  filter(Date > ymd("2015-11-09"),
         Date < ymd("2022-03-01"))

#look at rmse 
kd_rmse <- kd_joinC %>% 
  filter(!is.na(Kd_PARorCTD),
         !is.na(Kd_secchi)) %>% 
  mutate(Diff = Kd_PARorCTD - Kd_secchi)

plot(kd_rmse$Diff)
median(kd_rmse$Diff)
median(kd_rmse$Kd_PARorCTD)
median(kd_rmse$Kd_secchi)

round(rmse(kd_rmse$Kd_PARorCTD, kd_rmse$Kd_secchi), digits = 3)

cor(kd_rmse$Kd_PARorCTD, kd_rmse$Kd_secchi, method = "pearson")


#### E24 #### 
source('./Scripts/01_Data_Comp/E24_functions.R')

E24 <- kd_FIN

E0datetimes <- seq(ISOdate(2015, 11, 09,0,0,0,tz="UTC"),
                   ISOdate(2022,03,01,0,0,0,tz="UTC"), by="10 min")

PARvector <- incident(date = E0datetimes, latitude = 37.30, longitude = -79.84, elevation = (1663 * 0.3048),
                      timezone = -6, reflectance = TRUE)

E0_df <- data.frame(E0datetimes, PARvector)

E0_df_daily <- E0_df %>% 
  mutate(Date = as.Date(E0datetimes)) %>% 
  group_by(Date) %>% 
  summarise(E0 = mean(PARvector))

# get top of metalim 
temp_prof_15_18 <- read_csv("./Data/Model_Input/2015_18/FCR_2015_18_TempProfiles_hobos_and_interped2018.csv")

temp_prof_18_22 <- read_csv("./Data/Model_Input/2018_22/FCR_2018_22_TempProfiles.csv") %>% 
  dplyr::filter(dateTime > ymd_hms("2019-01-02 00:00:00"))

profiles <- bind_rows(temp_prof_15_18, temp_prof_18_22)

profiles <- profiles %>% 
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
         wtr_09.3 = temp9.3)

metalayers <- rLakeAnalyzer::ts.meta.depths(profiles, na.rm = T)

top_meta <- metalayers %>% 
  mutate(Date = as.Date(datetime)) %>% 
  group_by(Date) %>% 
  summarise(top_meta = mean(top, na.rm = T))


E24_fin <- left_join(E24, top_meta, by = "Date")
E24_fin <- left_join(E24_fin, E0_df_daily, by = "Date")

E24_fin <- E24_fin %>% 
  mutate(E24_umolpm2s = E0*(1-exp(-1*Kd*top_meta))*((Kd*top_meta)^-1))

E24_tojoin <- E24_fin %>% 
  select(Date, E24_umolpm2s) %>% 
  filter(!is.na(E24_umolpm2s))


#### schmidt ####
schmidt <- read_csv("./Data/Generated_Data/schmidt_timeseries.csv")

#### filtered chla ####
# filtchla <- read_csv("./Data/EDI2023/manual_chlorophyll_2014_2022.csv")
# 
# filtchla_A <- filtchla %>% 
#   filter(Reservoir == "FCR",
#          Site == 50,
#          Depth_m == 1.6) %>% 
#   mutate(Date = as.Date(DateTime)) %>% 
#   select(Date, Chla_ugL)

#### Combine ####
# head(exo_daily) 
# head(MET_FIN) #met that corrects for gaps w/ NDLAS interp 
# head(inf_daily_fin) #mean daily inflow in cms
# head(chemB)
# head(kd_FIN)
# head(E24_tojoin)
# head(schmidt)
####head(filtchla_A)

## join data together 
joinA <- full_join(MET_FIN, exo_daily, by = "Date")
joinB <- full_join(joinA, inf_daily_fin, by = "Date")
joinC <- full_join(joinB, chemB, by = "Date")
joinD <- full_join(joinC, kd_FIN, by = "Date")
joinE <- full_join(joinD, E24_tojoin, by = "Date")
join_fin <- full_join(joinE, schmidt, by = "Date")
# join_fin <- full_join(join_F, filtchla_A, by = "Date")

join_fin <- join_fin[order(join_fin$Date),]

join_fin <- join_fin %>% 
  filter(Date <= ymd("2022-03-01")) %>% 
  mutate(Rain_lag1 = lag(Rain_fin, 1),
         Rain_lag2 = lag(Rain_fin, 2))

#write.csv(join_fin, "./Data/EnviDrivers_daily_compiled_15_22_22nov23.csv", row.names = F)

drivers_ZT <- join_fin %>% 
  mutate(airtemp_ZT = scale(AirTemp_fin, center = TRUE, scale = TRUE),
         rain_ZT = scale(Rain_fin, center = TRUE, scale = TRUE),
         rain_lag1_ZT = scale(Rain_lag1, center = TRUE, scale = TRUE),
         rain_lag2_ZT = scale(Rain_lag2, center = TRUE, scale = TRUE),
         EXO_fdom_ZT = scale(daily_EXO_fdom, center = TRUE, scale = TRUE),
         EXO_chla_ZT = scale(daily_EXO_chla_ugL, center = TRUE, scale = TRUE),
         flow_ZT = scale(Fin_flow_cms, center = TRUE, scale = TRUE),
         TN_ZT = scale(TN_ugL, center = TRUE, scale = TRUE),
         TP_ZT = scale(TP_ugL, center = TRUE, scale = TRUE),
         DOC_ZT = scale(DOC_mgL, center = TRUE, scale = TRUE),
         SRP_ZT = scale(SRP_ugL, center = TRUE, scale = TRUE),
         NO3NO2_ZT = scale(NO3NO2_ugL, center = TRUE, scale = TRUE),
         NH4_ZT = scale(NH4_ugL, center = TRUE, scale = TRUE),
         Kd_ZT = scale(Kd, center = TRUE, scale = TRUE),
         E24_ZT = scale(E24_umolpm2s, center = TRUE, scale = TRUE),
         schmidt_ZT = scale(schmidt_daily, center = TRUE, scale = TRUE)
         #filt_chla_ZT = scale(Chla_ugL, center = TRUE, scale = TRUE)
         ) 

write.csv(drivers_ZT, "./Data/Generated_Data/EnviDrivers_daily_ZT.csv", row.names = F)


