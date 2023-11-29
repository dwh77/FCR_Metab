### Compile FCR EXO and met data from metabolism calculations 


### load packages needed ####
library(tidyverse) #for data filtering and compiling
library(LakeMetabolizer) #using to convert Shortwave to PAR for metabolism model input 
library(ggpmisc) # for lm for par comparison
library(Metrics) #for rmse
library(zoo) #for na.approx 

dir.create('Data/Model_Input/')
dir.create('Data/Model_Input/2018_22')
dir.create('Data/Model_Input/2015_18')

### Compile EXO and thermistor data ####

#read in data from EDI
catwalk_EDI <- read_csv("./Data/EDI2023/FCR_Catwalk_2018_2022.csv")


##Geting EXO DO and Sensor temp for each year and writing csvs to use in metab model

#filter catwalk data to desired timeframe and average to 10 minute data (theres a short periods with 1 minute data that mess up metab model)
catwalk <- catwalk_EDI %>% 
  dplyr::filter(DateTime >= ymd_hms("2018-07-04 00:00:00"),
         DateTime <= ymd_hms("2022-03-01 23:50:00")) %>%  #getting rid of data before consistent DO mg/L was deployed. EXO was deployed 2018-08-07 
  dplyr::select(DateTime, EXODO_mgL_1, EXOTemp_C_1, ThermistorTemp_C_surface, ThermistorTemp_C_1, ThermistorTemp_C_2, ThermistorTemp_C_3, 
         ThermistorTemp_C_4, ThermistorTemp_C_5, ThermistorTemp_C_6, ThermistorTemp_C_7, 
         ThermistorTemp_C_8, ThermistorTemp_C_9) %>% 
  dplyr::mutate(Date = as.Date(DateTime)) %>% 
  dplyr::mutate(Hour = hour(DateTime),
         minutes = minute(DateTime),
         minutes = as.numeric(minutes)) %>% 
  dplyr::mutate(Minute_tenminutebreaks = ifelse(minutes >= 56 | minutes <= 05,  00, minutes),
         Minute_tenminutebreaks = ifelse(minutes >= 06 & minutes <= 15,  10, Minute_tenminutebreaks),
         Minute_tenminutebreaks = ifelse(minutes >= 16 & minutes <= 25,  20, Minute_tenminutebreaks),
         Minute_tenminutebreaks = ifelse(minutes >= 26 & minutes <= 35,  30, Minute_tenminutebreaks),
         Minute_tenminutebreaks = ifelse(minutes >= 36 & minutes <= 45,  40, Minute_tenminutebreaks),
         Minute_tenminutebreaks = ifelse(minutes >= 46 & minutes <= 55,  50, Minute_tenminutebreaks),
         Hour = ifelse(minutes >= 56, Hour+1, Hour)) %>% #This is so that a time like 12:56 doesn't get classified as 12:00 when it should be 13:00
  dplyr::rename(Minute = Minute_tenminutebreaks) %>% 
  dplyr::group_by(Date, Hour, Minute) %>% 
  dplyr::summarize(DO = mean(EXODO_mgL_1, na.rm = T),
            sensorTemp = mean(EXOTemp_C_1, na.rm = T),
            temp0.1 = mean(ThermistorTemp_C_surface, na.rm  = T),
            temp1.0 = mean(ThermistorTemp_C_1, na.rm  = T),
            temp2.0 = mean(ThermistorTemp_C_2, na.rm  = T),
            temp3.0 = mean(ThermistorTemp_C_3, na.rm  = T),
            temp4.0 = mean(ThermistorTemp_C_4, na.rm  = T),
            temp5.0 = mean(ThermistorTemp_C_5, na.rm  = T),
            temp6.0 = mean(ThermistorTemp_C_6, na.rm  = T),
            temp7.0 = mean(ThermistorTemp_C_7, na.rm  = T),
            temp8.0 = mean(ThermistorTemp_C_8, na.rm  = T),
            temp9.0 = mean(ThermistorTemp_C_9, na.rm  = T),
            .groups = "drop") %>% 
  dplyr::mutate(dateTime = make_datetime(year = year(Date), month = month(Date), day = day(Date), hour = Hour, min = Minute, tz ="UTC")) %>% 
  dplyr::select(dateTime, DO, sensorTemp, temp0.1, temp1.0, temp2.0, temp3.0, temp4.0, temp5.0, temp6.0, temp7.0, temp8.0, temp9.0)

#check that there's no duplicated data
catwalk$dateTime[duplicated(catwalk$dateTime)]


##DO and Temp data 2018 to 2022
do_18_22 <- catwalk %>%
  dplyr::filter(dateTime >= ymd_hms("2018-08-29 00:00:00"),
         dateTime <= ymd_hms("2022-03-01 23:50:00")) %>%
  dplyr::select(dateTime, DO) 
write.csv(do_18_22, file = "./Data/Model_Input/2018_22/FCR_2018_22_DO.csv", row.names = FALSE)


temp_18_22 <- catwalk %>%
  dplyr::filter(dateTime >= ymd_hms("2018-08-29 00:00:00"),
         dateTime <= ymd_hms("2022-03-01 23:50:00")) %>%
  dplyr::select(dateTime, sensorTemp) 
write.csv(temp_18_22, file = "./Data/Model_Input/2018_22/FCR_2018_22_sensorTemp.csv", row.names = FALSE)


#Temp profiles 2018 to 2022
temp_profiles_2018_22 <- catwalk %>%
  dplyr::select(-DO, -sensorTemp) %>% 
  dplyr::filter(dateTime >= ymd_hms("2018-01-01 00:00:00"),
         dateTime <= ymd_hms("2022-03-01 23:50:00"))
write.csv(temp_profiles_2018_22, file = "./Data/Model_Input/2018_22/FCR_2018_22_TempProfiles.csv", row.names = FALSE)





### compile met variables and summarize to 10 or 15 for exo or wvwa #### 

#read in data from EDI
met_EDI <- read_csv("./Data/EDI2023/FCR_Met_final_2015_2022.csv")

##first get 10 min par and windspeed for EXO data 

#filter to desired variables and summarize data to 10 min from 1 minute
met_tenmin <- met_EDI %>% 
  dplyr::filter(DateTime >= ymd_hms("2018-08-29 00:00:00"),
         DateTime <= ymd_hms("2022-03-01 23:50:00")) %>% 
  dplyr::select(DateTime, ShortwaveRadiationUp_Average_W_m2, WindSpeed_Average_m_s, PAR_umolm2s_Average) %>% 
  dplyr::mutate(Date = as.Date(DateTime)) %>% 
  dplyr::mutate(Hour = hour(DateTime),
         minutes = minute(DateTime)) %>% 
  dplyr::mutate(Minute_tenminutebreaks = ifelse(minutes >= 56 | minutes <= 05,  00, minutes),
         Minute_tenminutebreaks = ifelse(minutes >= 06 & minutes <= 15,  10, Minute_tenminutebreaks),
         Minute_tenminutebreaks = ifelse(minutes >= 16 & minutes <= 25,  20, Minute_tenminutebreaks),
         Minute_tenminutebreaks = ifelse(minutes >= 26 & minutes <= 35,  30, Minute_tenminutebreaks),
         Minute_tenminutebreaks = ifelse(minutes >= 36 & minutes <= 45,  40, Minute_tenminutebreaks),
         Minute_tenminutebreaks = ifelse(minutes >= 46 & minutes <= 55,  50, Minute_tenminutebreaks),
         Hour = ifelse(minutes >= 56, Hour+1, Hour)) %>% 
  dplyr::rename(Minute = Minute_tenminutebreaks) %>% 
  dplyr::mutate(DateTimeA = make_datetime(year = year(Date), month = month(Date), day = day(Date), hour = Hour, min = Minute, tz ="UTC")) %>% 
  dplyr::group_by(DateTimeA) %>% 
  dplyr::summarize(ShortwaveRadiationUp_Average_W_m2_tenminmean = mean(ShortwaveRadiationUp_Average_W_m2, na.rm = TRUE),
                   WindSpeed_Average_m_s_tenminmean = mean(WindSpeed_Average_m_s, na.rm = TRUE),
                   PAR_Average_umol_s_m2_tenminmean = mean(PAR_umolm2s_Average, na.rm  = T),
                   .groups = "drop") %>% 
  dplyr::rename(DateTime = DateTimeA)

met_tenmin$DateTime[duplicated(met_tenmin$DateTime)]


#converting our SW data to PAR for use in the metabfunction 
met_tenmin_PAR <- sw.to.par(met_tenmin, sw.col = "ShortwaveRadiationUp_Average_W_m2_tenminmean", coeff = 2.114)

#compare converted PAR to actual PAR 
# met_tenmin_PAR %>% 
#   ggplot(aes(x = PAR_Average_umol_s_m2_tenminmean, y = par ))+
#   geom_point()+
#   geom_abline(slope = 1, intercept = 0, size=2, linetype =2, col = "red")+ #1:1 line
#   geom_smooth(method = "lm", inherit.aes = T)+ 
#   stat_poly_eq(formula=y~x,label.x = "left",label.y="top",parse=TRUE,inherit.aes = T,
#                aes(x = PAR_Average_umol_s_m2_tenminmean, y = par, label=paste(..adj.rr.label..,..p.value.label..,sep="~~~"),size=3))


#2018 to 2022
PAR_from_SW_2018_22 <- met_tenmin_PAR %>%
  dplyr::select(DateTime, par) %>%
  dplyr::rename("dateTime" = DateTime,
         "PAR" = "par")

write.csv(PAR_from_SW_2018_22, file = "./Data/Model_Input/2018_22/FCR_2018_22_PAR_EXO.csv", row.names = FALSE)


windspeed_2018_22 <- met_tenmin %>%
  dplyr::select(DateTime, WindSpeed_Average_m_s_tenminmean) %>%
  dplyr::rename("dateTime" = DateTime,
         "windSpeed" = "WindSpeed_Average_m_s_tenminmean")

write.csv(windspeed_2018_22, file = "./Data/Model_Input/2018_22/FCR_2018_22_WindSpeed_EXO.csv", row.names = FALSE)




##### Next get 15 min par and windspeed for WVWA data 

#filter to desired variables and summarize data to 10 min from 1 minute
met_15min <- met_EDI %>% 
  dplyr::filter(DateTime >= ymd_hms("2015-11-08 00:00:00"),
         DateTime <= ymd_hms("2019-01-02 00:00:00")) %>% 
  dplyr::select(DateTime, ShortwaveRadiationUp_Average_W_m2, WindSpeed_Average_m_s) %>% 
  dplyr::mutate(Date = as.Date(DateTime)) %>% 
  dplyr::mutate(Hour = hour(DateTime),
         minutes = minute(DateTime)) %>% 
  dplyr::mutate(Minute_15minutebreaks = ifelse(minutes >= 53 | minutes <= 07,  00, minutes),
         Minute_15minutebreaks = ifelse(minutes >= 08 & minutes <= 22,  15, Minute_15minutebreaks),
         Minute_15minutebreaks = ifelse(minutes >= 23 & minutes <= 37,  30, Minute_15minutebreaks),
         Minute_15minutebreaks = ifelse(minutes >= 38 & minutes <= 52,  45, Minute_15minutebreaks),
         Hour = ifelse(minutes >= 53, Hour+1, Hour)) %>%  
  dplyr::rename(Minute = Minute_15minutebreaks) %>% 
  dplyr::mutate(DateTimeA = make_datetime(year = year(Date), month = month(Date), day = day(Date), hour = Hour, min = Minute, tz ="UTC")) %>% 
  dplyr::group_by(DateTimeA) %>% 
  dplyr::summarize(ShortwaveRadiationUp_Average_W_m2_15minmean = mean(ShortwaveRadiationUp_Average_W_m2, na.rm = TRUE),
                   WindSpeed_Average_m_s_15minmean = mean(WindSpeed_Average_m_s, na.rm = TRUE),
                   .groups = "drop")  %>% 
  dplyr::rename(DateTime = DateTimeA)

met_15min$DateTime[duplicated(met_15min$DateTime)]



#converting our SW data to PAR for use in the metabfunction 
met_15min_PAR <- sw.to.par(met_15min, sw.col = "ShortwaveRadiationUp_Average_W_m2_15minmean", coeff = 2.114)

met_15min_PAR <- met_15min_PAR %>% 
  dplyr::rename(PAR = par,
                windSpeed = WindSpeed_Average_m_s_15minmean)


#bring in NLDAS regression data to gap fill missing time periods (gaps; 140ct2016 - 19dec2016, 6feb - 24apr2017, 19mar - 5apr2018)
nldas <- read.csv("./Data/Generated_Data/Met_input_from_NLDAS.csv")


missing_met_dates <- c(seq(ymd_hms("2016-10-14 16:00:00"), ymd_hms("2016-12-19 14:00:00"), by = "1 hour"),
                        seq(ymd_hms("2017-02-06 10:00:00"), ymd_hms("2017-04-24 13:00:00"), by = "1 hour"),
                        seq(ymd_hms("2018-03-19 11:00:00"), ymd_hms("2018-04-05 13:00:00"), by = "1 hour"))

nldas_gapfill <- nldas %>% 
  dplyr::mutate(DateTime = ymd_hms(time)) %>% 
  dplyr::filter(DateTime %in% missing_met_dates) %>%   
  dplyr::select(DateTime, WindSpeed_regress, PAR_from_SWregress) %>%
  dplyr::rename(PAR = PAR_from_SWregress,
                windSpeed = WindSpeed_regress)

met_15min_PAR_nldasfilled <- dplyr::full_join(met_15min_PAR, nldas_gapfill, by = c("DateTime", "PAR", "windSpeed")) %>% 
  dplyr::arrange(DateTime)

met_15min_PAR_nldasfilled$DateTime[duplicated(met_15min_PAR_nldasfilled$DateTime)]

#PAR and WindSpeed 2015 to 2018
PAR_from_SW_2015_18 <- met_15min_PAR_nldasfilled %>%
  dplyr::select(DateTime, PAR) %>%
  dplyr::rename(dateTime = DateTime)

write.csv(PAR_from_SW_2015_18, file = "./Data/Model_Input/2015_18/FCR_2015_18_PAR_15min_withNLDAS.csv", row.names = FALSE)


windspeed_2015_18 <- met_15min_PAR_nldasfilled %>%
  dplyr::select(DateTime, windSpeed) %>%
  dplyr::rename(dateTime = DateTime)

write.csv(windspeed_2015_18, file = "./Data/Model_Input/2015_18/FCR_2015_18_WindSpeed_15min_withNLDAS.csv", row.names = FALSE)


####  Writing temp profiles from HOBOs for WVWA period ####
hobo_edi <- read.csv("./Data/EDI2023/FCR_hobos_2015_2018.csv")

hobo <- hobo_edi %>% 
  dplyr::select(3:12) %>% 
  dplyr::mutate(DateTime = ymd_hms(DateTime)) %>% 
  dplyr::mutate(Date = as.Date(DateTime),
         Hour = hour(DateTime),
         minutes = minute(DateTime)) %>% 
  dplyr::mutate(Minute_15minutebreaks = ifelse(minutes >= 53 | minutes <= 07,  00, minutes),
         Minute_15minutebreaks = ifelse(minutes >= 08 & minutes <= 22,  15, Minute_15minutebreaks),
         Minute_15minutebreaks = ifelse(minutes >= 23 & minutes <= 37,  30, Minute_15minutebreaks),
         Minute_15minutebreaks = ifelse(minutes >= 38 & minutes <= 52,  45, Minute_15minutebreaks),
         Hour = ifelse(minutes >= 53, Hour+1, Hour)) %>%  #This is so that a time like 12:56 doesn't get classified as 12:00 when it should be 13:00
  dplyr::rename(Minute = Minute_15minutebreaks) %>% 
  dplyr::group_by(Date, Hour, Minute) %>% 
  dplyr::summarize(temp1.0 = mean(wtr_1, na.rm = T),
            temp2.0 = mean(wtr_2, na.rm = T),
            temp3.0 = mean(wtr_3, na.rm = T),
            temp4.0 = mean(wtr_4, na.rm = T),
            temp5.0 = mean(wtr_5, na.rm = T),
            temp6.0 = mean(wtr_6, na.rm = T),
            temp7.0 = mean(wtr_7, na.rm = T),
            temp8.0 = mean(wtr_8, na.rm = T),
            temp9.3 = mean(wtr_9.3, na.rm = T),
            .groups = "drop") %>% 
  dplyr::mutate(dateTime = make_datetime(year = year(Date), month = month(Date), day = day(Date), hour = Hour, min = Minute, tz ="UTC")) %>% 
  dplyr::select(dateTime, temp1.0, temp2.0, temp3.0, temp4.0, temp5.0, temp6.0, temp7.0, temp8.0, temp9.3)

profiles2015_18 <- hobo %>%
  dplyr::filter(dateTime >= ymd_hms("2015-11-08 00:00:00"))

#check for duplicated datetimes
profiles2015_18$dateTime[duplicated(profiles2015_18$dateTime)]

#summarize across duplicated datetimes
profiles2015_18 <- profiles2015_18 %>% 
  dplyr::group_by(dateTime) %>% 
  dplyr::summarise(across(.cols = where(is.numeric), .fns = ~mean(.x, na.rm = TRUE)), .groups = "drop")

#write csv
#write.csv(profiles2015_18, file = "./Data/Model_Input/2015_18/FCR_2015_18_TempProfiles_hobos.csv", row.names = FALSE)



#### Making input files w/ GLM and interpolated 2018 #####

##read in data 
glmtemp_hourly <- read_csv("./Data/Carey2022_glm/ModeledWaterTempFCRForDexter_22Aug2022.csv")

#change column names 
temp_profiles_glm_hourly <- glmtemp_hourly %>% 
  dplyr::select(-1) %>% 
  dplyr::rename(dateTime = DateTime,
         temp0.1 = temp_0.1,
         temp1.0 = temp_1,
         temp2.0 = temp_2,
         temp3.0 = temp_3,
         temp4.0 = temp_4,
         temp5.0 = temp_5,
         temp6.0 = temp_6,
         temp7.0 = temp_7,
         temp8.0 = temp_8,
         temp9.0 = temp_9)

#2018 to 2019
temp_profiles_2018_19_glm_hourly <- temp_profiles_glm_hourly %>% 
  dplyr::filter(dateTime >= ymd_hms("2018-01-01 00:00:00"),
         dateTime <= ymd_hms("2020-01-01 23:50:00")) 


### Compare glm hourly to catwalk thermistors to get RMSE 

#use 2019-2020 for comp
temp_prof_18_22_long <- temp_profiles_2018_22 %>% 
  pivot_longer(c(-1), names_to = "name", values_to = "Catwalk")

temp_glm_long <- temp_profiles_glm_hourly %>% 
  filter(dateTime >= ymd("2019-01-01"),
         dateTime <= ymd("2020-12-31")) %>% 
  pivot_longer(c(-1), names_to = "name", values_to = "GLM_hourly")

joined_glm_cat_rmse <- left_join(temp_glm_long, temp_prof_18_22_long, by = c("dateTime", "name")) %>% 
  dplyr::mutate(Diff = Catwalk - GLM_hourly) %>% 
  dplyr::filter(!is.na(Catwalk)) #getting rid of NAs to run rmse 

plot(joined_glm_cat_rmse$dateTime, joined_glm_cat_rmse$Diff)

round(rmse(joined_glm_cat_rmse$GLM_hourly, joined_glm_cat_rmse$Catwalk), digits = 3)


### getting january hobos
jan18_hobos <- profiles2015_18 %>% 
  dplyr::filter(dateTime >= ymd("2017-12-31"))

### getting july and onward SCC Sensors 

scc_thermistors <- temp_profiles_2018_22 %>% 
  dplyr::filter(dateTime <= ymd("2019-01-02"))

### Find data for in between these two 

## YSI 
ysi <- read_csv("./Data/Generated_Data/FCR_YSI_do_temp.csv")

ysi18_temp_profiles_cleaned <- ysi %>% 
  dplyr::filter(DateTime > ymd("2017-12-31"),
         DateTime < ymd("2019-01-01")) %>%
  dplyr::filter(!is.na(Temp_C)) %>% 
  dplyr::select(DateTime, Depth_m, Temp_C) %>% 
  pivot_wider(names_from = Depth_m, values_from = Temp_C) %>% 
  dplyr::select(DateTime, "0.1", "0.8", "1.6", "2.8", "3.8", "5", "6.2", "8", "9") %>% 
  dplyr::rename(dateTime = DateTime,
         temp0.1 = "0.1",
         temp0.8 = "0.8", 
         temp1.6 = "1.6", 
         temp2.8 = "2.8", 
         temp3.8 = "3.8", 
         temp5.0 = "5", 
         temp6.2 = "6.2", 
         temp8.0 = "8", 
         temp9.0 = "9")

## CTD 
ctd <- read_csv("./Data/Generated_Data/FCR_CTD_temp_do.csv")
head(ctd)

ctd18_temp_profiles_cleaned <- ctd %>% 
  dplyr::filter(Date > ymd("2017-12-31"),
         Date < ymd("2019-01-01")) %>% 
  dplyr::filter(!is.na(Temp_C)) %>% 
  dplyr::select(Date, Depth_m, Temp_C) %>% 
  dplyr::filter(Depth_m <= 9) %>% 
  pivot_wider(names_from = Depth_m, values_from = Temp_C) %>%   
  dplyr::select(Date, "0.1", "0.7", "1.6", "2.9", "3.8", "5", "6.2", "8", "9" ) %>% 
  dplyr::rename(dateTime = Date,
         temp0.1 = "0.1",
         temp0.8 = "0.7", 
         temp1.6 = "1.6", 
         temp2.8 = "2.9", 
         temp3.8 = "3.8", 
         temp5.0 = "5", 
         temp6.2 = "6.2", 
         temp8.0 = "8",
         temp9.0 = "9")

## Put CTD and YSI together 
ctd_ysi_binded_18 <- rbind(ysi18_temp_profiles_cleaned, ctd18_temp_profiles_cleaned)
ctd_ysi_binded_18 <- dplyr::arrange(ctd_ysi_binded_18, dateTime)

ctd_ysi_binded_18_daily <- ctd_ysi_binded_18 %>% 
  pivot_longer(cols = c(2:10)) %>% 
  dplyr::mutate(Date = as.Date(dateTime)   ,
         dateTimeTest = ymd_hms(paste(Date, "12:00:00", sep = ""))
  ) %>% 
  dplyr::group_by(dateTimeTest, name) %>% 
  dplyr::summarize(meantemp = mean(value, na.rm = T),
            .groups = "drop") %>% 
  pivot_wider(names_from = name, values_from = meantemp) %>% 
  dplyr::rename(dateTime = dateTimeTest)


## Interpolating for CTD and YSI timeframe 
# wvwa_dates <- as.data.frame(seq(ymd_hms("2018-01-15 12:00:00"), ymd_hms("2018-07-05 12:00:00"), by="days"))
wvwa_dates <- as.data.frame(seq(ymd_hms("2018-01-15 12:00:00"), ymd_hms("2018-12-17 12:00:00"), by="days"))
colnames(wvwa_dates)[1] <- "dateTime"


ctd_ysi_binded_18_daily_interpolated <- left_join(wvwa_dates, ctd_ysi_binded_18_daily, by = "dateTime")

ctd_ysi_binded_18_daily_interpolated$temp0.1 <- na.approx(ctd_ysi_binded_18_daily_interpolated$temp0.1)
ctd_ysi_binded_18_daily_interpolated$temp0.8 <- na.approx(ctd_ysi_binded_18_daily_interpolated$temp0.8)
ctd_ysi_binded_18_daily_interpolated$temp1.6 <- na.approx(ctd_ysi_binded_18_daily_interpolated$temp1.6)
ctd_ysi_binded_18_daily_interpolated$temp2.8 <- na.approx(ctd_ysi_binded_18_daily_interpolated$temp2.8)
ctd_ysi_binded_18_daily_interpolated$temp3.8 <- na.approx(ctd_ysi_binded_18_daily_interpolated$temp3.8)
ctd_ysi_binded_18_daily_interpolated$temp5.0 <- na.approx(ctd_ysi_binded_18_daily_interpolated$temp5.0)
ctd_ysi_binded_18_daily_interpolated$temp6.2 <- na.approx(ctd_ysi_binded_18_daily_interpolated$temp6.2)
ctd_ysi_binded_18_daily_interpolated$temp8.0 <- na.approx(ctd_ysi_binded_18_daily_interpolated$temp8.0)
ctd_ysi_binded_18_daily_interpolated$temp9.0 <- na.approx(ctd_ysi_binded_18_daily_interpolated$temp9.0)

##compare CTD/YSI to 
binded_dailyinterp_forComp <- ctd_ysi_binded_18_daily_interpolated %>% 
  dplyr::select(dateTime, temp0.1, temp5.0, temp8.0, temp9.0) %>% 
  pivot_longer(-1)

scc_thermistors_long <- scc_thermistors %>% 
  dplyr::select(dateTime, temp0.1, temp5.0, temp8.0, temp9.0) %>% 
  pivot_longer(-1)

dailyinterp_toCatwalk_comp <- left_join(binded_dailyinterp_forComp, scc_thermistors_long, by = c("dateTime", "name")) %>% 
  dplyr::rename(ctd_ysi_interp = value.x,
         Catwalk = value.y) %>% 
  dplyr::mutate(Diff = Catwalk - ctd_ysi_interp) %>% 
  dplyr::filter(!is.na(Catwalk)) #getting rid of NAs to run rmse 

plot(dailyinterp_toCatwalk_comp$dateTime, dailyinterp_toCatwalk_comp$Diff)

round(rmse(dailyinterp_toCatwalk_comp$ctd_ysi_interp, dailyinterp_toCatwalk_comp$Catwalk), digits = 3)

## Put all data together 

temp_profiles_18gap_glm_hourly <- temp_profiles_glm_hourly %>% 
  filter(dateTime >= ymd_hms("2018-01-15 00:00:00"),
         dateTime <= ymd_hms("2018-04-22 23:50:00")) 


ctd_ysi_binded_18_daily_interpolated_forjoin <- ctd_ysi_binded_18_daily_interpolated %>% 
  filter(dateTime >= ymd("2018-04-23"),
         dateTime <= ymd("2018-07-06"))

joined <- dplyr::bind_rows(jan18_hobos, temp_profiles_18gap_glm_hourly, ctd_ysi_binded_18_daily_interpolated_forjoin, scc_thermistors)

profiles2018_final <- joined %>% 
  select(dateTime, temp0.1, temp0.8, temp1.0, temp1.6, temp2.0, temp2.8, temp3.0, temp3.8,
         temp4.0, temp5.0, temp6.0, temp6.2, temp7.0, temp8.0, temp9.0, temp9.3)


#write temp profile for 2018 for wvwa sondes 
#write.csv(profiles2018_final, file = "./Data/Model_Input/2015_18/FCR_2018_TempProfiles_hobo_glm_ctdysi_scc.csv", row.names = FALSE) 

#make a csv that has hobos 2015 - jan2018 the interpolated values from csv above 

profiles2015_18_forbind <- profiles2015_18 %>% 
  dplyr::filter(dateTime < ymd_hms("2017-12-31 00:00:00")) %>% 
  dplyr::mutate(temp0.1 = NA, temp0.8 = NA, temp1.6 = NA, temp2.8 = NA, temp3.8 = NA, temp6.2=NA, temp9.0=NA)

profiles2015_18_withInterpolated2018 <- rbind(profiles2015_18_forbind, profiles2018_final)

profiles2015_18_withInterpolated2018 <- profiles2015_18_withInterpolated2018 %>% 
  select(dateTime, temp0.1, temp0.8, temp1.0, temp1.6, temp2.0, temp2.8, temp3.0, temp3.8,
       temp4.0, temp5.0, temp6.0, temp6.2, temp7.0, temp8.0, temp9.0, temp9.3)


write.csv(profiles2015_18_withInterpolated2018, file = "./Data/Model_Input/2015_18/FCR_2015_18_TempProfiles_hobos_and_interped2018.csv", row.names = FALSE) 
