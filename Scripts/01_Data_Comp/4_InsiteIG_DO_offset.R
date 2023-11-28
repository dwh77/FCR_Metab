### Compiling WVWA DO sensor data, QAQC, and comparing to EXO and YSI for offset 
### Dexter Howard 
### 29Jan2021
### updated again by DWH in September 2022 


library(tidyverse)
library(patchwork)

#### DEVELOPING OFFSET FOR WVWA DATA FROM EXO #### 

#### read in wvwa DO sonde data from edi
dosonde_edi <- read_csv("./Data/EDI2023/FCR_DOsondes_2012_2018.csv")


WVWA_QAQC <- dosonde_edi %>% 
  dplyr::select(DateTime, Temp_C_1m, DO_mgL_1m) %>% 
  dplyr::filter(DateTime > ymd("2015-11-09"))

#### Pulling in EXO data 

catwalk_EDI <- read_csv("./Data/EDI2023/FCR_Catwalk_2018_2022.csv", col_types = cols(.default = "d", Reservoir = "c", DateTime = "T"))

exo <- catwalk_EDI %>% 
  dplyr::select(DateTime, EXODO_mgL_1) %>% 
  dplyr::filter(DateTime >= ymd("2018-08-29"),
         DateTime <= ymd("2019-01-01")) 


#### pulling in YSI data 
ysi <- read_csv("./Data/FCR_YSI_do_temp_for_metabolism.csv")

ysi_1_16 <- ysi %>% 
  dplyr::filter(Depth_m %in% c(1.0, 1.6)) %>% 
  dplyr::filter(DateTime >= ymd_hms("2018-08-01 00:00:00"),
         DateTime <= ymd_hms("2019-01-01 00:00:00"))


#### Comparing overlap of WVWA and EXO 

wvwa_comp <- WVWA_QAQC %>% 
  dplyr::filter(DateTime >= ymd_hms("2018-08-01 00:00:00"),
         DateTime <= ymd_hms("2019-01-01 00:00:00")) %>% 
  dplyr::rename(dateTime = DateTime,
                DO1m_Dissolved_Oxygen_ppm = DO_mgL_1m) 


exo_comp <- exo %>% 
  dplyr::rename(dateTime = DateTime,
                DO = EXODO_mgL_1) %>% 
  dplyr::filter(dateTime >= ymd_hms("2018-08-01 00:00:00"),
         dateTime <= ymd_hms("2019-01-01 00:00:00")) 


overlap_nocorrection <- ggplot()+
  geom_line(data = wvwa_comp, mapping = aes(x = dateTime, y = DO1m_Dissolved_Oxygen_ppm, color = "InsiteIG"))+
  geom_line(data = exo_comp, mapping = aes(x = dateTime, y = DO, color = "EXO"))+
  # geom_point(data = ysi_1_16, mapping = aes(x = DateTime, y = DO_mgL, color = "YSI"))+
  labs(x = "Date", y = "DO (mg/L)")+
  ggtitle("DO sonde comparison - no correction")+
  theme_classic()

overlap_nocorrection

## Calculating offset 

diff <- left_join(wvwa_comp, exo_comp, by = "dateTime")

diff <- diff %>% 
  dplyr::select(dateTime, DO, DO1m_Dissolved_Oxygen_ppm) %>% 
  dplyr::mutate(DO_diff = DO - DO1m_Dissolved_Oxygen_ppm)

plot(diff$dateTime, diff$DO_diff, main = "EXO DO - WVWA DO")

summary(diff$DO_diff)

overlap_correction <- ggplot()+
  geom_line(data = wvwa_comp, mapping = aes(x = dateTime, y = DO1m_Dissolved_Oxygen_ppm + 1.687, color = "InsiteIG"))+
  geom_line(data = exo_comp, mapping = aes(x = dateTime, y = DO, color = "EXO"))+
  # geom_point(data = ysi_1_16, mapping = aes(x = DateTime, y = DO_mgL, color = "YSI"))+
  labs(x = "Date", y = "DO (mg/L)")+
  ggtitle("DO sonde comparison - correction \n Insite_corr = Insite + 1.687")+
  theme_classic()

overlap_correction

overlap_comparison <- overlap_nocorrection / overlap_correction

overlap_comparison

ggsave(filename = "./Figures/MS_final/DO_offset_overlap_comparison.png", 
       overlap_comparison, device = "png", width = 230, height = 150, units = "mm")
  


#### Writing WVWA DO data for model based on EXO offset ####
#Using offset of WVWA + 1.687 based on difference between EXO and WVWA from sensor deployment to 2019-01-01


wvwa_DO_final <- WVWA_QAQC %>% 
  dplyr::mutate(DO_corrected = DO_mgL_1m + 1.687)

wvwa_DO_final <- wvwa_DO_final %>% 
  dplyr::select(DateTime, Temp_C_1m, DO_corrected) %>% 
  dplyr::rename(dateTime = DateTime,
         sensorTemp = Temp_C_1m,
         DO = DO_corrected) 
  
# DO and Temp 2015 to 2018

do_2015_18 <- wvwa_DO_final %>% 
  dplyr::filter(dateTime >= ymd_hms("2015-11-09 16:00:00")) %>% 
  dplyr::select(dateTime, DO) 

write.csv(do_2015_18, file = "./Data/Model_Input/2015_18/FCR_2015_18_wvwa_DO.csv", row.names = FALSE)

temp_2015_18 <- wvwa_DO_final %>% 
  dplyr::filter(dateTime >= ymd_hms("2015-11-09 16:00:00")) %>% 
  dplyr::select(dateTime, sensorTemp)

write.csv(temp_2015_18, file = "./Data/Model_Input/2015_18/FCR_2015_18_wvwa_sensorTemp.csv", row.names = FALSE)



