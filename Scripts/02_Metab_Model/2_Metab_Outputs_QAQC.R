#### QAQC metab outputs 
#DWH september 2022

#load packages
library(tidyverse)

#check wd after running models
getwd()
setwd(here::here())
getwd()


#### read in data ####
## EXO data 
metaboutput2018_22 <- read.table("./Data/Model_Output/FCR2018_22/FCR2018_22 outputDCR.txt",header=T,sep='\t')
head(metaboutput2018_22)
metaboutput2018_22 <- metaboutput2018_22 %>% 
  mutate(solarDay = as.Date(solarDay))

##InsiteIG (wvwa) data 
metaboutput2015_18 <- read.table("./Data/Model_Output/FCR2015_18/FCR2015_18 outputDCR.txt",header=T,sep='\t')
metaboutput2015_18 <- metaboutput2015_18 %>% 
  mutate(solarDay = as.Date(solarDay)) %>% 
  filter(solarDay < ymd("2018-08-29"))

#bind together
metab_all <- rbind(metaboutput2015_18, metaboutput2018_22)


#### QAQC ####
#check initial plots 
plot(metab_all$solarDay, metab_all$GPP_mgO2perLperday, xlab="Solar Day", ylab="GPP, (mg O2 L-1 d-1)")
plot(metab_all$solarDay, metab_all$R_mgO2perLperday, xlab="Solar Day", ylab="R, (mg O2 L-1 d-1)")
plot(metab_all$solarDay, metab_all$NEM_mgO2perLperday, xlab="Solar Day", ylab="NEM, (mg O2 L-1 d-1)")

#Set values < 0.001 to NA and remove NEM when R or GPP is an NA 
metab_all_QAQC <- metab_all %>% 
  mutate(GPP_QAQC = ifelse(GPP_mgO2perLperday < 0.001, NA, GPP_mgO2perLperday),
         R_QAQC = ifelse(R_mgO2perLperday < 0.001, NA, R_mgO2perLperday),
         NEM_QAQC = ifelse( is.na(GPP_QAQC) | is.na(R_QAQC), NA, NEM_mgO2perLperday) ) 

#number of NAs 
R_na <- metab_all_QAQC %>% 
  filter(is.na(R_QAQC))

GPP_na <- metab_all_QAQC %>% 
  filter(is.na(GPP_QAQC))

NEM_na <- metab_all_QAQC %>% 
  filter(is.na(NEM_QAQC))

#removing days w/ NEM NAs 
metab_QAQCfin <- metab_all_QAQC %>% 
  filter(!is.na(NEM_QAQC)) 

#plot intial QAQC
metab_QAQCfin %>% 
  select(solarDay, GPP_QAQC, NEM_QAQC, R_QAQC) %>% 
  pivot_longer(cols = c(2:4)) %>% 
  ggplot(aes(x = solarDay, y = value))+
  geom_point()+
  facet_wrap(~name, scales = "free_y", nrow = 3)+
  theme_classic()

#check on extreme R outliers,  
R_outlier_check_greater15 <- metab_QAQCfin %>% 
  filter(R_QAQC > 15)
head(R_outlier_check_greater15)

#remvoign four values over 19 and after checking daily model fits, and these had poor model fits 

metab_QAQCfin2 <- metab_QAQCfin %>% 
  filter(R_QAQC < 19)

metab_QAQCfin2 %>% 
  select(solarDay, GPP_QAQC, NEM_QAQC, R_QAQC) %>% 
  pivot_longer(cols = c(2:4)) %>% 
  ggplot(aes(x = solarDay, y = value))+
  geom_point()+
  facet_wrap(~name, scales = "free_y", nrow = 3)+
  theme_classic()


#writing csv for export
write.csv(metab_QAQCfin2, "./Data/Model_Output/MetabOutput_QAQC_15_22.csv", row.names = FALSE)











