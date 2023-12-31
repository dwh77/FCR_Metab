### Compiling CTD and YSI data for FCR metabolism

### load packages needed ####
library(tidyverse)


#### YSI processing #### 
#read in ysi data
ysi_edi <- read.csv("./Data/EDI2023/YSI_PAR_profiles_2013-2022.csv")

#check YSI for missing 2018 timeframe 
ysi2018 <- ysi_edi %>% 
  filter(Reservoir == "FCR",
         Site == 50,
         Depth_m > 0,
         DateTime >= ymd("2018-01-10"), #these dates are to check if theres YSI casts in timeframe where we're missing thermistor temp profiles
         DateTime <= ymd("2018-07-10"),
         !is.na(Temp_C)) 

ysi_2018dates <- unique(ysi2018$DateTime)
ysi_2018dates
ysi_2018dates <- ymd_hms(ysi_2018dates)

#filer YSI to desired reservoir, site, timeframe, and variables
ysi_fcrmetab <- ysi_edi %>% 
  filter(Reservoir == "FCR",
         Site == 50,
         Depth_m > 0,
         DateTime >= ymd_hms("2015-11-01 00:00:00")) %>% 
  select(Reservoir, Site, DateTime, Depth_m, Temp_C, DO_mgL, DOsat_percent, Flag_Temp_C, Flag_DO_mgL, Flag_DOsat_percent)


#filter desired data and write csv. this csv will be used in data comp script to get by year data
ysi_fcrmetab_write <- ysi_fcrmetab %>% 
  select(DateTime, Depth_m, Temp_C, DO_mgL, DOsat_percent)

#create directory to write data and make csv
dir.create("./Data/Generated_Data")
write.csv(ysi_fcrmetab_write, "./Data/Generated_Data/FCR_YSI_DO_temp.csv", row.names = FALSE)



#### CTD processing ####
#read in data
ctd_edi <- read.csv("./Data/EDI2023/CTD_final_2013_2022.csv")

#change Date from character to datetime
ctd_edi$Date <- ymd_hms(ctd_edi$Date)

#check for casts on dates we're missing thermistors in 2018
ctd_2018 <- ctd_edi %>% 
  filter(Reservoir == "FCR",
         Site == 50,
         Depth_m > 0,
         Date >= ymd("2018-01-10"), 
         Date <= ymd("2018-07-10")) 

ctd_2018dates <- unique(ctd_2018$Date)


#filter to desired data 
ctd_fcrmetab <- ctd_edi %>% 
  filter(Reservoir == "FCR",
         Site == 50,
         Depth_m > 0,
         Date >= ymd_hms("2015-10-22 00:00:00")) %>%
  select(Reservoir, Site, Date, Depth_m, Temp_C, DO_mgL, DOsat_percent, Flag_Temp_C, Flag_DO_mgL, Flag_DOsat_percent)


#Process CTD to make it easier to work with depths

#rename ctd data to fit with prior code for processing
ctd <- ctd_fcrmetab


# filter out depths in the CTD cast that are closest to these specified values.
df.final<-data.frame()
ctd1<-ctd %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 0.1)))
ctd2<-ctd %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 0.4)))
ctd3<-ctd %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 0.7)))
ctd4<-ctd %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 1)))
ctd5<-ctd %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 1.3)))
ctd6<-ctd %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 1.6)))
ctd7<-ctd %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 1.9)))
ctd8<-ctd %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 2.3)))
ctd9<-ctd %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 2.6)))
ctd10<-ctd %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 2.9)))
ctd11<-ctd %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 3.2)))
ctd12<-ctd %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 3.5)))
ctd13<-ctd %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 3.8)))
ctd14<-ctd %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 4.1)))
ctd15<-ctd %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 4.4)))
ctd16<-ctd %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 4.7)))
ctd17<-ctd %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 5)))
ctd18<-ctd %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 5.3)))
ctd19<-ctd %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 5.6)))
ctd20<-ctd %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 5.9)))
ctd21<-ctd %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 6.2)))
ctd22<-ctd %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 6.5)))
ctd23<-ctd %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 6.8)))
ctd24<-ctd %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 7.1)))
ctd25<-ctd %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 7.4)))
ctd26<-ctd %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 7.7)))
ctd27<-ctd %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 8)))
ctd28<-ctd %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 8.3)))
ctd29<-ctd %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 8.7)))
ctd30<-ctd %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 9)))
ctd31<-ctd %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 9.3)))
ctd32<-ctd %>% group_by(Date) %>% slice(which.min(abs(as.numeric(Depth_m) - 9.6)))


# Bind each of the data layers together.
df.final = rbind(ctd1,ctd2,ctd3,ctd4,ctd5,ctd6,ctd7,ctd8,ctd9,ctd10,ctd11,ctd12,ctd13,ctd14,ctd15,ctd16,ctd17,ctd18,ctd19,
                 ctd20,ctd21,ctd22,ctd23,ctd24,ctd25,ctd26,ctd27,ctd28,ctd29,ctd30,ctd31, ctd32)

# Re-arrange the data frame by date
ctd <- arrange(df.final, Date)
ctd$Depth_m <- round(as.numeric(ctd$Depth_m), digits = 1) 
ctd <- ctd[!duplicated(ctd),] #removed duplicated rows 


#drop flags since there's no flag issues besides days sample wasn't collected for DO 
ctd <- ctd %>% 
  select(-Flag_Temp_C, -Flag_DO_mgL, -Flag_DOsat_percent)


#write filtered ctd csv

write.csv(ctd, "./Data/Generated_Data/FCR_CTD_temp_DO.csv", row.names = FALSE)

  



