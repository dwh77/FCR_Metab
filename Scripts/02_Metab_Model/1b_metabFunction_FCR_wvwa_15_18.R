#metabFunc_30Jul
#Script for running metabolism model
#CTS 14 May 2009
#Updated 30 Jul 2009
#Updated 07Nov2011 for the Ashokan reservoir to deal with lots of holes in DO and other data
#Updated 07Feb2013 DCR to including 
#extending out the timing of the each day +- timing is specified in the script file
#Adding in an additional variable that will esitmate the offset of the sensor
#this will be reflected in the new file: metabLoss_v5_DCR.R
#Updated 09Apr2013
#Includes a vector called ModelVersion<-c(1,0,3,0)
#Which includes variables to determine if there is DO initial fitting, simple/photoinhib for GPP, hours extended, and atm flux for Ice periods
#now associated with metabLoss_v7 and metabPredix_v7
#Updated 17Jun2013
#Includes a vector called ModelVersion<-c(1,0,3,0,0)
#Which includes variables to determine if there is DO initial fitting, simple/photoinhib for GPP, hours extended, and atm flux for Ice periods
#now associated with metabLoss_v8 and metabPredix_v8
#Updated 22Feb2018 to run 2007-2009 continuous under ice data
#Runs complete data set from August 2007 to Jan 2009 including under ice
#Updated 14May2018 to incorporate the time series DO correction

#This model version array tells the model what to run for parameter fits
#ModelVersion[1]
#If 1, then fit DO initial
#If 2, then DO initial is first DO value of the day
#If 3, then DO initial is the average of the first 1 hour (6 time periods)
#ModelVersion[2]
#If 0, then simple model for GPP
#If 1, then include photoinhibition model
#ModelVersion[3]
#The number is the number of hours to extend the day into the preceeding day and the following sunset
#This accounts for the fact that PAR starts to go up before the sunrise day that is calculated
#ModelVersion[4]
#If 0, then no atm flux during ice periods, requires the iceOut.vector and iceOn.vector with dates
#If 1, then treat ice periods like any other with respect to atmospheric flux
#ModelVersion[5]
#If 0, then read in the old DO and thermister data from SunapeeXXXX data/Proofed Data/
#If 1, then run the script Sunapee 2012 MergeDataDOproofing.R and read in the combed DO and thermister data from SunapeeXXXX data/CombedDeep/
#Turn off the DO offset!!
#In CombedDeep, the thermister data is now FullWtr_2007.wtr
#In CombedDeep, the DO data is now Sunapee_XXXX_CombedDO.txt
#If 2, then run the script Sunapee 2012 DOproofing MEAN SAT.R and read in the combed DO and thermister data from SunapeeXXXX data/CombedDeep/
#Turn off the DO offset!!
#In CombedDeep, the thermister data is now FullWtr_2007.wtr
#In CombedDeep, the DO data is now Sunapee_XXXX_CombedDO.txt
#If 3, then just use the usual DO and do not estimate an offset
#ModelVersion[6]
#This is the number of hours after sunrise that we should remove
#This is to eliminate that odd blip at the beginning of many days

### Modfied March 2021 by DWH for FCR ###

## check working directory, if its not at base of project ./FCR_Metabolism, manually reset to this point
getwd()

##check this addition from FEO works 
setwd(here::here())

#load packages
#ad in if from line 88
if(!require(tidyverse)){install.packages("tidyverse")}

library(tidyverse)


#Inputs are: #All of these should be loaded to the workspace before running this script
ModelVersion<-c(3,0,0,0,0,0) #pulling from line 15, adding 6 as 0 based on line ~644 in for loop. Was (1,0,3,0,0,0) 

lat <- 37.30       # Latitude of lake, decimal degrees, north positive
long <- -79.84      # Longitude of lake: wasn't in script before, adding just have
elev  <- 1663 * 0.3048       # Elevation of lake surface, m above sea level. Coverting from ft to m based on Water authortiy elevation values
#windHeight: ignoring for now since this is corrected for in FCR met data  # Height above lake surface at which wind speed is measured, m
timeStep  <- 15   # Time interval between DO measurements, minutes
sensorDepth <- 1.0  # Depth of DO sensor, m
outName <- "FCR2015_18"    # Text to use in labeling outputs, e.g. 'Acton2008'. Character.
dirData  <- "./Data/Model_Input/2015_18"    # Directory where data file are located, e.g. 'C:/GLEON/Acton'. Character.
dataIn  <- c("FCR_2015_18_wvwa_DO.csv", "FCR_2015_18_wvwa_sensorTemp.csv", 
             "FCR_2015_18_TempProfiles_hobos_and_interped2018.csv",
             "FCR_2015_18_PAR_15min_withNLDAS.csv", "FCR_2015_18_WindSpeed_15min_withNLDAS.csv")     # Names of files that contain the data. This should be a character vector of
#   length 5, e.g. c('Acton_2008_DO.txt','Acton_2008_PAR_5min.txt','Acton_2008_windSpeed.txt',
#   'Acton_2008_sensorTemp.txt','Acton_2008_tempProfile.txt')
dirFxns  <- "./Scripts/Metab_Model"    # Directory where functions are located, e.g. 'C:/GLEON/functions'
dirDump  <- "./Data/Model_Output/FCR2015_18a"   # Directory where outputs should be dumped, e.g. 'C:/GLEON/Acton/Results'
dir.create(dirDump)

# ice <- read.csv("./Data/Ice_Data_2013_2022.csv")
# ice$Date <- as.Date(ice$Date)
ice <- read_csv("./Data/Seasons_Ice_corrected_groupings.csv")


##DCR: here we are expanding out the day to include points before and after sunrises
#Hours before and after sunrise. We're using just 24 hours for now
hours.BeforeAfter<-ModelVersion[3]

#Initialized dataframe each year to export data for analysis
DF.residuals.Analysis.Year<-as.data.frame(matrix(nrow=0,ncol=12))

#Load in the suncalc package for sunrise and sunsets
if(!require(suncalc)){install.packages("suncalc")}
library(suncalc)

########################################
#Set up

#Load in required functions
#setwd(dirFxns), not setting wd here so that script can be run on any computer
source('./Scripts/Metab_Model/metabLoss_v8.R')
source('./Scripts/Metab_Model/metabPredix_v8.R')
source('./Scripts/Metab_Model/fillHoles.R')
source('./Scripts/Metab_Model/calcZMixDens.R')

#Set environment tz variable to GMT
Sys.setenv(tz="EST")


########################################
#Read and organize data

##
#Read data
# setwd(dirData) skipping changing directory so each of us can work of off exact same script
# View(dataIn) #view data in to ensure you're reading in right data set each time from subsetting functions 

#Divide PAR by 1000 to convert from measured units (umol m-2 s-1) to model units (mmol m-2 s-1)
#Note change for Acton here - na.strings="na"
# dataPAR <- read.table(dataIn[4],header=T,sep='\t',colClasses=c(dateTime="POSIXct"))
dataPAR <- read.csv("./Data/Model_Input/2015_18/FCR_2015_18_PAR_15min_withNLDAS.csv")
dataPAR$dateTime <- as.POSIXct(dataPAR$dateTime)
dataPAR$PAR[dataPAR$PAR<2]<-0
dataPAR$PAR <- dataPAR$PAR/1000

#Wind speed
# dataWind <- read.table(dataIn[5],header=T,sep='\t',colClasses=c(dateTime="POSIXct"))
dataWind <- read.csv("./Data/Model_Input/2015_18/FCR_2015_18_WindSpeed_15min_withNLDAS.csv")
dataWind$dateTime <- as.POSIXct(dataWind$dateTime)


#Water temp at depth of DO sensor
# dataSensorTemp <- read.table(dataIn[2],header=T,sep='\t',colClasses=c(dateTime="POSIXct"))
dataSensorTemp <- read.csv("./Data/Model_Input/2015_18/FCR_2015_18_wvwa_sensorTemp.csv")
dataSensorTemp$dateTime <- as.POSIXct(dataSensorTemp$dateTime)



#Read in either the original DO and thermister data if the ModelVersion[5]=1 
#Or read in the combed versions from CombedDeep file
if(ModelVersion[5]==0|ModelVersion[5]==3|ModelVersion[5]==4){
  #DO
  # dataDO <- read.table(dataIn[1],header=T,sep='\t',colClasses=c(dateTime="POSIXct"))
  dataDO <- read.csv("./Data/Model_Input/2015_18/FCR_2015_18_wvwa_DO.csv")
  dataDO$dateTime <- as.POSIXct(dataDO$dateTime)
  
  
  #Temp profile
  # dataTempProfile <- read.table(dataIn[3],header=T,sep='\t',colClasses=c(dateTime="POSIXct"))
  dataTempProfile <- read.csv("./Data/Model_Input/2015_18/FCR_2015_18_TempProfiles_hobos_and_interped2018.csv")
  dataTempProfile$dateTime <- as.POSIXct(dataTempProfile$dateTime)
  
  
}else if(ModelVersion[5]==1|ModelVersion[5]==2){
  setwd(dirCombedData)
  #Read combed DO data
  dataDO <- read.table(dataInXtra[3],header=T,sep='\t',colClasses=c(dateTime="POSIXct"))
  
  #Read in Combed Temp profile
  dataTempProfile <- read.table(dataInXtra[4],header=T,sep='\t',colClasses=c(DateTime="POSIXct"))
  names(dataTempProfile)[1]<-"dateTime"
  
}else{}


##
#Display some info about time grain of measurements

#Print first five time readings for each variable
firstFiveTimes <- data.frame(DO=dataDO$dateTime[1:5],PAR=dataPAR$dateTime[1:5],windSpeed=dataWind$dateTime[1:5],sensorTemp=dataSensorTemp$dateTime[1:5],tempProfile=dataTempProfile$dateTime[1:5])
print('First five time readings for each variable'); print(firstFiveTimes)

#Calculate first differences of time readings, display unique values
difTimesDO <- diff(dataDO$dateTime); print(table(difTimesDO))
#these lines were to identify weird gaps that are now fixed in lines 110 - 117
 # which(difTimesDO != 15)
 # test <- as.numeric(difTimesDO)
 # test <- test[test != 10]
 # print(test[!is.na(test)])
difTimesPAR <- diff(dataPAR$dateTime); print(table(difTimesPAR))
difTimesWindSpeed <- diff(dataWind$dateTime); print(table(difTimesWindSpeed))
difTimesSensorTemp <- diff(dataSensorTemp$dateTime); print(table(difTimesSensorTemp))
difTimesTempProfile <- diff(dataTempProfile$dateTime); print(table(difTimesTempProfile))


##
#Remove rows with duplicate dateTime stamps (and warn)

#Function to find duplicate dateTime stamps
findNotDupRows <- function(dataInName)
{
  #This function returns the indexes of rows where the dateTime is NOT a duplicate
  #of the dateTime in a previous row
  #dataInName is character, e.g. "dataPAR"
  dataIn <- eval(parse(text=dataInName))
  #Find duplicated time stamps
  dups <- duplicated(dataIn$dateTime)
  #If no duplicated time stamps, notDupRows=all rows in dataIn
  if (all(dups==FALSE))
  {
    notDupRows <- c(1:dim(dataIn)[1])
  } else
    #If at least one time stamp is duplicated, warn, and notDupRows is indexes
    #of rows where dateTime is not duplicate of the dateTime in a previous row
  {
    notDupRows <- which(dups==FALSE)
    nDups <- dim(dataIn)[1]-length(notDupRows)
    print(paste("Warning:",nDups,"rows with duplicate time stamps in",dataInName,"will be removed"))
  }
  #Return dupRows
  return(notDupRows)
}

notDupRows <- findNotDupRows("dataDO")
dataDO <- dataDO[notDupRows,]

notDupRows <- findNotDupRows("dataPAR")
dataPAR <- dataPAR[notDupRows,]

notDupRows <- findNotDupRows("dataWind")
dataWind <- dataWind[notDupRows,]

notDupRows <- findNotDupRows("dataSensorTemp")
dataSensorTemp <- dataSensorTemp[notDupRows,]

notDupRows <- findNotDupRows("dataTempProfile")
dataTempProfile <- dataTempProfile[notDupRows,]


##
#Make all data sets extend from startTime to endTime by timeStep
#Note that for some lakes it may be necessary to aggregate some variables to coarser time scale to get match up

#Round all time down to nearest timeStep (e.g. if timeStep is 5, round 00:07 to 00:05)
floorMins <- function(dataIn)
{
  #Pull out dateTime column and name it x
  x <- dataIn$dateTime
  nRows <- length(x)
  #Truncate each dateTime to hour; convert to class numeric
  floorHour <- as.POSIXct(trunc(x[1:nRows],"hour"))
  floorNumeric <- as.numeric(floorHour)
  #Create sequence from floorNumeric to next hour by timeStep (in seconds)
  seqSec <- seq(0,3600,60*timeStep)
  #Create matrix where each row is floorNumeric + the elements of seqSec
  matSec <- matrix(rep(seqSec,nRows),nrow=nRows,byrow=T)
  matSec <- floorNumeric + matSec
  #Calculate abs(time difference) between each element of x and the timeStep intervals
  difs <- abs(as.numeric(x) - matSec)
  #Find the minimum absolute difference in each row and select the corresponding time from matSec
  whichMin <- apply(difs,1,which.min)
  rowNames <- as.numeric(rownames(data.frame(whichMin)))
  matIndex <- (whichMin-1)*nRows + rowNames
  matSecFlat <- matrix(matSec,ncol=1)
  outTime <- as.POSIXct(matSecFlat[matIndex],origin="1970-01-01")
  #Return outTime
  return(outTime)
}

dataDO$dateTime <- floorMins(dataDO)
dataPAR$dateTime <- floorMins(dataPAR)
dataWind$dateTime <- floorMins(dataWind)
dataSensorTemp$dateTime <- floorMins(dataSensorTemp)
dataTempProfile$dateTime <- floorMins(dataTempProfile)

#Find the latest first time point and the earliest last time point of all the data
startTime <- max(min(dataDO$dateTime),min(dataPAR$dateTime),min(dataWind$dateTime),min(dataSensorTemp$dateTime),min(dataTempProfile$dateTime))
endTime <- min(max(dataDO$dateTime),max(dataPAR$dateTime),max(dataWind$dateTime),max(dataSensorTemp$dateTime),max(dataTempProfile$dateTime))

#Data.frame with one column "dateTime" which is sequence of times at time interval of timeStep, from startTrim to endTrim
completeTimes <- data.frame(dateTime=seq(startTime,endTime,paste(timeStep,"mins")))

#Merge all of input data.frames with completeTimes, so that they now all extend from startTime to endTime by timeStep
dataDO <- merge(completeTimes,dataDO,by="dateTime",all.x=T)
dataPAR <- merge(completeTimes,dataPAR,by="dateTime",all.x=T)
dataWind <- merge(completeTimes,dataWind,by="dateTime",all.x=T)
dataSensorTemp <- merge(completeTimes,dataSensorTemp,by="dateTime",all.x=T)
dataTempProfile <- merge(completeTimes,dataTempProfile,by="dateTime",all.x=T)

########################################
#Calculate sunrise, sunset

#Days of year for which to calculate sunrise and sunset
daysVec <- seq.POSIXt(trunc(startTime,"day"),trunc(endTime,"day"),"1 day")
#Day of year
day <- as.numeric(format(daysVec,format="%j"))


##New code here for sunrise/sunset####
#Get the sunrise and sunset times  
SunriseSunsetTimes<-getSunlightTimes(date=as.Date(daysVec),keep=c("sunrise","sunset"),lat=lat,lon=long,tz="EST")  

#Create data.frame with sunrise, sunset times for each day
sun <- data.frame(day=daysVec, sunrise=SunriseSunsetTimes$sunrise, sunset=SunriseSunsetTimes$sunset)

###################################################################################
##DCR: here we are expanding out the day to include points before and after sunrises
sun$sunrise<-sun$sunrise-(60*60*hours.BeforeAfter)
sun$sunset<-sun$sunset+(60*60*hours.BeforeAfter) 

########################################
#Trim data sets so that they end at last time before last sunrise
# i.e. lop off partial day at end

#Trim
endTrim <- max(sun$sunrise)
dataDO <- dataDO[dataDO$dateTime < endTrim,]
dataPAR <- dataPAR[dataPAR$dateTime < endTrim,]
dataWind<- dataWind[dataWind$dateTime < endTrim,]
dataSensorTemp <- dataSensorTemp[dataSensorTemp$dateTime < endTrim,]
dataTempProfile <- dataTempProfile[dataTempProfile$dateTime < endTrim,]
completeTimes <- data.frame(dateTime=completeTimes[completeTimes$dateTime < endTrim,])

#(Useful later) Vector giving which solar day each time in completeTimes belongs to
solarDaysBreaks <- sun$sunrise[sun$sunrise <= endTrim]
solarDaysVec <- cut.POSIXt(completeTimes$dateTime,breaks=solarDaysBreaks)


########################################
#Fill gaps in data

##
#DO - fill gaps up to 1 hours long (60 minutes)
dataDO <- fillHoles(dataDO,maxLength=60,timeStep=timeStep)
##Export filled data DO for trouble-shooting
#write.csv(dataDO,"DOfilled.csv")

##
#PAR - linearly interpolate gaps up to 60 min long
dataPAR <- fillHoles(dataPAR,maxLength=60,timeStep=timeStep)

##
#sensorTemp - linearly interpolate gaps up to 60 min long
dataSensorTemp <- fillHoles(dataSensorTemp,maxLength=366,timeStep=timeStep)

##
#windSpeed - fill with daily average as long as at least 80% of data are available

#Loop over days
for (i in 1:length(unique(solarDaysVec))-1){
  
  #Extract data between sunrise on day i and sunrise on day i+1
  timeSlice <- c(sun$sunrise[i], sun$sunrise[i+1]+(2*60*60*hours.BeforeAfter)) 
  #timeSlice <- c(sun$sunrise[i], sun$sunrise[i+1])  #DWH changed line above since we're not using 26 hour days
  dataTemp <- dataWind[dataWind$dateTime>=timeSlice[1] & dataWind$dateTime<timeSlice[2],]
  
  #Determine total number of observations, and number that are NA
  nTot <- length(dataTemp$windSpeed)
  nNA <- length(which(is.na(dataTemp$windSpeed)))
  
  #If >20% of obs are NA, skip to next i
  if (nNA/nTot > 0.20) next else
    
  {
    #Calculate mean windSpeed and sub in for NA values
    meanSpeed <- mean(dataTemp$windSpeed,na.rm=T)
    naRows <- as.numeric(row.names(dataTemp[is.na(dataTemp$windSpeed),]))
    dataWind$windSpeed[naRows] <- meanSpeed
  }
}

##
#tempProfile - linearly interpolate gaps up to 60 min long 

nCols <- dim(dataTempProfile)[2]

#Loop over the columns of dataTempProfile
for (i in 2:nCols){
  dataTemp <- dataTempProfile[,c(1,i)]
  dataTempFilled <- fillHoles(dataTemp,maxLength=366,timeStep=timeStep)
  dataTempProfile[,i] <- dataTempFilled[,2]
}

########################################
#Calculate zMix and fluxDummy

#Calc zMix

##Calculate the density profile for each column
##Create new file of densities with the same column headers as the tempProfile
dataDensProfile<-dataTempProfile

##Loop through the columns
for (j in 2:ncol(dataDensProfile)){
  dataDensProfile[,j]=1000*(1 - (dataTempProfile[,j]+288.9414)/(508929.2*(dataTempProfile[,j]+68.12963))*(dataTempProfile[,j]-3.9863)^2)
  ##Calculation for converting temperature to density:
  ##rho = 1000(1 - (T+288.9414)/(508929.2*(T+68.12963))*(T-3.9863)^2)
  # rho_t_24 = 1000*(24 - (24+288.9414)/(508929.2*(24+68.12963))*(24-3.9863)^2)
  
  ##End of column for loop
}


dataZMix <- calcZMixDens(dataDensProfile)

#Plot zMix
setwd(dirDump)
maxDepth <- max(as.numeric(substr(colnames(dataTempProfile)[2:nCols],5,10)))
pdf(file=paste(outName,'zMix.pdf'))
plot(zMix~dateTime,data=dataZMix,ylim=c(maxDepth,0))
dev.off()

#Identify when to shut off atmospheric flux
#If zMix > sensorDepth, then sensor is in mixed layer and fluxDummy = 1
#If zMix <= sensorDepth, then there is stratification at or above sensor and fluxDummy = 0 -- shut off atmosphere at this time step
fluxDummy <- as.numeric(dataZMix$zMix>sensorDepth)

flux_zmix_dataframe <- cbind(dataZMix, fluxDummy)

#Will = 0 if the date is before ice out in the spring
##as.numeric(as.Date(dataZMix$dateTime)>iceOut)

#will = 0 if the date is after ice On in the fall/early winter
##as.numeric(as.Date(dataZMix$dateTime)<iceOn)

###adding in ice on/off vectors 
#getting 2015 ice cover dates
head(ice)
head(flux_zmix_dataframe)
#chaging flux dummy to 0 for days w/ ice cover
fluxtest <- flux_zmix_dataframe %>% 
  mutate(Date = as.Date(dateTime))

fluxtest <- left_join(fluxtest, ice, by = "Date")
# fluxtest <- fluxtest %>% 
#   mutate(fluxdummynew = ifelse(Date %in% icedates, 0, fluxDummy))



fluxtest <- fluxtest %>% 
  mutate(Date = as.Date(dateTime)) %>% 
  mutate(icevector = ifelse(is.na(icevector), 0, icevector)) %>% 
  mutate(fluxdummynew = ifelse(icevector == 1, 0, fluxDummy))

write.csv(fluxtest, "./zmix_fluxdummy_2015_18.csv")

fluxDummy <- fluxtest$fluxdummynew



#If the ModelVersion[4]=0 then shut off atm flux during ice off periods
#Else leave flux dummy as is
# if(ModelVersion[4]==0){  #modify in here for our ice
#   #Use ice off dates to determine when to shut off flux dummy because of ice
#   #Shut off flux between 02Dec2007 and 24Apr2008 
#   fluxDummy<-fluxDummy*as.numeric(as.Date(dataZMix$dateTime)<iceOn.vector[1]|as.Date(dataZMix$dateTime)>iceOut.vector[2])
#   
#   #after 01Dec2008 ice on
#   fluxDummy<-fluxDummy*as.numeric(as.Date(dataZMix$dateTime)<iceOn.vector[2])
#   
# } else {#Do not change flux dummy
# }

########################################
#Merge data for convenience

#Merge
data1 <- dplyr::left_join(dataDO,dataPAR,by="dateTime")
data1 <- dplyr::left_join(data1,dataWind,by="dateTime")
data1 <- dplyr::left_join(data1,dataSensorTemp,by="dateTime")
data1 <- dplyr::left_join(data1,dataZMix,by="dateTime")

########################################
#Report on lengths of NA strings in data

#For each variable in data1, find the length of each run of NA in the data

#Have to temporarily sub in -99999 for NA to get rle() to work as desired
data1Temp <- as.matrix(data1[,2:6])
whichNA <- which(is.na(data1Temp))
data1Temp[whichNA] <- -99999
#DO
rleOut <- rle(data1Temp[,"DO"])
whichNA <- which(rleOut$values==-99999)
print('Lengths of DO NA strings')
print(sort(rleOut$lengths[whichNA]))
#PAR
rleOut <- rle(data1Temp[,"PAR"])
whichNA <- which(rleOut$values==-99999)
print('Lengths of PAR NA strings')
print(sort(rleOut$lengths[whichNA]))
#windSpeed
rleOut <- rle(data1Temp[,"windSpeed"])
whichNA <- which(rleOut$values==-99999)
print('Lengths of windSpeed NA strings')
print(sort(rleOut$lengths[whichNA]))
#sensorTemp
rleOut <- rle(data1Temp[,"sensorTemp"])
whichNA <- which(rleOut$values==-99999)
print('Lengths of sensorTemp NA strings')
print(sort(rleOut$lengths[whichNA]))
#zMix
rleOut <- rle(data1Temp[,"zMix"])
whichNA <- which(rleOut$values==-99999)
print('Lengths of zMix NA strings')
print(sort(rleOut$lengths[whichNA]))

rm(data1Temp)


########################################
#Calculate DOSat and kO2 at each time step

##
#Calculate average atmospheric pressure at elevation of lake
#Using the 'barometric formula' from Wikipedia - should double check
#Values of Rstar, g0, M are according to US Standard Atmosphere 1976; use ISO or SI instead?

#Constants
Pb <- 101325        #static pressure, pascals
Tb <- 288.15        #standard temp, K
Lb <- -0.0065       #standard temp lapse rate, K m-1
h <- elev           #elevation above sea level, m
hb <- 0             #elevation at bottom of atmospheric layer 0, m (note layer 0 extends to 11000 masl)
Rstar <-  8.31432   #universal gas constant, N m mol-1 K-1 (equiv to J K-1 mol-1)  SI: 8.314472
g0 <- 9.80665       #acceleration of gravity, m s-1
M <- 0.0289644      #molar mass of Earth's air, kg mol-1

#Pressure, in Pa (pascals)
P <- Pb * (Tb/(Tb+Lb*(h-hb)))^(g0*M/(Rstar*Lb))
# In mmHg
atmPres <- P*0.00750061683


##
#Calculate DO saturation
#Use eqn from Weiss 1970 Deep Sea Res. 17:721-735; simplified since salinity=0
# ln DO = A1 + A2 100/T + A3 ln T/100 + A4 T/100

attach(data1)

#Convert sensorTemp to Kelvin
sensorTempK <- sensorTemp + 273.15

#Weiss equation
A1 <- -173.4292;  A2 <- 249.6339;  A3 <- 143.3483;  A4 <- -21.8492
DOSat <- exp(((A1 + (A2*100/sensorTempK) + A3*log(sensorTempK/100) + A4*(sensorTempK/100))))

#Correction for local average atmospheric pressure
u <- 10^(8.10765 - (1750.286/(235+sensorTemp)))
DOSat <- (DOSat*((atmPres-u)/(760-u)))   #ml/L
DOSat <- DOSat/1000                      #L/L

#Convert using standard temperature and pressure. 
#Similar to calculating saturation DO at STP in ml/L, converting to mg?L (at STP),
#and then doing the above temperature and pressure conversions.
R <- 0.082057  #L atm deg-1 mol-1
O2molWt <- 15.999*2
convFactor <- O2molWt*(1/R)*(1/273.15)*(760/760) #g/L
DOSat <- DOSat*convFactor*1000                   #mg/L

##
#Calculate kO2
wp <- 0.15                       #exponent of wind profile power relationship, Smith 1985 Plant, Cell & Environment 8:387-398
#wind10 <- (10/windHeight)^wp * windSpeed
wind10 <- windSpeed #changing this to just equal wind speed since wind Height was normalized in published met data

k600 <- 2.07 + 0.215*wind10^1.7  #k600 in cm hr-1 per Cole and Caraco 1998;
k600 <- k600*24/100              #k600 in m day-1
schmidt <- 1800.6 - 120.1*sensorTemp + 3.7818*sensorTemp^2 - 0.047608*sensorTemp^3
kO2 <- k600*(schmidt/600)^-0.5   #Jahne et al. 87. exp could also be -.67
kO2 <- kO2*(timeStep/1440)       #change kO2 to units of m/(timeStep*min)

detach(data1)


########################################
#Fit model to find parameter estimates for each day

##
#Set up

#Organize input data
data2 <- data.frame(dateTime=data1$dateTime, DOObs=data1$DO, DOSat=DOSat, irr=data1$PAR, kO2=kO2, zMix=data1$zMix, fluxDummy=fluxDummy)


#Set up data.frame to store output of optimizations
nDays <- dim(sun)[1] - 1  #Not sure if this indexing works appropriately for all lakes
dateRange <- c(sun$day[1],sun$day[nDays])
outDays <- seq(dateRange[1],dateRange[2],"1 day")

#Modify this to include a column for the new variable: DO.offset
optimOut <- data.frame(solarDay=outDays,nll=rep(NA,nDays), iotaPEst=rep(NA,nDays), rhoEst=rep(NA,nDays), DOInitEst=rep(NA,nDays), DO.offset=rep(NA,nDays), Pmax=rep(NA,nDays),optimCode=rep(NA,nDays), R2=rep(NA,nDays), AIC=rep(NA,nDays))

#GPPFit calculated within the loop
GPPFitOut <- data.frame(solarDay=outDays,GPPFit=rep(NA,nDays))

#Calculate appropriate starting guesses for iota and rho parameters
#  This takes as reasonable guesses
#    iotaP = 0.5 (mg L-1 d-1)/(mmol m-2 s-1)
#    Pmax = 25 (mg L-1 d-1)/(mmol m-2 s-1)
#    rho  = 0.5 mg L-1 d-1
#  And converts them to appropriate units given time step of model

PmaxGuess <- 25*timeStep/1440
iotaPGuess <- 0.5*timeStep/1440
rhoGuess <- 0.5*timeStep/1440
DO.offsetGuess<-1.5

#Remove some stuff
rm(DOSat, kO2, fluxDummy)

#Change directory in prep for dumping results
# setwd(dirDump) #this is already set from above 



##
#Useful stuff for plotting

#Limits for y-axis for drivers
irrLims <- range(data2$irr,na.rm=T)
zMixLims <- c(max(data2$zMix,na.rm=T),0)
atmFluxLims <- c(-10*max(data2$kO2,na.rm=T)*min(data2$zMix,na.rm=T),10*max(data2$kO2,na.rm=T)*min(data2$zMix,na.rm=T))



#Set up pdf device
pdf(file=paste(outName,'daily fits.pdf'),width=11,height=8.5)
layout(rbind(matrix(c(1:14),nrow=2,byrow=F),matrix(c(15:28),nrow=2,byrow=F)),heights=c(1,1.5,1,1.5))
par(mar=c(1,2,0,0)+0.1)


##
#Run optimization for each day

for (i in 1:nDays){
  
  #Print occasional progress report on i
  if (i %in% seq(1,nDays,10)) (print(paste("Starting day",i)))
  
  #Extract data between sunrise on day i and sunrise on day i+1
  #Extends the temporary data for fitting to hours.BeforeAfter after the last sunrise
  timeSlice <- c(sun$sunrise[i], sun$sunrise[i+1]+(2*60*60*hours.BeforeAfter))
  #timeSlice <- c(sun$sunrise[i], sun$sunrise[i+1])
  dataTemp <- data2[data2$dateTime>=timeSlice[1] & data2$dateTime<timeSlice[2],]
  
  #If more than 20% of DOObs is missing, or if any NA in DOSat, irr, kO2, or zMix, 
  # return NA for optimization results and plot blank plots
  nTot <- length(dataTemp$DOObs)
  nNA <- length(which(is.na(dataTemp$DOObs)))
  if ((nNA/nTot > 0.20) |  any(is.na(dataTemp[,3:6])))
  {
    optimOut[i,2:7] <- NA
    frame(); frame()
    next
  } else
    
    #Otherwise, fit model and make plots
  {
    
    #For guess of initial DOHat, use first obs unless that is NA, in which case use min obs
    if (is.na(dataTemp$DOObs[1])==F) {(DOInit <- dataTemp$DOObs[1])} else {
      DOInit <- mean(dataTemp$DOObs[1:40],na.rm=T)
      #If the first value is false, we need a first DO value so set that equal to the mean of the first 40 values
      dataTemp$DOObs[1]<-mean(dataTemp$DOObs[1:40],na.rm=T)      
    }
    
    ##Code to excise the early morning hours depending on ModelVersion[6]
    ##To turn it off, set ModelVersion[6]=0
    
    #First find the dateTime of sunrise for this day
    #Or from the PAR values where you transition from 0s to >0 numbers
    PAR.sunrise<-dataTemp$dateTime[max(2,min(which(dataTemp$irr!=0,arr.ind=TRUE)))]
    
    
    #Take subset of data that is only the dateTime and Observed DO 
    keep.DF.temp<-c("dateTime","DOObs")
    DF.temp<-dataTemp[keep.DF.temp]
    
    #for all values from dataTemp$DOObs from sunrise to SR+ModelVersion[6]
    DF.temp$DOObs[DF.temp$dateTime>=PAR.sunrise & DF.temp$dateTime<(PAR.sunrise+60*60*ModelVersion[6])]<-NA
    
    
    ###Optional fill holes - this might work without filling the holes
    ###This will linearly interpolate the removed data
    
    #DF.temp<-fillHoles(DF.temp,maxLength=(ModelVersion[6]*60+10),timeStep=timeStep)
    
    #Update the dataTemp data frame with the excised data
    dataTemp$DOObs<-DF.temp$DOObs
    
    
    #Find parameter values by minimizing nll
    #parameter 1 is always iotaGuess, 2 is always rhoGuess, 3 is always DOoffset
    
    #4 will be DO initial if that is turned on, or Pmax if that is turned on
    
    #If 5 will be Pmax if photoinhibation and DO init is turned on
    
    
    parGuess <- log(c(iotaPGuess,rhoGuess))
    
    if(ModelVersion[5]==0){parGuess<-c(parGuess,log(c(DO.offsetGuess)))}
    
    if(ModelVersion[1]==1){parGuess<-c(parGuess,log(c(DOInit)))}
    
    if(ModelVersion[2]==1){parGuess<-c(parGuess,log(c(PmaxGuess)))}
    
    #This is the model fitting type
    op.Method<-"BFGS"
    
    #Lower and upper bounds for the L-BFGS-B method
    #lower.bound=log(c(0.000001,0.000001,0.000001))
    #upper.bound=log(c(50*timeStep/1440,50*timeStep/1440,5))
    
    optimTemp <- optim(parGuess,metabLoss_v8,dataIn=dataTemp,method=op.Method)
    
    #Save min nll
    optimOut[i,2] <- optimTemp$value
    #Save parameter estimates
    #  Multiply by 1440/timeStep to get from units of timeStep^-1 to units of day^-1
    #for iotaP and rho
    optimOut[i,3:4] <- exp(optimTemp$par[1:2])*(1440/timeStep)
    
    
    #Save estimate of parameter for the DO.offset parameter if it is estimated
    if(ModelVersion[5]==0){
      optimOut[i,6] <- exp(optimTemp$par[3])  
    } else {optimOut[i,6] <- NA}
    
    #Unpack Do initial, if it is estimated, store the estimate, otherwise, store NA
    if(ModelVersion[1]>1){
      optimOut[i,5] <- NA  
    } else if(ModelVersion[5]==0&ModelVersion[1]==1) {optimOut[i,5] <- exp(optimTemp$par[4])
    } else{optimOut[i,5] <- exp(optimTemp$par[3])}
    
    
    #Here is the Pmax parameter, multiple out to get from units of timeStep^-1 to units of day^-1 
    if(ModelVersion[2]==0){
      #Not fitting Photoinhibition
      optimOut[i,7]<-NA
      
    } else if(ModelVersion[5]==1|ModelVersion[5]==3|ModelVersion[5]==4 & ModelVersion[1]>1){
      #Fitting only Photoinhibition, not DO initial or DO offset
      optimOut[i,7] <- exp(optimTemp$par[3])*(1440/timeStep)
    } else if(ModelVersion[5]==0|ModelVersion[5]==3|ModelVersion[5]==4 & ModelVersion[1]>1){
      #Fitting Photoinhibition and DO initial
      optimOut[i,7] <- exp(optimTemp$par[4])*(1440/timeStep)
    } else if(ModelVersion[5]==1 & ModelVersion[1]==1){
      #Fitting only Photoinhibition, not DO initial
      optimOut[i,7] <- exp(optimTemp$par[4])*(1440/timeStep)
    } else {
      #Fitting only Photoinhibition, not DO initial
      optimOut[i,7] <- exp(optimTemp$par[5])*(1440/timeStep)
    }
    
    
    
    #Save code indicating whether nlm completed successfully
    #0 indicates completion
    #1 indicates the interation limit has been reached
    optimOut[i,8] <- optimTemp$convergence
    
    #Calculate atmFlux and DOHat given max likelihood parameter estimates
    predix <- metabPredix_v8(optimTemp$par,dataTemp)
    DOHat <- predix$DOHat
    atmFlux <- predix$atmFlux
    res <- predix$res
    
    #Calcualate SST, etc...
    obs <- na.exclude(dataTemp$DOObs) #drop DO NAs first so that observation and residuals line up
    dat <- data.frame(obs = obs, res = res)
    #dat <- na.exclude(dat) # remove rows w/ NA's
    SST <- sum(dat$obs^2, na.rm=T)
    SSE <- sum(dat$res^2, na.rm=T)
    chi <- sum(((dat$res^2)/dat$obs), na.rm=T)
    AIC <- chi - 2*length(dat$obs)
    SSM <- (SST-SSE)
    R2 <- (SST-SSE)/SST

    
    optimOut[i,9] <- R2
    optimOut[i,10] <- AIC	
    
    
    ##Calculate GPP in the correct units for each day and store it in GPPFitOut
    #GPP is in mg L-1 day-1
    
    #First calculate total mmol photons m-2 for each time step
    solarFlux <- dataTemp$irr*timeStep*60  #mmol m-2 timeStep-1
    
    #Convert iotaP and Pmax units
    iotaPNewUnits <- optimOut$iotaPEst[i]*10/(10*60*1440) #(mg L-1 time step-1) / (mmol m-2 time step-1)
    PmaxNewUnits <- optimOut$Pmax[i]*(10/1440) #(mg L-1 time step-1)
    
    #GPP for that time
    if(ModelVersion[2]==1){
      GPP_timePeriod<-PmaxNewUnits*(1-exp(-iotaPNewUnits*solarFlux/PmaxNewUnits))  
    } else {GPP_timePeriod<-iotaPNewUnits*solarFlux}
    
    
    GPPFitOut$GPPFit[i]<-sum(GPP_timePeriod)
    
    #############################################################
    #Debug and analysis
    #print(cbind(optimOut[i,],GPPFitOut$GPPFit[i],op.Method))
    #############################################################
    
    #Plot irradiance (orange points), zMix (dashed line), atmFlux (hollow black points)
    #y-axis tick labels are for atmFlux; positive values are flux into lake and negative values are flux out of lake
    par(mar=c(1,2,0,0)+0.1)
    plot(dataTemp$irr~dataTemp$dateTime, ylim=irrLims, axes=F, xlab="", ylab="", pch=18, col="dark orange")
    axis.POSIXct(1,dataTemp$dateTime,labels=F); box()
    text(x=min(dataTemp$dateTime),y=irrLims[2],labels=format(dataTemp$dateTime[1],format="%d-%b"),adj=c(0,1))
    par(new=T); plot(dataTemp$zMix~dataTemp$dateTime, ylim=zMixLims, type="l", lty=2, axes=F, xlab="", ylab="")
    par(new=T); plot(atmFlux~dataTemp$dateTime, axes=F, xlab="", ylab=""); axis(2)
    
    #Plot observed and predicted DO
    yLims <- range(c(DOHat,dataTemp$DOObs),na.rm=T)
    par(mar=c(2,2,0,0)+0.1)
    plot(DOHat ~ dataTemp$dateTime, ylim=yLims, type="l", axes=F, xlab="", ylab="")
    axis.POSIXct(1,dataTemp$dateTime,format="%H:%M")
    axis(2)
    box()
    points(dataTemp$DOObs ~ dataTemp$dateTime)
    meanDOSat <- round(mean(dataTemp$DOSat,na.rm=T),1)
    text(x=min(dataTemp$dateTime),y=yLims[2],labels=paste('DOSat',meanDOSat),adj=c(0,1))
    
    
    
    ###Export data for debugging
    #dataTemp is going to get you dateTime, DOObs, DOSat, irr, kO2, zmix, fluxdummy
    
    #Need wind
    dataTemp.windSpeed <- data1[data1$dateTime>=timeSlice[1] & data1$dateTime<timeSlice[2],]$windSpeed
    
    #Need DO predicted (DOHat), atmospheric flux (atmFlux)
    
    #This data frame will merge all the different columns together for a residual analysis
    DF.residuals.debug<-cbind(dataTemp,dataTemp.windSpeed,DOHat,atmFlux)
    
    #Calculate the actual residuals
    DF.residuals.debug$residuals<-DF.residuals.debug$DOObs-DF.residuals.debug$DOHat
    
    #Give the day for the entire column
    DF.residuals.debug$DAY<-optimOut$solarDay[i]
    
    #Merge with intialized dataFrame to keep track for the entire year
    DF.residuals.Analysis.Year<-rbind(DF.residuals.Analysis.Year,DF.residuals.debug)
    
  }
  
}  #end loop over nDays


#Close pdf graphics device
dev.off()

#Display optimOut
print(optimOut)

#Dump optimOut to dirDump
write.table(optimOut,paste(outName,'optimOut.txt'))

#Plot estimates of iota and rho
pdf(file=paste(outName,'time series of rho and iota.pdf'),width=8.5,height=11)
par(mfrow=c(2,1),mar=c(3,4,1,1))
plot(iotaPEst ~ solarDay, data=optimOut, xlab="", ylab="iotaPEst, (mg L-1 d-1) / (mmol m-2 s-1)")
plot(rhoEst ~ solarDay, data=optimOut, xlab="", ylab="rhoEst, (mg L-1 d-1)")
dev.off()


#######################################################
# Bookkeeping Calculations below, MCV 6/11/2009
attach(data2)

# calculate gasflux for each datapoint
gasflux <- kO2/timeStep *(DOSat - DOObs) / zMix  # grams O2 m^-3 minute^-1

# Find delta O2 / delta time
dDO<-diff(DOObs)
dT=diff(dateTime)
dDOdT <- dDO / as.numeric(dT,units="mins") # grams O2 m^-3 minute^-1

# calcuate "metabolism" value for each timestep, where m = R in dark, NEP in light
# since metabolism is calculated for a time interval (opposed to a point), use 
#   mean of zMix, and gasflux from start and end of each interval
m <- dDOdT - gasflux[1:length(gasflux)-1]+diff(gasflux)/2
# grams O2 m^-3 minute^-1


# For each night, find areal rate of respiration
R_nightly <- rep(NA,nDays)
for (i in 1:nDays)
{
  
  r <- which(dateTime>sun$sunset[i] & dateTime<sun$sunrise[i+1])
  R_nightly[i] <- mean(m[r],na.rm=T)*1440    # grams O2 m^-3 day^-1
  rm(r)
}

# For each daylight period, average R from night before and after
Rmean<-rep(NA,length(R_nightly)-1)
for (i in 1:length(R_nightly)-1){
  Rmean[i]<-(R_nightly[i]+R_nightly[i+1])/2 # grams O2 m^-3 day^-1
}

# Fill in R values for first and last days of deployment based on single nights
R <- c(R_nightly[1], Rmean)  # grams O2 m^-3 day^-1

# Calculate daylight NEP (this is NOT 24hr NEP)
Nd <- rep(NA,nDays)
for (i in 1:length(Nd)){
  r <- which(dateTime>sun$sunrise[i] & dateTime<sun$sunset[i])
  Nd[i] <- mean(m[r],na.rm=T) * (as.double(sun$sunset[i]-sun$sunrise[i],units="mins"))
  # grams O2 m^-3 daylight period^-1
  # line above multiplies the rate per minute by the number of daylight minutes
  rm(r)
}

# Calculate GPP & true NEP (24hr NEP)  # grams O2 m^-3 day^-1    
GPP <- Nd + (-R * as.double(sun$sunset[1:nDays]-sun$sunrise[1:nDays],units="days") )
# line above multiplies R (per day) by fraction of day that is daylight.
NEP <- GPP + R

#Redefine R so that positive respiration rates are displayed as positive
R <- -R

#Detach data2
detach(data2)

#Write out results of bookkeeping estimates
bookOut <- data.frame(solarDay=optimOut$solarDay,GPP,R)
write.table(bookOut,paste(outName,'bookOut.txt'))



########################################
#Plots to compare bookkeeping and fitting estimates

#Calculations for GPP fit are moved above into the for loop

#Write out GPPFit
#GPPFitOut <- data.frame(solarDay=optimOut$solarDay,GPPFit)
write.table(GPPFitOut,paste(outName,'GPPFitOut.txt'))

#Plot R vs rhoEst and GPP vs. GPPFit
pdf(file=paste(outName,'model comparison.pdf'),width=10,height=5)
par(mfrow=c(1,2))
plot(R ~ optimOut$rhoEst, xlab="R (max likelihood est)", ylab="R (bookkeeping est)"); abline(0,1)
plot(GPP ~ GPPFitOut$GPPFit, xlab="GPP (max likelihood est)", ylab="GPP (bookkeeping est)"); abline(0,1)
dev.off()


##Here I will output ultimate results in units that want (mg O2 L-1 day-1)

outputDCR<-optimOut
outputDCR<-merge(optimOut,GPPFitOut)

names(outputDCR)[grep("rhoEst", colnames(outputDCR))]<-"R_mgO2perLperday"
names(outputDCR)[grep("GPPFit", colnames(outputDCR))]<-"GPP_mgO2perLperday"
outputDCR$NEM_mgO2perLperday<-outputDCR$GPP_mgO2perLperday-outputDCR$R_mgO2perLperday
outputDCR$PtoR<-outputDCR$GPP_mgO2perLperday/outputDCR$R_mgO2perLperday

#Delete unwanted data frames
outputDCR$nll<-NULL
outputDCR$iotaPEst<-NULL
outputDCR$optimCode<-NULL
outputDCR$R2<-NULL
outputDCR$AIC<-NULL

##output as a tab delimited file
write.table(outputDCR,paste(outName,'outputDCR.txt'),sep="\t")

#Plot estimates of GPP and R
pdf(file=paste(outName,'Time series of GPP R NEM PtoR DCR.pdf'),width=8.5,height=11)
par(mfrow=c(2,1),mar=c(3,4,1,1))
plot(GPP_mgO2perLperday ~ solarDay, data=outputDCR, xlab="", ylab="GPP, (mg O2 L-1 d-1)")
plot(R_mgO2perLperday ~ solarDay, data=outputDCR, xlab="", ylab="R, (mg O2 L-1 d-1)")
plot(NEM_mgO2perLperday ~ solarDay, data=outputDCR, xlab="", ylab="NEM, (mg O2 L-1 d-1)")
plot(PtoR ~ solarDay, data=outputDCR, xlab="", ylab="P:R", ylim=c(0,2))

dev.off()


##END



