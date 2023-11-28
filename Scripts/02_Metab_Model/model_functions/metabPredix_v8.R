#Return DOHat, residuals, atmFlux, and NLL given par estimates
#CTS 18 Aug 2009
#Modified DCR 07Feb2013
#Includes an extra variable to estimate the drift of the DO sensor
#Modified DCR 17Jun2013
#Includes Photoinhibition and DO initial estimates

metabPredix_v8 <- function (parsIn,dataIn) {

##
#Inputs

#Unpack parameters and exponentiate to force positive
  iotaP <- exp(parsIn[1])
  rho <- exp(parsIn[2])
  
  
  #Save estimate of parameter for the DO.offset parameter if it is estimated
  if(ModelVersion[5]==0){
    DO.offset <- exp(optimTemp$par[3])  
  } else {DO.offset <- 0}

  
  #This is for the model fitting of initial DO
  #If the ModelVersion[3]=1 then we are fitting initial DO, otherwise, it is the first value or average of first 1 hour
  #If ModelVersion[5]==0 then we are fitting the DO offset and DO initial
  #If ModelVersion[5]>0 then we are just fitting DO initial
  if(ModelVersion[5]==0&ModelVersion[1]==1){
    #We are fitting DO initial
    DOInit <- exp(parsIn[4])
  } else if(ModelVersion[5]>0&ModelVersion[1]==1){
    DOInit <- exp(parsIn[3])
  } else if(ModelVersion[1]==2){
    #DOinitial is the first DO value of the day
    DOInit<-dataIn$DOObs[1]
  } else {
    #DO initial is the average of the first 1 hour 
    DOInit<-mean(dataIn$DOObs[1:6],na.rm=TRUE)
  }
  
  
  #Photoinhibition
  if(ModelVersion[2]==0){
    #No fitting photoihibition
    Pmax<-NA 
  } else if(ModelVersion[5]>0 & ModelVersion[1]>1){
    #Fitting only Photoinhibition, not DO initial
    Pmax <- exp(parsIn[3])
  } else if(ModelVersion[5]==0 & ModelVersion[1]>1){
    #Fitting only Photoinhibition and DO offset not DO initial
    Pmax <- exp(parsIn[4])
  } else if(ModelVersion[5]>0 & ModelVersion[1]==1){
    #Fitting only Photoinhibition and DO initial, not DO offset
    Pmax <- exp(parsIn[4])
  } else {
    #Fitting only Photoinhibition, DO initial, and DO offset
    Pmax <- exp(parsIn[5])
  }
  
  

#Label useful things
nObs <- dim(dataIn)[1]
kO2 <- dataIn$kO2
DOSat <- dataIn$DOSat
zMix <- dataIn$zMix
irr <- dataIn$irr
DOObs <- dataIn$DOObs
fluxDummy <- dataIn$fluxDummy

##
#Calculate predictions and residuals

#Set up output
DOHat <- rep(NA,nObs)
atmFlux <- rep(NA,nObs)  
#Initialize DOHat
DOHat[1] <- DOInit

  #Calculate atmFlux and predicted DO for each time point
  #Fluxes out of lake have negative sign
  if(ModelVersion[2]==1){
    #With Photoinhibition
    for (i in 1:(nObs-1)) {
      atmFlux[i] <- fluxDummy[i] * -kO2[i] * ((DOHat[i]-DO.offset) - DOSat[i]) / zMix[i]  
      DOHat[i+1] <- DOHat[i] + Pmax*(1-exp(-iotaP*irr[i]/Pmax)) - rho + atmFlux[i]
    }
  } else {
    for (i in 1:(nObs-1)) {
      atmFlux[i] <- fluxDummy[i] * -kO2[i] * ((DOHat[i]-DO.offset) - DOSat[i]) / zMix[i]  
      DOHat[i+1] <- DOHat[i] + iotaP*irr[i] - rho + atmFlux[i]
    }
    
  }

#Compare observed and predicted DO; calculate residuals and NLL
#Exclude from calculation any cases where DOObs=NA
if (any(is.na(DOObs)))
  {
  NAObs <- which(is.na(DOObs))
  res <- DOObs[-NAObs] - DOHat[-NAObs]
  } else
  
  {
  res <- DOObs - DOHat
  }

nRes <- length(res)
SSE <- sum(res^2)
sigma2 <- SSE/nRes
NLL <- 0.5*((SSE/sigma2) + nRes*log(2*pi*sigma2))

#Set up output structure
dataOut <- list(atmFlux=atmFlux,DOHat=DOHat,res=res,NLL=NLL)

#Return dataOut
return(dataOut)

}
