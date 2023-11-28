#### EDI data downloads 
dir.create('Data/EDI2023')


###CTD
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/200/13/27ceda6bc7fdec2e7d79a6e4fe16ffdf"
infile1 <- paste0(getwd(),"/Data/EDI2023/CTD_final_2013_2022.csv")
download.file(inUrl1,infile1)

###YSI
inUrl2  <- "https://pasta.lternet.edu/package/data/eml/edi/198/11/6e5a0344231de7fcebbe6dc2bed0a1c3"
infile2 <- paste0(getwd(),"/Data/EDI2023/YSI_PAR_profiles_2013-2022.csv")
download.file(inUrl2,infile2)


###Met station
inUrl3  <- "https://pasta.lternet.edu/package/data/eml/edi/389/7/02d36541de9088f2dd99d79dc3a7a853"
infile3 <- paste0(getwd(),"/Data/EDI2023/FCR_Met_final_2015_2022.csv")
download.file(inUrl3,infile3)


### EXO data
inUrl4 <- "https://pasta.lternet.edu/package/data/eml/edi/271/7/71e6b946b751aa1b966ab5653b01077f"
infile4 <- paste0(getwd(),"/Data/EDI2023/FCR_Catwalk_2018_2022.csv")
download.file(inUrl4,infile4)


### dowload WVWA DO sonde data from EDI
inUrl5  <- "https://pasta.lternet.edu/package/data/eml/edi/1357/1/25b7a46649a57545ff6c34d8d474464f"
infile5 <- paste0(getwd(),"/Data/EDI2023/FCR_DOsondes_2012_2018.csv")
download.file(inUrl5,infile5)

##Download FCR HOBOs from EDI
inUrl6  <- "https://pasta.lternet.edu/package/data/eml/edi/1357/1/461a43f5bc61a213365defb035962289"
infile6 <- paste0(getwd(),"/Data/EDI2023/FCR_hobos_2015_2018.csv")
download.file(inUrl6,infile6)

### Bathymetry
inUrl7  <- "https://pasta.lternet.edu/package/data/eml/edi/1254/1/f7fa2a06e1229ee75ea39eb586577184"
infile7 <- paste0(getwd(),"/Data/EDI2023/Reservoir_bathy_EDI.csv")
download.file(inUrl7,infile7)


### download weir data 
inUrl8  <- "https://pasta.lternet.edu/package/data/eml/edi/202/10/c065ff822e73c747f378efe47f5af12b"
infile8 <- paste0(getwd(),"/Data/EDI2023/FCR_weir_2022.csv")
download.file(inUrl8,infile8)


###download chem data 
inUrl9  <- "https://pasta.lternet.edu/package/data/eml/edi/199/11/509f39850b6f95628d10889d66885b76"
infile9 <- paste0(getwd(),"/Data/EDI2023/FCR_chem_2022.csv")
download.file(inUrl9,infile9)

### download secchi
inUrl10  <- "https://pasta.lternet.edu/package/data/eml/edi/198/11/81f396b3e910d3359907b7264e689052"
infile10 <- paste0(getwd(),"/Data/EDI2023/FCR_secchi_2022.csv")
download.file(inUrl10,infile10)


### ice data 
#https://portal.edirepository.org/nis/mapbrowse?packageid=edi.456.4
inUrl11  <- "https://pasta.lternet.edu/package/data/eml/edi/456/4/8454a757d0203bd8913d8adc8607f6c4"
infile11 <- paste0(getwd(),"/Data/EDI2023/Ice_Data_2013_2022.csv")
download.file(inUrl11,infile11)


### download filtered chla
# inUrl12  <- "https://pasta.lternet.edu/package/data/eml/edi/555/3/2f670c8316af634c76effdd205623912"
# infile12 <- paste0(getwd(),"/Data/EDI2023/manual_chlorophyll_2014_2022.csv")
# download.file(inUrl12,infile12)



