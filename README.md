# FCR_Metabolism
## Summary
Analysis of lake metabolism rates in Falling Creek Reservoir over multiple years 

## Recreating all analysis steps 

Run Install_packages.R and Download_data.R located in './Scripts/'. These scripts will download all need packages and data to rerun the entire analysis. 

Next run the Scripts in the folder './Scripts/01_Data_Comp/'in the following order 

- 1_CTD_YSI_Data_Compilation.R is used to compile temperature and light data from CTD and YSI profiles for the time frame before high frequency thermistors were deployed. Also used to compile CTD and YSI DO data for offset check between WVWA and EXO sensors.
- 2_NLDAS_Met_data_fill.R is used to develop regression between NLDAS and observed meteorology at FCR for metabolism model and environmental driver analysis 
- 3_EDI_Data_Comp.R is used to compile EXO DO data, thermistor, and meteorological data published on EDI
- 4_InsiteIG_DO_offset.R is used to compile and QAQC InsiteIG sensor data from EDI and develop offset between the InsiteIG and EXO DO sensors 
- 5_Seasons_Grouping_Testing.R is to develop seasons definitions based on thermal stratification thresholds calculated for FCR
- 6_Seasons_definitions_Ice_Vectors.Rmd is used to created csv files that contain season groupings and ice cover classification for the different season definitions outlined in Methods 2.5.2 and Supplementary Text S3. 
- 7_EnviDrivers_Comp.R is used to compile EDI data for the environmental driver analysis 
- E24_functions.R contains two functions that are called in script 7 to calculated E24

Next run the following two scripts located in './Scripts/02_Metab_Model/' to calculate and QAQC metabolism rates 

- 1a_metabFunction_FCR_EXO_18_22.R and 1b_metabFunction_FCR_wvwa_15_18.R are scripts that run the metabolism model for each time frame that the EXO or InsiteIG DO sensor was deployed. These scripts can be run in either order
- folder 'model_functions' contains four scripts: calcZMixDens.R, fillHoles.R, metabLoss_v8.R, metabPredix_v.R, that are functions used to run the metabFunction scripts 
- 2_Metab_Outputs_QAQC.R compiles and QAQCs metab model outputs to be used in data analysis scripts


Next run the following scripts located in './Scripts/03_Data_Analysis/' in order to QAQC metabolism data and run the analyses presented in the manuscript

- 1_Metab_Figures.Rmd is primary data analysis and visualization script 
- 2_Metab_Envi_Drivers.Rmd conducts environmental driver analysis

## Recreating figures or reruning metabolism models 

The workflow above is designed to be able to recreate the analysis at multiple entry points. All data products needed to run the metabolism model are created in './Scripts/01_Data_Comp/' and are stored in './Data/Model_Input/'. 

To just rerun the data analysis scripts, scripts 1-3 in './Scripts/03_Data_Analysis/' can be run, as all output data from the metabolism model are store in './Data/Model_Output/' along with a QAQCd csv file generated from script 2_Metab_Outputs_QAQC.R.


## Data folder 
This folders holds data that is generated from data compilation scripts and output data from metabolism models (Model_Output). This is also where downloaded EDI data will be held in './Data/EDI2023/' (which will be created once Scripts/Download_data.R is run). The folder Carey2022_glm holds modeled water temperature data and NLDAS-2 meteorological data from Carey et al. 2022 (add doi) that is used for generating interpolated water temperature profiles and developing offset of observed and modeled meteorology.


