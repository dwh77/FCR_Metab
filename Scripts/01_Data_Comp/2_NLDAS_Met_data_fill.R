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

#### Weibul correction of windspeed #### 

##pdf plots 
joinedSW |> 
  select(time, WindSpeed, Wind_hourly) |> 
  rename(WindSpeed_NLDAS = WindSpeed) |> 
  pivot_longer(-1) |> 
  ggplot()+
  geom_density(aes(value, fill = name))


##calculate weibul for local wind speed
#make data frame to hold weibul curve data
data <- data.frame(x = seq(0, 15, length.out = 100))

#remove NAs
localSW_noNA <- joinedSW |> filter(!is.na(Wind_hourly))
#estimate weibull paramaters 
local_params <- EnvStats::eweibull(localSW_noNA$Wind_hourly, method = "mle")
local_shape <- local_params$parameters[1]
local_scale <- local_params$parameters[2]

data$WD_local_fullYear <-
  dweibull(data$x, shape = local_shape, scale = local_scale)

##calculate weibul for NLDAS wind speed
#remove NAs
NLDASSW_noNA <- joinedSW |> filter(!is.na(WindSpeed))
#estimate weibull paramaters 
nldas_params <- EnvStats::eweibull(NLDASSW_noNA$WindSpeed, method = "mle")
nldas_shape <- nldas_params$parameters[1]
nldas_scale <- nldas_params$parameters[2]

data$WD_nldas_fullYear <-
  dweibull(data$x, shape = nldas_shape, scale = nldas_scale)

##weibul for leaf on 
joinedSW_leafon <- joinedSW |> 
  mutate(month = month(time)) |> 
  filter(month %in% c(5,6,7,8,9,10)) |> 
  filter(!is.na(Wind_hourly), !is.na(WindSpeed))

nldas_params_leafon <- EnvStats::eweibull(joinedSW_leafon$WindSpeed, method = "mle")
nldas_shape_leafon <- nldas_params_leafon$parameters[1]
nldas_scale_leafon <- nldas_params_leafon$parameters[2]

data$WD_nldas_leafon <-
  dweibull(data$x, shape = nldas_shape_leafon, scale = nldas_scale_leafon)


local_params_leafon <- EnvStats::eweibull(joinedSW_leafon$Wind_hourly, method = "mle")
local_shape_leafon <- local_params_leafon$parameters[1]
local_scale_leafon <- local_params_leafon$parameters[2]

data$WD_local_leafon <-
  dweibull(data$x, shape = local_shape_leafon, scale = local_scale_leafon)


#weibul for leaf off
joinedSW_leafoff <- joinedSW |> 
  mutate(month = month(time)) |> 
  filter(month %in% c(12,1,2,3))|> 
  filter(!is.na(Wind_hourly), !is.na(WindSpeed))

nldas_params_leafoff <- EnvStats::eweibull(joinedSW_leafoff$WindSpeed, method = "mle")
nldas_shape_leafoff <- nldas_params_leafoff$parameters[1]
nldas_scale_leafoff <- nldas_params_leafoff$parameters[2]

data$WD_nldas_leafoff <-
  dweibull(data$x, shape = nldas_shape_leafoff, scale = nldas_scale_leafoff)

local_params_leafoff <- EnvStats::eweibull(joinedSW_leafoff$Wind_hourly, method = "mle")
local_shape_leafoff <- local_params_leafoff$parameters[1]
local_scale_leafoff <- local_params_leafoff$parameters[2]

data$WD_local_leafoff <-
  dweibull(data$x, shape = local_shape_leafoff, scale = local_scale_leafoff)

data |> 
  pivot_longer(-1) |> 
  ggplot()+
  geom_line(aes(x = x, y = value, color = name), size = 1.5)+
  labs(x = "Wind Speed (m/s)", y = "Density", color = "Data and Grouping")+
  theme_classic()

#### corrected NLDAS data using weibul distributions 

#vectors to make data frames
percentiles <- seq(0.001, 0.999, by = 0.001)
x <- c(1:999)

#local weibul percentiles 
weibul_local_df <- data.frame("Percentile" = percentiles, "x_local" = x)

for(i in (weibul_local_df$Percentile)) {
   z <- EnvStats::eqweibull(localSW_noNA$Wind_hourly, p = i)
   y <- z$quantiles[1]
   
   weibul_local_df <- weibul_local_df |> 
     mutate(x_local = ifelse(Percentile == i, y, x_local))
 }


#nldas weibul percentiles
weibul_nldas_df <- data.frame("Percentile" = percentiles, "x_nldas" = x)

for(i in (weibul_nldas_df$Percentile)) {
  z <- EnvStats::eqweibull(NLDASSW_noNA$WindSpeed, p = i)
  y <- z$quantiles[1]
  
  weibul_nldas_df <- weibul_nldas_df |> 
    mutate(x_nldas = ifelse(Percentile == i, y, x_nldas))
  }


#bind weibul percentiles into one data frame
weibul_percentil_bind <- left_join(weibul_local_df, weibul_nldas_df, by = "Percentile")



### bind wind data to webiul percentiles

weibul_percentil_bind_local <- weibul_percentil_bind |>
  select(Percentile, x_local) |>
  mutate(x_local = round(x_local, digits = 3))

weibul_percentil_bind_nldas <- weibul_percentil_bind |>
  select(Percentile, x_nldas) |>
  mutate(x_nldas = round(x_nldas, digits = 3))

# 
# head(joinedSW)
# 
# local_wind <- joinedSW |> 
#   select(time, Wind_hourly) |> 
#   mutate(Wind_hourly = round(Wind_hourly, digits = 2))
# 
# local_wind_percent <- left_join(local_wind, weibul_percentil_bind_local, by = c("Wind_hourly" = "x_local")) |> 
#   filter(!is.na(Wind_hourly))

nldas_wind <- joinedSW |>
  select(time, WindSpeed) |>
  mutate(WindSpeed = round(WindSpeed, digits = 3))


##put nldas wind speed to percentile 
ldt <- data.table::data.table( nldas_wind , key = "WindSpeed" )
dt <- data.table::data.table( weibul_percentil_bind_nldas , key = "x_nldas" )

zz <- dt[ ldt , list(time, WindSpeed, Percentile , x_nldas  ) , roll = "nearest" ]


## now we have each nldas time lined up to a percentile, we'll join the local wind percentile data 
corrected_nldasWIND <- left_join(zz, weibul_percentil_bind_local, by = "Percentile")

corrected_nldasWIND |> 
  ggplot()+
  geom_density(aes(x_local))+ xlim(0,15)

nldas_wind_weibulcorrected <- corrected_nldasWIND |> 
  select(time, x_local) |> 
  rename(WindSpeed_nldas_weibulcor = x_local)

joinedSW <- left_join(joinedSW, nldas_wind_weibulcorrected, by = 'time')

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

ggsave(filename = "./Figures/Fig_S2_Shortwave_NLDAS_FCRmet.jpg",
       SW_nldas_fcrmet, device = "jpg", width = 180, height = 150, units = "mm")

#### Comparing NLDAS to FCR WindSpeed data ####

Wind_nldas_fcrmet <- joinedSW %>% 
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

Wind_nldasWeibul_fcrmet <- joinedSW %>% 
  select(time, Wind_hourly, WindSpeed_nldas_weibulcor) %>% 
  rename(FCRmet = Wind_hourly,
         NLDAS = WindSpeed_nldas_weibulcor) %>% 
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

Wind_nldasWeibul_fcrmet

ggsave(filename = "./Figures/Fig_S3_WindSpeed_NLDASweibul_FCRmet.jpg",
       Wind_nldasWeibul_fcrmet, device = "jpg", width = 180, height = 150, units = "mm")



#### Make input files for SW and WindSpeed based on NLDAS ####

sw_lm <- summary(lm(joinedSW$SW_hourly~joinedSW$ShortWave))
# wind_lm <- summary(lm(joinedSW$Wind_hourly~joinedSW$WindSpeed))


model_input <- joinedSW %>% 
  select(time, SW_hourly, Wind_hourly, ShortWave, WindSpeed_nldas_weibulcor) %>% 
  mutate(SW_regress = ( coefficients(sw_lm)[1] + (coefficients(sw_lm)[2]* ShortWave) ) ,
         #WindSpeed_regress = ( coefficients(wind_lm)[1] + (coefficients(wind_lm)[2] * WindSpeed) )
         ) %>% 
  mutate(PAR_from_SWregress = sw.to.par.base(SW_regress, coeff = 2.114)) %>% 
  rename(SW_hourly_metstation = SW_hourly,
         Wind_hourly_metstation = Wind_hourly,
         ShortWave_nldas = ShortWave)

#write csv
write.csv(model_input, "./Data/Generated_Data/Met_input_from_NLDAS.csv", row.names = F)


#check rmse and model fit for SW and Wind
model_input_tests <- model_input %>% 
  mutate(Diff_SW = SW_hourly_metstation - ShortWave_nldas
         #Diff_Wind = Wind_hourly_metstation - WindSpeed_nldas
         ) %>% 
  filter(!is.na(SW_hourly_metstation),
         !is.na(Wind_hourly_metstation))
  
plot(model_input_tests$time, model_input_tests$Diff_SW)

round(rmse(model_input_tests$SW_hourly_metstation, model_input_tests$ShortWave_nldas), digits = 3)


# plot(model_input_tests$time, model_input_tests$Diff_Wind)
# round(rmse(model_input_tests$Wind_hourly_metstation, model_input_tests$WindSpeed_nldas), digits = 3)




