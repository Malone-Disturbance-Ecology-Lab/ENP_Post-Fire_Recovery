# This script creates the NDVI Model Dataframe:

# Workflow
# 1. Imports and formats the baseline points.
# 2. spectral ("BL_Spec.RDATA"), 
# 3. seasonal condition ("BL_DAYMET.RDATA"), 
# 4. fire history ("BL_FireHist.RDATA") data to Baseline sample points.
# 5. Integrates all data into a Master dataframe ("BL_Master_df.RDATA")

rm(list=ls())

library(sf)
library(dplyr)
library(tidyverse)
library(terra)
library(ggplot2)
library(readr)
library(spatialEco)

library(tidyr)
library(reshape2) 
library(reshape) 

library(lubridate)
library(datasets)
library(caret)
library(sqldf)
library(splitstackshape)


# Set working directory
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")

# load upland sample points
BL_smpl_pts <- sf::st_read(dsn = "/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Sampling/BL_smpl_pts.shp")  %>% st_transform(crs = 4326) # sf file

# Add lat and long into the file:
BL_smpl_pts$lat <- sf::st_coordinates(BL_smpl_pts)[,1]
BL_smpl_pts$lon <- sf::st_coordinates(BL_smpl_pts)[,2]

BL_df <- as.data.frame(BL_smpl_pts)

##########################################################################################################################################################
# 1. SPECTRAL OBSERVATIONS 
##########################################################################################################################################################

# load spectral master dataframe
load("./Image_Processing/Master_Spec/Spec_Master.RDATA")

# SUBSET DATA TO JUST BASELINE SAMPLE POINTS and SAMPLE PERIOD.........................................................................................
# Format date to get year

BL_Spec_org  <- Spec_Master %>% mutate(Obs_Date = as.Date(Obs_Date, format= "%Y%m%d"),
                                      Obs_Year = format(Obs_Date, "%Y") %>% as.numeric(),
                                      Obs_month = format(Obs_Date, "%m")) %>% filter(Obs_Year >= 2010) %>% filter( ptID %in% BL_df$ptID) %>% mutate(  B1= (B1* 0.0000275) + (-0.2), 
                                                                                                                                                      B2= (B2* 0.0000275) + (-0.2), 
                                                                                                                                                      B3= (B3* 0.0000275) + (-0.2),
                                                                                                                                                      B4= (B4* 0.0000275) + (-0.2), 
                                                                                                                                                      B5= (B5* 0.0000275) + (-0.2),
                                                                                                                                                      B7= (B7* 0.0000275) + (-0.2))


# BAND MIN AND MAX CALUCULATIONS.........................................................................................................
BL_Spec_org  %>% names

BL_Spec_Summary <- BL_Spec_org  %>% reframe(.by= ptID,
                                       Pt.B1max = max(B1, na.rm=T),
                                       Pt.B2max = max(B2, na.rm=T),
                                       Pt.B3max = max(B3, na.rm=T),
                                       Pt.B4max = max(B4, na.rm=T),
                                       Pt.B5max = max(B5, na.rm=T),
                                       Pt.B7max = max(B7, na.rm=T),
                                       
                                       Pt.B1min = min(B1, na.rm=T),
                                       Pt.B2min = min(B2, na.rm=T),
                                       Pt.B3min = min(B3, na.rm=T),
                                       Pt.B4min = min(B4, na.rm=T),
                                       Pt.B5min = min(B5, na.rm=T),
                                       Pt.B7min = min(B7, na.rm=T),
                                       NIR.SWIR1 = Pt.B4max/ Pt.B5max,
                                       SWIR1.SWIR2 = Pt.B5max / Pt.B7max)
                                       

# Merge min max infor back into the main file by ptID:
BL_Spec <- BL_Spec_org %>% full_join(BL_Spec_Summary, by='ptID')

# save
save(BL_Spec, file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/BL_Spec_08082025.RDATA")

##########################################################################################################################################################
# 2. SEASONAL CONDITION DATA
##########################################################################################################################################################
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")

# CREATE DATAFRAMES FOR EACH VARIABLE: EXTRACT DATA TO POINTS, MELT, JOIN.
# load DAYMET stacks (created in DAYMET script)
tmin_EVG <- terra::rast("./Seasonal_Cond/tmin_EVG.tif")
tmax_EVG <- terra::rast("./Seasonal_Cond/tmax_EVG.tif")
precip_EVG <- terra::rast("./Seasonal_Cond/precip_EVG.tif")

# SEASONAL CONDITION SAMPLE POINTS
# make a copy of sample points for seasonal conditions
SeasonalCond_sp <- BL_smpl_pts

SeasonalCond_sp <- st_transform(SeasonalCond_sp, "+proj=lcc +lat_0=42.5 +lon_0=-100 +lat_1=25 +lat_2=60 +x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs")

# DAYMET DATES
# DAYMET layers are stacked in order but dates are not assigned. 
# list of months and years generated in Excel. Use to make df of dates. 
setwd("./Seasonal_Cond")
Date <-read.csv(file="DAYMET_dates.csv") # list of dates has been pre-generated
Date$Yr_Mo <- as.Date(paste0(as.character(Date$Date), '01'), format='%Y%m%d')  


# MIN TEMPERATURE.........................................................................................................................................
tmin_EVG <- terra::rast("tmin_EVG.tif")

# assign layer names or z dimension specifying month and year to DAYMENT 
names(tmin_EVG) <- paste("tmin", Date$Yr_Mo, sep=".")

terra::time(tmin_EVG) <- Date$Yr_Mo

# extract data for each variable to each observation based on observation month
BL_tmin_df <- terra::extract(tmin_EVG, SeasonalCond_sp, xy=F, ID=F) %>% as.data.frame() %>% mutate_at(vars(1:252), funs(round(., 1))) %>% cbind(SeasonalCond_sp)

# turn DF longways
BL_tmin_long <- reshape::melt(BL_tmin_df, na.rm=FALSE, value.name= "tmin", id=c("ptID","geometry","lat","lon")) 

# format
BL_tmin <- BL_tmin_long %>% mutate( var = as.character(variable)) %>% 
  transform(variable=stringr::str_replace(var,"tmin.","")) %>% mutate(tmin = value,
                              date = as.Date(variable, format= "%Y-%m-%d"),
                              Obs_Year = format(date, "%Y"),
                              Obs_month = format(date, "%m")) %>% dplyr::select(ptID, lat, lon, Obs_Year, Obs_month, tmin)
                              

# MAX TEMPERATURES.........................................................................................................................................
tmax_EVG <- terra::rast("tmax_EVG.tif")
# assign layer names or z dimension specifying month and year to DAYMENT 
names(tmax_EVG ) <- paste("tmax", Date$Yr_Mo, sep=".")
terra::time(tmax_EVG ) <- Date$Yr_Mo

# extract data for each variable to each observation based on observation month
BL_tmax_df <- terra::extract(tmax_EVG, SeasonalCond_sp, xy=F, ID=F) %>% as.data.frame() %>% mutate_at(vars(1:252), funs(round(., 1))) %>% cbind(SeasonalCond_sp)

# turn DF longways
BL_tmax_long <- reshape::melt(BL_tmax_df, na.rm=FALSE, value.name= "tmax", id=c("ptID","geometry","lat","lon")) 



BL_tmax <- BL_tmax_long %>% mutate( var = as.character(variable)) %>% 
  transform(variable= stringr::str_replace(var,"tmax.","")) %>% mutate(tmax = value,
                                                                       date = as.Date(variable, format= "%Y-%m-%d"),
                                                                       Obs_Year = format(date, "%Y"),
                                                                       Obs_month = format(date, "%m")) %>% dplyr::select(ptID, lat, lon, Obs_Year, Obs_month, tmax)

# PRECIPITATION ..........................................................................................................................................
precip_EVG <- terra::rast("precip_EVG.tif")
# assign layer names or z dimension specifying month and year to DAYMENT 
names(precip_EVG) <- paste("precip", Date$Yr_Mo, sep=".")
terra::time(precip_EVG) <- Date$Yr_Mo

# extract data for each variable to each observation based on observation month
BL_precip_df <- terra::extract(precip_EVG, SeasonalCond_sp, xy=F, ID=F) %>% as.data.frame() %>% mutate_at(vars(1:252), funs(round(., 1))) %>% cbind(SeasonalCond_sp)

# turn DF longways
BL_precip_long <- reshape::melt(BL_precip_df, na.rm=FALSE, value.name= "precip", id=c("ptID","geometry","lat","lon")) 

# format
BL_precip <- BL_precip_long %>% mutate( var = as.character(variable)) %>% 
  transform(variable= stringr::str_replace(var,"precip.","")) %>% mutate(precip = value,
                                                                       date = as.Date(variable, format= "%Y-%m-%d"),
                                                                       Obs_Year = format(date, "%Y"),
                                                                       Obs_month = format(date, "%m")) %>% dplyr::select(ptID, lat, lon, Obs_Year, Obs_month, precip)

# MERGE DAYMET DATAFRAMES............................................................................................................................
# could use any join type for this because number of observations are the same for all dfs. 

BL_DAYMET = full_join(BL_tmin, BL_tmax, by=c("ptID",  "lat", "lon", "Obs_Year", "Obs_month"), keep=FALSE) %>% 
  full_join(BL_precip, by=c("ptID",  "lat", "lon", "Obs_Year", "Obs_month"), keep=FALSE) 
# save DF

save(BL_DAYMET, file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/BL_DAYMET_08082025.RDATA")

##########################################################################################################################################################
# 3. FIRE HISTORY DATA
##########################################################################################################################################################
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")

# FIRE HISTORY METRICS........................................................................................................................................
# make a copy for fire history
Fire_Histry_sp <- BL_smpl_pts

# TOTAL FIRES
# total number of fires experienced by the point prior to AND including the recovery fire
load("./Fire_History/FireHistory_df.RDATA")
# Keep only rows with point ID that matches BL_Smpl_pts
Fire_Histry_sp <- st_transform(Fire_Histry_sp, crs = 32617)  # UTM Zone 17N (WGS84). "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"

BL_Frequency  <- filter(FireHistory_df, ptID %in% Fire_Histry_sp$ptID)

# make fire history dataframe
BL_FireHistory <- BL_Frequency %>% dplyr::select(ptID)

BL_Frequency %>% names
# calculate total fires
BL_FireHistory$TotalFires <- rowSums(BL_Frequency[ , c(2:44)]) # 1978-2020
summary(BL_FireHistory)

# PREVIOUS INTERVAL 
# number of years since fire prior to the recovery fire
load("./Fire_History/FireYears_df.RDATA")
names(FireYears_df)
# Keep only rows with point ID that matches Recov_Smpl_pts
BL_FireYear  <- filter(FireYears_df, ptID %in% Fire_Histry_sp$ptID)
# identify last fire year
names(BL_FireYear)
BL_FireYear$LFY <- apply(BL_FireYear[2:34], 1, max, na.rm=TRUE) # all prior to 2010
# identify penultimate fire year
# make last fire year NA...there must be a more efficient way to do this...but it works
BL_FireNA <- BL_FireYear
BL_FireNA %>% names
BL_FireNA$year_1978[BL_FireNA$LFY == 1978] <- NA
BL_FireNA$year_1979[BL_FireNA$LFY == 1979] <- NA
BL_FireNA$year_1980[BL_FireNA$LFY == 1980] <- NA
BL_FireNA$year_1981[BL_FireNA$LFY == 1981] <- NA
BL_FireNA$year_1982[BL_FireNA$LFY == 1982] <- NA
BL_FireNA$year_1983[BL_FireNA$LFY == 1983] <- NA
BL_FireNA$year_1984[BL_FireNA$LFY == 1984] <- NA
BL_FireNA$year_1985[BL_FireNA$LFY == 1985] <- NA
BL_FireNA$year_1986[BL_FireNA$LFY == 1986] <- NA
BL_FireNA$year_1987[BL_FireNA$LFY == 1987] <- NA
BL_FireNA$year_1988[BL_FireNA$LFY == 1988] <- NA
BL_FireNA$year_1989[BL_FireNA$LFY == 1989] <- NA
BL_FireNA$year_1990[BL_FireNA$LFY == 1990] <- NA
BL_FireNA$year_1991[BL_FireNA$LFY == 1991] <- NA
BL_FireNA$year_1992[BL_FireNA$LFY == 1992] <- NA
BL_FireNA$year_1993[BL_FireNA$LFY == 1993] <- NA
BL_FireNA$year_1994[BL_FireNA$LFY == 1994] <- NA
BL_FireNA$year_1995[BL_FireNA$LFY == 1995] <- NA
BL_FireNA$year_1996[BL_FireNA$LFY == 1996] <- NA
BL_FireNA$year_1997[BL_FireNA$LFY == 1997] <- NA
BL_FireNA$year_1998[BL_FireNA$LFY == 1998] <- NA
BL_FireNA$year_1999[BL_FireNA$LFY == 1999] <- NA
BL_FireNA$year_2000[BL_FireNA$LFY == 2000] <- NA
BL_FireNA$year_2001[BL_FireNA$LFY == 2001] <- NA
BL_FireNA$year_2002[BL_FireNA$LFY == 2002] <- NA
BL_FireNA$year_2003[BL_FireNA$LFY == 2003] <- NA
BL_FireNA$year_2004[BL_FireNA$LFY == 2004] <- NA
BL_FireNA$year_2005[BL_FireNA$LFY == 2005] <- NA
BL_FireNA$year_2006[BL_FireNA$LFY == 2006] <- NA
BL_FireNA$year_2007[BL_FireNA$LFY == 2007] <- NA
BL_FireNA$year_2008[BL_FireNA$LFY == 2008] <- NA
BL_FireNA$year_2009[BL_FireNA$LFY == 2009] <- NA

# Penultimate fire year is new max
names(BL_FireNA)
BL_FireYear$PenUltFY <- apply(BL_FireNA[2:34], 1, max, na.rm=TRUE)

# add to BL_FireHistory
BL_FireHistory$LFY <- BL_FireYear$LFY
BL_FireHistory$PenUltFY <- BL_FireYear$PenUltFY
BL_FireHistory$Prev.Int <- BL_FireHistory$LFY - BL_FireHistory$PenUltFY
# convert inf to NA 
BL_FireHistory[sapply(BL_FireHistory, is.infinite)] <- NA

# save DF
setwd("./Baseline")
save(BL_FireHistory, file="BL_FireHist.RDATA")
write_csv(BL_FireHistory, file="BL_FireHist.csv")

##########################################################################################################################################################
# 4. Baseline DATAFRAME INTEGRATION  
##########################################################################################################################################################
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")

# merge fire history and seasonal conditions with spectral observations to generate Baseline master dataframe
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/BL_DAYMET_08082025.RDATA")
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/BL_FireHist.RDATA")
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/BL_Spec_08082025.RDATA")

# Merge DAYMET with spectral observations 
# note: X and Y projections do not currently match across dataframes but it's ok because point IDs will be used.

BL_DAYMET$Obs_Year <- as.numeric(BL_DAYMET$Obs_Year)
BL_DAYMET$Obs_month <- as.numeric(BL_DAYMET$Obs_month )

BL_Spec$Obs_Year <- as.numeric(BL_Spec$Obs_Year )
BL_Spec$Obs_month <- as.numeric(BL_Spec$Obs_month)

BL_Master_df1 <- BL_Spec %>% full_join(BL_DAYMET, by=c("ptID", "Obs_Year", "Obs_month"))


# Merge FireHistory wiht Master
BL_Master_df <- BL_Master_df1 %>% left_join( BL_FireHistory, by=c("ptID"))

# Additional processing of dataframe: ####

length(unique(BL_Master_df$ptID)) # 14,981 pts

# filter to only obs with TSF >= 7 years
BL_Master_df <- BL_Master_df %>% mutate(TSF = Obs_Year - LFY,
                                        Prev.Int = replace_na( Prev.Int, 80),
                                        PenUltFY = replace_na( PenUltFY, 1977),
                                        Obs_month = as.numeric(Obs_month),
                                        Obs_Year = as.numeric(Obs_Year)) %>% filter(TSF >= 7)

BL_Master_df %>% summary

# TEST FOR VARIABILITY OF NDVI
# look at the pts with low NDVI, is there a pattern in their variance?
summary(BL_Master_df$NDVI) 
# Min.        1st Qu.    Median      Mean     3rd Qu.      Max. 
# -0.001093  0.202053  0.252070  0.250824   0.297974   0.502491  
# we want points that are not experiencing other disturbances or anything else weird. "stable pts"

BL_NDVI_var <- BL_Master_df %>% 
  group_by(ptID) %>%
  summarise(NDVIvar=var(NDVI, na.rm = T), NDVImax=max(NDVI, na.rm = T))

summary(BL_NDVI_var$NDVIvar)
#  Min.        1st Qu.       Median     Mean      3rd Qu.        Max. 
# 0.0002287  0.0011790.  0.0017590   0.0020267   0.0026056   0.0140301
summary(BL_NDVI_var$NDVImax)
# Min.     1st Qu.   Median    Mean   3rd Qu.    Max. 
# 0.1318  0.2813    0.3309   0.3300   0.3786   0.5025 
BL_Master_df <- BL_Master_df %>% left_join( BL_NDVI_var, by="ptID")
BL_Master_df$NDVImax

# Check out NDVI distribution
lowNDVI <- BL_Master_df[which(BL_Master_df$NDVI < .2),]
length(unique(lowNDVI$ptID))  # 13,907 out of 14,981 pts. so almost everybody hits a low NDVI. Not necessarily a location issue then.
summary(lowNDVI$NDVIvar) 

# look at pts that never get high
neverHigh <-  BL_Master_df[which(BL_Master_df$NDVImax < .2),]
length(unique(neverHigh$ptID)) #  173

# filter for just the obs with NDVImax greater than .2
highMAX <- BL_Master_df %>% filter(NDVImax >= .2)
length(unique(highMAX$ptID)) # 14,808

# now filter to remove all low obs
HIGH <- highMAX %>% filter(NDVI >= .2) 
length(unique(HIGH$ptID)) # no loss of pts 

# Filter df by number of observations per point...........................................................................................
# for a generalizable model, want to use the BEST data so points that have observations from all times of year.
# start with perfection. How does that limit the variables? how closely does the distribution match the original dataset?
# Frequency Table: How many points of each EcoType burned how many times?
# Total  
HIGH$ptID <- as.factor(HIGH$ptID) 
ptObsTotal <- dplyr::count(HIGH,ptID)
summary(ptObsTotal$n)
# Min.    1st Qu.    Median   Mean    3rd Qu.     Max. 
# 1.00   34.00     50.00     52.38   68.00     150.00

# remove pts with less than top quartile of obs
ptObsTotal$n <- as.numeric(ptObsTotal$n)
# mean
ptObsTotalRM68 <- ptObsTotal[which(ptObsTotal$n >= 68),] 
BL_Total_sample <- HIGH[HIGH$ptID %in% ptObsTotalRM68$ptID, ]

# How many points are there?
length(unique(BL_Total_sample$ptID)) # 3,736
# How many observations are there?
length(BL_Total_sample$ptID) # 333,290

# save
BL_Master_filtered <- BL_Total_sample

# SUBSET TESTING AND TRAINING DATA
# use stratified sampling so that data has same distribution (real, train, test)
BL_Master_filtered.nona.pineland <- BL_Master_filtered %>% filter(EcoType == 'Pineland') %>% na.omit 

BL_Master_filtered.nona <- BL_Master_filtered %>% na.omit

BL_Master_filtered.nona$ptID %>% unique %>% length
# Save
save(BL_Master_filtered, file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/BL_Master_filtered_08082025.RDATA")

save(BL_Master_filtered.nona,
     BL_Master_filtered.nona.pineland, file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/BL_Master_df_08082025.RDATA")
