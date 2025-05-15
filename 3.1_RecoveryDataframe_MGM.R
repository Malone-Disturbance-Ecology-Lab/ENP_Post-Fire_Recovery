# Recovery MASTER
# M.Grace McLeod (2023)

# This script extracts spectral, seasonal condition, and fire history data to Recovery sample points


rm(list=ls())

library(sf)
library(rgdal)
library(sp)
library(dplyr)
library(tidyverse)
library(raster)
library(ggplot2)
library(readr)
library(spatialEco)
library(tidyr)
library(MASS) 
library(reshape2) 
library(reshape) 
library(terra)
library(lubridate)
library(sqldf)
library(splitstackshape)


# load Recovery sample points
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")
Recov_smpl_pts <- rgdal::readOGR(dsn = "./Sampling", layer = "Recov_smpl_pts")
Recov_df <- as.data.frame(Recov_smpl_pts)


##########################################################################################################################################################
# SPECTRAL DATA 
##########################################################################################################################################################

# load Spec_Master
load("./Image_Processing/Master_Spec/Spec_Master.RDATA")

# SUBSET SPEC DATA TO JUST RECOVERY SAMPLE POINTS and SAMPLE PERIOD
# Format date to get year
Spec_Master$Obs_Date <- as.Date(Spec_Master$Obs_Date, format= "%Y%m%d")
Spec_Master$Obs_Year <- format(Spec_Master$Obs_Date, "%Y")
Spec_Master$Obs_month <- format(Spec_Master$Obs_Date, "%m")

# Keep only rows with point ID that matches Recov_Smpl_pts
Recov_Spec  <- filter(Spec_Master, ptID %in% Recov_df$ptID)

# BAND MIN AND MAX CALUCULATIONS (only due bands necessary for baseline)
Recov_Spec$Pt.B4max<- NA
Recov_Spec$Pt.B5max<- NA
Recov_Spec$Pt.B7max<- NA
Recov_Spec$Pt.B4min<- NA

# loop to calc min and max per band per point
for (i in unique(Recov_Spec$ptID)){
  print(i)
  # subset by point id
  SubID <- Recov_Spec[which(Recov_Spec$ptID == i),]
  # get min and max
  Recov_Spec$Pt.B4min[Recov_Spec$ptID ==i] <-min(SubID$B4) #B4
  Recov_Spec$Pt.B4max[Recov_Spec$ptID ==i] <-max(SubID$B4)
  Recov_Spec$Pt.B5max[Recov_Spec$ptID ==i] <-max(SubID$B5)
  Recov_Spec$Pt.B7max[Recov_Spec$ptID ==i] <-max(SubID$B7)
}

# Band Ratios
Recov_Spec$NIR.SWIR1 <- Recov_Spec$Pt.B4max / Recov_Spec$Pt.B5max
Recov_Spec$SWIR1.SWIR2 <- Recov_Spec$Pt.B5max / Recov_Spec$Pt.B7max

#save 
setwd("/Users/gracemcleod/Documents/EvergladesFire/Recovery")
save(Recov_Spec, file="Recov_Spec.RDATA")

rm(list=ls())

##########################################################################################################################################################
# SEASONAL CONDITION DATA
##########################################################################################################################################################

# load Recovery sample points
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")
Recov_smpl_pts <- rgdal::readOGR(dsn = "./Sampling", layer = "Recov_smpl_pts")
# make a copy for seasonal conditions
SeasonalCond_sp <- Recov_smpl_pts
# check crs
crs(SeasonalCond_sp)
SeasonalCond_sp <- spTransform(SeasonalCond_sp, crs( "+proj=lcc +lat_0=42.5 +lon_0=-100 +lat_1=25 +lat_2=60 +x_0=0 +y_0=0 +ellps=WGS84
+units=km +no_defs"))

# DAYMET DATES
# DAYMET layers are stacked in order but dates are not assigned. 
# list of months and years generated in Excel. Use to make df of dates. 
setwd("./Seasonal_Cond")
Date <-read.csv(file="DAYMET_dates.csv") # list of dates has been pre-generated
Date$Yr_Mo <- as.Date(paste0(as.character(Date$Date), '01'), format='%Y%m%d')  

# MIN TEMPERATURE.........................................................................................................................................
tmin_EVG <- stack("tmin_EVG.tif")
# assign layer names or z dimension specifying month and year to DAYMENT 
names(tmin_EVG) <- paste("tmin", Date$Yr_Mo, sep=".")
# extract data for each variable to each observation based on observation month
Recov_tmin_sp <- raster::extract(tmin_EVG, SeasonalCond_sp, method= "bilinear", buffer=NULL, df=TRUE, sp=TRUE, factors=TRUE)
Recov_tmin_df <- as.data.frame (Recov_tmin_sp)
# round to just 1 decimal
Recov_tmin_df <- Recov_tmin_df%>% mutate_at(vars(2:253), funs(round(., 1)))
# turn DF longways
Recov_tmin <- melt(Recov_tmin_df, na.rm=FALSE, value.name="tmin", id=c("ptID", "coords.x1", "coords.x2"))
# format columns
colnames(Recov_tmin)[5] <- "tmin" 
Recov_tmin$variable <- as.character(Recov_tmin$variable) 
Recov_tmin <- Recov_tmin %>%
  transform(variable=str_replace(variable,"tmin.",""))
Recov_tmin$variable <- as.Date(Recov_tmin$variable, format= "%Y.%m.%d")
Recov_tmin$Obs_Year <- format(Recov_tmin$variable, "%Y")
Recov_tmin$Obs_month <- format(Recov_tmin$variable, "%m")
Recov_tmin$variable <- NULL


# MAX TEMPERATURES.........................................................................................................................................
tmax_EVG <- stack("tmax_EVG.tif")
# assign layer names or z dimension specifying month and year to DAYMENT 
names(tmax_EVG) <- paste("tmax", Date$Yr_Mo, sep=".")
# extract data for each variable to each observation based on observation month
Recov_tmax_sp <- raster::extract(tmax_EVG, SeasonalCond_sp, method= "bilinear", buffer=NULL, df=TRUE, sp=TRUE, factors=TRUE)
Recov_tmax_df <- as.data.frame (Recov_tmax_sp)
# round to just 1 decimal
Recov_tmax_df <- Recov_tmax_df%>% mutate_at(vars(2:253), funs(round(., 1)))
# turn DF longways
Recov_tmax <- melt(Recov_tmax_df, na.rm=FALSE, value.name="tmax", id=c("ptID",  "coords.x1", "coords.x2"))
# format columns
colnames(Recov_tmax)[5] <- "tmax" 
Recov_tmax$variable <- as.character(Recov_tmax$variable) 
Recov_tmax <- Recov_tmax %>%
  transform(variable=str_replace(variable,"tmax.",""))
Recov_tmax$variable <- as.Date(Recov_tmax$variable, format= "%Y.%m.%d")
Recov_tmax$Obs_Year <- format(Recov_tmax$variable, "%Y")
Recov_tmax$Obs_month <- format(Recov_tmax$variable, "%m")
Recov_tmax$variable <- NULL


# PRECIPITATION.........................................................................................................................................
precip_EVG <- stack("precip_EVG.tif")
# assign layer names or z dimension specifying month and year to DAYMENT 
names(precip_EVG) <- paste("precip", Date$Yr_Mo, sep=".")
# extract data for each variable to each observation based on observation month
Recov_precip_sp <- raster::extract(precip_EVG, SeasonalCond_sp, method= "bilinear", buffer=NULL, df=TRUE, sp=TRUE, factors=TRUE)
Recov_precip_df <- as.data.frame (Recov_precip_sp)
# round to just 1 decimal
Recov_precip_df <- Recov_precip_df%>% mutate_at(vars(2:253), funs(round(., 1)))
# turn DF longways
Recov_precip <- melt(Recov_precip_df, na.rm=FALSE, value.name="precip", id=c("ptID",  "coords.x1", "coords.x2"))
# format columns
colnames(Recov_precip)[5] <- "precip" 
Recov_precip$variable <- as.character(Recov_precip$variable) 
Recov_precip <- Recov_precip %>%
  transform(variable=str_replace(variable,"precip.",""))
Recov_precip$variable <- as.Date(Recov_precip$variable, format= "%Y.%m.%d")
Recov_precip$Obs_Year <- format(Recov_precip$variable, "%Y")
Recov_precip$Obs_month <- format(Recov_precip$variable, "%m")
Recov_precip$variable <- NULL


# MERGE DAYMET DATAFRAMES.........................................................................................................................................
# use join() instead of merge() because it's much faster with large dfs
# could use any join type for this because number of observations are the same for all dfs. 
Recov_DAYMET = full_join(Recov_tmin, Recov_tmax, by=c("ptID",  "coords.x1", "coords.x2", "Obs_Year", "Obs_month"), keep=FALSE) 
Recov_DAYMET = full_join(Recov_DAYMET, Recov_precip, by=c("ptID",  "coords.x1", "coords.x2", "Obs_Year", "Obs_month"), keep=FALSE) 
# save DF
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery")
save(Recov_DAYMET, file="Recov_DAYMET.RDATA")
write.csv(Recov_DAYMET, file="Recov_DAYMET.csv")

rm(list=ls())

##########################################################################################################################################################
# FIRE HISTORY DATA
##########################################################################################################################################################

# FIRE-SPECIFIC INFO.......................................................................................................................................
# Extract from shapefiles
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")

# load RECOV.master.RDATA or RECOV.master.shp (See Fire_Data_Edits Script)
load("./Fire_History/RECOV.master.RDATA")
# load Recovery sample points
Recov_smpl_pts <- rgdal::readOGR(dsn = "./Sampling", layer = "Recov_smpl_pts")
# make copy and rename Recovery points
Recov_Fire_info <- Recov_smpl_pts

# check crs
crs(RECOV.master.sp)
crs(Recov_Fire_info)
Recov_Fire_info <- spTransform(Recov_Fire_info, crs("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))

# extract shapefile info for each point
TABLE <- sp::over(Recov_Fire_info, RECOV.master.sp)
Recov_Fire_info$FireName <- TABLE$FireName 
Recov_Fire_info$FireNumber <- TABLE$FireNumber
Recov_Fire_info$StartDate <- TABLE$StartDate  
Recov_Fire_info$EndDate <- TABLE$EndDate 
Recov_Fire_info$FireYear<- TABLE$Year
Recov_Fire_info$FireType <- TABLE$FireType
# view in df
Recov_Fire_info_df <- as.data.frame(Recov_Fire_info)
head(Recov_Fire_info_df)

# SAVE 
writeOGR(Recov_Fire_info, dsn="./Recovery", layer="Recov_Fire_Info", driver="ESRI Shapefile", overwrite=TRUE)
save(Recov_Fire_info, Recov_Fire_info_df, RECOV.master.sp,  file="Recov_Fire_Info.RDATA")


# FIRE HISTORY METRICS........................................................................................................................................

# load Recovery sample points
Recov_smpl_pts <- rgdal::readOGR(dsn = "./Sampling", layer = "Recov_smpl_pts")
# make a copy for fire history
Fire_Histry_sp <- Recov_smpl_pts

# TOTAL FIRES
# total number of fires experienced by the point prior to AND including the recovery fire
load("./Fire_History/FireHistory_df.RDATA")
names(FireHistory_df)
# Keep only rows with point ID that matches Recov_Smpl_pts
Fire_Histry_sp <- spTransform(Fire_Histry_sp, crs("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))
Fire_Histry_sp <- as.data.frame(Fire_Histry_sp)
summary(Fire_Histry_sp$coords.x1)
summary(FireHistory_df$coords.x1)
Recov_Frequency  <- filter(FireHistory_df, ptID %in% Fire_Histry_sp$ptID)
# make fire history dataframe
Recov_FireHistory <- Recov_Frequency %>% dplyr::select(ptID, coords.x1, coords.x2)
# calculate total fires
names(Recov_Frequency)
Recov_FireHistory$TotalFires <- rowSums(Recov_Frequency[ , c(3:32)]) # 1978-2007
summary(Recov_FireHistory)

# PREVIOUS INTERVAL 
# number of years since fire prior to the recovery fire
load("./Fire_History/FireYears_df.RDATA")
names(FireYears_df)
# Keep only rows with point ID that matches Recov_Smpl_pts
Recov_FireYear  <- filter(FireYears_df, ptID %in% Fire_Histry_sp$ptID)
# remove all fires/years after 2007
Recov_FireYear <- Recov_FireYear[1:32] 
# identify recovery fire year
names(Recov_FireYear)
Recov_FireYear$RecFY <- apply(Recov_FireYear[3:32], 1, max, na.rm=TRUE)
# identify penultimate fire year
# make rec fire year NA
Recov_FireNA <- Recov_FireYear
Recov_FireNA$EVG_.2001[Recov_FireNA$RecFY == 2001] <- NA
Recov_FireNA$EVG_.2002[Recov_FireNA$RecFY == 2002] <- NA
Recov_FireNA$EVG_.2003[Recov_FireNA$RecFY == 2003] <- NA
Recov_FireNA$EVG_.2004[Recov_FireNA$RecFY == 2004] <- NA
Recov_FireNA$EVG_.2005[Recov_FireNA$RecFY == 2005] <- NA
Recov_FireNA$EVG_.2006[Recov_FireNA$RecFY == 2006] <- NA
Recov_FireNA$EVG_.2007[Recov_FireNA$RecFY == 2007] <- NA
# Penultimate fire year is new max
names(Recov_FireNA)
Recov_FireYear$PenUltFY <- apply(Recov_FireNA[3:32], 1, max, na.rm=TRUE)
# check it
view(Recov_FireYear)
# add to Recov_FireHistory
Recov_FireHistory$RecFY <- Recov_FireYear$RecFY
Recov_FireHistory$PenUltFY <- Recov_FireYear$PenUltFY
Recov_FireHistory$Prev.Int <- Recov_FireHistory$RecFY - Recov_FireHistory$PenUltFY
# convert inf to NA 
Recov_FireHistory[sapply(Recov_FireHistory, is.infinite)] <- NA


# MERGE THE TWO DATAFRAMES..........................................................................................................................................
summary(Recov_Fire_info_df)
summary(Recov_FireHistory)
Recov_FireHistory <- merge(Recov_Fire_info_df, Recov_FireHistory, by=c("ptID", "coords.x1", "coords.x2"))
# check that RecFY and FireYear match
head(Recov_FireHistory)
# remove RecFY
Recov_FireHistory$RecFY <- NULL

# save
setwd("./Recovery")
save(Recov_FireHistory, file="Recov_FireHist.RDATA")
write.csv(Recov_FireHistory, file="Recov_FireHist.csv")

rm(list=ls())


##########################################################################################################################################################
# DATAFRAME INTEGRATION
##########################################################################################################################################################
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery")

# load all variable dataframes
load("./Recov_FireHist.RDATA")
load("./Recov_DAYMET.RDATA")
load("./Recov_Spec.RDATA")

# Merge DAYMET with spectral observations 
# note: X and Y projections do not currently match across dataframes but it's ok because point IDs will be used.
Recov_Spec$Obs_month <- as.numeric(Recov_Spec$Obs_month)
Recov_DAYMET$Obs_month <- as.numeric(Recov_DAYMET$Obs_month)
Recov_DAYMET$Obs_Year <- as.numeric(Recov_DAYMET$Obs_Year)
Recov_Master <- merge(Recov_Spec, Recov_DAYMET, by=c("ptID",  "Obs_Year", "Obs_month")) 
# remove coords from DAYMET
Recov_Master$coords.x1.y <- NULL; Recov_Master$coords.x2.y <- NULL 
Recov_Master$coords.x1 <- Recov_Master$coords.x1.x  ;  Recov_Master$coords.x1.x <-NULL
Recov_Master$coords.x2 <- Recov_Master$coords.x2.x  ;  Recov_Master$coords.x2.x <-NULL
# Merge FireHistory with Master
Recov_Master <- merge(Recov_Master, Recov_FireHistory, by=c("ptID"))
Recov_Master$coords.x1.y <- NULL; Recov_Master$coords.x2.y <- NULL
Recov_Master$coords.x1 <- Recov_Master$coords.x1.x  ;  Recov_Master$coords.x1.x <-NULL
Recov_Master$coords.x2 <- Recov_Master$coords.x2.x  ;  Recov_Master$coords.x2.x <-NULL
# correct number of pts?
length(unique(Recov_Master$ptID))


# Save DF
save(Recov_Master, file="Recov_Master.RDATA")




