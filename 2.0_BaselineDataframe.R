# Baseline Dataframe
# M.Grace McLeod (2023)

# This script merges
# 1. spectral ("BL_Spec.RDATA"), 
# 2. seasonal condition ("BL_DAYMET.RDATA"), 
# 3. fire history ("BL_FireHist.RDATA") data to Baseline sample points.
# 4. Integrates all data into a Master dataframe ("BL_Master_df.RDATA")

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
library(randomForest)
library(datasets)
library(caret)
library(sqldf)
library(splitstackshape)
#install.packages(“corrplot”)
library(corrplot)
library(cowplot)
library(viridis)
library(terra)

# Set working directory
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")

# load upland sample points
BL_smpl_pts <- rgdal::readOGR(dsn = "./Sampling", layer = "BL_smpl_pts")
BL_df <- as.data.frame(BL_smpl_pts)


##########################################################################################################################################################
# 1. SPECTRAL OBSERVATIONS 
##########################################################################################################################################################

# load spectral master dataframe
load("./Image_Processing/Master_Spec/Spec_Master.RDATA")

# SUBSET DATA TO JUST BASELINE SAMPLE POINTS and SAMPLE PERIOD.........................................................................................
# Format date to get year
Spec_Master$Obs_Date <- as.Date(Spec_Master$Obs_Date, format= "%Y%m%d")
Spec_Master$Obs_Year <- format(Spec_Master$Obs_Date, "%Y")
Spec_Master$Obs_month <- format(Spec_Master$Obs_Date, "%m")
# Subset to 2010-2020
Spec_Master$Obs_Year <- as.numeric(Spec_Master$Obs_Year)
BL_Spec <- Spec_Master[which(Spec_Master$Obs_Year >= 2010),]
# Keep only rows with point ID that matches BL_Smpl_pts
BL_Spec  <- filter(BL_Spec, ptID %in% BL_df$ptID)

# BAND MIN AND MAX CALUCULATIONS.........................................................................................................................
# Band Max per pt
BL_Spec$Pt.B1max<- NA
BL_Spec$Pt.B2max<- NA
BL_Spec$Pt.B3max<- NA
BL_Spec$Pt.B4max<- NA
BL_Spec$Pt.B5max<- NA
BL_Spec$Pt.B7max<- NA
# Band min per pt
BL_Spec$Pt.B1min<- NA
BL_Spec$Pt.B2min<- NA
BL_Spec$Pt.B3min<- NA
BL_Spec$Pt.B4min<- NA
BL_Spec$Pt.B5min<- NA
BL_Spec$Pt.B7min<- NA
# loop to calc min and max per band per point
for (i in unique(BL_Spec$ptID)){
  print(i)
  # subset by point id
  SubID <- BL_Spec[which(BL_Spec$ptID == i),]
  # get min and max
  BL_Spec$Pt.B1min[BL_Spec$ptID ==i] <-min(SubID$B1) #B1
  BL_Spec$Pt.B1max[BL_Spec$ptID ==i] <-max(SubID$B1)
  BL_Spec$Pt.B2min[BL_Spec$ptID ==i] <-min(SubID$B2) #B2
  BL_Spec$Pt.B2max[BL_Spec$ptID ==i] <-max(SubID$B2)
  BL_Spec$Pt.B3min[BL_Spec$ptID ==i] <-min(SubID$B3) #B3
  BL_Spec$Pt.B3max[BL_Spec$ptID ==i] <-max(SubID$B3)
  BL_Spec$Pt.B4min[BL_Spec$ptID ==i] <-min(SubID$B4) #B4
  BL_Spec$Pt.B4max[BL_Spec$ptID ==i] <-max(SubID$B4)
  BL_Spec$Pt.B5min[BL_Spec$ptID ==i] <-min(SubID$B5) #B5
  BL_Spec$Pt.B5max[BL_Spec$ptID ==i] <-max(SubID$B5)
  BL_Spec$Pt.B7min[BL_Spec$ptID ==i] <-min(SubID$B7) #B7
  BL_Spec$Pt.B7max[BL_Spec$ptID ==i] <-max(SubID$B7)
}


# BAND RATIOS...........................................................................................................................................
# Calculate band ratios from max values, and additional indicies for improving model fit
# landscape
BL_Spec$NIR.SWIR1 <-   BL_Spec$Pt.B4max /   BL_Spec$Pt.B5max
BL_Spec$SWIR1.SWIR2 <-   BL_Spec$Pt.B5max /   BL_Spec$Pt.B7max

# save DF
setwd("./Baseline")
save(BL_Spec, file="BL_Spec.RDATA")

##########################################################################################################################################################
# 2. SEASONAL CONDITION DATA
##########################################################################################################################################################
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")

# CREATE DATAFRAMES FOR EACH VARIABLE: EXTRACT DATA TO POINTS, MELT, JOIN.
# load DAYMET stacks (created in DAYMET script)
tmin_EVG <- stack("./Seasonal_Cond/tmin_EVG.tif")
tmax_EVG <- stack("./Seasonal_Cond/tmax_EVG.tif")
precip_EVG <- stack("./Seasonal_Cond/precip_EVG.tif")

# SEASONAL CONDITION SAMPLE POINTS
# make a copy of sample points for seasonal conditions
SeasonalCond_sp <- BL_smpl_pts
# check crs'
crs(SeasonalCond_sp)
crs(tmin_EVG)
SeasonalCond_sp <- spTransform(SeasonalCond_sp, crs( "+proj=lcc +lat_0=42.5 +lon_0=-100 +lat_1=25 +lat_2=60 +x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs"))

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
BL_tmin_sp <- raster::extract(tmin_EVG, SeasonalCond_sp, method= "bilinear", buffer=NULL, df=TRUE, sp=TRUE, factors=TRUE)
BL_tmin_df <- as.data.frame (BL_tmin_sp)
# round to just 1 decimal
BL_tmin_df <- BL_tmin_df%>% mutate_at(vars(3:255), funs(round(., 1)))
# turn DF longways
BL_tmin <- melt(BL_tmin_df, na.rm=FALSE, value.name="tmin", id=c("ptID", "coords.x1", "coords.x2"))
# format columns
colnames(BL_tmin)[5] <- "tmin" 
BL_tmin$variable <- as.character(BL_tmin$variable) 
BL_tmin <- BL_tmin %>%
  transform(variable=str_replace(variable,"tmin.",""))
BL_tmin$variable <- as.Date(BL_tmin$variable, format= "%Y.%m.%d")
BL_tmin$Obs_Year <- format(BL_tmin$variable, "%Y")
BL_tmin$Obs_month <- format(BL_tmin$variable, "%m")
BL_tmin$variable <- NULL

# MAX TEMPERATURES.........................................................................................................................................
tmax_EVG <- stack("tmax_EVG.tif")
# assign layer names or z dimension specifying month and year to DAYMENT 
names(tmax_EVG) <- paste("tmax", Date$Yr_Mo, sep=".")
# extract data for each variable to each observation based on observation month
BL_tmax_sp <- raster::extract(tmax_EVG, SeasonalCond_sp, method= "bilinear", buffer=NULL, df=TRUE, sp=TRUE, factors=TRUE)
BL_tmax_df <- as.data.frame (BL_tmax_sp)
# round to just 1 decimal
BL_tmax_df <- BL_tmax_df%>% mutate_at(vars(3:255), funs(round(., 1)))
# turn DF longways
BL_tmax <- melt(BL_tmax_df, na.rm=FALSE, value.name="tmax", id=c("ptID",  "coords.x1", "coords.x2"))
# format columns
colnames(BL_tmax)[5] <- "tmax" 
BL_tmax$variable <- as.character(BL_tmax$variable) 
BL_tmax <- BL_tmax %>%
  transform(variable=str_replace(variable,"tmax.",""))
BL_tmax$variable <- as.Date(BL_tmax$variable, format= "%Y.%m.%d")
BL_tmax$Obs_Year <- format(BL_tmax$variable, "%Y")
BL_tmax$Obs_month <- format(BL_tmax$variable, "%m")
BL_tmax$variable <- NULL


# PRECIPITATION ..........................................................................................................................................
precip_EVG <- stack("precip_EVG.tif")
# assign layer names or z dimension specifying month and year to DAYMENT 
names(precip_EVG) <- paste("precip", Date$Yr_Mo, sep=".")
# extract data for each variable to each observation based on observation month
BL_precip_sp <- raster::extract(precip_EVG, SeasonalCond_sp, method= "bilinear", buffer=NULL, df=TRUE, sp=TRUE, factors=TRUE)
BL_precip_df <- as.data.frame (BL_precip_sp)
# round to just 1 decimal
BL_precip_df <- BL_precip_df%>% mutate_at(vars(3:255), funs(round(., 1)))
# turn DF longways
BL_precip <- melt(BL_precip_df, na.rm=FALSE, value.name="precip", id=c("ptID", "coords.x1", "coords.x2"))
# format columns
colnames(BL_precip)[5] <- "precip" 
BL_precip$variable <- as.character(BL_precip$variable) 
BL_precip <- BL_precip %>%
  transform(variable=str_replace(variable,"precip.",""))
BL_precip$variable <- as.Date(BL_precip$variable, format= "%Y.%m.%d")
BL_precip$Obs_Year <- format(BL_precip$variable, "%Y")
BL_precip$Obs_month <- format(BL_precip$variable, "%m")
BL_precip$variable <- NULL


# MERGE DAYMET DATAFRAMES............................................................................................................................
# could use any join type for this because number of observations are the same for all dfs. 
BL_DAYMET = full_join(BL_tmin, BL_tmax, by=c("ptID",  "coords.x1", "coords.x2", "Obs_Year", "Obs_month"), keep=FALSE) 
BL_DAYMET = full_join(BL_DAYMET, BL_precip, by=c("ptID", "coords.x1", "coords.x2", "Obs_Year", "Obs_month"), keep=FALSE) 
# save DF
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline")
save(BL_DAYMET, file="BL_DAYMET.RDATA")

##########################################################################################################################################################
# 3. FIRE HISTORY DATA
##########################################################################################################################################################
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")


# FIRE HISTORY METRICS........................................................................................................................................

# load Recovery sample points
BL_smpl_pts <- st_read("./Sampling/BL_smpl_pts.shp")
# make a copy for fire history
Fire_Histry_sp <- BL_smpl_pts

# TOTAL FIRES
# total number of fires experienced by the point prior to AND including the recovery fire
load("./Fire_History/FireHistory_df.RDATA")
# Keep only rows with point ID that matches BL_Smpl_pts
Fire_Histry_sp <- st_transform(Fire_Histry_sp, crs = 32617)  # UTM Zone 17N (WGS84). "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"
Fire_Histry_sp <- as.data.frame(Fire_Histry_sp)
BL_Frequency  <- filter(FireHistory_df, ptID %in% Fire_Histry_sp$ptID)
# make fire history dataframe
BL_FireHistory <- BL_Frequency %>% dplyr::select(ptID, coords.x1, coords.x2)
# calculate total fires
names(BL_Frequency)
BL_FireHistory$TotalFires <- rowSums(BL_Frequency[ , c(3:45)]) # 1978-2020
summary(BL_FireHistory)

# PREVIOUS INTERVAL 
# number of years since fire prior to the recovery fire
load("./Fire_History/FireYears_df.RDATA")
names(FireYears_df)
# Keep only rows with point ID that matches Recov_Smpl_pts
BL_FireYear  <- filter(FireYears_df, ptID %in% Fire_Histry_sp$ptID)
# identify last fire year
names(BL_FireYear)
BL_FireYear$LFY <- apply(BL_FireYear[3:35], 1, max, na.rm=TRUE) # all prior to 2010
# identify penultimate fire year
# make last fire year NA...there must be a more efficient way to do this...but it works
BL_FireNA <- BL_FireYear
BL_FireNA$EVG_.1978[BL_FireNA$LFY == 1978] <- NA
BL_FireNA$EVG_.1979[BL_FireNA$LFY == 1979] <- NA
BL_FireNA$EVG_.1980[BL_FireNA$LFY == 1980] <- NA
BL_FireNA$EVG_.1981[BL_FireNA$LFY == 1981] <- NA
BL_FireNA$EVG_.1982[BL_FireNA$LFY == 1982] <- NA
BL_FireNA$EVG_.1983[BL_FireNA$LFY == 1983] <- NA
BL_FireNA$EVG_.1984[BL_FireNA$LFY == 1984] <- NA
BL_FireNA$EVG_.1985[BL_FireNA$LFY == 1985] <- NA
BL_FireNA$EVG_.1986[BL_FireNA$LFY == 1986] <- NA
BL_FireNA$EVG_.1987[BL_FireNA$LFY == 1987] <- NA
BL_FireNA$EVG_.1988[BL_FireNA$LFY == 1988] <- NA
BL_FireNA$EVG_.1989[BL_FireNA$LFY == 1989] <- NA
BL_FireNA$EVG_.1990[BL_FireNA$LFY == 1990] <- NA
BL_FireNA$EVG_.1991[BL_FireNA$LFY == 1991] <- NA
BL_FireNA$EVG_.1992[BL_FireNA$LFY == 1992] <- NA
BL_FireNA$EVG_.1993[BL_FireNA$LFY == 1993] <- NA
BL_FireNA$EVG_.1994[BL_FireNA$LFY == 1994] <- NA
BL_FireNA$EVG_.1995[BL_FireNA$LFY == 1995] <- NA
BL_FireNA$EVG_.1996[BL_FireNA$LFY == 1996] <- NA
BL_FireNA$EVG_.1997[BL_FireNA$LFY == 1997] <- NA
BL_FireNA$EVG_.1998[BL_FireNA$LFY == 1998] <- NA
BL_FireNA$EVG_.1999[BL_FireNA$LFY == 1999] <- NA
BL_FireNA$EVG_.2000[BL_FireNA$LFY == 2000] <- NA
BL_FireNA$EVG_.2001[BL_FireNA$LFY == 2001] <- NA
BL_FireNA$EVG_.2002[BL_FireNA$LFY == 2002] <- NA
BL_FireNA$EVG_.2003[BL_FireNA$LFY == 2003] <- NA
BL_FireNA$EVG_.2004[BL_FireNA$LFY == 2004] <- NA
BL_FireNA$EVG_.2005[BL_FireNA$LFY == 2005] <- NA
BL_FireNA$EVG_.2006[BL_FireNA$LFY == 2006] <- NA
BL_FireNA$EVG_.2007[BL_FireNA$LFY == 2007] <- NA
BL_FireNA$EVG_.2008[BL_FireNA$LFY == 2008] <- NA
BL_FireNA$EVG_.2009[BL_FireNA$LFY == 2009] <- NA
# Penultimate fire year is new max
names(BL_FireNA)
BL_FireYear$PenUltFY <- apply(BL_FireNA[3:35], 1, max, na.rm=TRUE)
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

rm(list=ls())


##########################################################################################################################################################
# 4. Baseline DATAFRAME INTEGRATION  
##########################################################################################################################################################
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")

# merge fire history and seasonal conditions with spectral observations to generate Baseline master dataframe
load("./Baseline/BL_DAYMET.RDATA")
load("./Baseline/BL_FireHist.RDATA")
load("./Baseline/BL_Spec.RDATA")

# Merge DAYMET with spectral observations 
# note: X and Y projections do not currently match across dataframes but it's ok because point IDs will be used.
BL_DAYMET$coords.x1 <- NULL ; BL_DAYMET$coords.x2 <- NULL # remove coords from DAYMET 
BL_Master_df <- merge(BL_Spec, BL_DAYMET, by=c("ptID", "Obs_Year", "Obs_month"))
# Merge FireHistory wiht Master
BL_FireHistory$coords.x1 <- NULL ; BL_FireHistory$coords.x2 <- NULL # remove coords from FireHistory
BL_Master_df <- merge(BL_Master_df, BL_FireHistory, by=c("ptID"))

# plot trends
BL_Master_df$TotalFires <- as.factor(BL_Master_df$TotalFires)
ggplot(BL_Master_df) +
  geom_jitter(aes(x=Obs_Year, y=NDVI, group=TotalFires, color=TotalFires))+
  geom_smooth(aes(x=Obs_Year, y=NDVI, group=TotalFires, color=TotalFires)) 

# Save DF
setwd("./Baseline")
save(BL_Master_df, file="BL_Master_df.RDATA")

