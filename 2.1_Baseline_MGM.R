# Baseline MODEL
# M.Grace McLeod (2023)

# This script 
# 1. Extracts spectral, fire history, and seasonal condition data to Baseline sample points
# 2. Runs Random Forest models to determine variable importance
# 3. Fits a Baseline model 

# Uses OG sampling design, but updated filtering and only for pinelands
# START AT "FORMAT DATA FOR RANDOM FORESTS" STEP
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")

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


# load upland sample points
BL_smpl_pts <- rgdal::readOGR(dsn = "./Sampling", layer = "BL_smpl_pts")
BL_df <- as.data.frame(BL_smpl_pts)


##########################################################################################################################################################
# SPECTRAL OBSERVATIONS 
##########################################################################################################################################################

# load Spec_Master
load("./Image_Processing/Master_Spec/Spec_Master.RDATA")
#Spec_Master <- read_csv("./Image_Processing/Master_Spec/Spec_Master.csv")

# SUBSET SPEC DATA TO JUST Baseline SAMPLE POINTS and SAMPLE PERIOD.........................................................................................
# Format date to get year
Spec_Master$Obs_Date <- as.Date(Spec_Master$Obs_Date, format= "%Y%m%d")
Spec_Master$Obs_Year <- format(Spec_Master$Obs_Date, "%Y")
Spec_Master$Obs_month <- format(Spec_Master$Obs_Date, "%m")
# Subset to 2010-2020
Spec_Master$Obs_Year <- as.numeric(Spec_Master$Obs_Year)
BL_Spec <- Spec_Master[which(Spec_Master$Obs_Year >= 2010),]
# Keep only rows with point ID that matches BL_Smpl_pts
BL_Spec  <- filter(BL_Spec, ptID %in% BL_df$ptID)
length(unique(BL_Spec$ptID)) # check it

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


# ADDITIONAL INDICIES ..........................................................................................................................................
# EVI = 2.5 * ((Band 4 – Band 3) / (Band 4 + 6 * Band 3 – 7.5 * Band 1 + 1))
BL_Spec$EVI <- 2.5 * ((BL_Spec$B4 - BL_Spec$B3) / (BL_Spec$B4 + 6 * BL_Spec$B3 - 7.5 * BL_Spec$B1 + 1))
# SAVI = ((Band 4 – Band 3) / (Band 4 + Band 3 + 0.5)) * (1.5)
BL_Spec$SAVI <- ((BL_Spec$B4 - BL_Spec$B3) / (BL_Spec$B4 + BL_Spec$B3 + 0.5)) * (1.5)


# save DF
setwd("./Baseline")
save(BL_Spec, file="BL_Spec.RDATA")
#write_scv(BL_Spec, "BL_Spec.csv")

rm(list=ls())

##########################################################################################################################################################
# SEASONAL CONDITION DATA
##########################################################################################################################################################
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")

# CREATE DATAFRAMES FOR EACH VARIABLE: EXTRACT DATA TO POINTS, MELT, JOIN.
# load DAYMET stacks (created in DAYMET script)
tmin_EVG <- stack("./Seasonal_Cond/tmin_EVG.tif")
tmax_EVG <- stack("./Seasonal_Cond/tmax_EVG.tif")
precip_EVG <- stack("./Seasonal_Cond/precip_EVG.tif")

# SEASONAL CONDITION SAMPLE POINTS
# load Baseline sample points
setwd("./Sampling")
BL_smpl_pts <- rgdal::readOGR(dsn = ".", layer = "BL_smpl_pts")
# make a copy for seasonal conditions
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
# use join() instead of merge() because it's much faster with large dfs
# could use any join type for this because number of observations are the same for all dfs. 
BL_DAYMET = full_join(BL_tmin, BL_tmax, by=c("ptID",  "coords.x1", "coords.x2", "Obs_Year", "Obs_month"), keep=FALSE) 
BL_DAYMET = full_join(BL_DAYMET, BL_precip, by=c("ptID", "coords.x1", "coords.x2", "Obs_Year", "Obs_month"), keep=FALSE) 
# save DF
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline")
save(BL_DAYMET, file="BL_DAYMET.RDATA")
#write.csv(BL_DAYMET, file="BL_DAYMET.csv")
rm(list=ls())

##########################################################################################################################################################
# FIRE HISTORY DATA
##########################################################################################################################################################
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")


# FIRE HISTORY METRICS........................................................................................................................................

# load Recovery sample points
BL_smpl_pts <- rgdal::readOGR(dsn = "./Sampling", layer = "BL_smpl_pts")
# make a copy for fire history
Fire_Histry_sp <- BL_smpl_pts

# TOTAL FIRES
# total number of fires experienced by the point prior to AND including the recovery fire
load("./Fire_History/FireHistory_df.RDATA")
names(FireHistory_df)
# Keep only rows with point ID that matches BL_Smpl_pts
Fire_Histry_sp <- spTransform(Fire_Histry_sp, crs("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))
Fire_Histry_sp <- as.data.frame(Fire_Histry_sp)
summary(Fire_Histry_sp$coords.x1)
summary(FireHistory_df$coords.x1)
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
# check it
view(BL_FireYear)
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











# FIRE HISTORY SAMPLE POINTS
# load Baseline sample points
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")
BL_smpl_pts <- rgdal::readOGR(dsn = "./Sampling", layer = "BL_smpl_pts")
BL_df <- as.data.frame(BL_smpl_pts)
# make a copy for seasonal conditions
FireHistory_sp <- BL_smpl_pts

# LOAD FIRE HISTORY DATA
# time since fire
tsf <- "./Fire_History/FH_layers/TSF_BL.tif"
TSF_BL <- stack(tsf)
# freq 1978-2009
freq <- "./Fire_History/FH_layers/Freq_1978_2009.tif"
Freq1978_2009 <- raster(freq)



# EXTRACT FIRE HISTORY DATA TO SAMPLE POINTS...................................................................................................................
# TSF stack 
# format layer names for TSF
yrs <- 2010:2020
yr <- as.Date(as.character(yrs), format = "%Y")
y <- year(yr)
names(TSF_BL) <- paste("tsf", y)

# extract raster data to BL sample points
FireHistory_sp <- raster::extract(Freq1978_2009, FireHistory_sp, method= "simple", buffer=NULL, df=TRUE, sp=TRUE, factors=TRUE)
FireHistory_sp <- raster::extract(TSF_BL,  FireHistory_sp, method= "simple", buffer=NULL, df=TRUE, sp=TRUE, factors=TRUE)
# convert to dataframe
BL_FireHistory <- as.data.frame(FireHistory_sp)

# format columns
colnames(BL_FireHistory)[3] <- "TotalFires" 
# melt TSF columns and create Obs_Year variable
BL_FireHistory <- melt(BL_FireHistory, na.rm=FALSE, value.name="TSF", id=c("ptID", "TotalFires", "coords.x1", "coords.x2"))
colnames(BL_FireHistory)[colnames(BL_FireHistory) == "value"] <- "TSF" 
BL_FireHistory$variable <- as.character(BL_FireHistory$variable) 
BL_FireHistory<- BL_FireHistory%>%
  transform(variable=str_replace(variable,"tsf.",""))
colnames(BL_FireHistory)[colnames(BL_FireHistory) == "variable"] <- "Obs_Year" 


# PREVIOUS INTERVAL
# load BL_FireYears
BL_FireYears <- read_csv("./Fire_History/FH_layers/BL_FireYears.csv")
# match points to BL_FireHistory
BL_FireYears <- BL_FireYears[BL_FireYears$ptID %in% BL_FireHistory$ptID,]
# merge
head(BL_FireHistory)
head(BL_FireYears)
BL_FireYears <- as.data.frame(BL_FireYears)
BL_FireYears[sapply(BL_FireYears, is.infinite)] <- NA
BL_FireHistory <- merge(BL_FireHistory, BL_FireYears, by=c("ptID"))
# calculate Prev.Int
BL_FireHistory$Prev.Int <- BL_FireHistory$LFYear - BL_FireHistory$PFYear

# save DF
setwd("./Baseline")
save(BL_FireHistory, file="BL_FireHist.RDATA")
#write_csv(BL_FireHistory, file="BL_FireHist.csv")

rm(list=ls())

##########################################################################################################################################################
# Baseline DATAFRAME INTEGRATION  
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
# check it
head(BL_Master_df)

# plot trends
BL_Master_df$TotalFires <- as.factor(BL_Master_df$TotalFires)
head(BL_Master_df)
length(unique(BL_Master_df$ptID[BL_Master_df$TotalFires == 8]))
ggplot(BL_Master_df) +
  geom_jitter(aes(x=Obs_Year, y=NDVI, group=TotalFires, color=TotalFires))+
  geom_smooth(aes(x=Obs_Year, y=NDVI, group=TotalFires, color=TotalFires)) 

# Save DF
setwd("./Baseline")
save(BL_Master_df, file="BL_Master_df.RDATA")


##########################################################################################################################################################
# FORMAT DATA FOR RANDOM FOREST
##########################################################################################################################################################
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")

# load BL_master_df
load("./Baseline/BL_Master_df.RDATA")
length(unique(BL_Master_df$ptID)) # 14,981 pts

# filter to only obs with TSF >= 7 years
BL_Master_df$TSF <- BL_Master_df$Obs_Year - BL_Master_df$LFY
BL_Master_df <- BL_Master_df[which(BL_Master_df$TSF >= 7),] 
length(unique(BL_Master_df$ptID)) # 14,981 pts (interestingly, same number of pts just more obs. so the points that dont burn foever doesnt change) 
summary(BL_Master_df$TSF)
# Set Prev.Int with NAs to 80
BL_Master_df$Prev.Int [is.na(BL_Master_df$Prev.Int)] =  80
# Make Obs_month numeric
BL_Master_df$Obs_month <- as.numeric(BL_Master_df$Obs_month)


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
BL_Master_df <- merge(BL_Master_df, BL_NDVI_var, by="ptID")

# Check out NDVI distribution
lowNDVI <- BL_Master_df[which(BL_Master_df$NDVI < .2),]
length(unique(lowNDVI$ptID))  # 13,907 out of 14,981 pts. so almost everybody hits a low NDVI. Not necessarily a location issue then.
summary(lowNDVI$NDVIvar) 
# look at pts that never get high
neverHigh <-  BL_Master_df[which(BL_Master_df$NDVImax < .2),]
length(unique(neverHigh$ptID)) #  173

# filter for just the obs with NDVImax greater than .2
highMAX <- BL_Master_df[which(BL_Master_df$NDVImax >= .2),]
length(unique(highMAX$ptID)) # 14,808
# now filter that to remove all low obs
HIGH <- highMAX[which(highMAX$NDVI >=.2),] # 200,000 fewer obs (23% of obs under .2)
length(unique(HIGH$ptID)) # no loss of pts 


# look at relationship between NDVI and other variables 
# Seasonality
ggplot(BL_Master_df, aes(x=Obs_month, y=NDVI)) +
  geom_smooth() +
  ggtitle("all obs")
ggplot(lowNDVI, aes(x=Obs_month, y=NDVI)) +
  geom_smooth() +
  ggtitle("Obs with NDVI < 0.2")
ggplot(HIGH, aes(x=Obs_month, y=NDVI)) +
  geom_jitter() +
  geom_smooth()+
  ggtitle("pts with max above .2")
ggplot(HIGH, aes(x=precip, y=NDVI)) +
  geom_jitter() +
  geom_smooth()+
  ggtitle("precip")
ggplot(HIGH, aes(x=tmax, y=NDVI)) +
  geom_jitter() +
  geom_smooth()+
  ggtitle("tmax")
ggplot(HIGH, aes(x=tmax, y=NDVI)) +
  geom_jitter() +
  geom_smooth()+
  ggtitle("tmin")
# fire history
HIGH$TotalFires <- as.factor(HIGH$TotalFires)
ggplot(HIGH, aes(x=TSF, y=NDVI)) +
  geom_jitter() +
  geom_smooth() +
  ggtitle("Time since fire")
ggplot(HIGH, aes(x=TotalFires, y=NDVI)) +
  geom_jitter() +
  ggtitle("Total Fires")
ggplot(HIGH, aes(x=MFRI, y=NDVI)) +
  geom_jitter() +
  ggtitle("MFRI")
# characterization spec
ggplot(HIGH, aes(x=Pt.B4max, y=NDVI)) +
  geom_jitter() +
  geom_smooth()+
  ggtitle("B4max")
ggplot(HIGH, aes(x=Pt.B4min, y=NDVI)) +
  geom_jitter() +
  geom_smooth()+
  ggtitle("B4min")
ggplot(HIGH, aes(x=NIR.SWIR1, y=NDVI)) +
  geom_jitter() +
  geom_smooth()+
  ggtitle("NIR.SWIR1")
ggplot(HIGH, aes(x=SWIR1.SWIR2, y=NDVI)) +
  geom_jitter() +
  geom_smooth()+
  ggtitle("SWIR1.SWIR2")


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


# CHECK DATA DISTRIBUTION ACROSS CONDITIONS
OG.TSF <- ggplot(BL_Master_df, aes(x=TSF)) +
  geom_density() +
  ggtitle("OG.TSF")
OG.TF <- ggplot(BL_Master_df, aes(x=TotalFires)) +
  geom_density()+
  ggtitle("OG.TF")
OG.precip <- ggplot(BL_Master_df, aes(x=precip)) +
  geom_density()+
  ggtitle("OG.precip")
OG.tmax <- ggplot(BL_Master_df, aes(x=tmax)) +
  geom_density()+
  ggtitle("OG.tmax")
OG.B4 <- ggplot(BL_Master_df, aes(x=B4)) +
  geom_density()+
  ggtitle("OG.B4")
OG.NIR.SWIR <- ggplot(BL_Master_df, aes(x=NIR.SWIR1)) +
  geom_density() +
  ggtitle("OG.NIR.SWIR")

new.TSF <- ggplot(BL_Total_sample, aes(x=TSF)) +
  geom_density() +
  ggtitle("new.TSF")
new.TF <- ggplot(BL_Total_sample, aes(x=TotalFires)) +
  geom_density()+
  ggtitle("new.TF")
new.precip <- ggplot(BL_Total_sample, aes(x=precip)) +
  geom_density()+
  ggtitle("new.precip")
new.tmax <- ggplot(BL_Total_sample, aes(x=tmax)) +
  geom_density()+
  ggtitle("new.tmax")
new.B4 <- ggplot(BL_Total_sample, aes(x=B4)) +
  geom_density()+
  ggtitle("new.B4")
new.NIR.SWIR <- ggplot(BL_Total_sample, aes(x=NIR.SWIR1)) +
  geom_density() +
  ggtitle("new.NIR.SWIR")
# plot all together
library(patchwork)
(OG.TSF + new.TSF) / (OG.TF + new.TF) /(OG.precip + new.precip) / (OG.tmax + new.tmax) / (OG.B4 + new.B4) / (OG.NIR.SWIR + new.NIR.SWIR)

# save
setwd("./Baseline")
BL_Master_filtered <- BL_Total_sample
save(BL_Master_filtered, file="BL_Master_filtered.RDATA")


# SUBSET TESTING AND TRAINING DATA
# use stratified sampling so that data has same distribution (real, train, test)
train <- stratified(BL_Master_filtered, c("Obs_month", "tmin", "tmax", "precip", "TotalFires", "Prev.Int"), .6)
test <- anti_join(BL_Master_filtered, train)
# save datasets to use for all models
save(train, test, file="BL_train_test.RDATA")


##########################################################################################################################################################
# RANDOM FOREST 
##########################################################################################################################################################
# Random Forests is a "black box" machine learning algorithm for classification (response is a factor) and regression (response is numeric). 
# advantageous because it avoids overfitting
# can deal with a large number of variables
# tells importance of each variable as a predictor of the response
# very few assumptions so widely applicable with out much data manipulation

# Some resources:
# How Random Forests work:  https://www.listendata.com/2014/11/random-forest-with-r.html
# Tutorial: https://www.r-bloggers.com/2021/04/random-forest-in-r/
# https://hackernoon.com/random-forest-regression-in-r-code-and-interpretation
# https://towardsdatascience.com/random-forests-an-ensemble-of-decision-trees-37a003084c6c

# load train and test data
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")
load("./Baseline/BL_train_test.RDATA")

# STRATIFY DATASET(S)
# start with a small sample that maintains a distribution across all the data. 
# fit the model. 
# Evaluate the model against the testing data.  
# then increase the sample size to see if your model improves with more data
train.1 <- stratified(train, c( "Obs_month", "tmin", "tmax", "precip", "TotalFires",  "Prev.Int"), .1)
train.25 <- stratified(train, c( "Obs_month", "tmin", "tmax", "precip", "TotalFires", "Prev.Int"), .25)

# NDVI ............................................................................................................................................................... 


# VARIABLE SELECTION / MODEL FITTING 

# STEP 1: Manipulating model parameters
# all variables, default settings
system.time(NDVI_rf_1.1  <- randomForest( NDVI ~   precip + tmax + tmin + Obs_month + 
                                            Prev.Int + TotalFires  + 
                                            NIR.SWIR1 + SWIR1.SWIR2 + 
                                            Pt.B1max + Pt.B2max + Pt.B3max + Pt.B4max + Pt.B5max + Pt.B7max +
                                            Pt.B1min + Pt.B2min + Pt.B3min + Pt.B4min + Pt.B5min + Pt.B7min,
                                        data= train.1,
                                        importance=TRUE,
                                        verbose=TRUE,
                                        predicted=TRUE,
                                        keep.inbag=TRUE)) # VE = 60.95 time= 40 sec
# all variables, change ntree
system.time(NDVI_rf_1.2  <- randomForest( NDVI ~   precip + tmax + tmin + Obs_month + 
                                         Prev.Int + TotalFires  + 
                                          NIR.SWIR1 + SWIR1.SWIR2 + 
                                          Pt.B1max + Pt.B2max + Pt.B3max + Pt.B4max + Pt.B5max + Pt.B7max +
                                          Pt.B1min + Pt.B2min + Pt.B3min + Pt.B4min + Pt.B5min + Pt.B7min,
                                        data= train.1,
                                        ntree= 1000,
                                        importance=TRUE,
                                        verbose=TRUE,
                                        predicted=TRUE,
                                        keep.inbag=TRUE)) # VE = 61.27 # time= 77 sec (does not improve VE, much slower, keep default)
# all variables, change mtry
system.time(NDVI_rf_1.3  <- randomForest( NDVI ~   precip + tmax + tmin + Obs_month + 
                                         Prev.Int + TotalFires  + 
                                          NIR.SWIR1 + SWIR1.SWIR2 + 
                                          Pt.B1max + Pt.B2max + Pt.B3max + Pt.B4max + Pt.B5max + Pt.B7max +
                                          Pt.B1min + Pt.B2min + Pt.B3min + Pt.B4min + Pt.B5min + Pt.B7min,
                                        data= train.1,
                                        ntree= 500,
                                        mtry= 5,
                                        importance=TRUE,
                                        verbose=TRUE,
                                        predicted=TRUE,
                                        keep.inbag=TRUE)) # VE = 60.97 (does not improve VE, keep default)
# all variables, change nodes
# nodes = 10: time= 23 sec, VE= 59.4
# nodes = 15: time= 20 sec, VE= 59.26
# nodes = 20: time= 18 sec, VE= 58.99
# set to 10
system.time(NDVI_rf_1.4  <- randomForest( NDVI ~   precip + tmax + tmin + Obs_month + 
                                            Prev.Int + TotalFires  + 
                                            NIR.SWIR1 + SWIR1.SWIR2 + 
                                            Pt.B1max + Pt.B2max + Pt.B3max + Pt.B4max + Pt.B5max + Pt.B7max +
                                            Pt.B1min + Pt.B2min + Pt.B3min + Pt.B4min + Pt.B5min + Pt.B7min,
                                          data= train.1,
                                          ntree= 500,
                                          mtry= 3,
                                          nodesize = 10,  
                                          importance=TRUE,
                                          verbose=TRUE,
                                          predicted=TRUE,
                                          keep.inbag=TRUE)) # VE = 61.49 
varImpPlot(NDVI_rf_1.4)


# STEP 2:
# To start removing variables, first look at corrilation and remove the highly correlated ones
# Corrilation plots
load("./Baseline/BL_Master_filtered.RDATA")
BL_Master_filtered$Obs_month <- as.numeric(BL_Master_filtered$Obs_month)
BL_Master_noNA <- na.omit(BL_Master_filtered)
head(BL_Master_filtered)
# all considered variables
M <-cor(BL_Master_filtered[, c("Pt.B1max" ,"Pt.B2max" , "Pt.B3max", "Pt.B4max" , "Pt.B5max", "Pt.B7max" , 
                         "Pt.B1min", "Pt.B2min" ,"Pt.B3min", "Pt.B4min", "Pt.B5min", "Pt.B7min",
                         "NIR.SWIR1", "SWIR1.SWIR2",
                         "Obs_month", "tmin", "tmax", "precip", 
                         "TotalFires", "Prev.Int",
                         "coords.x1", "coords.x2")]) 
corrplot(M, method="circle") 
corrplot(M, method="color")
corrplot(M, method="number")

# decide between B1-3max
# no B1max: VE= 59.78
# no B2max: VE= 59.75
# no B3max: VE= 59.8
# without any: VE= 60.08 (none very important)
# Remove all 3
system.time(NDVI_rf_2.1  <- randomForest( NDVI ~  precip + tmax + tmin + Obs_month + 
                                            Prev.Int + TotalFires  + 
                                            NIR.SWIR1 + SWIR1.SWIR2 + 
                                            Pt.B4max + Pt.B5max + Pt.B7max +
                                            Pt.B1min + Pt.B2min + Pt.B3min + Pt.B4min + Pt.B5min + Pt.B7min,
                                          data= train.1,
                                          ntree= 500,
                                          mtry= 3,
                                          nodesize = 10,  
                                          importance=TRUE,
                                          verbose=TRUE,
                                          predicted=TRUE,
                                          keep.inbag=TRUE)) # VE= 60.08
# decide between B5max and B7max
# no B5: VE= 60.49
# no B7: VE= 60.45
# without either: VE= 60.72
# does not matter for VE, B5 more important by a little
# remove B7max
system.time(NDVI_rf_2.2  <- randomForest(NDVI ~   precip + tmax + tmin + Obs_month + 
                                           Prev.Int + TotalFires  + 
                                           NIR.SWIR1 + SWIR1.SWIR2 + 
                                           Pt.B4max + Pt.B5max + 
                                           Pt.B1min + Pt.B2min + Pt.B3min + Pt.B4min + Pt.B5min + Pt.B7min,
                                          data= train.1,
                                          ntree= 500,
                                          mtry= 3,
                                          nodesize = 10,  
                                          importance=TRUE,
                                          verbose=TRUE,
                                          predicted=TRUE,
                                          keep.inbag=TRUE)) # VE= 61.3
# decide between B2-7min
# no B2min: VE= 60.45
# no B3min: VE= 60.44
# no B4min: VE= 60.11
# no B5min: VE= 60.6
# no B7min: VE= 60.54
# remove 2 and 7: VE= 60.78
# keep B4min only
system.time(NDVI_rf_2.3  <- randomForest(NDVI ~  precip + tmax + tmin + Obs_month + 
                                           Prev.Int + TotalFires  + 
                                           NIR.SWIR1 + SWIR1.SWIR2 + 
                                           Pt.B4max + Pt.B5max + 
                                           Pt.B4min,
                                          data= train.1,
                                          ntree= 500,
                                          mtry= 3,
                                          nodesize = 10,  
                                          importance=TRUE,
                                          verbose=TRUE,
                                          predicted=TRUE,
                                          keep.inbag=TRUE)) # VE= 60.28
# decide between tmin and tmax
# no tmin: VE=  60.23
# no tmax: VE= 59.43
# similar VE, tmax more important
# remove tmin
system.time(NDVI_rf_2.4  <- randomForest(NDVI ~ precip + tmax + Obs_month + 
                                           Prev.Int + TotalFires  + 
                                           NIR.SWIR1 + SWIR1.SWIR2 + 
                                           Pt.B4max + Pt.B5max + 
                                           Pt.B4min,
                                          data= train.1,
                                          ntree= 500,
                                          mtry= 3,
                                          nodesize = 10,  
                                          importance=TRUE,
                                          verbose=TRUE,
                                          predicted=TRUE,
                                          keep.inbag=TRUE)) # VE=  60.23
# re-do cor plot with remaining variables 
M <-cor(BL_Master_filtered[, c( "Pt.B4max" , "Pt.B5max",
                        "Pt.B4min", 
                         "NIR.SWIR1", "SWIR1.SWIR2",
                         "Obs_month",  "tmax", "precip", 
                         "TotalFires", "Prev.Int",
                         "coords.x1", "coords.x2")]) 
corrplot(M, method="circle") 
corrplot(M, method="color")
corrplot(M, method="number")
# does it imporve to decide between B4max, B5max and NIR:SWIR?
# NIR:SWIR way more important
# no B4max: VE= 59.38
# no B5max: VE= 60.33
# no NIR.SWIR: VE=  60.17
system.time(NDVI_rf_2.5  <- randomForest(NDVI ~ precip + tmax + Obs_month + 
                                           Prev.Int + TotalFires  + 
                                           NIR.SWIR1 + SWIR1.SWIR2 + 
                                           Pt.B4max + Pt.B5max + 
                                           Pt.B4min,
                                          data= train.1,
                                          ntree= 500,
                                          mtry= 3,
                                          nodesize = 10,  
                                          importance=TRUE,
                                          verbose=TRUE,
                                          predicted=TRUE,
                                          keep.inbag=TRUE)) # VE= 61.08
# NIR:SWIR and B4max are pretty close
# but both important. RF can handle it
# keep both. remove B5max
varImpPlot(NDVI_rf_2.4)

# STEP 3: 
# Anything never come up as important? rerun NDVI_rf_2.7 a couple times
# Does it imporve fit to remove it?
# no B4min: VE=  59.65
# no SWIR1:SWIR2: VE= 59.66
# no Obs_month: VE= 53.94 (keep!)
# keep all. 
system.time(NDVI_rf_3.1  <- randomForest(NDVI ~ precip + tmax + Obs_month + 
                                           Prev.Int + TotalFires  + 
                                           SWIR1.SWIR2 + 
                                           Pt.B4max + Pt.B5max + 
                                           Pt.B4min,
                                          data= train.1,
                                          ntree= 500,
                                          mtry= 3,
                                          nodesize = 10,  
                                          importance=TRUE,
                                          verbose=TRUE,
                                          predicted=TRUE,
                                          keep.inbag=TRUE)) # VE= no change to 2.5
# save models
setwd("./Baseline")
save(NDVI_rf_1.1, NDVI_rf_1.2, NDVI_rf_1.3, NDVI_rf_1.4,
     NDVI_rf_2.1, NDVI_rf_2.2, NDVI_rf_2.3, NDVI_rf_2.4, NDVI_rf_2.5,
     file="NDVI_rf_fitting.RDATA")

# STEP 4:
# try adding more vegetation levels 
# LOAD DATA
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")
# load veg data subset to uplands
load("./Veg_layers/Uplands_all.RDATA")
head(Uplands_df)
# load Baseline sample pts
BL_smpl_pts <- rgdal::readOGR(dsn = "./Sampling", layer = "BL_smpl_pts")
BL_df <- as.data.frame(BL_smpl_pts)
# load Baseline data
load("./Baseline/BL_train_test.RDATA")
# EXTRACT OTHER VEGETATION LEVELS
# make copy and rename recovery points
BL_veg_info <- BL_smpl_pts
# check crs
crs(Uplands)     # +proj=utm +zone=17 +datum=NAD83 +units=m +no_defs 
crs(BL_veg_info) # +proj=utm +zone=17 +datum=WGS84 +units=m +no_defs
Uplands <- spTransform(Uplands, "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")
# extract shapefile info for each point
TABLE <- sp::over(BL_veg_info, Uplands)
head(TABLE)
BL_veg_info$VegCode_L <- TABLE$VegCode_L 
BL_veg_info$VegCode_N  <- TABLE$VegCode_N 
BL_veg_info$L1_name  <- TABLE$L1_name   
BL_veg_info$L2_name  <- TABLE$L2_name  
BL_veg_info$L3_name <- TABLE$L3_name 
BL_veg_info$L4_name <- TABLE$L4_name 
BL_veg_info$L5_name <- TABLE$L5_name 
BL_veg_info$L6_name <- TABLE$L6_name
BL_veg_info$L7_name <- TABLE$L7_name 
# view in df
BL_veg_info_df <- as.data.frame(BL_veg_info)
head(BL_veg_info_df)
BL_veg_info_df$coords.x1 <- NULL ; BL_veg_info_df$coords.x2 <- NULL
# INTEGRATE WITH BL DATA
# train
class(train$ptID)
class(BL_veg_info_df$ptID)
BL_veg_info_df$ptID <- as.factor(BL_veg_info_df$ptID)
train <- merge(train, BL_veg_info_df, by="ptID")
train$L3_name <- as.factor(train$L3_name)
train$L4_name <- as.factor(train$L4_name)
train$L5_name <- as.factor(train$L5_name)
train$L6_name <- as.factor(train$L6_name)
head(train)
# test
test<- merge(test, BL_veg_info_df, by="ptID")
train$L7_name <- as.factor(train$L7_name)
test$L3_name <- as.factor(test$L3_name)
test$L4_name <- as.factor(test$L4_name)
test$L5_name <- as.factor(test$L5_name)
test$L6_name <- as.factor(test$L6_name)
test$L7_name <- as.factor(test$L7_name)
head(test)
# save datasets to use for all models
setwd("./Baseline")
save(train, test, file="BL_train_test.RDATA")

# all veg levels
train.1 <- stratified(train, c( "Obs_month", "tmin", "tmax", "precip", "TotalFires", "Prev.Int"), .1)
names(train.1)
system.time(NDVI_rf_4.1  <- randomForest(NDVI ~  precip + tmax + Obs_month + 
                                           Prev.Int + TotalFires  + 
                                           NIR.SWIR1 + SWIR1.SWIR2 + 
                                           Pt.B4max +  
                                           Pt.B3min + Pt.B4min +
                                           L1_name + L2_name + L3_name + L4_name + L5_name + L6_name + L7_name,
                                         data= train.1,
                                         ntree= 500,
                                         mtry= 3,
                                         nodesize = 10,  
                                         importance=TRUE,
                                         verbose=TRUE,
                                         predicted=TRUE,
                                         keep.inbag=TRUE)) # VE= 58.89
varImpPlot(NDVI_rf_4.1)
# not worth adding, doesn't really help



# FULL TRAINING DATA (same combo as 2.5)
system.time(NDVI_rf_5.1  <- randomForest(NDVI ~  precip + tmax + Obs_month + 
                                           Prev.Int + TotalFires  + 
                                           NIR.SWIR1 + SWIR1.SWIR2 + 
                                           Pt.B4max + Pt.B5max + 
                                           Pt.B4min,
                                          data= train,
                                          ntree= 500,
                                          mtry= 3,
                                          nodesize = 10,  
                                          importance=TRUE,
                                          verbose=TRUE,
                                          predicted=TRUE,
                                          keep.inbag=TRUE)) # VE=  71.21  time= 1:40hr
varImpPlot(NDVI_rf_5.1)
# apply to testing data
# prediction and confusion matricies
p1 <- predict(NDVI_rf_5.1, train)
summary(p1)
# TESTING DATA
test$NDVI_rf <- predict(NDVI_rf_5.1, test)
# see corrilation for prediction
summary(lm(test$NDVI ~ test$NDVI_rf))       # R2= 0.7287  

# save the model
#setwd("./Baseline")
NDVI_rf <- NDVI_rf_5.1
save(NDVI_rf, file="NDVI_rf.RDATA")



# pretty plots.......................................................................
# make variable importance into a dataframe
ImpData <- as.data.frame(importance(NDVI_rf))
ImpData$Var.Names <- row.names(ImpData)
ImpData$Var.Names <- factor(ImpData$Var.Names, levels=c("NIR.SWIR1",
                                                        "Pt.B4max",
                                                        "Pt.B5max",
                                                        "precip",
                                                        "TotalFires",
                                                        "Obs_month",
                                                        "tmax",
                                                        "Pt.B4min",
                                                        "Prev.Int",
                                                        "SWIR1.SWIR2") )
# ggplot
ndvi <- ggplot(ImpData, aes(x=Var.Names, y=`%IncMSE`)) +
  geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=`%IncMSE`), color="darkgreen") +
  geom_point(aes(size = IncNodePurity), color="green", alpha=0.6) +
  theme_light() +
  coord_flip() +
  theme(
    legend.position="bottom",
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    text= element_text(size=15),
    axis.ticks.y = element_blank()) +
  xlab("Predictor Variable") +
  ylab("Importance (%IncMSE)") +
  labs(title="NDVI", size="NodePurity")


# TARGET CONDITIONS
summary(BL_Master_df$TotalFires)
# pine = historical interval ~ 5yrs
pine <- BL_Master_df[which(BL_Master_df$EcoType == "Pineland"),]
mean(pine$NDVI[pine$TotalFires == 9]) # NDVI = 0.2503734
mean(pine$NIR.SWIR1[pine$TotalFires == 9]) # NIR:SWIR = 1.178277
mean(pine$SWIR1.SWIR2[pine$TotalFires == 9]) # SWIR1:SWIR2 = 1.278579
mean(pine$Pt.B4max[pine$TotalFires == 9]) # B4max = 18242.24
mean(pine$MFRI[pine$TotalFires == 9]) # MFRI = 3.626805
summary(pine$NDVI[pine$TotalFires == 9])


##########################################################################################################################################################
# VARIABLE MULTICOLINEARITY 
##########################################################################################################################################################
#https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline")

# load BL_master_df
load("./BL_Master_filtered.RDATA")

# Corrilation plots
names(BL_Master_df) # which columns do you want to compare?
BL_Master_df$Obs_month <- as.numeric(BL_Master_df$Obs_month)
BL_Master_noNA <- na.omit(BL_Master_df)

# all considered variables
M <-cor(BL_Master_df[, c("Pt.B1max" ,"Pt.B2max" , "Pt.B3max", "Pt.B4max" , "Pt.B5max", "Pt.B7max" , 
                         "Pt.B1min", "Pt.B2min" ,"Pt.B3min", "Pt.B4min", "Pt.B5min", "Pt.B7min",
                         "NIR.SWIR1", "SWIR1.SWIR2",
                         "Obs_month", "tmin", "tmax", "precip", 
                         "TotalFires", "MFRI", "Prev.Int",
                         "coords.x1", "coords.x2")]) 
corrplot(M, method="circle") 
corrplot(M, method="color")
corrplot(M, method="number")
# varaibles used in model 
M <-cor(BL_Master_df[, c( "Pt.B4max" ,"NIR.SWIR1", "SWIR1.SWIR2",
                          "Obs_month", "tmax", "precip", "tmin",
                          "TotalFires", "MFRI","Prev.Int",
                          "coords.x1", "coords.x2")]) 
corrplot(M, method="circle") 
corrplot(M, method="color")
corrplot(M, method="color", tl.col="black", addCoef.col="grey50")

# Luke's code for cowplot
# Arrange #
legend <- get_legend(gg2 + theme(legend.position="top"))
gg <- align_plots(gg1, gg2,align = "hv", axis = "tblr")

tiff(filename = "BL_corplot.tiff",  
     width = 8, height = 8, res = 600, units = "in",
     compression = "lzw")
corrplot(M, method="color", tl.col="black", addCoef.col="grey50")
dev.off() 



# POST-FIRE NDVI TRENDS...........................................................................................
load("./BL_Master_filtered.RDATA")
# assign frequency categories
summary(BL_Master_filtered$TotalFires)
class(BL_Master_filtered$TotalFires)
summary(BL_Master_filtered$TSF)
class(BL_Master_filtered$TSF)
BL_Master_filtered$TSF <- as.factor(BL_Master_filtered$TSF)
# plot distribution
ggplot(BL_Master_filtered, aes(x=TotalFires)) +
  geom_histogram() +
  scale_x_continuous(breaks=1:9) +
  theme(text = element_text(size = 20))
# plot trend
BL_Master_filtered$TotalFires <- as.factor(BL_Master_filtered$TotalFires)
cols <- c("1" = "#cfccba", "2" = "#e8e5cb", "3" = "#c3d5a3", "4" = "#9bc184",
          "5"= "#669f60", "6"= "#3d7c3c", "7"="#1e5b24", "8"="#1f3d13", "9"="#224313")

ggplot(BL_Master_filtered, aes(x=TSF, 
                      y=NDVI,
                      group=TotalFires,
                      color=TotalFires)) +
  geom_smooth(formula = y ~ s(x, bs = "cs", k=7)) +
  labs(
    x="Years Post-Fire",
    y= "Mean NDVI",
    title = "Baseline Locations",
    color="Total Fires", 
    fill= "Fire Frequency") + 
  #scale_fill_manual(values=cols) +
  #scale_color_manual(values=cols)+
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE) +
  theme_bw()+
  theme(
    text = element_text(size = 20),
    #panel.background = element_rect(fill='transparent'), #transparent panel bg
    #plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    #panel.grid.major = element_blank(), #remove major gridlines
    #panel.grid.minor = element_blank(), #remove minor gridlines
    #legend.background = element_rect(fill='transparent'), #transparent legend bg
    #legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )
# save plot with transparent background
ggsave('FCEprez.png', px, bg='transparent')


ggplot(BL_Master_filtered, aes(x=TotalFires, y=NDVI, color=TotalFires)) +
  geom_boxplot()


