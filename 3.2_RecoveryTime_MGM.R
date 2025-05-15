# RECOVERY TIME
# M.Grace McLeod (2023)

# This script 
# 1. formats recovery point data for the baseline model ("Recov_BL.RDATA")
# 2. applies the baseline NDVI model to the recovery points and calculates confidence intervals ("Recov_Combo.RDATA")
# 3. calculates the time to each threshold, even if max threshold is higher ("Recov_Rates.RDATA")
# 4. calculates recovery time to maximum threshold ("Max_Thresh.RDATA")
# 5. generates a single dataframe with all recovery time information ("Recov_Drivers.RDATA")


# FIX LA PROBLEMA
# list variables in datasets, make them MATCH 
# keep track of number of NAs, min and max, classes
# test: filter out NAs and see if thats the only filter that happens
# looks like there is a filtering step somewhere in the new df


rm(list=ls())

library(randomForest)
library(stats)
library(dplyr)
library(car)
library(devtools)
library(randomForestCI)
library(gtools)
library(ggplot2)
library(viridis)
library(splitstackshape)
library(MetBrewer)
library(cowplot)
library(lubridate)
library(patchwork)
#install_github("swager/randomForestCI")

# function to calculate performance error for random forests (Author: Sparkle Malone)
rfpred <- function( df, model){
  m.varhat.sg <- randomForestInfJack(model, df, calibrate = TRUE)
  df$model <- m.varhat.sg$y.hat
  df$var <- sqrt(m.varhat.sg$var.hat)
  
  return(df)
}


setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod") 

##########################################################################################################################################################
# 1. FORMAT DATA
##########################################################################################################################################################

# Copy of Recov_Master 
load("./Recovery/Recov_Master.RDATA")
Recov_BL <- Recov_Master
rm(Recov_Master)
length(unique(Recov_BL$ptID)) # 157,323

# REMOVE UNBURNED POINTS
# or those with no severity/EndDate
# load Severity dataframe 
load("./Severity/Sev_df.RDATA")
# match classes
summary(Sev_df)
summary(Recov_BL)
Recov_BL$EndDate <- as.Date(Recov_BL$EndDate)
Recov_BL$StartDate <- as.Date(Recov_BL$StartDate)
# merge dataframes
Recov_BL <- merge(Recov_BL, Sev_df, by=c("ptID", "StartDate", "EndDate", "FireYear"))
head(Recov_BL)
rm(Sev_df)
# filter severity to a minimum
length(unique(Sev_df$ptID[is.na(Sev_df$Severity)]))
sev_No_NA <- Sev_df %>% drop_na(Severity) ; length(unique(sev_No_NA$ptID)) # 153,057
sev_No_NA.1 <- sev_No_NA[which(sev_No_NA$Severity > .01),] ; length(unique(sev_No_NA.1$ptID)) # 114,294 (25% reduction, so 25% had no observable change)
Recov_BL <- Recov_BL[which(Recov_BL$Severity >.01),] 
# number of remaining pts
length(unique(Recov_BL$ptID)) # 114,294


# REMOVE OBSERVATIONS PRIOR TO RECOVERY FIRE
# count them first and record number of obs per pt
test1 <- Recov_Master %>% count(ptID)
colnames(test1)[which(names(test1) == "n")] <- "n_obs_total"
Recov_BL <- merge(Recov_BL, test1, by="ptID") 
# remove observations prior to EndDate
Recov_BL <- Recov_BL[which(Recov_BL$Obs_Date > Recov_BL$EndDate),]
# count number of obs post-fire
test2 <- Recov_BL  %>% count(ptID)
colnames(test2)[which(names(test2) == "n")] <- "n_obs_post"
Recov_BL <- merge(Recov_BL, test2, by="ptID") 

# Make FireYear numeric
Recov_BL$FireYear <- as.numeric(Recov_BL$FireYear)

# Set Prev.Int with NAs to 80
Recov_BL$Prev.Int [is.na(Recov_BL$Prev.Int)] =  80
summary(Recov_BL$Prev.Int)

# save 
setwd("./Recovery") 
save(Recov_BL, file="Recov_BL.RDATA")

##########################################################################################################################################################
# 2. BASELINE MODEL 
##########################################################################################################################################################

# SET VALUES ....................................................................................................................................
# for fire history and spectral variables so baseline represents "ideal" conditions closer to spectral saturation
# use baseline points 
load("./Baseline/BL_Master_filtered.RDATA")
#load("./Recovery/Recov_BL.RDATA")

# Total Fires
Recov_BL$TotalFires.sub <- NA
# assign sub TotalFires values
Recov_BL$TotalFires.sub <- 9
head(Recov_BL)
# has to replace TotalFires for model to work on it
Recov_BL$TotalFires.og <- Recov_BL$TotalFires # make a copy
Recov_BL$TotalFires <- Recov_BL$TotalFires.sub # sub in the sub values

# Previous Interval 
# use the values associated with the selected fire history values
mean(BL_Master_filtered$Prev.Int[BL_Master_filtered$TotalFires == 9]) # 6
# assign sub values
Recov_BL$Prev.Int.sub <- 6
head(Recov_BL)
Recov_BL$Prev.Int.og <- Recov_BL$Prev.Int 
Recov_BL$Prev.Int <- Recov_BL$Prev.Int.sub 

# Pt.B4max
# use the values associated with the selected fire history values
mean(BL_Master_filtered$Pt.B4max[BL_Master_filtered$TotalFires == 9]) # 21526.95
# assign substitutes
Recov_BL$Pt.B4max.sub <- 21526.95
head(Recov_BL)
Recov_BL$Pt.B4max.og <- Recov_BL$Pt.B4max
Recov_BL$Pt.B4max <- Recov_BL$Pt.B4max.sub

# Pt.B4min
# use the values associated with the selected fire history values
mean(BL_Master_filtered$Pt.B4min[BL_Master_filtered$TotalFires == 9]) # 10637.73
# assign substitutes
Recov_BL$Pt.B4min.sub <- 10637.73
head(Recov_BL)
Recov_BL$Pt.B4min.og <- Recov_BL$Pt.B4min
Recov_BL$Pt.B4min <- Recov_BL$Pt.B4min.sub

# Pt.B5max
# use the values associated with the selected fire history values
mean(BL_Master_filtered$Pt.B5max[BL_Master_filtered$TotalFires == 9]) # 15097.51
# assign substitutes
Recov_BL$Pt.B5max.sub <-15097.51
head(Recov_BL)
Recov_BL$Pt.B5max.og <- Recov_BL$Pt.B5max
Recov_BL$Pt.B5max <- Recov_BL$Pt.B5max.sub

# SWIR1.SWIR2
# use the values associated with the selected fire history values
mean(BL_Master_filtered$SWIR1.SWIR2[BL_Master_filtered$TotalFires == 9]) # 1.259455
# assign substitutes
Recov_BL$SWIR1.SWIR2.sub <- 1.259455
head(Recov_BL)
Recov_BL$SWIR1.SWIR2.og <- Recov_BL$SWIR1.SWIR2
Recov_BL$SWIR1.SWIR2 <- Recov_BL$SWIR1.SWIR2.sub

# NIR.SWIR1
# use the values associated with the selected fire history values
mean(BL_Master_filtered$NIR.SWIR1[BL_Master_filtered$TotalFires == 9]) # 1.429776
# assign substitutes
Recov_BL$NIR.SWIR1.sub <- 1.429776
head(Recov_BL)
Recov_BL$NIR.SWIR1.og <- Recov_BL$NIR.SWIR1
Recov_BL$NIR.SWIR1 <- Recov_BL$NIR.SWIR1.sub


# save
#setwd("I:/Malone Lab/ENP Fire/Grace_McLeod/Recovery") 
#save(Recov_BL, file="Recov_BL_sub.RDATA")


# APPLY BASELINE MODEL ...............................................................................................................................

# load baseline model
load("./Baseline/NDVI_rf.RDATA")

# RUN MODELS ON RECOVERY POINTS
Recov_BL$NDVI_rf <- predict(NDVI_rf, Recov_BL)  # r2=   
summary(lm(Recov_BL$NDVI ~ Recov_BL$NDVI_rf))

# save 
save(Recov_BL, file="Recov_BL_final.RDATA")


# CALCULATE CONFIDENCE INTERVALS

# subset dataframe (to avoid exhausting vector memory)
# create a blank dataframe that chunk dfs can be bound to
Recov_Combo <- data.frame(matrix(ncol=55, nrow=0)) 
colnames(Recov_Combo) <- colnames(Recov_BL)
# it can handle 1,000,000 at a time for sure so do first 15,000,000 then remainder
# Make list of low ends of the chunk 
a <- seq(1, 150000001, 10000)
a <- seq(14240001, 15500001, 10000)


# Loop through sub1 
for (i in a) {
  print(i)
  # subset to the chunk
  b <- i+ 9999 
  chunk <-  Recov_BL[c(i:b),] 
  # use Jacknife to get variance out of model (apply funciton)  
  # NDVI
  Recov_sub <- rfpred(df=chunk, model=NDVI_rf) 
  colnames(Recov_sub)[which(names(Recov_sub) == "model")] <- "model.NDVI"
  colnames(Recov_sub)[which(names(Recov_sub) == "var")] <- "var.NDVI"
  # append to big dataframe
  Recov_Combo <- smartbind(Recov_Combo, Recov_sub)
  
}

# run on remaining data 
Recov_sub <- rfpred (df=Recov_BL[15510001:15560817,], model=NDVI_rf)
colnames(Recov_sub)[which(names(Recov_sub) == "model")] <- "model.NDVI"
colnames(Recov_sub)[which(names(Recov_sub) == "var")] <- "var.NDVI"
Recov_Combo <- smartbind(Recov_Combo, Recov_sub)


# get lower threshold
Recov_Combo$lwrNDVI <- Recov_Combo$model.NDVI -  Recov_Combo$var.NDVI 
# get upper threshold
Recov_Combo$uprNDVI <- Recov_Combo$model.NDVI + Recov_Combo$var.NDVI

# get immediate post-fire NDVI
IMpre <- Recov_Combo[which(Recov_Combo$PostDate == Recov_Combo$Obs_Date),]
IMpre <- IMpre %>% select(ptID, NDVI)
colnames(IMpre)[which(names(IMpre) == "NDVI")] <- "PostNDVI"
Recov_Combo <- merge(Recov_Combo, IMpre, by="ptID")
rm(IMpre)
length(unique(Recov_Combo$ptID)) # 114,294

# REMOVE OBS AFTER ANTICEDENT FIRE
load("./Fire_History/FireYears_df.RDATA")
# filter FH to recov_combo pts
Recov_FH <- subset(FireYears_df, ptID %in% Recov_Combo$ptID)
# only care about 2008-2020 bc we know they only burned once before that (recFire)
names(Recov_FH)
Recov_FH <- Recov_FH %>% dplyr::select(ptID,EVG_.2008,EVG_.2009,EVG_.2010,EVG_.2011,
                                       EVG_.2012,EVG_.2013,EVG_.2014,EVG_.2015,EVG_.2016,EVG_.2017,EVG_.2018,EVG_.2019,EVG_.2020, coords.x1, coords.x2)
# Year of antecedent fire
Recov_FH$AFyear <- NA
Recov_FH$AFyear <- apply(Recov_FH[,2:14], 1, FUN=min,  na.rm = TRUE)
Recov_FH[sapply(Recov_FH, is.infinite)] <- NA
Recov_FH <- Recov_FH %>% dplyr::select(ptID, AFyear)
# merge with Recov_Combo
Recov_Combo2 <- merge(Recov_Combo, Recov_FH, by="ptID", all.x = T)
# remove observations after antecedent fire
Recov_Combo3 <- Recov_Combo2[which(Recov_Combo2$Obs_Year < Recov_Combo2$AFyear | is.na(Recov_Combo2$AFyear)),]
summary(Recov_Combo3$AFyear) # should still have NAs
length(unique(Recov_Combo3$ptID)) # 114,058

# rename
Recov_Combo <- Recov_Combo3
rm(Recov_Combo2, Recov_Combo3)

# Re-assign original values to substituted columns
Recov_Combo$Pt.B4min <- Recov_Combo$Pt.B4min.og ; Recov_Combo$Pt.B4min.og <- NULL ;  Recov_Combo$Pt.B4min.sub <- NULL 
Recov_Combo$Pt.B4max <- Recov_Combo$Pt.B4max.og ; Recov_Combo$Pt.B4max.og <- NULL ;  Recov_Combo$Pt.B4max.sub <- NULL 
Recov_Combo$Pt.B5max <- Recov_Combo$Pt.B5max.og ; Recov_Combo$Pt.B5max.og <- NULL ;  Recov_Combo$Pt.B5max.sub <- NULL 
Recov_Combo$Pt.B7max <- Recov_Combo$Pt.B7max.og ; Recov_Combo$Pt.B7max.og <- NULL ;  Recov_Combo$Pt.B7max.sub <- NULL 
Recov_Combo$SWIR1.SWIR2 <- Recov_Combo$SWIR1.SWIR2.og ; Recov_Combo$SWIR1.SWIR2.og <- NULL ; Recov_Combo$SWIR1.SWIR2.sub <- NULL
Recov_Combo$NIR.SWIR1 <- Recov_Combo$NIR.SWIR1.og ; Recov_Combo$NIR.SWIR1.og <- NULL ; Recov_Combo$NIR.SWIR1.sub <- NULL
Recov_Combo$NDVI_rf <- NULL


# save
setwd("./Recovery") 
save(Recov_Combo, file="Recov_Combo.RDATA")



##########################################################################################################################################################
# 3. RECOVERY RATE 
##########################################################################################################################################################
# Who REACHES the threshold AT ALL 
# Time required to reach each threshold, even for points with a higher max. 

load("./Recovery/Recov_Combo.RDATA")

# 100% (reach lower limit of baseline)...............................................................................................
# pull out all obs with NDVI > threshold  
Recov_100 <- Recov_Combo[which(Recov_Combo$NDVI >= (Recov_Combo$lwrNDVI * 1)),]
# summarize by pt to get min date (first date with NDVI above threshold)
Rec_Date <- Recov_100 %>% 
  group_by(ptID) %>%
  summarise(Rec100_Date = min(Obs_Date, na.rm = T)) 
# join min_date with Recov_threshold info by ptID
Recov_100 <- left_join(Recov_100, Rec_Date, by="ptID")
rm(Rec_Date)

# CALCULATE RECOVERY TIME
# between fire end date and rec date
class(Recov_100$Rec100_Date)
Recov_100$Rec100_Date <- as.Date(Recov_100$Rec100_Date)
Recov_100$EndDate <- as.Date(Recov_100$EndDate)
Recov_100$Rec100_Yrs <- as.numeric((Recov_100$Rec100_Date - Recov_100$EndDate) /356)
# Keep only the recovery observation for each point
Recov_100 <- Recov_100[which(Recov_100$Obs_Date == Recov_100$Rec100_Date),]
Recov_100$Rec100_NDVI <- Recov_100$NDVI


# 90% ...............................................................................................
Recov_90 <- Recov_Combo[which(Recov_Combo$NDVI >= (Recov_Combo$lwrNDVI * .9)),]
Rec_Date <- Recov_90 %>% 
  group_by(ptID) %>%
  summarise(Rec90_Date = min(Obs_Date, na.rm = T)) 
Recov_90 <- left_join(Recov_90, Rec_Date, by="ptID")
rm(Rec_Date)
Recov_90$Rec90_Date <- as.Date(Recov_90$Rec90_Date)
Recov_90$EndDate <- as.Date(Recov_90$EndDate)
Recov_90$Rec90_Yrs <- as.numeric((Recov_90$Rec90_Date - Recov_90$EndDate) /356)
Recov_90 <- Recov_90[which(Recov_90$Obs_Date == Recov_90$Rec90_Date),]
Recov_90$Rec90_NDVI <- Recov_90$NDVI


# 80% ...............................................................................................
Recov_80 <- Recov_Combo[which(Recov_Combo$NDVI >= (Recov_Combo$lwrNDVI * .8)),]
Rec_Date <- Recov_80 %>% 
  group_by(ptID) %>%
  summarise(Rec80_Date = min(Obs_Date, na.rm = T)) 
Recov_80 <- left_join(Recov_80, Rec_Date, by="ptID")
rm(Rec_Date)
Recov_80$Rec80_Date <- as.Date(Recov_80$Rec80_Date)
Recov_80$EndDate <- as.Date(Recov_80$EndDate)
Recov_80$Rec80_Yrs <- as.numeric((Recov_80$Rec80_Date - Recov_80$EndDate) /356)
Recov_80 <- Recov_80[which(Recov_80$Obs_Date == Recov_80$Rec80_Date),]
Recov_80$Rec80_NDVI <- Recov_80$NDVI


# 70% ...............................................................................................
Recov_70 <- Recov_Combo[which(Recov_Combo$NDVI >= (Recov_Combo$lwrNDVI * .7)),]
Rec_Date <- Recov_70 %>% 
  group_by(ptID) %>%
  summarise(Rec70_Date = min(Obs_Date, na.rm = T)) 
Recov_70 <- left_join(Recov_70, Rec_Date, by="ptID")
rm(Rec_Date)
Recov_70$Rec70_Date <- as.Date(Recov_70$Rec70_Date)
Recov_70$EndDate <- as.Date(Recov_70$EndDate)
Recov_70$Rec70_Yrs <- as.numeric((Recov_70$Rec70_Date - Recov_70$EndDate) /356)
Recov_70 <- Recov_70[which(Recov_70$Obs_Date == Recov_70$Rec70_Date),]
Recov_70$Rec70_NDVI <- Recov_70$NDVI


# 60% ...............................................................................................
Recov_60 <- Recov_Combo[which(Recov_Combo$NDVI >= (Recov_Combo$lwrNDVI * .6)),]
Rec_Date <- Recov_60 %>% 
  group_by(ptID) %>%
  summarise(Rec60_Date = min(Obs_Date, na.rm = T)) 
Recov_60 <- left_join(Recov_60, Rec_Date, by="ptID")
rm(Rec_Date)
Recov_60$Rec60_Date <- as.Date(Recov_60$Rec60_Date)
Recov_60$EndDate <- as.Date(Recov_60$EndDate)
Recov_60$Rec60_Yrs <- as.numeric((Recov_60$Rec60_Date - Recov_60$EndDate) /356)
Recov_60 <- Recov_60[which(Recov_60$Obs_Date == Recov_60$Rec60_Date),]
Recov_60$Rec60_NDVI <- Recov_60$NDVI


# 50% ...............................................................................................
Recov_50 <- Recov_Combo[which(Recov_Combo$NDVI >= (Recov_Combo$lwrNDVI * .5)),]
Rec_Date <- Recov_50 %>% 
  group_by(ptID) %>%
  summarise(Rec50_Date = min(Obs_Date, na.rm = T)) 
Recov_50 <- left_join(Recov_50, Rec_Date, by="ptID")
rm(Rec_Date)
Recov_50$Rec50_Date <- as.Date(Recov_50$Rec50_Date)
Recov_50$EndDate <- as.Date(Recov_50$EndDate)
Recov_50$Rec50_Yrs <- as.numeric((Recov_50$Rec50_Date - Recov_50$EndDate) /356)
Recov_50 <- Recov_50[which(Recov_50$Obs_Date == Recov_50$Rec50_Date),]
Recov_50$Rec50_NDVI <- Recov_50$NDVI


# 40% ...............................................................................................
Recov_40 <- Recov_Combo[which(Recov_Combo$NDVI >= (Recov_Combo$lwrNDVI * .4)),]
Rec_Date <- Recov_40 %>% 
  group_by(ptID) %>%
  summarise(Rec40_Date = min(Obs_Date, na.rm = T)) 
Recov_40 <- left_join(Recov_40, Rec_Date, by="ptID")
rm(Rec_Date)
Recov_40$Rec40_Date <- as.Date(Recov_40$Rec40_Date)
Recov_40$EndDate <- as.Date(Recov_40$EndDate)
Recov_40$Rec40_Yrs <- as.numeric((Recov_40$Rec40_Date - Recov_40$EndDate) /356)
Recov_40 <- Recov_40[which(Recov_40$Obs_Date == Recov_40$Rec40_Date),]
Recov_40$Rec40_NDVI <- Recov_40$NDVI


# 30% ...............................................................................................
Recov_30 <- Recov_Combo[which(Recov_Combo$NDVI >= (Recov_Combo$lwrNDVI * .3)),]
Rec_Date <- Recov_30 %>% 
  group_by(ptID) %>%
  summarise(Rec30_Date = min(Obs_Date, na.rm = T)) 
Recov_30 <- left_join(Recov_30, Rec_Date, by="ptID")
rm(Rec_Date)
Recov_30$Rec30_Date <- as.Date(Recov_30$Rec30_Date)
Recov_30$EndDate <- as.Date(Recov_30$EndDate)
Recov_30$Rec30_Yrs <- as.numeric((Recov_30$Rec30_Date - Recov_30$EndDate) /356)
Recov_30 <- Recov_30[which(Recov_30$Obs_Date == Recov_30$Rec30_Date),]
Recov_30$Rec30_NDVI <- Recov_30$NDVI

# 20% ...............................................................................................
Recov_20 <- Recov_Combo[which(Recov_Combo$NDVI >= (Recov_Combo$lwrNDVI * .2)),]
Rec_Date <- Recov_20 %>% 
  group_by(ptID) %>%
  summarise(Rec20_Date = min(Obs_Date, na.rm = T)) 
Recov_20 <- left_join(Recov_20, Rec_Date, by="ptID")
rm(Rec_Date)
Recov_20$Rec20_Date <- as.Date(Recov_20$Rec20_Date)
Recov_20$EndDate <- as.Date(Recov_20$EndDate)
Recov_20$Rec20_Yrs <- as.numeric((Recov_20$Rec20_Date - Recov_20$EndDate) /356)
Recov_20 <- Recov_20[which(Recov_20$Obs_Date == Recov_20$Rec20_Date),]
Recov_20$Rec20_NDVI <- Recov_20$NDVI

# 10% ...............................................................................................
Recov_10 <- Recov_Combo[which(Recov_Combo$NDVI >= (Recov_Combo$lwrNDVI * .1)),]
Rec_Date <- Recov_10 %>% 
  group_by(ptID) %>%
  summarise(Rec10_Date = min(Obs_Date, na.rm = T)) 
Recov_10 <- left_join(Recov_10, Rec_Date, by="ptID")
rm(Rec_Date)
Recov_10$Rec10_Date <- as.Date(Recov_10$Rec10_Date)
Recov_10$EndDate <- as.Date(Recov_10$EndDate)
Recov_10$Rec10_Yrs <- as.numeric((Recov_10$Rec10_Date - Recov_10$EndDate) /356)
Recov_10 <- Recov_10[which(Recov_10$Obs_Date == Recov_10$Rec10_Date),]
Recov_10$Rec10_NDVI <- Recov_10$NDVI


# Make Recovery Rate Dataframe.....................................................................................

# remove columns that corespond to only one observation
Recov_10 <- Recov_10 %>%
  select(ptID, coords.x1, coords.x2,
         StartDate, EndDate, FireYear, PenUltFY, AFyear, FireName, FireNumber, 
         FireType, TotalFires, Prev.Int, PreNBR, PreNDVI, 
         PostNBR, PostNDVI, PreDate, PostDate, DateDif, Severity, PreDateDif, PostDateDif,
         Rec10_Date, Rec10_Yrs, Rec10_NDVI)
Recov_20 <- Recov_20 %>%
  select(ptID, coords.x1, coords.x2,
         Rec20_Date, Rec20_Yrs, Rec20_NDVI)
Recov_30 <- Recov_30 %>%
  select(ptID, coords.x1, coords.x2,
         Rec30_Date, Rec30_Yrs, Rec30_NDVI)
Recov_40 <- Recov_40 %>%
  select(ptID, coords.x1, coords.x2,  
         Rec40_Date, Rec40_Yrs, Rec40_NDVI)
Recov_50 <- Recov_50 %>%
  select(ptID, coords.x1, coords.x2,  
         Rec50_Date, Rec50_Yrs, Rec50_NDVI)
Recov_60 <- Recov_60 %>%
  select(ptID, coords.x1, coords.x2,  
         Rec60_Date, Rec60_Yrs, Rec60_NDVI)
Recov_70 <- Recov_70 %>%
  select(ptID, coords.x1, coords.x2, 
         Rec70_Date, Rec70_Yrs, Rec70_NDVI)
Recov_80 <- Recov_80 %>%
  select(ptID, coords.x1, coords.x2,  
         Rec80_Date, Rec80_Yrs, Rec80_NDVI)
Recov_90 <- Recov_90 %>%
  select(ptID, coords.x1, coords.x2, 
         Rec90_Date, Rec90_Yrs, Rec90_NDVI)
Recov_100 <- Recov_100 %>%
  select(ptID, coords.x1, coords.x2, 
         Rec100_Date, Rec100_Yrs, Rec100_NDVI)

# merge threhold dataframes
Recov_Rate <- merge(Recov_10, Recov_20, all.x=T)
Recov_Rate <- merge(Recov_Rate, Recov_30, all.x=T)
Recov_Rate <- merge(Recov_Rate, Recov_40, all.x=T)
Recov_Rate <- merge(Recov_Rate, Recov_50, all.x=T)
Recov_Rate <- merge(Recov_Rate, Recov_60, all.x=T)
Recov_Rate <- merge(Recov_Rate, Recov_70, all.x=T)
Recov_Rate <- merge(Recov_Rate, Recov_80, all.x=T)
Recov_Rate <- merge(Recov_Rate, Recov_90, all.x=T)
Recov_Rate <- merge(Recov_Rate, Recov_100, all.x=T)
head(Recov_Rate)

# Save
setwd("./Recovery")
save(Recov_Rate, file="Recov_Rates.RDATA")

##########################################################################################################################################################
# 4. MAX THRESHOLD
##########################################################################################################################################################

load("./Recovery/Recov_Combo.RDATA")

# 100% .........................................................................................................
# reach lower limit of model confidence interval 
# pull out all obs with NDVI > lowNDVI  
Recov_100 <- Recov_Combo[which(Recov_Combo$NDVI >= Recov_Combo$lwrNDVI),]
# summarize by pt to get min date (first date with NDVI above threshold)
Rec_Date <- Recov_100 %>% 
  group_by(ptID) %>%
  summarise(Rec_Date = min(Obs_Date, na.rm = T)) 
# join min_date with Recov_Combo by ptID
Recov_100 <- left_join(Recov_100, Rec_Date, by="ptID")
rm(Rec_Date)
# how many pts reach target?
length(unique(Recov_100$ptID)) 
# CALCULATE RECOVERY TIME
# pull out distinct pt info
Recov_Time_100 <- Recov_100 %>% 
  filter(Rec_Date == Obs_Date) %>%
  dplyr::select(ptID, Rec_Date, EndDate, FireYear, Severity,
                Pt.B5max, uprNDVI, lwrNDVI, model.NDVI, var.NDVI, NDVI) %>%
  distinct()
# calculate recovery time
Recov_Time_100$Rec_Date <- as.Date(Recov_Time_100$Rec_Date)
Recov_Time_100$EndDate <- as.Date(Recov_Time_100$EndDate)
Recov_Time_100$Rec_Days <- Recov_Time_100$Rec_Date - Recov_Time_100$EndDate
Recov_Time_100$Rec_Yrs <- as.numeric(Recov_Time_100$Rec_Days) /356


# 90%...............................................................................................
# remove pts that recovered to 100
Recov_90 <- anti_join(Recov_Combo, Recov_Time_100, by="ptID")
# calculate 90% of lower limit
Recov_90$lwr90 <- Recov_90$lwrNDVI * 0.9
# filter for obs at or above threshold
Recov_90 <- Recov_90[which(Recov_90$NDVI >= Recov_90$lwr90),]
length(unique(Recov_90$ptID)) #19,553
# summarize by pt to get min date (first date with NDVI above threshold)
Rec_Date <- Recov_90 %>% 
  group_by(ptID) %>%
  summarise(Rec_Date = min(Obs_Date, na.rm = T)) 
# join min_date with Recov_Combo by ptID
Recov_90 <- left_join(Recov_90, Rec_Date, by="ptID")
rm(Rec_Date)
# how many pts reach target?
length(unique(Recov_90$ptID)) 
# CALCULATE RECOVERY TIME
# pull out distinct pt info
Recov_Time_90 <- Recov_90 %>% 
  filter(Rec_Date == Obs_Date) %>%
  dplyr::select(ptID, Rec_Date, EndDate, FireYear, Severity,
                Pt.B5max, uprNDVI, lwrNDVI, model.NDVI, var.NDVI, NDVI) %>%
  distinct()
# calculate recovery time
class(Recov_Time_90$EndDate)
Recov_Time_90$Rec_Date <- as.Date(Recov_Time_90$Rec_Date)
Recov_Time_90$EndDate <- as.Date(Recov_Time_90$EndDate)
Recov_Time_90$Rec_Days <- Recov_Time_90$Rec_Date - Recov_Time_90$EndDate
Recov_Time_90$Rec_Yrs <- as.numeric(Recov_Time_90$Rec_Days) /356

# 80%...............................................................................................
# remove pts that recovered beyond threshold
Recov_80 <- anti_join(Recov_Combo, Recov_Time_100, by="ptID")
Recov_80 <- anti_join(Recov_80, Recov_Time_90, by="ptID")
# calculate 80% of lower limit
Recov_80$lwr80 <- Recov_80$lwrNDVI * 0.8
# filter for obs at or above threshold
Recov_80 <- Recov_80[which(Recov_80$NDVI >= Recov_80$lwr80),]
# summarize by pt to get min date (first date with NDVI above threshold)
Rec_Date <- Recov_80 %>% 
  group_by(ptID) %>%
  summarise(Rec_Date = min(Obs_Date, na.rm = T)) 
# join min_date with Recov_Combo by ptID
Recov_80 <- left_join(Recov_80, Rec_Date, by="ptID")
rm(Rec_Date)
# how many pts reach target?
length(unique(Recov_80$ptID)) 
# CALCULATE RECOVERY TIME
# pull out distinct pt info
Recov_Time_80 <- Recov_80 %>% 
  filter(Rec_Date == Obs_Date) %>%
  dplyr::select(ptID, Rec_Date, EndDate, FireYear, Severity,
                Pt.B5max, uprNDVI, lwrNDVI, model.NDVI, var.NDVI, NDVI) %>%
  distinct()
# calculate recovery time
class(Recov_Time_80$EndDate)
Recov_Time_80$Rec_Date <- as.Date(Recov_Time_80$Rec_Date)
Recov_Time_80$EndDate <- as.Date(Recov_Time_80$EndDate)
Recov_Time_80$Rec_Days <- Recov_Time_80$Rec_Date - Recov_Time_80$EndDate
Recov_Time_80$Rec_Yrs <- as.numeric(Recov_Time_80$Rec_Days) /356

# 70%...............................................................................................
# remove pts that recovered beyond threshold
Recov_70 <- anti_join(Recov_Combo, Recov_Time_100, by="ptID")
Recov_70 <- anti_join(Recov_70, Recov_Time_90, by="ptID")
Recov_70 <- anti_join(Recov_70, Recov_Time_80, by="ptID")
# calculate 70% of lower limit
Recov_70$lwr70 <- Recov_70$lwrNDVI * 0.7
# filter for obs at or above threshold
Recov_70 <- Recov_70[which(Recov_70$NDVI >= Recov_70$lwr70),]
# summarize by pt to get min date (first date with NDVI above threshold)
Rec_Date <- Recov_70 %>% 
  group_by(ptID) %>%
  summarise(Rec_Date = min(Obs_Date, na.rm = T)) 
# join min_date with Recov_Combo by ptID
Recov_70 <- left_join(Recov_70, Rec_Date, by="ptID")
rm(Rec_Date)
# how many pts reach target?
length(unique(Recov_70$ptID)) 
# CALCULATE RECOVERY TIME
# pull out distinct pt info
Recov_Time_70 <- Recov_70 %>% 
  filter(Rec_Date == Obs_Date) %>%
  dplyr::select(ptID, Rec_Date, EndDate, FireYear, Severity,
                Pt.B5max, uprNDVI, lwrNDVI, model.NDVI, var.NDVI, NDVI) %>%
  distinct()
# calculate recovery time
class(Recov_Time_70$EndDate)
Recov_Time_70$Rec_Date <- as.Date(Recov_Time_70$Rec_Date)
Recov_Time_70$EndDate <- as.Date(Recov_Time_70$EndDate)
Recov_Time_70$Rec_Days <- Recov_Time_70$Rec_Date - Recov_Time_70$EndDate
Recov_Time_70$Rec_Yrs <- as.numeric(Recov_Time_70$Rec_Days) /356



# 60%...............................................................................................
# remove pts that recovered beyond threshold
Recov_60 <- anti_join(Recov_Combo, Recov_Time_100, by="ptID")
Recov_60 <- anti_join(Recov_60, Recov_Time_90, by="ptID")
Recov_60 <- anti_join(Recov_60, Recov_Time_80, by="ptID")
Recov_60 <- anti_join(Recov_60, Recov_Time_70, by="ptID")
# calculate 60% of lower limit
Recov_60$lwr60 <- Recov_60$lwrNDVI * 0.6
# filter for obs at or above threshold
Recov_60 <- Recov_60[which(Recov_60$NDVI >= Recov_60$lwr60),]
# summarize by pt to get min date (first date with NDVI above threshold)
Rec_Date <- Recov_60 %>% 
  group_by(ptID) %>%
  summarise(Rec_Date = min(Obs_Date, na.rm = T)) 
# join min_date with Recov_Combo by ptID
Recov_60 <- left_join(Recov_60, Rec_Date, by="ptID")
rm(Rec_Date)
# how many pts reach target?
length(unique(Recov_60$ptID)) 
# CALCULATE RECOVERY TIME
# pull out distinct pt info
Recov_Time_60 <- Recov_60 %>% 
  filter(Rec_Date == Obs_Date) %>%
  dplyr::select(ptID, Rec_Date, EndDate, FireYear, Severity,
                Pt.B5max, uprNDVI, lwrNDVI, model.NDVI, var.NDVI, NDVI) %>%
  distinct()
# calculate recovery time
class(Recov_Time_60$EndDate)
Recov_Time_60$Rec_Date <- as.Date(Recov_Time_60$Rec_Date)
Recov_Time_60$EndDate <- as.Date(Recov_Time_60$EndDate)
Recov_Time_60$Rec_Days <- Recov_Time_60$Rec_Date - Recov_Time_60$EndDate
Recov_Time_60$Rec_Yrs <- as.numeric(Recov_Time_60$Rec_Days) /356


# 50%...............................................................................................
# remove pts that recovered beyond threshold
Recov_50 <- anti_join(Recov_Combo, Recov_Time_100, by="ptID")
Recov_50 <- anti_join(Recov_50, Recov_Time_90, by="ptID")
Recov_50 <- anti_join(Recov_50, Recov_Time_80, by="ptID")
Recov_50 <- anti_join(Recov_50, Recov_Time_70, by="ptID")
Recov_50 <- anti_join(Recov_50, Recov_Time_60, by="ptID")
# calculate 50% of lower limit
Recov_50$lwr50 <- Recov_50$lwrNDVI * 0.5
# filter for obs at or above threshold
Recov_50 <- Recov_50[which(Recov_50$NDVI >= Recov_50$lwr50),]
# summarize by pt to get min date (first date with NDVI above threshold)
Rec_Date <- Recov_50 %>% 
  group_by(ptID) %>%
  summarise(Rec_Date = min(Obs_Date, na.rm = T)) 
# join min_date with Recov_Combo by ptID
Recov_50 <- left_join(Recov_50, Rec_Date, by="ptID")
rm(Rec_Date)
# how many pts reach target?
length(unique(Recov_50$ptID)) 
# CALCULATE RECOVERY TIME
# pull out distinct pt info
Recov_Time_50 <- Recov_50 %>% 
  filter(Rec_Date == Obs_Date) %>%
  dplyr::select(ptID, Rec_Date, EndDate, FireYear, Severity,
                Pt.B5max, uprNDVI, lwrNDVI, model.NDVI, var.NDVI, NDVI) %>%
  distinct()
# calculate recovery time
class(Recov_Time_50$EndDate)
Recov_Time_50$Rec_Date <- as.Date(Recov_Time_50$Rec_Date)
Recov_Time_50$EndDate <- as.Date(Recov_Time_50$EndDate)
Recov_Time_50$Rec_Days <- Recov_Time_50$Rec_Date - Recov_Time_50$EndDate
Recov_Time_50$Rec_Yrs <- as.numeric(Recov_Time_50$Rec_Days) /356


# <50%...............................................................................................
# remove pts that recovered beyond threshold
Recov_undr50 <- anti_join(Recov_Combo, Recov_Time_100, by="ptID")
Recov_undr50 <- anti_join(Recov_undr50, Recov_Time_90, by="ptID")
Recov_undr50 <- anti_join(Recov_undr50, Recov_Time_80, by="ptID")
Recov_undr50 <- anti_join(Recov_undr50, Recov_Time_70, by="ptID")
Recov_undr50 <- anti_join(Recov_undr50, Recov_Time_60, by="ptID")
Recov_undr50 <- anti_join(Recov_undr50, Recov_Time_50, by="ptID")
# calculate undr50% of lower limit
Recov_undr50$lwrundr50 <- Recov_undr50$lwrNDVI * 0.5
# filter for obs at or above threshold
Recov_undr50 <- Recov_undr50[which(Recov_undr50$NDVI < Recov_undr50$lwrundr50),]
# summarize by pt to get min date (first date with NDVI above threshold)
Rec_Date <- Recov_undr50 %>% 
  group_by(ptID) %>%
  summarise(Rec_Date = min(Obs_Date, na.rm = T)) 
# join min_date with Recov_Combo by ptID
Recov_undr50 <- left_join(Recov_undr50, Rec_Date, by="ptID")
rm(Rec_Date)
# how many pts reach target?
length(unique(Recov_undr50$ptID)) 
# CALCULATE RECOVERY TIME
# pull out distinct pt info
Recov_Time_undr50 <- Recov_undr50 %>% 
  filter(Rec_Date == Obs_Date) %>%
  dplyr::select(ptID, Rec_Date, EndDate, FireYear, Severity,
                Pt.B5max, uprNDVI, lwrNDVI, model.NDVI, var.NDVI, NDVI) %>%
  distinct()
# calculate recovery time
class(Recov_Time_undr50$EndDate)
Recov_Time_undr50$Rec_Date <- as.Date(Recov_Time_undr50$Rec_Date)
Recov_Time_undr50$EndDate <- as.Date(Recov_Time_undr50$EndDate)
Recov_Time_undr50$Rec_Days <- Recov_Time_undr50$Rec_Date - Recov_Time_undr50$EndDate
Recov_Time_undr50$Rec_Yrs <- as.numeric(Recov_Time_undr50$Rec_Days) /356


# MERGE.................................................................................................
# Select desired columns
Recov_Time_100 <- Recov_Time_100 %>%
  select(ptID, Pt.B5max, uprNDVI, lwrNDVI, model.NDVI, var.NDVI, NDVI, Rec_Date, Rec_Yrs)
Recov_Time_90 <- Recov_Time_90 %>%
  select(ptID,  Pt.B5max, uprNDVI, lwrNDVI, model.NDVI, var.NDVI, NDVI, Rec_Date, Rec_Yrs)
Recov_Time_80 <- Recov_Time_80 %>%
  select(ptID,  Pt.B5max, uprNDVI, lwrNDVI, model.NDVI, var.NDVI, NDVI, Rec_Date, Rec_Yrs)
Recov_Time_70 <- Recov_Time_70 %>%
  select(ptID,  Pt.B5max, uprNDVI, lwrNDVI, model.NDVI, var.NDVI, NDVI, Rec_Date, Rec_Yrs)
Recov_Time_60 <- Recov_Time_60 %>%
  select(ptID,  Pt.B5max, uprNDVI, lwrNDVI, model.NDVI, var.NDVI, NDVI, Rec_Date, Rec_Yrs)
Recov_Time_50 <- Recov_Time_50 %>%
  select(ptID,  Pt.B5max, uprNDVI, lwrNDVI, model.NDVI, var.NDVI, NDVI, Rec_Date, Rec_Yrs)
Recov_Time_undr50 <- Recov_Time_undr50 %>%
  select(ptID,  Pt.B5max, uprNDVI, lwrNDVI, model.NDVI, var.NDVI, NDVI, Rec_Date, Rec_Yrs)

# add coumn for threshold
Recov_Time_100$threshold <- "100"
Recov_Time_90$threshold <- "90"
Recov_Time_80$threshold <- "80"
Recov_Time_70$threshold <- "70"
Recov_Time_60$threshold <- "60"
Recov_Time_50$threshold <- "50"
Recov_Time_undr50$threshold <- "<50"
# merge together
Max_Thresh <- rbind(Recov_Time_100, Recov_Time_90)
Max_Thresh <- rbind(Max_Thresh, Recov_Time_80)
Max_Thresh <- rbind(Max_Thresh, Recov_Time_70)
Max_Thresh <- rbind(Max_Thresh, Recov_Time_60)
Max_Thresh <- rbind(Max_Thresh, Recov_Time_50)
Max_Thresh <- rbind(Max_Thresh, Recov_Time_undr50)
head(Max_Thresh)
unique(Max_Thresh$threshold)
# make threshold a factor
Max_Thresh$threshold <- as.factor(Max_Thresh$threshold)

# SAVE
setwd("./Recovery")
save(Max_Thresh, file="Max_Thresh.RDATA")

##########################################################################################################################################################
# 5. INTEGRATE ALL TIME RECOVERY DATA
##########################################################################################################################################################

Recov_Drivers <- merge(Recov_Rate, Max_Thresh, by="ptID") 

# Save 
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery")
save(Recov_Drivers, file="Recov_Drivers.RDATA")









