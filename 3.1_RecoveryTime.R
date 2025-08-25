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

library(stats)
library(dplyr)
library(car)
library(devtools)

library(gtools)
library(ggplot2)
library(viridis)
library(splitstackshape)
library(MetBrewer)
library(cowplot)
library(lubridate)
library(patchwork)

library(randomForestCI)
#devtools::install_github("swager/randomForestCI")


setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod") 

##########################################################################################################################################################
# 1. FORMAT DATA
##########################################################################################################################################################

# Copy of Recov_Master 
rm(list=ls())

load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/Recov_Master_082025.RDATA")

Recov_prep <- Recov_Master2 %>% mutate(EndDate = as.Date(EndDate),
                                       StartDate = as.Date(StartDate))

length(unique(Recov_prep$ptID)) # 157,323

# REMOVE UNBURNED POINTS
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Severity/Sev_df.RDATA")

# Evaluate cutoff for severity
length(unique(Sev_df$ptID[is.na(Sev_df$Severity)]))
sev_No_NA <- Sev_df %>% drop_na(Severity) ; length(unique(sev_No_NA$ptID)) # 153,057
sev_No_NA.1 <- sev_No_NA[which(sev_No_NA$Severity > .01),] ; length(unique(sev_No_NA.1$ptID)) # 114,294 (25% reduction, so 25% had no observable change)
rm( sev_No_NA , sev_No_NA.1 )

Sev_df <- Sev_df %>% mutate( FireYear = as.numeric(FireYear ))

Recov_prep_sev <- Recov_prep %>% left_join( Sev_df, by=c("ptID", "StartDate", "EndDate", "FireYear")) %>% filter( Severity >.01) 


# REMOVE OBSERVATIONS PRIOR TO RECOVERY FIRE
# or those with no severity/EndDate
# count them first and record number of obs per pt
counts.total <- Recov_Master2 %>% count(ptID) %>% rename(n_obs_total = n )
counts.post <- Recov_prep_sev  %>% count(ptID) %>% rename(n_obs_post = n )
counts.inx <- counts.total %>% full_join(counts.post, by="ptID")


Recov_prep_final <- Recov_prep_sev %>%  full_join(counts.inx, by="ptID") %>% filter(Obs_Date > EndDate) 
rm( Recov_prep)
rm( counts.total , counts.post, Recov_Master2)
rm( Recov_prep_sev)

# number of remaining pts
length(unique(Recov_prep_final$ptID)) # 114,294

# Check unique ecotypes!
Recov_BL <- Recov_prep_final

##########################################################################################################################################################
# 2. BASELINE MODEL 
##########################################################################################################################################################

# SET VALUES ....................................................................................................................................
# for fire history and spectral variables so baseline represents "ideal" conditions closer to spectral saturation

load(file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/BL_Master_filtered_08082025.RDATA")

# Total Fires
Recov_BL$TotalFires.og <- Recov_BL$TotalFires
# assign sub TotalFires values
Recov_BL$TotalFires <- 9

# Previous Interval 
Recov_BL$Prev.Int.og <- Recov_BL$Prev.Int 
Recov_BL$Prev.Int <- mean(BL_Master_filtered$Prev.Int[BL_Master_filtered$TotalFires == 9])

# Pt.B4max
Recov_BL$Pt.B4max.og <- Recov_BL$Pt.B4max
Recov_BL$Pt.B4max <- mean(BL_Master_filtered$Pt.B4max[BL_Master_filtered$TotalFires == 9])

# Pt.B4min
Recov_BL$Pt.B4min.og <- Recov_BL$Pt.B4min
Recov_BL$Pt.B4min <- mean(BL_Master_filtered$Pt.B4min[BL_Master_filtered$TotalFires == 9]) 

# Pt.B5max
Recov_BL$Pt.B5max.og <- Recov_BL$Pt.B5max
Recov_BL$Pt.B5max <- mean(BL_Master_filtered$Pt.B5max[BL_Master_filtered$TotalFires == 9]) 

# SWIR1.SWIR2
Recov_BL$SWIR1.SWIR2.og <- Recov_BL$SWIR1.SWIR2
Recov_BL$SWIR1.SWIR2 <- mean(BL_Master_filtered$SWIR1.SWIR2[BL_Master_filtered$TotalFires == 9])

# NIR.SWIR1
Recov_BL$NIR.SWIR1.og <- Recov_BL$NIR.SWIR1
Recov_BL$NIR.SWIR1 <- mean(BL_Master_filtered$NIR.SWIR1[BL_Master_filtered$TotalFires == 9]) 

# save
save(Recov_BL, file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recov_BL_sub.RDATA")


# APPLY BASELINE MODEL ...............................................................................................................................

# load baseline model
load( file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/NDVI_rf_08082025.RDATA")
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/BL_train_test_08082025.RDATA")

library(randomForest)
Recov_BL_pinelands <- Recov_BL %>% filter( EcoType == "Pineland")

Recov_BL_pinelands$ptID %>% unique %>% length # 86,303

# RUN MODELS ON RECOVERY POINTS
Recov_BL_pinelands$NDVI_rf <- predict(NDVI_rf, Recov_BL_pinelands)  

# save 
save(Recov_BL, Recov_BL_pinelands, file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/Recov_BL_final.RDATA")


# CALCULATE CONFIDENCE INTERVALS

# subset dataframe (to avoid exhausting vector memory)
# create a blank dataframe that chunk dfs can be bound to
Recov_Combo <- data.frame(matrix(ncol=55, nrow=0)) 
colnames(Recov_Combo) <- colnames(Recov_BL)
# it can handle 1,000,000 at a time for sure so do first 15,000,000 then remainder
# Make list of low ends of the chunk 
Recov_BL_pinelands$ptID %>% length /2

rfpred <- function( df, model){
  m.varhat.sg <- randomForestInfJack(model, df, calibrate = TRUE)
  df$model <- m.varhat.sg$y.hat
  df$var <- sqrt(m.varhat.sg$var.hat)
  
  return(df)
}

# One at a time
a <- seq(1, 6418692, 10000)

# Loop through sub1 
for (i in a) {
  print(i)
  # subset to the chunk
  b <- i+ 9999 
  chunk <-  Recov_BL_pinelands[c(i:b),] 
  # use Jacknife to get variance out of model (apply funciton)  
  # NDVI
  Recov_sub <- rfpred(df=chunk, model=NDVI_rf) 
  colnames(Recov_sub)[which(names(Recov_sub) == "model")] <- "model.NDVI"
  colnames(Recov_sub)[which(names(Recov_sub) == "var")] <- "var.NDVI"
  # append to big dataframe
  Recov_Combo <- gtools::smartbind(Recov_Combo, Recov_sub)
  
}

a <- c(seq(6418693, 12837385, 10000), 12837384)

# Loop through sub1 
for (i in a) {
  print(i)
  # subset to the chunk
  b <- i+ 9999 
  chunk <-  Recov_BL_pinelands[c(i:b),] 
  # use Jacknife to get variance out of model (apply funciton)  
  # NDVI
  Recov_sub <- rfpred(df=chunk, model=NDVI_rf) 
  colnames(Recov_sub)[which(names(Recov_sub) == "model")] <- "model.NDVI"
  colnames(Recov_sub)[which(names(Recov_sub) == "var")] <- "var.NDVI"
  # append to big dataframe
  Recov_Combo <- gtools::smartbind(Recov_Combo, Recov_sub)
  
}

# get lower threshold
Recov_Combo$lwrNDVI <- Recov_Combo$model.NDVI -  Recov_Combo$var.NDVI 
# get upper threshold
Recov_Combo$uprNDVI <- Recov_Combo$model.NDVI + Recov_Combo$var.NDVI

# get immediate post-fire NDVI
IMpre <- Recov_Combo[which(Recov_Combo$PostDate == Recov_Combo$Obs_Date),]
IMpre <- IMpre %>% select(ptID, NDVI)
colnames(IMpre)[which(names(IMpre) == "NDVI")] <- "PostNDVI"
Recov_Combo <- full_join(Recov_Combo, IMpre, by="ptID")
rm(IMpre)
length(unique(Recov_Combo$ptID)) # 114,294

# REMOVE OBS AFTER ANTICEDENT FIRE
load("./Fire_History/FireYears_df.RDATA")
# filter FH to recov_combo pts
Recov_FH <- subset(FireYears_df, ptID %in% Recov_Combo$ptID)
# only care about 2008-2020 bc we know they only burned once before that (recFire)
names(Recov_FH)
Recov_FH <- Recov_FH %>% dplyr::select(ptID,year_2008,year_2009,year_2010,year_2011,
                                       year_2012,year_2013,year_2014,year_2015,year_2016,year_2017,year_2018,year_2019,year_2020, x, y)
# Year of antecedent fire
Recov_FH$AFyear <- NA
Recov_FH$AFyear <- apply(Recov_FH[,2:14], 1, FUN=min,  na.rm = TRUE)
Recov_FH[sapply(Recov_FH, is.infinite)] <- NA
Recov_FH <- Recov_FH %>% dplyr::select(ptID, AFyear)
# full_join with Recov_Combo

Recov_Combo2 <- left_join(Recov_Combo, Recov_FH, by="ptID")

# remove observations after antecedent fire
Recov_Combo3 <- Recov_Combo2[which(Recov_Combo2$Obs_Year < Recov_Combo2$AFyear | is.na(Recov_Combo2$AFyear)),]
summary(Recov_Combo3$AFyear) # should still have NAs
length(unique(Recov_Combo3$ptID)) # 114,058/ 86135

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

# Additional adjustments:
Recov_Combo.final <- Recov_Combo %>% filter(PostDateDif <= 180, Severity > 0.01, PreNDVI >= 0.2 )
Recov_Combo <- Recov_Combo.final

# save
save(Recov_Combo, file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/Recov_Combo_082025.RDATA")

##########################################################################################################################################################
# 3. RECOVERY RATE 
##########################################################################################################################################################
# Who REACHES the threshold AT ALL 
# Time required to reach each threshold, even for points with a higher max. 

load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/Recov_Combo_082025.RDATA")

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

Recov_100$ptID %>% unique %>% length()
Recov_Combo$ptID %>% unique %>% length()

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

Recov_90$ptID %>% unique %>% length()
Recov_Combo$ptID %>% unique %>% length()

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

Recov_80$ptID %>% unique %>% length()
Recov_Combo$ptID %>% unique %>% length()

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

Recov_70$ptID %>% unique %>% length()
Recov_Combo$ptID %>% unique %>% length()

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

Recov_60$ptID %>% unique %>% length()
Recov_Combo$ptID %>% unique %>% length()

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

Recov_50$ptID %>% unique %>% length()
Recov_Combo$ptID %>% unique %>% length()

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

Recov_40$ptID %>% unique %>% length()
Recov_Combo$ptID %>% unique %>% length()

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

Recov_30$ptID %>% unique %>% length()
Recov_Combo$ptID %>% unique %>% length()

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

Recov_20$ptID %>% unique %>% length()
Recov_Combo$ptID %>% unique %>% length()

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

Recov_10$ptID %>% unique %>% length()
Recov_Combo$ptID %>% unique %>% length()


# Make Recovery Rate Dataframe.....................................................................................

# remove columns that corespond to only one observation

Recov_10n <- Recov_10 %>%
  select(ptID, x, y,
         StartDate, EndDate, FireYear, PenUltFY, AFyear, FireName, FireNumber, 
         FireType, TotalFires, Prev.Int, PreNBR, PreNDVI, 
         PostNBR, PostNDVI, PreDate, PostDate, DateDif, Severity, PreDateDif, PostDateDif,
         Rec10_Date, Rec10_Yrs, Rec10_NDVI)
Recov_20n <- Recov_20 %>%
  select(ptID,Rec20_Date, Rec20_Yrs, Rec20_NDVI)
Recov_30n <- Recov_30 %>%
  select(ptID,Rec30_Date, Rec30_Yrs, Rec30_NDVI)
Recov_40n <- Recov_40 %>%
  select(ptID, Rec40_Date, Rec40_Yrs, Rec40_NDVI)
Recov_50n <- Recov_50 %>%
  select(ptID,Rec50_Date, Rec50_Yrs, Rec50_NDVI)
Recov_60n <- Recov_60 %>%
  select(ptID, Rec60_Date, Rec60_Yrs, Rec60_NDVI)
Recov_70n <- Recov_70 %>%
  select(ptID,Rec70_Date, Rec70_Yrs, Rec70_NDVI)
Recov_80n <- Recov_80 %>%
  select(ptID,Rec80_Date, Rec80_Yrs, Rec80_NDVI)
Recov_90n <- Recov_90 %>%
  select(ptID, Rec90_Date, Rec90_Yrs, Rec90_NDVI)
Recov_100n <- Recov_100 %>%
  select(ptID,Rec100_Date, Rec100_Yrs, Rec100_NDVI)

# full_join threhold dataframes
Recov_Rate1 <- full_join(Recov_10n, Recov_20n)
Recov_Rate2 <- full_join(Recov_Rate1, Recov_30n)
Recov_Rate3 <- full_join(Recov_Rate2, Recov_40n)
Recov_Rate4 <- full_join(Recov_Rate3, Recov_50n)
Recov_Rate5 <- full_join(Recov_Rate4, Recov_60n)
Recov_Rate6 <- full_join(Recov_Rate5, Recov_70n)
Recov_Rate7 <- full_join(Recov_Rate6, Recov_80n)
Recov_Rate8 <- full_join(Recov_Rate7, Recov_90n)
Recov_Rate <- full_join(Recov_Rate8, Recov_100n)
head(Recov_Rate)


# Save
save(Recov_Rate, file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/Recov_Rates_082025.RDATA")

##########################################################################################################################################################
# 4. MAX THRESHOLD
##########################################################################################################################################################

load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/Recov_Combo_082025.RDATA")

# 100% .........................................................................................................
# reach lower limit of model confidence interval 
# pull out all obs with NDVI > lowNDVI  


recov.time <- function(df, threshold, adj){
  
  Recov <- df[which(df$NDVI >= df$lwrNDVI * threshold),] %>% anti_join( adj, by="ptID")
  # summarize by pt to get min date (first date with NDVI above threshold)
  Rec_Date <- Recov  %>% 
    group_by(ptID) %>%
    summarise(Rec_Date = min(Obs_Date, na.rm = T)) 
  # join min_date with Recov_Combo by ptID
  Recov_final <- left_join(Recov, Rec_Date, by="ptID")
  rm(Rec_Date)
  
  Recov_final %>% names
  Recov_Time <- Recov_final %>% 
    filter(Rec_Date == Obs_Date) %>%
    dplyr::select(ptID, Rec_Date, EndDate, FireYear, Severity,
                  B5, uprNDVI, lwrNDVI, model.NDVI, var.NDVI, NDVI, Severity, PreNDVI ) %>%
    distinct()
  # calculate recovery time
  Recov_Time$Rec_Date <- as.Date(Recov_Time$Rec_Date)
  Recov_Time$EndDate <- as.Date(Recov_Time$EndDate)
  Recov_Time$Rec_Days <- Recov_Time$Rec_Date - Recov_Time$EndDate
  Recov_Time$Rec_Yrs <- as.numeric(Recov_Time$Rec_Days) /356
  Recov_Time_final <- Recov_Time%>%
    select(ptID, B5, uprNDVI, lwrNDVI, model.NDVI, var.NDVI, NDVI, Rec_Date, Rec_Yrs,Severity, PreNDVI )
  
  return( Recov_Time_final)
}

Recov_Time_100 <- recov.time(df = Recov_Combo , threshold=1 , adj= data.frame(ptID=0) ) %>% mutate(threshold = "100")

Recov_Time_90 <- recov.time(df = Recov_Combo, threshold=0.9 , adj=Recov_Time_100) %>% mutate(threshold = "90")

Recov_Time_80 <- recov.time(df = Recov_Combo, threshold=0.8 , adj=rbind(Recov_Time_100, Recov_Time_90)) %>% mutate(threshold = "80")

Recov_Time_70 <- recov.time(df = Recov_Combo, threshold=0.7 , adj=rbind(Recov_Time_100, Recov_Time_90, Recov_Time_80)) %>% mutate(threshold = "70")

Recov_Time_60 <- recov.time(df = Recov_Combo, threshold=0.6 , adj=rbind(Recov_Time_100, Recov_Time_90, Recov_Time_80, Recov_Time_70)) %>% mutate(threshold = "60")

Recov_Time_50 <- recov.time(df = Recov_Combo, threshold=0.5 , adj=rbind(Recov_Time_100, Recov_Time_90, Recov_Time_80, Recov_Time_70, Recov_Time_60))%>% mutate(threshold = "50")

Recov_Time_40 <- recov.time(df = Recov_Combo, threshold=0.4 , adj=rbind(Recov_Time_100, Recov_Time_90, Recov_Time_80, Recov_Time_70, Recov_Time_60, Recov_Time_50)) %>% mutate(threshold = "<50")

Recov_Time_30 <- recov.time(df = Recov_Combo, threshold=0.3 , adj=rbind(Recov_Time_100, Recov_Time_90, Recov_Time_80, Recov_Time_70, Recov_Time_60, Recov_Time_50, Recov_Time_40))  %>% mutate(threshold = "<50")

Recov_Time_20 <- recov.time(df = Recov_Combo, threshold=0.2 , adj=rbind(Recov_Time_100, Recov_Time_90, Recov_Time_80, Recov_Time_70, Recov_Time_60, Recov_Time_50, Recov_Time_40, Recov_Time_30))  %>% mutate(threshold = "<50")
 
Recov_Time_10 <- recov.time(df = Recov_Combo, threshold=0.1 , adj=rbind(Recov_Time_100, Recov_Time_90, Recov_Time_80, Recov_Time_70, Recov_Time_60, Recov_Time_50, Recov_Time_40, Recov_Time_30, Recov_Time_20))  %>% mutate(threshold = "<50")


# Data Prep
# add threshold column
Recov_Time_100$thrshold <- "100%"
Recov_Time_90$thrshold <- "90%"
Recov_Time_80$thrshold <- "80%"
Recov_Time_70$thrshold <- "70%"
Recov_Time_60$thrshold <- "60%"
Recov_Time_50$thrshold <- "50%"
Recov_Time_40$thrshold <- "<50%"
Recov_Time_30$thrshold <- "<50%"
Recov_Time_20$thrshold <- "<50%"

Recov_Time_100$R_percent <- 100
Recov_Time_90$R_percent <- 90
Recov_Time_80$R_percent <- 80
Recov_Time_70$R_percent <- 70
Recov_Time_60$R_percent <- 60
Recov_Time_50$R_percent <- 50
Recov_Time_40$R_percent <- 40
Recov_Time_30$R_percent <- 30
Recov_Time_20$R_percent <- 20



# merge 
thresholds <- rbind(Recov_Time_100, Recov_Time_90, Recov_Time_80, Recov_Time_70, Recov_Time_60, Recov_Time_50, Recov_Time_40, Recov_Time_30, Recov_Time_20)

thresholds %>% names
# full_join together
Max_Thresh <- rbind(Recov_Time_100, Recov_Time_90, Recov_Time_80, Recov_Time_70,
                    Recov_Time_60, Recov_Time_50, Recov_Time_40, Recov_Time_30, Recov_Time_20)

head(Max_Thresh)
unique(Max_Thresh$threshold)
# make threshold a factor
Max_Thresh$threshold <- as.factor(Max_Thresh$threshold)

# Threshold DF:

thresholds$thrshold <- as.factor(thresholds$thrshold)
thresholds$thrshold <- factor(thresholds$thrshold, levels=c("<50%", "50%", "60%", "70%", "80%", "90%", "100%") )
levels(thresholds$thrshold)

thresholds$Rec_Yrs %>% summary


rec_tbl_new <- thresholds %>% reframe( .by = thrshold,
                                       R_percent =  R_percent %>% as.numeric,
                                       locations = length(ptID %>% unique),
                                       mean_time =  mean(Rec_Yrs, na.rm=T), 
                                       max_time =  max(Rec_Yrs, na.rm=T), 
                                       min_time=  min(Rec_Yrs, na.rm=T),
                                       model.NDVI = mean(model.NDVI),
                                       PreNDVI = mean(PreNDVI)) %>% distinct %>% 
  mutate( 
    locations_percent = locations/sum(locations)*100,
    locations_percent = cumsum(locations_percent),
    sample_size =  cumsum(locations),
    rec_time_range = max_time - min_time)

# SAVE
save(Max_Thresh, 
     Recov_Time_100, Recov_Time_90, Recov_Time_80, Recov_Time_70,
     Recov_Time_60, Recov_Time_50, Recov_Time_40, Recov_Time_30, Recov_Time_20,
     rec_tbl_new,thresholds,
     
     file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/Max_Thresh_082025.RDATA")

##########################################################################################################################################################
# 5. INTEGRATE ALL TIME RECOVERY DATA
##########################################################################################################################################################


load(file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/Max_Thresh_082025.RDATA")
load(file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/Recov_Rates_082025.RDATA")

Recov_Drivers <- full_join(Recov_Rate, Max_Thresh, by="ptID") 

# Save 
save(Recov_Drivers, file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/Recov_Drivers_082025.RDATA")
