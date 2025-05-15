# LANDSCAPE SUMMARY
# M. Grace McLeod (2024)


# This script summarizes the (1) fire history and (2) climate patterns for the landscape
# across pineland communitites during the study period. 


library(sf)
library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)
library(psych)
library(sp)
library(terra)
library(reshape2)
library(stringr)
library(patchwork)


rm(list=ls())
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")

##############################################################################################################################
# 1. FIRE HISTORY BY VEG TYPE
##############################################################################################################################

# projcrs="+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"


# MASTER DATAFRAME........................................................................................................
# All vegetation classes
veg_info <- st_read("./Veg_layers/veg_info/Veg_info_intersect.shp")
projcrs="+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"
veg_info <- st_as_sf(x = veg_info,                         
                     coords = c("lon", "lat"),
                     crs = projcrs)

# Fire history 
load("./Fire_History/FireYears_df.RDATA")
load("./Fire_History/FireHistory_df.RDATA")

# Fire characteristics
# shp from EVG_AllFires_1978-2020 (see Fire_Data_Edits) extracted to upland_pts_sf.shp (see VegLayers)
# in GIS using the intersect tool
fire_character <- st_read("./Landscape_Summary/Upland_FireCharacter_sp/Upland_FireCharacter.shp")
projcrs="+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"
fire_character <- st_as_sf(x = fire_character,                         
                     coords = c("lon", "lat"),
                     crs = projcrs)

# FIRE HISTORY 
# pull out desired columns from fire history (TF)
names(FireHistory_df)
Fire_by_veg <- FireHistory_df %>% 
  select(ptID, coords.x1, coords.x2, freq_1978_2020)
# pair veg info with fire history
Fire_by_veg <- merge(veg_info, Fire_by_veg, by="ptID")

# calcualte Previous Interval for all locations
# identify last fire year
names(FireYears_df)
FireYears_df$LFY <- apply(FireYears_df[3:45], 1, max, na.rm=TRUE)
# identify penultimate fire year
# make last fire year NA...there must be a more efficient way to do this...but it works
FY_NA <- FireYears_df
FY_NA$EVG_.1978[FY_NA$LFY == 1978] <- NA
FY_NA$EVG_.1979[FY_NA$LFY == 1979] <- NA
FY_NA$EVG_.1980[FY_NA$LFY == 1980] <- NA
FY_NA$EVG_.1981[FY_NA$LFY == 1981] <- NA
FY_NA$EVG_.1982[FY_NA$LFY == 1982] <- NA
FY_NA$EVG_.1983[FY_NA$LFY == 1983] <- NA
FY_NA$EVG_.1984[FY_NA$LFY == 1984] <- NA
FY_NA$EVG_.1985[FY_NA$LFY == 1985] <- NA
FY_NA$EVG_.1986[FY_NA$LFY == 1986] <- NA
FY_NA$EVG_.1987[FY_NA$LFY == 1987] <- NA
FY_NA$EVG_.1988[FY_NA$LFY == 1988] <- NA
FY_NA$EVG_.1989[FY_NA$LFY == 1989] <- NA
FY_NA$EVG_.1990[FY_NA$LFY == 1990] <- NA
FY_NA$EVG_.1991[FY_NA$LFY == 1991] <- NA
FY_NA$EVG_.1992[FY_NA$LFY == 1992] <- NA
FY_NA$EVG_.1993[FY_NA$LFY == 1993] <- NA
FY_NA$EVG_.1994[FY_NA$LFY == 1994] <- NA
FY_NA$EVG_.1995[FY_NA$LFY == 1995] <- NA
FY_NA$EVG_.1996[FY_NA$LFY == 1996] <- NA
FY_NA$EVG_.1997[FY_NA$LFY == 1997] <- NA
FY_NA$EVG_.1998[FY_NA$LFY == 1998] <- NA
FY_NA$EVG_.1999[FY_NA$LFY == 1999] <- NA
FY_NA$EVG_.2000[FY_NA$LFY == 2000] <- NA
FY_NA$EVG_.2001[FY_NA$LFY == 2001] <- NA
FY_NA$EVG_.2002[FY_NA$LFY == 2002] <- NA
FY_NA$EVG_.2003[FY_NA$LFY == 2003] <- NA
FY_NA$EVG_.2004[FY_NA$LFY == 2004] <- NA
FY_NA$EVG_.2005[FY_NA$LFY == 2005] <- NA
FY_NA$EVG_.2006[FY_NA$LFY == 2006] <- NA
FY_NA$EVG_.2007[FY_NA$LFY == 2007] <- NA
FY_NA$EVG_.2008[FY_NA$LFY == 2008] <- NA
FY_NA$EVG_.2009[FY_NA$LFY == 2009] <- NA
FY_NA$EVG_.2010[FY_NA$LFY == 2010] <- NA
FY_NA$EVG_.2011[FY_NA$LFY == 2011] <- NA
FY_NA$EVG_.2012[FY_NA$LFY == 2012] <- NA
FY_NA$EVG_.2013[FY_NA$LFY == 2013] <- NA
FY_NA$EVG_.2014[FY_NA$LFY == 2014] <- NA
FY_NA$EVG_.2015[FY_NA$LFY == 2015] <- NA
FY_NA$EVG_.2016[FY_NA$LFY == 2016] <- NA
FY_NA$EVG_.2017[FY_NA$LFY == 2017] <- NA
FY_NA$EVG_.2018[FY_NA$LFY == 2018] <- NA
FY_NA$EVG_.2019[FY_NA$LFY == 2019] <- NA
FY_NA$EVG_.2020[FY_NA$LFY == 2020] <- NA
# Penultimate fire year is new max
names(FY_NA)
FireYears_df$PenUltFY <- apply(FY_NA[3:45], 1, max, na.rm=TRUE)
# replace inf with NA
FireYears_df$LFY[is.infinite(FireYears_df$LFY)] <- NA
FireYears_df$PenUltFY[is.infinite(FireYears_df$PenUltFY)] <- NA
# calculate Prev.Int
FireYears_df$Prev.Int <- FireYears_df$LFY - FireYears_df$PenUltFY
# calculate TSF
FireYears_df$TSF <- 2020 - FireYears_df$LFY
# pair with fire_by_veg
FireYears_sub <- FireYears_df %>% 
  select(ptID, LFY, PenUltFY, Prev.Int)
# calculate mean fire return interval
Fire_by_veg$mfri <- (2020-1978)/(Fire_by_veg$freq_1978_2020 +1)
# pair veg info with fire history
Fire_by_veg <- merge(FireYears_sub, Fire_by_veg, by="ptID")

# FIRE CHARACTERISTICS
# pull out desired columns 
names(fire_character)
fire_character <- fire_character %>%
  select(ptID, FIRE_ID, FireNumber, FireName, StartDate, EndDate, Year, FireType)
# fix errors
unique(fire_character$FireType)
fire_character$FireType[fire_character$FireType == "Rx"] <- "RX"
# merge with fire_by_veg
Fire_by_veg <- merge(Fire_by_veg, fire_character, by="ptID")

# VEG CLASSES
# clean up veg classes
# replace undetermined values in Level 7
unique(Fire_by_veg$L7_name)
Fire_by_veg$L7_name[Fire_by_veg$L7_name == "Pine Flatwoods-Shrubs"] <- "Pine Flatwoods-Shrub" 
Fire_by_veg$L7_name[Fire_by_veg$L7_name == "Pine Upland" ] <- "Pine Upland-Undetermined" 
Fire_by_veg$L7_name[Fire_by_veg$L7_name == "Undetermined"] <- Fire_by_veg$L6_name[Fire_by_veg$L7_name == "Undetermined"]
Fire_by_veg$L7_name[Fire_by_veg$L7_name == "Pine Flatwood-Saw Palmetto" ] <- "Pine Flatwoods-Saw Palmetto"
# differentiate ecosystem type and dominant community for better visualization
Fire_by_veg$var <- Fire_by_veg$L7_name
Fire_by_veg <- Fire_by_veg %>%
  separate(var, c("EcoSys", "DomCom"), "-")
unique(Fire_by_veg$DomCom)
Fire_by_veg$DomCom[Fire_by_veg$DomCom == "Graminoids"] <- "Graminoid"
Fire_by_veg$DomCom[Fire_by_veg$DomCom == "Shrubs"] <- "Shrub"
Fire_by_veg$DomCom[Fire_by_veg$DomCom == "Undetermined"] <- Fire_by_veg$L6_name[Fire_by_veg$DomCom == "Undetermined"]
# level 4
Fire_by_veg$L4_name[Fire_by_veg$L4_name == "Pine Upland-Shrub"] <- "Pine Upland"
# regional
Fire_by_veg$Regional_VegCat <- "Pine Flatwoods"
Fire_by_veg$Regional_VegCat[Fire_by_veg$L4_name == "Pine Rockland"] <- "Pine Rockland"
Fire_by_veg$Regional_VegCat <- as.factor(Fire_by_veg$Regional_VegCat)

# clean up 
names(Fire_by_veg)
Fire_by_veg$FID_upland <- NULL
Fire_by_veg$coords.x2.y <- NULL
Fire_by_veg$coords.x1.y <- NULL
Fire_by_veg$Uplnds_ <- NULL
Fire_by_veg$FID_pinela <- NULL

# save
save(Fire_by_veg, file="Landscape_Summary.RDATA")
# spatial file
Fire_by_veg_sp <- st_as_sf(x=Fire_by_veg,
                           coords=c("coords.x1.x", "coords.x2.x"),
                           crs=projcrs)
st_write(Fire_by_veg_sp, "./Landscape_Summery", layer="Fire_by_veg_sp.shp", driver="ESRI Shapefile")


# FIRE HISTORY BY VEG CLASS................................................................................................

mean(Fire_by_veg$freq_1978_2020)
sd(Fire_by_veg$freq_1978_2020)
sd(Fire_by_veg$freq_1978_2020) / sqrt(length(Fire_by_veg$freq_1978_2020))

# STATS 
# means table
means_table <- Fire_by_veg %>%
  group_by(EcoSys) %>%
  summarize(mean_freq = mean(freq_1978_2020))
# data table
data_table <- describeBy(Fire_by_veg$freq_1978_2020, group=Fire_by_veg$EcoSys, fast=TRUE, digits = 2)
data_table <- describeBy(Fire_by_veg$freq_1978_2020, group=Fire_by_veg$DomCom, fast=TRUE, digits = 2)

# PREVIOUS INTERVAL
# means table
means_table <- Fire_by_veg %>%
  group_by(EcoSys) %>%
  summarize(mean_freq = mean(Prev.Int))
# data table
describeBy(Fire_by_veg$Prev.Int, group=Fire_by_veg$EcoSys, fast=TRUE, digits = 2)
describeBy(Fire_by_veg$Prev.Int, group=Fire_by_veg$DomCom, fast=TRUE, digits = 2)

# TIME SINCE FIRE 
summary(Fire_by_veg$TSF)
# means table
means_table <- Fire_by_veg %>%
  group_by(EcoSys) %>%
  summarize(mean_freq = mean(TSF))
# data table
describeBy(Fire_by_veg$TSF, group=Fire_by_veg$EcoSys, fast=TRUE, digits = 2)
describeBy(Fire_by_veg$TSF, group=Fire_by_veg$DomCom, fast=TRUE, digits = 2)


# UNBURNED 
load("./Landscape_Summary/Landscape_Summary.RDATA")
# Select unburned areas
unburned <- Fire_by_veg[which(Fire_by_veg$freq_1978_2020 ==0),]
0.0009 * 2859 # = 2.5731 km2


# FIRE TYPE 
# load fire perimeters clipped to pinelands 
Pine_fires_shp <- st_read("./Fire_History/Ecosyst_Clipped_Fires/Pine_fires_Km2.shp")
head(Pine_fires_shp)
# number of fires
WF <- Pine_fires_shp[which(Pine_fires_shp$FireType =="WF"),]
RX <- Pine_fires_shp[which(Pine_fires_shp$FireType =="RX"),]
length(unique(WF$FireNumber))
length(unique(RX$FireNumber))
# annual burned area 
sum(WF$Area_km2) # WF: sum= 693.4558
mean(WF$Area_km2) # WF: mean= 0.606698
sum(RX$Area_km2) # RX: sum= 1418.603
mean(RX$Area_km2) # RX: mean= 1.61941
# data table (for means)
describeBy(Fire_by_veg$freq_1978_2020, group=Fire_by_veg$EcoSys, fast=TRUE, digits = 2)
describeBy(Fire_by_veg$freq_1978_2020, group=Fire_by_veg$DomCom, fast=TRUE, digits = 2)


############################################################################################################################
# 2. CLIMATE SUMMARIES
############################################################################################################################

# DAYMET..................................................................................................................

projcrs="+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"
DAYMETcrs <- "+proj=lcc +lat_0=42.5 +lon_0=-100 +lat_1=25 +lat_2=60 +x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs"
# load DAYMET
tmin <- rast("./Seasonal_Cond/tmin_EVG.TIF")
tmax <- rast("./Seasonal_Cond/tmax_EVG.TIF")
precip <- rast("./Seasonal_Cond/precip_EVG.TIF")
# Load sample pts
smpl_pts <- st_read("./Sample_pts_upland.shp")
# match crs to rasters
crs(smpl_pts)
smpl_pts_vect <- vect("./Sampling/Sample_pts_upland.shp")
smpl_pts_DAYMETcrs <- project(smpl_pts_vect, precip)
smpl_pts_DAYMETcrs <- project(smpl_pts_vect, precip)
# assign dates to raster layers
Date <- read.csv(file="./Seasonal_Cond/DAYMET_dates.csv")
Date$Yr_Mo <- as.Date(paste0(as.character(Date$Date), '01'), format="%Y%m%d")


# tmin
# assign layer names to z dimensions
names(tmin) <- paste("tmin", Date$Yr_Mo, sep=".")
# extract data to each observation based on month
smpl_pts_tmin <- terra::extract(tmin, smpl_pts_DAYMETcrs, method="bilinear", xy=T, bind=T)
# round to 1 decimal
smpl_pts_tmin_df <- as.data.frame(smpl_pts_tmin)
names(smpl_pts_tmin_df)
smpl_pts_tmin_df <- smpl_pts_tmin_df %>% mutate_at(vars(3:254), funs(round(., 1)))
head(smpl_pts_tmin_df)
smpl_pts_tmin_df$Uplnds_ <- NULL
# rotate df longways
names(smpl_pts_tmin_df)
smpl_pts_tmin_df <- reshape2::melt(smpl_pts_tmin_df, na.rm=F, value.name="tmin", id=c("ptID", "x", "y"))
# format columns
smpl_pts_tmin_df$variable <- as.character(smpl_pts_tmin_df$variable)
smpl_pts_tmin_df <- smpl_pts_tmin_df %>%
  transform(variable=str_replace(variable, "tmin.", ""))
smpl_pts_tmin_df$variable <- as.Date(smpl_pts_tmin_df$variable, format="%Y-%m-%d")
smpl_pts_tmin_df$Obs_Year <- format(smpl_pts_tmin_df$variable, "%Y")
smpl_pts_tmin_df$Obs_month <- format(smpl_pts_tmin_df$variable, "%m")
head(smpl_pts_tmin_df)
smpl_pts_tmin_df$variable <- NULL


# tmax
# assign layer names to z dimensions
names(tmax) <- paste("tmax", Date$Yr_Mo, sep=".")
# extract data to each observation based on month
smpl_pts_tmax <- terra::extract(tmax, smpl_pts_DAYMETcrs, method="bilinear", xy=T, bind=T)
# round to 1 decimal
smpl_pts_tmax_df <- as.data.frame(smpl_pts_tmax)
names(smpl_pts_tmax_df)
smpl_pts_tmax_df <- smpl_pts_tmax_df %>% mutate_at(vars(3:254), funs(round(., 1)))
head(smpl_pts_tmax_df)
smpl_pts_tmax_df$Uplnds_ <- NULL
# rotate df longways
names(smpl_pts_tmax_df)
smpl_pts_tmax_df <- reshape2::melt(smpl_pts_tmax_df, na.rm=F, value.name="tmax", id=c("ptID", "x", "y"))
# format columns
smpl_pts_tmax_df$variable <- as.character(smpl_pts_tmax_df$variable)
smpl_pts_tmax_df <- smpl_pts_tmax_df %>%
  transform(variable=str_replace(variable, "tmax.", ""))
smpl_pts_tmax_df$variable <- as.Date(smpl_pts_tmax_df$variable, format="%Y-%m-%d")
smpl_pts_tmax_df$Obs_Year <- format(smpl_pts_tmax_df$variable, "%Y")
smpl_pts_tmax_df$Obs_month <- format(smpl_pts_tmax_df$variable, "%m")
head(smpl_pts_tmax_df)
smpl_pts_tmax_df$variable <- NULL


# precip
# assign dates to raster layers
Date <- read.csv(file="./Seasonal_Cond/DAYMET_dates.csv")
Date$Yr_Mo <- as.Date(paste0(as.character(Date$Date), '01'), format="%Y%m%d")
# assign layer names to z dimensions
names(precip) <- paste("precip", Date$Yr_Mo, sep=".")
# extract data to each observation based on month
smpl_pts_precip <- terra::extract(precip, smpl_pts_DAYMETcrs, method="bilinear", xy=T, bind=T)
# round to 1 decimal
smpl_pts_precip_df <- as.data.frame(smpl_pts_precip)
names(smpl_pts_precip_df)
smpl_pts_precip_df <- smpl_pts_precip_df %>% mutate_at(vars(3:254), funs(round(., 1)))
head(smpl_pts_precip_df)
smpl_pts_precip_df$Uplnds_ <- NULL
# rotate df longways
names(smpl_pts_precip_df)
smpl_pts_precip_df <- reshape2::melt(smpl_pts_precip_df, na.rm=F, value.name="precip", id=c("ptID", "x", "y"))
# format columns
smpl_pts_precip_df$variable <- as.character(smpl_pts_precip_df$variable)
smpl_pts_precip_df <- smpl_pts_precip_df %>%
  transform(variable=str_replace(variable, "precip.", ""))
smpl_pts_precip_df$variable <- as.Date(smpl_pts_precip_df$variable, format="%Y-%m-%d")
smpl_pts_precip_df$Obs_Year <- format(smpl_pts_precip_df$variable, "%Y")
smpl_pts_precip_df$Obs_month <- format(smpl_pts_precip_df$variable, "%m")
head(smpl_pts_precip_df)
smpl_pts_precip_df$variable <- NULL


# DAYMET Master
# join 
Upland_DAYMET <- full_join(smpl_pts_precip_df, smpl_pts_tmax_df, by=c("ptID", "x", "y", "Obs_Year", "Obs_month"), keep=F)
Upland_DAYMET <- full_join(Upland_DAYMET, smpl_pts_tmin_df, by=c("ptID", "x", "y", "Obs_Year", "Obs_month"), keep=F)
# save
setwd("./Seasonal_Cond")
save(Upland_DAYMET, file="Upland_DAYMET.RDATA")



# PDSI ............................................................................................................

# check/rest wd
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")

# LOAD DATA
projcrs="+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"
# load pdsi
pdsi <- rast("./Climate/PDSI_stack.TIF")
# Load sample pts
smpl_pts <- st_read("./Sampling/Sample_pts_upland.shp")
# match crs
smpl_pts_vect <- vect("./Sample_pts_upland.shp")
smpl_pts_PDSIcrs <- project(smpl_pts_vect, pdsi)
smpl_pts_PDSIcrs$Uplnds_ <- NULL

# EXTRACT DATA
# files are too large, run in chunks

# chunk 1
pdsi1 <- pdsi[[1:700]]
smpl_pts_pdsi1 <- terra::extract(pdsi1, smpl_pts_PDSIcrs, method="bilinear", xy=T, bind=T)
# format data
smpl_pts_pdsi1 <- as.data.frame(smpl_pts_pdsi1)
names(smpl_pts_pdsi1)
smpl_pts_pdsi1 <- smpl_pts_pdsi1 %>% mutate_at(vars(2:701), funs(round(., 1)))
head(smpl_pts_pdsi1)
names(smpl_pts_pdsi1)
smpl_pts_pdsi1 <- reshape2::melt(smpl_pts_pdsi1, na.rm=F, value.name="pdsi", id=c("ptID", "x", "y"))
colnames(smpl_pts_pdsi1)[which(names(smpl_pts_pdsi1) == "variable")] <- "Obs_Date"
rm(pdsi1)
smpl_pts_pdsi1$Obs_Date <- as.Date(smpl_pts_pdsi1$Obs_Date, format="%Y-%m-%d")
smpl_pts_pdsi1$Obs_Year <- format(smpl_pts_pdsi1$Obs_Date, "%Y")
#smpl_pts_pdsi1$Obs_month <- format(smpl_pts_pdsi1$Obs_Date, "%m")
# save 
setwd("./Climate")
save(smpl_pts_pdsi1, file="Uplands_PDSI.RDATA")
rm(smpl_pts_pdsi1, smpl_pts_vect, means_pdsi1, pdsi1)
gc()
# calculate means
means_pdsi1 <- smpl_pts_pdsi1 %>%
  group_by(Obs_Year) %>%
  summarize(pdsi=mean(pdsi))
# save
save(means_pdsi1, file="means_pdsi1.RDATA")


# chunk 2
pdsi2 <- pdsi[[701:1100]]
smpl_pts_pdsi2 <- terra::extract(pdsi2, smpl_pts_PDSIcrs, method="bilinear", xy=T, bind=T)
# format data
smpl_pts_pdsi2 <- as.data.frame(smpl_pts_pdsi2)
head(smpl_pts_pdsi2)
names(smpl_pts_pdsi2)
smpl_pts_pdsi2 <- smpl_pts_pdsi2 %>% mutate_at(vars(2:401), funs(round(., 1)))
smpl_pts_pdsi2 <- reshape2::melt(smpl_pts_pdsi2, na.rm=F, value.name="pdsi", id=c("ptID", "x", "y"))
colnames(smpl_pts_pdsi2)[which(names(smpl_pts_pdsi2) == "variable")] <- "Obs_Date"
rm(pdsi2, pdsi,smpl_pts, smpl_pts_vect)
gc()
smpl_pts_pdsi2$Obs_Date <- as.Date(smpl_pts_pdsi2$Obs_Date, format="%Y-%m-%d")
smpl_pts_pdsi2$Obs_Year <- format(smpl_pts_pdsi2$Obs_Date, "%Y")
# save
save(smpl_pts_pdsi2, file="Uplands_PDSI2.RDATA")
# calculate means
means_pdsi2 <- smpl_pts_pdsi2 %>%
  group_by(Obs_Year) %>%
  summarize(pdsi=mean(pdsi))
# save
save(means_pdsi2, file="means_pdsi2.RDATA")


# chunk 3
pdsi3 <- pdsi[[1101:1533]]
smpl_pts_pdsi3 <- terra::extract(pdsi3, smpl_pts_PDSIcrs, method="bilinear", xy=T, bind=T)
# format data
smpl_pts_pdsi3 <- as.data.frame(smpl_pts_pdsi3)
head(smpl_pts_pdsi3)
names(smpl_pts_pdsi3)
smpl_pts_pdsi3 <- smpl_pts_pdsi3 %>% mutate_at(vars(2:434), funs(round(., 1)))
smpl_pts_pdsi3 <- reshape2::melt(smpl_pts_pdsi3, na.rm=F, value.name="pdsi", id=c("ptID", "x", "y"))
colnames(smpl_pts_pdsi3)[which(names(smpl_pts_pdsi3) == "variable")] <- "Obs_Date"
rm(pdsi3, pdsi)
gc()
smpl_pts_pdsi3$Obs_Date <- as.Date(smpl_pts_pdsi3$Obs_Date, format="%Y-%m-%d")
smpl_pts_pdsi3$Obs_Year <- format(smpl_pts_pdsi3$Obs_Date, "%Y")
# save
save(smpl_pts_pdsi3, file="Uplands_pdsi3.RDATA")
# calculate means
means_pdsi3 <- smpl_pts_pdsi3 %>%
  group_by(Obs_Year) %>%
  summarize(pdsi=mean(pdsi))
# save
save(means_pdsi3, file="means_pdsi3.RDATA")


# Combine chunks
# would be great to combine all PDSI data but it's just too big. 
# combine means for figure
load("./means_pdsi1.RDATA")
load("./means_pdsi2.RDATA")
load("./means_pdsi3.RDATA")
# join
mean_pdsi <- rbind(means_pdsi1, means_pdsi2, means_pdsi3)
mean_pdsi <- mean_pdsi %>% mutate_at(vars(2), funs(round(., 1)))
# save
save(mean_pdsi, file="mean_pdsi_combo.RDATA")


# CLIMATE TIME SERIES............................................................................................................

setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")
load("./Seasonal_Cond/Upland_DAYMET.RDATA")
load("./Climate/mean_pdsi_combo.RDATA")

# DAYMET: calculate annual means
means_daymet <- Upland_DAYMET %>%
  group_by(Obs_Year) %>%
  summarize(precip=mean(precip),
            tmax=mean(tmax),
            tmin=mean(tmin))
# round to one decimal
means_daymet <- means_daymet %>% mutate_at(vars(2:4), funs(round(., 1)))
means_daymet$Obs_Year <- as.numeric(means_daymet$Obs_Year)
mean_pdsi$Obs_Year <- as.numeric(mean_pdsi$Obs_Year)

# STATS
tmax <- lm(means_daymet$tmax ~ means_daymet$Obs_Year)
summary(tmax)
tmin <- lm(means_daymet$tmin ~ means_daymet$Obs_Year)
summary(tmin)
precip <- lm(means_daymet$precip ~ means_daymet$Obs_Year)
summary(precip)


