# Recovery MASTER
# Sparkle Malone, updated and corrected Graces code.

# This script extracts spectral, seasonal condition, and fire history data to Recovery sample points

rm(list=ls())

library(sf)
library(tidyverse)
library(ggplot2)


# load Recovery sample points
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")
Recov_smpl_pts <- sf::st_read('/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Sampling/Recov_smpl_pts.shp') %>% st_transform(crs = 4326) # sf file

# Add lat and long into the file:
Recov_smpl_pts$lat <- sf::st_coordinates(Recov_smpl_pts)[,1]
Recov_smpl_pts$lon <- sf::st_coordinates(Recov_smpl_pts)[,2]

Recov_df <- as.data.frame(Recov_smpl_pts)


##########################################################################################################################################################
# SPECTRAL DATA 
##########################################################################################################################################################

# load Spec_Master
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Image_Processing/Master_Spec/Spec_Master.RDATA")


Recov_Spec  <- Spec_Master %>% mutate(Obs_Date = as.Date(Obs_Date, format= "%Y%m%d"),
                                       Obs_Year = format(Obs_Date, "%Y") %>% as.numeric(),
                                       Obs_month = format(Obs_Date, "%m")) %>% filter( ptID %in% Recov_df$ptID) %>% mutate(  B1= (B1* 0.0000275) + (-0.2), 
                                                                                                                                                       B2= (B2* 0.0000275) + (-0.2), 
                                                                                                                                                       B3= (B3* 0.0000275) + (-0.2),
                                                                                                                                                       B4= (B4* 0.0000275) + (-0.2), 
                                                                                                                                                       B5= (B5* 0.0000275) + (-0.2),
                                                                                                                                                       B7= (B7* 0.0000275) + (-0.2))


#save 
save(Recov_Spec, file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/Recov_Spec_082025.RDATA")

##########################################################################################################################################################
# SEASONAL CONDITION DATA
##########################################################################################################################################################

SeasonalCond_sp <- Recov_smpl_pts

setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")
tmin_EVG <- terra::rast("./Seasonal_Cond/tmin_EVG.tif") %>% terra::project( 'epsg:4326') 
tmax_EVG <- terra::rast("./Seasonal_Cond/tmax_EVG.tif") %>% terra::project( 'epsg:4326') 
precip_EVG <- terra::rast("./Seasonal_Cond/precip_EVG.tif") %>% terra::project( 'epsg:4326') 

# DAYMET DATES
Date <-read.csv(file="./Seasonal_Cond/DAYMET_dates.csv") # list of dates has been pre-generated
Date$Yr_Mo <- as.Date(paste0(as.character(Date$Date), '01'), format='%Y%m%d')

# Adjust Rasters:

names(tmin_EVG) <- paste("tmin", Date$Yr_Mo, sep=".")
terra::time(tmin_EVG) <- Date$Yr_Mo

names(tmax_EVG ) <- paste("tmax", Date$Yr_Mo, sep=".")
terra::time(tmax_EVG ) <- Date$Yr_Mo

names(precip_EVG) <- paste("precip", Date$Yr_Mo, sep=".")
terra::time(precip_EVG) <- Date$Yr_Mo

# MIN TEMPERATURE.........................................................................................................................................
  
Recov_tmin_df <- terra::extract(tmin_EVG, SeasonalCond_sp, exact=TRUE) %>% as.data.frame() %>% mutate_at(vars(1:252), funs(round(., 1))) %>% cbind(SeasonalCond_sp) 

# turn DF longways
Recov_tmin_long <- tidyr::pivot_longer(Recov_tmin_df, cols=c(2:253), names_to = "variable", values_to = "value")
  
# format
Recov_tmin <- Recov_tmin_long %>% mutate( var = as.character(variable)) %>% 
  transform(variable=stringr::str_replace(var,"tmin.","")) %>% mutate(tmin = value,
                                                                      date = as.Date(variable, format= "%Y-%m-%d"),
                                                                      Obs_Year = format(date, "%Y"),
                                                                      Obs_month = format(date, "%m")) %>% dplyr::select(ptID, lat, lon, Obs_Year, Obs_month, tmin)

# MAX TEMPERATURES.........................................................................................................................................

# extract data for each variable to each observation based on observation month
Recov_tmax_df <- terra::extract(tmax_EVG, SeasonalCond_sp, exact=TRUE) %>% as.data.frame() %>% mutate_at(vars(1:252), funs(round(., 1))) %>% cbind(SeasonalCond_sp)

# turn DF longways
Recov_tmax_long <- tidyr::pivot_longer(Recov_tmax_df, cols=c(2:253), names_to = "variable", values_to = "value")

Recov_tmax <- Recov_tmax_long %>% mutate( var = as.character(variable)) %>% 
  transform(variable= stringr::str_replace(var,"tmax.","")) %>% mutate(tmax = value,
                                                                       date = as.Date(variable, format= "%Y-%m-%d"),
                                                                       Obs_Year = format(date, "%Y"),
                                                                       Obs_month = format(date, "%m")) %>% dplyr::select(ptID, lat, lon, Obs_Year, Obs_month, tmax)

# PRECIPITATION.........................................................................................................................................

# extract data for each variable to each observation based on observation month
Recov_precip_df <- terra::extract(precip_EVG, SeasonalCond_sp, exact=TRUE) %>% as.data.frame() %>% mutate_at(vars(1:252), funs(round(., 1))) %>% cbind(SeasonalCond_sp)

# turn DF longways
Recov_precip_long <- tidyr::pivot_longer(Recov_precip_df, cols=c(2:253), names_to = "variable", values_to = "value")

# format
Recov_precip <- Recov_precip_long %>% mutate( var = as.character(variable)) %>% 
  transform(variable= stringr::str_replace(var,"precip.","")) %>% mutate(precip = value,
                                                                         date = as.Date(variable, format= "%Y-%m-%d"),
                                                                         Obs_Year = format(date, "%Y"),
                                                                         Obs_month = format(date, "%m")) %>% dplyr::select(ptID, lat, lon, Obs_Year, Obs_month, precip)


# MERGE DAYMET DATAFRAMES.........................................................................................................................................
# use join() instead of merge() because it's much faster with large dfs
# could use any join type for this because number of observations are the same for all dfs. 
Recov_DAYMET1 = full_join(Recov_tmin, Recov_tmax, by=c("ptID",  "lat", "lon", "Obs_Year", "Obs_month"), keep=FALSE) 

Recov_DAYMET <- Recov_DAYMET1 %>% full_join( Recov_precip, by=c("ptID",  "lat", "lon", "Obs_Year", "Obs_month"), keep=FALSE) 

# save DF
save(Recov_DAYMET, file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/Recov_DAYMET_082025.RDATA")
write.csv(Recov_DAYMET, file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/Recov_DAYMET.csv")

rm(list=ls())

##########################################################################################################################################################
# FIRE HISTORY DATA
##########################################################################################################################################################

# FIRE-SPECIFIC INFO.......................................................................................................................................
# Extract from shapefiles
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")
Recov_smpl_pts <- sf::st_read('/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Sampling/Recov_smpl_pts.shp') %>% st_transform(crs = 4326) # sf file

# Add lat and long into the file:
Recov_smpl_pts$lat <- sf::st_coordinates(Recov_smpl_pts)[,1]
Recov_smpl_pts$lon <- sf::st_coordinates(Recov_smpl_pts)[,2]


# load RECOV.master.RDATA or RECOV.master.shp (See Fire_Data_Edits Script)
load(file.path("/", "Volumes", "malonelab", "Research", "ENP", "ENP Fire", "FireHistory", "RECOV.master_082025.RDATA"))

# make copy and rename Recovery points
Recov_Fire_info <- Recov_smpl_pts 

# check crs
Recov_Fire_info <- sf::st_transform(Recov_Fire_info, "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")
RECOV.master.sp <- RECOV.master.sp %>%  sf::st_transform( "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")

# extract shapefile info for each point
TABLE <- sf::st_intersection(Recov_Fire_info, RECOV.master.sp)

# view in df
Recov_Fire_info_df <- as.data.frame(TABLE) %>% rename(
  FireName = Fire_Nm, 
  FireNumber = Fr_Nmbr,
  FireName = Fire_Nm, 
  FireNumber = Fr_Nmbr,
  StartDate = Disc_Dt , 
  EndDate = Dcld_Dt ,
  FireYear= Year,
  FireType = Fir_Typ,
  StartDate = Disc_Dt,  
  EndDate = Dcld_Dt, 
  FireYear= Year,
  FireType = Fir_Typ) %>% select(-File_Nm.1 ,-Fire_ID.1,, -Fr_Nmbr.1 ,
                                    -Fire_Nm.1 ,-Year.1 ,   -Disc_Dt.1,-Dcld_Dt.1, 
                                    -Dat_Flg.1,-Fir_Typ.1)


Recov_Fire_info <- TABLE %>% rename(
  FireName = Fire_Nm, 
  FireNumber = Fr_Nmbr,
  FireName = Fire_Nm, 
  FireNumber = Fr_Nmbr,
  StartDate = Disc_Dt , 
  EndDate = Dcld_Dt ,
  FireYear= Year,
  FireType = Fir_Typ,
  StartDate = Disc_Dt,  
  EndDate = Dcld_Dt, 
  FireYear= Year,
  FireType = Fir_Typ) %>% select(-File_Nm.1 ,-Fire_ID.1,, -Fr_Nmbr.1 ,
                                 -Fire_Nm.1 ,-Year.1 ,   -Disc_Dt.1,-Dcld_Dt.1, 
                                 -Dat_Flg.1,-Fir_Typ.1)

# SAVE 
st_write(Recov_Fire_info, dsn="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/Recov_Fire_Info_082025.shp", append=FALSE)

save(Recov_Fire_info, Recov_Fire_info_df, RECOV.master.sp,  file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/Recov_Fire_Info_082025.RDATA")


# FIRE HISTORY METRICS........................................................................................................................................
# make a copy for fire history
Fire_Histry_sp <- Recov_smpl_pts

# TOTAL FIRES
# total number of fires experienced by the point prior to AND including the recovery fire
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Fire_History/FireHistory_df.RDATA")

# Keep only rows with point ID that matches Recov_Smpl_pts
Fire_Histry_sp <- sf::st_transform(Fire_Histry_sp,"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")
Fire_Histry_sp <- as.data.frame(Fire_Histry_sp)

Recov_Frequency  <- filter(FireHistory_df, ptID %in% Fire_Histry_sp$ptID)
Recov_Frequency %>% names()

# make fire history dataframe
Recov_FireHistory <- Recov_Frequency %>% dplyr::select(ptID, x, y)

# calculate total fires
Recov_Frequency %>% names
Recov_FireHistory$TotalFires <- rowSums(Recov_Frequency[ , c(2:31)]) # 1978-2007
summary(Recov_FireHistory)

# PREVIOUS INTERVAL 
# number of years since fire prior to the recovery fire
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Fire_History/FireYears_df.RDATA")
names(FireYears_df)
# Keep only rows with point ID that matches Recov_Smpl_pts
Recov_FireYear  <- filter(FireYears_df, ptID %in% Fire_Histry_sp$ptID)

Recov_FireYear %>% summary
# remove all fires/years after 2007
Recov_FireYear %>% names
Recov_FireYear.2007 <- Recov_FireYear[2:31] 
Recov_FireYear.2007  %>% summary

# identify recovery fire year
names(Recov_FireYear.2007)
Recov_FireYear.2007$RecFY <- apply(Recov_FireYear.2007,1, max, na.rm=TRUE)
# identify penultimate fire year
# make rec fire year NA
Recov_FireNA <- Recov_FireYear.2007
names(Recov_FireNA)
Recov_FireNA$year_2001[Recov_FireNA$RecFY == 2001] <- NA
Recov_FireNA$year_2002[Recov_FireNA$RecFY == 2002] <- NA
Recov_FireNA$year_2003[Recov_FireNA$RecFY == 2003] <- NA
Recov_FireNA$year_2004[Recov_FireNA$RecFY == 2004] <- NA
Recov_FireNA$year_2005[Recov_FireNA$RecFY == 2005] <- NA
Recov_FireNA$year_2006[Recov_FireNA$RecFY == 2006] <- NA
Recov_FireNA$year_2007[Recov_FireNA$RecFY == 2007] <- NA

# Penultimate fire year is new max
names(Recov_FireNA)
Recov_FireYear.2007$PenUltFY <- apply(Recov_FireNA[1:30], 1, max, na.rm=TRUE)
# check it
view(Recov_FireYear.2007)

# add to Recov_FireHistory
Recov_FireYear.2007[sapply(Recov_FireYear.2007, is.infinite)] <- NA

Recov_FireHistory$RecFY <- Recov_FireYear.2007$RecFY
Recov_FireHistory$PenUltFY <-Recov_FireYear.2007$PenUltFY
Recov_FireHistory$Prev.Int <- Recov_FireHistory$RecFY - Recov_FireHistory$PenUltFY
Recov_FireHistory$Prev.Int [is.na(Recov_FireHistory$Prev.Int)] =  80

# convert inf to NA 
Recov_FireHistory[sapply(Recov_FireHistory, is.infinite)] <- NA

# MERGE THE TWO DATAFRAMES..........................................................................................................................................

Recov_FireHistory$FireYear <- Recov_FireHistory$RecFY
# Not the right amount of points. Examine the join more carefully!
Recov_FireHistory %>% summary()
Recov_Fire_info_df %>% summary()

Recov_FireHistory_final <- left_join( Recov_FireHistory,
                                      Recov_Fire_info_df, 
                                      by=c( 'ptID', 'FireYear')) %>% distinct()

Recov_FireHistory_final$RecFY <- NULL




# save
save(Recov_FireHistory_final, file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/Recov_FireHist.RDATA")
write.csv(Recov_FireHistory_final, file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/Recov_FireHist.csv")

##########################################################################################################################################################
# DATAFRAME INTEGRATION
##########################################################################################################################################################
rm(list=ls())

setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery")

# load all variable dataframes
load("./Recov_FireHist.RDATA")
load("./Recov_DAYMET_082025.RDATA")
load("./Recov_Spec_082025.RDATA")
load(file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/Recov_Fire_Info_082025.RDATA")

# Merge DAYMET with spectral observations 
# note: X and Y projections do not currently match across dataframes but it's ok because point IDs will be used.
Recov_Spec$Obs_month <- as.numeric(Recov_Spec$Obs_month)
Recov_DAYMET$Obs_month <- as.numeric(Recov_DAYMET$Obs_month)
Recov_DAYMET$Obs_Year <- as.numeric(Recov_DAYMET$Obs_Year)

Recov_Master <- full_join(Recov_Spec, Recov_DAYMET, by=c("ptID",  "Obs_Year", "Obs_month")) %>% distinct()


# Merge FireHistory with Master
Recov_Master2 <- full_join(Recov_Master, Recov_FireHistory_final, by=c("ptID")) %>% select(-x, -y,-coords.x1  ,-coords.x2, -lat.x,-lon.x,  -lat.y,-lon.y) %>% distinct()

# correct number of pts?
length(unique(Recov_Master2$ptID))

# Save DF
save(Recov_Master2, file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/Recov_Master_082025.RDATA")