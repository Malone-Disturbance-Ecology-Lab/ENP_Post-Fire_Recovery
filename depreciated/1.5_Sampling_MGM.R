#  SAMPLING DESIGN
# M.Grace McLeod (2023)

# Sample points are filtered to remove unburned locations
# Points are selected for Baseline and Recovery sets

# * "Bonus Pts" should be deleted prior to publishing


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
library(varhandle)
library(tools)


setwd("/Volumes/inwedata/Malone Lab/ENP Fire/Grace_McLeod")

# load fire history data 
load("./Fire_History/FireHistory_df.RDATA")
FireHistory_df <- read_csv("./Fire_History/FireHistory_df.csv")

##########################################################################################################################################################
# SUBSET SAMPLE POINTS
##########################################################################################################################################################

# BURNED POINTS: points that have burned at least once on record (1978-2020)
# remove points that never burned (1978-2020)
Burned_pts_df <- FireHistory_df[which(FireHistory_df$freq_1978_2020 != 0),]        #  411,265 pts (2,859 unburned / 414,124 origninal = >99% have burned on reccord)
# select only columns of interest
Burned_smpl_pts <- Burned_pts_df[,c("ptID", "coords.x1", "coords.x2" )]
# convert back to spatial
coordinates(Burned_smpl_pts) <- ~ coords.x1 + coords.x2
# re-assign projection
proj4string(Burned_smpl_pts) <- CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")
crs(Burned_smpl_pts)
# save
writeOGR(Burned_smpl_pts, dsn="./Sampling", layer="Burned_smpl_pts", driver="ESRI Shapefile", overwrite=TRUE)


# RECOVERY POINTS: pts that burned once during recovery window (2001-2007)
Recov_pts_df <- Burned_pts_df[which(Burned_pts_df$freq_2001_2007 == 1),]          # 157,323 pts  
# convert to spatial
Recov_smpl_pts <- Recov_pts_df[,c("ptID",  "coords.x1", "coords.x2" )]
coordinates(Recov_smpl_pts) <- ~ coords.x1 + coords.x2
# re-assign projection
proj4string(Recov_smpl_pts) <- CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")
crs(Recov_smpl_pts)
# save 
writeOGR(Recov_smpl_pts, dsn="./Sampling", layer="Recov_smpl_pts", driver="ESRI Shapefile", overwrite=TRUE)


# BASELINE POINTS: pts that remain unburned from 2010-2020
leftovers <- anti_join(Burned_pts_df, Recov_pts_df) # do not want to repeat recov pts
BL_pts_df <- leftovers[which(leftovers$freq_2010_2020 == 0),]                       # 14,981 
# convert to spatial
BL_smpl_pts <- BL_pts_df[,c("ptID", "coords.x1", "coords.x2" )]
coordinates(BL_smpl_pts) <- ~ coords.x1 + coords.x2
# re-assign projection
proj4string(BL_smpl_pts) <- CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")
crs(BL_smpl_pts)
# save 
writeOGR(BL_smpl_pts, dsn="./Sampling", layer="BL_smpl_pts", driver="ESRI Shapefile", overwrite=TRUE)






