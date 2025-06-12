# FIRE HISTORY LAYERS
# S.Malone

# Preprocessing:
# This fire history raster layers created in this workflow: https://github.com/Malone-Disturbance-Ecology-Lab/Everglades-Fire-History

# This script: 
# 1. Extracts raster data to upland sample pts to generate fire history dataframes


library(sf)
library(tidyverse)
library(terra)

rm(list=ls())
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")

# Point to FireHistory folder on Malone Lab server
firehist_folder <- file.path("/", "Volumes", "malonelab", "Research", "ENP", "ENP Fire", "FireHistory")

# Read in tidy annual burned/unburned rasters
burned_rasters <- terra::rast(file.path(firehist_folder, "EVER_BICY_1978_2023_burned.tif")) 

# Read in tidy year of fire occurrence rasters
year_rasters <- terra::rast(file.path(firehist_folder, "EVER_BICY_1978_2023_year_occurrence.tif"))

# load upland sample pts
Sample_pts_upland <- sf::st_read(dsn="./Sampling", layer="Sample_pts_upland_052925")


# FIRE FREQUENCY (for total fires)  ####
names(burned_rasters) 

# extract to sample pts
FireHistory <- terra::extract(burned_rasters, Sample_pts_upland, method= "simple",
                              ID=TRUE,xy=TRUE )

FireHistory_df <- as.data.frame(FireHistory)
# calculate total fires 
FireHistory_df %>% names
FireHistory_df <- FireHistory_df %>%
  # rename ID column to ptID
  dplyr::rename(ptID = ID)
FireHistory_df[is.na(FireHistory_df)] <- 0 # turn NAs to 0

# total fire history
FireHistory_df$freq_1978_2020 <- rowSums(FireHistory_df[ , c(2:44)])
summary(FireHistory_df$freq_1978_2020)
names(FireHistory_df)

# recovery window 
# leave 2000 so there is room to calculate severity/see pre-fire spectral values
FireHistory_df$freq_2001_2007 <- rowSums(FireHistory_df[ , c(25:31)])
summary(FireHistory_df$freq_2001_2007)

# BL sample period
FireHistory_df$freq_2010_2020 <- rowSums(FireHistory_df[ , c(34:44)])
summary(FireHistory_df$freq_2010_2020)

write.csv(FireHistory_df, "/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Fire_History/FireHistory_df.csv")
save(FireHistory_df, file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Fire_History/FireHistory_df.RDATA")


# FIRE YEAR (for time since fire ####
FireYears <- terra::extract(year_rasters, Sample_pts_upland, method= "simple", ID=TRUE,xy=TRUE )

FireYears_df <- as.data.frame(FireYears)
FireYears_df <- FireYears_df %>%
  # rename ID column to ptID
  dplyr::rename(ptID = ID)
summary(FireYears_df)


write.csv(FireYears_df, "/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Fire_History/FireYears_df.csv")
save(FireYears_df, file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Fire_History/FireYears_df.RDATA")

