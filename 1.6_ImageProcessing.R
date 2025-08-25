# LANDSAT IMAGE PROCESSING
# M.Grace McLeod (2022)


# this script does the following steps to process images for 4 ARD Landsat tiles covering the Everglades

# Images were downloaded manually from Earth Explorer:  https://earthexplorer.usgs.gov/

# 1. open tar files and saves Landsat .tif image 
# 2. extracts spectral data from tile-specific stacks to upland sample points and filters data for clouds and QAQC
# 3. merges all tile-specific dataframes to make a Spectral Master dataframe and calculates spectral indicies ("Spec_Master.RDATA")


rm(list=ls())

library(landsat)
library(sp)
library(fields)
library(RGISTools)
library(remotes)
library(tools)
library(stringr)
library(tidyverse)
library(terra)

setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Image_Processing")

##########################################################################################################################################################
# 1. OPENING TAR FILES. *** run once***
##########################################################################################################################################################

# BULK DOWNLOAD FROM EARTH EXPLORER
# download completed in chunks. Run on each 
bulk <- "./BulkDownload/BulkDownLoad_2000_2005"
#bulk <- "./BulkDownload/BulkDownLoad_2006_2009"
#bulk <- "./BulkDownload/BulkDownLoad_2010_2012"
#bulk <- "./BulkDownload/BulkDownLoad_2013_2015"
#bulk <- "./BulkDownload/BulkDownLoad_2016_2020"

# creat path to store tif files
LStif <- "./LS_tif"

# SELECT .TAR FILES
# make a list of all the .tar files
all_tars <- list.files(bulk, pattern="tar") 
# just the surface reflectance (SR) ones
SR.tar <- list.files(bulk, pattern = glob2rx("*_SR*.tar$"), full.names = TRUE)

# LOOP TO UNTAR FILES 

for (tar in SR.tar){
  setwd(bulk)
  print(tar)
  LStifs <- untar(tar, list=FALSE, exdir = LStif) 
}

##########################################################################################################################################################
# 2. IMAGE PROCESSING
##########################################################################################################################################################

# SET UP THE ENVIRONMENT  
# Pull out tifs from new directory 
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Image_Processing") 
LStif <- "./LS_tif"
# make a list of all the .tif files
all_tifs <- list.files(LStif, pattern = "TIF") 
# Load upland sample points 
#setwd("/Volumes/inwedata/Malone Lab/ENP Fire/Grace_McLeod/Sampling")
Smpl_pts <- sf::st_read(dsn = "/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Sampling", layer = "Sample_pts_upland")
# set crs to match landsat
Smpl_pts <- st_transform(Smpl_pts, crs("+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))


# TILE 26_18.............................................................................................................................................................
# Group by tile
tile_26_18 <- list.files(LStif, pattern = glob2rx("*_026018*.TIF$"), full.names = TRUE)
# create stack
system.time(stack_26_18 <- stack(tile_26_18))
# subset sample points to tile
test1 <- subset(stack_26_18, 1)
Sub_pts_2618 <- st_as_sf(crop(vect(Smpl_pts),  ext(test1)))
# extract to tile sample points
system.time(Ext_26_18_sp <- terra::extract(stack_26_18, Sub_pts_2618, method= "simple", buffer=NULL, df=TRUE, sp=TRUE, factors=TRUE))
Ext_26_18_df <- as.data.frame (Ext_26_18_sp)

# SUBSET DATAFRAME (to avoid exhausting vector memory)
# split into chunks *10,000 at a time seemed to go well* 
    Sub_26_18 <- Ext_26_18_df[1:13,]
    # Sub_26_18 <- Ext_26_18_df[7893:15785,]
# melt df longways
melt_26_18 <- melt(Sub_26_18, na.rm=FALSE, id=c("ptID", "EcoType", "coords.x1", "coords.x2"))
# split names
class( melt_26_18$variable)
melt_26_18$variable <- as.character( melt_26_18$variable)
melt_26_18 [c('LS7', 'CU', 'Tile', 'Obs_Date', "DLdate", 'colect', 'type', 'Band')] <- str_split_fixed(melt_26_18$variable, '_', 8)
# get rid of unnecessary columns 
melt_26_18$variable <- NULL ; melt_26_18$LS7 <- NULL ; melt_26_18$CU <- NULL ; melt_26_18$DLdate <- NULL ;  melt_26_18$colect <- NULL ;  melt_26_18$type <- NULL 
# turn band into columns
Spec_26_18 <- melt_26_18 %>%
  spread(Band, value)
# remove rows where band values =NA 
# if value is NA for one band, it will be for all.
Spec_26_18_noNA <- Spec_26_18 %>% drop_na(c('B1', 'B2', 'B3', 'B4', 'B5', 'B7'))
# QAQC
# make all flagged values = NA
Spec_26_18_qaqc <- Spec_26_18_noNA
Spec_26_18_qaqc$CLOUD_QA[Spec_26_18_qaqc$CLOUD_QA > 0] <- NA
Spec_26_18_qaqc$PIXEL[Spec_26_18_qaqc$PIXEL != 5440] <- NA  # 5440 is code for totally clear pixels
Spec_26_18_qaqc$RADSAT[Spec_26_18_qaqc$RADSAT > 0] <- NA
# remove rows with flagged values
Spec_26_18_clean <- Spec_26_18_qaqc %>% drop_na(c('CLOUD_QA', 'PIXEL', 'RADSAT'))
Spec_26_18_clean$RADSAT <- NULL ; Spec_26_18_clean$QA <- NULL ; Spec_26_18_clean$PIXEL <- NULL ; Spec_26_18_clean$LINEAGE <- NULL ; Spec_26_18_clean$CLOUD_QA <- NULL ; Spec_26_18_clean$ATMOS_OPACITY <- NULL
# rename this subset and clear environment before running next subset
    # Spec_26_18_a <- Spec_26_18_clean
    # Spec_26_18_b <- Spec_26_18_clean
rm(melt_26_18, Spec_26_18, Spec_26_18_clean, Spec_26_18_noNA, Spec_26_18_qaqc, Sub_26_18)
# merge abcd
Spec_26_18 <- rbind(Spec_26_18_a, Spec_26_18_b)

# save
setwd("./Master_Spec")
save (stack_26_18, Ext_26_18_df, Spec_26_18, file="Tile_26_18.RDATA")
# clean up envr before next tile
rm(stack_26_18, Ext_26_18_df, Sub_pts_2618, Ext_26_18_sp, Spec_26_18_a, Spec_26_18_b) 



# TILE 26_19.............................................................................................................................................................
# Group by tile
tile_26_19 <- list.files(LStif, pattern = glob2rx("*_026019*.TIF$"), full.names = TRUE)
# create stack
system.time(stack_26_19 <- stack(tile_26_19))
# subset sample points to tile
test1 <- subset(stack_26_19, 1)
Sub_pts_2619 <- st_as_sf(crop(vect(Smpl_pts),  ext(test1)))
# extract to tile sample points
system.time(Ext_26_19_sp <- raster::extract(stack_26_19, Sub_pts_2619, method= "simple", buffer=NULL, df=TRUE, sp=TRUE, factors=TRUE))
Ext_26_19_df <- as.data.frame (Ext_26_19_sp)

# SUBSET DATAFRAME (to avoid exhausting vector memory)
# split into chunks (42168/4= 10542)
      # Sub_26_19 <- Ext_26_19_df[1:10542,]
      # Sub_26_19 <- Ext_26_19_df[10543:21084,]
      # Sub_26_19 <- Ext_26_19_df[21085:31626,]
      # Sub_26_19 <- Ext_26_19_df[31627:42168,]
# melt df longways
melt_26_19 <- melt(Sub_26_19, na.rm=FALSE, id=c("ptID", "EcoType", "coords.x1", "coords.x2"))
# split names
class( melt_26_19$variable)
melt_26_19$variable <- as.character( melt_26_19$variable)
melt_26_19 [c('LS7', 'CU', 'Tile', 'Obs_Date', "DLdate", 'colect', 'type', 'Band')] <- str_split_fixed(melt_26_19$variable, '_', 8)
# get rid of unnecessary columns 
melt_26_19$variable <- NULL ; melt_26_19$LS7 <- NULL ; melt_26_19$CU <- NULL ; melt_26_19$DLdate <- NULL ;  melt_26_19$colect <- NULL ;  melt_26_19$type <- NULL 
# turn band into columns
Spec_26_19 <- melt_26_19 %>%
  spread(Band, value)
# remove rows where band values =NA 
# if value is NA for one band, it will be for all.
Spec_26_19_noNA <- Spec_26_19 %>% drop_na(c('B1', 'B2', 'B3', 'B4', 'B5', 'B7'))
# QAQC
# make all flagged values = NA
Spec_26_19_qaqc <- Spec_26_19_noNA
Spec_26_19_qaqc$CLOUD_QA[Spec_26_19_qaqc$CLOUD_QA > 0] <- NA
Spec_26_19_qaqc$PIXEL[Spec_26_19_qaqc$PIXEL != 5440] <- NA  # 5440 is code for totally clear pixels
Spec_26_19_qaqc$RADSAT[Spec_26_19_qaqc$RADSAT > 0] <- NA
# remove rows with flagged values
Spec_26_19_clean <- Spec_26_19_qaqc %>% drop_na(c('CLOUD_QA', 'PIXEL', 'RADSAT'))
Spec_26_19_clean$RADSAT <- NULL ; Spec_26_19_clean$QA <- NULL ; Spec_26_19_clean$PIXEL <- NULL ; Spec_26_19_clean$LINEAGE <- NULL ; Spec_26_19_clean$CLOUD_QA <- NULL ; Spec_26_19_clean$ATMOS_OPACITY <- NULL
# rename this subset and clear environment before running next subset
     # Spec_26_19_a <- Spec_26_19_clean
     # Spec_26_19_b <- Spec_26_19_clean
     # Spec_26_19_c <- Spec_26_19_clean
     # Spec_26_19_d <- Spec_26_19_clean
rm(melt_26_19, Spec_26_19, Spec_26_19_clean, Spec_26_19_noNA, Spec_26_19_qaqc, Sub_26_19)
# merge abcd
Spec_26_19 <- rbind(Spec_26_19_a, Spec_26_19_b)
Spec_26_19 <- rbind(Spec_26_19, Spec_26_19_c)
Spec_26_19 <- rbind(Spec_26_19, Spec_26_19_d)

# save
#setwd("/Volumes/inwedata/Malone Lab/ENP Fire/Grace_McLeod/Image_Processing/Master_Spec")
save (stack_26_19, Ext_26_19_df, Spec_26_19, file="Tile_26_19.RDATA")
# clean up envr before next tile
rm(stack_26_19, Sub_pts_2619, Ext_26_19_sp, Ext_26_19_df, Spec_26_19_a, Spec_26_19_b, Spec_26_19_c, Spec_26_19_d)


# TILE 27_18.............................................................................................................................................................
# Group by tile
tile_27_18 <- list.files(LStif, pattern = glob2rx("*_027018*.TIF$"), full.names = TRUE)
# create stack
system.time(stack_27_18 <- stack(tile_27_18))
# subset sample points to tile
test1 <- subset(stack_27_18, 1)
Sub_pts_2718 <- st_as_sf(crop(vect(Smpl_pts),  ext(test1)))
# extract to tile sample points
system.time(Ext_27_18_sp <- raster::extract(stack_27_18, Sub_pts_2718, method= "simple", buffer=NULL, df=TRUE, sp=TRUE, factors=TRUE))
Ext_27_18_df <- as.data.frame (Ext_27_18_sp)

# SUBSET DATAFRAME (to avoid exhausting vector memory)
# split into chunks (80131/8= 10016.38)
    # Sub_27_18 <- Ext_27_18_df[1:10016,]
    # Sub_27_18 <- Ext_27_18_df[10017:20032,]
    # Sub_27_18 <- Ext_27_18_df[20033:30048,]
    # Sub_27_18 <- Ext_27_18_df[30049:40064,]
    # Sub_27_18 <- Ext_27_18_df[40065:50080,]
    # Sub_27_18 <- Ext_27_18_df[50081:60096,]
    # Sub_27_18 <- Ext_27_18_df[60097:70112,]
    # Sub_27_18 <- Ext_27_18_df[70113:80131,]
# melt df longways
melt_27_18 <- melt(Sub_27_18, na.rm=FALSE, id=c("ptID", "EcoType", "coords.x1", "coords.x2"))
# split names
class( melt_27_18$variable)
melt_27_18$variable <- as.character( melt_27_18$variable)
melt_27_18 [c('LS7', 'CU', 'Tile', 'Obs_Date', "DLdate", 'colect', 'type', 'Band')] <- str_split_fixed(melt_27_18$variable, '_', 8)
# get rid of unnecessary columns 
melt_27_18$variable <- NULL ; melt_27_18$LS7 <- NULL ; melt_27_18$CU <- NULL ; melt_27_18$DLdate <- NULL ;  melt_27_18$colect <- NULL ;  melt_27_18$type <- NULL 
# turn band into columns
Spec_27_18 <- melt_27_18 %>%
  spread(Band, value)
# remove rows where band values =NA 
# if value is NA for one band, it will be for all.
Spec_27_18_noNA <- Spec_27_18 %>% drop_na(c('B1', 'B2', 'B3', 'B4', 'B5', 'B7'))
# QAQC
# make all flagged values = NA
Spec_27_18_qaqc <- Spec_27_18_noNA
Spec_27_18_qaqc$CLOUD_QA[Spec_27_18_qaqc$CLOUD_QA > 0] <- NA
Spec_27_18_qaqc$PIXEL[Spec_27_18_qaqc$PIXEL != 5440] <- NA  # 5440 is code for totally clear pixels
Spec_27_18_qaqc$RADSAT[Spec_27_18_qaqc$RADSAT > 0] <- NA
# remove rows with flagged values
Spec_27_18_clean <- Spec_27_18_qaqc %>% drop_na(c('CLOUD_QA', 'PIXEL', 'RADSAT'))
Spec_27_18_clean$RADSAT <- NULL ; Spec_27_18_clean$QA <- NULL ; Spec_27_18_clean$PIXEL <- NULL ; Spec_27_18_clean$LINEAGE <- NULL ; Spec_27_18_clean$CLOUD_QA <- NULL ; Spec_27_18_clean$ATMOS_OPACITY <- NULL
# rename this subset and clear environment before running next subset
    # Spec_27_18_a <- Spec_27_18_clean
    # Spec_27_18_b <- Spec_27_18_clean
    # Spec_27_18_c <- Spec_27_18_clean
    # Spec_27_18_d <- Spec_27_18_clean
    # Spec_27_18_e <- Spec_27_18_clean
    # Spec_27_18_f <- Spec_27_18_clean
    # Spec_27_18_g <- Spec_27_18_clean
    # Spec_27_18_h <- Spec_27_18_clean
rm(melt_27_18, Spec_27_18, Spec_27_18_clean, Spec_27_18_noNA, Spec_27_18_qaqc, Sub_27_18)
# merge abcd
Spec_27_18 <- rbind(Spec_27_18_a, Spec_27_18_b)
Spec_27_18 <- rbind(Spec_27_18, Spec_27_18_c)
Spec_27_18 <- rbind(Spec_27_18, Spec_27_18_d)
Spec_27_18 <- rbind(Spec_27_18, Spec_27_18_e)
Spec_27_18 <- rbind(Spec_27_18, Spec_27_18_f)
Spec_27_18 <- rbind(Spec_27_18, Spec_27_18_g)
Spec_27_18 <- rbind(Spec_27_18, Spec_27_18_h)

# save
#setwd("/Volumes/inwedata/Malone Lab/ENP Fire/Grace_McLeod/Image_Processing/Master_Spec")
save (stack_27_18, Ext_27_18_df, Spec_27_18,  file="Tile_27_18.RDATA") 
# clean up envr before next tile
rm(stack_27_18,  Ext_27_18_df, Sub_pts_2718, Ext_27_18_sp, Spec_27_18_a, Spec_27_18_b, Spec_27_18_c, Spec_27_18_d, Spec_27_18_e, Spec_27_18_f, Spec_27_18_g, Spec_27_18_h) 


# TILE 27_19.............................................................................................................................................................
# Group by tile
tile_27_19 <- list.files(LStif, pattern = glob2rx("*_027019*.TIF$"), full.names = TRUE)
# create stack
system.time(stack_27_19 <- stack(tile_27_19))
# subset sample points to tile
test1 <- subset(stack_27_19, 1)
Sub_pts_2719 <- st_as_sf(crop(vect(Smpl_pts),  ext(test1)))
# extract to tile sample points
system.time(Ext_27_19_sp <- raster::extract(stack_27_19, Sub_pts_2719, method= "simple", buffer=NULL, df=TRUE, sp=TRUE, factors=TRUE))
Ext_27_19_df <- as.data.frame (Ext_27_19_sp)

# SUBSET DATAFRAME (to avoid exhausting vector memory)
#  ....this one is huge...need a new appraoch...
# create a blank Spec_27_29 dataframe that chunk dfs can be appended to
Spec_27_19 <- data.frame(matrix(ncol=12, nrow=0))
colnames(Spec_27_19) <- colnames(Spec_26_18)
# Its a prime number....do the first 410000 in chunks, then the remaining 3686....stupid I know but it should work
sub1 <- Ext_27_19_df[1:410000,]
      # just to do it faster on 2 computers 
     # sub1 <- Ext_27_19_df[1:210000,]
        #sub1<- Ext_27_19_df[210001:410000,]
# Make list of low ends of the chunk 
a <- seq(1, 400001, 10000)
      # a <- seq(40001, 200001, 10000 )
        #a <- seq(210001, 400001, 10000)
# Loop through sub1 
for (i in a) {
  print(i)
  # subset to the chunk
  b <- i+9999 
  chunk <- sub1[i:b,] 
  #  RUN NORMAL PROCESS
  # melt df longways
  melt_27_19 <- melt(chunk, na.rm=FALSE, id=c("ptID", "EcoType", "coords.x1", "coords.x2"))
  # split names
  class( melt_27_19$variable)
  melt_27_19$variable <- as.character( melt_27_19$variable)
  melt_27_19 [c('LS7', 'CU', 'Tile', 'Obs_Date', "DLdate", 'colect', 'type', 'Band')] <- str_split_fixed(melt_27_19$variable, '_', 8)
  # get rid of unnecessary columns 
  melt_27_19$variable <- NULL ; melt_27_19$LS7 <- NULL ; melt_27_19$CU <- NULL ; melt_27_19$DLdate <- NULL ;  melt_27_19$colect <- NULL ;  melt_27_19$type <- NULL 
  # turn band into columns
  bands_27_19 <- melt_27_19 %>%
    spread(Band, value)
  # remove rows where band values =NA 
  # if value is NA for one band, it will be for all.
  Spec_27_19_noNA <- bands_27_19 %>% drop_na(c('B1', 'B2', 'B3', 'B4', 'B5', 'B7'))
  # QAQC
  # make all flagged values = NA
  Spec_27_19_qaqc <- Spec_27_19_noNA
  Spec_27_19_qaqc$CLOUD_QA[Spec_27_19_qaqc$CLOUD_QA > 0] <- NA
  Spec_27_19_qaqc$PIXEL[Spec_27_19_qaqc$PIXEL != 5440] <- NA  # 5440 is code for totally clear pixels
  Spec_27_19_qaqc$RADSAT[Spec_27_19_qaqc$RADSAT > 0] <- NA
  # remove rows with flagged values
  Spec_27_19_clean <- Spec_27_19_qaqc %>% drop_na(c('CLOUD_QA', 'PIXEL', 'RADSAT'))
  Spec_27_19_clean$RADSAT <- NULL ; Spec_27_19_clean$QA <- NULL ; Spec_27_19_clean$PIXEL <- NULL ; Spec_27_19_clean$LINEAGE <- NULL ; Spec_27_19_clean$CLOUD_QA <- NULL ; Spec_27_19_clean$ATMOS_OPACITY <- NULL
  # append to Spec_27_19 df
  Spec_27_19 <- rbind(Spec_27_19, Spec_27_19_clean)
  # clean up 
  rm(melt_27_19, Spec_27_19_clean, Spec_27_19_noNA, Spec_27_19_qaqc, bands_27_19)
  
}

# save
#setwd("/Volumes/inwedata/Malone Lab/ENP Fire/Grace_McLeod/Image_Processing/Master_Spec")
save (Spec_27_19, file="Spec_27_19.RDATA") 


      
# Run on ramaining lines 
sub2 <- Ext_27_19_df[410001:413686,]
melt_27_19 <- melt(sub2, na.rm=FALSE, id=c("ptID", "EcoType", "coords.x1", "coords.x2"))
class( melt_27_19$variable)
melt_27_19$variable <- as.character( melt_27_19$variable)
melt_27_19 [c('LS7', 'CU', 'Tile', 'Obs_Date', "DLdate", 'colect', 'type', 'Band')] <- str_split_fixed(melt_27_19$variable, '_', 8)
melt_27_19$variable <- NULL ; melt_27_19$LS7 <- NULL ; melt_27_19$CU <- NULL ; melt_27_19$DLdate <- NULL ;  melt_27_19$colect <- NULL ;  melt_27_19$type <- NULL 
bands_27_19 <- melt_27_19 %>%
  spread(Band, value)
Spec_27_19_noNA <- bands_27_19 %>% drop_na(c('B1', 'B2', 'B3', 'B4', 'B5', 'B7'))
Spec_27_19_qaqc <- Spec_27_19_noNA
Spec_27_19_qaqc$CLOUD_QA[Spec_27_19_qaqc$CLOUD_QA > 0] <- NA
Spec_27_19_qaqc$PIXEL[Spec_27_19_qaqc$PIXEL != 5440] <- NA  # 5440 is code for totally clear pixels
Spec_27_19_qaqc$RADSAT[Spec_27_19_qaqc$RADSAT > 0] <- NA
Spec_27_19_clean <- Spec_27_19_qaqc %>% drop_na(c('CLOUD_QA', 'PIXEL', 'RADSAT'))
Spec_27_19_clean$RADSAT <- NULL ; Spec_27_19_clean$QA <- NULL ; Spec_27_19_clean$PIXEL <- NULL ; Spec_27_19_clean$LINEAGE <- NULL ; Spec_27_19_clean$CLOUD_QA <- NULL ; Spec_27_19_clean$ATMOS_OPACITY <- NULL
# append to Spec_27_19
Spec_27_19 <- rbind(Spec_27_19, Spec_27_19_clean)
rm(melt_27_19, Spec_27_19_clean, Spec_27_19_noNA, Spec_27_19_qaqc, bands_27_19)

# save
#setwd("/Volumes/inwedata/Malone Lab/ENP Fire/Grace_McLeod/Image_Processing/Master_Spec")
save (stack_27_19, Ext_27_19_df, Spec_27_19, file="Tile_27_19.RDATA") 
# clean up envr 
rm(list=ls())


##########################################################################################################################################################
# 3. DATAFRAME INTEGRATION
##########################################################################################################################################################
# if wd not already set
setwd("/Volumes/inwedata/Malone Lab/ENP Fire/Grace_McLeod/Image_Processing/Master_Spec")
# load files
load("./Tile_26_18.RDATA")
load("./Tile_26_19.RDATA")
load("./Tile_27_18.RDATA")
load("./Tile_27_19.RDATA")

# bind tiles together
Spec_Master <- rbind(Spec_26_18, Spec_26_19)
Spec_Master <- rbind(Spec_Master , Spec_27_18)
Spec_Master <- rbind(Spec_Master , Spec_27_19)
# Save 
#setwd("/Volumes/inwedata/Malone Lab/ENP Fire/Grace_McLeod/Image_Processing/Master_Spec")
save(Spec_Master, file="Spec_Master.RDATA")

# CALCULATE SPECTRAL INDECIES
# Band combinations guide (https://www.esri.com/arcgis-blog/products/product/imagery/band-combinations-for-landsat-8/)
# turn off factors
options(stringsAsFactors = FALSE)

# NDVI 
# calculate NDVI using the red (band 3) and nir (band 4) bands  
# NVDI = (NIR - Red) / (NIR + Red) 
Spec_Master$NDVI <- (Spec_Master$B4 - Spec_Master$B3) / (Spec_Master$B4 + Spec_Master$B3)   
# NBR 
# using nir (band 4) and swir (band 7) *there is no band 6 so band 7 is in the 6th position
# NBR = (NIR - SWIR) / (NIR + SWIR)
Spec_Master$NBR <- (Spec_Master$B4 - Spec_Master$B7) / (Spec_Master$B4 + Spec_Master$B7)   
# NBR 2 
# using swir1 (band 5) and swir2 (band 7) *suposidly good for measuring recovery
# NBR2 = (SWIR1 - SWIR2) / (SWIR1 + SWIR2)
Spec_Master$NBR2 <- (Spec_Master$B5 - Spec_Master$B7) / (Spec_Master$B5 + Spec_Master$B7)  

# check it
head(Spec_Master)

# SAVE SPECTRAL MASTER DATAFRAME
#setwd("/Volumes/inwedata/Malone Lab/ENP Fire/Grace_McLeod/Image_Processing/Master_Spec")
write.csv(Spec_Master, file="Spec_Master.csv")
save(Spec_Master, file="Spec_Master.RDATA")










