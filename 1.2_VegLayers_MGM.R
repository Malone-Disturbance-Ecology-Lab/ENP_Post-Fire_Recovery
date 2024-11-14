# Upland Vegetation Layers 
# M. Grace McLeod (2022)


# This script:
# 1. generates the upland vegetation layer for Everglades National Park and Big Cypress National Preserve
# 2. generates upland sample points

# Data downloaded from the Vegetation Mapping Project of Everglades National Park and Big Cypress National Preserve
  # EVER: https://irma.nps.gov/DataStore/Reference/Profile/2286556
  # eBICY: https://irma.nps.gov/DataStore/Reference/Profile/2288126
  # wBICY: https://irma.nps.gov/DataStore/Reference/Profile/2278515
# Geospatial data from original download (.mbd and .gdb files) was imported into ArcMap 10.8.1 and "EVER/BICY_VegMap_Vegetation_Dissolve" files were converted to shapefiles for analysis.
  
##########################################################################################################################################################

library(sf)
library(rgdal)
library(sp)
library(dplyr)
library(tidyverse)
library(raster)
library(ggplot2)
library(fasterize)
library(readr)

rm(list=ls())
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Veg_layers")

##########################################################################################################################################################
# UPLAND VEGETATION LAYERS
##########################################################################################################################################################

# IMPORTING SHAPEFILES.................................................................................................

#Make a list of all shapeliness in that folder
shps <- dir("./EVG_shp_original", "*.shp$")

#import shapefiles
eBICY <- readOGR(dsn = ".", layer = "eBICY")
# make dataframe
eBICY_df <- as.data.frame(eBICY)
# wBICY
wBICY <- readOGR(dsn = ".", layer = "wBICY")
wBICY_df <- as.data.frame(wBICY)
# ENP
ENP <- readOGR(dsn = ".", layer = "EVER")
ENP_df <- as.data.frame(ENP)


# CLEAN UP DATA AND MAKE IT MATCH FOR ALL REGIONS.................................................................................................

# check column names
colnames(ENP_df)
colnames(eBICY_df)
colnames(wBICY_df)

# ENP: rename columns
ENP$L1_name <- ENP$L1_Name ;  ENP$L1_Name <- NULL
ENP$L2_name <- ENP$L2_Name ;  ENP$L2_Name <- NULL
ENP$L3_name <- ENP$L3_Name ;  ENP$L3_Name <- NULL
ENP$L4_name <- ENP$L4_Name ;  ENP$L4_Name <- NULL
ENP$L5_name <- ENP$L5_Name ;  ENP$L5_Name <- NULL
ENP$L6_name <- ENP$L6_Name ;  ENP$L6_Name <- NULL
ENP$L7_name <- ENP$L7_Name ;  ENP$L7_Name <- NULL
ENP$VegCode_N <- ENP$NAME ;  ENP$NAME <- NULL
ENP$VegCode_L <- ENP$VegCode_Le ;  ENP$VegCode_Le <- NULL
# recorder 
ENP <- ENP[, c(1,2,3,14,13,6,7,8,9,10,11,12,4,5)]
# update dataframe and check names
ENP_df <- as.data.frame(ENP)

# wBICY: add L7_name column 
L7_NA <- "Undetermined"
wBICY$L7_name <- L7_NA
# rename "SHAPE.." to "Shape..."
wBICY$Shape_Area <- wBICY$SHAPE_Area ;  wBICY$SHAPE_Area <- NULL
wBICY$Shape_Leng <- wBICY$SHAPE_Leng ;  wBICY$SHAPE_Leng <- NULL
# update dataframe and check names
wBICY_df <- as.data.frame(wBICY)

#check crs'
crs(eBICY) #+proj=utm +zone=17 +datum=NAD83 +units=m +no_defs
crs(wBICY)
crs(ENP)

# COMBINE FOR MASTER VEGETATION DATASET.................................................................................................

# add column to each specifying which region it originally belonged to 
wBICY_reg <- "wBICY"
wBICY$Map_region <- wBICY_reg
eBICY_reg <- "eBICY"
eBICY$Map_region <- eBICY_reg
ENP_reg <- "EVER"
ENP$Map_region <- ENP_reg

# Combine using rbind function
VegAll <- rbind(eBICY, wBICY, ENP)
VegAll_df <- as.data.frame(VegAll)

# save as .RDATA
writeOGR(VegAll, ".", layer = "VegAll", overwrite_layer=TRUE, driver= "ESRI Shapefile")


# SUBSET TO UPLAND ECOSYSTEMS.................................................................................................

# get all unique L2 values
unique(VegAll_df$L2_name)

# subset for just uplands
Uplands <- VegAll[which(VegAll$L2_name== "Upland Shrubland"| 
                          VegAll$L2_name== "Upland Scrub"|
                          VegAll$L2_name== "Upland Forest"|
                          VegAll$L2_name== "Upland Woodland"),]
#check it
Uplands_df <- as.data.frame(Uplands)
plot(Uplands)

# save uplands dataset
save(Uplands, file = 'Uplands_all.RDATA')

# Combine similar types identified as pineland 
  # subset to only desired L3 veg types
unique(Uplands$L3_name)
UpEco <- Uplands[which(Uplands$L3_name== "Pine Woodland"|
                         Uplands$L3_name== "Pine Upland"),]
# Assign a Pineland value
UpEco$EcoType[UpEco$L3_name =="Pine Woodland"] <- "Pineland"
UpEco$EcoType[UpEco$L3_name =="Pine Upland"] <- "Pineland"
# Create new column and assign 1 for pine
  # numeric value will be used for rasterizing
UpEco$value[UpEco$EcoType == "Pineland"] <- 1

# save polygons for clipping/masking
writeOGR(UpEco, ".", layer = "Uplands_poly", overwrite_layer=TRUE, driver= "ESRI Shapefile")


# RASTERIZE .................................................................................................

# UPLANDS
# load reference raster (Created from Landsat image)
load("/Volumes/inwedata/Malone Lab/ENP Fire/Grace_McLeod/AOI/Landsat_ref/Raster30x30.tif") 

# make sure crs of dataframe and reference raster match
crs(Raster30x30) #+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs 
crs(UpEco) #+proj=utm +zone=17 +datum=NAD83 +units=m +no_defs 
UpEco <- spTransform(UpEco, "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")

# Convert to simple feature (sf) to run faserize()
UpEco_sf <- st_as_sf(UpEco)

# Fasterize: a faster version of the rasterize function   
Uplands_raster <- fasterize(UpEco_sf, Raster30x30, field= "value", fun="sum", background= NA) 
plot(Uplands_raster)
unique(Uplands_raster$layer)
  
# save raster
writeRaster(Uplands_raster, "Uplands_raster.tif", overwrite=TRUE)


citation("raster")
citation("fasterize")


##########################################################################################################################################################
# UPLAND SAMPLE PONTS FROM VEGETATION RASTER
##########################################################################################################################################################

# CONVERT GRID USED FOR RASTERIZING INTO A SPATIAL POINT DATAFRAME.....................................................................................................
# Import grid (Upland Raster)
# Uplands
Uplands_raster <- raster("/Volumes/inwedata/Malone Lab/ENP Fire/Grace_McLeod/Veg_layers/Uplands_raster.tif")
# Convert from raster grid to spatial points dataframe  
Sample_pts <- rasterToPoints(Uplands_raster, fun=NULL, spatial=TRUE)  # 414,124 pts
# assign point IDs
Sample_pts$ptID <- 1:414124 
# Plot  
plot(Sample_pts)
# SAVE UPLAND POINTS
setwd("/Volumes/inwedata/Malone Lab/ENP Fire/Grace_McLeod/Sampling")
writeOGR(Sample_pts, dsn=".", layer="Sample_pts_upland", driver="ESRI Shapefile", overwrite=TRUE)



##########################################################################################################################################################
# SUMMARY 
##########################################################################################################################################################
setwd("/Volumes/inwedata/Malone Lab/ENP Fire/Grace_McLeod")

# load sample pts
smpl_pts <- st_read("/Sampling/Sample_pts_upland.shp")
crs(smpl_pts)
# load all uplands
load("./Veg_layers/Uplands_all.RDATA")


# subset to pinelands
UpEco <- Uplands[which(Uplands$L3_name== "Pine Woodland"| +
                         Uplands$L3_name== "Pine Upland"),]
# convert to sf
projcrs <- "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"
polygons <- st_as_sf(x = UpEco,                         
               coords = c("lon", "lat"),
               crs = projcrs)
points <- st_as_sf(x = smpl_pts,                         
                coords = c("lon", "lat"),
                crs = "NAD83 / UTM zone 17N")
# match crs
polygons <-  st_transform(df, crs = st_crs(points))
st_crs(polygons) <-"+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"
polygons <- st_as_sf(x = polygons,                         
               coords = c("lon", "lat"),
               crs = projcrs)
# check overlap
plot(st_geometry(polygons[1,]), axes=T)
plot(st_geometry(points), pch=1, col = 'blue', add = TRUE)

# write out simple features (join in GIS)
st_write(polygons, "./Veg_layers", layer ="pineland_polys_sf.shp", driver="ESRI Shapefile")
st_write(points, "./Sampling", layer ="upland_pts_sf.shp", driver="ESRI Shapefile")

# .........................................
# GIS: used Intersect tool to combine files
#..........................................

# load veg_all shapefiles
veg_info <- st_read("./Veg_layers/veg_info/Veg_info_intersect.shp")
veg_info <- st_as_sf(x = veg_info,                         
                     coords = c("lon", "lat"),
                     crs = projcrs)








