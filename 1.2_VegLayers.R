# Upland Vegetation Layers 
# M. Grace McLeod (2022)


# This script:
# 1. generates the upland vegetation layer for Everglades National Park and Big Cypress National Preserve
# 2. generates upland sample points

# Data downloaded from the Vegetation Mapping Project of Everglades National Park and Big Cypress National Preserve
  # EVER: https://irma.nps.gov/DataStore/Reference/Profile/2286556
  # eBICY: https://irma.nps.gov/DataStore/Reference/Profile/2288126
  # wBICY: https://irma.nps.gov/DataStore/Reference/Profile/2278515
# Geospatial data from original download (.mbd and .gdb files) was imported into ArcMap 10.8.1 and 
#     "EVER/BICY_VegMap_Vegetation_Dissolve" files were converted to shapefiles for analysis.
  

library(sf)
library(tidyverse)
library(ggplot2)
library(fasterize)
library(readr)
library(terra)

rm(list=ls())
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Veg_layers")

##########################################################################################################################################################
# 1. UPLAND VEGETATION LAYERS
##########################################################################################################################################################

# Products:
# "VegAll.shp" = Combined vegetation dataset for Everglades landscape
# "Uplands_all.RDATA"= Vegetation data for uplands only
# "Uplands_poly.shp" = Vegetation data for pinelends only
# "Uplands_raster.tif" = Raster of pinelands 


# COMBINE FOR MASTER VEGETATION DATASET.................................................................................................
# combine the two Big Cypress and the Everglades regions into one dataset

# IMPORTING SHAPEFILES
#Make a list of all shapeliness in that folder
shps <- dir("./EVG_shp_original", "*.shp$")
#import shapefiles
eBICY <- st_read(dsn = "./EVG_shp_original/", layer = "eBICY")
# make dataframe
eBICY_df <- as.data.frame(eBICY)
# wBICY
wBICY <- st_read(dsn = "./EVG_shp_original/", layer = "wBICY")
wBICY_df <- as.data.frame(wBICY)
# ENP
ENP <- st_read(dsn = "./EVG_shp_original", layer = "EVER")
ENP_df <- as.data.frame(ENP)


# CLEAN UP DATA AND MAKE IT MATCH FOR ALL REGIONS
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
ENP$Shape_Leng <- ENP$Shape_Area <- NA

ENP %>% names
eBICY %>% names
# recorder 
ENP <- ENP[, c(eBICY %>% names)]

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


# COMBINE FOR MASTER VEGETATION DATASET
# add column to each specifying which region it originally belonged to 
wBICY_reg <- "wBICY"
wBICY$Map_region <- wBICY_reg
eBICY_reg <- "eBICY"
eBICY$Map_region <- eBICY_reg
ENP_reg <- "EVER"
ENP$Map_region <- ENP_reg

eBICY %>% names
eBICY %>% names
ENP %>% names

# Combine using rbind function
VegAll <- rbind(eBICY, wBICY, ENP)
VegAll_df <- as.data.frame(VegAll)

# save
# Combined vegetation dataset for Everglades landscape
st_write(VegAll, ".", layer = "VegAll", overwrite_layer=TRUE, driver= "ESRI Shapefile")


# SUBSET TO UPLAND ECOSYSTEMS.................................................................................................

# ALL UPLANDS 
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

# save 
# Vegetation data for uplands only
save(Uplands, file = 'Uplands_all.RDATA')


# PINELANDS ONLY 
# Combine similar types identified as pineland 
  # subset to only desired L3 veg types
unique(Uplands$L3_name)
UpEco <- Uplands[which(Uplands$L3_name== "Pine Woodland"|
                         Uplands$L3_name== "Pine Upland"),]

UpEco %>% summary
# Assign a Pineland value
UpEco$EcoType[UpEco$L3_name =="Pine Woodland"] <- "Pineland"
UpEco$EcoType[UpEco$L3_name =="Pine Upland"] <- "Pineland"

# Create new column and assign 1 for pine
  # numeric value will be used for rasterizing
UpEco$value[UpEco$EcoType == "Pineland"] <- 1

# save as polygons for clipping/masking
# Vegetation data for pinelands only
writeOGR(UpEco, ".", layer = "Uplands_poly", overwrite_layer=TRUE, driver= "ESRI Shapefile")


# RASTERIZE .................................................................................................

# UPLANDS
# load reference raster (Created from Landsat image)
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/AOI/Landsat_ref/Raster30x30.tif") 

# make sure crs of dataframe and reference raster match
crs(Raster30x30) #+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs 
crs(UpEco) #+proj=utm +zone=17 +datum=NAD83 +units=m +no_defs 
UpEco <- st_transform(UpEco, "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")

# Convert to simple feature (sf) to run faserize()
UpEco_sf <- st_as_sf(UpEco)

# Fasterize: a faster version of the rasterize function   
Uplands_raster <- fasterize(UpEco_sf, Raster30x30, field= "value", fun="sum", background= NA) 
plot(Uplands_raster)
unique(Uplands_raster$layer)
  
# save raster
# pinelands only
writeRaster(Uplands_raster, "Uplands_raster.tif", overwrite=TRUE)


citation("raster")
citation("fasterize")


##########################################################################################################################################################
# 2. UPLAND SAMPLE PONTS FROM VEGETATION RASTER
##########################################################################################################################################################

# Products:
# "Sample_pts_upland" = sample points for pinelands


# CONVERT GRID USED FOR RASTERIZING INTO A SPATIAL POINT DATAFRAME.....................................................................................................

# Import grid (Upland Raster)
# Uplands
Uplands_raster <- rast("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Veg_layers/Uplands_raster.tif")

# Convert from raster grid to spatial points dataframe  
Sample_pts <- as.points(Uplands_raster)  %>% st_as_sf()# 414,124 pts
# assign point IDs
Sample_pts$ptID <- 1:414124 
# Plot  
plot(Sample_pts)

# Need to add Ecotype into the sample points:
UpEco.proj <- st_transform(UpEco, st_crs(Sample_pts) ) %>% select(EcoType) # Ecotype

UpEco.proj_raster <- rasterize(UpEco.proj, Uplands_raster, field= "EcoType", background= NA)

#Sample_pts_eco <- st_intersections( Sample_pts, UpEco.proj) # Extract information

EcoType.xy <- terra::extract( UpEco.proj_raster, Sample_pts, xy=T)

Sample_pts_eco <- cbind(Sample_pts,EcoType.xy[, c(2:4)]) %>% rename(coords.x1 = x , coords.x2 = y )

# SAVE UPLAND POINTS
# pinelands only
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Sampling")

st_write(Sample_pts_eco, dsn=".", layer="Sample_pts_upland_052925", driver="ESRI Shapefile", overwrite=TRUE)









