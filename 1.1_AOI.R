# AREA OF INTEREST, RASTER GRID, and LANDSCAPE BOUNDARY
# M. Grace McLeod (2023)


# This script:
# Produces the raster grid and Everglades boundary files used throughout this project.

# 1. Uses Landsat images downloaded from Earth Explorer and an "area of interest" (AOI) 
#     shapefile drawn in Google Earth Pro to generate the 30m raster grid ("Raster30x30.tif")

# 2. Uses the management boundaries of Everglades National Park and Big Cypress National Preserve 
#     to make a coumbined landscape boundary shapefile ("EVG_bound.shp")


# load libraries
library(terra)
library(sf)
library(cleangeo)

# clearn Global Environment 
rm(list=ls())
# set working directory
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/AOI")


# 1: AOI AND 30X30m GRID TO MATCH LANDSAT RESOLUTION ...........................................................................................................

  ## CREATING AN AREA OF INTEREST (AOI) ##
  # create AIO to cover entire extent of all shapefiels. Use Google Earth Pro.
  AOI <- st_read("EVG_AOI.kml")
  AOI <- st_transform(AOI, crs("+proj=longlat +datum=WGS84 +no_defs")) #must be in this crs for APPEARS
  
  st_write(AOI, dsn= ".", layer= "AOI", driver="ESRI Shapefile")
  save(AOI, file = "AOI.sp")
  
  ## CREATING 30X30m RASTER FROM LANDSAT ##
  # EarthExplorer-download 1 date of Landsat area within the AIO
  # import raster images (.tif)
  r1 <- rast("./Landsat_ref/LC08_L1TP_015042_20200307_20200822_02_T1_refl.tif")
  r2 <- rast("./Landsat_ref/LC08_L1TP_015043_20200307_20200822_02_T1_refl.tif")
  r3 <- rast("./Landsat_ref/LC08_L1TP_016042_20200314_20200822_02_T1_refl.tif")
  # Merge images together using Moasic()
  Raster30x30 <- mosaic(r1, r2, r3, fun = mean)
  # Save 
  save(Raster30x30, file="Raster30x30.tif", overwrite=T)
  

# 2: LANDSCAPE BOUNDARY SHAPE FILES...........................................................................................................

  # Open National Park boundary shapefiles downloaded from 
  # https://public-nps.opendata.arcgis.com/datasets/nps-boundary-1/explore?location=28.799880%2C-76.293142%2C6.60
  NPSbound <- st_read("./NPS_boundary/nps_boundary.shp")
  
  # Select only Everglades National Park (EVER) and Big Cypress National Preserve (BICY)
  EVGbound <- NPSbound[which(NPSbound$UNIT_CODE== "BICY" |
                               NPSbound$UNIT_CODE =="EVER"),]
  # check it out
  plot(EVGbound$geometry)
  # assign crs
  EVGbound <- st_transform(EVGbound, crs("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs "))
  # clean geometry...just incase
  EVGbound <- sf::st_make_valid(EVGbound, errors.only = NULL, strategy = "POLYGONATION", verbose = FALSE)
  # Save EVG bound as shp
  st_write(EVGbound, ".", layer = "EVG_bound", overwrite_layer=TRUE, driver= "ESRI Shapefile")
  
  
  
  