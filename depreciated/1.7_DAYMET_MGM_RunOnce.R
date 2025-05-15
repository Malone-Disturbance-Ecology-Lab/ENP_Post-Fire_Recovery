# DAYMENT METEROROLOGICAL DATA
  # AUTHOR: M.GRACE MCLEOD (2022)


# This script:
  # 1. downloads monthly climate data from DAYMET 
  # 2. stacks monthly layers
  # 3. masks stacks to the Everglades boundary


rm(list=ls())

library(ncdf4)
library(raster)
library(rgdal)
library(ggplot2)
library(dplyr)  
library(knitr) 
library(markdown)
library(rmarkdown)
library(covr) 
library(testthat)
library(tools)
library(FedData)
library(daymetr) 
library(terra)

# ADDITIONAL RESOURCES
# https://daac.ornl.gov/DAYMET/guides/Daymet_V4_Monthly_Climatology.html
# https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1855
# https://cran.r-project.org/web/packages/daymetr/daymetr.pdf
# https://tmieno2.github.io/R-as-GIS-for-Economists/daymet-with-daymetr-and-feddata.html
# https://cran.r-project.org/web/packages/daymetr/vignettes/daymetr-vignette.html

setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Seasonal_Cond")

##########################################################################################################################################################
# IMPORTING DATA
##########################################################################################################################################################

# DOWNLOAD DAYMET GRIDED CLIMATE DATA..........................................................................................................
# Monthly precip
download_daymet_ncss( c(26.5, -81.75, 25, -80.25), 
                      param= 'prcp', 
                      start = 2000, end = 2020, 
                      frequency = "monthly", 
                      path = "./DAYMET_monthly_prcp")    
# Monthly tmax
download_daymet_ncss( c(26.5, -81.75, 25, -80.25), 
                      param= 'tmax', 
                      start = 2000, end = 2020, 
                      frequency = "monthly", 
                      path = "./DAYMET_monthly_tmax")    
# Monthly tmin
download_daymet_ncss( c(26.5, -81.75, 25, -80.25), 
                      param= 'tmin', 
                      start = 2000, end = 2020, 
                      frequency = "monthly", 
                      path = "./DAYMET_monthly_tmin") 


##########################################################################################################################################################
# DATA PROCESSING 
##########################################################################################################################################################

# LoadEVG boundary for masking
EVGbound <- readOGR("/Volumes/inwedata/Malone Lab/ENP Fire/Grace_McLeod/AOI/EVG_bound.shp")
# Change projection to match DAYMET data
EVGbound <- spTransform(EVGbound, CRSobj = CRS("+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs"))

# PRECIPITATION .........................................................................................................
# load files
precip.files <- dir("./DAYMET_monthly_prcp", "*_ncss.nc$")
precip.files<- list_files_with_exts("./DAYMET_monthly_prcp", "nc")
# stack
precip_stack <- stack(precip.files)
# change projections
raster::projection(precip_stack) <- "+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs"
# mask to EVGbound
precip_EVG <- mask(precip_stack, EVGbound)
# check it
sub1 <- subset(precip_EVG, 1)
plot(sub1) 
plot(EVGbound, add=TRUE)
# save 
#setwd("/Volumes/inwedata/Malone Lab/ENP Fire/Grace_McLeod/Seasonal_Cond")
writeRaster(precip_EVG, "precip_EVG.tif", overwrite=TRUE)


# MAX TEMPERATURE.........................................................................................................
# load files
tmax.files <- dir("./DAYMET_monthly_tmax", "*_ncss.nc$")
tmax.files<- list_files_with_exts("./DAYMET_monthly_tmax", "nc")
# stack
tmax_stack <- stack(tmax.files)
# change projections
raster::projection(tmax_stack) <- "+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs"
# mask to EVGbound
tmax_EVG <- mask(tmax_stack, EVGbound)
# check it
sub1 <- subset(tmax_EVG, 1)
plot(sub1) 
plot(EVGbound, add=TRUE)
# save 
#setwd("/Volumes/inwedata/Malone Lab/ENP Fire/Grace_McLeod/Seasonal_Cond")
writeRaster(tmax_EVG, "tmax_EVG.tif", overwrite=TRUE)


# MIN TEMPERATURE.........................................................................................................
# load files
tmin.files <- dir("./DAYMET_monthly_tmin", "*_ncss.nc$")
tmin.files<- list_files_with_exts("./DAYMET_monthly_tmin", "nc")
# stack
tmin_stack <- stack(tmin.files)
# change projections
raster::projection(tmin_stack) <- "+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs"
# mask to EVGbound
tmin_EVG <- mask(tmin_stack, EVGbound)
# check it
sub1 <- subset(tmin_EVG, 1)
plot(sub1) 
plot(EVGbound, add=TRUE)
# save 
#setwd("/Volumes/inwedata/Malone Lab/ENP Fire/Grace_McLeod/Seasonal_Cond")
writeRaster(tmin_EVG, "tmin_EVG.tif", overwrite=TRUE)


# SAVE .RDATA..............................................................................................................................
# setwd("/Volumes/inwedata/Malone Lab/ENP Fire/Grace_McLeod/Seasonal_Cond")
# save(precip_EVG, tmax_EVG, tmin_EVG, file="DAYMET.RDATA")

citation("daymetr")
