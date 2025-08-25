# 1.9 PDSI Layers

rm(list=ls())

library(ncdf4)
library(terra)

setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")

##########################################################################################################################################################
# 1. IMPORTING PDSI DATA. *** run once***
##########################################################################################################################################################
# Data downloaded from: 
# http://thredds.northwestknowledge.net:8080/thredds/ncss/grid/agg_met_pdsi_1979_CurrentYear_CONUS.nc/dataset.html
# north: 26.332183, west: -81.955031,  south: 25.088531,  east: -80.214021
# 2000/01/01  -  2020/12/31

# Citation: Abatzoglou, J.T., 2013, Development of gridded surface meteorological data for ecological applications and modeling, International Journal of Climatology, DOI: 10.1002/joc.3413
# CRS: +proj=longlat +datum=WGS84 +no_defs

# Read nc file metadata
nc <- nc_open("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Climate/agg_met_pdsi_1979_CurrentYear_CONUS.nc") 

# load files
PDSI.files <- list_files_with_exts("./Climate", "nc")
# stack it
PDSI_stack <- rast(PDSI.files)
names(PDSI_stack)
# check it out (first and last layers)
sub1 <- subset(PDSI_stack, 1)
sub2 <- subset(PDSI_stack, 1533)
plot(sub1)
plot(sub2)

# FORMATTING DATES....................................................................................................................................................
# Assign dates to layer names 
# extract day values from layer names, you will need to convert them to dates using the fact that these integers are 'days since 1900-01-01' in the Gregorian calendar. 
Dates <- as.data.frame(names(PDSI_stack))
Dates$`names(PDSI_stack)`<- gsub("[A-z]","",Dates$`names(PDSI_stack)`) #remove everything other than date
Dates$`names(PDSI_stack)`<- gsub("=","",Dates$`names(PDSI_stack)`)
Dates$`names(PDSI_stack)`<- as.numeric(Dates$`names(PDSI_stack)`) #make numeric
Dates$`names(PDSI_stack)` <- as.Date(Dates$`names(PDSI_stack)`, "1900-01-01")
# re-assign dates to layer names
names(PDSI_stack) <- Dates$`names(PDSI_stack)`

# save 
setwd("./Climate")
save(PDSI_stack, file="PDSI_stack.RDATA")
writeRaster(PDSI_stack, "PDSI_stack.tif", overwrite=TRUE)
