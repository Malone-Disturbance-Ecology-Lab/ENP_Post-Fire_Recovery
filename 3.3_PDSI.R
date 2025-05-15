# PALMER DROUGHT SEVERITY INDEX (PDSI) DATA
  # M. Grace McLeod (2023)


# This script calculates mean pdsi values for recovery points at varing 
# intervals pre- and post-fire, using PDSI data on a 5 day temporal resolution.
# 1. Processes PDSI data, converting .nc files to a rastser stack
# 2. Extracts PDSI data to recovery sample points ("Recov_PDSI.RDATA")
# 3. Summarizes pdsi during the recovery period ("RecPrdPDSI_summary.RDATA")


rm(list=ls())

library(ncdf4)
library(terra)
library(tools)
library(reshape2)
library(readr)
library(ggplot2)
library(dplyr)  
library(knitr) 
library(plotly)
library(scPDSI)
library(sf)
library(sp)

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
nc <- nc_open("./Climate/agg_met_pdsi_1979_CurrentYear_CONUS.nc") 

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


##########################################################################################################################################################
# 2. EXTRACT TO RECOVERY SAMPLE POINTS 
##########################################################################################################################################################
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")

# EXTRACT PDSI DATA TO RECOVERY SAMPLE POINTS.......................................................................
# load PDSI 
PDSI_stack <- rast("./Climate/PDSI_stack.tif")
# load Recovery sample points (crs="+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")
Recov_smpl_pts <- vect("./Sampling/Recov_smpl_pts.shp")
rec.df <- as.data.frame(Recov_smpl_pts)

#  change crs of df to match stack
crs(PDSI_stack)
crs(Recov_smpl_pts) <- "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"
Recov_smpl_pts <- project(Recov_smpl_pts, PDSI_stack)
sub1 <- subset(PDSI_stack, 1)
plot(sub1)
plot(Recov_smpl_pts, add=T)

# extract PDSI to points
Recov_pdsi_sp <- terra::extract(PDSI_stack, Recov_smpl_pts, method= "bilinear", bind=T)
Recov_pdsi_df <- as.data.frame (Recov_pdsi_sp)
# round to just 1 decimal
Recov_pdsi_df <- Recov_pdsi_df%>% mutate_at(vars(2:1534), funs(round(., 1)))
# turn DF longways
Recov_pdsi <- reshape2::melt(Recov_pdsi_df, na.rm=FALSE, value.name="pdsi", id=c("ptID"))
# format columns
colnames(Recov_pdsi)[which(names(Recov_pdsi) == "variable")] <- "Week.pdsi" #rename
class(Recov_pdsi$Week.pdsi)
Recov_pdsi$Week.pdsi <- as.character(Recov_pdsi$Week.pdsi)
Recov_pdsi$Week.pdsi<- as.Date(Recov_pdsi$Week.pdsi, format= "%Y-%m-%d") #make date


# FILTER TO RECOVERED LOCATIONS ONLY........................................................................................................
# load Recov_Time
load("./Recovery/Recov_Drivers.RDATA")
# subset to just recovered points
Recov_pdsi <- subset(Recov_pdsi, ptID %in% Recov_Drivers$ptID)
length(unique(Recov_pdsi$ptID))

# Check it out
summary(Recov_pdsi)
#Min.       1st Qu.  Median    Mean    3rd Qu.    Max. 
# -5.3000 -1.9000   -1.0000  -0.3175  1.4000    8.1000

# save 
setwd("./Climate")
save(Recov_pdsi, file="Recov_PDSI.RDATA") 


##########################################################################################################################################################
# 3. RECOVERY PERIOD SUMMARY STATISTICS 
##########################################################################################################################################################
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")

# Load pdsi data
load("./Climate/Recov_PDSI.RDATA")
# Load Recovery time data
load("./Recovery/Recov_Drivers.RDATA")

# make reduced df with just ptID and time required to get to each threshold
pdsi_summary <- Recov_Drivers %>%
  select(ptID, StartDate,
         Rec_Date)

# calculate pdsi monthly instead of weekly
Recov_pdsi$YrMo.pdsi <- format(Recov_pdsi$Week.pdsi, "%Y-%m")

for (i in 1:length(pdsi_summary$ptID)) {
  print(i)
  # Extract fire start and recovery date
  StartDate <- as.Date(pdsi_summary$StartDate[i])
  RecDate <- as.Date(pdsi_summary$Obs_Date[i])
  # Filter by pdID and recovery period
  RecPrd_pdsi <- Recov_pdsi %>%
    filter(ptID == pdsi_summary$ptID[i], Week.pdsi > StartDate, Week.pdsi < RecDate)
  # Calculate summary statistics
  pdsi_summary$pdsi.mean[i] <- mean(RecPrd_pdsi$pdsi, na.rm = TRUE)
  pdsi_summary$pdsi.var[i] <- var(RecPrd_pdsi$pdsi, na.rm = TRUE)
  pdsi_summary$pdsi.min[i] <- min(RecPrd_pdsi$pdsi, na.rm = TRUE)
  pdsi_summary$pdsi.max[i] <- max(RecPrd_pdsi$pdsi, na.rm = TRUE)
  pdsi_summary$pdsi.sd[i] <- sd(RecPrd_pdsi$pdsi, na.rm = TRUE)
  # Count wet and dry weeks
  nwet <- RecPrd_pdsi %>% filter(pdsi > 4) %>% nrow() 
  ndry <- RecPrd_pdsi %>% filter(pdsi < -4) %>% nrow()
  ntot <- nrow(RecPrd_pdsi)
  # Calculate wet / dry percent of recovery period
  pdsi_summary$pdsi.wet[i] <- ifelse(ntot > 0, nwet / ntot * 100, NA)
  pdsi_summary$pdsi.dry[i] <- ifelse(ntot > 0, ndry / ntot * 100, NA)
}


setwd("./Climate")
save(pdsi_summary, file="RecPrdPDSI_summary.RDATA")




