# PALMER DROUGHT SEVERITY INDEX (PDSI) DATA
  # M. Grace McLeod (2023)


# This script calculates mean pdsi values for recovery points at varing 
# intervals pre- and post-fire, using PDSI data on a 5 day temporal resolution.
# 1. Processes PDSI data, converting .nc files to a rastser stack
# 2. Extracts PDSI data to recovery sample points
# 3. Summarizes pdsi.....
  # during the recovery period
  # 1 year pre-post 
  # 6 months pre-post
  # 1 month pre-post
  # 1 week pre-post


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
# IMPORTING PDSI DATA. *** run once***
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
# EXTRACT TO RECOVERY SAMPLE POINTS 
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
load("./Recovery/Recov_Rates.RDATA")
# subset to just recovered points
Recov_pdsi <- subset(Recov_pdsi, ptID %in% Recov_Rate$ptID)
length(unique(Recov_pdsi$ptID))

# Check it out
summary(Recov_pdsi)
#Min.       1st Qu.  Median    Mean    3rd Qu.    Max. 
# -5.3000 -1.9000   -1.0000  -0.3175  1.4000    8.1000

# save 
setwd("./Climate")
save(Recov_pdsi, file="Recov_PDSI.RDATA") 


##########################################################################################################################################################
# RECOVERY PERIOD SUMMARY STATISTICS
##########################################################################################################################################################
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")

# Load pdsi data
load("./Climate/Recov_PDSI.RDATA")
# Load Recovery Rates DF
load("./Recovery/Recov_Rates.RDATA")
# make reduced df with just ptID and time required to get to each threshold
rates_dates <- Recov_Rate %>%
  select(ptID, StartDate,
         Rec10_Date, 
         Rec20_Date,
         Rec30_Date,
         Rec40_Date,  
         Rec50_Date,
         Rec60_Date, 
         Rec70_Date,
         Rec80_Date,
         Rec90_Date,
         Rec100_Date)

# calculate pdsi monthly insead of weekly
Recov_pdsi$YrMo.pdsi <- format(Recov_pdsi$Week.pdsi, "%Y-%m")


# Filter pdsi obs to just recovery period
for (i in 1:length(rates_dates$ptID)) {
  print(i)
  # pull out fire StartDate 
  StartDate <- as.Date(rates_dates$StartDate[i])
  # subset pdsi data by threshold
  pdsi10 <- Recov_pdsi %>%
    filter(Week.pdsi > StartDate &  Week.pdsi < rates_dates$Rec10_Date[i])
  rates_dates$pdsi10.mean[i] <- mean(pdsi10$pdsi, na.rm=T)
  rates_dates$pdsi10.var[i] <- var(pdsi10$pdsi, na.rm=T)
  rates_dates$pdsi10.min[i] <- min(pdsi10$pdsi, na.rm=T)
  rates_dates$pdsi10.max[i] <- max(pdsi10$pdsi, na.rm=T)
  nwet10 <- pdsi10 %>% filter(pdsi > 4 ) %>% nrow() 
  ndry10 <- pdsi10 %>% filter(pdsi < -4 ) %>% nrow()
  ntot10 <- pdsi10 %>% nrow()
  rates_dates$pdsi10.wet[i] <- nwet10 / ntot10 *100
  rates_dates$pdsi10.dry[i] <- ndry10 / ntot10 *100
  rm(pdsi10, nwet10, ndry10, ntot10)
  
  pdsi20 <- Recov_pdsi %>%
    filter(Week.pdsi > StartDate &  Week.pdsi < rates_dates$Rec20_Date[i])
  rates_dates$pdsi20.mean[i] <- mean(pdsi20$pdsi, na.rm=T)
  rates_dates$pdsi20.var[i] <- var(pdsi20$pdsi, na.rm=T)
  rates_dates$pdsi20.min[i] <- min(pdsi20$pdsi, na.rm=T)
  rates_dates$pdsi20.max[i] <- max(pdsi20$pdsi, na.rm=T)
  nwet20 <- pdsi20 %>% filter(pdsi > 4 ) %>% nrow() 
  ndry20 <- pdsi20 %>% filter(pdsi < -4 ) %>% nrow()
  ntot20 <- pdsi20 %>% nrow()
  rates_dates$pdsi20.wet[i] <- nwet20 / ntot20 *100
  rates_dates$pdsi20.dry[i] <- ndry20 / ntot20 *100
  rm(pdsi20, nwet20, ndry20, ntot20)
  
  pdsi30 <- Recov_pdsi %>%
    filter(Week.pdsi > StartDate &  Week.pdsi < rates_dates$Rec30_Date[i])
  rates_dates$pdsi30.mean[i] <- mean(pdsi30$pdsi, na.rm=T)
  rates_dates$pdsi30.var[i] <- var(pdsi30$pdsi, na.rm=T)
  rates_dates$pdsi30.min[i] <- min(pdsi30$pdsi, na.rm=T)
  rates_dates$pdsi30.max[i] <- max(pdsi30$pdsi, na.rm=T)
  nwet30 <- pdsi30 %>% filter(pdsi > 4 ) %>% nrow() 
  ndry30 <- pdsi30 %>% filter(pdsi < -4 ) %>% nrow()
  ntot30 <- pdsi30 %>% nrow()
  rates_dates$pdsi30.wet[i] <- nwet30 / ntot30 *100
  rates_dates$pdsi30.dry[i] <- ndry30 / ntot30 *100
  rm(pdsi30, nwet30, ndry30, ntot30)
  
  pdsi40 <- Recov_pdsi %>%
    filter(Week.pdsi > StartDate &  Week.pdsi < rates_dates$Rec40_Date[i])
  rates_dates$pdsi40.mean[i] <- mean(pdsi40$pdsi, na.rm=T)
  rates_dates$pdsi40.var[i] <- var(pdsi40$pdsi, na.rm=T)
  rates_dates$pdsi40.min[i] <- min(pdsi40$pdsi, na.rm=T)
  rates_dates$pdsi40.max[i] <- max(pdsi40$pdsi, na.rm=T)
  nwet40 <- pdsi40 %>% filter(pdsi > 4 ) %>% nrow() 
  ndry40 <- pdsi40 %>% filter(pdsi < -4 ) %>% nrow()
  ntot40 <- pdsi40 %>% nrow()
  rates_dates$pdsi40.wet[i] <- nwet40 / ntot40 *100
  rates_dates$pdsi40.dry[i] <- ndry40 / ntot40 *100
  rm(pdsi40, nwet40, ndry40, ntot40)
  
  pdsi50 <- Recov_pdsi %>%
    filter(Week.pdsi > StartDate &  Week.pdsi < rates_dates$Rec50_Date[i])
  rates_dates$pdsi50.mean[i] <- mean(pdsi50$pdsi, na.rm=T)
  rates_dates$pdsi50.var[i] <- var(pdsi50$pdsi, na.rm=T)
  rates_dates$pdsi50.min[i] <- min(pdsi50$pdsi, na.rm=T)
  rates_dates$pdsi50.max[i] <- max(pdsi50$pdsi, na.rm=T)
  nwet50 <- pdsi50 %>% filter(pdsi > 4 ) %>% nrow() 
  ndry50 <- pdsi50 %>% filter(pdsi < -4 ) %>% nrow()
  ntot50 <- pdsi50 %>% nrow()
  rates_dates$pdsi50.wet[i] <- nwet50 / ntot50 *100
  rates_dates$pdsi50.dry[i] <- ndry50 / ntot50 *100
  rm(pdsi50, nwet50, ndry50, ntot50)
  
  pdsi60 <- Recov_pdsi %>%
    filter(Week.pdsi > StartDate &  Week.pdsi < rates_dates$Rec60_Date[i])
  rates_dates$pdsi60.mean[i] <- mean(pdsi60$pdsi, na.rm=T)
  rates_dates$pdsi60.var[i] <- var(pdsi60$pdsi, na.rm=T)
  rates_dates$pdsi60.min[i] <- min(pdsi60$pdsi, na.rm=T)
  rates_dates$pdsi60.max[i] <- max(pdsi60$pdsi, na.rm=T)
  nwet60 <- pdsi60 %>% filter(pdsi > 4 ) %>% nrow() 
  ndry60 <- pdsi60 %>% filter(pdsi < -4 ) %>% nrow()
  ntot60 <- pdsi60 %>% nrow()
  rates_dates$pdsi60.wet[i] <- nwet60 / ntot60 *100
  rates_dates$pdsi60.dry[i] <- ndry60 / ntot60 *100
  rm(pdsi60, nwet60, ndry60, ntot60)
  
  pdsi70 <- Recov_pdsi %>%
    filter(Week.pdsi > StartDate &  Week.pdsi < rates_dates$Rec70_Date[i])
  rates_dates$pdsi70.mean[i] <- mean(pdsi70$pdsi, na.rm=T)
  rates_dates$pdsi70.var[i] <- var(pdsi70$pdsi, na.rm=T)
  rates_dates$pdsi70.min[i] <- min(pdsi70$pdsi, na.rm=T)
  rates_dates$pdsi70.max[i] <- max(pdsi70$pdsi, na.rm=T)
  nwet70 <- pdsi70 %>% filter(pdsi > 4 ) %>% nrow() 
  ndry70 <- pdsi70 %>% filter(pdsi < -4 ) %>% nrow()
  ntot70 <- pdsi70 %>% nrow()
  rates_dates$pdsi70.wet[i] <- nwet70 / ntot70 *100
  rates_dates$pdsi70.dry[i] <- ndry70 / ntot70 *100
  rm(pdsi70, nwet70, ndry70, ntot70)
  
  pdsi80 <- Recov_pdsi %>%
    filter(Week.pdsi > StartDate &  Week.pdsi < rates_dates$Rec80_Date[i])
  rates_dates$pdsi80.mean[i] <- mean(pdsi80$pdsi, na.rm=T)
  rates_dates$pdsi80.var[i] <- var(pdsi80$pdsi, na.rm=T)
  rates_dates$pdsi80.min[i] <- min(pdsi80$pdsi, na.rm=T)
  rates_dates$pdsi80.max[i] <- max(pdsi80$pdsi, na.rm=T)
  nwet80 <- pdsi80 %>% filter(pdsi > 4 ) %>% nrow() 
  ndry80 <- pdsi80 %>% filter(pdsi < -4 ) %>% nrow()
  ntot80 <- pdsi80 %>% nrow()
  rates_dates$pdsi80.wet[i] <- nwet80 / ntot80 *100
  rates_dates$pdsi80.dry[i] <- ndry80 / ntot80 *100
  rm(pdsi80, nwet80, ndry80, ntot80)
  
  pdsi90 <- Recov_pdsi %>%
    filter(Week.pdsi > StartDate &  Week.pdsi < rates_dates$Rec90_Date[i])
  rates_dates$pdsi90.mean[i] <- mean(pdsi90$pdsi, na.rm=T)
  rates_dates$pdsi90.var[i] <- var(pdsi90$pdsi, na.rm=T)
  rates_dates$pdsi90.min[i] <- min(pdsi90$pdsi, na.rm=T)
  rates_dates$pdsi90.max[i] <- max(pdsi90$pdsi, na.rm=T)
  nwet90 <- pdsi90 %>% filter(pdsi > 4 ) %>% nrow() 
  ndry90 <- pdsi90 %>% filter(pdsi < -4 ) %>% nrow()
  ntot90 <- pdsi90 %>% nrow()
  rates_dates$pdsi90.wet[i] <- nwet90 / ntot90 *100
  rates_dates$pdsi90.dry[i] <- ndry90 / ntot90 *100
  rm(pdsi90, nwet90, ndry90, ntot90)
  
  pdsi100 <- Recov_pdsi %>%
    filter(Week.pdsi > StartDate &  Week.pdsi < rates_dates$Rec100_Date[i])
  rates_dates$pdsi100.mean[i] <- mean(pdsi100$pdsi, na.rm=T)
  rates_dates$pdsi100.var[i] <- var(pdsi100$pdsi, na.rm=T)
  rates_dates$pdsi100.min[i] <- min(pdsi100$pdsi, na.rm=T)
  rates_dates$pdsi100.max[i] <- max(pdsi100$pdsi, na.rm=T)
  nwet100 <- pdsi100 %>% filter(pdsi > 4 ) %>% nrow() 
  ndry100 <- pdsi100 %>% filter(pdsi < -4 ) %>% nrow()
  ntot100 <- pdsi100 %>% nrow()
  rates_dates$pdsi100.wet[i] <- nwet100 / ntot100 *100
  rates_dates$pdsi100.dry[i] <- ndry100 / ntot100 *100
  rm(pdsi100, nwet100, ndry100, ntot100)

}
setwd("./Climate")
save(rates_dates, file="Rec_pdsi_summary.RDATA")


##########################################################################################################################################################
# POST-FIRE CLIMATE
##########################################################################################################################################################

# add average pdsi columns
Recov_Time_NDVI$RecPrd.pdsi <- NA
Recov_Time_NDVI$Post.1yr.pdsi <- NA
Recov_Time_NDVI$Post.6mo.pdsi <- NA
Recov_Time_NDVI$Post.3mo.pdsi <- NA
Recov_Time_NDVI$Post.1mo.pdsi <- NA
Recov_Time_NDVI$Post.1wk.pdsi <- NA


# LOOP
# for each ptID
for (i in unique(pdsi_master$ptID)) {
  print(i)
  # subset to ptID
  subID <- pdsi_master[which(pdsi_master$ptID == i),]
  # pull out fire EndDate and Recovery Date
  EndDate <- Recov_Time_NDVI$EndDate[Recov_Time_NDVI$ptID == i]
  RecDate <- Recov_Time_NDVI$Obs_Date[Recov_Time_NDVI$ptID == i]
  # subset pdsi observations to recovery period (EndDate -> Rec_Date)
  subDate <- subID[which(subID$Week.pdsi >= EndDate & subID$Week.pdsi <= RecDate),]
  # take mean of all pdsi values and write into Recov_Time df
  Recov_Time_NDVI$RecPrd.pdsi[Recov_Time_NDVI$ptID == i] <- mean(subDate$pdsi)
  
  # try taking mean for 1yr post (assuming points took at least 1 year to recover)
  subDate <- subDate[order(subDate$Week.pdsi),] # make sure observations are ordered chronologically
  try(Recov_Time_NDVI$Post.1yr.pdsi[Recov_Time_NDVI$ptID == i] <- mean(subDate$pdsi[1:73])) # 73 five-day obs per year
  # try 6 mo post
  try(Recov_Time_NDVI$Post.6mo.pdsi[Recov_Time_NDVI$ptID == i] <- mean(subDate$pdsi[1:37])) # ~ 37 obs
  # try 3 mo post
  try(Recov_Time_NDVI$Post.3mo.pdsi[Recov_Time_NDVI$ptID == i] <- mean(subDate$pdsi[1:18])) # ~ 18 obs
  # try 1 mo post
  try(Recov_Time_NDVI$Post.1mo.pdsi[Recov_Time_NDVI$ptID == i] <- mean(subDate$pdsi[1:6])) # ~ 6 obs
  # try 1 week post
  try(Recov_Time_NDVI$Post.1wk.pdsi[Recov_Time_NDVI$ptID == i] <- mean(subDate$pdsi[1:2])) # ~ 2 obs (puts it anywhere from 5-10 days)
}



# Save
#setwd("/Volumes/inwedata/Malone Lab/ENP Fire/Grace_McLeod/Climate")
PostFire_pdsi <- Recov_Time_NDVI
save(PostFire_pdsi, file="PostFire_pdsi.RDATA")


##########################################################################################################################################################
# PRE-FIRE CLIMATE
#######################################################################################################################################################

# add average pdsi columns
Recov_Time_NDVI$Pre.1yr.pdsi <- NA
Recov_Time_NDVI$Pre.6mo.pdsi <- NA
Recov_Time_NDVI$Pre.3mo.pdsi <- NA
Recov_Time_NDVI$Pre.1mo.pdsi <- NA
Recov_Time_NDVI$Pre.1wk.pdsi <- NA

# LOOP
# for each ptID
for (i in unique(pdsi_master$ptID)) {
  print(i)
  # subset to ptID
  subID <- pdsi_master[which(pdsi_master$ptID == i),]
  # pull out fire EndDate and Recovery Date
  StartDate <- Recov_Time_NDVI$StartDate[Recov_Time_NDVI$ptID == i]
  # subset pdsi observations to recovery period (EndDate -> Rec_Date)
  subDate <- subID[which(subID$Week.pdsi <= StartDate),]
  
  # try taking mean for 1yr post (assuming points took at least 1 year to recover)
  subDate <- subDate[order(subDate$Week.pdsi, decreasing = T),] # make sure observations REVERSE chronological
  try(Recov_Time_NDVI$Pre.1yr.pdsi[Recov_Time_NDVI$ptID == i] <- mean(subDate$pdsi[1:73])) # 73 five-day obs per year
  # try 6 mo pre
  try(Recov_Time_NDVI$Pre.6mo.pdsi[Recov_Time_NDVI$ptID == i] <- mean(subDate$pdsi[1:37])) # ~ 37 obs
  # try 3 mo pre
  try(Recov_Time_NDVI$Pre.3mo.pdsi[Recov_Time_NDVI$ptID == i] <- mean(subDate$pdsi[1:18])) # ~ 18 obs
  # try 1 mo pre
  try(Recov_Time_NDVI$Pre.1mo.pdsi[Recov_Time_NDVI$ptID == i] <- mean(subDate$pdsi[1:6])) # ~ 6 obs
  # try 1 week pre
  try(Recov_Time_NDVI$Pre.1wk.pdsi[Recov_Time_NDVI$ptID == i] <- mean(subDate$pdsi[1:2])) # ~ 2 obs (puts it anywhere from 5-10 days)
}



# Save
#setwd("/Volumes/inwedata/Malone Lab/ENP Fire/Grace_McLeod/Climate")
PreFire_pdsi <- Recov_Time_NDVI
save(PreFire_pdsi, file="PreFire_pdsi.RDATA")



##########################################################################################################################################################
# TROUBLE SHOOTING
#######################################################################################################################################################

# ~2,000 points had NA's for post-fire pdsi values. 
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Climate/PostFire_pdsi.RDATA")

# pull out points with NAs
post.NAs <- PostFire_pdsi[apply(is.na(PostFire_pdsi[,13:17]), 1, all),]
RecPrd.NAs <- PostFire_pdsi[is.na(PostFire_pdsi$RecPrd.pdsi),]

# 1. The points that recovered in less than 10 days did not get a value for 1wk pdsi or above. 
# 2. The points that recovered in less than 5 days did not get a RecPrd pdsi

# RUN LOOP TO ASSIGN WEEK-OF PDSI TO 1wk AND RecPrd FOR THOSE POINTS

# subset pdsi_master to those points
pdsi_master <- subset(pdsi_master, ptID %in% RecPrd.NAs$ptID)

# LOOP
# for each ptID
for (i in unique(pdsi_master$ptID)) {
  print(i)
  # subset to ptID
  subID <- pdsi_master[which(pdsi_master$ptID == i),]
  # pull out fire EndDate and Recovery Date
  EndDate <- Recov_Time_NDVI$EndDate[Recov_Time_NDVI$ptID == i]
  RecDate <- Recov_Time_NDVI$Obs_Date[Recov_Time_NDVI$ptID == i]
  # subset pdsi observations to recovery period (EndDate -> Rec_Date)
  subDate <- subID[which(subID$Week.pdsi >= EndDate & subID$Week.pdsi <= RecDate),]
  # take mean of all pdsi values and write into Recov_Time df
  Recov_Time_NDVI$RecPrd.pdsi[Recov_Time_NDVI$ptID == i] <- mean(subDate$pdsi)
  
  # try taking mean for 1yr post (assuming points took at least 1 year to recover)
  subDate <- subDate[order(subDate$Week.pdsi),] # make sure observations are ordered chronologically
  try(Recov_Time_NDVI$Post.1yr.pdsi[Recov_Time_NDVI$ptID == i] <- mean(subDate$pdsi[1:73])) # 73 five-day obs per year
  # try 6 mo post
  try(Recov_Time_NDVI$Post.6mo.pdsi[Recov_Time_NDVI$ptID == i] <- mean(subDate$pdsi[1:37])) # ~ 37 obs
  # try 3 mo post
  try(Recov_Time_NDVI$Post.3mo.pdsi[Recov_Time_NDVI$ptID == i] <- mean(subDate$pdsi[1:18])) # ~ 18 obs
  # try 1 mo post
  try(Recov_Time_NDVI$Post.1mo.pdsi[Recov_Time_NDVI$ptID == i] <- mean(subDate$pdsi[1:6])) # ~ 6 obs
  # try 1 week post
  try(Recov_Time_NDVI$Post.1wk.pdsi[Recov_Time_NDVI$ptID == i] <- mean(subDate$pdsi[1:2])) # ~ 2 obs (puts it anywhere from 5-10 days)
}





