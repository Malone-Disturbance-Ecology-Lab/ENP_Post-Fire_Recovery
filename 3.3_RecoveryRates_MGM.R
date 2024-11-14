# RECOVERY RATE
# M. Grace McLeod
# 2023


# This script calculates recovery rates for each recovery sample point

library(dplyr)
library(ggplot2)

rm(list=ls())

setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery") 


##########################################################################################################################################################
# NDVI RECOVERY RATES
##########################################################################################################################################################

# Load all recovery point observations
load("./Recov_Combo_final.RDATA")
# load recovery points used in analysis
load("/Volumes/inwedata/Malone Lab/ENP Fire/Grace_McLeod/DRIVERS/Recov_Time_NDVI_analysis.RDATA")


# DATA PREP................................................................................

# remove unnecessary columns
head(Recov_Combo)
names(Recov_Combo)
Recov_obs <- Recov_Combo %>% dplyr::select(ptID, EndDate, Obs_Year, Obs_month, Obs_Date,
                                        B1, B2, B3, B4, B5, B7,
                                        NDVI, NBR,
                                        NDVI_rf, NBR_rf)

# go ahead and set a severity minimum (moved to Analysis script)
#summary(Recov_Time_NDVI$Severity)
#ggplot(Recov_Time_NDVI, aes(x=EcoType,y=Severity))+geom_violin()
# let's go with a minimum of 0.1 
Recov_Time_NDVI <- Recov_Time_NDVI[which(Recov_Time_NDVI$Severity >.01),] # - ~8,000 pts

# remove points not used in analysis
Recov_obs <- subset(Recov_obs, ptID %in% Recov_Time_NDVI$ptID)

#pull out NDVI at time of recovery for each point, and recovery date
names(Recov_Time_NDVI)
rec_data <- Recov_Time_NDVI %>% dplyr::select(ptID, Obs_Date, TotalFires, NDVI)
colnames(rec_data)[which(names(rec_data) == "Obs_Date")] <- "Rec_Date"
colnames(rec_data)[which(names(rec_data) == "NDVI")] <- "Rec_NDVI"
names(rec_data)
rm(Recov_Time_NDVI)


# merge with rec_data
# so each observation has the recovery date and NDVI info
Recov_obs <- merge(Recov_obs, rec_data, by="ptID")
head(Recov_obs)

# REMOVE OBSERVATIONS AFTER RECOVERY DATE
# new blank dataframe
Rec_rates <- data.frame(matrix(ncol = 18, nrow = 0))
names(Rec_rates) <- names(Recov_obs)
# loop
for ( i in unique(Recov_obs$ptID)){
  print(i)
  # subset to ptID
  subID <- Recov_obs[which(Recov_obs$ptID == i),]
  # filter out observations after recovery date
  subRecPrd <- subID[which(subID$Obs_Date <= subID$Rec_Date),]
  # append to Rec_rates df
  Rec_rates <- rbind(Rec_rates, subRecPrd)
}

length(unique(Rec_rates$ptID)) # check number of points
head(Rec_rates)

# CALCULATE RECOVERY RATES...............................................................................
# proportion of recover at each obs date (observed NDVI / target NDVI)
Rec_rates$Rec_rate <- Rec_rates$NDVI / Rec_rates$NDVI_rf
summary(Rec_rates$Rec_rate)

# can you pivot NDVI and Obs_Date columns so they are numbered/labled 1-x or a-x?

# Save 
setwd("/Volumes/inwedata/Malone Lab/ENP Fire/Grace_McLeod/RECOVERY") 
save(Rec_rates, file="Rec_rates.RDATA")

##########################################################################################################################################################
# PDSI RECOVERY RATES
##########################################################################################################################################################
setwd("/Volumes/inwedata/Malone Lab/ENP Fire/Grace_McLeod/Climate") 

# Load PDSI data
load("/Volumes/inwedata/Malone Lab/ENP Fire/Grace_McLeod/Climate/Recov_PDSI.RDATA")
head(Recov_pdsi)
load("/Volumes/inwedata/Malone Lab/ENP Fire/Grace_McLeod/RECOVERY/Rec_rates.RDATA")

# remove points not used in analysis
Recov_pdsi <- subset(Recov_pdsi, ptID %in% Rec_rates$ptID)

# convert date to year/month in both dataframes
# pdsi
head(Recov_pdsi)
Recov_pdsi$YrMo <- Recov_pdsi$Week.pdsi
class(Recov_pdsi$YrMo)
Recov_pdsi$YrMo <- format(Recov_pdsi$YrMo, "%Y-%m")
# ndvi
Rec_rates$YrMo <- format(as.Date(Rec_rates$Obs_Date, "%Y-%m-%d"), "%Y-%m")
head(Rec_rates)

# select only relevant columns
Recov_pdsi <- Recov_pdsi %>% dplyr::select(ptID, pdsi, YrMo)


# need to calculate monthly mean/min/max for each point.
# calculate mean/min/max pdsi for the month of the NDVI obs
# group by ptID and YearMonth
pdsi_rates <- Recov_pdsi %>%
  group_by(ptID, YrMo) %>%
  summarize(mo.mean=mean(pdsi, na.rm=T), mo.max=max(pdsi, na.rm = T), mo.min=min(pdsi, na.rm = T))
# check it
head(pdsi_rates)

# merge by year/month and ptID
Rec_rates <- merge(Rec_rates, pdsi_rates, by=c("ptID", "YrMo"))
head(Rec_rates)

# Save 
setwd("/Volumes/inwedata/Malone Lab/ENP Fire/Grace_McLeod/RECOVERY") 
save(Rec_rates, file="Rec_rates.RDATA")


##########################################################################################################################################################
# FILTER TO STAY RECOVERED
##########################################################################################################################################################
library(dplyr)
head(Rec_rates)

# observations available
load("/Volumes/inwedata/Malone Lab/ENP Fire/Grace_McLeod/RECOVERY/Recov_Master.RDATA")
load("/Users/gracemcleod/Documents/Thesis/Recov_Master.RDATA")
test1 <- Recov_Master %>% count(ptID)
colnames(test1)[which(names(test1) == "n")] <- "obs.total"

# observations since fire date
obs.PostFire <- Recov_Master[which(Recov_Master$Obs_Date > Recov_Master$EndDate),] 
test2 <- obs.PostFire  %>% count(ptID)
colnames(test2)[which(names(test2) == "n")] <- "obs.PostFire"

# observations since fire above threshold
load("/Volumes/inwedata/Malone Lab/ENP Fire/Grace_McLeod/RECOVERY/Recov_Combo_final.RDATA")
load("/Users/gracemcleod/Documents/Thesis/Recov_Combo_final.RDATA")
obs.Post.Above <- Recov_Combo[which(Recov_Combo$NDVI >=Recov_Combo$lwrNDVI),]
test3 <- obs.PostFire  %>% count(ptID)
colnames(test3)[which(names(test3) == "n")] <- "obs.Post.Above"

# Observations since fire below threshold
obs.Post.Below <- Recov_Combo[which(Recov_Combo$NDVI < Recov_Combo$lwrNDVI),]
test4 <- obs.PostFire  %>% count(ptID)
colnames(test4)[which(names(test4) == "n")] <- "obs.Post.Below"

# Minimum date of observation above threshold
test5 <- obs.Post.Above %>%
  group_by(ptID) %>%
  summarise(minDate=min(Obs_Date))

# Maximum date of observation since fire below threshold.
test6 <- obs.Post.Below %>%
  group_by(ptID) %>%
  summarise(minDate=min(Obs_Date))

# COMBINE 
Combo <- merge(test1, test2, by="ptID")
Combo <- merge(Combo, test3, by="ptID")
Combo <- merge(Combo, test4, by="ptID")
Combo <- merge(Combo, test5, by="ptID")
Combo <- merge(Combo, test6, by="ptID")
summary(Combo)















#######################################################
# TRASH
#######################################################

# subset to run on laptop
291743697/4 # = 72935924

pdsi_1 <- Recov_pdsi[1:72935924,]
save(pdsi_1, file="pdsi_1.RDATA")
rm(pdsi_1)
pdsi_2 <- Recov_pdsi[72935925:145871849,]
save(pdsi_2, file="pdsi_2.RDATA")
rm(pdsi_2)
pdsi_3 <- Recov_pdsi[145871850:218807774,]
save(pdsi_3, file="pdsi_3.RDATA")
rm(pdsi_3)
pdsi_4 <- Recov_pdsi[218807775:291743697,]
save(pdsi_4, file="pdsi_4.RDATA")
rm(pdsi_4)


# remove duplicates
Recov_pdsi2 <- unique(Recov_pdsi[c("ptID", "pdsi", "YrMo")])
head(Recov_pdsi)

# merge with ndvi recovery rates
test2 <- left_join(Rec_rates, Recov_pdsi, by=c("ptID", "YrMo"))
test <- merge(Rec_rates, Recov_pdsi, by=c("ptID", "YrMo"))
head(test2)

# pull out just ptID and mean/min/max
pdsi_calcs <- Recov_pdsi %>% dplyr::select(ptID, mo_mean, mo_min, mo_max)





