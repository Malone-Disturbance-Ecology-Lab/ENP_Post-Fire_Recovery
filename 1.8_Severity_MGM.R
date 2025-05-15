# FIRE SEVERITY 
# M. Grace McLeod (2023)


# This script calculates fire severity as the change in NBR pre- and post-fire for the Recovery Sample Points.
# 1. Identifies last available pre-fire image and first available post-fire image
# 2. Calculates the time between imaging dates and fire start and end dates

library(tidyverse)
library(tidyr)
library(ggplot2)
library(stats)
library(sp)
library(sf)
library(rgdal)
library(dplyr)
library(raster)
library(dichromat)

rm(list=ls())
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Severity")

##########################################################################################################################################################
# 1. CALCUALTE SEVERITY
##########################################################################################################################################################

# Load Recovery Master Dataframe
load("/Volumes/inwedata/Malone Lab/ENP Fire/Grace_McLeod/Recovery/Recov_Master.RDATA")

# Make a Severity dataframe
Sev_df <- Recov_Master %>% dplyr::select(ptID, StartDate, EndDate, FireYear) %>%
  distinct()
# add columns for pre and post values
Sev_df$PreNBR <- NA
Sev_df$PreNDVI <- NA
Sev_df$PostNBR <- NA


# IMMEDIATE PRE AND POST FIRE SPEC ................................................................................................................................................

# start off cycling through points. 
# loop thorugh points in sev_df
for (i in Sev_df$ptID){
    print(i)
  # subset obs for this point from master
  subID <- Recov_Master[which(Recov_Master$ptID == i),]
  # Immediate pre-fire
      # pull out pre-fire obs 
      StartDate <- as.Date(Sev_df$StartDate[Sev_df$ptID ==i])
      try(subPRE <- subID[which(subID$Obs_Date < StartDate),])
      rm(StartDate)
      # take max and write into sev_df 
      try(Sev_df$PreDate[Sev_df$ptID == i] <- as.character(subPRE$Obs_Date[subPRE$Obs_Date == max(subPRE$Obs_Date)]))
      try(Sev_df$PreNBR[Sev_df$ptID == i] <- subPRE$NBR[subPRE$Obs_Date == max(subPRE$Obs_Date)])
      try(Sev_df$PreNDVI[Sev_df$ptID == i] <- subPRE$NDVI[subPRE$Obs_Date == max(subPRE$Obs_Date)])
  # Immediate post-fire
      # pull out post-fire obs 
      EndDate <- as.Date(Sev_df$EndDate[Sev_df$ptID ==i])
      try(subPOST <- subID[which(subID$Obs_Date > EndDate),])
      rm(EndDate)
      # take min and write into sev_df 
      try(Sev_df$PostDate[Sev_df$ptID == i] <- as.character(subPOST$Obs_Date[subPOST$Obs_Date == min(subPOST$Obs_Date)]))
      try(Sev_df$PostNBR[Sev_df$ptID == i] <- subPOST$NBR[subPOST$Obs_Date == min(subPOST$Obs_Date)])

}

# check it: Are there any cases where pre date > post date? (data entry errors from the Park Service)
class(Sev_df$PreDate)  ;  class(Sev_df$StartDate)
Sev_df$PreDate <- as.Date(Sev_df$PreDate)  ;  Sev_df$StartDate <- as.Date(Sev_df$StartDate)
test <- Sev_df[which(Sev_df$StartDate > Sev_df$EndDate),]
length(unique(test$ptID))
# If so, remove from dataframe
# Sev_df <- anti_join(Sev_df, test)


# SEVERITY CALCULATION .................................................................................................................................................
Sev_df$DateDif <- as.Date(Sev_df$PostDate) - as.Date(Sev_df$PreDate)
Sev_df$Severity <- Sev_df$PreNBR -Sev_df$PostNBR

# check it out 
summary(Sev_df$Severity)
summary(Sev_df$PreNBR) 
summary(Sev_df$PostNBR) # NAs because EndDate close to end of sample period and no clear image after OR EndDate not recorded

# round to just two decimal places
Sev_df$Severity <- format(round(Sev_df$Severity, 2), nsmall = 2)
Sev_df$Severity <- as.numeric(Sev_df$Severity)

# plot
ggplot(Sev_df, aes(y=Severity)) +
  geom_boxplot() +
  scale_color_manual(values=c( "#ff6633"))+ 
  labs( y="Fire Severity" ) +
  theme_bw()+
  theme(legend.position = "none") +
  theme(axis.title = element_text(size = 20), axis.text = element_text(size=20))
Sev.1 <- Sev_df[which(Sev_df$Severity >= .1),]
pre <- ggplot(Sev.1, aes(y=PreNBR)) +
  geom_boxplot() +
  scale_color_manual(values=c( "#ff6633"))+ 
  labs( y="Pre-NBR" ) +
  theme_bw()+
  theme(legend.position = "none") +
  theme(axis.title = element_text(size = 20), axis.text = element_text(size=20))
post <- ggplot(Sev.1, aes(y=PostNBR)) +
  geom_boxplot() +
  scale_color_manual(values=c( "#ff6633"))+ 
  labs( y="Post-NBR" ) +
  theme_bw()+
  theme(legend.position = "none") +
  theme(axis.title = element_text(size = 20), axis.text = element_text(size=20))
pre + post

# save 
#setwd("/Volumes/inwedata/Malone Lab/ENP Fire/Grace_McLeod/Severity")
write.csv(Sev_df, file="Severity.csv")
save(Sev_df, file="Sev_df.RDATA")


##########################################################################################################################################################
#  2. PRE- POST- DATE DIFFERENCE
##########################################################################################################################################################
load("./Sev_df.RDATA")

# CALCUALTE RANGE BETWEEN START-DATE AND PRE-DATE...................................................................................................
class(Sev_df$PreDate)  ;  class(Sev_df$StartDate)
Sev_df$PreDate <- as.Date(Sev_df$PreDate)  ;  Sev_df$StartDate <- as.Date(Sev_df$StartDate)
Sev_df$PreDateDif <- as.numeric(Sev_df$StartDate - Sev_df$PreDate)
summary(Sev_df$PreDateDif)
#   Min.    1st Qu.    Median    Mean   3rd Qu.    Max.    
#  1.00     9.00      25.00     37.33   52.00    1028.00    
# Almost all points have a pre-fire obs within 100 days of fire, and top 3rd quartile have obs within 38 days. 
length(unique(Sev_df$ptID[Sev_df$PreDateDif <= 38])) #  99,577
length(unique(Sev_df$ptID[Sev_df$PreDateDif > 38])) # 57,746
# Reasonable to say you need pre-fire value within one month or two satalite passes! 

# PLOT to visualize distribution
ggplot(aes(x=PreDateDif, y=PreNBR), data=Sev_df) +
  geom_point(color="#ff6633") +
  geom_smooth(color="black") +
  labs(x= "Days Pre-Fire", y="Pre-Fire NBR") +
  ggtitle("Pineland")  + 
  theme(text = element_text(size = 20)) +
  ylim(-0.4, 0.4) +
  scale_x_reverse()+
  theme(panel.background = element_blank()) +
  theme(panel.grid.major = element_line(color = "black",
                                        linewidth = 0.05,
                                        linetype = 1)) 


# CALCUALTE RANGE BETWEEN END-DATE AND POST-DATE...................................................................................................
class(Sev_df$PostDate)  ;  class(Sev_df$EndDate)
Sev_df$PostDate <- as.Date(Sev_df$PostDate)  ;  Sev_df$EndDate <- as.Date(Sev_df$EndDate)
Sev_df$PostDateDif <- as.numeric(Sev_df$PostDate - Sev_df$EndDate)
summary(Sev_df$PostDateDif)
#    Min.   1st Qu.    Median    Mean   3rd Qu.    Max.    NA's 
#    1.00   13.00     31.00     46.12   66.00    602.00    4266 
# Lots of NAs due to fires with no recorded EndDate. 
no_end <- Sev_df[which(is.na(Sev_df$EndDate)),]
# Sev_df <- Sev_df %>% drop_na(PostDateDif)
# Almost all points have a post-fire obs within 100 days of fire. 
length(unique(Sev_df$ptID[Sev_df$PostDateDif <= 41])) # 112,942
length(unique(Sev_df$ptID[Sev_df$PostDateDif > 41])) # 36,886
# Reasonable to say you need post-fire value within one month or two satalite passes! 

# SAVE updates 
#setwd("/Volumes/inwedata/Malone Lab/ENP Fire/Grace_McLeod/Severity")
write.csv(Sev_df, file="Severity.csv")
save(Sev_df, file="Sev_df.RDATA")






