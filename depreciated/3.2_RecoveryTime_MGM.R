# RECOVERY TIME
# M.Grace McLeod (2023)

# This script 
# 1. formats recovery point data for the baseline model 
# 2. plots baseline variable relationships with post-fire NDVI to see required conditions for saturation
# 3. runs the baseline model on recovery points
# 4. calculates recovery time 


rm(list=ls())

library(randomForest)
library(stats)
library(dplyr)
library(car)
library(devtools)
library(randomForestCI)
library(gtools)
library(ggplot2)
library(viridis)
library(splitstackshape)
library(MetBrewer)
library(cowplot)
library(lubridate)
library(patchwork)
#install_github("swager/randomForestCI")

# function to calculate performance error for random forests (Author: Sparkle Malone)
rfpred <- function( df, model){
  m.varhat.sg <- randomForestInfJack(model, df, calibrate = TRUE)
  df$model <- m.varhat.sg$y.hat
  df$var <- sqrt(m.varhat.sg$var.hat)
  
  return(df)
}


setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod") 

##########################################################################################################################################################
# FORMAT DATA
##########################################################################################################################################################

# Copy of Recov_Master 
load("./Recovery/Recov_Master.RDATA")
Recov_BL <- Recov_Master
rm(Recov_Master)
length(unique(Recov_BL$ptID)) # 157,323

# REMOVE UNBURNED POINTS
# or those with no severity/EndDate
# load Severity dataframe 
load("./Severity/Sev_df.RDATA")
# match classes
summary(Sev_df)
summary(Recov_BL)
Recov_BL$EndDate <- as.Date(Recov_BL$EndDate)
Recov_BL$StartDate <- as.Date(Recov_BL$StartDate)
# merge dataframes
Recov_BL <- merge(Recov_BL, Sev_df, by=c("ptID", "StartDate", "EndDate", "FireYear"))
head(Recov_BL)
rm(Sev_df)
# filter severity to a minimum
length(unique(Sev_df$ptID[is.na(Sev_df$Severity)]))
sev_No_NA <- Sev_df %>% drop_na(Severity) ; length(unique(sev_No_NA$ptID)) # 153,057
sev_No_NA.1 <- sev_No_NA[which(sev_No_NA$Severity > .01),] ; length(unique(sev_No_NA.1$ptID)) # 114,294 (25% reduction, so 25% had no observable change)
Recov_BL <- Recov_BL[which(Recov_BL$Severity >.01),] 
# number of remaining pts
length(unique(Recov_BL$ptID)) # 114,294


# REMOVE OBSERVATIONS PRIOR TO RECOVERY FIRE
# count them first and record number of obs per pt
test1 <- Recov_Master %>% count(ptID)
colnames(test1)[which(names(test1) == "n")] <- "n_obs_total"
Recov_BL <- merge(Recov_BL, test1, by="ptID") 
# remove observations prior to EndDate
Recov_BL <- Recov_BL[which(Recov_BL$Obs_Date > Recov_BL$EndDate),]
# count number of obs post-fire
test2 <- Recov_BL  %>% count(ptID)
colnames(test2)[which(names(test2) == "n")] <- "n_obs_post"
Recov_BL <- merge(Recov_BL, test2, by="ptID") 

# Make FireYear numeric
Recov_BL$FireYear <- as.numeric(Recov_BL$FireYear)

# Set Prev.Int with NAs to 80
Recov_BL$Prev.Int [is.na(Recov_BL$Prev.Int)] =  80
summary(Recov_BL$Prev.Int)

# save 
#setwd("./Recovery") 
#save(Recov_BL, file="Recov_BL.RDATA")


##########################################################################################################################################################
# EXPLORATION OF POST-FIRE TRENDS.  *** do not need to run ****
##########################################################################################################################################################
# baseline values need to represent spectral saturation (the point at which change is no longer detected by the sensor)
# graph variable relationships to post-fire NDVI for setting substitute baseline values


# FREQUENCY................................................................................................................
# includes code for giving figures a transparant background!

# calcualte time since fire for each observation
Recov_BL$FireYear <- as.numeric(Recov_BL$FireYear)
Recov_BL$TSF <- Recov_BL$Obs_Year - Recov_BL$FireYear
summary(Recov_BL$TSF)

# assign frequency categories
summary(Recov_BL$TotalFires)
class(Recov_BL$TotalFires)
Recov_BL$TSF <- as.factor(Recov_BL$TSF)
# plot distribution
ggplot(Recov_BL, aes(x=TotalFires)) +
  geom_histogram() +
  scale_x_continuous(breaks=1:9) +
  theme(text = element_text(size = 20))
# plot trend
Recov_BL$TotalFires <- as.factor(Recov_BL$TotalFires)

px<- ggplot(Recov_BL, aes(x=TSF, 
                          y=NDVI,
                          group=TotalFires,
                          color=TotalFires)) +
  geom_smooth() +
  labs(
    x="Years Post-Fire",
    y= "Mean NDVI",
    title = "Recovery Locations",
    color="Total Fires", 
    fill= "Fire Frequency") + 
  #scale_fill_manual(values=cols) +
  #scale_color_manual(values=cols)+
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE) +
  theme_bw()+
  theme(
    text = element_text(size = 20),
    #panel.background = element_rect(fill='transparent'), #transparent panel bg
    #plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    #panel.grid.major = element_blank(), #remove major gridlines
    #panel.grid.minor = element_blank(), #remove minor gridlines
    #legend.background = element_rect(fill='transparent'), #transparent legend bg
    #legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )
# save plot with transparent background
#ggsave('FCEprez.png', px, bg='transparent')




# PREV.INT ................................................................................................................
# not actually a variable used in baseline but nice to see

# format data
summary(Recov_BL$Prev.Int)
class(Recov_BL$Prev.Int)
Recov_BL$TSF <- as.factor(Recov_BL$TSF)
# plot distribution
ggplot(Recov_BL, aes(x=Prev.Int)) +
  geom_histogram()  +
  theme(text = element_text(size = 20))
# plot trend
Recov_BL$Prev.Int <- as.factor(Recov_BL$Prev.Int)
py <- ggplot(Recov_BL, aes(x=TSF, 
                           y=NDVI,
                           group=Prev.Int,
                           color=Prev.Int)) +
  geom_smooth() +
  labs(
    x="Years Post-Fire",
    y= "Normalized Difference Vegetation Index (NDVI)",
    title = "Change in NDVI with time since fire in Everglades Recov_BLlands",
    color="Years Since Penultimate Fire", 
    fill= "Years Since Penultimate Fire") + 
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE) +
  theme(
    text = element_text(size = 18)
  )


# Pt.B4max ................................................................................................................

# format
class(Recov_BL$Pt.B4max )
summary(Recov_BL$Pt.B4max )
# plot distribution
ggplot(Recov_BL, aes(x=Pt.B4max )) +
  geom_histogram()  +
  theme(text = element_text(size = 20))
# plot trend
Recov_BL <- Recov_BL %>% mutate(new_bin = cut(Pt.B4max, breaks=9)) # bin data
summary(Recov_BL$new_bin)
Recov_BL$new_bin <- as.factor(Recov_BL$new_bin)
ggplot(Recov_BL, aes(x=TSF, 
                     y=NDVI,
                     group=new_bin ,
                     color=new_bin )) +
  geom_smooth() +
  labs(
    x="Years Post-Fire",
    y= "Normalized Difference Vegetation Index (NDVI)",
    title = "Change in NDVI with time since fire in Everglades Recov_BLlands",
    color="Pt.B4max", 
    fill= "Pt.B4max") + 
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE) +
  theme(
    text = element_text(size = 18)
  )



# SWIR1:SWIR2 ................................................................................................................

# format
class(Recov_BL$SWIR1.SWIR2)
summary(Recov_BL$SWIR1.SWIR2)
# plot distribution
ggplot(Recov_BL, aes(x=SWIR1.SWIR2)) +
  geom_histogram()  +
  theme(text = element_text(size = 20))
# plot trend
Recov_BL <- Recov_BL %>% mutate(new_bin = cut(SWIR1.SWIR2, breaks=9)) # bin data
summary(Recov_BL$new_bin)
Recov_BL$new_bin <- as.factor(Recov_BL$new_bin)
ggplot(Recov_BL, aes(x=TSF, 
                     y=NDVI,
                     group=new_bin ,
                     color=new_bin )) +
  geom_smooth() +
  labs(
    x="Years Post-Fire",
    y= "Normalized Difference Vegetation Index (NDVI)",
    title = "Change in NDVI with time since fire in Everglades Recov_BLlands",
    color="SWIR1.SWIR2", 
    fill= "SWIR1.SWIR2") + 
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE) +
  theme(
    text = element_text(size = 18)
  )



# NIR:SWIR1 ................................................................................................................

# format
class(Recov_BL$NIR.SWIR1)
summary(Recov_BL$NIR.SWIR1)
# plot distribution
ggplot(Recov_BL, aes(x=NIR.SWIR1)) +
  geom_histogram()  +
  theme(text = element_text(size = 20))
# plot trend
Recov_BL <- Recov_BL %>% mutate(new_bin = cut(NIR.SWIR1, breaks=8)) # bin data
summary(Recov_BL$new_bin)
mean(Recov_BL$NIR.SWIR1[Recov_BL$new_bin == "(1.53,1.65]"]) # 1.561387 (use for false value)
Recov_BL$new_bin <- as.factor(Recov_BL$new_bin)
ggplot(Recov_BL, aes(x=TSF, 
                     y=NDVI,
                     group=new_bin ,
                     color=new_bin )) +
  geom_smooth() +
  labs(
    x="Years Post-Fire",
    y= "Normalized Difference Vegetation Index (NDVI)",
    title = "Change in NDVI with time since fire in Everglades Recov_BLlands",
    color="NIR.SWIR1", 
    fill= "NIR.SWIR1") + 
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE) +
  theme(
    text = element_text(size = 18)
  )

##########################################################################################################################################################
# BASELINE MODEL 
##########################################################################################################################################################

# SET VALUES ....................................................................................................................................
# for fire history and spectral variables so baseline represents "ideal" conditions closer to spectral saturation
# use baseline points 
load("./Baseline/BL_Master_filtered.RDATA")
#load("./Recovery/Recov_BL.RDATA")

# Total Fires
Recov_BL$TotalFires.sub <- NA
# assign sub TotalFires values
Recov_BL$TotalFires.sub <- 9
head(Recov_BL)
# has to replace TotalFires for model to work on it
Recov_BL$TotalFires.og <- Recov_BL$TotalFires # make a copy
Recov_BL$TotalFires <- Recov_BL$TotalFires.sub # sub in the sub values

# Previous Interval 
# use the values associated with the selected fire history values
mean(BL_Master_filtered$Prev.Int[BL_Master_filtered$TotalFires == 9]) # 6
# assign sub values
Recov_BL$Prev.Int.sub <- 6
head(Recov_BL)
Recov_BL$Prev.Int.og <- Recov_BL$Prev.Int 
Recov_BL$Prev.Int <- Recov_BL$Prev.Int.sub 

# Pt.B4max
# use the values associated with the selected fire history values
mean(BL_Master_filtered$Pt.B4max[BL_Master_filtered$TotalFires == 9]) # 21526.95
# assign substitutes
Recov_BL$Pt.B4max.sub <- 21526.95
head(Recov_BL)
Recov_BL$Pt.B4max.og <- Recov_BL$Pt.B4max
Recov_BL$Pt.B4max <- Recov_BL$Pt.B4max.sub

# Pt.B4min
# use the values associated with the selected fire history values
mean(BL_Master_filtered$Pt.B4min[BL_Master_filtered$TotalFires == 9]) # 10637.73
# assign substitutes
Recov_BL$Pt.B4min.sub <- 10637.73
head(Recov_BL)
Recov_BL$Pt.B4min.og <- Recov_BL$Pt.B4min
Recov_BL$Pt.B4min <- Recov_BL$Pt.B4min.sub

# Pt.B5max
# use the values associated with the selected fire history values
mean(BL_Master_filtered$Pt.B5max[BL_Master_filtered$TotalFires == 9]) # 15097.51
# assign substitutes
Recov_BL$Pt.B5max.sub <-15097.51
head(Recov_BL)
Recov_BL$Pt.B5max.og <- Recov_BL$Pt.B5max
Recov_BL$Pt.B5max <- Recov_BL$Pt.B5max.sub

# SWIR1.SWIR2
# use the values associated with the selected fire history values
mean(BL_Master_filtered$SWIR1.SWIR2[BL_Master_filtered$TotalFires == 9]) # 1.259455
# assign substitutes
Recov_BL$SWIR1.SWIR2.sub <- 1.259455
head(Recov_BL)
Recov_BL$SWIR1.SWIR2.og <- Recov_BL$SWIR1.SWIR2
Recov_BL$SWIR1.SWIR2 <- Recov_BL$SWIR1.SWIR2.sub

# NIR.SWIR1
# use the values associated with the selected fire history values
mean(BL_Master_filtered$NIR.SWIR1[BL_Master_filtered$TotalFires == 9]) # 1.429776
# assign substitutes
Recov_BL$NIR.SWIR1.sub <- 1.429776
head(Recov_BL)
Recov_BL$NIR.SWIR1.og <- Recov_BL$NIR.SWIR1
Recov_BL$NIR.SWIR1 <- Recov_BL$NIR.SWIR1.sub


# save
#setwd("I:/Malone Lab/ENP Fire/Grace_McLeod/Recovery") 
#save(Recov_BL, file="Recov_BL_sub.RDATA")


# APPLY BASELINE MODEL ...............................................................................................................................

# load baseline model
load("./Baseline/NDVI_rf.RDATA")

# RUN MODELS ON RECOVERY POINTS
Recov_BL$NDVI_rf <- predict(NDVI_rf, Recov_BL)  # r2=   
summary(lm(Recov_BL$NDVI ~ Recov_BL$NDVI_rf))

# save 
#setwd("I:/Malone Lab/ENP Fire/Grace_McLeod/Recovery")
#save(Recov_BL, file="Recov_BL_final.RDATA")


# CALCULATE CONFIDENCE INTERVALS

# subset dataframe (to avoid exhausting vector memory)
# create a blank dataframe that chunk dfs can be bound to
Recov_Combo <- data.frame(matrix(ncol=55, nrow=0)) 
colnames(Recov_Combo) <- colnames(Recov_BL)
# it can handle 1,000,000 at a time for sure so do first 15,000,000 then remainder
# Make list of low ends of the chunk 
a <- seq(1, 150000001, 10000)
a <- seq(14240001, 15500001, 10000)


# Loop through sub1 
for (i in a) {
  print(i)
  # subset to the chunk
  b <- i+ 9999 
  chunk <-  Recov_BL[c(i:b),] 
  # use Jacknife to get variance out of model (apply funciton)  
  # NDVI
  Recov_sub <- rfpred(df=chunk, model=NDVI_rf) 
  colnames(Recov_sub)[which(names(Recov_sub) == "model")] <- "model.NDVI"
  colnames(Recov_sub)[which(names(Recov_sub) == "var")] <- "var.NDVI"
  # NBR
  #Recov_sub <- rfpred(df=Recov_sub, model=NBR_rf)  
  #colnames(Recov_sub)[which(names(Recov_sub) == "model")] <- "model.NBR"
  #colnames(Recov_sub)[which(names(Recov_sub) == "var")] <- "var.NBR"
  # append to big dataframe
  Recov_Combo <- smartbind(Recov_Combo, Recov_sub)
  
}

# run on remaining data 
Recov_sub <- rfpred (df=Recov_BL[15510001:15560817,], model=NDVI_rf)
colnames(Recov_sub)[which(names(Recov_sub) == "model")] <- "model.NDVI"
colnames(Recov_sub)[which(names(Recov_sub) == "var")] <- "var.NDVI"
#Recov_sub <- rfpred (df=Recov_sub, model=NBR_rf)
#colnames(Recov_sub)[which(names(Recov_sub) == "model")] <- "model.NBR"
#colnames(Recov_sub)[which(names(Recov_sub) == "var")] <- "var.NBR"
Recov_Combo <- smartbind(Recov_Combo, Recov_sub)


# get lower threshold
Recov_Combo$lwrNDVI <- Recov_Combo$model.NDVI -  Recov_Combo$var.NDVI 
#Recov_Combo$lwrNBR <- Recov_Combo$model.NBR -  Recov_Combo$var.NBR
# get upper threshold
Recov_Combo$uprNDVI <- Recov_Combo$model.NDVI + Recov_Combo$var.NDVI
#Recov_Combo$uprNBR <- Recov_Combo$model.NBR + Recov_Combo$var.NBR


# save updates 
#setwd("./Recovery")
#save(Recov_Combo, file="Recov_Combo.RDATA")


##########################################################################################################################################################
# RECOVERY TIME: DATA PREP
##########################################################################################################################################################
#load("./Recovery/Recov_Combo.RDATA")

# get immediate post-fire NDVI
IMpre <- Recov_Combo[which(Recov_Combo$PostDate == Recov_Combo$Obs_Date),]
IMpre <- IMpre %>% select(ptID, NDVI)
colnames(IMpre)[which(names(IMpre) == "NDVI")] <- "PostNDVI"
Recov_Combo <- merge(Recov_Combo, IMpre, by="ptID")
rm(IMpre)
length(unique(Recov_Combo$ptID)) # 114,294

# REMOVE OBS AFTER ANTICEDENT FIRE
load("./Fire_History/FireYears_df.RDATA")
# filter FH to recov_combo pts
Recov_FH <- subset(FireYears_df, ptID %in% Recov_Combo$ptID)
# only care about 2008-2020 bc we know they only burned once before that (recFire)
names(Recov_FH)
Recov_FH <- Recov_FH %>% dplyr::select(ptID,EVG_.2008,EVG_.2009,EVG_.2010,EVG_.2011,
                                       EVG_.2012,EVG_.2013,EVG_.2014,EVG_.2015,EVG_.2016,EVG_.2017,EVG_.2018,EVG_.2019,EVG_.2020, coords.x1, coords.x2)
# Year of antecedent fire
Recov_FH$AFyear <- NA
names(Recov_FH)
Recov_FH$AFyear <- apply(Recov_FH[,2:14], 1, FUN=min,  na.rm = TRUE)
Recov_FH[sapply(Recov_FH, is.infinite)] <- NA
summary(Recov_FH$AFyear)
names(Recov_FH)
Recov_FH <- Recov_FH %>% dplyr::select(ptID, AFyear)
# merge with Recov_Combo
Recov_Combo2 <- merge(Recov_Combo, Recov_FH, by="ptID", all.x = T)
# remove observations after antecedent fire
class(Recov_Combo2$AFyear)
summary(Recov_Combo2$AFyear)
summary(Recov_Combo2$Obs_Year)
Recov_Combo3 <- Recov_Combo2[which(Recov_Combo2$Obs_Year < Recov_Combo2$AFyear | is.na(Recov_Combo2$AFyear)),]
summary(Recov_Combo3$AFyear) # should still have NAs
length(unique(Recov_Combo3$ptID)) # 114,058

# losing pts....why?
#test <- anti_join(Recov_Combo2, Recov_Combo3)
#length(unique(test$ptID)) 
#test2 <- test[which(test$Obs_Year >= test$AFyear),] # ~ 200 pts have no observations prior to the anticedent fire
#test3 <- test[which(test$Obs_Year < test$AFyear),]

# rename
Recov_Combo <- Recov_Combo3
rm(Recov_Combo2, Recov_Combo3)

# Re-assign original values to substituted columns
Recov_Combo$Pt.B4min <- Recov_Combo$Pt.B4min.og ; Recov_Combo$Pt.B4min.og <- NULL ;  Recov_Combo$Pt.B4min.sub <- NULL 
Recov_Combo$Pt.B4max <- Recov_Combo$Pt.B4max.og ; Recov_Combo$Pt.B4max.og <- NULL ;  Recov_Combo$Pt.B4max.sub <- NULL 
Recov_Combo$Pt.B5max <- Recov_Combo$Pt.B5max.og ; Recov_Combo$Pt.B5max.og <- NULL ;  Recov_Combo$Pt.B5max.sub <- NULL 
Recov_Combo$Pt.B7max <- Recov_Combo$Pt.B7max.og ; Recov_Combo$Pt.B7max.og <- NULL ;  Recov_Combo$Pt.B7max.sub <- NULL 
Recov_Combo$SWIR1.SWIR2 <- Recov_Combo$SWIR1.SWIR2.og ; Recov_Combo$SWIR1.SWIR2.og <- NULL ; Recov_Combo$SWIR1.SWIR2.sub <- NULL
Recov_Combo$NIR.SWIR1 <- Recov_Combo$NIR.SWIR1.og ; Recov_Combo$NIR.SWIR1.og <- NULL ; Recov_Combo$NIR.SWIR1.sub <- NULL
Recov_Combo$NDVI_rf <- NULL


# save
setwd("./Recovery") 
save(Recov_Combo, file="Recov_Combo2.RDATA")


##########################################################################################################################################################
# RECOVERY TIME TO MAX THRESHOLD
##########################################################################################################################################################
# Who REACHES the threshold AT ALL 

load("./Recov_Combo2.RDATA")

# 100% .........................................................................................................
# reach lower limit of model confidence interval 
# pull out all obs with NDVI > lowNDVI  
Recov_100 <- Recov_Combo[which(Recov_Combo$NDVI >= Recov_Combo$lwrNDVI),]
# summarize by pt to get min date (first date with NDVI above threshold)
Rec_Date <- Recov_100 %>% 
  group_by(ptID) %>%
  summarise(Rec_Date = min(Obs_Date, na.rm = T)) 
# join min_date with Recov_Combo by ptID
Recov_100 <- left_join(Recov_100, Rec_Date, by="ptID")
rm(Rec_Date)
# how many pts reach target?
length(unique(Recov_100$ptID)) # 60,296
# CALCULATE RECOVERY TIME
# pull out distinct pt info
Recov_Time_100 <- Recov_100 %>% 
  dplyr::select(ptID, Rec_Date, EndDate, FireYear, Severity) %>%
  distinct()
# calculate recovery time
class(Recov_Time_100$EndDate)
Recov_Time_100$Rec_Date <- as.Date(Recov_Time_100$Rec_Date)
Recov_Time_100$EndDate <- as.Date(Recov_Time_100$EndDate)
Recov_Time_100$Rec_Days <- Recov_Time_100$Rec_Date - Recov_Time_100$EndDate
Recov_Time_100$Rec_Yrs <- as.numeric(Recov_Time_100$Rec_Days) /356
head(Recov_Time_100)
# add back other info from Recov_Combo
Recov_Time_100$Obs_Date <- Recov_Time_100$Rec_Date
names(Recov_Time_100)
Recov_Time_100 <- merge(Recov_Time_100, Recov_100, by=c("ptID","Rec_Date", "EndDate" ,"FireYear" ,"Severity", "Obs_Date"))

# check it out
summary(Recov_Time_100$Rec_Yrs) 
# Min.      1st Qu.      Median     Mean      3rd Qu.      Max. 
# 0.002809  0.741573  1.240169   2.096887  2.308989   20.087079 
Recov_Time_100$Rec_Days <- as.numeric(Recov_Time_100$Rec_Days)
ggplot(Recov_Time_100, aes(y=Rec_Yrs)) +
  geom_boxplot(col="orange") + 
  geom_density(col="purple") + 
  labs(title = "Recover to Target (100%)") +
  theme_bw()
# percentage that gets this far
60296 / 114058 * 100 # = 52.86% 
# number of fires
length(unique(Recov_Combo$FireName)) # 172
length(unique(Recov_Time_100$FireName)) # 162

# 90%...............................................................................................
# remove pts that recovered to 100
Recov_90 <- anti_join(Recov_Combo, Recov_Time_100, by="ptID")
# calculate 90% of lower limit
Recov_90$lwr90 <- Recov_90$lwrNDVI * 0.9
# filter for obs at or above threshold
Recov_90 <- Recov_90[which(Recov_90$NDVI >= Recov_90$lwr90),]
length(unique(Recov_90$ptID)) #19,553
# summarize by pt to get min date (first date with NDVI above threshold)
Rec_Date <- Recov_90 %>% 
  group_by(ptID) %>%
  summarise(Rec_Date = min(Obs_Date, na.rm = T)) 
# join min_date with Recov_Combo by ptID
Recov_90 <- left_join(Recov_90, Rec_Date, by="ptID")
rm(Rec_Date)
# how many pts reach target?
length(unique(Recov_90$ptID)) # 19,553
# CALCULATE RECOVERY TIME
# pull out distinct pt info
Recov_Time_90 <- Recov_90 %>% 
  dplyr::select(ptID, Rec_Date, EndDate, FireYear, Severity) %>%
  distinct()
# calculate recovery time
class(Recov_Time_90$EndDate)
Recov_Time_90$Rec_Date <- as.Date(Recov_Time_90$Rec_Date)
Recov_Time_90$EndDate <- as.Date(Recov_Time_90$EndDate)
Recov_Time_90$Rec_Days <- Recov_Time_90$Rec_Date - Recov_Time_90$EndDate
Recov_Time_90$Rec_Yrs <- as.numeric(Recov_Time_90$Rec_Days) /356
head(Recov_Time_90)
# add back other info from Recov_Combo
Recov_Time_90$Obs_Date <- Recov_Time_90$Rec_Date
names(Recov_Time_90)
Recov_Time_90 <- merge(Recov_Time_90, Recov_90, by=c("ptID","Rec_Date", "EndDate" ,"FireYear" ,"Severity", "Obs_Date"))
rm(Recov_90)

# check it out
summary(Recov_Time_90$Rec_Yrs) 
# Min.      1st Qu.      Median     Mean      3rd Qu.      Max. 
# 0.002809  1.126404  1.615169  2.776745   3.643258    19.862360 
Recov_Time_90$Rec_Days <- as.numeric(Recov_Time_90$Rec_Days)
ggplot(Recov_Time_90, aes(y=Rec_Yrs)) +
  geom_boxplot(col="orange") + 
  geom_density(col="purple") + 
  labs(title = "Reach 90% threshold") +
  theme_bw()
# percentage that gets this far
length(unique(Recov_Combo$ptID))
19553 / 114058 * 100 # = 17.14% of all pts (36.37 % of the 53,762 "unrecovered" pts)
# number of fires
length(unique(Recov_Time_90$FireName)) # 139

# 80%...............................................................................................
# remove pts that recovered beyond threshold
Recov_80 <- anti_join(Recov_Combo, Recov_Time_100, by="ptID")
Recov_80 <- anti_join(Recov_80, Recov_Time_90, by="ptID")
# calculate 80% of lower limit
Recov_80$lwr80 <- Recov_80$lwrNDVI * 0.8
# filter for obs at or above threshold
Recov_80 <- Recov_80[which(Recov_80$NDVI >= Recov_80$lwr80),]
# summarize by pt to get min date (first date with NDVI above threshold)
Rec_Date <- Recov_80 %>% 
  group_by(ptID) %>%
  summarise(Rec_Date = min(Obs_Date, na.rm = T)) 
# join min_date with Recov_Combo by ptID
Recov_80 <- left_join(Recov_80, Rec_Date, by="ptID")
rm(Rec_Date)
# how many pts reach target?
length(unique(Recov_80$ptID)) # 16,857
# CALCULATE RECOVERY TIME
# pull out distinct pt info
Recov_Time_80 <- Recov_80 %>% 
  dplyr::select(ptID, Rec_Date, EndDate, FireYear, Severity) %>%
  distinct()
# calculate recovery time
class(Recov_Time_80$EndDate)
Recov_Time_80$Rec_Date <- as.Date(Recov_Time_80$Rec_Date)
Recov_Time_80$EndDate <- as.Date(Recov_Time_80$EndDate)
Recov_Time_80$Rec_Days <- Recov_Time_80$Rec_Date - Recov_Time_80$EndDate
Recov_Time_80$Rec_Yrs <- as.numeric(Recov_Time_80$Rec_Days) /356
head(Recov_Time_80)
# add back other info from Recov_Combo
Recov_Time_80$Obs_Date <- Recov_Time_80$Rec_Date
names(Recov_Time_80)
Recov_Time_80 <- merge(Recov_Time_80, Recov_80, by=c("ptID","Rec_Date", "EndDate" ,"FireYear" ,"Severity", "Obs_Date"))
rm(Recov_80)

# check it out
summary(Recov_Time_80$Rec_Yrs) 
# Min.      1st Qu.      Median     Mean      3rd Qu.      Max. 
# 0.002809  0.867978   1.356742  2.289042   2.890449   19.342697 
Recov_Time_80$Rec_Days <- as.numeric(Recov_Time_80$Rec_Days)
ggplot(Recov_Time_80, aes(y=Rec_Yrs)) +
  geom_boxplot(col="orange") + 
  geom_density(col="purple") + 
  labs(title = "Reach 80% threshold") +
  theme_bw()
# percentage that gets this far
16857 /114058 * 100 # =  14.78%
# number of fires
length(unique(Recov_Time_80$FireName)) # 129

# 70%...............................................................................................
# remove pts that recovered beyond threshold
Recov_70 <- anti_join(Recov_Combo, Recov_Time_100, by="ptID")
Recov_70 <- anti_join(Recov_70, Recov_Time_90, by="ptID")
Recov_70 <- anti_join(Recov_70, Recov_Time_80, by="ptID")
# calculate 70% of lower limit
Recov_70$lwr70 <- Recov_70$lwrNDVI * 0.7
# filter for obs at or above threshold
Recov_70 <- Recov_70[which(Recov_70$NDVI >= Recov_70$lwr70),]
# summarize by pt to get min date (first date with NDVI above threshold)
Rec_Date <- Recov_70 %>% 
  group_by(ptID) %>%
  summarise(Rec_Date = min(Obs_Date, na.rm = T)) 
# join min_date with Recov_Combo by ptID
Recov_70 <- left_join(Recov_70, Rec_Date, by="ptID")
rm(Rec_Date)
# how many pts reach target?
length(unique(Recov_70$ptID)) # 12,833
# CALCULATE RECOVERY TIME
# pull out distinct pt info
Recov_Time_70 <- Recov_70 %>% 
  dplyr::select(ptID, Rec_Date, EndDate, FireYear, Severity) %>%
  distinct()
# calculate recovery time
class(Recov_Time_70$EndDate)
Recov_Time_70$Rec_Date <- as.Date(Recov_Time_70$Rec_Date)
Recov_Time_70$EndDate <- as.Date(Recov_Time_70$EndDate)
Recov_Time_70$Rec_Days <- Recov_Time_70$Rec_Date - Recov_Time_70$EndDate
Recov_Time_70$Rec_Yrs <- as.numeric(Recov_Time_70$Rec_Days) /356
head(Recov_Time_70)
# add back other info from Recov_Combo
Recov_Time_70$Obs_Date <- Recov_Time_70$Rec_Date
names(Recov_Time_70)
Recov_Time_70 <- merge(Recov_Time_70, Recov_70, by=c("ptID","Rec_Date", "EndDate" ,"FireYear" ,"Severity", "Obs_Date"))
rm(Recov_70)

# check it out
summary(Recov_Time_70$Rec_Yrs) 
# Min.      1st Qu.      Median     Mean      3rd Qu.      Max. 
# 0.002809  0.575843   1.005618  1.732911  1.941011    18.247191 
Recov_Time_70$Rec_Days <- as.numeric(Recov_Time_70$Rec_Days)
ggplot(Recov_Time_70, aes(y=Rec_Yrs)) +
  geom_boxplot(col="orange") + 
  geom_density(col="purple") + 
  labs(title = "Reach 70% threshold") +
  theme_bw()
# percentage that gets this farss
12833 /114058 * 100 # =  11.25%
# number of fires
length(unique(Recov_Time_70$FireName)) # 99


# 60%...............................................................................................
# remove pts that recovered beyond threshold
Recov_60 <- anti_join(Recov_Combo, Recov_Time_100, by="ptID")
Recov_60 <- anti_join(Recov_60, Recov_Time_90, by="ptID")
Recov_60 <- anti_join(Recov_60, Recov_Time_80, by="ptID")
Recov_60 <- anti_join(Recov_60, Recov_Time_70, by="ptID")
# calculate 60% of lower limit
Recov_60$lwr60 <- Recov_60$lwrNDVI * 0.6
# filter for obs at or above threshold
Recov_60 <- Recov_60[which(Recov_60$NDVI >= Recov_60$lwr60),]
# summarize by pt to get min date (first date with NDVI above threshold)
Rec_Date <- Recov_60 %>% 
  group_by(ptID) %>%
  summarise(Rec_Date = min(Obs_Date, na.rm = T)) 
# join min_date with Recov_Combo by ptID
Recov_60 <- left_join(Recov_60, Rec_Date, by="ptID")
rm(Rec_Date)
# how many pts reach target?
length(unique(Recov_60$ptID)) # 4,342
# CALCULATE RECOVERY TIME
# pull out distinct pt info
Recov_Time_60 <- Recov_60 %>% 
  dplyr::select(ptID, Rec_Date, EndDate, FireYear, Severity) %>%
  distinct()
# calculate recovery time
class(Recov_Time_60$EndDate)
Recov_Time_60$Rec_Date <- as.Date(Recov_Time_60$Rec_Date)
Recov_Time_60$EndDate <- as.Date(Recov_Time_60$EndDate)
Recov_Time_60$Rec_Days <- Recov_Time_60$Rec_Date - Recov_Time_60$EndDate
Recov_Time_60$Rec_Yrs <- as.numeric(Recov_Time_60$Rec_Days) /356
head(Recov_Time_60)
# add back other info from Recov_Combo
Recov_Time_60$Obs_Date <- Recov_Time_60$Rec_Date
names(Recov_Time_60)
Recov_Time_60 <- merge(Recov_Time_60, Recov_60, by=c("ptID","Rec_Date", "EndDate" ,"FireYear" ,"Severity", "Obs_Date"))
rm(Recov_60)

# check it out
summary(Recov_Time_60$Rec_Yrs) 
# Min.      1st Qu.      Median     Mean      3rd Qu.      Max. 
# 0.002809  0.575843  0.575843   1.039068  1.134129   13.632022 
Recov_Time_60$Rec_Days <- as.numeric(Recov_Time_60$Rec_Days)
ggplot(Recov_Time_60, aes(y=Rec_Yrs)) +
  geom_boxplot(col="orange") + 
  geom_density(col="purple") + 
  labs(title = "Reach 60% threshold") +
  theme_bw()
# percentage that gets this far
4342 /114058 * 100 # =  3.81%
# number of fires
length(unique(Recov_Time_60$FireName)) # 77

# 50%...............................................................................................
# remove pts that recovered beyond threshold
Recov_50 <- anti_join(Recov_Combo, Recov_Time_100, by="ptID")
Recov_50 <- anti_join(Recov_50, Recov_Time_90, by="ptID")
Recov_50 <- anti_join(Recov_50, Recov_Time_80, by="ptID")
Recov_50 <- anti_join(Recov_50, Recov_Time_70, by="ptID")
Recov_50 <- anti_join(Recov_50, Recov_Time_60, by="ptID")
# calculate 50% of lower limit
Recov_50$lwr50 <- Recov_50$lwrNDVI * 0.5
# filter for obs at or above threshold
Recov_50 <- Recov_50[which(Recov_50$NDVI >= Recov_50$lwr50),]
# summarize by pt to get min date (first date with NDVI above threshold)
Rec_Date <- Recov_50 %>% 
  group_by(ptID) %>%
  summarise(Rec_Date = min(Obs_Date, na.rm = T)) 
# join min_date with Recov_Combo by ptID
Recov_50 <- left_join(Recov_50, Rec_Date, by="ptID")
rm(Rec_Date)
# how many pts reach target?
length(unique(Recov_50$ptID)) # 174
# CALCULATE RECOVERY TIME
# pull out distinct pt info
Recov_Time_50 <- Recov_50 %>% 
  dplyr::select(ptID, Rec_Date, EndDate, FireYear, Severity) %>%
  distinct()
# calculate recovery time
class(Recov_Time_50$EndDate)
Recov_Time_50$Rec_Date <- as.Date(Recov_Time_50$Rec_Date)
Recov_Time_50$EndDate <- as.Date(Recov_Time_50$EndDate)
Recov_Time_50$Rec_Days <- Recov_Time_50$Rec_Date - Recov_Time_50$EndDate
Recov_Time_50$Rec_Yrs <- as.numeric(Recov_Time_50$Rec_Days) /356
head(Recov_Time_50)
# add back other info from Recov_Combo
Recov_Time_50$Obs_Date <- Recov_Time_50$Rec_Date
names(Recov_Time_50)
Recov_Time_50 <- merge(Recov_Time_50, Recov_50, by=c("ptID","Rec_Date", "EndDate" ,"FireYear" ,"Severity", "Obs_Date"))
rm(Recov_50)

# check it out
summary(Recov_Time_50$Rec_Yrs) 
# Min.      1st Qu.      Median     Mean      3rd Qu.      Max. 
# 0.002809 0.201545    0.533708   0.596442   0.917135.  2.558989
Recov_Time_50$Rec_Days <- as.numeric(Recov_Time_50$Rec_Days)
ggplot(Recov_Time_50, aes(y=Rec_Yrs)) +
  geom_boxplot(col="orange") + 
  geom_density(col="purple") + 
  labs(title = "Reach 50% threshold") +
  theme_bw()
# percentage that gets this far
174 /114058 * 100 # =  0.15%
# number of fires
length(unique(Recov_Time_50$FireName)) # 38

# <50%...............................................................................................
# remove pts that recovered beyond threshold
Recov_undr50 <- anti_join(Recov_Combo, Recov_Time_100, by="ptID")
Recov_undr50 <- anti_join(Recov_undr50, Recov_Time_90, by="ptID")
Recov_undr50 <- anti_join(Recov_undr50, Recov_Time_80, by="ptID")
Recov_undr50 <- anti_join(Recov_undr50, Recov_Time_70, by="ptID")
Recov_undr50 <- anti_join(Recov_undr50, Recov_Time_60, by="ptID")
Recov_undr50 <- anti_join(Recov_undr50, Recov_Time_undr50, by="ptID")
# calculate undr50% of lower limit
Recov_undr50$lwrundr50 <- Recov_undr50$lwrNDVI * 0.5
# filter for obs at or above threshold
Recov_undr50 <- Recov_undr50[which(Recov_undr50$NDVI < Recov_undr50$lwrundr50),]
# summarize by pt to get min date (first date with NDVI above threshold)
Rec_Date <- Recov_undr50 %>% 
  group_by(ptID) %>%
  summarise(Rec_Date = min(Obs_Date, na.rm = T)) 
# join min_date with Recov_Combo by ptID
Recov_undr50 <- left_join(Recov_undr50, Rec_Date, by="ptID")
rm(Rec_Date)
# how many pts reach target?
length(unique(Recov_undr50$ptID)) # 169
# CALCULATE RECOVERY TIME
# pull out distinct pt info
Recov_Time_undr50 <- Recov_undr50 %>% 
  dplyr::select(ptID, Rec_Date, EndDate, FireYear, Severity) %>%
  distinct()
# calculate recovery time
class(Recov_Time_undr50$EndDate)
Recov_Time_undr50$Rec_Date <- as.Date(Recov_Time_undr50$Rec_Date)
Recov_Time_undr50$EndDate <- as.Date(Recov_Time_undr50$EndDate)
Recov_Time_undr50$Rec_Days <- Recov_Time_undr50$Rec_Date - Recov_Time_undr50$EndDate
Recov_Time_undr50$Rec_Yrs <- as.numeric(Recov_Time_undr50$Rec_Days) /356
head(Recov_Time_undr50)
# add back other info from Recov_Combo
Recov_Time_undr50$Obs_Date <- Recov_Time_undr50$Rec_Date
names(Recov_Time_undr50)
Recov_Time_undr50 <- merge(Recov_Time_undr50, Recov_undr50, by=c("ptID","Rec_Date", "EndDate" ,"FireYear" ,"Severity", "Obs_Date"))
rm(Recov_undr50)

# check it out
summary(Recov_Time_undr50$Rec_Yrs) 
# Min.      1st Qu.      Median     Mean      3rd Qu.      Max. 
# 0.002809  0.036517  0.126405     0.205771  0.283708   1.497191 
Recov_Time_undr50$Rec_Days <- as.numeric(Recov_Time_undr50$Rec_Days)
ggplot(Recov_Time_undr50, aes(y=Rec_Yrs)) +
  geom_boxplot(col="orange") + 
  geom_density(col="purple") + 
  labs(title = "Reach <50% threshold") +
  theme_bw()
# percentage that gets this far
169 /114058 * 100 # =  0.15%
# number of fires
length(unique(Recov_Time_undr50$FireName)) # 37

# MERGE.................................................................................................
# Select desired columns
names(Recov_Time_100)
Recov_Time_100 <- Recov_Time_100 %>%
  select(ptID, coords.x1, coords.x2, StartDate, EndDate, FireYear, FireNumber, 
         FireType, TotalFires, Prev.Int, Severity, PostDateDif, model.NDVI, NDVI, Obs_Date, Rec_Yrs)
Recov_Time_90 <- Recov_Time_90 %>%
  select(ptID, coords.x1, coords.x2, StartDate, EndDate, FireYear, FireNumber, 
         FireType, TotalFires, Prev.Int, Severity, PostDateDif, model.NDVI, NDVI, Obs_Date, Rec_Yrs)
Recov_Time_80 <- Recov_Time_80 %>%
  select(ptID, coords.x1, coords.x2, StartDate, EndDate, FireYear, FireNumber, 
         FireType, TotalFires, Prev.Int, Severity, PostDateDif, model.NDVI, NDVI, Obs_Date, Rec_Yrs)
Recov_Time_70 <- Recov_Time_70 %>%
  select(ptID, coords.x1, coords.x2, StartDate, EndDate, FireYear, FireNumber, 
         FireType, TotalFires, Prev.Int, Severity, PostDateDif, model.NDVI, NDVI, Obs_Date, Rec_Yrs)
Recov_Time_60 <- Recov_Time_60 %>%
  select(ptID, coords.x1, coords.x2, StartDate, EndDate, FireYear, FireNumber, 
         FireType, TotalFires, Prev.Int, Severity, PostDateDif, model.NDVI, NDVI, Obs_Date, Rec_Yrs)
Recov_Time_50 <- Recov_Time_50 %>%
  select(ptID, coords.x1, coords.x2, StartDate, EndDate, FireYear, FireNumber, 
         FireType, TotalFires, Prev.Int, Severity, PostDateDif, model.NDVI, NDVI, Obs_Date, Rec_Yrs)
Recov_Time_undr50 <- Recov_Time_undr50 %>%
  select(ptID, coords.x1, coords.x2, StartDate, EndDate, FireYear, FireNumber, 
         FireType, TotalFires, Prev.Int, Severity, PostDateDif, model.NDVI, NDVI, Obs_Date, Rec_Yrs)

# add coumn for threshold
Recov_Time_100$threshold <- "100"
Recov_Time_90$threshold <- "90"
Recov_Time_80$threshold <- "80"
Recov_Time_70$threshold <- "70"
Recov_Time_60$threshold <- "60"
Recov_Time_50$threshold <- "50"
Recov_Time_undr50$threshold <- "<50"

# merge together
Rate_drivers <- rbind(Recov_Time_100, Recov_Time_90)
Rate_drivers <- rbind(Rate_drivers, Recov_Time_80)
Rate_drivers <- rbind(Rate_drivers, Recov_Time_70)
Rate_drivers <- rbind(Rate_drivers, Recov_Time_60)
Rate_drivers <- rbind(Rate_drivers, Recov_Time_50)
Rate_drivers <- rbind(Rate_drivers, Recov_Time_undr50)
head(Rate_drivers)
unique(Rate_drivers$threshold)

# make threshold a factor
Rate_drivers$threshold <- as.factor(Rate_drivers$threshold)

# SAVE..........................................................................................................
save(Recov_Time_100, Recov_Time_90, Recov_Time_80, Recov_Time_70, Recov_Time_60, Recov_Time_50, Recov_Time_undr50, file = "Recov_time_thresholds.RDATA")
save(Rate_drivers, file="Rate_drivers.RDATA")

##########################################################################################################################################################
# RECOVERY TIME TO X% OF THRESHOLD (and beyond)....use for recovery rate
##########################################################################################################################################################
# Who REACHES the threshold AT ALL 
load("./Recovery/Recov_Combo2.RDATA")

# 100% (reach lower limit of baseline)...............................................................................................
# pull out all obs with NDVI > threshold  
Recov_100 <- Recov_Combo[which(Recov_Combo$NDVI >= (Recov_Combo$lwrNDVI * 1)),]
# summarize by pt to get min date (first date with NDVI above threshold)
Rec_Date <- Recov_100 %>% 
  group_by(ptID) %>%
  summarise(Rec100_Date = min(Obs_Date, na.rm = T)) 
# join min_date with Recov_threshold info by ptID
Recov_100 <- left_join(Recov_100, Rec_Date, by="ptID")
rm(Rec_Date)

# CALCULATE RECOVERY TIME
# between fire end date and rec date
class(Recov_100$Rec100_Date)
Recov_100$Rec100_Date <- as.Date(Recov_100$Rec100_Date)
Recov_100$EndDate <- as.Date(Recov_100$EndDate)
Recov_100$Rec100_Yrs <- as.numeric((Recov_100$Rec100_Date - Recov_100$EndDate) /356)
# Keep only the recovery observation for each point
Recov_100 <- Recov_100[which(Recov_100$Obs_Date == Recov_100$Rec100_Date),]
Recov_100$Rec100_NDVI <- Recov_100$NDVI

# check it out
summary(Recov_100$Rec100_Yrs) 
# Min.      1st Qu.      Median     Mean      3rd Qu.      Max. 
# 0.002809  0.741573  1.240169   2.096887  2.308989   20.087079 
ggplot(Recov_100, aes(y=Rec100_Yrs)) +
  geom_boxplot(col="orange") + 
  geom_density(col="purple") + 
  labs(title = "Recover to Target (100%)") +
  theme_bw()
# how many pts reach target?
length(Recov_100$ptID) # 60,296
# percentage that gets this far
length(Recov_100$ptID) / 114058 * 100 # = 52.86% 
# number of fires
length(unique(Recov_Combo$FireName)) # 178
length(unique(Recov_100$FireName)) # 162

# 90% ...............................................................................................
Recov_90 <- Recov_Combo[which(Recov_Combo$NDVI >= (Recov_Combo$lwrNDVI * .9)),]
Rec_Date <- Recov_90 %>% 
  group_by(ptID) %>%
  summarise(Rec90_Date = min(Obs_Date, na.rm = T)) 
Recov_90 <- left_join(Recov_90, Rec_Date, by="ptID")
rm(Rec_Date)
Recov_90$Rec90_Date <- as.Date(Recov_90$Rec90_Date)
Recov_90$EndDate <- as.Date(Recov_90$EndDate)
Recov_90$Rec90_Yrs <- as.numeric((Recov_90$Rec90_Date - Recov_90$EndDate) /356)
Recov_90 <- Recov_90[which(Recov_90$Obs_Date == Recov_90$Rec90_Date),]
Recov_90$Rec90_NDVI <- Recov_90$NDVI
summary(Recov_90$Rec90_Yrs) 
# Min.      1st Qu.      Median     Mean      3rd Qu.      Max. 
# 0.002809  0.469101   0.941011  1.431211  1.429775    19.862360 
ggplot(Recov_90, aes(y=Rec90_Yrs)) +
  geom_boxplot(col="orange") + 
  geom_density(col="purple") + 
  labs(title = "Recover to Target (90%)") +
  theme_bw()
length(Recov_90$ptID) # 79,849
length(Recov_90$ptID) / 114058 * 100 # = 70.01%
length(unique(Recov_90$FireName)) # 168

# 80% ...............................................................................................
Recov_80 <- Recov_Combo[which(Recov_Combo$NDVI >= (Recov_Combo$lwrNDVI * .8)),]
Rec_Date <- Recov_80 %>% 
  group_by(ptID) %>%
  summarise(Rec80_Date = min(Obs_Date, na.rm = T)) 
Recov_80 <- left_join(Recov_80, Rec_Date, by="ptID")
rm(Rec_Date)
Recov_80$Rec80_Date <- as.Date(Recov_80$Rec80_Date)
Recov_80$EndDate <- as.Date(Recov_80$EndDate)
Recov_80$Rec80_Yrs <- as.numeric((Recov_80$Rec80_Date - Recov_80$EndDate) /356)
Recov_80 <- Recov_80[which(Recov_80$Obs_Date == Recov_80$Rec80_Date),]
Recov_80$Rec80_NDVI <- Recov_80$NDVI
summary(Recov_80$Rec80_Yrs) 
# Min.      1st Qu.      Median     Mean      3rd Qu.      Max. 
# 0.002809  0.306180  0.575843   0.951633  1.039326   19.342697 
ggplot(Recov_80, aes(y=Rec80_Yrs)) +
  geom_boxplot(col="orange") + 
  geom_density(col="purple") + 
  labs(title = "Recover to Target (80%)") +
  theme_bw()
length(Recov_80$ptID) # 96,706
length(Recov_80$ptID) / 114058 * 100 # = 84.79%
length(unique(Recov_80$FireName)) # 177

# 70% ...............................................................................................
Recov_70 <- Recov_Combo[which(Recov_Combo$NDVI >= (Recov_Combo$lwrNDVI * .7)),]
Rec_Date <- Recov_70 %>% 
  group_by(ptID) %>%
  summarise(Rec70_Date = min(Obs_Date, na.rm = T)) 
Recov_70 <- left_join(Recov_70, Rec_Date, by="ptID")
rm(Rec_Date)
Recov_70$Rec70_Date <- as.Date(Recov_70$Rec70_Date)
Recov_70$EndDate <- as.Date(Recov_70$EndDate)
Recov_70$Rec70_Yrs <- as.numeric((Recov_70$Rec70_Date - Recov_70$EndDate) /356)
Recov_70 <- Recov_70[which(Recov_70$Obs_Date == Recov_70$Rec70_Date),]
Recov_70$Rec70_NDVI <- Recov_70$NDVI
summary(Recov_70$Rec70_Yrs) 
# Min.      1st Qu.      Median     Mean      3rd Qu.      Max. 
# 0.002809  0.193820  0.370787   0.594367  0.643258    18.247191
ggplot(Recov_70, aes(y=Rec70_Yrs)) +
  geom_boxplot(col="orange") + 
  geom_density(col="purple") + 
  labs(title = "Recover to Target (70%)") +
  theme_bw()
length(Recov_70$ptID) # 109,539
length(Recov_70$ptID) / 114058 * 100 # = 96.04%
length(unique(Recov_70$FireName)) # 178

# 60% ...............................................................................................
Recov_60 <- Recov_Combo[which(Recov_Combo$NDVI >= (Recov_Combo$lwrNDVI * .6)),]
Rec_Date <- Recov_60 %>% 
  group_by(ptID) %>%
  summarise(Rec60_Date = min(Obs_Date, na.rm = T)) 
Recov_60 <- left_join(Recov_60, Rec_Date, by="ptID")
rm(Rec_Date)
Recov_60$Rec60_Date <- as.Date(Recov_60$Rec60_Date)
Recov_60$EndDate <- as.Date(Recov_60$EndDate)
Recov_60$Rec60_Yrs <- as.numeric((Recov_60$Rec60_Date - Recov_60$EndDate) /356)
Recov_60 <- Recov_60[which(Recov_60$Obs_Date == Recov_60$Rec60_Date),]
Recov_60$Rec60_NDVI <- Recov_60$NDVI
summary(Recov_60$Rec60_Yrs) 
# Min.      1st Qu.      Median     Mean      3rd Qu.      Max. 
# 0.002809  0.117978  0.230337  0.324307  0.441011      13.632022 
ggplot(Recov_60, aes(y=Rec60_Yrs)) +
  geom_boxplot(col="orange") + 
  geom_density(col="purple") + 
  labs(title = "Recover to Target (60%)") +
  theme_bw()
length(Recov_60$ptID) # 113,881
length(Recov_60$ptID) / 114058 * 100 # = 99.84%
length(unique(Recov_60$FireName)) # 178

# 50% ...............................................................................................
Recov_50 <- Recov_Combo[which(Recov_Combo$NDVI >= (Recov_Combo$lwrNDVI * .5)),]
Rec_Date <- Recov_50 %>% 
  group_by(ptID) %>%
  summarise(Rec50_Date = min(Obs_Date, na.rm = T)) 
Recov_50 <- left_join(Recov_50, Rec_Date, by="ptID")
rm(Rec_Date)
Recov_50$Rec50_Date <- as.Date(Recov_50$Rec50_Date)
Recov_50$EndDate <- as.Date(Recov_50$EndDate)
Recov_50$Rec50_Yrs <- as.numeric((Recov_50$Rec50_Date - Recov_50$EndDate) /356)
Recov_50 <- Recov_50[which(Recov_50$Obs_Date == Recov_50$Rec50_Date),]
Recov_50$Rec50_NDVI <- Recov_50$NDVI
summary(Recov_50$Rec50_Yrs) 
# Min.      1st Qu.      Median     Mean      3rd Qu.      Max. 
# 0.002809  0.058989  0.160112   0.196052   0.238764    4.603933 
ggplot(Recov_50, aes(y=Rec50_Yrs)) +
  geom_boxplot(col="orange") + 
  geom_density(col="purple") + 
  labs(title = "Recover to Target (50%)") +
  theme_bw()
length(Recov_50$ptID) # 114,055
length(Recov_50$ptID) / 114058 * 100 # = 99.99%
length(unique(Recov_50$FireName)) # 178

# 40% ...............................................................................................
Recov_40 <- Recov_Combo[which(Recov_Combo$NDVI >= (Recov_Combo$lwrNDVI * .4)),]
Rec_Date <- Recov_40 %>% 
  group_by(ptID) %>%
  summarise(Rec40_Date = min(Obs_Date, na.rm = T)) 
Recov_40 <- left_join(Recov_40, Rec_Date, by="ptID")
rm(Rec_Date)
Recov_40$Rec40_Date <- as.Date(Recov_40$Rec40_Date)
Recov_40$EndDate <- as.Date(Recov_40$EndDate)
Recov_40$Rec40_Yrs <- as.numeric((Recov_40$Rec40_Date - Recov_40$EndDate) /356)
Recov_40 <- Recov_40[which(Recov_40$Obs_Date == Recov_40$Rec40_Date),]
Recov_40$Rec40_NDVI <- Recov_40$NDVI
summary(Recov_40$Rec40_Yrs) 
# Min.      1st Qu.      Median     Mean      3rd Qu.      Max. 
# 0.002809  0.036517  0.092697  0.139926   0.193820     2.356742
ggplot(Recov_40, aes(y=Rec40_Yrs)) +
  geom_boxplot(col="orange") + 
  geom_density(col="purple") + 
  labs(title = "Recover to Target (40%)") +
  theme_bw()
length(Recov_40$ptID) # 114,058 
length(Recov_40$ptID) / 114058 * 100 # = 100%

# 30% ...............................................................................................
Recov_30 <- Recov_Combo[which(Recov_Combo$NDVI >= (Recov_Combo$lwrNDVI * .3)),]
Rec_Date <- Recov_30 %>% 
  group_by(ptID) %>%
  summarise(Rec30_Date = min(Obs_Date, na.rm = T)) 
Recov_30 <- left_join(Recov_30, Rec_Date, by="ptID")
rm(Rec_Date)
Recov_30$Rec30_Date <- as.Date(Recov_30$Rec30_Date)
Recov_30$EndDate <- as.Date(Recov_30$EndDate)
Recov_30$Rec30_Yrs <- as.numeric((Recov_30$Rec30_Date - Recov_30$EndDate) /356)
Recov_30 <- Recov_30[which(Recov_30$Obs_Date == Recov_30$Rec30_Date),]
Recov_30$Rec30_NDVI <- Recov_30$NDVI
summary(Recov_30$Rec30_Yrs) 
# Min.      1st Qu.      Median     Mean      3rd Qu.      Max. 
# 0.002809  0.036517  0.087079   0.121807   0.179775    1.691011 
ggplot(Recov_30, aes(y=Rec30_Yrs)) +
  geom_boxplot(col="orange") + 
  geom_density(col="purple") + 
  labs(title = "Recover to Target (30%)") +
  theme_bw()

# 20% ...............................................................................................
Recov_20 <- Recov_Combo[which(Recov_Combo$NDVI >= (Recov_Combo$lwrNDVI * .2)),]
Rec_Date <- Recov_20 %>% 
  group_by(ptID) %>%
  summarise(Rec20_Date = min(Obs_Date, na.rm = T)) 
Recov_20 <- left_join(Recov_20, Rec_Date, by="ptID")
rm(Rec_Date)
Recov_20$Rec20_Date <- as.Date(Recov_20$Rec20_Date)
Recov_20$EndDate <- as.Date(Recov_20$EndDate)
Recov_20$Rec20_Yrs <- as.numeric((Recov_20$Rec20_Date - Recov_20$EndDate) /356)
Recov_20 <- Recov_20[which(Recov_20$Obs_Date == Recov_20$Rec20_Date),]
Recov_20$Rec20_NDVI <- Recov_20$NDVI
summary(Recov_20$Rec20_Yrs) 
# Min.      1st Qu.      Median     Mean      3rd Qu.      Max. 
# 0.002809  0.036517   0.081461   0.117977   0.179775    1.691011
ggplot(Recov_20, aes(y=Rec20_Yrs)) +
  geom_boxplot(col="orange") + 
  geom_density(col="purple") + 
  labs(title = "Recover to Target (20%)") +
  theme_bw()

# 10% ...............................................................................................
Recov_10 <- Recov_Combo[which(Recov_Combo$NDVI >= (Recov_Combo$lwrNDVI * .1)),]
Rec_Date <- Recov_10 %>% 
  group_by(ptID) %>%
  summarise(Rec10_Date = min(Obs_Date, na.rm = T)) 
Recov_10 <- left_join(Recov_10, Rec_Date, by="ptID")
rm(Rec_Date)
Recov_10$Rec10_Date <- as.Date(Recov_10$Rec10_Date)
Recov_10$EndDate <- as.Date(Recov_10$EndDate)
Recov_10$Rec10_Yrs <- as.numeric((Recov_10$Rec10_Date - Recov_10$EndDate) /356)
Recov_10 <- Recov_10[which(Recov_10$Obs_Date == Recov_10$Rec10_Date),]
Recov_10$Rec10_NDVI <- Recov_10$NDVI
summary(Recov_10$Rec10_Yrs) 
# Min.      1st Qu.      Median     Mean      3rd Qu.      Max. 
# 0.002809 0.036517  0.081461    0.117687   0.179775    1.691011 
ggplot(Recov_10, aes(y=Rec10_Yrs)) +
  geom_boxplot(col="orange") + 
  geom_density(col="purple") + 
  labs(title = "Recover to Target (10%)") +
  theme_bw()

# Save 
setwd("./Recovery")
save(Recov_100, Recov_90, Recov_80, Recov_70, Recov_60, Recov_50, Recov_40,
     Recov_30, Recov_20, Recov_10, file = "Recov_thresholds.RDATA")


# Make Recovery Rate Dataframe.....................................................................................

# remove columns that correlate to only one observation
Recov_10 <- Recov_10 %>%
  select(ptID, coords.x1, coords.x2,
         StartDate, EndDate, FireYear, PenUltFY, AFyear, FireName, FireNumber, 
         FireType, TotalFires, Prev.Int, PreNBR, PreNDVI, 
         PostNBR, PostNDVI, PreDate, PostDate, DateDif, Severity, PreDateDif, PostDateDif,
         Pt.B5max, model.NDVI, var.NDVI, lwrNDVI, uprNDVI,  
         Rec10_Date, Rec10_Yrs, Rec10_NDVI)
Recov_20 <- Recov_20 %>%
  select(ptID, coords.x1, coords.x2,
         Rec20_Date, Rec20_Yrs, Rec20_NDVI)
Recov_30 <- Recov_30 %>%
  select(ptID, coords.x1, coords.x2,
         Rec30_Date, Rec30_Yrs, Rec30_NDVI)
Recov_40 <- Recov_40 %>%
  select(ptID, coords.x1, coords.x2,  
         Rec40_Date, Rec40_Yrs, Rec40_NDVI)
Recov_50 <- Recov_50 %>%
  select(ptID, coords.x1, coords.x2,  
         Rec50_Date, Rec50_Yrs, Rec50_NDVI)
Recov_60 <- Recov_60 %>%
  select(ptID, coords.x1, coords.x2,  
         Rec60_Date, Rec60_Yrs, Rec60_NDVI)
Recov_70 <- Recov_70 %>%
  select(ptID, coords.x1, coords.x2, 
         Rec70_Date, Rec70_Yrs, Rec70_NDVI)
Recov_80 <- Recov_80 %>%
  select(ptID, coords.x1, coords.x2,  
         Rec80_Date, Rec80_Yrs, Rec80_NDVI)
Recov_90 <- Recov_90 %>%
  select(ptID, coords.x1, coords.x2, 
         Rec90_Date, Rec90_Yrs, Rec90_NDVI)
Recov_100 <- Recov_100 %>%
  select(ptID, coords.x1, coords.x2, 
         Rec100_Date, Rec100_Yrs, Rec100_NDVI)

# merge threhold dataframes
Recov_Rate <- merge(Recov_10, Recov_20, all.x=T)
Recov_Rate <- merge(Recov_Rate, Recov_30, all.x=T)
Recov_Rate <- merge(Recov_Rate, Recov_40, all.x=T)
Recov_Rate <- merge(Recov_Rate, Recov_50, all.x=T)
Recov_Rate <- merge(Recov_Rate, Recov_60, all.x=T)
Recov_Rate <- merge(Recov_Rate, Recov_70, all.x=T)
Recov_Rate <- merge(Recov_Rate, Recov_80, all.x=T)
Recov_Rate <- merge(Recov_Rate, Recov_90, all.x=T)
Recov_Rate <- merge(Recov_Rate, Recov_100, all.x=T)
head(Recov_Rate)

# Save
save(Recov_Rate, file="Recov_Rates.RDATA")


##########################################################################################################################################################
# POINT-SPECIFIC MAX RECOVERY
##########################################################################################################################################################
# identify "max recovery"
# pull out "recovery date"
# make summary dataframe for analysis

# load all observation data
load("Y:/Recovery/Recov_Combo2.RDATA")

# identify "max recovery"
# minimum difference between baseline and observed
Recov_Combo$RecDif <- Recov_Combo$model.NDVI - Recov_Combo$NDVI # difference between baseline and obs NDVI
min.dif <- Recov_Combo %>%   
  group_by(ptID) %>%
  summarise(min.dif = min(model.NDVI - NDVI),
            MaxRec_Date = Obs_Date[RecDif == min.dif])
Recov_Combo <- merge(Recov_Combo, min.dif, by="ptID") 
# remove observations after max recovery date
Recov_Combo$Obs_Date <- as.Date(Recov_Combo$Obs_Date)
Recov_Combo$MaxRec_Date <- as.Date(Recov_Combo$MaxRec_Date)
Pre_Rec_obs <- Recov_Combo[which(Recov_Combo$Obs_Date <= Recov_Combo$MaxRec_Date),]

# calculate time to "max recovery"
Pre_Rec_obs$StartDate <- as.Date(Pre_Rec_obs$StartDate)
Pre_Rec_obs$Obs_Date <- as.Date(Pre_Rec_obs$Obs_Date)
Pre_Rec_obs$Rec_Yrs <- Pre_Rec_obs$Obs_Date - Pre_Rec_obs$StartDate
Pre_Rec_obs$Rec_Yrs <- as.numeric(Pre_Rec_obs$Rec_Yrs)
Pre_Rec_obs$Rec_Yrs <- Pre_Rec_obs$Rec_Yrs / 365

# keep only max rec obs
Rec_Time <- Recov_Combo[which(Recov_Combo$RecDif == Recov_Combo$min.dif),] 

# SAVE
setwd("./Recovery")
save(Pre_Rec_obs, Rec_Time, file="MaxRec.RDATA")

##########################################################################################################################################################
# STAY RECOVERED
##########################################################################################################################################################
load("I:/Malone Lab/ENP Fire/Grace_McLeod/Recovery/Recov_Combo2.RDATA")
load("H:/Recovery/Recov_Combo2.RDATA")

# SET THRESHOLD.........................................................................................................................................
# For pts that reach X% of target
# pull out all obs with NDVI > threshold 

# X %
Recov_Combo$lwrThresh <- Recov_Combo$lwrNDVI *.7
threshold <- Recov_Combo[which(Recov_Combo$NDVI >= Recov_Combo$lwrThresh),]

# how many pts reach target?
length(unique(threshold$ptID)) 
# summarize by pt to get min date (first date with NDVI above threshold)
threshold$Obs_Date <- as.Date(threshold$Obs_Date)
Min_Date <- threshold %>% 
  group_by(ptID) %>%
  summarise(Min_Date = min(Obs_Date, na.rm = T)) 
# subset all observations to just those pts
Recov_Combo_sub <- subset(Recov_Combo, ptID %in% threshold$ptID)
length(unique(Recov_Combo_sub$ptID)) # 60,296
rm(threshold)
# join min_date with Recov_Combo_sub by ptID
Recov_Combo_sub <- left_join(Recov_Combo_sub, Min_Date, by="ptID")
summary(Recov_Combo_sub$Min_Date)
rm(Min_Date)

# identify date of last observation
Recov_Combo_sub$Obs_Date <- as.Date(Recov_Combo_sub$Obs_Date)
Last_obs_Date <- Recov_Combo_sub %>% 
  group_by(ptID) %>%
  summarise(Last_obs_Date = max(Obs_Date, na.rm = T))
# join last obs date with Recov_Combo by ptID
Recov_Combo_sub <- left_join(Recov_Combo_sub, Last_obs_Date, by="ptID")
rm(Last_obs_Date)
# Subset to only obs after initial recovery date 
Recov_Combo_sub$Obs_Date <- as.Date(Recov_Combo_sub$Obs_Date)
Recov_Combo_sub$Min_Date <- as.Date(Recov_Combo_sub$Min_Date)
Recov_Combo_sub <- Recov_Combo_sub[which(Recov_Combo_sub$Obs_Date >= Recov_Combo_sub$Min_Date),]

# new columns
Recov_Combo_sub$Rec_Date <- NA
Recov_Combo_sub$n_rec_obs <- NA
Recov_Combo_sub$prcnt_rec_obs <- NA
Recov_Combo_sub$t_rec_obs <- NA

# New dataframe
Recov_Time_NDVI <- data.frame(matrix(nrow = 0, ncol = 68))
colnames(Recov_Time_NDVI) <- colnames(Recov_Combo_sub)


# STAY RECOVERED for x% of obs ..................................................................
# cross at x% threshold and stay above for x% of observations
# must stay recovered for at least x months

# Loop
for (i in unique(Recov_Combo_sub$ptID)) {
  print(i)
  # subset to i
  subID <- Recov_Combo_sub[which(Recov_Combo_sub$ptID == i),]
  # index all dates above threshold
  dates.above <- subID$Obs_Date[subID$NDVI >= subID$lwrThresh]
  # loop through possible recovery dates
  for (d in 1:length(dates.above)){
    # subset to obs after min_date[d]
    subDate <- subID[which(subID$Obs_Date >= d),]
    # identify total number of obs
    n.obs <- as.numeric(length(subDate$NDVI))
    # identify number of obs above threshold
    n.above <- as.numeric(length(subDate$NDVI[subDate$NDVI >= subDate$lwrThresh]))
    # identify possible time spent recovered
    t.above <- (interval((subID$Obs_Date[d]), (unique(subID$Last_obs_Date))) %/% months(1))
    # if desired % of obs over threshold....
    if ((n.above / n.obs *100) >= 80 & t.above >= 3){
      # record percentage of obs used to decide
      subID$prcnt_rec_obs <- as.numeric(n.above / n.obs *100)
      # record number of observation used
      subID$n_rec_obs <- n.obs
      # record time spent recovered
      subID$t_rec_obs <- t.above
      # write out d as recovery date
      subID$Rec_Date <- subID$Obs_Date[d]
      RecRow <- subID[which(subID$Rec_Date == subID$Obs_Date),]
      Recov_Time_NDVI <- smartbind(RecRow, Recov_Time_NDVI)
      rm(RecRow)
      break
    } else {
      print("not recovered")
    }
  }
}


# CHECK IT OUT

# 100% threshold, 100% of obs, min 3mo
length(unique(Recov_Combo$ptID)) # 114,058 total pts
length(unique(Recov_Combo_sub$ptID)) # 60,296 pts cross threshold (53%)
length(Recov_Time_NDVI$ptID) # 69 pts stay recovered (0.6%)
summary(Recov_Time_NDVI$n_rec_obs) # most had ~ 10 obs in that period
# save 
setwd("H:/Recovery")
save(Recov_Time_NDVI, file="Rec_100_100%_3mo.RDATA")
# 100% threshold, 80% of obs, min 3mo
length(unique(Recov_Combo$ptID)) # 114,058 total pts
length(unique(Recov_Combo_sub$ptID)) # 60,296 pts cross threshold (53%)
length(Recov_Time_NDVI$ptID) # 1,552 pts stay recovered (1.4%)
summary(Recov_Time_NDVI$n_rec_obs) # most had ~50  obs in that period
# save 
setwd("H:/Recovery")
save(Recov_Time_NDVI, file="Rec_100_80%_3mo.RDATA")

# 80% threshold, 100% of obs, min 3mo
length(unique(Recov_Combo$ptID)) # 114,058 total pts
length(unique(Recov_Combo_sub$ptID)) # 96,706 pts cross threshold (85%)
length(Recov_Time_NDVI$ptID) # 3,398 pts stay recovered (3%)
summary(Recov_Time_NDVI$n_rec_obs) # most had ~25 obs in that period
# save 
setwd("H:/Recovery")
save(Recov_Time_NDVI, file="Rec_80_100%_3mo.RDATA")
# 80% threshold, 80% of obs, min 3mo
length(unique(Recov_Combo$ptID)) # 114,058 total pts
length(unique(Recov_Combo_sub$ptID)) # 96,706 cross threshold (85%)
length(Recov_Time_NDVI$ptID) # 36,084 pts stay recovered (32%)
summary(Recov_Time_NDVI$n_rec_obs) # most had ~50 obs in that period
# save 
setwd("H:/Recovery")
save(Recov_Time_NDVI, file="Rec_80_80%_3mo.RDATA")

# 70% threshold, 100% of obs, min 3mo
length(unique(Recov_Combo$ptID)) # 114,058 total pts
length(unique(Recov_Combo_sub$ptID)) # 109,539 cross threshold (96%)
length(Recov_Time_NDVI$ptID) # 10,702 pts stay recovered (9%)
summary(Recov_Time_NDVI$n_rec_obs) # most had ~35 obs in that period
# save 
setwd("H:/Recovery")
save(Recov_Time_NDVI, file="Rec_70_100%_3mo.RDATA")
# 70% threshold, 80% of obs, min 3mo
length(unique(Recov_Combo$ptID)) # 114,058 total pts
length(unique(Recov_Combo_sub$ptID)) # 109,539 cross threshold (96%)
length(Recov_Time_NDVI$ptID) # 66,201 pts stay recovered (54%)
summary(Recov_Time_NDVI$n_rec_obs) # most had ~50 obs in that period
# save 
setwd("H:/Recovery")
save(Recov_Time_NDVI, file="Rec_70_80%_3mo.RDATA")


######################################################################################################################################
# EXPLORING TRENDS
#####################################################################################################################################

# IM POST-FIRE GREEN UP ISSUE???
# fresh grass is super green and can inflate NDVI values until other vegetation shades them out/scenes in fall.
# ex: Lacotoure et al. NDVI recovery in 10 days. 
# What if we identify the first two times the obs_ndvi crosses the threshold (from being below), 
# and see what the time lag generally is between the two?

# Or plot the data and look at the trends. Novel...................................................................................

# add column for threshold reached
Recov_Combo$lwrThresh <- Recov_Combo$lwrNDVI *1
Recov_Combo$uprThresh[Recov_Combo$NDVI >= Recov_Combo$lwrThresh] <- 100
Recov_Combo$uprThresh <- as.factor(Recov_Combo$uprThresh)

# calcualte time since fire for each observation
Recov_Combo$TSF <- Recov_Combo$Obs_Year - Recov_Combo$FireYear
# plot 
ggplot(Recov_Combo, aes(x=TSF, 
                          y=NDVI,
                          group=uprThresh,
                          color=uprThresh)) +
  geom_smooth() +
  labs(
    x="Years Post-Fire",
    y= "NDVI",
    title = "Recovery Locations",
    color="Upr Threshold", 
    fill= "Upr Threshold") + 
  #scale_fill_manual(values=cols) +
  #scale_color_manual(values=cols)+
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE) +
  theme_bw()+
  theme(
    text = element_text(size = 20)
  )


# Remove obs prior to initial rec date
getwd()
load("./Recov_time_thresholds.RDATA")
# pull out rec_date and ptID
Rec.Dates <- Recov_Time_80 %>%
  select(ptID, Rec_Date)
Recov_Combo$Rec_Date <- NA
Recov_Combo <- merge(Recov_Combo, Rec.Dates, by="ptID")


##########################################################################################################################################################
# UNRECOVERED
##########################################################################################################################################################

load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Manuscript/DriversData.RDATA")

# Compare mean 
means <- pdsi.total.summary.FH.maxThreshold %>%
  group_by(Threshold) %>%
  summarise(thresh.ndvi.mean = mean(NDVI))
mean(pdsi.total.summary.FH.maxThreshold$PreNDVI)
# 80% threshold = pre-fire NDVI 


# WHO DOES NOT MAKE IT TO 80?
# pull out locations with max threshold < 80
unrecovered <- pdsi.total.summary.FH.maxThreshold
unrecovered$rec.status <- "unrecovered"
unrecovered$rec.status[unrecovered$Threshold >= 80] <- "recovered"
length(unrecovered$ptID[unrecovered$rec.status == "unrecovered"]) # 17,524
length(unrecovered$ptID[unrecovered$rec.status == "recovered"]) # 96,700

# plot to see distributions
names(unrecovered)
# interesting
ggplot(unrecovered) +
  geom_density(aes(x=Prev.Int, color=rec.status)) # spike in unrecovered with low prev.int (~ 5yrs). So hard to recover if just burned
ggplot(unrecovered) +
  geom_density(aes(x=PreNDVI, color=rec.status)) # pre-NDVI is lower for unrecovered
ggplot(unrecovered) +
  geom_density(aes(x=pdsi.mean, color=rec.status)) # unrecovered has bigger hump for wet... but spread is similar
ggplot(unrecovered) +
  geom_density(aes(x=pdsi.min, color=rec.status)) # recovered has bigger hump for dry...but spread is similar
ggplot(unrecovered) +
  geom_density(aes(x=FireYear, color=rec.status)) # similar. More unrecovered fires in 2002. More recovered fires in 2007
# no pattern 
ggplot(unrecovered) +
  geom_density(aes(x=Obs_Mo, color=rec.status)) # both recovred at all times of year
ggplot(unrecovered) +
  geom_density(aes(x=pdsi.max, color=rec.status)) # same
ggplot(unrecovered) +
  geom_density(aes(x=Severity, color=rec.status)) # not much, maybe unrecovered had more locations with lower severity
ggplot(unrecovered) +
  geom_density(aes(x=TotalFires, color=rec.status)) # no real pattern

##########################################################################################################################################################
# TRASH / EXTRA
##########################################################################################################################################################

# Make dataframe for Sparkle
load("./Recov_Combo2.RDATA")
# OG values for baseline not included....add back in 
Recov_Combo$TotalFires <- Recov_Combo$TotalFires.og ; Recov_Combo$TotalFires.og <- NULL ;  Recov_Combo$TotalFires.sub <- NULL
Recov_Combo$Prev.Int <- Recov_Combo$Prev.Int.og ; Recov_Combo$Prev.Int.og <- NULL ; Recov_Combo$Prev.Int.sub <- NULL
Recov_Combo$Pt.B4min <- Recov_Combo$Pt.B4min.og ; Recov_Combo$Pt.B4min.og <- NULL ;  Recov_Combo$Pt.B4min.sub <- NULL 
Recov_Combo$Pt.B4max <- Recov_Combo$Pt.B4max.og ; Recov_Combo$Pt.B4max.og <- NULL ;  Recov_Combo$Pt.B4max.sub <- NULL 
Recov_Combo$Pt.B5max <- Recov_Combo$Pt.B5max.og ; Recov_Combo$Pt.B5max.og <- NULL ;  Recov_Combo$Pt.B5max.sub <- NULL 
Recov_Combo$Pt.B7max <- Recov_Combo$Pt.B7max.og ; Recov_Combo$Pt.B7max.og <- NULL ;  Recov_Combo$Pt.B7max.sub <- NULL 
Recov_Combo$SWIR1.SWIR2 <- Recov_Combo$SWIR1.SWIR2.og ; Recov_Combo$SWIR1.SWIR2.og <- NULL ; Recov_Combo$SWIR1.SWIR2.sub <- NULL
Recov_Combo$NIR.SWIR1 <- Recov_Combo$NIR.SWIR1.og ; Recov_Combo$NIR.SWIR1.og <- NULL ; Recov_Combo$NIR.SWIR1.sub <- NULL
# reduce to only necesary columns
names(Recov_Combo)
head(Recov_Combo)
Recov_Combo$NDVI_rf <- NULL
# 





# 80 % (Test 4)
Recov_Combo$lwr80 <- Recov_Combo$lwrNDVI *.8
threshold <- Recov_Combo[which(Recov_Combo$NDVI >= Recov_Combo$lwr80),]
# how many pts reach target?
length(unique(threshold$ptID)) # 96,706
# summarize by pt to get min date (first date with NDVI above threshold)
threshold$Obs_Date <- as.Date(threshold$Obs_Date)
Min_Date <- threshold %>% 
  group_by(ptID) %>%
  summarise(Min_Date = min(Obs_Date, na.rm = T)) 
# subset all observations to just those pts
Recov_Combo_sub <- subset(Recov_Combo, ptID %in% threshold$ptID)
length(unique(Recov_Combo_sub$ptID)) # 60,296
rm(threshold)
# join min_date with Recov_Combo_sub by ptID
Recov_Combo_sub <- left_join(Recov_Combo_sub, Min_Date, by="ptID")
summary(Recov_Combo_sub$Min_Date)
rm(Min_Date)
# identify date of last observation
Recov_Combo_sub$Obs_Date <- as.Date(Recov_Combo_sub$Obs_Date)
Last_obs_Date <- Recov_Combo_sub %>% 
  group_by(ptID) %>%
  summarise(Last_obs_Date = max(Obs_Date, na.rm = T))
# join last obs date with Recov_Combo by ptID
Recov_Combo_sub <- left_join(Recov_Combo_sub, Last_obs_Date, by="ptID")
rm(Last_obs_Date)
# Subset to only obs after initial recovery date 
Recov_Combo_sub$Obs_Date <- as.Date(Recov_Combo_sub$Obs_Date)
Recov_Combo_sub$Min_Date <- as.Date(Recov_Combo_sub$Min_Date)
Recov_Combo_sub <- Recov_Combo_sub[which(Recov_Combo_sub$Obs_Date >= Recov_Combo_sub$Min_Date),]
# new columns
Recov_Combo_sub$Rec_Date <- NA
Recov_Combo_sub$rec_obs <- NA
# New dataframe
Recov_Time_NDVI <- data.frame(matrix(nrow = 0, ncol = 66)) 
colnames(Recov_Time_NDVI) <- colnames(Recov_Combo_sub)

# STAY RECOVERED
for (i in unique(Recov_Combo_sub$ptID)) {
  print(i)
  # subset to i
  subID <- Recov_Combo_sub[which(Recov_Combo_sub$ptID == i),]
  # index all dates above threshold
  dates.above <- subID$Obs_Date[subID$NDVI >= subID$lwr80]
  # loop through possible recovery dates
  for (d in 1:length(dates.above)){
    # subset to obs between min_date[d] and last obs date
    subDate <- subID[which(subID$Obs_Date >= d & subID$Obs_Date <= subID$Last_obs_Date),]
    # identify total number of obs
    n.obs <- as.numeric(length(subDate$NDVI))
    # identify number of obs above threshold
    n.above <- as.numeric(length(subDate$NDVI[subDate$NDVI >= subDate$lwr80]))
    # if desired % of obs over threshold....
    if ((n.above / n.obs *100) >= 80 & (interval((subID$Obs_Date[d]), (unique(subID$Last_obs_Date))) %/% months(1)) >= 3){
      # record percentage of obs used to decide
      subID$rec_obs <- as.numeric(n.above / n.obs *100)
      # write out d as recovery date
      subID$Rec_Date <- subID$Obs_Date[d]
      RecRow <- subID[which(subID$Rec_Date == subID$Obs_Date),]
      Recov_Time_NDVI <- smartbind(RecRow, Recov_Time_NDVI)
      rm(RecRow)
      break
    } else {
      print("not recovered")
    }
  }
}


#70% threshold for 6 months (test 3)
# 70 %
Recov_Combo$lwr70 <- Recov_Combo$lwrNDVI *.7
threshold <- Recov_Combo[which(Recov_Combo$NDVI >= Recov_Combo$lwr70),]
# how many pts reach target?
length(unique(threshold$ptID)) # 109,539
# summarize by pt to get min date (first date with NDVI above threshold)
threshold$Obs_Date <- as.Date(threshold$Obs_Date)
Min_Date <- threshold %>% 
  group_by(ptID) %>%
  summarise(Min_Date = min(Obs_Date, na.rm = T)) 
# subset all observations to just those pts
Recov_Combo_sub <- subset(Recov_Combo, ptID %in% threshold$ptID)
length(unique(Recov_Combo_sub$ptID)) # 60,296
rm(threshold)
# join min_date with Recov_Combo_sub by ptID
Recov_Combo_sub <- left_join(Recov_Combo_sub, Min_Date, by="ptID")
summary(Recov_Combo_sub$Min_Date)
rm(Min_Date)
# identify date of last observation
Recov_Combo_sub$Obs_Date <- as.Date(Recov_Combo_sub$Obs_Date)
Last_obs_Date <- Recov_Combo_sub %>% 
  group_by(ptID) %>%
  summarise(Last_obs_Date = max(Obs_Date, na.rm = T))
# join last obs date with Recov_Combo by ptID
Recov_Combo_sub <- left_join(Recov_Combo_sub, Last_obs_Date, by="ptID")
rm(Last_obs_Date)
# Subset to only obs after initial recovery date 
Recov_Combo_sub$Obs_Date <- as.Date(Recov_Combo_sub$Obs_Date)
Recov_Combo_sub$Min_Date <- as.Date(Recov_Combo_sub$Min_Date)
Recov_Combo_sub <- Recov_Combo_sub[which(Recov_Combo_sub$Obs_Date >= Recov_Combo_sub$Min_Date),]
# new columns
Recov_Combo_sub$Rec_Date <- NA
Recov_Combo_sub$n_rec_obs <- NA
# New dataframe
Recov_Time_NDVI <- data.frame(matrix(nrow = 0, ncol = 66)) 
colnames(Recov_Time_NDVI) <- colnames(Recov_Combo_sub)


# STAY RECOVERED
for (i in unique(Recov_Combo_sub$ptID)) {
  print(i)
  # subset to i
  subID <- Recov_Combo_sub[which(Recov_Combo_sub$ptID == i),]
  # index all dates above threshold
  dates.above <- subID$Obs_Date[subID$NDVI >= subID$lwr70]
  # loop through possible recovery dates
  for (d in 1:length(dates.above)){
    # identify date 3 months post-min date
    Mo6_Date <- subID$Obs_Date[d] %m+% months(6)
    # remove all obs prior to min_date
    sub_Date <- subID[which(subID$Obs_Date >= subID$Obs_Date[d]),]
    # subset to obs within that period 
    sub6mo <- sub_Date[which(sub_Date$Obs_Date <= Mo6_Date),]
    # if all obs over threshold....
    if (length(sub6mo$NDVI[sub6mo$NDVI >= sub6mo$lwr70]) == length(sub6mo$NDVI)){
      # record how many obs were used to decide
      subID$n_rec_obs <- as.numeric(length(sub6mo$Obs_Date))
      # write out d as recovery date
      subID$Rec_Date <- subID$Obs_Date[d]
      RecRow <- subID[which(subID$Rec_Date == subID$Obs_Date),]
      Recov_Time_NDVI <- smartbind(RecRow, Recov_Time_NDVI)
      rm(RecRow)
      break
    } else {
      print("not recovered")
    }
  }
}


# How many pts do you get if you just have to get AND STAY within 80% of the baseline?...........................

# Using Recov_Combo with post-anticedent fire obs removed...
load("I:/Malone Lab/ENP Fire/Grace_McLeod/Recovery/Recov_Combo2.RDATA")

# MEASURE RECOVERY TIME TO 80% OF THRESHOLD
# calculate 80% of thrshold
Recov_Combo$lwr80 <- Recov_Combo$lwrNDVI * 0.8
# pull out all obs with NDVI > low80 
Recov_NDVI <- Recov_Combo[which(Recov_Combo$NDVI >= Recov_Combo$lwr80),] # recovered obs
# summarize by pt to get min date (first date with NDVI above threshold)
Rec_Date <- Recov_NDVI %>% 
  group_by(ptID) %>%
  summarise(Rec_Date = min(Obs_Date, na.rm = T)) 
# join min_date with Recov_Combo by ptID
length(unique(Recov_NDVI$ptID)) # 96,706
Recov_Combo <- left_join(Recov_Combo, Rec_Date, by="ptID")
rm(Rec_Date)

# CALCULATE RECOVERY TIME
# pull out distinct pt info
Recov_Time_NDVI <- Recov_Combo %>% 
  dplyr::select(ptID, Rec_Date, EndDate, FireYear, Severity) %>%
  distinct()
# only pts in Recov_NDVI (remove unrecovered pts)
Recov_Time_NDVI <- subset(Recov_Time_NDVI, ptID %in% Recov_NDVI$ptID)
length(unique(Recov_Time_NDVI$ptID)) 
# calculate recovery time
class(Recov_Time_NDVI$Rec_Date)
Recov_Time_NDVI$Rec_Date <- as.Date(Recov_Time_NDVI$Rec_Date)
Recov_Time_NDVI$EndDate <- as.Date(Recov_Time_NDVI$EndDate)
Recov_Time_NDVI$Rec_Days <- Recov_Time_NDVI$Rec_Date - Recov_Time_NDVI$EndDate
Recov_Time_NDVI$Rec_Yrs <- as.numeric(Recov_Time_NDVI$Rec_Days) /356
# add back other info from Recov_Combo
Recov_Time_NDVI$Obs_Date <- Recov_Time_NDVI$Rec_Date
names(Recov_Time_NDVI)
Recov_Time_NDVI <- merge(Recov_Time_NDVI, Recov_Combo, by=c("ptID","Rec_Date", "EndDate" ,"FireYear" ,"Severity", "Obs_Date"))
# save
setwd("I:/Malone Lab/ENP Fire/Grace_McLeod/Recovery")
save(Recov_Time_NDVI, file="Recov_Time_NDVI_80.RDATA")

# FILTER FOR PTS THAT STAY ABOVE 80% THRESHOLD
# only pts that do recover from first fire
Recov_Combo_80 <- subset(Recov_Combo, ptID %in% Recov_Time_NDVI$ptID)
length(unique(Recov_Combo_80$ptID))
# now you can do things like subset to where obs_date > rec_date for only obs after recovery date
Post_Rec_obs <- Recov_Combo_80[which(Recov_Combo_80$Obs_Date >= Recov_Combo_80$Rec_Date),]
# number of obs > lwr?
high_obs <- Post_Rec_obs[which(Post_Rec_obs$NDVI >= Post_Rec_obs$lwr80),]
high_obs_n <- high_obs %>% 
  group_by(ptID) %>%
  count()
colnames(high_obs_n)[which(names(high_obs_n) == "n")] <- "post_above_n"
Recov_Time_NDVI <- left_join(Recov_Time_NDVI, high_obs_n, by="ptID")
Recov_Time_NDVI$post_above_n[is.na(Recov_Time_NDVI$post_above_n)] <- 0
rm(high_obs_n)
# number of obs < lwr?
low_obs <- Post_Rec_obs[which(Post_Rec_obs$NDVI < Post_Rec_obs$lwr80),]
low_obs_n <- low_obs %>% 
  group_by(ptID) %>%
  count()
colnames(low_obs_n)[which(names(low_obs_n) == "n")] <- "post_below_n"
Recov_Time_NDVI <- left_join(Recov_Time_NDVI, low_obs_n, by="ptID")
Recov_Time_NDVI$post_below_n[is.na(Recov_Time_NDVI$post_below_n)] <- 0
rm(low_obs_n)
# total number of post-rec obs
Recov_Time_NDVI$post_total_n <- Recov_Time_NDVI$post_above_n + Recov_Time_NDVI$post_below_n
summary(Recov_Time_NDVI$post_total_n)
# proportion above threshold
Recov_Time_NDVI$prop.over.80 <- Recov_Time_NDVI$post_above_n / Recov_Time_NDVI$post_total_n * 100
summary(Recov_Time_NDVI$prop.over.80)
length(unique(Recov_Time_NDVI$ptID[Recov_Time_NDVI$prop.over.80 == 100]))



# OLD STUFF ............................................................................................................................

# only pts that do recover from first fire
Recov_Combo_REC <- subset(Recov_Combo, ptID %in% Recov_Time_NDVI$ptID)
length(unique(Recov_Combo_REC$ptID))
# now you can do things like subset to where obs_date > rec_date for only obs after recovery date
Post_Rec_obs <- Recov_Combo_REC[which(Recov_Combo_REC$Obs_Date >= Recov_Combo_REC$Rec_Date),]
# mean NDVI post rec_date? 
Post_stats <- Post_Rec_obs %>%
  group_by(ptID) %>%
  summarize(mean.post.ndvi=mean(NDVI, na.rm=T)) 
Recov_Time_NDVI<- left_join(Recov_Time_NDVI, Post_stats, by="ptID")
rm(Post_stats)
# number of obs > lwr?
high_obs <- Post_Rec_obs[which(Post_Rec_obs$NDVI >= Post_Rec_obs$lwrNDVI),]
high_obs_n <- high_obs %>% 
  group_by(ptID) %>%
  count()
colnames(high_obs_n)[which(names(high_obs_n) == "n")] <- "post_above_n"
Recov_Time_NDVI <- left_join(Recov_Time_NDVI, high_obs_n, by="ptID")
rm(high_obs_n)
# number of obs < lwr?
low_obs <- Post_Rec_obs[which(Post_Rec_obs$NDVI < Post_Rec_obs$lwrNDVI),]
low_obs_n <- low_obs %>% 
  group_by(ptID) %>%
  count()
colnames(low_obs_n)[which(names(low_obs_n) == "n")] <- "post_below_n"
Recov_Time_NDVI <- left_join(Recov_Time_NDVI, low_obs_n, by="ptID")
rm(low_obs_n)
# total number of post-rec obs
Recov_Time_NDVI$post_total_n <- Recov_Time_NDVI$post_above_n + Recov_Time_NDVI$post_below_n

# HOW TO DECIDE WHAT GOES?
summary(Recov_Time_NDVI$mean.post.ndvi)
summary(Recov_Time_NDVI$post_above_n)
summary(Recov_Time_NDVI$post_below_n) # 1,183  never drop back below min
# number of obs >70 %
Post_Rec_obs$lwr70 <- Post_Rec_obs$lwrNDVI * 0.7
Recov_Time_NDVI_70 <- Post_Rec_obs[which(Post_Rec_obs$NDVI >= Post_Rec_obs$lwr70),]
Rec_70 <- Recov_Time_NDVI_70 %>% 
  group_by(ptID) %>%
  count()
colnames(Rec_70)[which(names(Rec_70) == "n")] <- "lwr70_n"
Recov_Time_NDVI <- left_join(Recov_Time_NDVI, Rec_70, by="ptID")
summary(Recov_Time_NDVI$lwr70_n)
# number of obs >80 %
Post_Rec_obs$lwr80 <- Post_Rec_obs$lwrNDVI * 0.8
Recov_Time_NDVI_80 <- Post_Rec_obs[which(Post_Rec_obs$NDVI >= Post_Rec_obs$lwr80),]
Rec_80 <- Recov_Time_NDVI_80 %>% 
  group_by(ptID) %>%
  count()
colnames(Rec_80)[which(names(Rec_80) == "n")] <- "lwr80_n"
Recov_Time_NDVI <- left_join(Recov_Time_NDVI, Rec_80, by="ptID")
summary(Recov_Time_NDVI$lwr80_n)
# number of obs >90 %
Post_Rec_obs$lwr90 <- Post_Rec_obs$lwrNDVI * 0.9
Recov_Time_NDVI_90 <- Post_Rec_obs[which(Post_Rec_obs$NDVI >= Post_Rec_obs$lwr90),]
Rec_90 <- Recov_Time_NDVI_90 %>% 
  group_by(ptID) %>%
  count()
colnames(Rec_90)[which(names(Rec_90) == "n")] <- "lwr90_n"
Recov_Time_NDVI <- left_join(Recov_Time_NDVI, Rec_90, by="ptID")
summary(Recov_Time_NDVI$lwr90_n)
# number of obs >95 %
Post_Rec_obs$lwr95 <- Post_Rec_obs$lwrNDVI * 0.95
Recov_Time_NDVI_95 <- Post_Rec_obs[which(Post_Rec_obs$NDVI >= Post_Rec_obs$lwr95),]
Rec_95 <- Recov_Time_NDVI_95 %>% 
  group_by(ptID) %>%
  count()
colnames(Rec_95)[which(names(Rec_95) == "n")] <- "lwr95_n"
Recov_Time_NDVI <- left_join(Recov_Time_NDVI, Rec_95, by="ptID")
summary(Recov_Time_NDVI$lwr95_n)

# WHAT PROPORTION OF POST-FIRE OBS STAY ABOVE X% OF THRESHOLD?
Recov_Time_NDVI$prop.over.95 <- Recov_Time_NDVI$lwr95_n / Recov_Time_NDVI$post_total_n * 100
summary(Recov_Time_NDVI$prop.over.95) 
Recov_Time_NDVI$prop.over.90 <- Recov_Time_NDVI$lwr90_n / Recov_Time_NDVI$post_total_n * 100
summary(Recov_Time_NDVI$prop.over.90)
Recov_Time_NDVI$prop.over.80 <- Recov_Time_NDVI$lwr80_n / Recov_Time_NDVI$post_total_n * 100
summary(Recov_Time_NDVI$prop.over.80)
Recov_Time_NDVI$prop.over.70 <- Recov_Time_NDVI$lwr70_n / Recov_Time_NDVI$post_total_n * 100
summary(Recov_Time_NDVI$prop.over.70)
# plot
p95 <- ggplot(Recov_Time_NDVI, aes(y=prop.over.95)) + 
  geom_boxplot() + 
  labs(title="proportion of post-rec obs >95% of threshold")
p90 <- ggplot(Recov_Time_NDVI, aes(y=prop.over.90)) + 
  geom_boxplot() + 
  labs(title="proportion of post-rec obs >90% of threshold")
p80 <- ggplot(Recov_Time_NDVI, aes(y=prop.over.80)) + 
  geom_boxplot() + 
  labs(title="proportion of post-rec obs >80% of threshold")
p70 <- ggplot(Recov_Time_NDVI, aes(y=prop.over.70)) + 
  geom_boxplot() + 
  labs(title="proportion of post-rec obs >70% of threshold")
p95 + p90 + p80 + p70
# how many pts stay above each threshold 
length(unique(Recov_Time_NDVI$ptID[Recov_Time_NDVI$prop.over.95 == 100]))
length(unique(Recov_Time_NDVI$ptID[Recov_Time_NDVI$prop.over.90 == 100]))
length(unique(Recov_Time_NDVI$ptID[Recov_Time_NDVI$prop.over.80 == 100]))
length(unique(Recov_Time_NDVI$ptID[Recov_Time_NDVI$prop.over.70 >= 100]))
# look at recovery time for each thrshold
t_100 <- Recov_Time_NDVI[which(is.na(Recov_Time_NDVI$prop.over.95)),]
t_95 <- Recov_Time_NDVI[which(Recov_Time_NDVI$prop.over.95 == 100),]
t_90 <- Recov_Time_NDVI[which(Recov_Time_NDVI$prop.over.90 == 100),]
t_80 <- Recov_Time_NDVI[which(Recov_Time_NDVI$prop.over.80 == 100),]
t_70 <- Recov_Time_NDVI[which(Recov_Time_NDVI$prop.over.70 == 100),]
t100 <- ggplot(t_100, aes(y=Rec_Yrs)) + 
  geom_boxplot(col="orange") + 
  geom_density(col="purple") +   
  labs(title="100% Recovered")
t95 <- ggplot(t_95, aes(y=Rec_Yrs)) + 
  geom_boxplot(col="orange") + 
  geom_density(col="purple") +   
  labs(title="95% Recovered")
t90 <- ggplot(t_90, aes(y=Rec_Yrs)) + 
  geom_boxplot(col="orange") + 
  geom_density(col="purple") + 
  labs(title="90% Recovered")
t80 <- ggplot(t_80, aes(y=Rec_Yrs)) + 
  geom_boxplot(col="orange") + 
  geom_density(col="purple") + 
  labs(title="80% Recovered")
t70 <- ggplot(t_80, aes(y=Rec_Yrs)) + 
  geom_boxplot(col="orange") + 
  geom_density(col="purple") + 
  labs(title="70% Recovered")
(t100 + t95 + t90 ) / (t80 + t70)


# NOTES
# if "unrecovered" pts follow similar rec.time patterns as "recovered" within x % of threshold, 
# measure recovery time to x percent for all pts and see if its about the same. 
# if you lump together, is it about the same pattern?

# Loop through by ptID
for (i in unique(Recov_Combo_sub$ptID)) {
  print(i)
  # subset to i
  subID <- Recov_Combo_sub[which(Recov_Combo_sub$ptID == i),]
  # pull out last obs date
  Last_obs_Date <- unique(subID$Last_obs_Date)
  
  # List of all dates above threshold
  above.dates <- subID$Obs_Date[subID$NDVI >= subID$lwrNDVI]
  
  
  # loop through above dates
  for (d in 1:length(above.dates)) {
    print(above.dates[d])
    # reassign min date to be d
    subID$Min_Date <- above.dates[d]
    # identify date 3 months post-min date
    Mo3_Date <- d %m+% months(3)
    # subset to obs within that period (min-Date to 3mo date)
    sub3mo <- subID[which(subID$Obs_Date <= Mo3_Date),]
    # calc above and below
    # number of obs > lwr?
    high_obs <- sub3mo[which(sub3mo$NDVI >= sub3mo$lwrNDVI),]
    high_obs_n <- length(high_obs$ptID)
    sub3mo$post_above_n <- high_obs_n
    sub3mo$post_above_n[is.na(sub3mo$post_above_n)] <- 0
    rm(high_obs_n, high_obs)
    # number of obs < lwr?
    low_obs <- sub3mo[which(sub3mo$NDVI < sub3mo$lwrNDVI),]
    low_obs_n <- length(low_obs$ptID)
    sub3mo$post_below_n <- low_obs_n
    sub3mo$post_below_n[is.na(sub3mo$post_below_n)] <- 0
    rm(low_obs_n)
    # total number of post-rec obs
    sub3mo$post_total_n <- sub3mo$post_above_n + sub3mo$post_below_n
    # if all recovered, write out date as rec date
    if(min(sub3mo$post_total_n) == min(sub3mo$post_above_n)){
      print("all above")
      # set rec_date to min_date
      sub3mo$Rec_Date <- sub3mo$Min_Date
      # write out that date 
      RecRow <- sub3mo[which(sub3mo$Rec_Date == sub3mo$Obs_Date),]
      Recov_Time_NDVI <- smartbind(RecRow, Recov_Time_NDVI)
      rm(RecRow)
    } else {
      print("some below") 
    }
    if(Min_Date >= Last_obs_Date %m-% months(3)){
      print("insuf. time")
      # set rec_date to NA
      sub3mo$Rec_Date <- NA
      # write out to recovery file
      RecRow <- sub3mo[which(is.na(sub3mo$Rec_Date)),]
      Recov_Time_NDVI <- smartbind(RecRow, Recov_Time_NDVI)
      rm(RecRow)
      break
    } else {
      print("rerun")
    }
    
    
    
    
    # REASSIGN MIN DATE AND REAPEAT
    repeat {
      # List of all remaining dates above threshold
      above.dates <- subID$Obs_Date[subID$NDVI >= subID$lwrNDVI]
      # loop through above dates
      for (d in 1:length(above.dates)) {
        print(above.dates[d])
        # reassign min date to be d
        subID$Min_Date <- above.dates[d]
        # identify date 3 months post-min date
        subID$Mo3_Date <- subID$Min_Date %m+% months(3)
        # subset to obs within that period (min-Date to 3mo date)
        sub3mo <- subID[which(subID$Obs_Date <= subID$Mo3_Date),]
        # recalc above and below
        # number of obs > lwr?
        high_obs <- sub3mo[which(sub3mo$NDVI >= sub3mo$lwrNDVI),]
        high_obs_n <- length(high_obs$ptID)
        sub3mo$post_above_n <- high_obs_n
        sub3mo$post_above_n[is.na(sub3mo$post_above_n)] <- 0
        rm(high_obs_n, high_obs)
        # number of obs < lwr?
        low_obs <- sub3mo[which(sub3mo$NDVI < sub3mo$lwrNDVI),]
        low_obs_n <- length(low_obs$ptID)
        sub3mo$post_below_n <- low_obs_n
        sub3mo$post_below_n[is.na(sub3mo$post_below_n)] <- 0
        rm(low_obs_n)
        # total number of post-rec obs
        sub3mo$post_total_n <- sub3mo$post_above_n + sub3mo$post_below_n
        # if all recovered, write out date as rec date
        if(min(sub3mo$post_total_n) == min(sub3mo$post_above_n)){
          print("all above")
          # set rec_date to min_date
          sub3mo$Rec_Date <- sub3mo$Min_Date
          # write out that date 
          RecRow <- sub3mo[which(sub3mo$Rec_Date == sub3mo$Obs_Date),]
          Recov_Time_NDVI <- smartbind(RecRow, Recov_Time_NDVI)
          rm(RecRow)
          break
        } else {
          print("some below") 
        }
      }
    }
    
    
    # NEED TO FIND NEXT DATE (repeat until recovered)
    repeat {
      # remove row where min date = obs date
      subID <- subID[which(subID$Min_Date != subID$Obs_Date),]
      # List of all remaining dates above threshold
      above.dates <- subID$Obs_Date[subID$NDVI >= subID$lwrNDVI]
      # loop through above dates
      for (d in 1:length(above.dates)) {
        print(above.dates[d])
        # reassign min date to be d
        subID$Min_Date <- above.dates[d]
        # identify date 3 months post-min date
        subID$Mo3_Date <- subID$Min_Date %m+% months(3)
        # subset to obs within that period (min-Date to 3mo date)
        sub3mo <- subID[which(subID$Obs_Date <= subID$Mo3_Date),]
        # recalc above and below
        # number of obs > lwr?
        high_obs <- sub3mo[which(sub3mo$NDVI >= sub3mo$lwrNDVI),]
        high_obs_n <- length(high_obs$ptID)
        sub3mo$post_above_n <- high_obs_n
        rm(high_obs_n, high_obs)
        # number of obs < lwr?
        low_obs <- sub3mo[which(sub3mo$NDVI < sub3mo$lwrNDVI),]
        low_obs_n <- length(low_obs$ptID)
        sub3mo$post_below_n <- low_obs_n
        rm(low_obs_n)
        # total number of post-rec obs
        sub3mo$post_total_n <- sub3mo$post_above_n + sub3mo$post_below_n
        # if all recovered, write out date as rec date
        if(min(sub3mo$post_total_n) == min(sub3mo$post_above_n)){
          print("all above")
          # set rec_date to min_date
          sub3mo$Rec_Date <- sub3mo$Min_Date
          # write out that date 
          RecRow <- sub3mo[which(sub3mo$Rec_Date == sub3mo$Obs_Date),]
          Recov_Time_NDVI <- smartbind(RecRow, Recov_Time_NDVI)
          rm(RecRow)
          break
        } else {
          print("some below") 
        }
      }
      
    }
  }
  
  
  
  
 




# CALCUALTE RECOVERY TIME....................TO LOWER NDVI THRSHOLD.............................................................................................................

# pull out all obs with NDVI > lowNDVI (recovered only)
Recov_NDVI <- Recov_Combo[which(Recov_Combo$NDVI >= Recov_Combo$lwrNDVI),] 
# summarize by pt to get min date (first date with NDVI above threshold)
Rec_Date <- Recov_NDVI %>% 
  group_by(ptID) %>%
  summarise(Rec_Date = min(Obs_Date, na.rm = T)) 
# join min_date with Recov_Combo by ptID
length(unique(Recov_Combo_og$ptID)) # 114,294
length(unique(Recov_Combo$ptID)) # 114,058
length(unique(Recov_NDVI$ptID)) # 60,296 
Recov_Combo <- left_join(Recov_Combo, Rec_Date, by="ptID")
rm(Rec_Date)
head(Recov_Combo)

# WHY DID SOME PTS NEVER RECOVER??? - See "Unrecovered" script 
# pull out unrecovered pts with obs prior to antecedent fire
unrecovered <- anti_join(Recov_Combo, Recov_NDVI, by="ptID")
length(unique(unrecovered$ptID)) # 53,762 (47 %)
# save 
setwd("I:/Malone Lab/ENP Fire/Grace_McLeod/Recovery")
save(unrecovered, file = "Unrecovered.RDATA")


# CALCULATE RECOVERY TIME
# pull out distinct pt info
Recov_Time_NDVI <- Recov_Combo %>% 
  dplyr::select(ptID, Rec_Date, EndDate, FireYear, Severity) %>%
  distinct()
# only pts in Recov_NDVI
Recov_Time_NDVI <- subset(Recov_Time_NDVI, ptID %in% Recov_NDVI$ptID)
# calculate recovery time
Recov_Time_NDVI <- na.omit(Recov_Time_NDVI)
class(Recov_Time_NDVI$EndDate)
Recov_Time_NDVI$Rec_Date <- as.Date(Recov_Time_NDVI$Rec_Date)
Recov_Time_NDVI$EndDate <- as.Date(Recov_Time_NDVI$EndDate)
Recov_Time_NDVI$Rec_Days <- Recov_Time_NDVI$Rec_Date - Recov_Time_NDVI$EndDate
Recov_Time_NDVI$Rec_Yrs <- as.numeric(Recov_Time_NDVI$Rec_Days) /356
# add back other info from Recov_Combo
Recov_Time_NDVI$Obs_Date <- Recov_Time_NDVI$Rec_Date
names(Recov_Time_NDVI)
Recov_Time_NDVI <- merge(Recov_Time_NDVI, Recov_Combo, by=c("ptID","Rec_Date", "EndDate" ,"FireYear" ,"Severity", "Obs_Date"))
# save
setwd("I:/Malone Lab/ENP Fire/Grace_McLeod/Recovery")
save(Recov_Time_NDVI, file="Recov_Time_NDVI1.RDATA")









# CLEAN UP RECOVERY DATAFRAME
names(Recov_Time_NDVI)
# Re-assign original values to substituted columns
Recov_Time_NDVI$TotalFires <- Recov_Time_NDVI$TotalFires.og ; Recov_Time_NDVI$TotalFires.og <- NULL ;  Recov_Time_NDVI$TotalFires.sub <- NULL
Recov_Time_NDVI$Prev.Int <- Recov_Time_NDVI$Prev.Int.og ; Recov_Time_NDVI$Prev.Int.og <- NULL ; Recov_Time_NDVI$Prev.Int.sub <- NULL
Recov_Time_NDVI$Pt.B4min <- Recov_Time_NDVI$Pt.B4min.og ; Recov_Time_NDVI$Pt.B4min.og <- NULL ;  Recov_Time_NDVI$Pt.B4min.sub <- NULL 
Recov_Time_NDVI$Pt.B4max <- Recov_Time_NDVI$Pt.B4max.og ; Recov_Time_NDVI$Pt.B4max.og <- NULL ;  Recov_Time_NDVI$Pt.B4max.sub <- NULL 
Recov_Time_NDVI$Pt.B5max <- Recov_Time_NDVI$Pt.B5max.og ; Recov_Time_NDVI$Pt.B5max.og <- NULL ;  Recov_Time_NDVI$Pt.B5max.sub <- NULL 
Recov_Time_NDVI$Pt.B7max <- Recov_Time_NDVI$Pt.B7max.og ; Recov_Time_NDVI$Pt.B7max.og <- NULL ;  Recov_Time_NDVI$Pt.B7max.sub <- NULL 
Recov_Time_NDVI$SWIR1.SWIR2 <- Recov_Time_NDVI$SWIR1.SWIR2.og ; Recov_Time_NDVI$SWIR1.SWIR2.og <- NULL ; Recov_Time_NDVI$SWIR1.SWIR2.sub <- NULL
Recov_Time_NDVI$NIR.SWIR1 <- Recov_Time_NDVI$NIR.SWIR1.og ; Recov_Time_NDVI$NIR.SWIR1.og <- NULL ; Recov_Time_NDVI$NIR.SWIR1.sub <- NULL


# Percent recovered: 
61856/114294*100 # = 54 %

# save updates
setwd("I:/Malone Lab/ENP Fire/Grace_McLeod/Recovery")
save(Recov_Time_NDVI, file="Recov_Time_NDVI.RDATA")





# WHAT'S DEIVING VARIATION IN RECOVERY TIME?

# check it out
summary(Recov_Time_NDVI)
Recov_Time_NDVI$Rec_Days <- as.numeric(Recov_Time_NDVI$Rec_Days)
ggplot(Recov_Time_NDVI, aes(y=Rec_Yrs)) +
  geom_boxplot(col="orange") + 
  geom_density(col="purple") + 
  theme_bw()
# relationship between recovery time and immediate post-fire NDVI
library(mgcv)
library(voxel)
library(emmeans)
postNDVI <- gam(Rec_Yrs ~  s(PostNDVI),
                data = Recov_Time_NDVI)
plot(postNDVI, pages = 1, scheme=2)
plotGAM(postNDVI, smooth.cov="PostNDVI", plotCI=TRUE) + 
  ylim(0, 20) +
  ylab("Recovery Time") + 
  xlab("Immediate Post-Fire NDVI") +
  theme(text = element_text(size = 20))
ggplot(Recov_Time_NDVI, aes(x=PostNDVI, y=Rec_Yrs)) +
  geom_smooth()+
  #geom_point(size=3) +
  #geom_bar(position="dodge", stat="identity") +
  #geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.3, size=.75) +
  theme_bw() +
  labs(y="Recovery Time (years)", x="Immediate Post NDVI") +
  theme(text = element_text(size = 20), 
        legend.position = "bottom") 
# recovery categories 
Recov_Time_NDVI$RecCat <- NA
Recov_Time_NDVI$RecCat[Recov_Time_NDVI$Rec_Yrs <2] <- "<2 years"
Recov_Time_NDVI$RecCat[Recov_Time_NDVI$Rec_Yrs >=2 & Recov_Time_NDVI$Rec_Yrs <10] <- "2-10 years"
Recov_Time_NDVI$RecCat[Recov_Time_NDVI$Rec_Yrs >10] <- ">10 years"
Recov_Time_NDVI$RecCat <- as.factor(Recov_Time_NDVI$RecCat)
Recov_Time_NDVI$RecCat <- factor(Recov_Time_NDVI$RecCat, levels=c("<2 years","2-10 years", ">10 years") )
levels(Recov_Time_NDVI$RecCat)
postNDVI_cat <- gam(PostNDVI ~ RecCat,
                    data=Recov_Time_NDVI)
m1_emm <- emmeans(postNDVI_cat, specs = c("RecCat"))
df <- as.data.frame(m1_emm)
ggplot(df, aes(x=RecCat, y=emmean, color=RecCat)) +
  geom_point(size=5) +
  #geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.5, linewidth=.75) +
  theme_bw() +
  labs(y="Immediate Post-Fire NDVI", x=NULL, color="Recovery Category") +
  theme(text = element_text(size = 20), 
        legend.position = "bottom") +
  scale_color_manual(values=c("#c9af0a", "#3f5c93", "#6d973f"))
# severity
sev.gam <- gam(Rec_Yrs ~  s(Severity),
               data = Recov_Time_NDVI)
plotGAM(sev.gam, smooth.cov="Severity", plotCI=TRUE) + 
  ylim(0, 20) +
  ylab("Recovery Time") + 
  xlab("Severity") +
  theme(text = element_text(size = 20))










# REMOVE PTS THAT REBURN PRIOR TO RECOVERY 
# load full fire history dataframe
load("I:/Malone Lab/ENP Fire/Grace_McLeod/Fire_History/FireYears_df.RDATA")
head(FireYears_df)
# filter FH to recov_NDVI pts
Recov_FH <- subset(FireYears_df, ptID %in% Recov_Time_NDVI$ptID)
# only care about 2008-2020 bc we know they only burned once before that (recFire)
names(Recov_FH)
Recov_FH <- Recov_FH %>% dplyr::select(ptID,EVG_.2008,EVG_.2009,EVG_.2010,EVG_.2011,
                                       EVG_.2012,EVG_.2013,EVG_.2014,EVG_.2015,EVG_.2016,EVG_.2017,EVG_.2018,EVG_.2019,EVG_.2020, coords.x1, coords.x2)
# Pull out rec year
Recov_Time_NDVI$Obs_Year <- format(as.Date(Recov_Time_NDVI$Rec_Date, format="%Y-%m-%d"),"%Y")
Rec_Yrs <- Recov_Time_NDVI %>% dplyr::select(ptID, Obs_Year)
names(Rec_Yrs)
Recov_FH <- merge(Recov_FH, Rec_Yrs, by="ptID")
Recov_FH$Obs_Year <- as.numeric(Recov_FH$Obs_Year)
summary(Recov_FH$Obs_Year)
# Year of antecedent fire
Recov_FH$AFyear <- NA
names(Recov_FH)
Recov_FH$AFyear <- apply(Recov_FH[,2:14], 1, FUN=min,  na.rm = TRUE)
Recov_FH[sapply(Recov_FH, is.infinite)] <- NA
summary(Recov_FH$AFyear)
# is antecedent fire before recovery year?
Recov_FH$reburn <- "no"
Recov_FH$reburn[Recov_FH$Obs_Year > Recov_FH$AFyear] <- "yes"   #  8314 (8314/70170 =11% reburn)
Recov_FH$reburn <- as.factor(Recov_FH$reburn)
summary(Recov_FH$reburn)
# Filter out pts that reburn prior to recovery
NoReburn <- Recov_FH[which(Recov_FH$reburn == "no"),]
Recov_Time_NDVI <- subset(Recov_Time_NDVI, ptID %in% NoReburn$ptID)
# remove observations after anticedent fire
Recov_Time_NDVI <- Recov_Time_NDVI[which(Recov_Time_NDVI$Obs_Year < Recov_Time_NDVI$AFY),]




##########################################################################################################################################################
# FINAL STATS: Severty, Recovery Time, ect
##########################################################################################################################################################



# RECOVERY TIME..........................................................................................................................

# NDVI
test <- Recov_Time_NDVI[which(Recov_Time_NDVI$Rec_Yrs <= 0),]
summary(Recov_Time_NDVI$Rec_Yrs)
# Min.       1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.002809  0.176966   0.365169  0.568108  0.688202   19.000000 
summary(Recov_Time_NDVI$Rec_Yrs[Recov_Time_NDVI$EcoType=="Pineland"])
summary(Recov_Time_NDVI$Rec_Yrs[Recov_Time_NDVI$EcoType=="Hammock"])
# pts recovered
148623/149826*100 # = 99.2%

# NBR
summary(Recov_Time_NBR$Rec_Yrs)
#Min.      1st Qu.    Median      Mean    3rd Qu.      Max. 
#0.002809  0.202247  0.426966   0.624050  0.693820   21.365169
# pts recovered
149301/149826*100 # = 99.65 %

# plot
levels(Recov_Time_NDVI$EcoType)
Recov_Time_NDVI$EcoType <- factor(Recov_Time_NDVI$EcoType, levels=c("Pineland", "Hammock") )
ggplot(Recov_Time_NDVI, aes(x=EcoType, y=Rec_Yrs)) + 
  geom_violin() +
  #scale_color_manual(values=c("#660066", "#ff6633"))+ 
  #labs(title="NDVI Recovery Time") +
  labs(x= "Ecosystem Type", y="Recovery Time (years)" ) +
  theme_bw()+
  theme(legend.position = "none") +
  theme(axis.title = element_text(size = 15), axis.text = element_text(size=15))+
  geom_hline(yintercept=1, linetype="dashed", size=1, color = "#c9af0a") +
  geom_hline(yintercept=5, linetype="dashed", size=1, color = "#3f5c93") +
  geom_hline(yintercept=20, linetype="dashed", size=1, color = "#6d973f")
  
  

Recov_Time_NDVI.5 <- Recov_Time_NDVI[which(Recov_Time_NDVI$Rec_Yrs <= 5),]
ggplot(Recov_Time_NDVI.5, aes(x=EcoType, y=Rec_Yrs, color=EcoType)) +
  geom_jitter()+
  geom_violin() +
  scale_color_manual(values=c("#660066", "#ff6633"))+ 
  labs(title="NDVI Recovery Time") +
  labs(title="NDVI Recovery Time (< 5 years)") +
  labs(x= "Ecosystem Type", y="Recovery Time (years)" ) +
  theme_bw()+
  theme(legend.position = "none") +
  theme(axis.title = element_text(size = 20), axis.text = element_text(size=20))

ggplot(Recov_Time_NDVI, aes(x=Rec_Yrs, color=EcoType)) +
  geom_density()
ggplot(Recov_Time_NDVI.5, aes(x=Rec_Yrs, color=EcoType)) +
  geom_density()


# PERCENTAGE OF PTS FOR <1, 1-5, AND >5 YEARS REC TIME (NDVI)
pine <- Recov_Time_NDVI[which(Recov_Time_NDVI$EcoType == "Pineland"),]
hmk <- Recov_Time_NDVI[which(Recov_Time_NDVI$EcoType == "Hammock"),]

# < 1yr
length(pine$ptID[pine$Rec_Yrs < 1])
# 100103
100103/114803*100
# 87.19546 %
length(hmk$ptID[hmk$Rec_Yrs < 1])
# 30944
30944 / 33820 *100
# 91.49616 %

# 1-5 yrs
length(pine$ptID[pine$Rec_Yrs >= 1 & pine$Rec_Yrs <5])
# 13946
13946/114803*100
# 12.14777 %
length(hmk$ptID[hmk$Rec_Yrs >= 1 & hmk$Rec_Yrs < 5])
# 2668
2668 / 33820 *100
# 7.888823 %

# >5 yrs
length(pine$ptID[pine$Rec_Yrs > 5])
# 754
754/114803*100
# 0.6567773 %
length(hmk$ptID[hmk$Rec_Yrs > 5])
# 208
208/ 33820 *100
# 0.6150207 %





# SEVERITY ...................................................................................................................................

summary(Recov_Time_NDVI$Severity[Recov_Time_NDVI$EcoType =="Pineland"])
# Min.        1st Qu.     Median      Mean       3rd Qu.      Max. 
# 0.0000025   0.0436144  0.0864582   0.0986416   0.1413603   0.9589902 
summary(Recov_Time_NDVI$Severity[Recov_Time_NDVI$EcoType =="Hammock"])
# Min.       1st Qu.       Median      Mean      3rd Qu.       Max. 
# 0.0000068   0.0329970   0.0804740   0.1049274  0.1564542   0.4789629 

# plot
ggplot(Recov_Time_NDVI, aes(x=EcoType, y=Severity, color=EcoType)) +
  geom_boxplot() +
  scale_color_manual(values=c("#660066", "#ff6633"))+ 
  labs(x= "Ecosystem Type", y="Fire Severity" ) +
  theme_bw()+
  theme(legend.position = "none") +
  theme(axis.title = element_text(size = 20), axis.text = element_text(size=20))





