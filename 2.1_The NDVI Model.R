# The NDVI MODEL 
# M.Grace McLeod (2023)

# This script merges
# 1. formats data for Random Forest ("BL_train_test.RDATA")
# 2. builds the Random Forest Model for predicting NDVI ("NDVI_rf.RDATA")
# 3. performs a sensitivity analysis to evaluate drivers of NDVI 

rm(list=ls())

library(sf)
library(rgdal)
library(sp)
library(dplyr)
library(tidyverse)
library(raster)
library(ggplot2)
library(readr)
library(spatialEco)
library(tidyr)
library(MASS) 
library(reshape2) 
library(reshape) 
library(terra)
library(lubridate)
library(randomForest)
library(datasets)
library(caret)
library(sqldf)
library(splitstackshape)
library(corrplot)
library(terra)
library(ggpubr)


# Set working directory
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")


##########################################################################################################################################################
# 1. FORMAT DATA FOR RANDOM FOREST
##########################################################################################################################################################
# load BL_master_df
load("./Baseline/BL_Master_df.RDATA")
length(unique(BL_Master_df$ptID)) # 14,981 pts

# filter to only obs with TSF >= 7 years
BL_Master_df$TSF <- BL_Master_df$Obs_Year - BL_Master_df$LFY
BL_Master_df <- BL_Master_df[which(BL_Master_df$TSF >= 7),] 
length(unique(BL_Master_df$ptID)) # 14,981 pts (interestingly, same number of pts just more obs. so the points that dont burn foever doesnt change) 
summary(BL_Master_df$TSF)
# Set Prev.Int with NAs to 80
BL_Master_df$Prev.Int [is.na(BL_Master_df$Prev.Int)] =  80
# Make Obs_month numeric
BL_Master_df$Obs_month <- as.numeric(BL_Master_df$Obs_month)


# TEST FOR VARIABILITY OF NDVI
# look at the pts with low NDVI, is there a pattern in their variance?
summary(BL_Master_df$NDVI) 
# Min.        1st Qu.    Median      Mean     3rd Qu.      Max. 
# -0.001093  0.202053  0.252070  0.250824   0.297974   0.502491  
# we want points that are not experiencing other disturbances or anything else weird. "stable pts"
BL_NDVI_var <- BL_Master_df %>% 
  group_by(ptID) %>%
  summarise(NDVIvar=var(NDVI, na.rm = T), NDVImax=max(NDVI, na.rm = T))
summary(BL_NDVI_var$NDVIvar)
#  Min.        1st Qu.       Median     Mean      3rd Qu.        Max. 
# 0.0002287  0.0011790.  0.0017590   0.0020267   0.0026056   0.0140301
summary(BL_NDVI_var$NDVImax)
# Min.     1st Qu.   Median    Mean   3rd Qu.    Max. 
# 0.1318  0.2813    0.3309   0.3300   0.3786   0.5025 
BL_Master_df <- merge(BL_Master_df, BL_NDVI_var, by="ptID")

# Check out NDVI distribution
lowNDVI <- BL_Master_df[which(BL_Master_df$NDVI < .2),]
length(unique(lowNDVI$ptID))  # 13,907 out of 14,981 pts. so almost everybody hits a low NDVI. Not necessarily a location issue then.
summary(lowNDVI$NDVIvar) 
# look at pts that never get high
neverHigh <-  BL_Master_df[which(BL_Master_df$NDVImax < .2),]
length(unique(neverHigh$ptID)) #  173

# filter for just the obs with NDVImax greater than .2
highMAX <- BL_Master_df[which(BL_Master_df$NDVImax >= .2),]
length(unique(highMAX$ptID)) # 14,808
# now filter that to remove all low obs
HIGH <- highMAX[which(highMAX$NDVI >=.2),] # 200,000 fewer obs (23% of obs under .2)
length(unique(HIGH$ptID)) # no loss of pts 

# Filter df by number of observations per point...........................................................................................
# for a generalizable model, want to use the BEST data so points that have observations from all times of year.
# start with perfection. How does that limit the variables? how closely does the distribution match the original dataset?
# Frequency Table: How many points of each EcoType burned how many times?
# Total  
HIGH$ptID <- as.factor(HIGH$ptID) 
ptObsTotal <- dplyr::count(HIGH,ptID)
summary(ptObsTotal$n)
# Min.    1st Qu.    Median   Mean    3rd Qu.     Max. 
# 1.00   34.00     50.00     52.38   68.00     150.00

# remove pts with less than top quartile of obs
ptObsTotal$n <- as.numeric(ptObsTotal$n)
# mean
ptObsTotalRM68 <- ptObsTotal[which(ptObsTotal$n >= 68),] 
BL_Total_sample <- HIGH[HIGH$ptID %in% ptObsTotalRM68$ptID, ]
# How many points are there?
length(unique(BL_Total_sample$ptID)) # 3,736
# How many observations are there?
length(BL_Total_sample$ptID) # 333,290

# save
setwd("./Baseline")
BL_Master_filtered <- BL_Total_sample
save(BL_Master_filtered, file="BL_Master_filtered.RDATA")


# SUBSET TESTING AND TRAINING DATA
# use stratified sampling so that data has same distribution (real, train, test)
train <- stratified(BL_Master_filtered, c("Obs_month", "tmin", "tmax", "precip", "TotalFires", "Prev.Int"), .6)
test <- anti_join(BL_Master_filtered, train)
# save datasets to use for all models
save(train, test, file="BL_train_test.RDATA")



##########################################################################################################################################################
# 2. RANDOM FOREST NDVI MODEL
##########################################################################################################################################################

# load train and test data
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")
load("./Baseline/BL_train_test.RDATA")

# NDVI ............................................................................................................................................................... 

# Corrilation plots
load("./Baseline/BL_Master_filtered.RDATA")
BL_Master_filtered$Obs_month <- as.numeric(BL_Master_filtered$Obs_month)
BL_Master_noNA <- na.omit(BL_Master_filtered)
# all considered variables
M <-cor(BL_Master_filtered[, c("Pt.B1max" ,"Pt.B2max" , "Pt.B3max", "Pt.B4max" , "Pt.B5max", "Pt.B7max" , 
                               "Pt.B1min", "Pt.B2min" ,"Pt.B3min", "Pt.B4min", "Pt.B5min", "Pt.B7min",
                               "NIR.SWIR1", "SWIR1.SWIR2",
                               "Obs_month", "tmin", "tmax", "precip", 
                               "TotalFires", "Prev.Int",
                               "coords.x1", "coords.x2")]) 
corrplot(M, method="circle") 
corrplot(M, method="color")
corrplot(M, method="number")


# Run Random Forest Model
system.time(NDVI_rf_5.1  <- randomForest(NDVI ~  precip + tmax + Obs_month + 
                                           Prev.Int + TotalFires  + 
                                           NIR.SWIR1 + SWIR1.SWIR2 + 
                                           Pt.B4max + Pt.B5max + 
                                           Pt.B4min,
                                         data= train,
                                         ntree= 500,
                                         mtry= 3,
                                         nodesize = 10,  
                                         importance=TRUE,
                                         verbose=TRUE,
                                         predicted=TRUE,
                                         keep.inbag=TRUE)) # VE=  71.21  time= 1:40hr
varImpPlot(NDVI_rf_5.1)
# apply to testing data
# prediction and confusion matricies
p1 <- predict(NDVI_rf_5.1, train)
summary(p1)
# TESTING DATA
test$NDVI_rf <- predict(NDVI_rf_5.1, test)
# see corrilation for prediction
summary(lm(test$NDVI ~ test$NDVI_rf))       # R2= 0.7287  

# save the model
#setwd("./Baseline")
NDVI_rf <- NDVI_rf_5.1
save(NDVI_rf, file="NDVI_rf.RDATA")


# TARGET CONDITIONS.....................................................................................
# calculate target conditions based on ideal fire history 

summary(BL_Master_df$TotalFires)
# pine = historical interval ~ 5yrs
pine <- BL_Master_df[which(BL_Master_df$EcoType == "Pineland"),]
mean(pine$NDVI[pine$TotalFires == 9]) # NDVI = 0.2503734
mean(pine$NIR.SWIR1[pine$TotalFires == 9]) # NIR:SWIR = 1.178277
mean(pine$SWIR1.SWIR2[pine$TotalFires == 9]) # SWIR1:SWIR2 = 1.278579
mean(pine$Pt.B4max[pine$TotalFires == 9]) # B4max = 18242.24
mean(pine$MFRI[pine$TotalFires == 9]) # MFRI = 3.626805
summary(pine$NDVI[pine$TotalFires == 9]) 


##########################################################################################################################################################
# 3. SENSITIVITY ANALYSIS 
##########################################################################################################################################################

# DATASET(S)
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")
setwd("./Baseline")
load(file="BL_train_test.RDATA")

load("./BL_train_test.RDATA")
load(file="NDVI_rf.RDATA")

train.tf1 <- train %>% filter( TotalFires == 1) 
train.tf1$TSF %>% hist

length(train.tf1$TSF[ train.tf1$TSF > 15])/ length(train.tf1$TSF)

# Get a list of variables: 

data_name = getCall(NDVI_rf)$data %>% eval

library(randomForest)

data_name$predicted <- predict(NDVI_rf,data_name )

sensitivity_var.df <- function( rf.model, var.name){
  data_name = getCall(rf.model)$data
  data_again = eval(data_name)
  data_again %>% summary
  
  var.df <- varImpPlot(rf.model) %>% as.data.frame %>% rownames
  
  summary.df <- data_again[, var.df] %>% summarise_all(mean) %>% select(!var.name) %>% mutate( summary='mean')
  
  VOI <- data.frame(VOI =seq( min(data_again[, var.name]) , max(data_again[, var.name]), 0.1))
  VOI[, var.name] <- VOI$VOI
  VOI <- VOI %>% select(  var.name)
  totalfires.summary <- summary.df %>% merge(VOI)
  totalfires.summary$NDVI.model <- predict(rf.model, totalfires.summary)
  
  return(totalfires.summary )
}

totalfires.summary <- sensitivity_var.df(rf.model = NDVI_rf, var.name =  'TotalFires')
Prev.Int.summary <- sensitivity_var.df(rf.model = NDVI_rf, var.name =  'Prev.Int')
SWIR1.SWIR2.summary <- sensitivity_var.df(rf.model = NDVI_rf, var.name =  'SWIR1.SWIR2')
Obs_month.summary <- sensitivity_var.df(rf.model = NDVI_rf, var.name =  'Obs_month')
Pt.B4max.summary <- sensitivity_var.df(rf.model = NDVI_rf, var.name =  'Pt.B4max')
Pt.B4min.summary <- sensitivity_var.df(rf.model = NDVI_rf, var.name =  'Pt.B4min')
NIR.SWIR1.summary <- sensitivity_var.df(rf.model = NDVI_rf, var.name =  'NIR.SWIR1')
Pt.B5max.summary <- sensitivity_var.df(rf.model = NDVI_rf, var.name =  'Pt.B5max')
tmax.summary <- sensitivity_var.df(rf.model = NDVI_rf, var.name =  'tmax')
precip.summary <- sensitivity_var.df(rf.model = NDVI_rf, var.name =  'precip')


# Save baseline summaries:

save(totalfires.summary ,Prev.Int.summary ,  SWIR1.SWIR2.summary,
     Obs_month.summary,Pt.B4max.summary, Pt.B4min.summary,
     NIR.SWIR1.summary, Pt.B5max.summary,tmax.summary,precip.summary, file = "/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline_Sensitivity.Rdata" )

load(file = "/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline_Sensitivity.Rdata" )

library(tidyverse)
library(ggplot2)
library(ggpubr)

Pt.B4min.summary <- Pt.B4min.summary %>% mutate(Pt.B4min.adj =  (Pt.B4min * 0.0000275) + (-0.2) )

# plot 

sensitivity_plot_NDVI <- function( df, VOI, label) {
  
  ggplot(data = df)  + 
    geom_smooth(aes( x= VOI, y = NDVI.model), col='black') + theme_bw() +
    ylab('NDVI') + xlab(label) +theme(text = element_text(size = 8))
}

p.1 <- sensitivity_plot_NDVI( df=totalfires.summary , VOI = totalfires.summary$TotalFires , label='Total Fires')
p.2 <- sensitivity_plot_NDVI( df=Prev.Int.summary , VOI = Prev.Int.summary$Prev.Int , label='Time Since Fire')
p.3 <- sensitivity_plot_NDVI( df=SWIR1.SWIR2.summary , VOI = SWIR1.SWIR2.summary$SWIR1.SWIR2 , label='SWIR1:SWIR2')
p.4 <- sensitivity_plot_NDVI( df=Obs_month.summary , VOI = Obs_month.summary$Obs_month , label='Month')
p.5 <- sensitivity_plot_NDVI( df=Pt.B4max.summary , VOI = Pt.B4max.summary$Pt.B4max , label='Maximum Band 4')
p.6 <- sensitivity_plot_NDVI( df=Pt.B4min.summary , VOI = Pt.B4min.summary$Pt.B4min.adj , label='Minimum Band 4')
p.7 <- sensitivity_plot_NDVI( df=Pt.B5max.summary , VOI = Pt.B5max.summary$Pt.B5max , label='Maximum Band 5')
p.8 <- sensitivity_plot_NDVI( df=NIR.SWIR1.summary , VOI = NIR.SWIR1.summary$NIR.SWIR1 , label='NIR.SWIR1')
p.9 <- sensitivity_plot_NDVI( df=tmax.summary , VOI = tmax.summary$tmax , label='Maximum Air Temperature')
p.10 <- sensitivity_plot_NDVI( df=precip.summary , VOI = precip.summary$precip , label='Precipitation')


# One to One plot

one2one <- ggplot(data = data_name)  + geom_point(aes( x= NDVI, y = predicted)) +
  geom_smooth(aes( x= NDVI, y = predicted), method= "lm",col="goldenrod") + theme_bw() +
  ylab('Predicted') + xlab('Observed') +theme(text = element_text(size = 8)) + 
  geom_abline(intercept = 0, slope = 1, col="red", linetype="dashed")

setwd('/Users/sm3466/YSE Dropbox/Sparkle Malone/Research/ENP_Post-Fire_Recovery/Figures')
png(filename="Baseline_Sensitivity.png",
    width = 2000, height=1400, res=400)
ggarrange( one2one, p.1, p.2, p.4,
           p.3, p.6, p.8, 
           p.9, p.10, 
           ncol=3, nrow=3, 
           labels =c("A", "B", "C", "D", "E",
                     "F", "G", "H", "I"),
           font.label=list(color="black",size=10))

dev.off()








