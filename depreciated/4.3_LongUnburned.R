
library(tidyverse)
library(randomForest)
library(ggplot2)
library(ggpubr)


setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")
setwd("./Baseline")
load( file="BL_train_test.RDATA")

setwd('/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline')
load(file="NDVI_rf.RDATA")

# NDVI directly after fire
train <- train %>% mutate( NDVI.Model = predict(NDVI_rf,  cur_data()))
test <- test %>% mutate( NDVI.Model = predict(NDVI_rf,  cur_data()))


# Based on the training data:
TF.2 <- train %>% filter( TotalFires < 2  )
TF.2$ptID %>% length/ train$ptID %>% length *100
long.unburned <- TF.2 %>% filter( TSF >= 15  ) 
long.unburned$ptID %>% length/ train$ptID %>% length *100

train$NDVI %>% hist
long.unburned$NDVI %>% hist
long.unburned$NDVI %>% median
  
# Based on the test data
TF.2 <- test %>% filter( TotalFires < 2  )
TF.2$ptID %>% length/ test$ptID %>% length *100
long.unburned <- TF.2 %>% filter( TSF > 15  ) 
long.unburned$ptID %>% length/ test$ptID %>% length *100

ideal <- train %>% filter( TotalFires == 9) %>% reframe(NIR.SWIR1 = mean(NIR.SWIR1), 
                                                        SWIR1.SWIR2 = mean(SWIR1.SWIR2))
# Set the target for the long unburned:
long.unburned.target <-long.unburned %>% filter( TotalFires < 2, TSF > 15) %>% mutate(
  TotalFires = 9,
  NIR.SWIR1 = ideal$NIR.SWIR1,
  SWIR1.SWIR2= ideal$SWIR1.SWIR2,
  NDVI.Model = predict(NDVI_rf,  cur_data()))


long.unburned
lm( NDVI.Model ~ NDVI, data =train) %>% summary
lm( NDVI.Model ~ NDVI, data =long.unburned) %>% summary
lm( NDVI.Model ~ NDVI, data =long.unburned.target) %>% summary

long.unburned %>% ggplot(aes( x= NDVI, y= NDVI.Model)) + 
  geom_point( aes( x= NDVI, y= NDVI.Model)) + 
  stat_smooth(method = "lm", se=FALSE, color="red", formula = y ~ x) + 
  #geom_point(alpha=0.1) +
  stat_cor(aes(label = paste(after_stat(rr.label), after_stat(..p.label..), sep = "~`,`~")), # adds R^2 and p-value
           r.accuracy = 0.01,
           p.accuracy = 0.001,
           label.x = 0.2, label.y = 0.8, size = 3) +
  stat_regline_equation(aes(label = ..eq.label..), # adds equation to linear regression
                        label.x = 0.2, label.y = 0.5, size =3) +
  geom_abline(intercept = 0, slope = 1, col = 'grey50',linetype="dashed")

long.unburned.target %>% ggplot(aes( x= NDVI, y= NDVI.Model)) + 
  geom_point( ) + 
  stat_smooth(method = "lm", se=FALSE, color="red", formula = y ~ x) + 
  #geom_point(alpha=0.1) +
  stat_cor(aes(label = paste(after_stat(rr.label), after_stat(..p.label..), sep = "~`,`~")), # adds R^2 and p-value
           r.accuracy = 0.01,
           p.accuracy = 0.001,
           label.x = 0.2, label.y = 0.8, size = 3) +
  stat_regline_equation(aes(label = ..eq.label..), # adds equation to linear regression
                        label.x = 0.2, label.y = 0.5, size =3) +
  geom_abline(intercept = 0, slope = 1, col = 'grey50',linetype="dashed")


rm(list=ls())
setwd('/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod')
# Load data
load("./Recovery/Recov_Rates.RDATA")


summary(Recov_Rate)

long.unburned.rec <- Recov_Rate %>% filter(Prev.Int > 15, 
                                                     TotalFires <2)

# How much of the recovery data is long unburned?
long.unburned.rec$ptID %>% length /Recov_Rate$ptID %>% length*100


long.unburned.rec.high <- long.unburned.rec %>% filter(PreNDVI >= model.NDVI)

# How much of the long unburned had a higher pre-NDVI than target?
long.unburned.rec.high$ptID %>% length/long.unburned.rec$ptID %>% length*100

# How much of the long unburned had a higher pre-NDVI than target?
long.unburned.rec.high$ptID %>% length/Recov_Rate$ptID %>% length *100

mean(long.unburned.rec.high$Rec100_Yrs, na.rm=T)
mean(long.unburned.rec.high$Rec80_Yrs, na.rm=T)
mean(long.unburned.rec.high$Rec70_Yrs, na.rm=T)

mean(Recov_Rate$Rec100_Yrs, na.rm=T)
mean(Recov_Rate$Rec80_Yrs, na.rm=T)
mean(Recov_Rate$Rec70_Yrs, na.rm=T)

# Resulted in recovery rates much faster than what is observed on the landscape. 
