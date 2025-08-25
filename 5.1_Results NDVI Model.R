# Results for the NDVI Model 

# Results

# This script generates figures for the Everglades Post-Fire Recovery Manuscript
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(randomForest)
library(ggridges)
library(sf)
library(sp)
library(terra)
library(leaflet)
library(ggspatial)
library(ggmap)
library(MetBrewer)
library(grid)
library(ggpubr)
library(corrplot)
library(ggcorrplot)

rm(list=ls())

# 3.1. The NDVI model  ####
rm(list= ls())

load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/BL_train_test_08082025.RDATA")

load( file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/NDVI_rf_08082025.RDATA")

NDVI_rf

train$pred <- predict(NDVI_rf, newdata=train)
test$pred <- predict(NDVI_rf, newdata=test)

lm( test$pred ~ test$NDVI) %>% summary

load(file='/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/Spatial_AutoCorrelation.RDATA'  )

library(spdep)
moran.results

library(gstats)


png(filename="/Users/sm3466/YSE Dropbox/Sparkle Malone/Research/ENP_Post-Fire_Recovery/Figures/Baseline_Variogram.png",
    width = 500, height=500, res=400)

dev.off()

# Summary statistics for the results section
test$NDVI %>% mean
test$pred %>% mean
test$pred %>% sd/sqrt(length(test$pred))

test$NDVI[test$TotalFires ==9] %>% mean
train$NDVI[train$TotalFires ==9] %>% mean

train$NDVI[train$EcoType == 'Pineland'] %>% mean(na.rm=T)
test$NDVI[train$EcoType == 'Pineland'] %>% mean(na.rm=T)
  
train$NDVI[train$EcoType == 'Hammock'] %>% mean

train$pred[train$EcoType == 'Pineland'] %>% mean
test$pred[train$EcoType == 'Pineland'] %>% mean(na.rm=T)


# Fire dependent

train.pineland$NDVI[ train.pineland$TSF >= 15 & train.pineland$TotalFires  < 2] %>% mean # Long unburned
train.pineland$NDVI[ train.pineland$TotalFires  >=7 ] %>% mean # 


length(train.pineland$NDVI[ train.pineland$TSF >= 15 & train.pineland$TotalFires  < 2])/ length( train.pineland$NDVI)
length(test.pineland$NDVI[ test.pineland$TSF >= 15 & test.pineland$TotalFires  < 2])/ length( test.pineland$NDVI)


# Figure 3: The NDVI Model Sensitivity Analysis: 

# BASELINE SENSITIVITY PLOT:  ####

load(file = "/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline_Sensitivity.Rdata" )

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
p.6 <- sensitivity_plot_NDVI( df=Pt.B4min.summary , VOI = Pt.B4min.summary$Pt.B4min , label='Minimum Band 4')
p.7 <- sensitivity_plot_NDVI( df=Pt.B5max.summary , VOI = Pt.B5max.summary$Pt.B5max , label='Maximum Band 5')
p.8 <- sensitivity_plot_NDVI( df=NIR.SWIR1.summary , VOI = NIR.SWIR1.summary$NIR.SWIR1 , label='NIR.SWIR1')
p.9 <- sensitivity_plot_NDVI( df=tmax.summary , VOI = tmax.summary$tmax , label='Air Temperature (°C)')
p.10 <- sensitivity_plot_NDVI( df=precip.summary , VOI = precip.summary$precip , label='Precipitation (mm)')


# One to One plot
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/BL_train_test_08082025.RDATA")
load( file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/NDVI_rf_08082025.RDATA")

data_name = getCall(NDVI_rf)$data %>% eval
data_name$predicted <- predict(NDVI_rf,data_name )

one2one <- ggplot(data = data_name)  + geom_point(aes( x= NDVI, y = predicted)) +
  geom_smooth(aes( x= NDVI, y = predicted), method= "lm",col="grey") + theme_bw() +
  ylab('Predicted') + xlab('Observed') +theme(text = element_text(size = 8)) + 
  geom_abline(intercept = 0, slope = 1, col="red", linetype="dashed")

png(filename="/Users/sm3466/YSE Dropbox/Sparkle Malone/Research/ENP_Post-Fire_Recovery/Figures/Baseline_Sensitivity.png",
    width = 2000, height=1400, res=400)
ggarrange( one2one, p.1, p.2, p.4,
           p.3, p.6, p.8, 
           p.9, p.10, 
           ncol=3, nrow=3, 
           labels =c("A", "B", "C", "D", "E",
                     "F", "G", "H", "I"),
           font.label=list(color="black",size=10))

dev.off()

# Observed versus predicted supplemental Figure:
load(file = "/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Figure3_Data_082025.RDATA")

test$pred <- predict(NDVI_rf, newdata= test)
summary.lm.val <- summary(lm(data=test, pred~ NDVI))
summary.lm.val$coefficients
summary.lm.val$r.squared

plot_BL_obsVSpred <- ggplot( data= test, aes(x= NDVI, y=pred)) + 
  geom_point() +
  geom_smooth(method='lm', se=TRUE, col='grey') + 
  labs(x="Observed", y="Predicted")+ 
  ylim(0.2, 0.6) + xlim(0.2, 0.5) + 
  geom_text( x=0.3, y=0.5, label="Y = 0.76*Observed + 0.068", size=3)+ geom_abline(intercept = 0, slope = 1, col="red", linetype="dashed") + 
  theme(text = element_text(size=20)) +theme_bw()


# Long-unburden:
total.baseline <- rbind( test, train)
long_unburned <- total.baseline %>% filter( TotalFires < 2 , Prev.Int > 15)
long_unburned$NDVI %>% mean(na.rm=T)
(long_unburned$ptID %>% unique %>% length / total.baseline$ptID %>% unique %>% length)*100

