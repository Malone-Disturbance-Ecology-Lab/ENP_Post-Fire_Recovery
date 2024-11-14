# DRIVERS ANALYSIS
# SPARKLE MALONE, GRACE McLEOD, 2024

# This script evaluates the influence of fire history and climate 
# on recovery threshold and recovery time. 

rm(list=ls())
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")

library(tidyverse)
library(corrplot)
library(mgcv)
library(voxel)
library(randomForest)
library(splitstackshape)
library(caret)
library(dplyr)
library(emmeans)
library(MetBrewer)
library(ggplot2)
library(caret)

# Load data
load("./Climate/Rec_pdsi_summary.RDATA")

#############################################################################################################################
# DATA PREP 
#############################################################################################################################

# data prep
names(pdsi.total.summary.FH.maxThreshold )
pdsi.total.summary.FH.maxThreshold$threshold
# normalize nwet and ndry by ntotal
pdsi.total.summary.FH.maxThreshold$nwet.frac = pdsi.total.summary.FH.maxThreshold$nwet/pdsi.total.summary.FH.maxThreshold$ntot
pdsi.total.summary.FH.maxThreshold$ndry.frac = pdsi.total.summary.FH.maxThreshold$ndry/pdsi.total.summary.FH.maxThreshold$ntot
# make Prev.Int more reasonable
pdsi.total.summary.FH$Prev.Int[pdsi.total.summary.FH$Prev.Int == 80] <- 45
pdsi.total.summary.FH.maxThreshold$Prev.Int[pdsi.total.summary.FH.maxThreshold$Prev.Int == 80] <- 45
# factorize Prev.Int (within or outside of historical interval ~ 10 yrs)
pdsi.total.summary.FH.maxThreshold$hist.Int[pdsi.total.summary.FH.maxThreshold$Prev.Int <= 10] <- 0
pdsi.total.summary.FH.maxThreshold$hist.Int[pdsi.total.summary.FH.maxThreshold$Prev.Int > 10] <- 1
pdsi.total.summary.FH.maxThreshold$hist.Int <- as.factor(pdsi.total.summary.FH.maxThreshold$hist.Int)
# factorize post date difference by 6mo intervals
pdsi.total.summary.FH.maxThreshold$PDD.6mo[pdsi.total.summary.FH.maxThreshold$PostDateDif <= 180] <- 0.5 
pdsi.total.summary.FH.maxThreshold$PDD.6mo[pdsi.total.summary.FH.maxThreshold$PostDateDif > 180 & pdsi.total.summary.FH.maxThreshold$PostDateDif <= 360] <- 2 
pdsi.total.summary.FH.maxThreshold$PDD.6mo[pdsi.total.summary.FH.maxThreshold$PostDateDif > 360 & pdsi.total.summary.FH.maxThreshold$PostDateDif <= 540] <- 1.5 
pdsi.total.summary.FH.maxThreshold$PDD.6mo[pdsi.total.summary.FH.maxThreshold$PostDateDif > 540] <- 2 
pdsi.total.summary.FH.maxThreshold$PDD.6mo <- as.factor(pdsi.total.summary.FH.maxThreshold$PDD.6mo)
pdsi.total.summary.FH.maxThreshold$FireYear.fac <- as.factor(pdsi.total.summary.FH.maxThreshold$FireYear)
# factorize observation month for pdsi
pdsi.total.summary.FH.maxThreshold$Obs_Mo <- format(pdsi.total.summary.FH.maxThreshold$Obs_Date,"%m")
pdsi.total.summary.FH.maxThreshold$Obs_Mo <- as.factor(pdsi.total.summary.FH.maxThreshold$Obs_Mo)
# bin and factorize preFire NDVI
load("./Recovery/Recov_Rates.RDATA")
PreNDVI <- Recov_Rate %>%
  select(ptID, PreNDVI)
pdsi.total.summary.FH.maxThreshold <- merge(pdsi.total.summary.FH.maxThreshold, PreNDVI, by="ptID")
summary(pdsi.total.summary.FH.maxThreshold$PreNDVI)
ggplot(pdsi.total.summary.FH.maxThreshold, aes(x=PreNDVI)) +
  geom_density()
pdsi.total.summary.FH.maxThreshold$PreNDVI.cat[pdsi.total.summary.FH.maxThreshold$PreNDVI <= 0.2] <- 1
pdsi.total.summary.FH.maxThreshold$PreNDVI.cat[pdsi.total.summary.FH.maxThreshold$PreNDVI > 0.2 & pdsi.total.summary.FH.maxThreshold$PreNDVI < 0.3] <- 2
pdsi.total.summary.FH.maxThreshold$PreNDVI.cat[pdsi.total.summary.FH.maxThreshold$PreNDVI > 0.3] <- 3
pdsi.total.summary.FH.maxThreshold$PreNDVI.cat <- as.factor(pdsi.total.summary.FH.maxThreshold$PreNDVI.cat)
# bin and factorize model_NDVI 
ggplot(pdsi.total.summary.FH.maxThreshold, aes(x=model.NDVI)) +
  geom_density()
pdsi.total.summary.FH.maxThreshold$modelNDVI.cat[pdsi.total.summary.FH.maxThreshold$model.NDVI <= 0.275] <- 1
pdsi.total.summary.FH.maxThreshold$modelNDVI.cat[pdsi.total.summary.FH.maxThreshold$model.NDVI > 0.275 & pdsi.total.summary.FH.maxThreshold$PreNDVI < 0.312] <- 2
pdsi.total.summary.FH.maxThreshold$modelNDVI.cat[pdsi.total.summary.FH.maxThreshold$model.NDVI > 0.312] <- 3
pdsi.total.summary.FH.maxThreshold$modelNDVI.cat <- as.factor(pdsi.total.summary.FH.maxThreshold$modelNDVI.cat)

# save edits
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Manuscript")
save(pdsi.total.summary.FH.maxThreshold, pdsi.total.summary.FH, pdsi.total.summary.FH, file="DriversData.RDATA")

# Correlation Matrix
names(pdsi.total.summary.FH.maxThreshold)
pdsi.total.summary.FH.maxThreshold$RecoveryTime.days <- pdsi.total.summary.FH.maxThreshold$RecoveryTime.days %>% as.numeric()
sum.MT <- pdsi.total.summary.FH.maxThreshold[, c(2:8,10, 13, 16:18,21,23)] %>% na.omit()
summary( sum.MT)
M <- cor(sum.MT)
corrplot(M, method="circle")

#############################################################################################################################
# Y = MaxThreshold 
#############################################################################################################################

# Data Manipulation
# Add indecies for fire frequency and pdsi conditions.............................................................................
# pdsi
summary(pdsi.total.summary.FH.maxThreshold$pdsi.mean)
pdsi.total.summary.FH.maxThreshold$pdsi.index <- "normal"
pdsi.total.summary.FH.maxThreshold$pdsi.index[pdsi.total.summary.FH.maxThreshold$pdsi.mean < -2] <- "dry"
pdsi.total.summary.FH.maxThreshold$pdsi.index[pdsi.total.summary.FH.maxThreshold$pdsi.mean > 2] <- "wet"
# fire frequency
# baased on historical MFRI over 29 years
# Total Fires: min= <4, med=4-6, high= >6
pdsi.total.summary.FH.maxThreshold$freq.index[pdsi.total.summary.FH.maxThreshold$TotalFires < 4] <- "<4"
pdsi.total.summary.FH.maxThreshold$freq.index[pdsi.total.summary.FH.maxThreshold$TotalFires >= 4 & pdsi.total.summary.FH.maxThreshold$TotalFires <=6] <- "4-6"
pdsi.total.summary.FH.maxThreshold$freq.index[pdsi.total.summary.FH.maxThreshold$TotalFires > 6 ] <- ">6"


# RANDOM FOREST......................................................................................................... 

# the point here is about suggesting improvement fo yo self. 
# Reaching a high threshold means you improved.
pdsi.total.summary.FH.maxThreshold[pdsi.total.summary.FH.maxThreshold == -Inf] <- NA
pdsi.total.summary.FH.maxThreshold.noNA <- pdsi.total.summary.FH.maxThreshold %>%
  na.omit 
pdsi.total.summary.FH.maxThreshold.noNA$Threshold <- as.factor(pdsi.total.summary.FH.maxThreshold.noNA$Threshold)
pdsi.total.summary.FH.maxThreshold.noNA$pdsi.index <- as.factor(pdsi.total.summary.FH.maxThreshold.noNA$pdsi.index)
pdsi.total.summary.FH.maxThreshold.noNA$freq.index <- as.factor(pdsi.total.summary.FH.maxThreshold.noNA$freq.index)


# Train and Test data
train <- stratified(pdsi.total.summary.FH.maxThreshold.noNA, c("Obs_Mo", "FireYear", "Severity", "PostDateDif", 
                                                               "pdsi.index",
                                                               "freq.index",
                                                               "model.NDVI", "PreNDVI"), .6)
test <- anti_join(pdsi.total.summary.FH.maxThreshold.noNA, train)

# Indecies RF
system.time(threshold_rf_index <- randomForest( Threshold ~  
                                                   PreNDVI +
                                                  freq.index +
                                                  pdsi.index,
                                                data= train,
                                                importance=TRUE,
                                                predicted=TRUE,
                                                keep.inbag=TRUE)) 
varImpPlot(threshold_rf_index)
p1 <- predict(threshold_rf_index, train)
confusionMatrix(p1, train$Threshold)
# try to just predict above or below 80
# include confusion matrix results in manuscript
train$rec.status <- "<80%"
train$rec.status[train$Threshold == "90" | train$Threshold == "100"] <- ">80%"
train$rec.status <- as.factor(train$rec.status)
system.time(t80_rf_index <- randomForest( rec.status ~  
                                                   PreNDVI +
                                                   freq.index +
                                                   pdsi.index,
                                                 data= train,
                                                 importance=TRUE,
                                                 predicted=TRUE,
                                                 keep.inbag=TRUE)) 
varImpPlot(t80_rf_index)
p1 <- predict(t80_rf_index, train)
confusionMatrix(p1, train$rec.status)

# save model 
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Manuscript")
save(threshold_rf_index, t80_rf_index, train, test, pdsi.total.summary.FH.maxThreshold.noNA, file="RF_threshold_index.RDATA")


# Sensitivity Analysis.............................................................................................................

# calculate mean for each index and adjust severity and PreNDVI values
mean(pdsi.total.summary.FH.maxThreshold$Severity[pdsi.total.summary.FH.maxThreshold$freq.index == "<4"], na.rm=T) # 0.1044326
mean(pdsi.total.summary.FH.maxThreshold$Severity[pdsi.total.summary.FH.maxThreshold$freq.index == "4-6"], na.rm=T) #  0.1007905
mean(pdsi.total.summary.FH.maxThreshold$Severity[pdsi.total.summary.FH.maxThreshold$freq.index == ">6"], na.rm=T) # 0.1066974

# make a synthetic dataframe 
sensitivity <- rbind(data.frame(Severity= 0.1044326, PreNDVI= seq(0.04, 0.4, 0.01), freq.index="<4", pdsi.index = "dry"), 
                                data.frame(Severity= 0.1007905, PreNDVI= seq(0.04, 0.4, 0.01), freq.index="4-6", pdsi.index = "dry"),
                                data.frame(Severity= 0.1066974, PreNDVI= seq(0.04, 0.4, 0.01), freq.index=">6", pdsi.index = "dry"),
                     data.frame(Severity= 0.1044326, PreNDVI= seq(0.04, 0.4, 0.01), freq.index="<4", pdsi.index = "normal"), 
                     data.frame(Severity= 0.1007905, PreNDVI= seq(0.04, 0.4, 0.01), freq.index="4-6", pdsi.index = "normal"),
                     data.frame(Severity= 0.1066974, PreNDVI= seq(0.04, 0.4, 0.01), freq.index=">6", pdsi.index = "normal"),
                     data.frame(Severity= 0.1044326, PreNDVI= seq(0.04, 0.4, 0.01), freq.index="<4", pdsi.index = "wet"), 
                     data.frame(Severity= 0.1007905, PreNDVI= seq(0.04, 0.4, 0.01), freq.index="4-6", pdsi.index = "wet"),
                     data.frame(Severity= 0.1066974, PreNDVI= seq(0.04, 0.4, 0.01), freq.index=">6", pdsi.index = "wet")) 
sensitivity$PostDateDif <- 600
sensitivity$freq.index <- as.factor(sensitivity$freq.index)
sensitivity$pdsi.index <- as.factor(sensitivity$pdsi.index)
# project model into dataframe
sensitivity$threshold.predict <- predict(threshold_rf_index, sensitivity)
sensitivity$threshold.predict80 <- predict(t80_rf_index, sensitivity)
# save
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Manuscript")
save(sensitivity, file="Sensitivity_data.RDATA")

# plot
sensitivity$freq.index <- factor(sensitivity$freq.index, levels=  c("<4", "4-6", ">6"))
ggplot(sensitivity) +
  geom_line(aes(y=PreNDVI, x=threshold.predict, col=pdsi.index), position= position_dodge(0.8), alpha= 0.5 , size=5) +
  facet_wrap(~freq.index) +
  theme_bw()

# compare training to testing
test$threshold.predict <- predict(threshold_rf_index, test)
# see corrilation for prediction
summary(lm(test$NDVI ~ test$threshold.predict))       # R2= 0.6447  






# TRASH........................................
# find high, med, and lows based on mean
summary(pdsi.total.summary.FH.maxThreshold$pdsi.mean) 
mean(pdsi.total.summary.FH.maxThreshold$pdsi.min[pdsi.total.summary.FH.maxThreshold$pdsi.mean < -3], na.rm=T) # -3.5 min
mean(pdsi.total.summary.FH.maxThreshold$pdsi.max[pdsi.total.summary.FH.maxThreshold$pdsi.mean < -3], na.rm=T) # -2.8 max
mean(pdsi.total.summary.FH.maxThreshold$pdsi.min[pdsi.total.summary.FH.maxThreshold$pdsi.mean < -0.4 & pdsi.total.summary.FH.maxThreshold$pdsi.mean > -0.5], na.rm=T) # -2.9 min
mean(pdsi.total.summary.FH.maxThreshold$pdsi.max[pdsi.total.summary.FH.maxThreshold$pdsi.mean < -0.4 & pdsi.total.summary.FH.maxThreshold$pdsi.mean > -0.5], na.rm=T) # 2.7 max
mean(pdsi.total.summary.FH.maxThreshold$pdsi.min[pdsi.total.summary.FH.maxThreshold$pdsi.mean > 2.5], na.rm=T) # 1.9 min
mean(pdsi.total.summary.FH.maxThreshold$pdsi.max[pdsi.total.summary.FH.maxThreshold$pdsi.mean > 2.5], na.rm=T) # 2.9 max
# make severity categories
summary(pdsi.total.summary.FH.maxThreshold$Severity)
ggplot(pdsi.total.summary.FH.maxThreshold) +
  geom_density(aes(x=Severity))
sev.df <- data.frame(Severity = c(0.02, 0.1, 0.96))
# make a synthetic dataframe
summary(pdsi.total.summary.FH.maxThreshold$PreNDVI)

sensitivity_preNDVI.pdsi <- rbind(data.frame(pdsi.mean = -3, pdsi.min= -3.5, pdsi.max= -2.8, PreNDVI= seq(0.04, 0.4, 0.01), index="dry"), 
                             data.frame(pdsi.mean = 0.45, pdsi.min= -2.9, pdsi.max= 2.7, PreNDVI= seq(0.04, 0.4, 0.01), index="normal"),
                             data.frame(pdsi.mean = 2.75, pdsi.min= 1.9, pdsi.max= 2.9, PreNDVI= seq(0.04, 0.4, 0.01), index="wet"))

sensitivity_preNDVI.freq <- rbind(data.frame(pdsi.mean = -3, pdsi.min= -3.5, pdsi.max= -2.8, PreNDVI= seq(0.04, 0.4, 0.01), index="<4"), 
                             data.frame(pdsi.mean = 0.45, pdsi.min= -2.9, pdsi.max= 2.7, PreNDVI= seq(0.04, 0.4, 0.01), index="4-6"),
                             data.frame(pdsi.mean = 2.75, pdsi.min= 1.9, pdsi.max= 2.9, PreNDVI= seq(0.04, 0.4, 0.01), index=">6"))

sensitivity_preNDVI <- merge(sensitivity_preNDVI, sev.df)

# project model into dataframe
sensitivity_preNDVI$threshold <- predict(threshold_rf_index, sensitivity_preNDVI)
# plot
ggplot(sensitivity_preNDVI) +
  geom_line(aes(x=PreNDVI, y=threshold, col=Severity), size=5) +
  facet_wrap(~Severity)
# you can make it to 100 under any conditions. 
# Any preNDVI can make it to 100 when wet. 
# for normal conditions, only really hight preNDVI make it to 100, and low make it to 70 or 80, but only high make it to 100
# for dry, only high ndvi make it to 100, low pre ndvi make it to 70 or 80
# so a low prendvi can be inhibiting to systems when there is not enough moisture. 
# ....
# If you make it to 100 under med conditions you had a high preNDVI. 
# getting to 100 under high conditions is associated with a low preNDVI


# Fire History 
summary(pdsi.total.summary.FH.maxThreshold$TotalFires) 
# Total Fires: min= 1-2, med=3-4, high=6-8, max= 10-11

# Prev.Int
mean(pdsi.total.summary.FH.maxThreshold$Prev.Int[pdsi.total.summary.FH.maxThreshold$TotalFires < 2], na.rm=T) # min 45
mean(pdsi.total.summary.FH.maxThreshold$Prev.Int[pdsi.total.summary.FH.maxThreshold$TotalFires >=3 & pdsi.total.summary.FH.maxThreshold$TotalFires <= 4], na.rm=T) # med =9
mean(pdsi.total.summary.FH.maxThreshold$Prev.Int[pdsi.total.summary.FH.maxThreshold$TotalFires >=6 & pdsi.total.summary.FH.maxThreshold$TotalFires <= 8], na.rm=T) # high =6
mean(pdsi.total.summary.FH.maxThreshold$Prev.Int[pdsi.total.summary.FH.maxThreshold$TotalFires >=10 & pdsi.total.summary.FH.maxThreshold$TotalFires <= 11], na.rm=T) # max=7

# Severity
mean(pdsi.total.summary.FH.maxThreshold$Severity[pdsi.total.summary.FH.maxThreshold$TotalFires < 2], na.rm=T) # min 0.1
mean(pdsi.total.summary.FH.maxThreshold$Severity[pdsi.total.summary.FH.maxThreshold$TotalFires >=3 & pdsi.total.summary.FH.maxThreshold$TotalFires <= 4], na.rm=T) # med = 0.09 
mean(pdsi.total.summary.FH.maxThreshold$Severity[pdsi.total.summary.FH.maxThreshold$TotalFires >=6 & pdsi.total.summary.FH.maxThreshold$TotalFires <= 8], na.rm=T) # high = 0.1
mean(pdsi.total.summary.FH.maxThreshold$Severity[pdsi.total.summary.FH.maxThreshold$TotalFires >=10 & pdsi.total.summary.FH.maxThreshold$TotalFires <= 11], na.rm=T) # max= 0.07

# make a synthetic dataframe 
summary(pdsi.total.summary.FH.maxThreshold$PreNDVI)
sensitivity_preNDVI_FH <- rbind(data.frame(TotalFires = 1.5, Prev.Int = 45, Severity= 0.1, PreNDVI= seq(0.04, 0.4, 0.01), index="min"), 
                             data.frame(TotalFires = 3.5, Prev.Int = 9, Severity= 0.09, PreNDVI= seq(0.04, 0.4, 0.01), index="med"),
                             data.frame(TotalFires = 7, Prev.Int = 6, Severity= 0.1, PreNDVI= seq(0.04, 0.4, 0.01), index="high"),
                             data.frame(TotalFires = 10.5, Prev.Int = 7, Severity= 0.07, PreNDVI= seq(0.04, 0.4, 0.01), index="max"))
# project model into dataframe
sensitivity_preNDVI_FH$PostDateDif <- 600
sensitivity_preNDVI_FH$threshold <- predict(threshold_rf_index, sensitivity_preNDVI_FH)
# plot
ggplot(sensitivity_preNDVI_FH) +
  geom_line(aes(x=PreNDVI, y=threshold, col=index), size=5) +
  facet_wrap(~index)

# when preNDVI less than .3, you have a lower chance of making it to 100
# at lower freq, you need prendvi > .3 in order to make it to 100
# more likely to reach 70 or higher if you have high or max freq
# for lower ndvi, you are likely to recover to 70. 


# GAMs..................not used................................................................................................................
set.seed(0)
M.maxT.FireHist <- gam(Threshold ~ hist.Int + PDD.6mo +
                         s(FireYear, k=3, by=hist.Int) + 
                         s(TotalFires, by=hist.Int) +
                         s(Severity, by=PDD.6mo) +
                         #s(coords.x1, coords.x2, bs='re') +
                         s(model.NDVI),
                       data =pdsi.total.summary.FH.maxThreshold, method = 'REML')
M.maxT.FireHist2 <- gam(Threshold ~ 
                          s(Prev.Int) +
                          s(coords.x1, coords.x2, bs='re') +
                          s(model.NDVI),
                       data =pdsi.total.summary.FH.maxThreshold, method = 'REML')
summary(M.maxT.FireHist) # 17%
set.seed(1)
gam.check(M.maxT.FireHist)
plot.gam(M.maxT.FireHist, pages=1)

# tried....FireYear, ObsMo, PreNDVI, modelNDVI, hist.Int, PDD.6mo (as by=)

ggplot(data=pdsi.total.summary.FH.maxThreshold, aes(x=Threshold, y=pdsi.min)) +
  geom_point()



set.seed(2)
M.maxT.Climate <- gam(Threshold ~ 
                        ti(pdsi.mean, Obs_Mo) + 
                        #pdsi.max+
                        #pdsi.min + 
                        #model.NDVI +
                        s(coords.x1, coords.x2, bs='re') ,
                      n.rep = 120,
                      data =pdsi.total.summary.FH.maxThreshold, method = 'REML')
summary(M.maxT.Climate)
set.seed(4)
gam.check(M.maxT.Climate)
plot.gam(M.maxT.Climate, pages=1)

set.seed(5)
M.maxT.Total <- gam(Threshold ~ hist.Int + PDD.6mo +
                      s(FireYear, k=3, by=hist.Int) + 
                      s(TotalFires, by=hist.Int) +
                      s(Severity, by=PDD.6mo) +
                      s(model.NDVI) +
                      s(pdsi.mean) +
                      s(pdsi.min) +
                      s(pdsi.max) +
                      s(coords.x1, coords.x2, bs='re'),
                       data =pdsi.total.summary.FH.maxThreshold, method = 'REML')



plotGAM(M.maxT.FireHist, smooth.cov="TotalFires", plotCI=TRUE, groupCovs="hist.Int") 
 # ylim(0, 20) +
  #labs(title = NULL) + 
  #labs(tag = ">5 yrs") +
  #ylab("Recovery Time") + 
  #xlab("Time Since Fire") +
  #scale_x_continuous(n.breaks=10)+
  #theme(text = element_text(size = 20))
plotGAM(M.maxT.FireHist2, smooth.cov="Prev.Int", plotCI=TRUE)
plotGAM(M.maxT.FireHist, smooth.cov="FireYear", plotCI=TRUE, groupCovs="hist.Int")
plotGAM(M.maxT.FireHist, smooth.cov="Severity", plotCI=TRUE, groupCovs="PDD.6mo")
plotGAM(M.maxT.FireHist, smooth.cov="model.NDVI", plotCI=TRUE)

plotGAM(M.maxT.Climate, smooth.cov="pdsi.mean", plotCI=TRUE)
plotGAM(M.maxT.Climate, smooth.cov="pdsi.max", plotCI=TRUE)
plotGAM(M.maxT.Climate, smooth.cov="pdsi.min", plotCI=TRUE)

# save
save(M.maxT.FireHist, file="GAMs.RDATA")





#############################################################################################################################
# Y = RecYears_MaxT 
#############################################################################################################################
# DATA PREP
# make threshold a factor
pdsi.total.summary.FH.maxThreshold$Threshold <- as.factor(pdsi.total.summary.FH.maxThreshold$Threshold)
# filter for post date dif < 6mo
pdsi.total.summary.PPD6mo <- pdsi.total.summary.FH.maxThreshold[which(pdsi.total.summary.FH.maxThreshold$PostDateDif <= 180),]
pdsi.total.summary.PPD6mo$Threshold <- as.factor(pdsi.total.summary.PPD6mo$Threshold)
# pdsi categories
pdsi.total.summary.PPD6mo$pdsiCAT <- "normal"
pdsi.total.summary.PPD6mo$pdsiCAT[pdsi.total.summary.PPD6mo$pdsi.mean <= -2] <- "dry"
pdsi.total.summary.PPD6mo$pdsiCAT[pdsi.total.summary.PPD6mo$pdsi.mean >= 2] <- "wet"
pdsi.total.summary.PPD6mo$pdsiCAT <- as.factor(pdsi.total.summary.PPD6mo$pdsiCAT)
levels(pdsi.total.summary.PPD6mo$pdsiCAT)
# threshold factor order
levels(pdsi.total.summary.PPD6mo$Threshold)
levels(pdsi.total.summary.PPD6mo$pdsiCAT)


# GAMs 
M.RecY_maxT.all <- gam(Rec_Yrs ~ 
                         Threshold + 
                         pdsiCAT +
                         #s(FireYear,k=3) + 
                         ti(TotalFires, ndry.frac)+ 
                         s(Prev.Int, by=threshold)+ 
                         s(Severity, by=Threshold) +
                         s(pdsi.min, by=Threshold) + 
                         s(pdsi.max, by=Threshold) +
                         s(pdsi.mean, by=Threshold) +
                         s(PreNDVI, by=Threshold) +
                         s(coords.x1, coords.x2, bs='re'),
                       data =pdsi.total.summary.PPD6mo, method = 'REML')
summary(M.RecY_maxT.all)
gam.check(M.RecY_maxT.all)

M.RecY.FireHist <- gam(Rec_Yrs ~ 
                         Threshold + 
                         #pdsiCAT +
                         #s(FireYear,k=3) + 
                         ti(TotalFires, ndry.frac)+ 
                         s(Prev.Int, by=threshold)+ 
                         s(Severity, by=Threshold) +
                         #s(pdsi.min, by=Threshold) + 
                         #s(pdsi.max, by=Threshold) +
                         #s(pdsi.mean, by=pdsiCAT) +
                         s(PreNDVI, by=Threshold) +
                         s(coords.x1, coords.x2, bs='re'),
                       data =pdsi.total.summary.PPD6mo, method = 'REML')
summary(M.RecY.FireHist)
gam.check(M.RecY.FireHist)


M.RecY.FireHist.TEST4 <- gam(Rec_Yrs ~ 
                         Threshold + 
                         #pdsiCAT +
                         #s(FireYear,k=3) + 
                         te(TotalFires, ndry.frac)+ 
                         s(Prev.Int, by=threshold)+ 
                         s(Severity, by=Threshold) +
                         #s(pdsi.min, by=Threshold) + 
                         #s(pdsi.max, by=Threshold) +
                         #s(pdsi.mean, by=pdsiCAT) +
                         s(PreNDVI, by=Threshold) +
                         #s(pdsi.mean, bs="re") +
                         s(coords.x1, coords.x2, bs='re'),
                       data =pdsi.total.summary.PPD6mo, method = 'REML')
summary(M.RecY.FireHist.TEST)
gam.check(M.RecY.FireHist.TEST)


ggplot(pdsi.total.summary.PPD6mo) +
  geom_jitter(aes(x=PreNDVI, y=Rec_Yrs, color=Threshold)) +
  facet_wrap(vars(Threshold))


# plots
TF <- plotGAM(M.RecY.FireHist.TEST, smooth.cov="TotalFires", plotCI=TRUE, groupCovs="Threshold") +
  labs(x="Total Fires", y="Recovery Time (yrs)", title="", tag="a", color="Threshold") +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 20)) +
  scale_x_continuous(limits = c(1, 9), n.breaks=6) +
  scale_color_met_d("Tam") 
TSF <- plotGAM(M.RecY.FireHist.TEST, smooth.cov="Prev.Int", plotCI=TRUE, groupCovs="Threshold") +
  labs(x="Time Since Fire (yrs)", y="Recovery Time (yrs)", title="", tag="b", color="Threshold") +
  theme_bw() +
  theme(legend.position = "right",
        text = element_text(size = 20)) +
  scale_color_met_d("Tam") +
  annotate("rect", xmin = 28, xmax = 44, ymin = -0.5, ymax = 7.5,
           alpha = 1,fill = "white")
SEV <- plotGAM(M.RecY.FireHist.TEST, smooth.cov="Severity", plotCI=TRUE, groupCovs="Threshold") +
  labs(x="Severity", y="Recovery Time (yrs)", title="", tag="c", color="Threshold") +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 20)) +  
  scale_x_continuous(limits = c(0, 0.5), n.breaks=4) +
  scale_color_met_d("Tam") 
PreNDVI <- plotGAM(M.RecY.FireHist.TEST, smooth.cov="PreNDVI", plotCI=TRUE, groupCovs="Threshold") +
  labs(x="Pre-Fire NDVI", y="Recovery Time (yrs)", title="", tag="d", color="Threshold") +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 20)) +
  scale_x_continuous(limits = c(0, 0.35), n.breaks=6) +
  scale_color_met_d("Tam") 

# save good models
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Manuscript/figures")
save(M.RecY.FireHist, M.RecY.FireHist.TEST, file="GAM_FireHist.RDATA")





# CLIMATE
# FOR CLIAMTE, TRY GROUPING RECOVERY TIME INTO SLOW,MED,FAST
# RUN GAM WITH REC_TIME AS CATEGORICAL. Keep it general. 
ggplot(pdsi.total.summary.PPD6mo, aes(x=Rec_Yrs))+
  geom_density()
pdsi.total.summary.PPD6mo$RecCat <- "medium"
pdsi.total.summary.PPD6mo$RecCat[pdsi.total.summary.PPD6mo$Rec_Yrs <= 3] <- "fast"
pdsi.total.summary.PPD6mo$RecCat[pdsi.total.summary.PPD6mo$Rec_Yrs > 10] <- "slow"
pdsi.total.summary.PPD6mo$RecCat <- as.factor(pdsi.total.summary.PPD6mo$RecCat)
levels(pdsi.total.summary.PPD6mo$RecCat)
# gam
M.Climate.pdsi.recCat2 <- gam(Rec_Yrs ~ 
                         pdsiCAT +
                        s(coords.x1, coords.x2, bs='re'),
                      data =pdsi.total.summary.PPD6mo, method = 'REML')
summary(M.Climate.pdsi.recCat2)
gam.check(M.Climate.pdsi.recCat2)
# estimated marginal means
m1_emm <- emmeans(M.Climate.pdsi.recCat2, specs = c("pdsiCAT"))
df <- as.data.frame(m1_emm)
#df$RecCat <- factor(df$RecCat, levels=c("<3 year", "3-10 years", ">10 years"))
df$pdsiCAT <- factor(df$pdsiCAT, levels=c("dry", "norm", "wet"))

# plot
ggplot(df, aes(x=pdsiCAT, y=emmean, color=pdsiCAT)) +
  geom_point(size=5) +
  #geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.5, linewidth=.75) +
  scale_color_manual(values=c("#a9845a", "#677853", "#738e8e"))+
  theme_bw() +
  labs(y="Recovery Time (years)", x="Post-Fire Climate", color="Recovery Category") +
  theme(text = element_text(size = 20), 
        legend.position = "none") 


# save climate model
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Manuscript")
save(M.Climate.pdsi.recCat2, file="GAM_Climate2.RDATA")








test <- getViz(M.RecY_maxT.FireHist)
plot(test)

M.RecY_maxT.Climate <- gam(Rec_Yrs ~
                             s(pdsi.mean) + s(pdsi.max)+ s(pdsi.min) + s(nwet)+ s(coords.x1, coords.x2, bs='re') ,
                           data =pdsi.total.summary.FH.maxThreshold, method = 'REML')

summary(M.RecY_maxT.Climate)

library(randomForest)

pdsi.total.summary.FH.maxThreshold %>% summary

M.RecY.rf <- randomForest(Rec_Yrs ~ 
                         Threshold + FireYear + 
                   TotalFires + ndry.frac +
                   Prev.Int + Severity +
                   pdsi.min + pdsi.max +
                   pdsi.mean+PreNDVI,  data = pdsi.total.summary.FH.maxThreshold %>% na.omit )







