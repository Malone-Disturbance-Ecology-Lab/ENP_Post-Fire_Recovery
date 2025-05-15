# DRIVERS ANALYSIS
# SPARKLE MALONE, GRACE McLEOD, 2024

# not most up to date

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

# NEW IDEAD:
# classification random forest
# sensitivity analysis for the variables

# RANDOM FOREST.........................................................................................................

# the point here is about suggesting improvement fo yo self. 
# Reaching a high threshold means you improved.
pdsi.total.summary.FH.maxThreshold[pdsi.total.summary.FH.maxThreshold == -Inf] <- NA
pdsi.total.summary.FH.maxThreshold.noNA <- pdsi.total.summary.FH.maxThreshold %>%
  na.omit 
pdsi.total.summary.FH.maxThreshold.noNA$Threshold <- as.factor(pdsi.total.summary.FH.maxThreshold.noNA$Threshold)

# Train and Test data
train <- stratified(pdsi.total.summary.FH.maxThreshold.noNA, c("Obs_Mo", "FireYear", "Severity", "PostDateDif", "TotalFires", "Prev.Int",
                                                          "pdsi.mean", "pdsi.max", "pdsi.min",
                                                          "model.NDVI", "PreNDVI"), .6)
test <- anti_join(pdsi.total.summary.FH.maxThreshold.noNA, train)


# Fire History RF
system.time(threshold_rf_FireHist  <- randomForest( Threshold ~   
                                                      PreNDVI+
                                                      Severity +
                                                      PostDateDif + 
                                                      TotalFires +
                                                      Prev.Int,
                                          data= train,
                                          importance=TRUE,
                                          verbose=TRUE,
                                          predicted=TRUE,
                                          keep.inbag=TRUE)) 
# Climate RF
system.time(threshold_rf_Climate  <- randomForest( Threshold ~   
                                                      pdsi.mean+ 
                                                      pdsi.max + 
                                                      pdsi.min +
                                                     #model.NDVI +
                                                     PreNDVI ,
                                                    data= train,
                                                    importance=TRUE,
                                                    verbose=TRUE,
                                                    predicted=TRUE,
                                                    keep.inbag=TRUE)) 
varImpPlot(threshold_rf_FireHist)



# Sensitivity Analysis.............................................................................................................
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
sensitivity_preNDVI <- rbind(data.frame(pdsi.mean = -3, pdsi.min= -3.5, pdsi.max= -2.8, PreNDVI= seq(0.04, 0.4, 0.01), index="dry"), 
                             data.frame(pdsi.mean = 0.45, pdsi.min= -2.9, pdsi.max= 2.7, PreNDVI= seq(0.04, 0.4, 0.01), index="normal"),
                             data.frame(pdsi.mean = 2.75, pdsi.min= 1.9, pdsi.max= 2.9, PreNDVI= seq(0.04, 0.4, 0.01), index="wet"))
sensitivity_preNDVI <- merge(sensitivity_preNDVI, sev.df)

# project model into dataframe
sensitivity_preNDVI$threshold <- predict(threshold_rf_Climate, sensitivity_preNDVI)
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
sensitivity_preNDVI_FH$threshold <- predict(threshold_rf_FireHist, sensitivity_preNDVI_FH)
# plot
ggplot(sensitivity_preNDVI_FH) +
  geom_line(aes(x=PreNDVI, y=threshold, col=index), size=5) +
  facet_wrap(~index)

# when preNDVI less than .3, you have a lower chance of making it to 100
# at lower freq, you need prendvi > .3 in order to make it to 100
# more likely to reach 70 or higher if you have high or max freq
# for lower ndvi, you are likely to recover to 70. 









# GAMs..................................................................................................................................
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
                         s(pdsi.mean, by=pdsiCAT) +
                         s(PreNDVI, by=Threshold) +
                         s(coords.x1, coords.x2, bs='re'),
                       data =pdsi.total.summary.PPD6mo, method = 'REML')

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

M.RecY.Climate <- gam(Rec_Yrs ~ 
                        Threshold + 
                        pdsiCAT +
                        #s(FireYear,k=3) + 
                        #ti(TotalFires, ndry.frac)+ 
                        #s(Prev.Int, by=threshold)+ 
                        #s(Severity, by=Threshold) +
                        s(pdsi.min, by=Threshold) + 
                        s(pdsi.max, by=Threshold) +
                        s(pdsi.mean, by=pdsiCAT) +
                        #(PreNDVI, by=Threshold) +
                        s(coords.x1, coords.x2, bs='re'),
                      data =pdsi.total.summary.PPD6mo, method = 'REML')




summary(M.RecY_maxT.FireHist)
gam.check(M.RecY_maxT.FireHist)

plotGAM(M.RecY_maxT.FireHist, smooth.cov="TotalFires", plotCI=TRUE, groupCovs="pdsiCAT")
plotGAM(M.RecY_maxT.FireHist, smooth.cov="Prev.Int", plotCI=TRUE, groupCovs="Threshold")
plotGAM(M.RecY_maxT.FireHist, smooth.cov="Severity", plotCI=TRUE, groupCovs="Threshold")
plotGAM(M.RecY_maxT.FireHist, smooth.cov="pdsi.mean", plotCI=TRUE, groupCovs="Threshold")
plotGAM(M.RecY_maxT.FireHist, smooth.cov="PreNDVI", plotCI=TRUE, groupCovs="Threshold")
test <- getViz(M.RecY_maxT.FireHist)
plot(test)


M.RecY_maxT.Climate <- gam(Rec_Yrs ~
                             s(pdsi.mean) + s(pdsi.max)+ s(pdsi.min) + s(nwet)+ s(coords.x1, coords.x2, bs='re') ,
                           data =pdsi.total.summary.FH.maxThreshold, method = 'REML')

summary(M.RecY_maxT.Climate)
