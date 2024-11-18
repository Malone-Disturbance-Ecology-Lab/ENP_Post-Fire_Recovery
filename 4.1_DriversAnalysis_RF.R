# DRIVERS ANALYSIS
# SPARKLE MALONE, GRACE McLEOD, 2024

# This script evaluates the influence of fire history and climate 
# on recovery threshold and recovery time. 

rm(list=ls())
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")

library(tidyverse)
library(corrplot)
library(randomForest)
library(splitstackshape)
library(ggplot2)


# Load data
load("./Climate/Rec_pdsi_summary.RDATA")
load("./Recovery/Recov_Rates.RDATA")

#############################################################################################################################
# DATA PREP 
#############################################################################################################################
PreNDVI <- Recov_Rate %>% select(ptID, PreNDVI)
pdsi.total.summary.FH.maxThreshold <- merge(pdsi.total.summary.FH.maxThreshold, PreNDVI, by="ptID")

# normalize nwet and ndry by ntotal
# make Prev.Int more reasonable
# factorize Prev.Int (within or outside of historical interval ~ 10 yrs)
# factorize post date difference by 6mo intervals, make fire year a factor, and format pdsi observation date
# bin and factorize preFire NDVI

driver.analysis <- pdsi.total.summary.FH.maxThreshold %>% filter(Threshold != 40, PostDateDif <= 180 ) %>% 
  mutate(nwet.frac = nwet/ntot,
         ndry.frac = ndry/ntot,
         Prev.Int = Prev.Int %>% as.numeric,
         hist.Int = cut( Prev.Int, breaks = c(1, 10, 80), labels=c("0-10", "11-50")) %>% as.factor,
         PDD.6mo = cut( PostDateDif, breaks = c( 0, 180, 360, 540, 602), labels=c("0-180", "181 - 360", "361 - 540", "+541")) %>% as.factor,
         FireYear.fac = FireYear %>%  as.factor(),
         Obs_Mo = format(Obs_Date,"%m") %>% as.factor(),
         RecoveryTime.days = RecoveryTime.days %>% as.numeric() ,
         PreNDVI.cat = cut(PreNDVI, breaks= c(0, 0.2, 0.3, 0.5), labels= c("0-0.2", "0.21-0.3", "0.31-0.43" )) %>% as.factor,
         modelNDVI.cat = cut(model.NDVI, breaks= c(0, 0.275, 0.312, 0.4), labels= c("0-0.275", "0.276-0.312", "0.313-0.4" )) %>% as.factor,
         pdsi.index =  cut(pdsi.mean, breaks = c(-4, -2 , 2, 3), labels=c("dry", "normal","wet" )) %>% as.factor,
         TotalFires = TotalFires %>%  as.numeric,
         freq.index =  cut(TotalFires, breaks = c(0, 4, 6, 11), labels=c("<4", "4-6",">6" )) %>% as.factor,
         Threshold = Threshold %>% as.factor %>% droplevels,
         rec.status = recode_factor( Threshold, '50'="<80%", '60'="<80%",'70'="<80%",'80'="<80%",'90'=">80%",'100'=">80%")) %>% na.omit %>% distinct 

driver.analysis$Threshold %>% levels



setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Manuscript")
save(driver.analysis, file="DriversData.RDATA")

# Train and Test data

train <- stratified(driver.analysis, c("Severity","PostDateDif",
                                       "pdsi.index","freq.index",
                                       "model.NDVI", "PreNDVI", "Threshold"), 0.5)

train <- driver.analysis %>%
  sample_frac(0.5) 

test <- anti_join(driver.analysis, train, by= 'ptID')


driver.analysis %>% summary
train  %>% summary
test  %>% summary

# Correlation Matrix
train %>% names
sum.MT <- train[, c(3:6, 8:10,16:21,23, 25:27)] %>% na.omit()
summary( sum.MT)
M <- cor(sum.MT)
corrplot::corrplot(M, method="circle")

#############################################################################################################################
# Y = MaxThreshold 
#############################################################################################################################
# RANDOM FOREST......................................................................................................... 

# We need models of recovery threshold reached and recovery time to maximum threshold:

# We focus on the the difference between >80 % and lessthan 80% maximum threshold reached.

# Variable Selection:
library(VSURF)
library(parallel)



T80_rf_index.vsurf <- VSURF(train[, c(3:6, 8:9, 14:19, 25, 35)], train[, 36], ntree = 500,
                       RFimplem = "randomForest", 
                       clusterType = "PSOCK", 
                       verbose = TRUE,
                       ncores = detectCores() - 2, parallel= TRUE)

T80_rf_index.vars <-names( train[, c(3:6, 8:9, 14:19, 25, 35)]) [T80_rf_index.vsurf$varselect.pred] 

T80_rf_index <- randomForest( rec.status ~ .,
                              data= train %>% select( c('rec.status', T80_rf_index.vars) ),
                              importance=TRUE,
                              predicted=TRUE,
                              keep.inbag=TRUE)

T80_rf_index 

varImpPlot(T80_rf_index)

train$T80_rf_index <- predict(T80_rf_index , train)
confusionMatrix(train$T80_rf_index, train$rec.status )

# Save Model 
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Manuscript")
save(T80_rf_index, test, train,T80_rf_index.vars , file="RF_threshold_index.RDATA")

# Recovery Time model: ####
Rc_Yrs_rf_index.vsurf <- VSURF(train[, c(3:6, 8:9, 14:19, 25, 35)], train[, 23], ntree = 500,
                            RFimplem = "randomForest", 
                            clusterType = "PSOCK", 
                            verbose = TRUE,
                            ncores = detectCores() - 2, parallel= TRUE)

Rc_Yrs_rf_index.vars <-names( train[, c(3:6, 8:9, 14:19, 25, 35)]) [Rc_Yrs_rf_index.vsurf$varselect.pred]


Rc_Yrs_rf_index <- randomForest( Rec_Yrs ~ .,
                              data= train %>% select( c('Rec_Yrs', T80_rf_index.vars) ),
                              importance=TRUE,
                              predicted=TRUE,
                              keep.inbag=TRUE)

Rc_Yrs_rf_index 

varImpPlot(Rc_Yrs_rf_index)

train$Rc_Yrs_rf_index <- predict(Rc_Yrs_rf_index , train)
confusionMatrix(train$Rc_Yrs_rf_index, train$Rec_Yrs )

# Save Model 
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Manuscript")
save(Rc_Yrs_rf_index, test, train,Rc_Yrs_rf_index.vars , file="RF_Rec_Yrs_index.RDATA")

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
  geom_line(aes(x=PreNDVI, y=Threshold, col=index), size=5) +
  facet_wrap(~index)

# when preNDVI less than .3, you have a lower chance of making it to 100
# at lower freq, you need prendvi > .3 in order to make it to 100
# more likely to reach 70 or higher if you have high or max freq
# for lower ndvi, you are likely to recover to 70. 


#############################################################################################################################
# Y = RecYears_MaxT 
#############################################################################################################################


test <- getViz(M.RecY_maxT.FireHist)
plot(test)

M.RecY_maxT.Climate <- gam(Rec_Yrs ~
                             s(pdsi.mean) + s(pdsi.max)+ s(pdsi.min) + s(nwet)+ s(coords.x1, coords.x2, bs='re') ,
                           data =pdsi.total.summary.FH.maxThreshold, method = 'REML')

summary(M.RecY_maxT.Climate)

library(randomForest)

#create ID column


M.RecY.rf <- randomForest(Rec_Yrs ~ 
                         Threshold + FireYear + 
                   TotalFires + ndry.frac +
                   Prev.Int + Severity +
                   pdsi.min + pdsi.max +
                   pdsi.mean + PreNDVI,  data = train )


setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Manuscript" )
save(train, test, M.RecY.rf, file="RF_Models.RDATA")

varImpPlot(M.RecY.rf )


# Sensitivity Analysis:

sensitivity_Rec_Yrs.Threshold <- data.frame(
  Threshold = unique(train$Threshold)[-7] %>% as.factor)

sensitivity_Rec_Yrs.FireYear <- data.frame(
  FireYear = seq( min(train$FireYear), max(train$FireYear), by=1))
  

sensitivity.model.vars.lower <- train %>% group_by(Threshold, FireYear) %>% summarise(TotalFires = quantile(TotalFires, 0.25),
                                                                                      ndry.frac = quantile(ndry.frac, 0.25 ),
                                                                                      Prev.Int=quantile(Prev.Int, 0.25),
                                                                                      Severity = quantile(Severity, 0.25),
                                                                                      pdsi.min = quantile(pdsi.min, 0.25),
                                                                                      pdsi.max = quantile(pdsi.max, 0.25),
                                                                                      pdsi.mean = quantile(pdsi.mean, 0.25),
                                                                                      PreNDVI = quantile( PreNDVI, 0.25),
                                                                                      Quantile = as.factor(0.25)) %>% as.data.frame()


sensitivity.model.vars.mean <- train %>% group_by(Threshold, FireYear) %>% summarise(TotalFires = quantile(TotalFires, 0.5),
                                                                                      ndry.frac = quantile(ndry.frac, 0.5 ),
                                                                                      Prev.Int=quantile(Prev.Int, 0.5),
                                                                                      Severity = quantile(Severity, 0.5),
                                                                                      pdsi.min = quantile(pdsi.min, 0.5),
                                                                                      pdsi.max = quantile(pdsi.max, 0.5),
                                                                                      pdsi.mean = quantile(pdsi.mean, 0.5),
                                                                                      PreNDVI = quantile( PreNDVI, 0.5),
                                                                                      Quantile = as.factor(0.5)) %>% as.data.frame()



sensitivity.model.vars.high <- train %>% group_by(Threshold, FireYear) %>% summarise(TotalFires = quantile(TotalFires, 0.75),
                                                                                      ndry.frac = quantile(ndry.frac, 0.75 ),
                                                                                      Prev.Int=quantile(Prev.Int, 0.75),
                                                                                      Severity = quantile(Severity, 0.75),
                                                                                      pdsi.min = quantile(pdsi.min, 0.75),
                                                                                      pdsi.max = quantile(pdsi.max, 0.75),
                                                                                      pdsi.mean = quantile(pdsi.mean, 0.75),
                                                                                      PreNDVI = quantile( PreNDVI, 0.75),
                                                                                      Quantile = as.factor(0.75)) %>% as.data.frame()

# This file has the range of conditions of interest:
sensitivity.model.vars <- rbind( sensitivity.model.vars.lower, sensitivity.model.vars.mean, sensitivity.model.vars.high)


# Sevsitivity files by topic
sensitivity.model.vars %>% select(-c(ndry.frac)) %>% cross_join( data.frame(
  ndry.frac = seq( min(train$ndry.frac), max(train$ndry.frac), by=0.05 ))))

sensitivity_Rec_Yrs.TotalFires <- sensitivity.model.vars %>% select( -c(TotalFires)) %>% cross_join( data.frame(
  TotalFires = seq( min(train$TotalFires), max(train$TotalFires), by=3) ))


sensitivity_Rec_Yrs.ndry.frac <- sensitivity.model.vars %>% select(-c(ndry.frac)) %>% cross_join( data.frame(
  ndry.frac = seq( min(train$ndry.frac), max(train$ndry.frac), by=0.05) ))

sensitivity_Rec_Yrs.Prev.Int <- data.frame(
  Prev.Int = seq( min(train$Prev.Int), max(train$Prev.Int) , by=10))
  
sensitivity_Rec_Yrs.Severity <- data.frame(
  Severity = seq( min(train$Severity), max(train$Severity) , by=0.15))

sensitivity_Rec_Yrs.pdsi.min <- data.frame(
  pdsi.min = seq( min(train$pdsi.min), max(train$pdsi.min) , by=1.7))

sensitivity_Rec_Yrs.pdsi.max <- data.frame(
  pdsi.max = seq( min(train$pdsi.max), max(train$pdsi.max) , by=1.7))

sensitivity_Rec_Yrs.pdsi.mean <- data.frame(
  pdsi.mean = seq( min(train$pdsi.mean), max(train$pdsi.mean) , by=1.7))

sensitivity_Rec_Yrs.PreNDVI <- data.frame(
  PreNDVI = seq( min(train$PreNDVI), max(train$PreNDVI) , by=0.1))

sensitivity <- cross_join(sensitivity_Rec_Yrs.FireYear,
                             sensitivity_Rec_Yrs.TotalFires) %>% 
  cross_join(sensitivity_Rec_Yrs.ndry.frac) %>% 
  cross_join(sensitivity_Rec_Yrs.Prev.Int) %>%  
  cross_join(sensitivity_Rec_Yrs.Severity ) %>% 
  cross_join( sensitivity_Rec_Yrs.pdsi.min) %>% 
  cross_join(sensitivity_Rec_Yrs.pdsi.max ) %>% 
  cross_join(sensitivity_Rec_Yrs.pdsi.mean ) %>% 
  cross_join(sensitivity_Rec_Yrs.PreNDVI) %>%  cross_join(sensitivity_Rec_Yrs.Threshold)

sensitivity$Threshold %>% levels

library(randomForest)
train %>% select( Threshold , FireYear , 
                    TotalFires , ndry.frac ,
                    Prev.Int , Severity ,
                    pdsi.min , pdsi.max ,
                    pdsi.mean , PreNDVI) %>% summary

sensitivity %>% select( Threshold , FireYear , 
                    TotalFires , ndry.frac ,
                    Prev.Int , Severity ,
                    pdsi.min , pdsi.max ,
                    pdsi.mean , PreNDVI) %>% summary

load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Manuscript/RF_Models.RDATA" )


sensitivity_Rec_Yrs.TotalFires$M.RecY.rf <- predict( M.RecY.rf, sensitivity_Rec_Yrs.TotalFires )


setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Manuscript" )
save(train, test, M.RecY.rf, sensitivity, file="RF_Models.RDATA")

ggplot(data = sensitivity_Rec_Yrs.TotalFires, aes( x= TotalFires , y=M.RecY.rf, color=Threshold )) + geom_smooth()

ggplot(data = sensitivity_Rec_Yrs.TotalFires, aes( x= TotalFires , y=M.RecY.rf, color=factor(FireYear) )) + geom_smooth()



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