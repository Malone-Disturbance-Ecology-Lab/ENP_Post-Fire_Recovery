# DRIVERS ANALYSIS
# SPARKLE MALONE, 2024

# This script evaluates the influence of fire history and climate on recovery thresholdand recovery time. 

# 1. train and test datasets
# 2. runs variable selection and Random Forest model for T80 and Recovery Time ("RF_threshold_index.RDATA")
# 3. runs sensitivity analysis to determine drivers ("Sensitivity_data.RDATA")

rm(list=ls())
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Manuscript")

library(tidyverse)
library(corrplot)
library(randomForest)
library(splitstackshape)
library(ggplot2)
library(caret)
library(ggridges)
library(ggpubr)
library(VSURF)
library(parallel)
library(ggpubr)

# Load data
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Climate/RecPrdPDSI_summary.RDATA")
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/Recov_Drivers.RDATA")

#############################################################################################################################
# 1. DATA PREP 
#############################################################################################################################

# add pdsi data for analysis-ready dataframe
# make Prev.Int more reasonable
# factorize Prev.Int (within or outside of historical interval ~ 10 yrs)
# factorize post date difference by 6mo intervals, make fire year a factor, and format pdsi observation date
# bin and factorize preFire NDVI

driver.analysis1 <- merge(Recov_Drivers, pdsi_summary, by=c("ptID", "StartDate", "Rec_Date")) %>%
  filter(PostDateDif <= 180 ) %>% 
  mutate(
    Prev.Int = Prev.Int %>% as.numeric,
    hist.Int = cut( Prev.Int, breaks = c(1, 10, 80), labels=c("0-10", "11-50")) %>% as.factor,
    PDD.6mo = cut( PostDateDif, breaks = c( 0, 180, 360, 540, 602), labels=c("0-180", "181 - 360", "361 - 540", "+541")) %>% as.factor,
    FireYear.fac = FireYear %>%  as.factor(),
    Obs_Mo = format(Rec_Date,"%m") %>% as.factor(),
    PreNDVI.cat = cut(PreNDVI, breaks= c(0, 0.2, 0.3, 0.5), labels= c("0-0.2", "0.21-0.3", "0.31-0.43" )) %>% as.factor,
    modelNDVI.cat = cut(model.NDVI, breaks= c(0, 0.275, 0.312, 0.4), labels= c("0-0.275", "0.276-0.312", "0.313-0.4" )) %>% as.factor,
    pdsi.index =  cut(pdsi.mean, breaks = c(-4, -2 , 2, 3), labels=c("dry", "normal","wet" )) %>% as.factor,
    TotalFires = TotalFires %>%  as.numeric,
    freq.index =  cut(TotalFires, breaks = c(0, 4, 6, 11), labels=c("<4", "4-6",">6" )) %>% as.factor,
    threshold= threshold%>% as.factor %>% droplevels,
    rec.status = recode_factor( threshold, '<50'="<80%", '50'="<80%", '60'="<80%",'70'="<80%",'80'="<80%",'90'=">80%",'100'=">80%")) %>% 
  distinct %>%
  mutate(threshold = factor(threshold, levels = c("100", "90", "80", "70", "60", "50", "<50")))

driver.analysis1$FireType[driver.analysis1$FireType == 'Rx'] <- "RX"

# select only desired variables
vars <- c("ptID","coords.x1","coords.x2",
          "threshold", "Rec_Date", "Rec_Yrs", "rec.status",
          "pdsi.mean", "pdsi.min", "pdsi.max", "pdsi.sd", "pdsi.wet", "pdsi.dry",
          "Prev.Int", "TotalFires", "hist.Int", "FireType", "FireYear",
          "Severity", "PostDateDif", "NDVI", "Obs_Mo", "model.NDVI", "PreNDVI",
          "FireYear.fac","PreNDVI.cat","modelNDVI.cat", "PDD.6mo", "freq.index", "pdsi.index")
driver.analysis2 <- driver.analysis1 %>% 
  select(all_of(vars)) %>%
  na.omit()

# Train and Test data

# stratify by variables
train1 <- stratified(driver.analysis, c("Severity","PostDateDif",
                                       "pdsi.index","freq.index",
                                       "model.NDVI", "PreNDVI", "threshold"), 0.5)

train <- train1 %>%
  sample_frac(0.5) 

test <- anti_join(driver.analysis, train, by= 'ptID')

driver.analysis %>% summary
train  %>% summary
test  %>% summary

# save data
save(train, test, driver.analysis, file="DriversData.RDATA")

#############################################################################################################################
# 2. DRIVER MODELS
#############################################################################################################################

# We focus on the the difference between >80 % and < 80% maximum threshold reached.

# Recovery Status (T80) model .............................................................................................

# Run VSURF
T80_rf_index.vsurf <- VSURF(train[, c("pdsi.mean", "pdsi.min", "pdsi.max", "pdsi.sd", "pdsi.wet", "pdsi.dry",
                                        "Prev.Int", "TotalFires", "FireType",
                                        "Severity", "PostDateDif","PreNDVI","freq.index")], 
                            train[["rec.status"]],
                            ntree = 500,
                            RFimplem = "randomForest", 
                            clusterType = "PSOCK", 
                            verbose = TRUE,
                            ncores = detectCores() - 2, parallel= TRUE)


# Get names of selected variables
T80_rf_index.vars <-names( train[, c("pdsi.mean", "pdsi.min", "pdsi.max", "pdsi.sd", "pdsi.wet", "pdsi.dry",
                                     "Prev.Int", "TotalFires", "FireType",
                                     "Severity", "PostDateDif","PreNDVI","freq.index")]) [T80_rf_index.vsurf$varselect.pred] 

# Random Forest
T80_rf_index <- randomForest( rec.status ~ .,
                              data= train %>% select(c(rec.status, all_of(T80_rf_index.vars))),
                              importance=TRUE,
                              predicted=TRUE,
                              keep.inbag=TRUE)

T80_rf_index 

varImpPlot(T80_rf_index)

train$T80_rf_index <- predict(T80_rf_index , train)
confusionMatrix(train$T80_rf_index, train$rec.status )

# Save Model 
save(T80_rf_index, test, train,T80_rf_index.vars , file="RF_threshold_index.RDATA")


# Recovery Time model...............................................................................................

# Run VSURF
Rc_Yrs_rf_index.vsurf <- VSURF(train[, c("pdsi.mean", "pdsi.min", "pdsi.max", "pdsi.sd", "pdsi.wet", "pdsi.dry",
                                         "Prev.Int", "TotalFires", "FireType",
                                         "Severity", "PostDateDif","PreNDVI","freq.index")], 
                               train[["Rec_Yrs"]], 
                               ntree = 500,
                               RFimplem = "randomForest", 
                               clusterType = "PSOCK", 
                               verbose = TRUE,
                               ncores = detectCores() - 2, parallel= TRUE)

# Get names of selected variables
Rc_Yrs_rf_index.vars <-names( train[, c("pdsi.mean", "pdsi.min", "pdsi.max", "pdsi.sd", "pdsi.wet", "pdsi.dry",
                                        "Prev.Int", "TotalFires", "FireType",
                                        "Severity", "PostDateDif","PreNDVI","freq.index")]) [Rc_Yrs_rf_index.vsurf$varselect.pred]

# Random Forest
Rc_Yrs_rf_index <- randomForest( Rec_Yrs ~ .,
                              data= train %>% select(c(Rec_Yrs, all_of(Rc_Yrs_rf_index.vars))),
                              importance=TRUE,
                              predicted=TRUE,
                              keep.inbag=TRUE)


Rc_Yrs_rf_index 

varImpPlot(Rc_Yrs_rf_index)

train$Rc_Yrs_rf_index <- predict(Rc_Yrs_rf_index , train)
test$Rc_Yrs_rf_index <- predict(Rc_Yrs_rf_index , test)

lm( train$Rec_Yrs ~ train$Rc_Yrs_rf_index) %>% summary
lm( test$Rec_Yrs ~ test$Rc_Yrs_rf_index) %>% summary

# Save Model 
save(Rc_Yrs_rf_index, test, train, Rc_Yrs_rf_index.vars , file="RF_Rec_Yrs_index.RDATA")



#############################################################################################################################
# 3. SENSITIVITY ANALYSIS
#############################################################################################################################

load("./DriversData.RDATA")
load("./RF_threshold_index.RDATA")
load("./RF_Rec_Yrs_index.RDATA")


# Sensitivity Analysis...................................................

sensitivity.mean <- driver.analysis %>% summarise(
  pdsi.sd = mean(pdsi.sd, na.rm=T),
  pdsi.max = mean(pdsi.max , na.rm=T), 
  pdsi.min = mean(pdsi.min , na.rm=T), 
  pdsi.mean = mean(pdsi.mean , na.rm=T),
  PreNDVI = mean(PreNDVI , na.rm=T),
  Prev.Int = mean(PreNDVI , na.rm=T),
  PostDateDif = mean(PostDateDif , na.rm=T),
  Severity = mean(Severity , na.rm=T)) %>% mutate(level="mean")

sensitivity.max <- driver.analysis %>% summarise(
  pdsi.sd = max(pdsi.sd, na.rm=T),
  pdsi.max = max(pdsi.max , na.rm=T), 
  pdsi.mean = max(pdsi.mean , na.rm=T), 
  pdsi.min = max(pdsi.min , na.rm=T),
  PreNDVI = max(PreNDVI , na.rm=T),
  Prev.Int = max(PreNDVI , na.rm=T),
  PostDateDif = max(PostDateDif , na.rm=T),
  Severity = max(Severity , na.rm=T))%>% mutate(level="max")

sensitivity.min <- driver.analysis %>% summarise(
  pdsi.sd = min(pdsi.sd, na.rm=T),
  pdsi.max = min(pdsi.max , na.rm=T), 
  pdsi.mean = min(pdsi.mean , na.rm=T), 
  pdsi.min = min(pdsi.min , na.rm=T),
  PreNDVI = min(PreNDVI , na.rm=T),
  Prev.Int = min(PreNDVI , na.rm=T),
  PostDateDif = min(PostDateDif , na.rm=T),
  Severity = min(Severity , na.rm=T))%>% mutate(level="min")

sensitivity.mean <- sensitivity.mean %>% rbind( sensitivity.max, sensitivity.min)


summary.var <- function( var){
  var.df = data.frame( var =seq(min(var), max(var), ((max(var)- min(var))/10))) %>% as.data.frame()
  return(var.df)
}
summary.df <- function( mean.conditions, target.var, df){
  
  target.var.df = df %>% select(target.var) %>% summary.var() %>% mutate(target.var = var ) %>% select(target.var)
  summary.df <- mean.conditions %>% select(-c(target.var)) %>%cross_join(target.var.df)
  return(summary.df)
}

sensitivity.mean.pdsi.sd <- summary.df(mean.conditions= sensitivity.mean, 
                                       target.var = 'pdsi.sd', df=driver.analysis ) %>% mutate(pdsi.sd=target.var ,
                                                                                               target.var = 'pdsi.sd') 

sensitivity.mean.pdsi.max <- summary.df(mean.conditions= sensitivity.mean, 
                                        target.var = 'pdsi.max', df=driver.analysis ) %>% mutate(pdsi.max=target.var ,
                                                                                                 target.var = 'pdsi.max') 

sensitivity.mean.pdsi.mean <- summary.df(mean.conditions= sensitivity.mean, 
                                         target.var = 'pdsi.mean', df=driver.analysis )%>% mutate(pdsi.mean=target.var ,
                                                                                                  target.var = 'pdsi.mean') 

sensitivity.mean.pdsi.min <- summary.df(mean.conditions= sensitivity.mean, 
                                        target.var = 'pdsi.min', df=driver.analysis ) %>% mutate(pdsi.min=target.var,
                                                                                                 target.var = 'pdsi.min') 



sensitivity.mean.PreNDVI <- summary.df(mean.conditions= sensitivity.mean, 
                                       target.var = 'PreNDVI', df=driver.analysis ) %>% mutate(PreNDVI=target.var,
                                                                                               target.var = 'PreNDVI') 

sensitivity.mean.Prev.Int<- summary.df(mean.conditions= sensitivity.mean, 
                                       target.var = 'Prev.Int', df=driver.analysis ) %>% mutate(Prev.Int=target.var,
                                                                                                target.var = 'Prev.Int') 

sensitivity.mean.PostDateDif<- summary.df(mean.conditions= sensitivity.mean, 
                                          target.var = 'PostDateDif', df=driver.analysis ) %>% mutate(PostDateDif=target.var,
                                                                                                      target.var = 'PostDateDif') 


sensitivity.mean.Severity <- summary.df(mean.conditions= sensitivity.mean, 
                                        target.var = 'Severity', df=driver.analysis ) %>% mutate(Severity=target.var,
                                                                                                 target.var = 'Severity') 


sensitivity.df <- rbind( sensitivity.mean.pdsi.sd,sensitivity.mean.pdsi.max,
                         sensitivity.mean.pdsi.mean, sensitivity.mean.pdsi.min, 
                         sensitivity.mean.Severity, sensitivity.mean.PostDateDif,
                         sensitivity.mean.Prev.Int, sensitivity.mean.PreNDVI)

sensitivity.df$Rc_Yrs_rf_index <-  predict(Rc_Yrs_rf_index ,sensitivity.df )
sensitivity.df$T80_rf_index <-  predict(T80_rf_index ,sensitivity.df)


rm( sensitivity.mean.pdsi.sd,sensitivity.mean.pdsi.max,
    sensitivity.mean.pdsi.mean, sensitivity.mean.pdsi.min, 
    sensitivity.mean.Severity,sensitivity.mean.PostDateDif,
    sensitivity.mean.Prev.Int,sensitivity.mean.PreNDVI)


# save
save(sensitivity.df, file="Sensitivity_data.RDATA")


