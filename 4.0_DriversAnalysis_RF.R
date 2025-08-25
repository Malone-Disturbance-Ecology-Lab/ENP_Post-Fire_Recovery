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
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/Recov_Drivers_082025.RDATA")

# Train and Test data [Run Once] ####

# stratify by variables
train <- stratified(driver.analysis2, c("Severity","PostDateDif",
                                       "pdsi.index","freq.index",
                                       "model.NDVI", "PreNDVI", "threshold"), 0.5)

#train <- train1 %>% sample_frac(0.75) 

test <- anti_join(driver.analysis2, train, by= 'ptID')

driver.analysis2 %>% summary
train  %>% summary
test  %>% summary

# save data
save(train, test, driver.analysis2, file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/DriversData_082025.RDATA")

#############################################################################################################################
# 2. DRIVER MODELS
#############################################################################################################################

load( file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/DriversData_082025.RDATA")

# Sample Size:
driver.analysis2$ptID %>% unique %>% length
# Recovery Status (T80) model 

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

T80_rf_index.vars <- c("pdsi.min", "pdsi.max","PreNDVI")
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
save(T80_rf_index, test, train,T80_rf_index.vars , file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/RF_threshold_index_082025.RDATA")


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

Rc_Yrs_rf_index.vars <- c("pdsi.min", "pdsi.max","pdsi.mean", "pdsi.sd")

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
save(Rc_Yrs_rf_index, test, train, Rc_Yrs_rf_index.vars , file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/RF_Rec_Yrs_index_082025.RDATA")



#############################################################################################################################
# 3. SENSITIVITY ANALYSIS
#############################################################################################################################

load(file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/DriversData_082025.RDATA")

# Sensitivity Analysis...................................................

sensitivity.mean <- driver.analysis2 %>% summarise(
  pdsi.sd = mean(pdsi.sd, na.rm=T),
  pdsi.max = mean(pdsi.max , na.rm=T), 
  pdsi.min = mean(pdsi.min , na.rm=T), 
  pdsi.mean = mean(pdsi.mean , na.rm=T),
  PreNDVI = mean(PreNDVI , na.rm=T),
  Prev.Int = mean(PreNDVI , na.rm=T),
  PostDateDif = mean(PostDateDif , na.rm=T),
  Severity = mean(Severity , na.rm=T)) %>% mutate(level="mean")

sensitivity.max <- driver.analysis2 %>% summarise(
  pdsi.sd = max(pdsi.sd, na.rm=T),
  pdsi.max = max(pdsi.max , na.rm=T), 
  pdsi.mean = max(pdsi.mean , na.rm=T), 
  pdsi.min = max(pdsi.min , na.rm=T),
  PreNDVI = max(PreNDVI , na.rm=T),
  Prev.Int = max(PreNDVI , na.rm=T),
  PostDateDif = max(PostDateDif , na.rm=T),
  Severity = max(Severity , na.rm=T))%>% mutate(level="max")

sensitivity.min <- driver.analysis2 %>% summarise(
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
                                       target.var = 'pdsi.sd', df=driver.analysis2 ) %>% mutate(pdsi.sd=target.var ,
                                                                                               target.var = 'pdsi.sd') 

sensitivity.mean.pdsi.max <- summary.df(mean.conditions= sensitivity.mean, 
                                        target.var = 'pdsi.max', df=driver.analysis2 ) %>% mutate(pdsi.max=target.var ,
                                                                                                 target.var = 'pdsi.max') 

sensitivity.mean.pdsi.mean <- summary.df(mean.conditions= sensitivity.mean, 
                                         target.var = 'pdsi.mean', df=driver.analysis2 )%>% mutate(pdsi.mean=target.var ,
                                                                                                  target.var = 'pdsi.mean') 

sensitivity.mean.pdsi.min <- summary.df(mean.conditions= sensitivity.mean, 
                                        target.var = 'pdsi.min', df=driver.analysis2 ) %>% mutate(pdsi.min=target.var,
                                                                                                 target.var = 'pdsi.min') 



sensitivity.mean.PreNDVI <- summary.df(mean.conditions= sensitivity.mean, 
                                       target.var = 'PreNDVI', df=driver.analysis2 ) %>% mutate(PreNDVI=target.var,
                                                                                               target.var = 'PreNDVI') 

sensitivity.mean.Prev.Int<- summary.df(mean.conditions= sensitivity.mean, 
                                       target.var = 'Prev.Int', df=driver.analysis2 ) %>% mutate(Prev.Int=target.var,
                                                                                                target.var = 'Prev.Int') 

sensitivity.mean.PostDateDif<- summary.df(mean.conditions= sensitivity.mean, 
                                          target.var = 'PostDateDif', df=driver.analysis2 ) %>% mutate(PostDateDif=target.var,
                                                                                                      target.var = 'PostDateDif') 


sensitivity.mean.Severity <- summary.df(mean.conditions= sensitivity.mean, 
                                        target.var = 'Severity', df=driver.analysis2 ) %>% mutate(Severity=target.var,
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
save(sensitivity.df, file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/Sensitivity_data.RDATA")


