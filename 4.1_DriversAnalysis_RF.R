# DRIVERS ANALYSIS
# SPARKLE MALONE, 2024

# This script evaluates the influence of fire history and climate 
# on recovery threshold and recovery time. 

rm(list=ls())
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")

library(tidyverse)
library(corrplot)
library(randomForest)
library(splitstackshape)
library(ggplot2)
library(caret)
library(ggridges)
library(ggpubr)


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

driver.analysis$FireType[driver.analysis$FireType == 'Rx'] <- "RX"

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


train[, c(3:6, 8:9, 15:19, 25, 35)] %>% names()

T80_rf_index.vsurf <- VSURF(train[, c(3:6, 8:9, 15:19, 25, 35)], train[, 36], ntree = 500,
                       RFimplem = "randomForest", 
                       clusterType = "PSOCK", 
                       verbose = TRUE,
                       ncores = detectCores() - 2, parallel= TRUE)

T80_rf_index.vars <-names( train[, c(3:6, 8:9, 15:19, 25, 35)]) [T80_rf_index.vsurf$varselect.pred] 

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
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Manuscript")
save(T80_rf_index, test, train,T80_rf_index.vars , file="RF_threshold_index.RDATA")

# Recovery Time model: ####
Rc_Yrs_rf_index.vsurf <- VSURF(train[, c(3:6, 8:9, 15:19, 25, 35)], train[, 23], ntree = 500,
                            RFimplem = "randomForest", 
                            clusterType = "PSOCK", 
                            verbose = TRUE,
                            ncores = detectCores() - 2, parallel= TRUE)

Rc_Yrs_rf_index.vars <-names( train[, c(3:6, 8:9, 15:19, 25, 35)]) [Rc_Yrs_rf_index.vsurf$varselect.pred]

train %>% names
Rc_Yrs_rf_index <- randomForest( Rec_Yrs ~ .,
                              data= train %>% select(all_of(c('Rec_Yrs', Rc_Yrs_rf_index.vars) )),
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
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Manuscript")
save(Rc_Yrs_rf_index, test, train,Rc_Yrs_rf_index.vars , file="RF_Rec_Yrs_index.RDATA")

# Data Visualization ####

library(ggplot2)
library(tidyverse)
library(randomForest)

load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Manuscript/DriversData.RDATA")
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Manuscript/RF_threshold_index.RDATA")
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Manuscript/RF_Rec_Yrs_index.RDATA")

#MODELS Summary: https://r-graph-gallery.com/table.html
Rc_Yrs_rf_index.vars 
T80_rf_index.vars

Rc_Yrs_rf_index
T80_rf_index




# Sensitivity Analysis......................................................####

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


library(randomForest)
# save
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Manuscript")
save(sensitivity.df,
     file="Sensitivity_data.RDATA")

# Plots for Rc_Yrs_rf_index:

plot.Rc_Yrs_rf_inde.pdsi.max <-  sensitivity.df %>% filter(target.var == 'pdsi.max') %>%  ggplot(aes(x = pdsi.max, y =Rc_Yrs_rf_index)) + geom_point(color= '#6f948c', alpha=0.5) +  geom_smooth( color= '#6f948c') + theme_bw() + ylab("Recovery (Years)")

plot.Rc_Yrs_rf_inde.pdsi.mean <- sensitivity.df %>% filter(target.var == 'pdsi.mean') %>%  ggplot(aes(x = pdsi.mean, y =Rc_Yrs_rf_index)) + geom_point(color= '#6f948c', alpha=0.5) +  geom_smooth( color= '#6f948c') + theme_bw()+ ylab("Recovery (Years)")

plot.Rc_Yrs_rf_inde.pdsi.min <- sensitivity.df %>% filter(target.var == 'pdsi.min') %>%  ggplot(aes(x = pdsi.min, y =Rc_Yrs_rf_index)) + geom_point(color= '#6f948c', alpha=0.5) +  geom_smooth( color= '#6f948c') + theme_bw() + ylab("Recovery (Years)")

plot.Rc_Yrs_rf_inde.pdsi.sd <- sensitivity.df %>% filter(target.var == 'pdsi.sd') %>%  ggplot(aes(x = pdsi.sd, y =Rc_Yrs_rf_index)) + geom_point(color= '#6f948c', alpha=0.5) +  geom_smooth( color= '#6f948c') + theme_bw() + ylab("Recovery (Years)")

library(ggpubr)

pdsi.grid <- ggarrange( 
  plot.Rc_Yrs_rf_inde.pdsi.max ,
  plot.Rc_Yrs_rf_inde.pdsi.mean,
  plot.Rc_Yrs_rf_inde.pdsi.min, 
  plot.Rc_Yrs_rf_inde.pdsi.sd, 
  nrow=2, ncol=2, labels = c("B", "C", "D", "E")
  
  )

colour_breaks <- c(1, 2, 3, 5, 7, 10, 15, 20)
colours <- c( "darkblue",'#6f948c', "cyan", "goldenrod",'red')

map.recYears <- driver.analysis %>% ggplot( ) + geom_point( aes( x=coords.x1, y = coords.x2, colour =  Rec_Yrs), size=0.3) +
  scale_colour_gradientn(
    limits  = range(driver.analysis$Rec_Yrs),
    colours = colours[c(1, seq_along(colours), length(colours))],
    values  = c(0, scales::rescale(colour_breaks, from = range(driver.analysis$Rec_Yrs)), 1),
  ) 



ggarrange(  ggarrange(map.recYears, labels="A"), pdsi.grid , nrow=1)


# Second Model
varImpPlot(T80_rf_index)

plot.T80_rf_index.PreNDVI <-  sensitivity.df %>% filter(target.var == 'PreNDVI') %>%  ggplot(aes(x = PreNDVI)) + geom_density(aes(colour =T80_rf_index), alpha=0.5)  + theme_bw() + ylab("Density") 

plot.T80_rf_index.pdsi.sd <- sensitivity.df %>% filter(target.var == 'pdsi.sd') %>%  ggplot(aes(x = pdsi.sd)) + geom_density(aes(colour =T80_rf_index), alpha=0.5)  + theme_bw() + ylab("Density") 

plot.T80_rf_index.pdsi.max <-sensitivity.df %>% filter(target.var == 'pdsi.max') %>%  ggplot(aes(x = pdsi.max)) + geom_density(aes(colour =T80_rf_index), alpha=0.5)  + theme_bw() + ylab("Density") 

plot.T80_rf_index.pdsi.min <-sensitivity.df %>% filter(target.var == 'pdsi.min') %>%  ggplot( aes(x = pdsi.min, y = level, fill = T80_rf_index)) + geom_density_ridges_gradient( ) +  theme_ridges() + theme(legend.position = "none")

plot.T80_rf_index.pdsi.mean <- sensitivity.df %>% filter(target.var == 'pdsi.mean') %>%  ggplot( aes(x = pdsi.mean, y = level, fill = T80_rf_index)) + geom_density_ridges_gradient( ) +  theme_ridges() + theme(legend.position = "none")

plot.T80_rf_index.PostDateDif <- sensitivity.df %>% filter(target.var == 'PostDateDif') %>%  ggplot( aes(x = PostDateDif, y = level, fill = T80_rf_index)) + geom_density_ridges_gradient( ) +  theme_ridges() + theme(legend.position = "none")

plot.T80_rf_index.Severity <-sensitivity.df %>% filter(target.var == 'Severity') %>%  ggplot(aes(x = Severity)) + geom_density(aes(colour =T80_rf_index))  + theme_bw() + ylab("Density") 

plot.T80_rf_index.Prev.Int <- sensitivity.df %>% filter(target.var == 'Prev.Int') %>%  ggplot( aes(x = Prev.Int, y = level, fill = T80_rf_index)) + geom_density_ridges_gradient( ) +  theme_ridges() + theme(legend.position = "none")

map.rec.status <- driver.analysis %>% ggplot( ) + geom_point( aes( x=coords.x1, y = coords.x2, colour =  rec.status), size=0.3) 


map.plot <- ggarrange( map.rec.status, labels= c( "A"))

driver.plots <- ggarrange( 
plot.T80_rf_index.PreNDVI,
plot.T80_rf_index.pdsi.sd,
plot.T80_rf_index.pdsi.max,
plot.T80_rf_index.pdsi.mean,
plot.T80_rf_index.pdsi.min,
plot.T80_rf_index.PostDateDif,
plot.T80_rf_index.Severity,
plot.T80_rf_index.Prev.Int, nrow=2, ncol=4, labels= c( "B", "C", "D", "E","F", "G", "H", "I"),common.legend = TRUE)


final.plot.T80 <- ggarrange(map.plot , driver.plots, ncol=1)
