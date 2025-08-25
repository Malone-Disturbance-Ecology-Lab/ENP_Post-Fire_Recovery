# The NDVI MODEL 
# Sparkle L. Malone

# This script merges
# 1. formats data for Random Forest ("BL_train_test.RDATA")
# 2. builds the Random Forest Model for predicting NDVI ("NDVI_rf.RDATA")
# 3. performs a sensitivity analysis to evaluate drivers of NDVI 

rm(list=ls())

library(sf)
library(tidyverse)
library(ggplot2)
library(grDevices)
library(randomForest)
library(stats)
library(corrplot)

# Testung and Training data [RUN ONCE] ####
load(file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/BL_Master_df_08082025.RDATA")

train <- splitstackshape::stratified(BL_Master_filtered.nona , c("Obs_month", "tmin", "tmax", "precip", "TotalFires", "Prev.Int"), .6)

test <- anti_join(BL_Master_filtered.nona , train)

train.pineland <- splitstackshape::stratified(BL_Master_filtered.nona.pineland , c("Obs_month", "tmin", "tmax", "precip", "TotalFires", "Prev.Int"), .6)

test.pineland <- anti_join(BL_Master_filtered.nona.pineland , train)

# save datasets to use for all models
save(train, test,
     train.pineland , test.pineland , file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/BL_train_test_08082025.RDATA")

##########################################################################################################################################################
# 2. RANDOM FOREST NDVI MODEL
##########################################################################################################################################################

# load train and test data

load(file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/BL_Master_df_08082025.RDATA")

load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/BL_train_test_08082025.RDATA")

# Sample size:
BL_Master_filtered.nona$ptID %>% unique %>% length

BL_Master_filtered.nona %>% filter (EcoType == "Pineland") %>% select(ptID) %>% unique %>% count

BL_Master_filtered.nona %>% filter (EcoType == "Hammock") %>% select(ptID) %>% unique %>% count

# Correlation plot Viz
M <-cor(train[, c("Pt.B1max" ,"Pt.B2max" , "Pt.B3max", "Pt.B4max" , "Pt.B5max", "Pt.B7max" , 
                  "Pt.B1min", "Pt.B2min" ,"Pt.B3min", "Pt.B4min", "Pt.B5min", "Pt.B7min",
                  "NIR.SWIR1", "SWIR1.SWIR2", "tmin", "tmax", "precip", 
                  "TotalFires", "Prev.Int")] %>% na.omit ) 
corrplot::corrplot(M, method="circle") 
corrplot::corrplot(M, method="color")
corrplot::corrplot(M, method="number")

# Run Random Forest Model
system.time(NDVI_rf <-  randomForest(NDVI ~  precip + tmax + Obs_month + 
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
                                     keep.inbag=TRUE)) 

system.time(NDVI_rf_pinelands <- randomForest(NDVI ~  precip + tmax + Obs_month + 
                                                Prev.Int + TotalFires  + 
                                                NIR.SWIR1 + SWIR1.SWIR2 + 
                                                Pt.B4max + Pt.B5max + 
                                                Pt.B4min,
                                              data= train.pineland,
                                              ntree= 500,
                                              mtry= 3,
                                              nodesize = 10,  
                                              importance=TRUE,
                                              verbose=TRUE,
                                              predicted=TRUE,
                                              keep.inbag=TRUE)) 

varImpPlot(NDVI_rf)
varImpPlot(NDVI_rf_pinelands)

# apply to testing data
# prediction and confusion matricies
# TESTING DATA
test$NDVI_rf <- predict(NDVI_rf, test)
test$NDVI_rf.pinelands <- predict(NDVI_rf_pinelands, test)
# see corrilation for prediction
summary(lm(test$NDVI ~ test$NDVI_rf))       # R2= 0.7287  
summary(lm(test$NDVI ~ test$NDVI_rf.pinelands)) # R2= 0.7287 

# save the model
#setwd("./Baseline")


save( NDVI_rf, NDVI_rf_pinelands, file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/NDVI_rf_08082025.RDATA")

test %>% ggplot() + geom_point(aes( x=NDVI, y = NDVI_rf, col=EcoType), size=0.1)

test %>% filter(EcoType == "Pineland") %>%  ggplot() + geom_point(aes( x=NDVI, y = NDVI_rf), size=0.1)

test %>% filter(EcoType == "Hammock") %>%  ggplot() + geom_point(aes( x=NDVI, y = NDVI_rf), size=0.1)

# TARGET CONDITIONS.....................................................................................
# calculate target conditions based on ideal fire history 

load(file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/BL_Master_df_08082025.RDATA")

summary(BL_Master_filtered.nona$TotalFires)
BL_Master_filtered.nona$EcoType %>% unique
# pine = historical interval ~ 5yrs
pine <-BL_Master_filtered.nona %>% filter( EcoType == "Pineland")
hammock <- BL_Master_filtered.nona %>% filter( EcoType == "Hammock")

pine$TSF %>% summary
mean(pine$NDVI[pine$TSF <= 10], na.rm=T) # NDVI = 0.2503734
mean(pine$NIR.SWIR1[pine$TSF <= 10]) # NIR:SWIR = 1.357873
mean(pine$SWIR1.SWIR2[pine$TSF <= 10]) # SWIR1:SWIR2 = 1.548661
mean(pine$Pt.B4max[pine$TSF <= 10]) # B4max = 0.3230611

hammock$TSF %>% summary
mean(hammock$NDVI, na.rm=T) # NDVI = 0.2503734
mean(hammock$NIR.SWIR1) # NIR:SWIR = 1.178277
mean(hammock$SWIR1.SWIR2) # SWIR1:SWIR2 = 1.278579
mean(hammock$Pt.B4max) # B4max = 0.383066

# Compare the two systems

##########################################################################################################################################################
# 3. SENSITIVITY ANALYSIS 
##########################################################################################################################################################

# DATASET(S)

load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/BL_train_test_08082025.RDATA")

load( file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/NDVI_rf_08082025.RDATA")

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
  
  vars <- varImpPlot(rf.model) %>% as.data.frame %>% rownames
  
  summary.df <- data_again %>% dplyr::select(any_of(vars)) %>% summarise_all(mean) %>% dplyr::select(!var.name) %>% mutate( summary='mean')
  
  target <- data_again %>% dplyr::select(any_of(var.name))
  
  VOI <- data.frame(VOI =seq( min(target[,1]) , max(target[,1]), 0.01))
  VOI[, var.name] <- VOI$VOI
  VOI <- VOI %>% dplyr::select(  var.name)
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
p.6 <- sensitivity_plot_NDVI( df=Pt.B4min.summary , VOI = Pt.B4min.summary$Pt.B4min , label='Minimum Band 4')
p.7 <- sensitivity_plot_NDVI( df=Pt.B5max.summary , VOI = Pt.B5max.summary$Pt.B5max , label='Maximum Band 5')
p.8 <- sensitivity_plot_NDVI( df=NIR.SWIR1.summary , VOI = NIR.SWIR1.summary$NIR.SWIR1 , label='NIR.SWIR1')
p.9 <- sensitivity_plot_NDVI( df=tmax.summary , VOI = tmax.summary$tmax , label= 'Maximum Air Temperature (°C)')
p.10 <- sensitivity_plot_NDVI( df=precip.summary , VOI = precip.summary$precip , label='Precipitation (mm)')


# One to One plot

one2one <- ggplot(data = data_name)  + geom_point(aes( x= NDVI, y = predicted)) +
  geom_smooth(aes( x= NDVI, y = predicted), method= "lm",col="grey") + theme_bw() +
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
