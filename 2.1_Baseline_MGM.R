
##########################################################################################################################################################
# FORMAT DATA FOR RANDOM FOREST
##########################################################################################################################################################
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")

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


# look at relationship between NDVI and other variables 
# Seasonality
ggplot(BL_Master_df, aes(x=Obs_month, y=NDVI)) +
  geom_smooth() +
  ggtitle("all obs")
ggplot(lowNDVI, aes(x=Obs_month, y=NDVI)) +
  geom_smooth() +
  ggtitle("Obs with NDVI < 0.2")
ggplot(HIGH, aes(x=Obs_month, y=NDVI)) +
  geom_jitter() +
  geom_smooth()+
  ggtitle("pts with max above .2")
ggplot(HIGH, aes(x=precip, y=NDVI)) +
  geom_jitter() +
  geom_smooth()+
  ggtitle("precip")
ggplot(HIGH, aes(x=tmax, y=NDVI)) +
  geom_jitter() +
  geom_smooth()+
  ggtitle("tmax")
ggplot(HIGH, aes(x=tmax, y=NDVI)) +
  geom_jitter() +
  geom_smooth()+
  ggtitle("tmin")
# fire history
HIGH$TotalFires <- as.factor(HIGH$TotalFires)
ggplot(HIGH, aes(x=TSF, y=NDVI)) +
  geom_jitter() +
  geom_smooth() +
  ggtitle("Time since fire")
ggplot(HIGH, aes(x=TotalFires, y=NDVI)) +
  geom_jitter() +
  ggtitle("Total Fires")
ggplot(HIGH, aes(x=MFRI, y=NDVI)) +
  geom_jitter() +
  ggtitle("MFRI")
# characterization spec
ggplot(HIGH, aes(x=Pt.B4max, y=NDVI)) +
  geom_jitter() +
  geom_smooth()+
  ggtitle("B4max")
ggplot(HIGH, aes(x=Pt.B4min, y=NDVI)) +
  geom_jitter() +
  geom_smooth()+
  ggtitle("B4min")
ggplot(HIGH, aes(x=NIR.SWIR1, y=NDVI)) +
  geom_jitter() +
  geom_smooth()+
  ggtitle("NIR.SWIR1")
ggplot(HIGH, aes(x=SWIR1.SWIR2, y=NDVI)) +
  geom_jitter() +
  geom_smooth()+
  ggtitle("SWIR1.SWIR2")


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


# CHECK DATA DISTRIBUTION ACROSS CONDITIONS
OG.TSF <- ggplot(BL_Master_df, aes(x=TSF)) +
  geom_density() +
  ggtitle("OG.TSF")
OG.TF <- ggplot(BL_Master_df, aes(x=TotalFires)) +
  geom_density()+
  ggtitle("OG.TF")
OG.precip <- ggplot(BL_Master_df, aes(x=precip)) +
  geom_density()+
  ggtitle("OG.precip")
OG.tmax <- ggplot(BL_Master_df, aes(x=tmax)) +
  geom_density()+
  ggtitle("OG.tmax")
OG.B4 <- ggplot(BL_Master_df, aes(x=B4)) +
  geom_density()+
  ggtitle("OG.B4")
OG.NIR.SWIR <- ggplot(BL_Master_df, aes(x=NIR.SWIR1)) +
  geom_density() +
  ggtitle("OG.NIR.SWIR")

new.TSF <- ggplot(BL_Total_sample, aes(x=TSF)) +
  geom_density() +
  ggtitle("new.TSF")
new.TF <- ggplot(BL_Total_sample, aes(x=TotalFires)) +
  geom_density()+
  ggtitle("new.TF")
new.precip <- ggplot(BL_Total_sample, aes(x=precip)) +
  geom_density()+
  ggtitle("new.precip")
new.tmax <- ggplot(BL_Total_sample, aes(x=tmax)) +
  geom_density()+
  ggtitle("new.tmax")
new.B4 <- ggplot(BL_Total_sample, aes(x=B4)) +
  geom_density()+
  ggtitle("new.B4")
new.NIR.SWIR <- ggplot(BL_Total_sample, aes(x=NIR.SWIR1)) +
  geom_density() +
  ggtitle("new.NIR.SWIR")
# plot all together
library(patchwork)
(OG.TSF + new.TSF) / (OG.TF + new.TF) /(OG.precip + new.precip) / (OG.tmax + new.tmax) / (OG.B4 + new.B4) / (OG.NIR.SWIR + new.NIR.SWIR)

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
# RANDOM FOREST 
##########################################################################################################################################################
# Random Forests is a "black box" machine learning algorithm for classification (response is a factor) and regression (response is numeric). 
# advantageous because it avoids overfitting
# can deal with a large number of variables
# tells importance of each variable as a predictor of the response
# very few assumptions so widely applicable with out much data manipulation

# Some resources:
# How Random Forests work:  https://www.listendata.com/2014/11/random-forest-with-r.html
# Tutorial: https://www.r-bloggers.com/2021/04/random-forest-in-r/
# https://hackernoon.com/random-forest-regression-in-r-code-and-interpretation
# https://towardsdatascience.com/random-forests-an-ensemble-of-decision-trees-37a003084c6c

# load train and test data
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")
load("./Baseline/BL_train_test.RDATA")

# STRATIFY DATASET(S)
# start with a small sample that maintains a distribution across all the data. 
# fit the model. 
# Evaluate the model against the testing data.  
# then increase the sample size to see if your model improves with more data
train.1 <- stratified(train, c( "Obs_month", "tmin", "tmax", "precip", "TotalFires",  "Prev.Int"), .1)
train.25 <- stratified(train, c( "Obs_month", "tmin", "tmax", "precip", "TotalFires", "Prev.Int"), .25)

# NDVI ............................................................................................................................................................... 


# VARIABLE SELECTION / MODEL FITTING 

# STEP 1: Manipulating model parameters
# all variables, default settings
system.time(NDVI_rf_1.1  <- randomForest( NDVI ~   precip + tmax + tmin + Obs_month + 
                                            Prev.Int + TotalFires  + 
                                            NIR.SWIR1 + SWIR1.SWIR2 + 
                                            Pt.B1max + Pt.B2max + Pt.B3max + Pt.B4max + Pt.B5max + Pt.B7max +
                                            Pt.B1min + Pt.B2min + Pt.B3min + Pt.B4min + Pt.B5min + Pt.B7min,
                                        data= train.1,
                                        importance=TRUE,
                                        verbose=TRUE,
                                        predicted=TRUE,
                                        keep.inbag=TRUE)) # VE = 60.95 time= 40 sec
# all variables, change ntree
system.time(NDVI_rf_1.2  <- randomForest( NDVI ~   precip + tmax + tmin + Obs_month + 
                                         Prev.Int + TotalFires  + 
                                          NIR.SWIR1 + SWIR1.SWIR2 + 
                                          Pt.B1max + Pt.B2max + Pt.B3max + Pt.B4max + Pt.B5max + Pt.B7max +
                                          Pt.B1min + Pt.B2min + Pt.B3min + Pt.B4min + Pt.B5min + Pt.B7min,
                                        data= train.1,
                                        ntree= 1000,
                                        importance=TRUE,
                                        verbose=TRUE,
                                        predicted=TRUE,
                                        keep.inbag=TRUE)) # VE = 61.27 # time= 77 sec (does not improve VE, much slower, keep default)
# all variables, change mtry
system.time(NDVI_rf_1.3  <- randomForest( NDVI ~   precip + tmax + tmin + Obs_month + 
                                         Prev.Int + TotalFires  + 
                                          NIR.SWIR1 + SWIR1.SWIR2 + 
                                          Pt.B1max + Pt.B2max + Pt.B3max + Pt.B4max + Pt.B5max + Pt.B7max +
                                          Pt.B1min + Pt.B2min + Pt.B3min + Pt.B4min + Pt.B5min + Pt.B7min,
                                        data= train.1,
                                        ntree= 500,
                                        mtry= 5,
                                        importance=TRUE,
                                        verbose=TRUE,
                                        predicted=TRUE,
                                        keep.inbag=TRUE)) # VE = 60.97 (does not improve VE, keep default)
# all variables, change nodes
# nodes = 10: time= 23 sec, VE= 59.4
# nodes = 15: time= 20 sec, VE= 59.26
# nodes = 20: time= 18 sec, VE= 58.99
# set to 10
system.time(NDVI_rf_1.4  <- randomForest( NDVI ~   precip + tmax + tmin + Obs_month + 
                                            Prev.Int + TotalFires  + 
                                            NIR.SWIR1 + SWIR1.SWIR2 + 
                                            Pt.B1max + Pt.B2max + Pt.B3max + Pt.B4max + Pt.B5max + Pt.B7max +
                                            Pt.B1min + Pt.B2min + Pt.B3min + Pt.B4min + Pt.B5min + Pt.B7min,
                                          data= train.1,
                                          ntree= 500,
                                          mtry= 3,
                                          nodesize = 10,  
                                          importance=TRUE,
                                          verbose=TRUE,
                                          predicted=TRUE,
                                          keep.inbag=TRUE)) # VE = 61.49 
varImpPlot(NDVI_rf_1.4)


# STEP 2:
# To start removing variables, first look at corrilation and remove the highly correlated ones
# Corrilation plots
load("./Baseline/BL_Master_filtered.RDATA")
BL_Master_filtered$Obs_month <- as.numeric(BL_Master_filtered$Obs_month)
BL_Master_noNA <- na.omit(BL_Master_filtered)
head(BL_Master_filtered)
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

# decide between B1-3max
# no B1max: VE= 59.78
# no B2max: VE= 59.75
# no B3max: VE= 59.8
# without any: VE= 60.08 (none very important)
# Remove all 3
system.time(NDVI_rf_2.1  <- randomForest( NDVI ~  precip + tmax + tmin + Obs_month + 
                                            Prev.Int + TotalFires  + 
                                            NIR.SWIR1 + SWIR1.SWIR2 + 
                                            Pt.B4max + Pt.B5max + Pt.B7max +
                                            Pt.B1min + Pt.B2min + Pt.B3min + Pt.B4min + Pt.B5min + Pt.B7min,
                                          data= train.1,
                                          ntree= 500,
                                          mtry= 3,
                                          nodesize = 10,  
                                          importance=TRUE,
                                          verbose=TRUE,
                                          predicted=TRUE,
                                          keep.inbag=TRUE)) # VE= 60.08
# decide between B5max and B7max
# no B5: VE= 60.49
# no B7: VE= 60.45
# without either: VE= 60.72
# does not matter for VE, B5 more important by a little
# remove B7max
system.time(NDVI_rf_2.2  <- randomForest(NDVI ~   precip + tmax + tmin + Obs_month + 
                                           Prev.Int + TotalFires  + 
                                           NIR.SWIR1 + SWIR1.SWIR2 + 
                                           Pt.B4max + Pt.B5max + 
                                           Pt.B1min + Pt.B2min + Pt.B3min + Pt.B4min + Pt.B5min + Pt.B7min,
                                          data= train.1,
                                          ntree= 500,
                                          mtry= 3,
                                          nodesize = 10,  
                                          importance=TRUE,
                                          verbose=TRUE,
                                          predicted=TRUE,
                                          keep.inbag=TRUE)) # VE= 61.3
# decide between B2-7min
# no B2min: VE= 60.45
# no B3min: VE= 60.44
# no B4min: VE= 60.11
# no B5min: VE= 60.6
# no B7min: VE= 60.54
# remove 2 and 7: VE= 60.78
# keep B4min only
system.time(NDVI_rf_2.3  <- randomForest(NDVI ~  precip + tmax + tmin + Obs_month + 
                                           Prev.Int + TotalFires  + 
                                           NIR.SWIR1 + SWIR1.SWIR2 + 
                                           Pt.B4max + Pt.B5max + 
                                           Pt.B4min,
                                          data= train.1,
                                          ntree= 500,
                                          mtry= 3,
                                          nodesize = 10,  
                                          importance=TRUE,
                                          verbose=TRUE,
                                          predicted=TRUE,
                                          keep.inbag=TRUE)) # VE= 60.28
# decide between tmin and tmax
# no tmin: VE=  60.23
# no tmax: VE= 59.43
# similar VE, tmax more important
# remove tmin
system.time(NDVI_rf_2.4  <- randomForest(NDVI ~ precip + tmax + Obs_month + 
                                           Prev.Int + TotalFires  + 
                                           NIR.SWIR1 + SWIR1.SWIR2 + 
                                           Pt.B4max + Pt.B5max + 
                                           Pt.B4min,
                                          data= train.1,
                                          ntree= 500,
                                          mtry= 3,
                                          nodesize = 10,  
                                          importance=TRUE,
                                          verbose=TRUE,
                                          predicted=TRUE,
                                          keep.inbag=TRUE)) # VE=  60.23
# re-do cor plot with remaining variables 
M <-cor(BL_Master_filtered[, c( "Pt.B4max" , "Pt.B5max",
                        "Pt.B4min", 
                         "NIR.SWIR1", "SWIR1.SWIR2",
                         "Obs_month",  "tmax", "precip", 
                         "TotalFires", "Prev.Int",
                         "coords.x1", "coords.x2")]) 
corrplot(M, method="circle") 
corrplot(M, method="color")
corrplot(M, method="number")
# does it imporve to decide between B4max, B5max and NIR:SWIR?
# NIR:SWIR way more important
# no B4max: VE= 59.38
# no B5max: VE= 60.33
# no NIR.SWIR: VE=  60.17
system.time(NDVI_rf_2.5  <- randomForest(NDVI ~ precip + tmax + Obs_month + 
                                           Prev.Int + TotalFires  + 
                                           NIR.SWIR1 + SWIR1.SWIR2 + 
                                           Pt.B4max + Pt.B5max + 
                                           Pt.B4min,
                                          data= train.1,
                                          ntree= 500,
                                          mtry= 3,
                                          nodesize = 10,  
                                          importance=TRUE,
                                          verbose=TRUE,
                                          predicted=TRUE,
                                          keep.inbag=TRUE)) # VE= 61.08
# NIR:SWIR and B4max are pretty close
# but both important. RF can handle it
# keep both. remove B5max
varImpPlot(NDVI_rf_2.4)

# STEP 3: 
# Anything never come up as important? rerun NDVI_rf_2.7 a couple times
# Does it imporve fit to remove it?
# no B4min: VE=  59.65
# no SWIR1:SWIR2: VE= 59.66
# no Obs_month: VE= 53.94 (keep!)
# keep all. 
system.time(NDVI_rf_3.1  <- randomForest(NDVI ~ precip + tmax + Obs_month + 
                                           Prev.Int + TotalFires  + 
                                           SWIR1.SWIR2 + 
                                           Pt.B4max + Pt.B5max + 
                                           Pt.B4min,
                                          data= train.1,
                                          ntree= 500,
                                          mtry= 3,
                                          nodesize = 10,  
                                          importance=TRUE,
                                          verbose=TRUE,
                                          predicted=TRUE,
                                          keep.inbag=TRUE)) # VE= no change to 2.5
# save models
setwd("./Baseline")
save(NDVI_rf_1.1, NDVI_rf_1.2, NDVI_rf_1.3, NDVI_rf_1.4,
     NDVI_rf_2.1, NDVI_rf_2.2, NDVI_rf_2.3, NDVI_rf_2.4, NDVI_rf_2.5,
     file="NDVI_rf_fitting.RDATA")

# STEP 4:
# try adding more vegetation levels 
# LOAD DATA
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")

# load veg data subset to uplands
load("./Veg_layers/Uplands_all.RDATA")
head(Uplands_df)
# load Baseline sample pts
BL_smpl_pts <- rgdal::readOGR(dsn = "./Sampling", layer = "BL_smpl_pts")
BL_df <- as.data.frame(BL_smpl_pts)
# load Baseline data
load("./Baseline/BL_train_test.RDATA")
# EXTRACT OTHER VEGETATION LEVELS
# make copy and rename recovery points
BL_veg_info <- BL_smpl_pts
# check crs
crs(Uplands)     # +proj=utm +zone=17 +datum=NAD83 +units=m +no_defs 
crs(BL_veg_info) # +proj=utm +zone=17 +datum=WGS84 +units=m +no_defs
Uplands <- spTransform(Uplands, "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")
# extract shapefile info for each point
TABLE <- sp::over(BL_veg_info, Uplands)
head(TABLE)
BL_veg_info$VegCode_L <- TABLE$VegCode_L 
BL_veg_info$VegCode_N  <- TABLE$VegCode_N 
BL_veg_info$L1_name  <- TABLE$L1_name   
BL_veg_info$L2_name  <- TABLE$L2_name  
BL_veg_info$L3_name <- TABLE$L3_name 
BL_veg_info$L4_name <- TABLE$L4_name 
BL_veg_info$L5_name <- TABLE$L5_name 
BL_veg_info$L6_name <- TABLE$L6_name
BL_veg_info$L7_name <- TABLE$L7_name 
# view in df
BL_veg_info_df <- as.data.frame(BL_veg_info)
head(BL_veg_info_df)
BL_veg_info_df$coords.x1 <- NULL ; BL_veg_info_df$coords.x2 <- NULL
# INTEGRATE WITH BL DATA
# train
class(train$ptID)
class(BL_veg_info_df$ptID)
BL_veg_info_df$ptID <- as.factor(BL_veg_info_df$ptID)
train <- merge(train, BL_veg_info_df, by="ptID")
train$L3_name <- as.factor(train$L3_name)
train$L4_name <- as.factor(train$L4_name)
train$L5_name <- as.factor(train$L5_name)
train$L6_name <- as.factor(train$L6_name)
head(train)
# test
test<- merge(test, BL_veg_info_df, by="ptID")
train$L7_name <- as.factor(train$L7_name)
test$L3_name <- as.factor(test$L3_name)
test$L4_name <- as.factor(test$L4_name)
test$L5_name <- as.factor(test$L5_name)
test$L6_name <- as.factor(test$L6_name)
test$L7_name <- as.factor(test$L7_name)
head(test)
# save datasets to use for all models
setwd("./Baseline")
save(train, test, file="BL_train_test.RDATA")

# all veg levels
train.1 <- stratified(train, c( "Obs_month", "tmin", "tmax", "precip", "TotalFires", "Prev.Int"), .1)
names(train.1)
system.time(NDVI_rf_4.1  <- randomForest(NDVI ~  precip + tmax + Obs_month + 
                                           Prev.Int + TotalFires  + 
                                           NIR.SWIR1 + SWIR1.SWIR2 + 
                                           Pt.B4max +  
                                           Pt.B3min + Pt.B4min +
                                           L1_name + L2_name + L3_name + L4_name + L5_name + L6_name + L7_name,
                                         data= train.1,
                                         ntree= 500,
                                         mtry= 3,
                                         nodesize = 10,  
                                         importance=TRUE,
                                         verbose=TRUE,
                                         predicted=TRUE,
                                         keep.inbag=TRUE)) # VE= 58.89
varImpPlot(NDVI_rf_4.1)
# not worth adding, doesn't really help



# FULL TRAINING DATA (same combo as 2.5)
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



# pretty plots.......................................................................
# make variable importance into a dataframe
ImpData <- as.data.frame(importance(NDVI_rf))
ImpData$Var.Names <- row.names(ImpData)
ImpData$Var.Names <- factor(ImpData$Var.Names, levels=c("NIR.SWIR1",
                                                        "Pt.B4max",
                                                        "Pt.B5max",
                                                        "precip",
                                                        "TotalFires",
                                                        "Obs_month",
                                                        "tmax",
                                                        "Pt.B4min",
                                                        "Prev.Int",
                                                        "SWIR1.SWIR2") )
# ggplot
ndvi <- ggplot(ImpData, aes(x=Var.Names, y=`%IncMSE`)) +
  geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=`%IncMSE`), color="darkgreen") +
  geom_point(aes(size = IncNodePurity), color="green", alpha=0.6) +
  theme_light() +
  coord_flip() +
  theme(
    legend.position="bottom",
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    text= element_text(size=15),
    axis.ticks.y = element_blank()) +
  xlab("Predictor Variable") +
  ylab("Importance (%IncMSE)") +
  labs(title="NDVI", size="NodePurity")


# TARGET CONDITIONS
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
# VARIABLE MULTICOLINEARITY 
##########################################################################################################################################################
#https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline")

library( corrplot )
library(ggplot2)
# load BL_master_df
load("./BL_Master_filtered.RDATA")

# Correlation plots
names(BL_Master_filtered) # which columns do you want to compare?
BL_Master_filtered$Obs_month <- as.numeric(BL_Master_filtered$Obs_month)
BL_Master_noNA <- na.omit(BL_Master_filtered)

BL_Master_df <- BL_Master_noNA
# all considered variables
M <-cor(BL_Master_df[, c("Pt.B1max" ,"Pt.B2max" , "Pt.B3max", "Pt.B4max" , "Pt.B5max", "Pt.B7max" , 
                         "Pt.B1min", "Pt.B2min" ,"Pt.B3min", "Pt.B4min", "Pt.B5min", "Pt.B7min",
                         "NIR.SWIR1", "SWIR1.SWIR2",
                         "Obs_month", "tmin", "tmax", "precip", 
                         "TotalFires",  "Prev.Int",
                         "coords.x1", "coords.x2")]) 
corrplot(M, method="circle") 
corrplot(M, method="color")
corrplot(M, method="number")
# varaibles used in model 
M <-cor(BL_Master_df[, c( "Pt.B4max" ,"NIR.SWIR1", "SWIR1.SWIR2",
                          "Obs_month", "tmax", "precip", "tmin",
                          "TotalFires","Prev.Int",
                          "coords.x1", "coords.x2")]) 
corrplot(M, method="circle") 
corrplot(M, method="color")
corrplot(M, method="color", tl.col="black", addCoef.col="grey50")

tiff(filename = "BL_corplot.tiff",  
     width = 8, height = 8, res = 600, units = "in",
     compression = "lzw")
corrplot(M, method="color", tl.col="black", addCoef.col="grey50")
dev.off() 



# POST-FIRE NDVI TRENDS...........................................................................................
load("./BL_Master_filtered.RDATA")
# assign frequency categories
summary(BL_Master_filtered$TotalFires)
class(BL_Master_filtered$TotalFires)
summary(BL_Master_filtered$TSF)
class(BL_Master_filtered$TSF)
BL_Master_filtered$TSF <- as.factor(BL_Master_filtered$TSF)
# plot distribution
ggplot(BL_Master_filtered, aes(x=TotalFires)) +
  geom_histogram() +
  scale_x_continuous(breaks=1:9) +
  theme(text = element_text(size = 20))
# plot trend
BL_Master_filtered$TotalFires <- as.factor(BL_Master_filtered$TotalFires)
cols <- c("1" = "#cfccba", "2" = "#e8e5cb", "3" = "#c3d5a3", "4" = "#9bc184",
          "5"= "#669f60", "6"= "#3d7c3c", "7"="#1e5b24", "8"="#1f3d13", "9"="#224313")

ggplot(BL_Master_filtered, aes(x=TSF, 
                      y=NDVI,
                      group=TotalFires,
                      color=TotalFires)) +
  geom_smooth(formula = y ~ s(x, bs = "cs", k=7)) +
  labs(
    x="Years Post-Fire",
    y= "Mean NDVI",
    title = "Baseline Locations",
    color="Total Fires", 
    fill= "Fire Frequency") + 
  #scale_fill_manual(values=cols) +
  #scale_color_manual(values=cols)+
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE) +
  theme_bw()+
  theme(
    text = element_text(size = 20),
    #panel.background = element_rect(fill='transparent'), #transparent panel bg
    #plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    #panel.grid.major = element_blank(), #remove major gridlines
    #panel.grid.minor = element_blank(), #remove minor gridlines
    #legend.background = element_rect(fill='transparent'), #transparent legend bg
    #legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )
# save plot with transparent background
ggsave('FCEprez.png', px, bg='transparent')


ggplot(BL_Master_filtered, aes(x=TotalFires, y=NDVI, color=TotalFires)) +
  geom_boxplot()


# Sensitivity Analysis: ####
library(randomForest)
library(dplyr)

rm(list=ls())

# DATASET(S)
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")
setwd("./Baseline")
load( file="BL_train_test.RDATA")

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



library(ggplot2)

sensitivity_plot_NDVI <- function( df, VOI, label) {
  
  ggplot(data = df)  + 
    geom_smooth(aes( x= VOI, y = NDVI.model), col='black') + theme_bw() +
    ylab('NDVI') + xlab(label) +theme(text = element_text(size = 8))
}

p.1 <- sensitivity_plot_NDVI( df=totalfires.summary , VOI = totalfires.summary$TotalFires , label='Total Fires')
p.2 <- sensitivity_plot_NDVI( df=Prev.Int.summary , VOI = Prev.Int.summary$Prev.Int , label='Previous Fire Interval')
p.3 <- sensitivity_plot_NDVI( df=SWIR1.SWIR2.summary , VOI = SWIR1.SWIR2.summary$SWIR1.SWIR2 , label='SWIR1:SWIR2')
p.4 <- sensitivity_plot_NDVI( df=Obs_month.summary , VOI = Obs_month.summary$Obs_month , label='Month')
p.5 <- sensitivity_plot_NDVI( df=Pt.B4max.summary , VOI = Pt.B4max.summary$Pt.B4max , label='Maximum Band 4')
p.6 <- sensitivity_plot_NDVI( df=Pt.B4min.summary , VOI = Pt.B4min.summary$Pt.B4min , label='Minimum Band 4')
p.7 <- sensitivity_plot_NDVI( df=Pt.B5max.summary , VOI = Pt.B5max.summary$Pt.B5max , label='Maximum Band 5')
p.8 <- sensitivity_plot_NDVI( df=NIR.SWIR1.summary , VOI = NIR.SWIR1.summary$NIR.SWIR1 , label='NIR.SWIR1')
p.9 <- sensitivity_plot_NDVI( df=tmax.summary , VOI = tmax.summary$tmax , label='Maximum Air Temperature')
p.10 <- sensitivity_plot_NDVI( df=precip.summary , VOI = precip.summary$precip , label='Precipitation')

library(ggpubr)

# One to One plot

one2one <- ggplot(data = data_name)  + geom_point(aes( x= NDVI, y = predicted)) +
  geom_smooth(aes( x= NDVI, y = predicted), method= "lm",col="goldenrod") + theme_bw() +
  ylab('Predicted') + xlab('Observed') +theme(text = element_text(size = 8)) + 
  geom_abline(intercept = 0, slope = 1, col="red", linetype="dashed")

setwd('/Users/sm3466/YSE Dropbox/Sparkle Malone/Research/ENP_Post-Fire_Recovery')
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








