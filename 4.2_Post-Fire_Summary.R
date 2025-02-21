# Summary of Fire impacts relative to ppre-fire NDVI

rm(list=ls())

setwd('/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod')
load("./Recovery/Recov_Master.RDATA")

setwd('/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline')
load(file="NDVI_rf.RDATA")

library(dplyr)
library(randomForest)

Recov_Master %>% names
Recov_Master$EndDate <- as.Date(Recov_Master$EndDate)

# Average PreFire NDVI:
summary.preNDVI <- Recov_Master %>% filter( Obs_Date < EndDate) %>% reframe( .by= ptID, NDVI.pre.mean = mean(NDVI))

# Get the first observation date after the fire:
Recov_Master <- Recov_Master %>% mutate( dateDIFF = Obs_Date - EndDate )

Min.Post.Date <- Recov_Master %>% filter( Obs_Date > EndDate)  %>% 
  reframe( .by= ptID, MinDateDIFF = min(dateDIFF)) # Find the first measurement after the fire

# NDVI directly after fire
Recov_Master.2 <- Recov_Master %>% left_join( Min.Post.Date , by='ptID') %>% mutate( NDVI.Model = predict(NDVI_rf,  cur_data()))

NDVI_Post <- Recov_Master.2 %>% filter(dateDIFF == MinDateDIFF ) %>% select(ptID, NDVI, NDVI.Model) %>%  mutate( post.NDVI = NDVI)

summary.NDVI.Fire <- summary.preNDVI %>% left_join(NDVI_Post , by='ptID') %>%  
  mutate(DIFF.NDVI= NDVI-NDVI.pre.mean, 
    Delta.NDVI = (DIFF.NDVI/NDVI.pre.mean)*100 )


save(summary.NDVI.Fire, file ="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/Recov_Summary.RDATA")

# See Summary Statistics of pre- and Postfire NDVI
summary.NDVI.Fire$Delta.NDVI %>% mean
(summary.NDVI.Fire$Delta.NDVI %>% sd)/sqrt(length(summary.NDVI.Fire$post.NDVI))

summary.NDVI.Fire$NDVI.pre.mean %>% range
summary.NDVI.Fire$NDVI.pre.mean %>% mean
summary.NDVI.Fire$NDVI.pre.mean %>% sd / sqrt(length(summary.NDVI.Fire$NDVI.pre.mean ))

summary.NDVI.Fire$post.NDVI %>% mean # Post fire NDVI
(summary.NDVI.Fire$post.NDVI %>% sd)/sqrt(length(summary.NDVI.Fire$post.NDVI))

summary.NDVI.Fire$Delta.NDVI %>% hist
summary.NDVI.Fire$Delta.NDVI %>% range

summary.NDVI.Fire$NDVI.Model %>% range(na.rm=T)
summary.NDVI.Fire$NDVI.Model %>% mean(na.rm=T)
summary.NDVI.Fire$NDVI.Model %>% sd(na.rm=T)/ sqrt(length(summary.NDVI.Fire$NDVI.Model))

summary.NDVI.Fire <- summary.NDVI.Fire %>% mutate( recovery = NDVI.Model - NDVI.pre.mean)

summary.NDVI.Fire$recovery %>% mean(na.rm=T)
summary.NDVI.Fire$recovery %>% sd(na.rm=T)/sqrt(length(summary.NDVI.Fire$post.NDVI))


# Recover difference caomparing pre and model:
100-0.06/0.08*100
