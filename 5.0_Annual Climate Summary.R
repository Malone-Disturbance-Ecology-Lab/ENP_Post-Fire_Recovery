# Annual CLimate Summary:

library(tidyverse)
library(sf)
library(terra)
# Summarise PDSI................................................................................................................
load('/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Climate/Recov_PDSI.RDATA' ) 

summary.pdsi <- Recov_pdsi %>% 
  mutate(Year = format(Week.pdsi, "%Y")) %>%  
  group_by(Year) %>% 
  summarise( PDSI.mean = mean(pdsi, na.rom=T))

load( '/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Seasonal_Cond/Upland_DAYMET.RDATA')

names(Upland_DAYMET)
summary.climate.prcp <- Upland_DAYMET %>%  
  group_by(Obs_Year, ptID) %>% 
  summarise(PRCP.tot = sum(precip, na.rom=T)) %>% 
  group_by(Obs_Year) %>% 
  summarise(PRCP = mean(PRCP.tot, na.rom=T)) %>% 
  mutate( Year = Obs_Year)
summary.climate.temp <- Upland_DAYMET %>%  
  group_by(Obs_Year) %>% 
  summarise(TMAX.mean = mean(tmax, na.rom=T), TMIN.mean = mean(tmin, na.rom=T)) %>% 
  mutate( Year = Obs_Year)
summary.tot <- summary.climate.temp %>% 
  left_join(summary.pdsi, by = 'Year') %>% 
  left_join(summary.climate.prcp, by='Year')

save(summary.tot, file='/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Climate/Annual_Climate_Summary_ENP.RDATA' )