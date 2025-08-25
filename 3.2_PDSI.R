# PALMER DROUGHT SEVERITY INDEX (PDSI) DATA

# This script calculates mean pdsi values for recovery points at varing 
# intervals pre- and post-fire, using PDSI data on a 5 day temporal resolution.
# 1. Processes PDSI data, converting .nc files to a rastser stack
# 2. Extracts PDSI data to recovery sample points ("Recov_PDSI.RDATA")
# 3. Summarizes pdsi during the recovery period ("RecPrdPDSI_summary.RDATA")


rm(list=ls())

library(ncdf4)
library(terra)
library(tidyverse)

# 1. EXTRACT TO RECOVERY SAMPLE POINTS 
##########################################################################################################################################################
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod")

# EXTRACT PDSI DATA TO RECOVERY SAMPLE POINTS.......................................................................
# load PDSI 
PDSI_stack <- terra::rast("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Climate/PDSI_stack.tif")
# load Recovery sample points (crs="+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")
Recov_smpl_pts <- sf::st_read('/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Sampling/Recov_smpl_pts.shp') %>% sf::st_transform(crs = 4326) # sf file

#  change crs of df to match stack
PDSI_stack <- terra::project(PDSI_stack, 'epsg:4326')


# extract PDSI to points
Recov_pdsi_sp <- terra::extract(PDSI_stack, Recov_smpl_pts, method= "bilinear", bind=T)
Recov_pdsi_df <- as.data.frame (Recov_pdsi_sp)

# round to just 1 decimal
Recov_pdsi_df <- Recov_pdsi_df %>% mutate_at(vars(2:1534), funs(round(., 1)))

# turn DF longways in groups

recov.longer <- function(df, start.col, end.col ){
  
  df1 <- df[, c(1, start.col:end.col )]
  
 
  df2 <-  tidyr::pivot_longer(df1 , cols=c(2: length(df1)), values_to ="Week.pdsi", 
                          names_to = "variable")  
  df3 <- df2 %>% 
    transform(d= stringr::str_replace(variable,"X","")) %>% mutate(
      date=as.Date( d, '%Y.%m.%d' ) ) %>% select(ptID, Week.pdsi, date )
  
  return(df3)
  
}

Recov_pdsi.100 <-  recov.longer(Recov_pdsi_df, 2, 100)
Recov_pdsi.200 <-  recov.longer(Recov_pdsi_df, 101, 200)
Recov_pdsi.300 <-  recov.longer(Recov_pdsi_df, 201, 300)
Recov_pdsi.400 <-  recov.longer(Recov_pdsi_df, 301, 400)
Recov_pdsi.500 <-  recov.longer(Recov_pdsi_df, 401, 500)
Recov_pdsi.600 <-  recov.longer(Recov_pdsi_df, 501, 600)
Recov_pdsi.700 <-  recov.longer(Recov_pdsi_df, 601, 700)
Recov_pdsi.800 <-  recov.longer(Recov_pdsi_df, 701, 800)
Recov_pdsi.900 <-  recov.longer(Recov_pdsi_df, 801, 900)
Recov_pdsi.1000 <-  recov.longer(Recov_pdsi_df, 901, 1000)

Recov_pdsi <- rbind(Recov_pdsi.100,
                    Recov_pdsi.200,
                    Recov_pdsi.300,
                    Recov_pdsi.400,
                    Recov_pdsi.500,
                    Recov_pdsi.600,
                    Recov_pdsi.700,
                    Recov_pdsi.800,
                    Recov_pdsi.900,
                    Recov_pdsi.1000)

Recov_pdsi <- Recov_pdsi %>% rename(PDSI = Week.pdsi) %>% mutate( YrMo.pdsi = date %>% format("%Y-%m"))

head(Recov_pdsi)

rm(Recov_pdsi.100, Recov_pdsi.200,Recov_pdsi.300,
  Recov_pdsi.400, Recov_pdsi.500,Recov_pdsi.600,
  Recov_pdsi.700,Recov_pdsi.800,Recov_pdsi.900,
  Recov_pdsi.1000)

# FILTER TO RECOVERED LOCATIONS ONLY
# load Recov_Time
load(file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/Recov_Drivers_082025.RDATA")

# subset to just recovered points
Recov_pdsi <- subset(Recov_pdsi, ptID %in% Recov_Drivers$ptID)

# Check it out
summary(Recov_pdsi$PDSI)

# save 
save(Recov_pdsi, file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Climate/Recov_PDSI.RDATA") 


##########################################################################################################################################################
# 2. RECOVERY PERIOD SUMMARY STATISTICS 
##########################################################################################################################################################

rm(list=ls())
# Load pdsi data
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Climate/Recov_PDSI.RDATA")
# Load Recovery time data
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/Recov_Drivers_082025.RDATA")

# make reduced df with just ptID and time required to get to each threshold
library(tidyverse)

Recov_Drivers_sub <- Recov_Drivers %>% select(ptID, StartDate, Rec_Date )
Recov_Drivers_sub$StartDate  <- Recov_Drivers_sub$StartDate %>% as.Date

Recov_Drivers_pdsi <- Recov_Drivers_sub %>% full_join( Recov_pdsi , by =c('ptID')) %>% 
  na.omit %>% filter(date > StartDate, date < Rec_Date)

pdsi_summary <- Recov_Drivers_pdsi %>%
  select(ptID, StartDate, Rec_Date, PDSI)  %>%
  reframe(.by=c(ptID),
         pdsi.mean = mean(PDSI, na.rm = TRUE),
         pdsi.var = var(PDSI, na.rm = TRUE),
         pdsi.min = min(PDSI, na.rm = TRUE),
         pdsi.max = max(PDSI, na.rm = TRUE),
         pdsi.sd =sd(PDSI, na.rm = TRUE),
         nwet = case_when( PDSI > 4 ~ length(PDSI), .default = 0),
         ndry = case_when( PDSI < -4 ~ length(PDSI), .default = 0), 
         ntot = length(PDSI),
         pdsi.wet = nwet / ntot * 100, 
         pdsi.dry = ndry / ntot * 100)

save(pdsi_summary, file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Climate/RecPrdPDSI_summary_082025.RDATA")

driver.analysis1 <- merge(Recov_Drivers, pdsi_summary, by="ptID") %>%
  filter(PostDateDif <= 180 ) %>% 
  mutate(
    Prev.Int = Prev.Int %>% as.numeric,
    hist.Int = cut( Prev.Int, breaks = c(1, 10, 80), labels=c("0-10", "11-50")) %>% as.factor,
    PDD.6mo = cut( PostDateDif, breaks = c( 0, 180, 360, 540, 602), labels=c("0-180", "181 - 360", "361 - 540", "+541")) %>% as.factor,
    FireYear.fac = FireYear %>%  as.factor(),
    Obs_Mo = format(Rec_Date,"%m") %>% as.factor(),
    #PreNDVI.cat = cut(PreNDVI, breaks= c(0, 0.2, 0.3, 0.5), labels= c("0-0.2", "0.21-0.3", "0.31-0.43" )) %>% as.factor,
    #modelNDVI.cat = cut(model.NDVI, breaks= c(0, 0.275, 0.312, 0.4), labels= c("0-0.275", "0.276-0.312", "0.313-0.4" )) %>% as.factor,
    pdsi.index =  cut(pdsi.mean, breaks = c(-4, -2 , 2, 3), labels=c("dry", "normal","wet" )) %>% as.factor,
    TotalFires = TotalFires %>%  as.numeric,
    freq.index =  cut(TotalFires, breaks = c(0, 4, 6, 11), labels=c("<4", "4-6",">6" )) %>% as.factor,
    threshold= threshold%>% as.factor %>% droplevels,
    rec.status = recode_factor( threshold, '<50'="<80%", '50'="<80%", '60'="<80%",'70'="<80%",'80'="<80%",'90'=">80%",'100'=">80%")) %>% 
  distinct %>%
  mutate(threshold = factor(threshold, levels = c("100", "90", "80", "70", "60", "50", "<50")))

driver.analysis1$FireType[driver.analysis1$FireType == 'Rx'] <- "RX"
driver.analysis1 %>% names()
# select only desired variables
vars <- c("ptID","x","y",
          "threshold", "Rec_Date", "Rec_Yrs", "rec.status",
          "pdsi.mean", "pdsi.min", "pdsi.max", "pdsi.sd", "pdsi.wet", "pdsi.dry",
          "Prev.Int", "TotalFires", "hist.Int", "FireType", "FireYear",
          "Severity", "PostDateDif", "NDVI", "Obs_Mo", "model.NDVI", "PreNDVI",
          "FireYear.fac","PreNDVI.cat","modelNDVI.cat", "PDD.6mo", "freq.index", "pdsi.index")

driver.analysis2 <- driver.analysis1 %>% 
  select(all_of(vars)) %>%
  na.omit()

save(driver.analysis2, file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/DriverAnalysis_DF_082025.RDATA")