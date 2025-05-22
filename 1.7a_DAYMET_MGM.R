## --------------------------------------------- ##
#       Daymet Meterological Data (Part A)
## --------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script:
## 1. downloads monthly climate data from DAYMET (RUN ONCE)

# ADDITIONAL RESOURCES
# https://daac.ornl.gov/DAYMET/guides/Daymet_V4_Monthly_Climatology.html
# https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1855
# https://cran.r-project.org/web/packages/daymetr/daymetr.pdf
# https://tmieno2.github.io/R-as-GIS-for-Economists/daymet-with-daymetr-and-feddata.html
# https://cran.r-project.org/web/packages/daymetr/vignettes/daymetr-vignette.html

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

rm(list=ls())

# Load necessary libraries
library(daymetr)

# Point to the seasonal conditions folder on server
season_dir <- file.path("/", "Volumes", "MaloneLab", "Research", "ENP", "ENP Fire", "Grace_McLeod", "Seasonal_Cond") 

## --------------------------------------------- ##
#       Importing Data ***RUN ONCE*** -----
## --------------------------------------------- ##

# Download Daymet gridded climate data
# Monthly precip
daymetr::download_daymet_ncss(c(26.5, -81.75, 25, -80.25), 
                              param= 'prcp', 
                              start = 2000, end = 2020, 
                              frequency = "monthly", 
                              path = file.path(season_dir, "DAYMET_monthly_prcp"))    
# Monthly tmax
daymetr::download_daymet_ncss(c(26.5, -81.75, 25, -80.25), 
                              param= 'tmax', 
                              start = 2000, end = 2020, 
                              frequency = "monthly", 
                              path = file.path(season_dir, "DAYMET_monthly_tmax"))    
# Monthly tmin
daymetr::download_daymet_ncss(c(26.5, -81.75, 25, -80.25), 
                              param= 'tmin', 
                              start = 2000, end = 2020, 
                              frequency = "monthly", 
                              path = file.path(season_dir, "DAYMET_monthly_tmin")) 
