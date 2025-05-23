## --------------------------------------------- ##
#       Daymet Meterological Data (Part B)
## --------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script:
## 2. stacks monthly layers and masks stacks to the Everglades boundary

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
# If you don't have the "librarian" package, uncomment the next line and run it to install the package
# install.packages("librarian")
librarian::shelf(tidyverse, sf, terra)

# Point to the seasonal conditions folder on server
season_dir <- file.path("/", "Volumes", "MaloneLab", "Research", "ENP", "ENP Fire", "Grace_McLeod", "Seasonal_Cond") 

## --------------------------------------------- ##
#              Data Processing -----
## --------------------------------------------- ##

# LoadEVG boundary for masking
EVGbound <- sf::st_read(file.path("/", "Volumes", "MaloneLab", "Research", "ENP", "ENP Fire", "Grace_McLeod", "AOI", "EVG_bound.shp")) %>%
  # Change projection to match Daymet data
  sf::st_transform("+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs")

# Precipitation ------------------------------------

# Load files
precip.files <- as.list(dir(file.path(season_dir, "DAYMET_monthly_prcp"), "*_ncss.nc$", full.names = T))

# Read in rasters
precip_rasters <- precip.files %>%
  purrr::map(.f = ~terra::rast(.))

# Stack rasters
precip_stack <- terra::rast(precip_rasters) %>%
  terra::project("+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs")

# Mask to EVGbound
precip_EVG <- terra::mask(precip_stack, EVGbound)

# Check it
sub1 <- terra::subset(precip_EVG, 1)
plot(sub1) 
plot(EVGbound, add=TRUE)

# Save 
terra::writeRaster(precip_EVG, "precip_EVG.tif", overwrite = TRUE)

# Max Temperature ------------------------------------
# Load files
tmax.files <- as.list(dir(file.path(season_dir, "DAYMET_monthly_tmax"), "*_ncss.nc$", full.names = T))

# Read in rasters
tmax_rasters <- tmax.files %>%
  purrr::map(.f = ~terra::rast(.))

# Stack rasters
tmax_stack <- terra::rast(tmax_rasters) %>%
  # Change projection
  terra::project("+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs")

# Mask to EVGbound
tmax_EVG <- terra::mask(tmax_stack, EVGbound)

# Check it
sub1 <- terra::subset(tmax_EVG, 1)
plot(sub1) 
plot(EVGbound, add=TRUE)

# Save 
terra::writeRaster(tmax_EVG, "tmax_EVG.tif", overwrite = TRUE)

# Min Temperature ------------------------------------

# Load files
tmin.files <- as.list(dir(file.path(season_dir, "DAYMET_monthly_tmin"), "*_ncss.nc$", full.names = T))

# Read in rasters
tmin_rasters <- tmin.files %>%
  purrr::map(.f = ~terra::rast(.))

# Stack rasters
tmin_stack <- terra::rast(tmin_rasters) %>%
  # Change projection
  terra::project("+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs")

# Mask to EVGbound
tmin_EVG <- terra::mask(tmin_stack, EVGbound)

# Check it
sub1 <- terra::subset(tmin_EVG, 1)
plot(sub1) 
plot(EVGbound, add=TRUE)

# Save 
terra::writeRaster(tmin_EVG, "tmin_EVG.tif", overwrite = TRUE)

