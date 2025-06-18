## --------------------------------------------- ##
#                Sampling Design
## --------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script:
## 1. removes unburned locations ("Burned_smpl_pts.shp")
## 2. selects recovery sample points based on fire history ("Recov_smpl_pts.shp")
## 3. selects baseline sample points based on fire history ("BL_smpl_pts.shp")

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

rm(list=ls())

# Load necessary libraries
# If you don't have the "librarian" package, uncomment the next line and run it to install the package
# install.packages("librarian")
librarian::shelf(tidyverse, sf)

# Point to the path in the server where the data is
dir <- file.path("/", "Volumes", "MaloneLab", "Research", "ENP", "ENP Fire", "Grace_McLeod") 

# Load fire history data 
FireHistory_df <- read.csv(file.path(dir, "Fire_History", "FireHistory_df.csv"))

## --------------------------------------------- ##
#             Subset Sample Points -----
## --------------------------------------------- ##

# 1. BURNED POINTS: points that have burned at least once on record (1978-2020)
Burned_pts_df <- FireHistory_df %>%
  # Remove points that never burned (1978-2020)
  dplyr::filter(freq_1978_2020 != 0)  %>% # 411,265 pts (2,859 unburned / 414,124 origninal = >99% have burned on reccord)
  # Rename columns
  dplyr::rename(coords.x1 = x) %>% 
  dplyr::rename(coords.x2 = y)

Burned_smpl_pts <- Burned_pts_df %>%
  # Select only columns of interest
  dplyr::select(ptID, coords.x1, coords.x2) %>%
  # Convert back to spatial
  sf::st_as_sf(coords = c("coords.x1", "coords.x2"),
               crs = "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")

sf::st_crs(Burned_smpl_pts)

# Save
sf::st_write(Burned_smpl_pts, dsn = file.path(dir, "Sampling"), layer = "Burned_smpl_pts_updated", driver = "ESRI Shapefile", append = FALSE)


# 2. RECOVERY POINTS: pts that burned once during recovery window (2001-2007)
Recov_pts_df <- Burned_pts_df %>%
  # Get the points that only burned once (2001-2007)
  dplyr::filter(freq_2001_2007 == 1) # 157,323 pts  

Recov_smpl_pts <- Recov_pts_df %>%
  # Select only columns of interest
  dplyr::select(ptID, coords.x1, coords.x2) %>%
  # Convert back to spatial
  sf::st_as_sf(coords = c("coords.x1", "coords.x2"),
               crs = "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")

sf::st_crs(Recov_smpl_pts)

# Save 
sf::st_write(Recov_smpl_pts, dsn = file.path(dir, "Sampling"), layer = "Recov_smpl_pts_updated", driver = "ESRI Shapefile", append = FALSE)


# 3. BASELINE POINTS: pts that remain unburned from 2010-2020
leftovers <- dplyr::anti_join(Burned_pts_df, Recov_pts_df) # do not want to repeat recov pts
BL_pts_df <- leftovers %>%
  # Get the points that were unburned from 2010-2020
  dplyr::filter(freq_2010_2020 == 0) # 14,981 pts

# convert to spatial
BL_smpl_pts <- BL_pts_df %>%
  # Select only columns of interest
  dplyr::select(ptID, coords.x1, coords.x2) %>%
  # Convert back to spatial
  sf::st_as_sf(coords = c("coords.x1", "coords.x2"),
               crs = "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs")

sf::st_crs(BL_smpl_pts)

# Save 
sf::st_write(BL_smpl_pts, dsn = file.path(dir, "Sampling"), layer = "BL_smpl_pts_updated", driver = "ESRI Shapefile", append = FALSE)

