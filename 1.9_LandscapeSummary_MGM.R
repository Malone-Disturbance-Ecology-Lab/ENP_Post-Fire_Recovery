## --------------------------------------------- ##
#                Landscape Summary
## --------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script summarizes the (1) fire history and (2) climate patterns for the landscape
## across pineland communitites during the study period. 

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

# Load necessary libraries
# If you don't have the "librarian" package, uncomment the next line and run it to install the package
# install.packages("librarian")
librarian::shelf(sf, tidyverse, psych, terra, tidyterra)

# Point to the path in the server where the data is
dir <- file.path("/", "Volumes", "MaloneLab", "Research", "ENP", "ENP Fire", "Grace_McLeod") 

## --------------------------------------------- ##
#         Fire History by Veg Type -----
## --------------------------------------------- ##

# Master dataframe --------------------------------

projcrs <- "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"

# All vegetation classes
veg_info <- sf::st_read(file.path(dir, "Veg_layers", "veg_info", "Veg_info_intersect.shp")) %>%
  # Make sure projection is correct
  sf::st_transform(projcrs)

# Fire characteristics
# shp from EVG_AllFires_1978-2020 (see Fire_Data_Edits) extracted to upland_pts_sf.shp (see VegLayers)
# in GIS using the intersect tool
fire_character <- sf::st_read(file.path(dir, "Landscape_Summary", "Upland_FireCharacter_sp", "Upland_FireCharacter.shp")) %>%
  # Make sure projection is correct
  sf::st_transform(projcrs)

# Fire history 
load(file.path(dir, "Fire_History", "FireYears_df.RDATA"))
load(file.path(dir, "Fire_History", "FireHistory_df.RDATA"))

# Pull out desired columns from fire history (TF)
names(FireHistory_df)
Fire_by_veg <- FireHistory_df %>% 
  dplyr::select(ptID, x, y, freq_1978_2020) 

# Pair veg info with fire history
Fire_by_veg_v2 <- dplyr::inner_join(veg_info, Fire_by_veg, by = "ptID") %>%
  # Calculate mean fire return interval
  dplyr::mutate(mfri = (2020-1978)/(freq_1978_2020+1))

# Calculate previous interval for all locations
names(FireYears_df)

FireYears_df_v2 <- FireYears_df %>% 
  # Reverse column order
  dplyr::select(rev(order(colnames(.)))) %>%
  # Combine the values of all the year columns together into one string 
  tidyr::unite(year_2023:year_1978, col = "combined", sep = ",", na.rm = T) %>%
  # Separate the first 2 parts into last fire year and penultimate fire year columns
  tidyr::separate_wider_delim(combined, delim = ",",
                              # Identify last fire year
                              # Identify penultimate fire year
                              names = c("LFY", "PenUltFY"),
                              too_few = "align_start",
                              too_many = "drop") %>%
  dplyr::mutate(LFY = as.numeric(LFY),
                PenUltFY = as.numeric(PenUltFY)) %>%
  # Fill out the empty values in LFY with NAs
  dplyr::mutate(LFY = dplyr::case_when(
    nchar(LFY) == 0 ~ NA,
    T ~ LFY
  )) %>%
  # Calculate Prev.Int
  dplyr::mutate(Prev.Int = LFY - PenUltFY,
                # Calculate TSF
                TSF = 2020 - LFY) 

FireYears_sub <- FireYears_df_v2 %>% 
  # Pair with fire_by_veg
  dplyr::select(ptID, LFY, PenUltFY, Prev.Int, TSF) 

# Pair veg info with fire history
Fire_by_veg_v3 <- dplyr::inner_join(FireYears_sub, Fire_by_veg_v2, by = "ptID")

# Fire characteristics
# Pull out desired columns 
names(fire_character)

fire_character_v2 <- fire_character %>%
  dplyr::select(ptID, FIRE_ID, FireNumber, FireName, StartDate, EndDate, Year, FireType) %>%
  # Fix errors
  dplyr::mutate(FireType = dplyr::case_when(
    FireType == "Rx" ~ "RX",
    T ~ FireType
  ))

# merge with fire_by_veg
Fire_by_veg_v4 <- dplyr::inner_join(Fire_by_veg_v3, fire_character_v2, by = "ptID")

# Veg classes
# Clean up veg classes
# Replace undetermined values in Level 7
unique(Fire_by_veg_v4$L7_name)
Fire_by_veg_v5 <- Fire_by_veg_v4 %>%
  dplyr::mutate(L7_name = dplyr::case_when(
    L7_name == "Pine Flatwoods-Shrubs" ~ "Pine Flatwoods-Shrub",
    L7_name == "Pine Upland" ~ "Pine Upland-Undetermined",
    L7_name == "Undetermined" ~ L6_name,
    L7_name == "Pine Flatwood-Saw Palmetto" ~ "Pine Flatwoods-Saw Palmetto",
    T ~ L7_name
  )) %>%
  # Differentiate ecosystem type and dominant community for better visualization
  dplyr::mutate(var = L7_name) %>%
  tidyr::separate_wider_delim(cols = var, delim = "-", names = c("EcoSys", "DomCom"),
                              too_few = "align_start") %>%
  dplyr::mutate(DomCom = dplyr::case_when(
    DomCom == "Graminoids" ~ "Graminoid",
    DomCom == "Shrubs" ~ "Shrub",
    DomCom == "Undetermined" ~ L6_name,
    T ~ DomCom
  )) %>%
  # Level 4
  dplyr::mutate(L4_name = dplyr::case_when(
    L4_name == "Pine Upland-Shrub" ~ "Pine Upland",
    T ~ L4_name
  )) %>%
  # Regional
  dplyr::mutate(Regional_VegCat = "Pine Flatwoods") %>%
  dplyr::mutate(Regional_VegCat = dplyr::case_when(
    L4_name == "Pine Rockland" ~ "Pine Rockland",
    T ~ Regional_VegCat
  )) %>%
  dplyr::mutate(Regional_VegCat = as.factor(Regional_VegCat)) %>%
  # Clean up 
  dplyr::select(-FID_upland, -Uplnds_, -FID_pinela, -geometry.y) %>%
  dplyr::rename(geometry = geometry.x)

names(Fire_by_veg_v5)

# save
save(Fire_by_veg_v5, file = file.path(dir, "Landscape_Summary", "Landscape_Summary_updated.RDATA"))

# Spatial file
Fire_by_veg_sp <- sf::st_as_sf(x = Fire_by_veg_v5,
                               coords = c("x", "y"),
                               crs = projcrs)

sf::st_write(Fire_by_veg_sp, file.path(dir, "Landscape_Summary"), layer = "Fire_by_veg_sp_updated.shp", driver = "ESRI Shapefile",
             append = FALSE)

# Fire history by veg class --------------------------------
mean(Fire_by_veg_v5$freq_1978_2020)
sd(Fire_by_veg_v5$freq_1978_2020)
sd(Fire_by_veg_v5$freq_1978_2020) / sqrt(length(Fire_by_veg_v5$freq_1978_2020))

# Stats 
# Means table
means_table <- Fire_by_veg_v5 %>%
  dplyr::group_by(EcoSys) %>%
  dplyr::summarize(mean_freq = mean(freq_1978_2020))
# Data table
data_table <- psych::describeBy(Fire_by_veg_v5$freq_1978_2020, group=Fire_by_veg_v5$EcoSys, fast=TRUE, digits = 2)
data_table <- psych::describeBy(Fire_by_veg_v5$freq_1978_2020, group=Fire_by_veg_v5$DomCom, fast=TRUE, digits = 2)

# Previous interval
# Means table
means_table <- Fire_by_veg_v5 %>%
  dplyr::group_by(EcoSys) %>%
  dplyr::summarize(mean_freq = mean(Prev.Int))
# Data table
psych::describeBy(Fire_by_veg_v5$Prev.Int, group=Fire_by_veg_v5$EcoSys, fast=TRUE, digits = 2)
psych::describeBy(Fire_by_veg_v5$Prev.Int, group=Fire_by_veg_v5$DomCom, fast=TRUE, digits = 2)

# Time since fire
summary(Fire_by_veg_v5$TSF)
# Means table
means_table <- Fire_by_veg_v5 %>%
  dplyr::group_by(EcoSys) %>%
  dplyr::summarize(mean_freq = mean(TSF))
# Data table
psych::describeBy(Fire_by_veg_v5$TSF, group=Fire_by_veg_v5$EcoSys, fast=TRUE, digits = 2)
psych::describeBy(Fire_by_veg_v5$TSF, group=Fire_by_veg_v5$DomCom, fast=TRUE, digits = 2)

# Unburned 
# Select unburned areas
unburned <- Fire_by_veg_v3 %>% dplyr::filter(freq_1978_2020 ==00)
0.0009 * 2859 # = 2.5731 km2

# Fire type
# Load fire perimeters clipped to pinelands 
Pine_fires_shp <- sf::st_read(file.path(dir, "Fire_History", "Ecosyst_Clipped_Fires", "Pine_fires_Km2.shp"))
head(Pine_fires_shp)
# number of fires
WF <- Pine_fires_shp %>% dplyr::filter(FireType =="WF")
RX <- Pine_fires_shp %>% dplyr::filter(FireType =="RX")
length(unique(WF$FireNumber))
length(unique(RX$FireNumber))
# annual burned area 
sum(WF$Area_km2) # WF: sum= 693.4558
mean(WF$Area_km2) # WF: mean= 0.606698
sum(RX$Area_km2) # RX: sum= 1418.603
mean(RX$Area_km2) # RX: mean= 1.61941
# data table (for means)
psych::describeBy(Fire_by_veg_v5$freq_1978_2020, group=Fire_by_veg_v5$EcoSys, fast=TRUE, digits = 2)
psych::describeBy(Fire_by_veg_v5$freq_1978_2020, group=Fire_by_veg_v5$DomCom, fast=TRUE, digits = 2)

## --------------------------------------------- ##
#              Climate Summaries -----
## --------------------------------------------- ##
# Clean up envr
rm(list = setdiff(ls(), c("dir", "projcrs")))

# Daymet ------------------------------------------

DAYMETcrs <- "+proj=lcc +lat_0=42.5 +lon_0=-100 +lat_1=25 +lat_2=60 +x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs"

# Load Daymet
tmin <- terra::rast(file.path(dir, "Seasonal_Cond", "tmin_EVG_updated.tif"))
tmax <- terra::rast(file.path(dir, "Seasonal_Cond", "tmax_EVG_updated.tif"))
precip <- terra::rast(file.path(dir,"Seasonal_Cond", "precip_EVG_updated.tif"))

# Create a list with all the rasters
temp_precip_list <- list(tmin, tmax, precip)
# Create a list with all the data types
temp_precip_names <- c("tmin", "tmax", "precip")
# Create an empty list to store data for later
df_list <- vector(mode = "list", length = 3)

# Load sample pts
smpl_pts_vect <- terra::vect(file.path(dir, "Sampling", "Sample_pts_upland_052925.shp"))

# Match crs to rasters
terra::crs(smpl_pts_vect)
smpl_pts_DAYMETcrs <- terra::project(smpl_pts_vect, terra::crs(precip))

# Assign dates to raster layers
Date <- read.csv(file = file.path(dir, "Seasonal_Cond", "DAYMET_dates.csv")) %>%
  dplyr::mutate(Yr_Mo = as.Date(paste0(as.character(Date), '01'), format="%Y%m%d"))

# Takes around 18 mins
for (i in 1:length(temp_precip_list)){
  # Get the name of the data we're working with
  data_name <- temp_precip_names[[i]]
  
  message("Currently on: ", data_name)
  
  # Assign layer names to z dimensions
  names(temp_precip_list[[i]]) <- paste(data_name, Date$Yr_Mo, sep=".")
  
  # Extract data to each observation based on month
  smpl_pts_extracted <- terra::extract(temp_precip_list[[i]], smpl_pts_DAYMETcrs, method="bilinear", xy=T, bind=T)
  
  # Convert to dataframe
  smpl_pts_df <- as.data.frame(smpl_pts_extracted) %>%
    # Round to 1 decimal
    dplyr::mutate(dplyr::across(.cols = dplyr::starts_with(data_name),
                                .fns = ~round(., digits = 1))) %>%
    # Remove Uplnds_ column and outdated crds_x1, crds_x2 columns
    dplyr::select(-c(Uplnds_, crds_x1, crds_x2))
  
  # Rotate df longways
  smpl_pts_df_v2 <- smpl_pts_df %>%
    tidyr::pivot_longer(cols = -c(ptID, EcoType, x, y),
                        names_to = "variable",
                        values_to = data_name) %>%
    # Format columns
    # Remove the data name character string from the values in variable
    dplyr::mutate(variable = stringr::str_replace(variable, paste0(data_name, "."), "")) %>%
    # Convert to date format
    dplyr::mutate(variable = as.Date(variable, format = "%Y.%m.%d")) %>%
    # Create columns for the year and month
    dplyr::mutate(Obs_Year = format(variable, "%Y"),
                  Obs_month = format(variable, "%m")) %>%
    # Remove the variable column
    dplyr::select(-variable)
  
  #head(smpl_pts_df_v2)
  
  df_list[[i]] <- smpl_pts_df_v2
  
}

# DAYMET Master
# Join 
Upland_DAYMET <- df_list %>%
  purrr::reduce(dplyr::full_join, by=c("ptID", "EcoType", "x", "y", "Obs_Year", "Obs_month"), keep = F) %>%
  dplyr::relocate(tmin, .before = tmax)

# Save
save(Upland_DAYMET, file = file.path(dir, "Seasonal_Cond", "Upland_DAYMET_updated.RDATA"))

# PDSI ------------------------------------------

# Clean up envr
rm(list = setdiff(ls(), c("dir", "projcrs")))

# Load pdsi
pdsi <- terra::rast(file.path(dir, "Climate", "PDSI_stack.tif"))

# Load sample pts
smpl_pts_vect <- terra::vect(file.path(dir, "Sampling", "Sample_pts_upland_052925.shp"))

# Match crs
smpl_pts_PDSIcrs <- terra::project(smpl_pts_vect, terra::crs(pdsi)) %>%
  # Remove Uplnds_, crds_x1, crds_x2 columns
  dplyr::select(-c(Uplnds_, crds_x1, crds_x2))

# Extract data
# Files are too large, run in chunks

# Divide into chunks
pdsi1 <- pdsi[[1:700]]
pdsi2 <- pdsi[[701:1100]]
pdsi3 <- pdsi[[1101:1533]]

# Put chunks in a list
pdsi_list <- list(pdsi1, pdsi2, pdsi3)

# Takes around 33 mins
# For each chunk...
for (i in 1:length(pdsi_list)){
  
  message("Currently on chunk: ", i, " of ", length(pdsi_list))
  
  # Extract data to sample points
  smpl_pts_pdsi <- terra::extract(pdsi_list[[i]], smpl_pts_PDSIcrs, method = "bilinear", xy = T, bind = T)
  
  # Format data
  smpl_pts_pdsi_df <- as.data.frame(smpl_pts_pdsi) %>%
    # Round to 1 decimal
    dplyr::mutate(dplyr::across(.cols = dplyr::starts_with("X2"),
                                .fns = ~round(., digits = 1)))
  
  smpl_pts_pdsi_df_v2 <- smpl_pts_pdsi_df %>%
    # Rotate df longways
    tidyr::pivot_longer(cols = -c(ptID, EcoType, x, y),
                        names_to = "Obs_Date",
                        values_to = "pdsi") %>%
    dplyr::mutate(Obs_Date = as.Date(Obs_Date, format="X%Y.%m.%d")) %>%
    dplyr::mutate(Obs_Year = format(Obs_Date, "%Y"), .after = Obs_Date)
  
  # calculate means
  means_pdsi <- smpl_pts_pdsi_df_v2 %>%
    dplyr::group_by(Obs_Year) %>%
    dplyr::summarize(pdsi = mean(pdsi))
  
  # Save
  save(smpl_pts_pdsi_df_v2, file = file.path(dir, "Climate", paste0("Uplands_PDSI", i, "_updated.RDATA")))
  save(means_pdsi, file = file.path(dir, "Climate", paste0("means_pdsi", i, "_updated.RDATA")))
  
}

# Combine chunks
# would be great to combine all PDSI data but it's just too big. 
# combine means for figure
load(file.path(dir, "Climate", "means_pdsi1_updated.RDATA"))
means_pdsi1 <- means_pdsi
rm(means_pdsi)

load(file.path(dir, "Climate", "means_pdsi2_updated.RDATA"))
means_pdsi2 <- means_pdsi
rm(means_pdsi)

load(file.path(dir, "Climate", "means_pdsi3_updated.RDATA"))
means_pdsi3 <- means_pdsi
rm(means_pdsi)

# Join
mean_pdsi <- rbind(means_pdsi1, means_pdsi2, means_pdsi3) %>%
  # Round to 1 decimal
  dplyr::mutate(pdsi = round(pdsi, digits = 1))

# Save
save(mean_pdsi, file = file.path(dir, "Climate", "mean_pdsi_combo_updated.RDATA"))

# Climate time series ------------------------------------------

# Clean up envr
rm(list = setdiff(ls(), c("dir")))

load(file.path(dir, "Seasonal_Cond", "Upland_DAYMET_updated.RDATA"))
load(file.path(dir, "Climate", "mean_pdsi_combo_updated.RDATA"))

# DAYMET: calculate annual means
means_daymet <- Upland_DAYMET %>%
  dplyr::group_by(Obs_Year) %>%
  dplyr::summarize(precip = mean(precip),
                   tmax = mean(tmax),
                   tmin = mean(tmin)) %>%
  # Round to one decimal
  dplyr::mutate(dplyr::across(.cols = -Obs_Year,
                              .fns = ~round(., digits = 1))) %>%
  # Convert Obs_Year column to numeric
  dplyr::mutate(Obs_Year = as.numeric(Obs_Year))

mean_pdsi <- mean_pdsi %>%
  # Convert Obs_Year column to numeric
  dplyr::mutate(Obs_Year = as.numeric(Obs_Year))

# Stats
tmax <- lm(means_daymet$tmax ~ means_daymet$Obs_Year)
summary(tmax)

tmin <- lm(means_daymet$tmin ~ means_daymet$Obs_Year)
summary(tmin)

precip <- lm(means_daymet$precip ~ means_daymet$Obs_Year)
summary(precip)
