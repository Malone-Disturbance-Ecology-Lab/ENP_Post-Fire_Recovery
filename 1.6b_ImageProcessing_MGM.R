## --------------------------------------------- ##
#       Landsat Image Processing (Part B)
## --------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script does the following steps to process images for 4 ARD Landsat tiles covering the Everglades

## Images were downloaded manually from Earth Explorer: https://earthexplorer.usgs.gov/

## 2. extracts spectral data from tile-specific stacks to upland sample points and filters data for clouds and QAQC
## 3. merges all tile-specific dataframes to make a Spectral Master dataframe and calculates spectral indices ("Spec_Master.RDATA")

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

# Load necessary libraries
# If you don't have the "librarian" package, uncomment the next line and run it to install the package
# install.packages("librarian")
librarian::shelf(sf, tidyverse, terra)

## --------------------------------------------- ##
#             Image Processing -----
## --------------------------------------------- ##
# Pull out tifs from new directory 
image_dir <- file.path("/", "Volumes", "MaloneLab", "Research", "ENP", "ENP Fire", "Grace_McLeod", "Image_Processing") 
LStif <-  file.path(image_dir, "LS_tif")

# Load upland sample points 
Smpl_pts <- terra::vect(file.path("/", "Volumes", "MaloneLab", "Research", "ENP", "ENP Fire", "Grace_McLeod", "Sampling", "Sample_pts_upland_052925.shp"))

# Tile 26_18 ------------------------------------------------------

# Group by tile
tile_26_18 <- as.list(list.files(LStif, pattern = glob2rx("*_026018*.TIF$"), full.names = TRUE))

# Read in rasters as a list
rast_tile_26_18 <- tile_26_18 %>%
  purrr::map(.f = ~terra::rast(.), .progress = T)

# Stack them
stack_26_18 <- terra::rast(rast_tile_26_18)

# Make sure sample points have the same crs as the rasters
Smpl_pts_26_18 <- Smpl_pts %>%
  terra::project(terra::crs(stack_26_18))

# Subset sample points to tile
Sub_pts_26_18 <- terra::crop(Smpl_pts_26_18, terra::ext(stack_26_18))

# Extract to tile sample points
# Takes a while (30 mins)
Ext_26_18_df <- terra::extract(stack_26_18, Sub_pts_26_18, xy = TRUE, ID = FALSE)

# Combine with uplands vector information
Ext_26_18_df_v2 <- as.data.frame(cbind(Sub_pts_26_18, Ext_26_18_df)) %>%
  # Remove the outdated coordinates columns (since we've reprojected the shapefile)
  dplyr::select(-crds_x1, -crds_x2) %>%
  # Rename the new coordinates columns as coords.x1, coords.x2
  dplyr::rename(coords.x1 = x,
                coords.x2 = y)

#glimpse((Ext_26_18_df_v2))

melt_26_18 <- Ext_26_18_df_v2 %>%
  # Pivot df longways
  tidyr::pivot_longer(cols = -c(Uplnds_, ptID, EcoType, coords.x1, coords.x2), 
                      names_to = "variable", 
                      values_to = "value") %>%
  # Convert variable column to character
  dplyr::mutate(variable = as.character(variable)) %>%
  # Create Tile, Obs_Date, Band columns
  dplyr::mutate(Tile = stringr::str_split_fixed(variable, '_', 8)[,3],
                Obs_Date = stringr::str_split_fixed(variable, '_', 8)[,4],
                Band = stringr::str_split_fixed(variable, '_', 8)[,8]) %>%
  # Get rid of unnecessary column 
  dplyr::select(-variable)

# Turn band into columns
Spec_26_18 <- melt_26_18 %>%
  tidyr::pivot_wider(names_from = "Band", values_from = "value")

Spec_26_18_clean <- Spec_26_18 %>% 
  # Remove rows where band values are NA 
  # If value is NA for one band, it will be for all
  tidyr::drop_na(B1, B2, B3, B4, B5, B7) %>%
  # Make all flagged values = NA
  # If CLOUD_QA > 0, make it NA
  dplyr::mutate(CLOUD_QA = dplyr::case_when(
    CLOUD_QA > 0 ~ NA,
    T ~ CLOUD_QA
  )) %>%
  # If PIXEL != 5440, make it NA
  # 5440 is code for totally clear pixels
  dplyr::mutate(PIXEL = dplyr::case_when(
    PIXEL != 5440 ~ NA,
    T ~ PIXEL
  )) %>%
  # If RADSET > 0, make it NA
  dplyr::mutate(RADSAT = dplyr::case_when(
    RADSAT > 0 ~ NA,
    T ~ RADSAT
  )) %>%
  # Remove rows with flagged values
  tidyr::drop_na(CLOUD_QA, PIXEL, RADSAT) %>%
  # Remove unnecessary columns
  dplyr::select(-c(RADSAT, PIXEL, LINEAGE, CLOUD_QA, ATMOS_OPACITY))

# Save  
save(Ext_26_18_df_v2, Spec_26_18_clean, file = file.path(image_dir, "Master_Spec", "Tile_26_18_updated.RDATA"))

# Clean up envr before next tile
rm(list = setdiff(ls(), c("Smpl_pts", "image_dir", "LStif")))
gc()

# Tile 26_19 ------------------------------------------------------

# Group by tile
tile_26_19 <- as.list(list.files(LStif, pattern = glob2rx("*_026019*.TIF$"), full.names = TRUE))

# Read in rasters as a list
rast_tile_26_19 <- tile_26_19 %>%
  purrr::map(.f = ~terra::rast(.), .progress = T)

# Stack them
stack_26_19 <- terra::rast(rast_tile_26_19)

# Make sure sample points have the same crs as the rasters
Smpl_pts_26_19 <- Smpl_pts %>%
  terra::project(terra::crs(stack_26_19))

# Subset sample points to tile
Sub_pts_26_19 <- terra::crop(Smpl_pts_26_19, terra::ext(stack_26_19))

# Extract to tile sample points
# Takes a while (30 mins)
Ext_26_19_df <- terra::extract(stack_26_19, Sub_pts_26_19, xy = TRUE, ID = FALSE)

# Combine with uplands vector information
Ext_26_19_df_v2 <- as.data.frame(cbind(Sub_pts_26_19, Ext_26_19_df)) %>%
  # Remove the outdated coordinates columns (since we've reprojected the shapefile)
  dplyr::select(-crds_x1, -crds_x2) %>%
  # Rename the new coordinates columns as coords.x1, coords.x2
  dplyr::rename(coords.x1 = x,
                coords.x2 = y)

#glimpse((Ext_26_19_df_v2))

melt_26_19 <- Ext_26_19_df_v2 %>%
  # Pivot df longways
  tidyr::pivot_longer(cols = -c(Uplnds_, ptID, EcoType, coords.x1, coords.x2), 
                      names_to = "variable", 
                      values_to = "value") %>%
  # Convert variable column to character
  dplyr::mutate(variable = as.character(variable)) %>%
  # Create Tile, Obs_Date, Band columns
  dplyr::mutate(Tile = stringr::str_split_fixed(variable, '_', 8)[,3],
                Obs_Date = stringr::str_split_fixed(variable, '_', 8)[,4],
                Band = stringr::str_split_fixed(variable, '_', 8)[,8]) %>%
  # Get rid of unnecessary column 
  dplyr::select(-variable)

# Turn band into columns
Spec_26_19 <- melt_26_19 %>%
  tidyr::pivot_wider(names_from = "Band", values_from = "value")

Spec_26_19_clean <- Spec_26_19 %>% 
  # Remove rows where band values are NA 
  # If value is NA for one band, it will be for all
  tidyr::drop_na(B1, B2, B3, B4, B5, B7) %>%
  # Make all flagged values = NA
  # If CLOUD_QA > 0, make it NA
  dplyr::mutate(CLOUD_QA = dplyr::case_when(
    CLOUD_QA > 0 ~ NA,
    T ~ CLOUD_QA
  )) %>%
  # If PIXEL != 5440, make it NA
  # 5440 is code for totally clear pixels
  dplyr::mutate(PIXEL = dplyr::case_when(
    PIXEL != 5440 ~ NA,
    T ~ PIXEL
  )) %>%
  # If RADSET > 0, make it NA
  dplyr::mutate(RADSAT = dplyr::case_when(
    RADSAT > 0 ~ NA,
    T ~ RADSAT
  )) %>%
  # Remove rows with flagged values
  tidyr::drop_na(CLOUD_QA, PIXEL, RADSAT) %>%
  # Remove unnecessary columns
  dplyr::select(-c(RADSAT, PIXEL, LINEAGE, CLOUD_QA, ATMOS_OPACITY, QA))

# Save  
save(Ext_26_19_df_v2, Spec_26_19_clean, file = file.path(image_dir, "Master_Spec", "Tile_26_19_updated.RDATA"))

# Clean up envr before next tile
rm(list = setdiff(ls(), c("Smpl_pts", "image_dir", "LStif")))
gc()

# Tile 27_18 ------------------------------------------------------

# Group by tile
tile_27_18 <- as.list(list.files(LStif, pattern = glob2rx("*_027018*.TIF$"), full.names = TRUE))

# Identify files that are either zero-bytes, can't be read in with terra::rast(), or plot with terra::plot()
# as these files will cause the R session to abort itself later on
problem_files <- c("LE07_CU_027018_20000603_20210425_02_SR_B7.TIF",
                   "LE07_CU_027018_20000612_20210425_02_SR_B7.TIF",
                   "LE07_CU_027018_20000612_20210425_02_SR_CLOUD_QA.TIF",
                   "LE07_CU_027018_20000619_20210425_02_SR_B5.TIF",
                   "LE07_CU_027018_20000619_20210425_02_SR_B7.TIF",
                   "LE07_CU_027018_20000628_20210425_02_SR_B1.TIF",
                   "LE07_CU_027018_20010428_20210426_02_SR_B3.TIF",
                   "LE07_CU_027018_20010724_20210426_02_SR_B4.TIF",
                   "LE07_CU_027018_20031103_20210428_02_SR_B5.TIF",
                   "LE07_CU_027018_20141203_20210502_02_SR_B5.TIF",
                   "LE07_CU_027018_20161122_20210502_02_QA_PIXEL.TIF",
                   "LE07_CU_027018_20171211_20210503_02_SR_CLOUD_QA.TIF",
                   "LE07_CU_027018_20190912_20210504_02_SR_B3.TIF",
                   "LE07_CU_027018_20190912_20210504_02_SR_B4.TIF",
                   "LE07_CU_027018_20191014_20210504_02_SR_ATMOS_OPACITY.TIF",
                   "LE07_CU_027018_20200923_20210504_02_SR_B7.TIF",
                   "LE07_CU_027018_20200923_20210504_02_SR_CLOUD_QA.TIF",
                   "LE07_CU_027018_20000628_20210425_02_SR_B2.TIF",
                   "LE07_CU_027018_20000628_20210425_02_SR_B3.TIF",
                   "LE07_CU_027018_20001205_20210426_02_SR_B7.TIF",
                   "LE07_CU_027018_20030925_20210428_02_SR_B4.TIF",
                   "LE07_CU_027018_20030925_20210428_02_SR_B5.TIF",
                   "LE07_CU_027018_20031128_20210428_02_SR_B4.TIF",
                   "LE07_CU_027018_20100428_20210430_02_SR_B1.TIF",
                   "LE07_CU_027018_20100428_20210430_02_SR_B2.TIF",
                   "LE07_CU_027018_20171220_20210503_02_SR_B2.TIF",
                   "LE07_CU_027018_20200907_20210504_02_SR_B5.TIF",
                   "LE07_CU_027018_20200907_20210504_02_SR_B7.TIF",
                   "LE07_CU_027018_20200923_20210504_02_SR_B5.TIF", 
                   "LE07_CU_027018_20010615_20210426_02_SR_B7.TIF",
                   "LE07_CU_027018_20010724_20210426_02_SR_B5.TIF", 
                   "LE07_CU_027018_20031128_20210428_02_SR_B3.TIF",
                   "LE07_CU_027018_20141212_20210502_02_SR_B1.TIF",
                   "LE07_CU_027018_20161115_20210502_02_SR_B5.TIF",
                   "LE07_CU_027018_20161115_20210502_02_SR_B7.TIF",
                   "LE07_CU_027018_20171211_20210503_02_SR_B5.TIF",
                   "LE07_CU_027018_20171220_20210503_02_QA_LINEAGE.TIF",
                   "LE07_CU_027018_20200930_20210504_02_SR_B2.TIF",
                   "LE07_CU_027018_20001205_20210426_02_SR_B5.TIF",
                   "LE07_CU_027018_20010428_20210426_02_SR_B1.TIF",
                   "LE07_CU_027018_20010428_20210426_02_SR_B2.TIF",
                   "LE07_CU_027018_20010615_20210426_02_SR_B4.TIF",
                   "LE07_CU_027018_20010615_20210426_02_SR_B5.TIF",
                   "LE07_CU_027018_20031103_20210428_02_SR_B7.TIF")

tile_27_18_fixed <- tile_27_18

# Remove problem files
for (i in problem_files){
  problem_index <- stringr::str_detect(tile_27_18_fixed, i)
  tile_27_18_fixed <- tile_27_18_fixed[-which(problem_index)]
}

# Read in rasters as a list
rast_tile_27_18 <- tile_27_18_fixed %>%
  purrr::map(.f = ~terra::rast(.), .progress = T)

# Stack them
stack_27_18 <- terra::rast(rast_tile_27_18)

# Make sure sample points have the same crs as the rasters
Smpl_pts_27_18 <- Smpl_pts %>%
  terra::project(terra::crs(stack_27_18))

# Subset sample points to tile
Sub_pts_27_18 <- terra::crop(Smpl_pts_27_18, terra::ext(stack_27_18))

# Extract to tile sample points
# Takes a while (30 mins)
Ext_27_18_df <- terra::extract(stack_27_18, Sub_pts_27_18, xy = TRUE, ID = FALSE)

# Combine with uplands vector information
Ext_27_18_df_v2 <- as.data.frame(cbind(Sub_pts_27_18, Ext_27_18_df)) %>%
  # Remove the outdated coordinates columns (since we've reprojected the shapefile)
  dplyr::select(-crds_x1, -crds_x2) %>%
  # Rename the new coordinates columns as coords.x1, coords.x2
  dplyr::rename(coords.x1 = x,
                coords.x2 = y)

#glimpse((Ext_27_18_df_v2))

# Divide into parts to avoid exceeding vector memory
Ext_27_18_df_v2_A <- Ext_27_18_df_v2[1:10000,]
Ext_27_18_df_v2_B <- Ext_27_18_df_v2[10001:20000,]
Ext_27_18_df_v2_C <- Ext_27_18_df_v2[20001:30000,]
Ext_27_18_df_v2_D <- Ext_27_18_df_v2[30001:40000,]
Ext_27_18_df_v2_E <- Ext_27_18_df_v2[40001:50000,]
Ext_27_18_df_v2_F <- Ext_27_18_df_v2[50001:nrow(Ext_27_18_df_v2),]

# Combine parts in a list
Ext_27_18_df_v2_parts <- list(Ext_27_18_df_v2_A,
                              Ext_27_18_df_v2_B,
                              Ext_27_18_df_v2_C,
                              Ext_27_18_df_v2_D,
                              Ext_27_18_df_v2_E,
                              Ext_27_18_df_v2_F)

# Create empty list to store the result for each part
Spec_27_18_clean_together <- list()

for (i in 1:length(Ext_27_18_df_v2_parts)){
  
  message("On part: ", i, " of ", length(Ext_27_18_df_v2_parts))
  
  melt_27_18 <- Ext_27_18_df_v2_parts[[i]] %>%
    # Pivot df longways
    tidyr::pivot_longer(cols = -c(Uplnds_, ptID, EcoType, coords.x1, coords.x2), 
                        names_to = "variable", 
                        values_to = "value") %>%
    # Convert variable column to character
    dplyr::mutate(variable = as.character(variable)) %>%
    # Create Tile, Obs_Date, Band columns
    dplyr::mutate(Tile = stringr::str_split_fixed(variable, '_', 8)[,3],
                  Obs_Date = stringr::str_split_fixed(variable, '_', 8)[,4],
                  Band = stringr::str_split_fixed(variable, '_', 8)[,8]) %>%
    # Get rid of unnecessary column 
    dplyr::select(-variable)
  
  # Turn band into columns
  Spec_27_18 <- melt_27_18 %>%
    tidyr::pivot_wider(names_from = "Band", values_from = "value")
  
  Spec_27_18_clean <- Spec_27_18 %>% 
    # Remove rows where band values are NA 
    # If value is NA for one band, it will be for all
    tidyr::drop_na(B1, B2, B3, B4, B5, B7) %>%
    # Make all flagged values = NA
    # If CLOUD_QA > 0, make it NA
    dplyr::mutate(CLOUD_QA = dplyr::case_when(
      CLOUD_QA > 0 ~ NA,
      T ~ CLOUD_QA
    )) %>%
    # If PIXEL != 5440, make it NA
    # 5440 is code for totally clear pixels
    dplyr::mutate(PIXEL = dplyr::case_when(
      PIXEL != 5440 ~ NA,
      T ~ PIXEL
    )) %>%
    # If RADSET > 0, make it NA
    dplyr::mutate(RADSAT = dplyr::case_when(
      RADSAT > 0 ~ NA,
      T ~ RADSAT
    )) %>%
    # Remove rows with flagged values
    tidyr::drop_na(CLOUD_QA, PIXEL, RADSAT) %>%
    # Remove unnecessary columns
    dplyr::select(-c(RADSAT, PIXEL, LINEAGE, CLOUD_QA, ATMOS_OPACITY))
  
  # Store result for this part
  Spec_27_18_clean_together[[i]] <- Spec_27_18_clean
}

# Turn into a dataframe
Spec_27_18_clean_together_df <- Spec_27_18_clean_together %>%
  purrr::map_dfr(.f = select, everything()) %>%
  dplyr::relocate(B2, .after = B1) %>%
  dplyr::relocate(B3, .after = B2) %>%
  dplyr::relocate(B4, .after = B3) %>%
  dplyr::relocate(B5, .after = B4) %>%
  dplyr::relocate(B7, .after = B5)

# Save  
save(Ext_27_18_df_v2, Spec_27_18_clean_together_df, file = file.path(image_dir, "Master_Spec", "Tile_27_18_updated.RDATA"))

# Clean up envr before next tile
rm(list = setdiff(ls(), c("Smpl_pts", "image_dir", "LStif")))
gc()

# Tile 27_19 ------------------------------------------------------

# Group by tile
tile_27_19 <- as.list(list.files(LStif, pattern = glob2rx("*_027019*.TIF$"), full.names = TRUE))

# Identify files that are either zero-bytes, can't be read in with terra::rast(), or plot with terra::plot()
# as these files will cause the R session to abort itself later on
problem_files <- c("LE07_CU_027019_20000315_20210425_02_QA_RADSAT.TIF",
                   "LE07_CU_027019_20000527_20210425_02_SR_B7.TIF",
                   "LE07_CU_027019_20000603_20210425_02_SR_B2.TIF",
                   "LE07_CU_027019_20020720_20210427_02_SR_B5.TIF",
                   "LE07_CU_027019_20020720_20210427_02_SR_B7.TIF",
                   "LE07_CU_027019_20020727_20210427_02_QA_LINEAGE.TIF",
                   "LE07_CU_027019_20020812_20210427_02_QA_LINEAGE.TIF",
                   "LE07_CU_027019_20020812_20210427_02_SR_B2.TIF",
                   "LE07_CU_027019_20061127_20210429_02_QA_LINEAGE.TIF",
                   "LE07_CU_027019_20061127_20210429_02_QA_PIXEL.TIF",
                   "LE07_CU_027019_20080618_20210430_02_QA_LINEAGE.TIF",
                   "LE07_CU_027019_20090204_20210430_02_SR_B3.TIF",
                   "LE07_CU_027019_20090204_20210430_02_SR_B4.TIF",
                   "LE07_CU_027019_20101021_20210430_02_QA_LINEAGE.TIF",
                   "LE07_CU_027019_20101021_20210430_02_QA_PIXEL.TIF",
                   "LE07_CU_027019_20131130_20210501_02_QA_RADSAT.TIF",
                   "LE07_CU_027019_20201228_20210504_02_SR_B1.TIF",
                   "LE07_CU_027019_20020821_20210427_02_SR_B7.TIF",
                   "LE07_CU_027019_20020821_20210427_02_SR_CLOUD_QA.TIF",
                   "LE07_CU_027019_20131123_20210501_02_SR_B7.TIF", #
                   "LE07_CU_027019_20000527_20210425_02_SR_CLOUD_QA.TIF",
                   "LE07_CU_027019_20020727_20210427_02_SR_B4.TIF",
                   "LE07_CU_027019_20020922_20210427_02_SR_B1.TIF",
                   "LE07_CU_027019_20020922_20210427_02_SR_B2.TIF",
                   "LE07_CU_027019_20131123_20210501_02_SR_B5.TIF",
                   "LE07_CU_027019_20140705_20210501_02_SR_B2.TIF",
                   "LE07_CU_027019_20140705_20210501_02_SR_B3.TIF",
                   "LE07_CU_027019_20140705_20210501_02_SR_B4.TIF",
                   "LE07_CU_027019_20140705_20210501_02_SR_B7.TIF",
                   "LE07_CU_027019_20201228_20210504_02_SR_B2.TIF")

tile_27_19_fixed <- tile_27_19

# Remove problem files
for (i in problem_files){
  problem_index <- stringr::str_detect(tile_27_19_fixed, i)
  tile_27_19_fixed <- tile_27_19_fixed[-which(problem_index)]
}

# Read in rasters as a list
rast_tile_27_19 <- tile_27_19_fixed %>%
  purrr::map(.f = ~terra::rast(.), .progress = T)

# Stack them
stack_27_19 <- terra::rast(rast_tile_27_19)

# Make sure sample points have the same crs as the rasters
Smpl_pts_27_19 <- Smpl_pts %>%
  terra::project(terra::crs(stack_27_19))

# Subset sample points to tile
Sub_pts_27_19 <- terra::crop(Smpl_pts_27_19, terra::ext(stack_27_19))

# Identify indices to split chunks on
chunks <- seq(1, length(Sub_pts_27_19), by = 10000)

# For each chunk...
for (i in 1:length(chunks)){
  
  message("On part: ", i, " of ", length(chunks))
  
  # Amount to increment the index by
  increment <- 9999
  
  # If we're on the last chunk, increment by 339831-330001=9830 instead
  if (i == length(chunks)) {
    increment <- length(Sub_pts_27_19)-chunks[i]
  }
  
  # Split into chunk
  Sub_pts_27_19_chunk <- Sub_pts_27_19[chunks[i]:(chunks[i]+increment)]
  
  # Extract to tile sample points
  # Takes a while (20 mins)
  Ext_27_19_df <- terra::extract(stack_27_19, Sub_pts_27_19_chunk, xy = TRUE, ID = FALSE)
  
  # Combine with uplands vector information
  Ext_27_19_df_v2 <- as.data.frame(cbind(Sub_pts_27_19_chunk, Ext_27_19_df)) %>%
    # Remove the outdated coordinates columns (since we've reprojected the shapefile)
    dplyr::select(-crds_x1, -crds_x2) %>%
    # Rename the new coordinates columns as coords.x1, coords.x2
    dplyr::rename(coords.x1 = x,
                  coords.x2 = y)
  
  # Create a file name for this chunk
  file_name1 <- paste0("Ext_27_19_df_", chunks[i], "_", chunks[i]+increment, ".csv")
  
  # Export extracted sample points
  write.csv(Ext_27_19_df_v2, 
            file = file.path(image_dir, "Master_Spec", "Tile_27_19_parts", file_name1), 
            row.names = F)
  
  #glimpse((Ext_27_19_df_v2))
  
  melt_27_19 <- Ext_27_19_df_v2 %>%
    # Pivot df longways
    tidyr::pivot_longer(cols = -c(Uplnds_, ptID, EcoType, coords.x1, coords.x2), 
                        names_to = "variable", 
                        values_to = "value") %>%
    # Convert variable column to character
    dplyr::mutate(variable = as.character(variable)) %>%
    # Create Tile, Obs_Date, Band columns
    dplyr::mutate(Tile = stringr::str_split_fixed(variable, '_', 8)[,3],
                  Obs_Date = stringr::str_split_fixed(variable, '_', 8)[,4],
                  Band = stringr::str_split_fixed(variable, '_', 8)[,8]) %>%
    # Get rid of unnecessary column 
    dplyr::select(-variable)
  
  # Turn band into columns
  Spec_27_19 <- melt_27_19 %>%
    tidyr::pivot_wider(names_from = "Band", values_from = "value")
  
  Spec_27_19_clean <- Spec_27_19 %>% 
    # Remove rows where band values are NA 
    # If value is NA for one band, it will be for all
    tidyr::drop_na(B1, B2, B3, B4, B5, B7) %>%
    # Make all flagged values = NA
    # If CLOUD_QA > 0, make it NA
    dplyr::mutate(CLOUD_QA = dplyr::case_when(
      CLOUD_QA > 0 ~ NA,
      T ~ CLOUD_QA
    )) %>%
    # If PIXEL != 5440, make it NA
    # 5440 is code for totally clear pixels
    dplyr::mutate(PIXEL = dplyr::case_when(
      PIXEL != 5440 ~ NA,
      T ~ PIXEL
    )) %>%
    # If RADSET > 0, make it NA
    dplyr::mutate(RADSAT = dplyr::case_when(
      RADSAT > 0 ~ NA,
      T ~ RADSAT
    )) %>%
    # Remove rows with flagged values
    tidyr::drop_na(CLOUD_QA, PIXEL, RADSAT) %>%
    # Remove unnecessary columns
    dplyr::select(-c(RADSAT, PIXEL, LINEAGE, CLOUD_QA, ATMOS_OPACITY))
  
  # Create another file name for this chunk
  file_name2 <- paste0("Spec_27_19_clean_", chunks[i], "_", chunks[i]+increment, ".csv")
  
  # Export clean extracted sample points
  write.csv(Spec_27_19_clean, 
            file = file.path(image_dir, "Master_Spec", "Tile_27_19_parts", file_name2), 
            row.names = F)
  
  rm(Sub_pts_27_19_chunk, Ext_27_19_df, Ext_27_19_df_v2, melt_27_19, Spec_27_19, Spec_27_19_clean)
  
}

# Point to the folder containing all the parts
parts <- file.path(image_dir, "Master_Spec", "Tile_27_19_parts")

# Identify all Ext_27_19_df_ files
Ext_27_19_together <- as.list(list.files(parts, pattern = "Ext_27_19_df_", full.names = TRUE))

# Combine files into a dataframe
Ext_27_19_together_df <- Ext_27_19_together %>%
  purrr::map(.f = read.csv) %>%
  purrr::map_dfr(.f = select, everything())

# Identify all Spec_27_19_clean_ files
Spec_27_19_clean_together <- as.list(list.files(parts, pattern = "Spec_27_19_clean_", full.names = TRUE))

# Combine files into a dataframe
Spec_27_19_clean_together_df <- Spec_27_19_clean_together %>%
  purrr::map(.f = read.csv) %>%
  purrr::map_dfr(.f = select, everything())

# Save  
save(Ext_27_19_together_df, Spec_27_19_clean_together_df, file = file.path(image_dir, "Master_Spec", "Tile_27_19_updated.RDATA"))

# Clean up envr
rm(list = setdiff(ls(), c("Smpl_pts", "image_dir", "LStif")))
gc()

## --------------------------------------------- ##
#           Dataframe Integration -----
## --------------------------------------------- ##

# Load files
load(file.path(image_dir, "Master_Spec", "Tile_26_18_updated.RDATA"))
load(file.path(image_dir, "Master_Spec", "Tile_26_19_updated.RDATA"))
load(file.path(image_dir, "Master_Spec", "Tile_27_18_updated.RDATA"))
load(file.path(image_dir, "Master_Spec", "Tile_27_19_updated.RDATA"))

Spec_27_19_clean_together_df <- Spec_27_19_clean_together_df %>%
  # Convert some columns to character first
  dplyr::mutate(Tile = as.character(Tile),
                Obs_Date = as.character(Obs_Date))

# List tiles together
Specs <- list(Spec_26_18_clean, Spec_26_19_clean, Spec_27_18_clean_together_df, Spec_27_19_clean_together_df)

# Calculate Spectral Indices
# Band combinations guide (https://www.esri.com/arcgis-blog/products/product/imagery/band-combinations-for-landsat-8/)
# Turn off factors
options(stringsAsFactors = FALSE)

# Bind tiles together
Spec_Master <- Specs %>%
  purrr::list_rbind(x = .) %>%
  # Calculate NDVI using the red (band 3) and nir (band 4) bands  
  # NVDI = (NIR - Red) / (NIR + Red) 
  dplyr::mutate(NDVI = (B4-B3)/(B4+B3),
                # using nir (band 4) and swir (band 7) *there is no band 6 so band 7 is in the 6th position
                # NBR = (NIR - SWIR) / (NIR + SWIR)
                NBR = (B4-B7)/(B4+B7),
                # using swir1 (band 5) and swir2 (band 7) *supposedly good for measuring recovery
                # NBR2 = (SWIR1 - SWIR2) / (SWIR1 + SWIR2)
                NBR2 = (B5-B7)/(B5+B7))

# Check it
head(Spec_Master)

# Save spectral master dataframe
save(Spec_Master, file = file.path(image_dir, "Master_Spec", "Spec_Master_updated.RDATA"))

# Export as csv
write.csv(Spec_Master, 
          file = file.path(image_dir, "Master_Spec", "Spec_Master_updated.csv"), 
          row.names = F)

