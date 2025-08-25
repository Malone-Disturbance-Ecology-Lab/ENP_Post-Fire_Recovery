## --------------------------------------------- ##
#       Landsat Image Processing (Part A)
## --------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script does the following steps to process images for 4 ARD Landsat tiles covering the Everglades

## Images were downloaded manually from Earth Explorer: https://earthexplorer.usgs.gov/

## 1. open tar files and saves Landsat .tif image (RUN ONCE)

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

rm(list=ls())

# Point to the image processing folder on server
image_dir <- file.path("/", "Volumes", "MaloneLab", "Research", "ENP", "ENP Fire", "Grace_McLeod", "Image_Processing")

## --------------------------------------------- ##
#       Opening TAR Files ***RUN ONCE*** -----
## --------------------------------------------- ##

# Bulk download from Earth Explorer
# Download completed in chunks. Run on each chunk.
bulk <- file.path(image_dir, "BulkDownload", "BulkDownLoad_2000_2005")
#bulk <- file.path(image_dir, "BulkDownload", "BulkDownLoad_2006_2009")
#bulk <- file.path(image_dir, "BulkDownload", "BulkDownLoad_2010_2012")
#bulk <- file.path(image_dir, "BulkDownload", "BulkDownLoad_2013_2015")
#bulk <- file.path(image_dir, "BulkDownload", "BulkDownLoad_2016_2020")

# Create path to store tif files
LStif <- file.path(image_dir, "LS_tif")

# Select .tar files
# Make a list of all the .tar files
all_tars <- list.files(bulk, pattern = "tar") 
# Just the surface reflectance (SR) ones
SR.tar <- list.files(bulk, pattern = glob2rx("*_SR*.tar$"), full.names = TRUE)

# Loop to untar files
for (tar in SR.tar){
  print(tar)
  LStifs <- untar(tar, list = FALSE, exdir = LStif) 
}