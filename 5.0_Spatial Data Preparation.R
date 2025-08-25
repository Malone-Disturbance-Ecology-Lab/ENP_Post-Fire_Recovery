# Spatial Data Preparation
# Creates layers needed in the map figures.
rm(list=ls())
# EVG boundary
EVG_bound <- read_sf(dsn = "/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/AOI", layer = "EVG_bound")
#st_crs(EVG_bound) <- 26917 # assign crs (NAD83 UTM17)
EVG_bound <- st_transform(EVG_bound, crs=4326) # transform to lat long

# Sample points for recovery
# use OG sample pts to fix coordinates in drivers df
smpl_pts <- read_sf(dsn = "/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Sampling", layer = "Recov_smpl_pts") %>% st_transform(crs=4326)

smpl_pts <- smpl_pts %>%
  select(ptID, geometry)
smpl_pts <- as.data.frame(smpl_pts)

load( file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/DriverAnalysis_DF_082025.RDATA")


drivers.sp <- left_join(driver.analysis2, smpl_pts, by="ptID")
drivers.sp <- st_as_sf(drivers.sp) %>%  st_transform(drivers.sp, crs=4326) # transform to lat long

# Upland sample points
# all pinelands, not just recovery pts
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Veg_layers/Uplands_all.RDATA")
upland_pts <- st_as_sf(Uplands)
st_crs(upland_pts)
upland_pts <- st_transform(upland_pts, crs=4326) 

# Florida Boundary
FL_bound <- read_sf(dsn = "/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/AOI/FL_boundary", 
                    layer = "Florida_State_Boundary")
st_crs(FL_bound)
FL_bound <- st_transform(FL_bound, crs=4326) # transform to lat long
# Extract the bounding box of `bound`
bbox_bound <- st_bbox(EVG_bound)
# Increase the bounding box by a certain factor, for example 0.01 degrees
# You can adjust this factor based on how much larger you want the bounding box to be
buffer_factor <- 0.05
# Create the new bounding box with adjusted coordinates
new_bbox <- st_sfc(
  st_polygon(list(matrix(c(
    bbox_bound["xmin"] - buffer_factor, bbox_bound["ymin"] - buffer_factor,  # bottom-left corner
    bbox_bound["xmin"] - buffer_factor, bbox_bound["ymax"] + buffer_factor,  # top-left corner
    bbox_bound["xmax"] + buffer_factor, bbox_bound["ymax"] + buffer_factor,  # top-right corner
    bbox_bound["xmax"] + buffer_factor, bbox_bound["ymin"] - buffer_factor,  # bottom-right corner
    bbox_bound["xmin"] - buffer_factor, bbox_bound["ymin"] - buffer_factor   # closing the polygon
  ), ncol = 2, byrow = TRUE))),
  crs = st_crs(EVG_bound)  # Use the CRS of `bound`
)
# clip
FL_clip <- st_intersection(FL_bound, new_bbox)
plot(FL_clip$geometry)

# FL inlay
inlay_plot <- ggplot() +
  geom_sf(data = FL_bound, fill = "lightgrey", color = "#3b3a3f") +
  geom_sf(data = EVG_bound, color = "#3b3a3f", fill = "#3b3a3f", size = 1.5) +    
  theme_classic() +
  theme(legend.position = "none", 
        axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.line = element_blank(),  
        axis.ticks = element_blank())  

# Plot the new bounding box
plot(new_bbox, col = "lightgray", border = "red", main = "Extended Bounding Box")
plot(FL_bound, add = TRUE, col = "blue", border = "black")  
plot(upland_pts, add=T, col= "black")

# Hexagonal Grid
# Get the bounding box of the data
hex_bbox <- st_bbox(drivers.sp)
if (!exists("hex_bbox") || !inherits(hex_bbox, "bbox")) {
  hex_bbox <- st_bbox(drivers.sp)
}
# Calculate the cell size to generate desired number of hexagons
cell_size <- sqrt(
  as.numeric(hex_bbox["xmax"] - hex_bbox["xmin"]) * 
    as.numeric(hex_bbox["ymax"] - hex_bbox["ymin"]) / 1000)
# Create a hexagonal grid
hex_grid <- st_make_grid(drivers.sp, cellsize = cell_size, square = FALSE)
hex_grid_sf <- st_sf(geometry = hex_grid)


# Save spatial files
save(EVG_bound, FL_bound, FL_clip, inlay_plot, new_bbox, 
     smpl_pts, upland_pts, drivers.sp, hex_grid_sf,
     file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Spatial_files.RDATA")
