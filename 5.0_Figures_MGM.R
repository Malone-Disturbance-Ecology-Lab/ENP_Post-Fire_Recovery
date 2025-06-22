# FIGURES
# GRACE MCLEOD, 2024

# This script generates figures for the Everglades Post-Fire Recovery Manuscript
# using common styles and color pallets.


library(ggplot2)
library(tidyverse)
library(ggpubr)
library(randomForest)
library(ggridges)
library(sf)
library(sp)
library(terra)
library(leaflet)
library(ggspatial)
library(ggmap)
library(MetBrewer)
library(grid)
library(ggpubr)
library(corrplot)
library(ggcorrplot)

rm(list=ls())

# Spatial Data 
load('/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Manuscript/Spatial_files.RDATA')
# Climate Data
load('/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Climate/Annual_Climate_Summary_ENP.RDATA')
#load("~/Documents/Thesis/RECOVERY/ENP_Post-Fire_Recovery/.RData")
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Manuscript/DriversData.RDATA")

### DATA PREP : ***Run Once*** ##############################################################################

# Spatial Fuckery.........................................................................

# EVG boundary
EVG_bound <- read_sf(dsn = "/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/AOI", layer = "EVG_bound")
#st_crs(EVG_bound) <- 26917 # assign crs (NAD83 UTM17)
EVG_bound <- st_transform(EVG_bound, crs=4326) # transform to lat long

# Sample points for recovery
# use OG sample pts to fix coordinates in drivers df
smpl_pts <- read_sf(dsn = "/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Sampling", layer = "Recov_smpl_pts")
st_crs(smpl_pts) <- 32617 # assign crs (WGS83 UTM17)
smpl_pts <- smpl_pts %>%
  select(ptID, geometry)
smpl_pts <- as.data.frame(smpl_pts)
drivers.sp <- left_join(driver.analysis, smpl_pts, by="ptID")
drivers.sp <- st_as_sf(drivers.sp)
st_crs(drivers.sp) <- 32617
drivers.sp <- st_transform(drivers.sp, crs=4326) # transform to lat long

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
     file="Spatial_files.RDATA")


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


### FIGURE 1: Pinelands Map ##################################################################

unique(upland_pts$L4_name)
class(upland_pts$L4_name)
upland_pts$L4_name[upland_pts$L4_name == "Pine Upland"] <- "Pine Flatwoods"

# plot
map.pinelands <- upland_pts %>%
  filter(L4_name %in% c("Pine Rockland", "Pine Flatwoods"))

map.pinelands_plot <- ggplot(map.pinelands) + 
  geom_sf(data = FL_clip, fill = "lightgrey", color = "lightgrey") + 
  geom_sf(data = EVG_bound, fill = NA, color = "black", size = 1.5) +    
  geom_sf(aes(colour = L4_name, fill = L4_name), size = 5) +
  theme_classic() +
  scale_color_manual(values = c("#777777", "#555555")) +
  scale_fill_manual(values = c("#777777", "#555555")) +
  theme(text = element_text(size = 35),
        legend.position = "bottom") +
  labs(fill = "Ecosystem Type", color = NULL, x=NULL, y=NULL) +
  guides(color = "none") +
  scale_x_continuous(
    breaks = seq(-81.5, -80, by = 0.5)) + 
  scale_y_continuous(
    breaks = seq(24, 27, by = 0.5))+ 
  # park lables
  #annotate("text", x = -80.6, y = 25.8, label = "Everglades National Park", size = 6.5, color = "black", fontface ="plain" ) +
 # annotate("text", x = -81.1, y = 26.3, label = "Big Cypress National Preserve", size = 6.5, color = "black", fontface = "plain") + 
  # FL inlay
  annotation_custom( 
    grob = ggplotGrob(inlay_plot), 
    xmin = -81.6, xmax = -81.2, ymin = 24.8, ymax = 25.2) +
  # north arrow 
  annotation_north_arrow(
    location = "tr", 
    width = unit(1, "cm"), 
    height = unit(1.5, "cm"),
    style = north_arrow_orienteering()) +
  # scale bar
  annotation_scale(
    location = "bl",
    width_hint = 0.2, 
    height = unit(0.4, "cm"), 
    style = "ticks", 
    bar_cols = c("black", "white"),
    text_cex = 2)


map.pinelands_plot_new <- ggplot(map.pinelands) + 
  geom_sf(data = FL_clip, fill = "lightgrey", color = "lightgrey") + 
  geom_sf(data = EVG_bound, fill = NA, color = "black", size = 1.5) +    
  geom_sf(aes(colour = L4_name, fill = L4_name), size = 5) +
  theme_classic() +
  scale_color_manual(values = c("black", "black")) +
  scale_fill_manual(values = c("black", "black")) +
  theme(text = element_text(size = 35),
        legend.position = "none") +
  labs(fill = "Ecosystem Type", color = NULL, x=NULL, y=NULL) +
  guides(color = "none") +
  scale_x_continuous(
    breaks = seq(-81.5, -80, by = 0.5)) + 
  scale_y_continuous(
    breaks = seq(24, 27, by = 0.5))+ 
  # park lables
  #annotate("text", x = -80.6, y = 25.8, label = "Everglades National Park", size = 6.5, color = "black", fontface ="plain" ) +
  # annotate("text", x = -81.1, y = 26.3, label = "Big Cypress National Preserve", size = 6.5, color = "black", fontface = "plain") + 
  # FL inlay
  annotation_custom( 
    grob = ggplotGrob(inlay_plot), 
    xmin = -81.6, xmax = -81.2, ymin = 24.8, ymax = 25.2) +
  # north arrow 
  annotation_north_arrow(
    location = "tr", 
    width = unit(1, "cm"), 
    height = unit(1.5, "cm"),
    style = north_arrow_orienteering()) +
  # scale bar
  annotation_scale(
    location = "bl",
    width_hint = 0.2, 
    height = unit(0.4, "cm"), 
    style = "ticks", 
    bar_cols = c("black", "white"),
    text_cex = 2)


# Save
setwd("Figures")
ggsave("Map_Pinelands_06212025.png", plot = map.pinelands_plot_new, width = 12, height = 10, dpi = 300)


### FIGURE 2: Historical Conditions- Climate Normals and Fire History  #################################################################

# Climate Normals ..........................................................................................
summary.tot$Year <- summary.tot$Year %>% as.numeric()

plot.temp <- ggplot(data =summary.tot) + 
  geom_line(aes(x = Year, y = TMAX.mean), linetype = "dashed", size = 2, alpha = 0.8) + 
  geom_line(aes(x = Year, y = TMIN.mean), linetype = "dotted", size = 2, alpha = 0.8) +
  #geom_point( aes(x=Year, y = TMAX.mean), col="#e76155") + 
  #geom_line( aes(x=Year, y = TMAX.mean), col="#e76155", size=2, alpha=.5) + 
  #geom_point( aes(x=Year, y = TMIN.mean), col="#376795") + 
  #geom_line( aes(x=Year, y = TMIN.mean), col="#376795", size=2, alpha=.5) +
  labs(x="", y="Temperature (C)", tag="B") + 
  ylim(15, 40) +
  theme_bw() + 
  theme(text = element_text(size = 25), 
        legend.position = "top") +
  annotate("text", x = 2010, y = 32, label = "Maximum", size = 8, fontface ="plain" ) +
  annotate("text", x = 2010, y = 21, label = "Minimum", size = 8,  fontface = "plain") 

plot.prcp <- ggplot(data =summary.tot) + 
  geom_line( aes(x=Year, y = PRCP), size=3, alpha=.5) +
  geom_point( aes(x=Year, y = PRCP)) + 
  labs(x="", y="Precipitation (mm)", tag="C") +
  ylim(1200, 2050) +
  theme_bw() + 
  theme(text = element_text(size = 25))

plot.pdsi <- ggplot(data =summary.tot) + 
  geom_line( aes(x=Year, y = PDSI.mean), size=3, alpha=.5) +
  geom_point( aes(x=Year, y = PDSI.mean)) + 
  labs(x="", y="PDSI", tag="D") +
  ylim(-4, 4) +
  theme_bw() + 
  theme(text = element_text(size = 25))

# combine 
plot_climate <- ggarrange(plot.temp,
                          plot.prcp,
                          plot.pdsi,
                          nrow=3, ncol=1,
                          #labels=c('B','C', 'D'),
                          font.label = list(size = 30))

# Fire History ..........................................................................................

# Map: fire frequency
# Spatial join to assign points to hexagons
hex_aggregated <- st_join(hex_grid_sf, drivers.sp, join = st_intersects) %>%
  group_by(geometry) %>%
  summarise(
    avg_TotalFires = mean(TotalFires, na.rm = TRUE),
    avg_PrevInt = mean(Prev.Int, na.rm = TRUE), 
    .groups = "drop")
# visualize
map_FireFreq <- hex_aggregated %>%
  arrange(avg_TotalFires) %>%  
  ggplot() + 
  geom_sf(data = FL_clip, fill = "lightgrey", color = "lightgrey") + 
  geom_sf(data = EVG_bound, fill = NA, color = "black", size = 1.5) +    
  geom_sf(
    aes(fill = avg_TotalFires), color = NA, size = 0.3) +  
  theme_classic() +
  scale_fill_gradientn(                                    
    colors = met.brewer("Tam", n = 6, direction = -1),  
    na.value = "transparent",                             
    breaks = c(1, 3, 5, 7, 9, 11),
    labels = c("1", "3", "5", "7", "9", "11")) +  
  scale_x_continuous(
    breaks = seq(-81.5, -80, by = 0.5)) + 
  scale_y_continuous(
    breaks = seq(24, 27, by = 0.5))+ 
  theme(
    text = element_text(size = 30),
    legend.position = "bottom") +
  labs(fill = "Total Fires",  
       tag = "A") +
  # FL inlay
  annotation_custom( 
    grob = ggplotGrob(inlay_plot), 
    xmin = -81.6, xmax = -81.2, ymin = 24.8, ymax = 25.2) +
  # north arrow 
  annotation_north_arrow(
    location = "tr", 
    width = unit(1.5, "cm"), 
    height = unit(2, "cm"),
    style = north_arrow_orienteering()) +
  # scale bar
  annotation_scale(
    location = "bl",
    width_hint = 0.2, 
    height = unit(0.4, "cm"), 
    style = "ticks", 
    bar_cols = c("black", "white"),
    text_cex = 2)

#map.plot <- ggarrange(map.FireFreq, labels= c("A"))
#final.plot.HistCond <- ggarrange(ggarrange(map.plot, labels="A"), climate_plots , nrow=1)
#final.plot.HistCond <- ggarrange(ggarrange(map.FireFreq), climate_plots , nrow=1)

# Final Figures
final.plot.HistCond <- ggarrange(map_FireFreq,
                                 plot_climate,
                                 nrow=1, ncol=2)


# Save
#getwd()
#setwd("/Users/gracemcleod/Documents/Thesis/RECOVERY/ENP_Post-Fire_Recovery/Figures")
ggsave("Historical_Conditions.png", plot = final.plot.HistCond, width = 18, height =12 , dpi = 300)



### FIGURE 3: Post-fire, Pre-fire, NDVI Enhanced #############################################################################

# Panel comparing:
# A: modeled vs observed NDVI
# B: wet vs dry season model values
# C: pre-fire vs post-fire vs baseline

# Data Prep: #####
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/NDVI_rf.RDATA")
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/BL_train_test.RDATA")
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/Recov_Combo2.RDATA")

# Observed versus Predicted
test$pred <- predict(NDVI_rf, newdata= test)
summary.lm.val <- summary(lm(data=test, pred~ NDVI))
summary.lm.val$coefficients
summary.lm.val$r.squared

plot_BL_obsVSpred <- ggplot( data= test,aes(x= NDVI, y=pred)) + 
  geom_point() +
  geom_smooth(method='lm', se=TRUE, col='grey71') + 
  labs(x="Observed", y="Predicted")+ 
  ylim(0.2, 0.5) + xlim(0.2, 0.5) + 
  geom_text( x=0.3, y=0.5, label="Y = 0.71*Observed + 0.082", size=6)+ 
  theme(text = element_text(size=20))


# RUN TO SEE WHAT ALL THIS IS ABOUT! 

# Subset the data needed:
Recovery.baseline <- Recov_Combo %>% select( ptID, StartDate, Obs_Date, model.NDVI, NDVI )
Recovery.baseline$Fire.date <- as.Date(Recovery.baseline$StartDate)
Recovery.baseline$NDVI.date <- as.Date(Recovery.baseline$Obs_Date)
Recovery.baseline$Time.diff <- as.numeric(Recovery.baseline$NDVI.date - Recovery.baseline$Fire.date )
summary(as.numeric(Recovery.baseline$Time.diff)) 
# There are no negative values. remove values within 
# Subset to the min date difference and measure the gap between model and observed ndvi:
Recovery.baseline.minimum <- Recovery.baseline %>% group_by(ptID) %>% summarise( Time.diff.min = min(Time.diff)) %>% full_join(Recovery.baseline, by = 'ptID')
Recovery.baseline.start <- Recovery.baseline.minimum %>%  filter( Time.diff.min == Time.diff)
Recovery.baseline.start$NDVI.Diff <- Recovery.baseline.start$model.NDVI - Recovery.baseline.start$NDVI
# For the MS we need to show how the baseline compares to NDVI right before the fire: 
time.diff.90.df <- Recovery.baseline.start %>% filter( Time.diff <= 90) 
mean( Recovery.baseline.start$model.NDVI) - mean(Recovery.baseline.start$NDVI)
mean( Recov_Combo$PreNDVI ) - mean(Recovery.baseline.start$NDVI)/ mean( Recovery.baseline.start$model.NDVI) - mean(Recovery.baseline.start$NDVI)*100
range(time.diff.90.df$NDVI)
mean(time.diff.90.df$NDVI)
mean(time.diff.90.df$model.NDVI) - mean(time.diff.90.df$NDVI)
mean(Recov_Combo$PreNDVI)- mean(time.diff.90.df$NDVI)
var(Recov_Combo$PreNDVI)/sqrt(length(Recov_Combo$PreNDVI))
time.diff.90.df$Month <- format( as.Date(time.diff.90.df $Obs_Date), '%m')
Recovery.baseline.start$Month <- format( as.Date(Recovery.baseline.start$Obs_Date), '%m')
Recovery.baseline.start.1 <- Recovery.baseline.start %>% filter(Month == "01")
Recovery.baseline.start.6 <- Recovery.baseline.start %>% filter(Month == "06")

# A:
# Figure un-burned (Train) density broken up into quantiles and baseline density for initial values.
train.predtions <- train
# set values for baseline
train.predtions$TotalFires <- 9
train.predtions$Prev.Int <- 6
train.predtions$Pt.B4max <- 21526.95
train.predtions$Pt.B4min <- 10637.73
train.predtions$Pt.B5max <- 15097.51
train.predtions$SWIR1.SWIR2 <- 1.259455
train.predtions$NIR.SWIR1 <- 1.429776
train.predtions$model.NDVI <-predict(NDVI_rf, data =train.predtions )

save(train.predtions, time.diff.90.df, Recovery.baseline.start,Recovery.baseline.start.1,
     Recovery.baseline.start.6, file = "/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Figure3_Data.RDATA")

# A: Modeled VS Observed ####

load(file = "/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Figure3_Data.RDATA" )
train.NDVI.plot <- ggplot(data =train.predtions) + 
  geom_density(aes(x= NDVI, y = ..scaled..), color="black" , size=1.5) + 
  geom_density(aes(x= model.NDVI, y = ..scaled..),linetype = "dashed", color="#1e5b24", size=1.5 ) + 
  annotate( geom="text", x= 0.33, y=1, label="Modeled", size=5 , fontface = "bold",color="#1e5b24") +
  annotate( geom="text", x= 0.35, y=0.8, label="Observed", size=5 , fontface = "bold",color="black") + 
  labs(y="Density", x=NULL, tag="B") +
  theme_bw() +
  theme(text = element_text(size = 20))


# B: Wet vs Dry Season NDVI Predictions
# THE BASELINE MODEL PREDICTIONS ARE MULTI-MODAL due to month/season
baseline.season.plot <- ggplot() + 
  geom_density(data =Recovery.baseline.start, aes(x= model.NDVI, y = ..scaled..), fill="#1e5b24", alpha=0.2 ,linetype = "dashed" , color="#1e5b24") + 
  geom_density(data= Recovery.baseline.start.1, aes(x= model.NDVI, y = ..scaled..),fill="#7e5421", alpha=0.5 ,linetype = "dashed" , color="#7e5421") +
  geom_density(data= Recovery.baseline.start.6, aes(x= model.NDVI, y = ..scaled..),fill="#2d77aa", alpha=0.5 ,linetype = "dashed" , color="#2d77aa")  + 
  labs(y="Density", x=NULL, tag="C") +
  annotate( geom="text", x= 0.28, y=1.3, label="Dry Season", size=5 , fontface = "bold", color="#7e5421") +
  annotate( geom="text", x= 0.32, y=1.15, label="Modeled", size=5 , fontface = "bold", color="#1e5b24")+
  annotate( geom="text", x= 0.35, y=1.3, label="Wet Season", size=5 , fontface = "bold",color="#2d77aa") +
  theme_bw() +
  theme(text = element_text(size = 20)) 


# C: Pre vs Post vs Baseline
recovery.density.plot <- ggplot(data =time.diff.90.df) + 
  geom_density(aes(x= NDVI, y = ..scaled..),fill="#a01202", alpha=0.2 ,linetype = "solid" , color="#a01202") + 
  geom_density(aes(x= model.NDVI, y = ..scaled..), fill="#1e5b24", alpha=0.4,linetype = "dashed" , color="#1e5b24") + 
  geom_density( data= Recov_Combo, aes( x= PreNDVI, y = ..scaled..), fill="goldenrod", alpha=0.2 ,linetype = "solid" , color="goldenrod") + 
  labs(y="Density", x="NDVI", tag="D") +
  annotate( geom="text", x= 0.1, y=1.13, label="Post-fire", size=5 , fontface = "bold", color="#a01202") +
  annotate( geom="text", x= 0.25, y=1.13, label="Pre-fire", size=5 , fontface = "bold", color="goldenrod") + 
  annotate( geom="text", x= 0.37, y=0.6, label="Enhanced", size=5 , fontface = "bold", color="#1e5b24", angle=-90) +
  theme_bw() +
  theme(text = element_text(size = 20))

recovery.density.plot.new <- ggplot(data =time.diff.90.df) + 
  geom_density(aes(x= NDVI, y = ..scaled..),fill="#a01202", alpha=0.2 ,linetype = "solid" , color="#a01202") + 
  geom_density(aes(x= model.NDVI, y = ..scaled..), fill="#1e5b24", alpha=0.4,linetype = "dashed" , color="#1e5b24") + 
  geom_density( data= Recov_Combo, aes( x= PreNDVI, y = ..scaled..), fill="goldenrod", alpha=0.2 ,linetype = "solid" , color="goldenrod") + 
  annotate( geom="text", x= 0.1, y=1.13, label="Post-fire", size=5 , fontface = "bold", color="#a01202") +
  annotate( geom="text", x= 0.25, y=1.13, label="Pre-fire", size=5 , fontface = "bold", color="goldenrod") + 
  annotate( geom="text", x= 0.37, y=0.6, label="Enhanced", size=5 , fontface = "bold", color="#1e5b24", angle=-90) +
  theme_bw() +
  theme(text = element_text(size = 20))

# Combine 
pannel_bl <- ggarrange(train.NDVI.plot,
                       baseline.season.plot,
                       recovery.density.plot ,
                       #labels = c("B", "C", "D"),
                       ncol = 1, nrow = 3)


# MAP: Pre-Fire NDVI
# Spatial join to assign points to hexagons
hex_aggregated <- st_join(hex_grid_sf, drivers.sp, join = st_intersects) %>%
  group_by(geometry) %>%
  summarise(
    avg_NDVI = mean(PreNDVI, na.rm= TRUE),
    .groups = "drop")
# map
map_PreNDVI<-hex_aggregated %>%
  arrange(avg_NDVI) %>%  # Reorder data so lower values are plotted first
  ggplot() + 
  geom_sf(data = FL_clip, fill = "lightgrey", color = "lightgrey") + 
  geom_sf(data = EVG_bound, fill = NA, color = "black", size = 1.5) +    
  geom_sf(aes(fill = avg_NDVI), color=NA) +
  theme_classic() +
  scale_fill_gradientn(
    colors = c("#FFFFE0", "gold", "goldenrod", "#D97706"),    
    #colors = met.brewer("VanGogh3", direction = 1), 
    na.value = "transparent",
    breaks = c(0.1534, 0.25, 0.35),
    labels = c("0.1", "0.2", "0.3")) +
  scale_x_continuous(
    breaks = seq(-81.5, -80, by = 0.5)) + 
  scale_y_continuous(
    breaks = seq(24, 27, by = 0.5))+ 
  theme(text = element_text(size = 20),
        legend.position = "bottom") +
  labs(fill = "Pre-Fire NDVI") +
  # FL inlay
  annotation_custom( 
    grob = ggplotGrob(inlay_plot), 
    xmin = -81.6, xmax = -81.2, ymin = 24.8, ymax = 25.2) +
  # north arrow 
  annotation_north_arrow(
    location = "tr", 
    width = unit(1.5, "cm"), 
    height = unit(2, "cm"),
    style = north_arrow_orienteering()) +
  # scale bar
  annotation_scale(
    location = "bl",
    width_hint = 0.2, 
    height = unit(0.4, "cm"), 
    style = "ticks", 
    bar_cols = c("black", "white"),
    text_cex = 2)

# Final Plot
map_bl <- ggarrange(map_PreNDVI,
                    labels="A",
                    font.label = list( size = 20))

baseline_plot <- ggarrange(map_bl,
                           pannel_bl,
                           ncol = 2, nrow = 1)


# Save
ggsave("baseline_plot.png", plot = baseline_plot, width = 12, height = 8, dpi = 300)
ggsave("Figure3_NDVI-pre-post-enhanced.png", plot = recovery.density.plot.new, width = 6, height = 4, dpi = 300)


### FIGURE 5: Recovery #################################################################

# Recovery Time by Threshold
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/Recov_Combo2.RDATA")
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/Recov_Rates.RDATA")
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/Recov_time_thresholds.RDATA")

# Data Prep
# add threshold column
Recov_Time_100$thrshold <- "100%"
Recov_Time_90$thrshold <- "90%"
Recov_Time_80$thrshold <- "80%"
Recov_Time_70$thrshold <- "70%"
Recov_Time_60$thrshold <- "60%"
Recov_Time_50$thrshold <- "50%"
Recov_Time_undr50$thrshold <- "<50%"
# merge 
thresholds <- full_join(Recov_Time_100, Recov_Time_90)
thresholds <- full_join(thresholds, Recov_Time_80)
thresholds <- full_join(thresholds, Recov_Time_70)
thresholds <- full_join(thresholds, Recov_Time_60)
thresholds <- full_join(thresholds, Recov_Time_50)
thresholds <- full_join(thresholds, Recov_Time_undr50)
rm(Recov_Time_100, Recov_Time_90, Recov_Time_80, Recov_Time_70, Recov_Time_60, 
   Recov_Time_50, Recov_Time_undr50)
# fix factor order
thresholds$thrshold <- as.factor(thresholds$thrshold)
thresholds$thrshold <- factor(thresholds$thrshold, levels=c("<50%", "50%", "60%", "70%", "80%", "90%", "100%") )
levels(thresholds$thrshold)


# RECOVERY RATE
# load recovery table (made by Sparkle)
rec_tbl <- read.csv("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/recoverytable.csv")
# calculate range of rec times per threshold
rec_tbl$rec_time_range <- rec_tbl$max_time - rec_tbl$min_time
# plot recovery to x threshold over time
plot_recovery_rate <- ggplot(rec_tbl) +
  geom_smooth(aes(x=mean_time, y=R_percent),color="#555555", linetype="dotted")+
  geom_point(aes(x=mean_time, y=R_percent, size=locations_percent, color=R_percent)) +
  scale_size(range = c(3, 10))+
  labs(y="Percent Recovered",
       x="Mean Recovery Time (yrs)", 
       size="Percent of Locations",
       #color="Recovery Time Range (yrs)",
       color="Threshold", 
       tag="B") +
  scale_color_met_c("VanGogh3") +
  theme_bw()+
  theme(
    text = element_text(size = 20),
    legend.position= "none") + 
  #legend.position = c(0.95, 0.05),  # Place legend in the lower right corner
  #legend.justification = c(1, 0)) +   # Adjust to anchor the legend to the corner
  guides(color = "none")


# RECOVERY TIME DISTIRIBUTION BY THRESHOLD
plot_recovery_densityRecTime <- ggplot(
  data = thresholds, 
  aes(x = Rec_Yrs, y = thrshold, fill = thrshold, height = after_stat(density))) +
  geom_density_ridges(
    scale = 7,
    stat = "density",
    adjust = 1.2,
    trim = TRUE, # Trim ridges to actual data range
    color = "#555555") + # Keep a consistent outline for all ridges
  scale_fill_brewer(palette = 5) +
  theme_bw() +  
  theme(
    legend.position = "none",
    text = element_text(size = 20),
    panel.grid.major.y = element_blank(),  
    panel.grid.minor.y = element_blank()) +  
  labs( x = "Recovery Time (yrs)", y = "Threshold", fill="Threshold", tag="C")

# RECOVERY NDVI DISTRIBUTION BY THREHOLD
plot_recovery_densityNDVI <- ggplot(data = thresholds) +
  geom_density( color="darkgreen", alpha=.7, aes(x = NDVI, fill= thrshold)) +
  geom_density(aes(x=PreNDVI), color="goldenrod", size=1.2) +
  scale_fill_brewer(palette = 5) +
  theme_bw() + 
  theme(legend.position = "none",
        text = element_text(size = 20),
        panel.grid.major.y = element_blank(),  
        panel.grid.minor.y = element_blank()) + 
  labs(x = "NDVI", y = "Density", fill="Threshold", tag="C") +
  geom_vline(xintercept = mean(thresholds$PreNDVI), color="goldenrod", linetype="dashed", size=1.2) +
  annotate(geom="text", x= 0.15, y=30, label="Pre-Fire NDVI", size=6 , color="goldenrod") 


# Combined Recovery Plots
plot_recovery.1 <- ggarrange(plot_recovery_rate,
                             plot_recovery_densityRecTime,
                             nrow=1, ncol=2,
                             #labels=c("B", "C"),
                             font.label = list(size = 20))
plot_recovery.2 <- ggarrange(plot_recovery.1,
                             plot_recovery_densityNDVI, 
                             nrow=2, ncol=1,
                             font.label = list(size = 20))

plot_recovery.3 <- ggarrange(plot_recovery_rate,
                             plot_recovery_densityNDVI, 
                             nrow=2, ncol=1,
                             font.label = list(size = 20))





# MAX RECOVERY THRESHOLD MAP
# make threhsold numeric
drivers.sp$thresh.num <- NA
drivers.sp$thresh.num[drivers.sp$threshold == "<50"] <- 40
drivers.sp$thresh.num[drivers.sp$threshold == "50"] <- 50
drivers.sp$thresh.num[drivers.sp$threshold == "60"] <- 60
drivers.sp$thresh.num[drivers.sp$threshold == "70"] <- 70
drivers.sp$thresh.num[drivers.sp$threshold == "80"] <- 80
drivers.sp$thresh.num[drivers.sp$threshold == "90"] <- 90
drivers.sp$thresh.num[drivers.sp$threshold == "100"] <- 100

# Map (threshold)
# apply hex grid
hex_aggregated <- st_join(hex_grid_sf, drivers.sp, join = st_intersects) %>%
  group_by(geometry) %>%
  summarise(
    thresh.num = mean(thresh.num, na.rm = TRUE),
    .groups = "drop")
# map
map_recovery_MapMaxThresh <- hex_aggregated %>%
  arrange(thresh.num) %>%  
  ggplot() + 
  geom_sf(data = FL_clip, fill = "lightgrey", color = "lightgrey") + 
  geom_sf(data = EVG_bound, fill = NA, color = "black", size = 1.5) +    
  geom_sf(aes(fill = thresh.num), color=NA) +
  theme_classic() +
  theme(
    text = element_text(size = 30))+
  scale_fill_gradientn(
    colors = met.brewer("VanGogh3", n = 7, direction = 1),     
    na.value = "transparent",                             
    breaks = c(40, 50, 60, 70, 80, 90, 100),
    labels = c("<50", "50", "60","70", "80", "90","100")) + 
  scale_x_continuous(
    breaks = seq(-81.5, -80, by = 0.5)) + 
  scale_y_continuous(
    breaks = seq(24, 27, by = 0.5))+ 
  theme(text = element_text(size = 25),
        legend.position = "bottom") +
  labs(color = "Max Threshold",
       tag = "A") +
  # FL inlay
  annotation_custom( 
    grob = ggplotGrob(inlay_plot), 
    xmin = -81.6, xmax = -81.2, ymin = 24.8, ymax = 25.2) +
  # north arrow 
  annotation_north_arrow(
    location = "tr", 
    width = unit(1.5, "cm"), 
    height = unit(2, "cm"),
    style = north_arrow_orienteering()) +
  # scale bar
  annotation_scale(
    location = "bl",
    width_hint = 0.2, 
    height = unit(0.4, "cm"), 
    style = "ticks", 
    bar_cols = c("black", "white"),
    text_cex = 2)+
  guides(fill = "none")


# Final Figure
plot_recovery_final <- ggarrange(map_recovery_MapMaxThresh,
                                 plot_recovery.2,
                                 nrow = 1, ncol=2)

# save
# getwd()
ggsave("Recovery_thresholds_part1.png", plot = map_recovery_MapMaxThresh, width = 5, height = 6, dpi = 300)
ggsave("Recovery_thresholds_part2.png", plot = plot_recovery.3, width = 5, height = 6, dpi = 300)

### FIGURE 6: Drivers of Enhancement (>T80) #################################################################
# Recivery Threshold Index (T80_rf_index) 
# uses an 80% threshold recovery index (T80_rf_index)
# to compare locations that exceeded 80% recovered to those that did not.

load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Manuscript/RF_threshold_index.RDATA")
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Manuscript/RF_Rec_Yrs_index.RDATA")
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Manuscript/Sensitivity_data.RDATA")

varImpPlot(T80_rf_index)

# rec.status by driver variable
plot.T80_rf_index.PreNDVI <-  sensitivity.df %>% filter(target.var == 'PreNDVI') %>%  
  ggplot(aes(x = PreNDVI)) + 
  geom_density(aes(colour =T80_rf_index, fill=T80_rf_index), alpha=0.5)  + 
  theme_bw() +
  labs(x= "Pre-Fire NDVI", y="Density", color="Recovery Status", fill="Recovery Status") +
  scale_color_manual(
    values = c("#9bc184","#1e5b24")) +
  scale_fill_manual(
    values = c("#9bc184","#1e5b24"))+
  guides(fill = "none") + 
  theme(
    text = element_text(size = 20), 
    legend.position = "none") + 
  guides(color = "none")

plot.T80_rf_index.Severity <-sensitivity.df %>% filter(target.var == 'Severity') %>%  
  ggplot(aes(x = Severity)) + 
  geom_density(aes(colour =T80_rf_index, fill =T80_rf_index), alpha=0.5)  + 
  theme_bw() + 
  labs(x= "Severity", y="Density", color="Recovery Status", fill="Recovery Status") +
  scale_color_manual(
    values = c("#9bc184","#1e5b24")) +
  scale_fill_manual(
    values = c("#9bc184","#1e5b24"))+
  guides(fill = "none") + 
  theme(
    text = element_text(size = 20), 
    legend.position = "none") + 
  guides(color = "none")

plot.T80_rf_index.pdsi.sd <- sensitivity.df %>% filter(target.var == 'pdsi.sd') %>%  
  ggplot(aes(x = pdsi.sd)) + 
  geom_density(aes(colour =T80_rf_index, fill=T80_rf_index), alpha=0.5)  + 
  theme_bw() + 
  labs(x= "PDSI sd", y="Density", color="Recovery Status" , fill="Recovery Status") +
  scale_color_manual(
    values = c("#9bc184","#1e5b24")) +
  scale_fill_manual(
    values = c("#9bc184","#1e5b24"))+
  guides(fill = "none") + 
  theme(
    text = element_text(size = 20), 
    legend.position = "none") + 
  guides(color = "none")

plot.T80_rf_index.pdsi.max <-sensitivity.df %>% filter(target.var == 'pdsi.max') %>%  
  ggplot(aes(x = pdsi.max)) + 
  geom_density(aes(colour =T80_rf_index, fill=T80_rf_index), alpha=0.5)  + 
  theme_bw() + 
  labs(x= "PDSI maximum", y="Density", color="Recovery Status" , fill="Recovery Status") +
  scale_color_manual(
    values = c("#9bc184","#1e5b24")) +
  scale_fill_manual(
    values = c("#9bc184","#1e5b24"))+
  guides(fill = "none") + 
  theme(
    text = element_text(size = 20), 
    legend.position = "none") + 
  guides(color = "none")

plot.T80_rf_index.pdsi.min <-sensitivity.df %>% filter(target.var == 'pdsi.min') %>%  
  ggplot( aes(x = pdsi.min, y = level, fill = T80_rf_index)) + 
  geom_density_ridges_gradient() +  
  theme_bw()+
  scale_color_manual(
    values = c("#9bc184","#1e5b24")) +
  scale_fill_manual(
    values = c("#9bc184","#1e5b24")) +
  theme(
    text = element_text(size = 20),
    legend.position = "none",
    axis.title.x = element_text(hjust = 0.5),
    axis.title.y = element_text(hjust = 0.5)) +
  labs(x="PDSI minimum", y="Index", fill="Recovery Status") + 
  guides(color = "none")

plot.T80_rf_index.pdsi.mean <- sensitivity.df %>% filter(target.var == 'pdsi.mean') %>%  
  ggplot( aes(x = pdsi.mean, y = level, fill = T80_rf_index)) + 
  geom_density_ridges_gradient( ) +  
  theme_bw()+
  scale_color_manual(
    values = c("#9bc184","#1e5b24")) +
  scale_fill_manual(
    values = c("#9bc184","#1e5b24")) +
  theme(
    text = element_text(size = 20),
    legend.position = "none",
    axis.title.x = element_text(hjust = 0.5),
    axis.title.y = element_text(hjust = 0.5)) +
  labs(x="PDSI mean", y="Index", fill="Recovery Status") + 
  guides(color = "none")

plot.T80_rf_index.PostDateDif <- sensitivity.df %>% filter(target.var == 'PostDateDif') %>%  
  ggplot( aes(x = PostDateDif, y = level, fill = T80_rf_index)) + 
  geom_density_ridges_gradient( ) +  
  theme_bw()+
  scale_color_manual(
    values = c("#9bc184","#1e5b24")) +
  scale_fill_manual(
    values = c("#9bc184","#1e5b24")) +
  theme(
    text = element_text(size = 20),
    legend.position = "none",
    axis.title.x = element_text(hjust = 0.5),
    axis.title.y = element_text(hjust = 0.5)) +
  labs(x="Post-Fire Date Difference", y="Index", fill="Recovery Status") + 
  guides(color = "none")

plot.T80_rf_index.Prev.Int <- sensitivity.df %>% filter(target.var == 'Prev.Int') %>%  
  ggplot( aes(x = Prev.Int, y = level, fill = T80_rf_index)) + 
  geom_density_ridges_gradient( ) +  
  theme_bw()+
  theme_bw()+
  scale_color_manual(
    values = c("#9bc184","#1e5b24")) +
  scale_fill_manual(
    values = c("#9bc184","#1e5b24")) +
  theme(
    text = element_text(size = 20),
    legend.position = "none",
    axis.title.x = element_text(hjust = 0.5),
    axis.title.y = element_text(hjust = 0.5)) +
  labs(x="Time Since Fire", y="Index", fill="Recovery Status") + 
  guides(color = "none")

# combine
plot_T80_drivers <- ggarrange( 
  plot.T80_rf_index.PreNDVI,
  #plot.T80_rf_index.Prev.Int,
  #plot.T80_rf_index.Severity,
  #plot.T80_rf_index.PostDateDif,  
  plot.T80_rf_index.pdsi.max,
  plot.T80_rf_index.pdsi.min,  
  #plot.T80_rf_index.pdsi.sd, 
  #plot.T80_rf_index.pdsi.mean,
  nrow=3, ncol=1, labels= c( "A","B", "C"),common.legend = TRUE,
  font.label = list(size = 16))


# save
setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Manuscript/Figures")
#ggsave("Drivers_T80_pannel.png", plot = plot_T80_drivers, width = 10, height = 14, dpi = 300)
ggsave("Drivers_T80_pannel.png", plot = plot_T80_drivers, width = 5, height = 10, dpi = 300)

### FIGURE 7: Drivers of Recovery Time #############################################################################################################

load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Manuscript/DriversData.RDATA")
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Manuscript/RF_threshold_index.RDATA")
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Manuscript/RF_Rec_Yrs_index.RDATA")
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Manuscript/Sensitivity_data.RDATA")

#MODELS:
Rc_Yrs_rf_index.vars 
T80_rf_index.vars

# Recovery Time Index (Rc_Yrs_rf_index) 
# shows relationships between recovery time (yrs) and PDSI variables
# also plots recovery time spatially

# rec_time vs. PDSI
plot.Rc_Yrs_rf_inde.pdsi.max <-  sensitivity.df %>% filter(target.var == 'pdsi.max') %>%  
  ggplot(aes(x = pdsi.max, y =Rc_Yrs_rf_index, color = pdsi.max)) + 
  geom_smooth(color = "#6f948c", fill="lightgrey", size=2) + 
  geom_point(size=4) +  
  theme_bw() +
  scale_color_gradientn(
    colors = met.brewer("Ingres", n = 4, direction = -1),  
    breaks = c(-4, 0, 4),
    labels = c("-4", "0", "4"),
    limits = c(-4, 4),                    
    oob = scales::squish) +
  labs(y="Recovery Time", x= "PDSI maximum", color="PDSI") +
  theme(
    text = element_text(size = 20), 
    legend.position = "none")

plot.Rc_Yrs_rf_inde.pdsi.mean <- sensitivity.df %>% filter(target.var == 'pdsi.mean') %>%  
  ggplot(aes(x = pdsi.mean, y =Rc_Yrs_rf_index, color = pdsi.mean)) + 
  geom_smooth(color = "#6f948c", fill="lightgrey", size=2) + 
  geom_point(size=4) +  
  theme_bw() +
  scale_color_gradientn(
    colors = met.brewer("Ingres", n = 4, direction = -1),  
    breaks = c(-4, 0, 4),
    labels = c("-4", "0", "4"),
    limits = c(-4, 4),                    
    oob = scales::squish) +
  labs(y="Recovery Time", x= "PDSI mean", color="PDSI") +
  theme(
    text = element_text(size = 20), 
    legend.position = "none")

plot.Rc_Yrs_rf_inde.pdsi.min <- sensitivity.df %>% filter(target.var == 'pdsi.min') %>%  
  ggplot(aes(x = pdsi.min, y =Rc_Yrs_rf_index, color = pdsi.min)) + 
  geom_smooth(color = "#6f948c", fill="lightgrey", size=2) +   
  geom_point(size=4) +  
  theme_bw() + 
  scale_color_gradientn(
    colors = met.brewer("Ingres", n = 4, direction = -1),  
    breaks = c(-4, 0, 4),
    labels = c("-4", "0", "4"),
    limits = c(-4, 4),                    
    oob = scales::squish) +
  labs(y="Recovery Time", x= "PDSI minimum", color="PDSI") +
  theme(
    text = element_text(size = 20), 
    legend.position = "none")

plot.Rc_Yrs_rf_inde.pdsi.sd <- sensitivity.df %>% filter(target.var == 'pdsi.sd') %>%  
  ggplot(aes(x = pdsi.sd, y = Rc_Yrs_rf_index, color = pdsi.sd)) +  
  geom_smooth(color = "#6f948c", fill="lightgrey", size=2) +  
  geom_point(size=4) +  
  theme_bw() + 
  scale_color_gradientn(
    colors = met.brewer("Ingres", n = 4, direction = -1),  
    breaks = c(-4, 0, 4),
    labels = c("-4", "0", "4"),
    limits = c(-4, 4),                    
    oob = scales::squish) +
  labs(y = "Recovery Time", x = "PDSI sd", color="PDSI") +
  theme(
    text = element_text(size = 20), 
    legend.position = "none")

plot_pdsi <- ggarrange( 
  plot.Rc_Yrs_rf_inde.pdsi.max ,
  plot.Rc_Yrs_rf_inde.pdsi.mean,
  plot.Rc_Yrs_rf_inde.pdsi.min, 
  plot.Rc_Yrs_rf_inde.pdsi.sd, 
  nrow=2, ncol=2, labels = c("A","B", "C", "D"),
  common.legend = TRUE,
  font.label = list(size = 20))


setwd("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Manuscript/Figures")
ggsave("Drivers_RecTime.png", plot = plot_pdsi, width = 12, height = 10, dpi = 300)



### SUPPLIMENTAL #################################################################

# Correlation Plot: Baseline NDVI Model.........................................................................................
# load BL_master_df
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/BL_Master_filtered.RDATA")
# varaibles used in model 
BL_Master_filtered$TotalFires <- as.numeric(BL_Master_filtered$TotalFires)
M <-cor(BL_Master_filtered[, c( "Pt.B1max" , "Pt.B2max" ,"Pt.B3max" ,"Pt.B4max" ,"Pt.B5max", "Pt.B7max",
                                "Pt.B1min" ,"Pt.B2min" ,"Pt.B3min" ,"Pt.B4min" ,"Pt.B5min", "Pt.B7min",
                                "NIR.SWIR1", "SWIR1.SWIR2",
                                "Obs_Year", "Obs_month", "tmax", "precip", "tmin",
                                "TotalFires","Prev.Int",
                                "coords.x1", "coords.x2")]) 
# plot 
corplot_bl <- ggcorrplot(M, 
                         method = "circle",        
                         type = "lower",          
                         lab = FALSE,              
                         lab_size = 3,            
                         colors = c("red", "white", "blue")) +
  theme(text = element_text(size = 20)) 

ggsave("SUP.BLcorplot.png", plot = corplot_bl, width = 8, height = 8, dpi = 300)


# POST-FIRE NDVI TRENDS...........................................................................................
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/BL_Master_filtered.RDATA")
# assign frequency categories
BL_Master_filtered$TSF <- as.factor(BL_Master_filtered$TSF)
BL_Master_filtered$TotalFires <- as.factor(BL_Master_filtered$TotalFires)
# plot distribution
ggplot(BL_Master_filtered, aes(x=TotalFires)) +
  geom_histogram() +
  scale_x_continuous(breaks=1:9) +
  theme(text = element_text(size = 20))
# plot trend
ggplot(BL_Master_filtered, aes(x=TSF, 
                               y=NDVI,
                               group=TotalFires,
                               color=TotalFires)) +
  geom_smooth(formula = y ~ s(x, bs = "cs", k=7)) +
  labs(
    x="Years Post-Fire",
    y= "Mean NDVI",
    title = "Baseline Locations",
    color="Total Fires", 
    fill= "Fire Frequency") + 
  #scale_fill_manual(values=cols) +
  #scale_color_manual(values=cols)+
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE) +
  theme_bw()+
  theme(
    text = element_text(size = 20),
    #panel.background = element_rect(fill='transparent'), #transparent panel bg
    #plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    #panel.grid.major = element_blank(), #remove major gridlines
    #panel.grid.minor = element_blank(), #remove minor gridlines
    #legend.background = element_rect(fill='transparent'), #transparent legend bg
    #legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )
# save plot with transparent background
ggsave('FCEprez.png', px, bg='transparent')

ggplot(BL_Master_filtered, aes(x=TotalFires, y=NDVI, color=TotalFires)) +
  geom_boxplot()


# Boxplot: Total Fires by Veg Class..........................................................
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Landscape_Summary.RDATA")

# data prep
Fire_by_veg$rockland_flatwoods <- "Pine Flatwoods"
Fire_by_veg$rockland_flatwoods[Fire_by_veg$L4_name == "Pine Rockland"] <- "Pine Rockland"
# plot
fire_by_veg <- ggplot(Fire_by_veg, aes(x=DomCom, y=freq_1978_2020, fill= rockland_flatwoods)) +
  geom_boxplot(alpha = 0.7, color="#333333") +
  scale_fill_manual(values = c("lightgrey","#333333"))+
  theme_bw() +
  theme(text=element_text(size=20), axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(
    x = "Dominant Community",
    y = "Total Fires",
    fill = "Ecosystem Type")

ggsave("SUP.fire_by_veg.png", plot = fire_by_veg, width = 10, height = 8, dpi = 300)



# Bar charts: Fire History..............................................................................

# load fire perimeters clipped to pinelands 
Pine_fires_shp <- st_read("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Fire_History/Ecosyst_Clipped_Fires/Pine_fires_Km2.shp")
head(Pine_fires_shp)
# number of fires
WF <- Pine_fires_shp[which(Pine_fires_shp$FireType =="WF"),]
RX <- Pine_fires_shp[which(Pine_fires_shp$FireType =="RX"),]
Pine_fires_shp$StartDate <- as.Date(Pine_fires_shp$StartDate, format= "%Y-%m-%d")
Pine_fires_shp$FireMo <- format(Pine_fires_shp$StartDate,"%m")
Pine_fires_shp$FireMo <- as.numeric(Pine_fires_shp$FireMo)
NoNA <- Pine_fires_shp[which(Pine_fires_shp$FireType == "RX" | Pine_fires_shp$FireType == "WF"),]

# annual burned area 
annual_area <- ggplot(NoNA, aes(x=Year, y=Area_km2, fill=FireType)) +
  geom_col(position="dodge")+
  scale_fill_manual(values = c("#efc86d","#6f9969"),
                    labels=c("Prescribed", "Wild"))+
  theme_bw() +
  theme(text=element_text(size=20), 
        legend.position="none") +
  scale_x_discrete(breaks = c(1980, 1990, 2000, 2010, 2020))+
  labs(
    x = "Year",
    y = "Burned Area (Km2)",
    fill=" Fire Type",
    tag="D")

# annual number of fires
annual_fires <- ggplot(NoNA, aes(x=Year, fill=FireType)) +
  geom_bar(position="dodge")+
  scale_fill_manual(values = c("#efc86d","#6f9969"),
                    labels=c("Prescribed", "Wild"))+
  theme_bw() +
  theme(text=element_text(size=20), 
        legend.position="none") +
  scale_x_discrete(breaks = c(1980, 1990, 2000, 2010, 2020))+
  labs(
    x = "Year",
    y = "Number of Fires",
    fill=" Fire Type",
    tag="B")

# mean monthly fires
fires_summary <- NoNA %>% # monthly means
  group_by(Year,FireMo, FireType) %>%
  summarise(num_fires = n())
average_fires <- fires_summary %>%
  group_by(FireMo, FireType) %>%
  summarise(avg_fires_per_month = mean(num_fires))
average_fires$FireMo <- as.factor(average_fires$FireMo)
average_fires <- average_fires[which(average_fires$FireMo != "NA"),]
mean_monthly_fires <- ggplot(average_fires, aes(x=FireMo, y= avg_fires_per_month, fill=FireType)) +
  geom_col(position="dodge") +
  scale_fill_manual(values = c("#efc86d","#6f9969"),
                    labels=c("Prescribed", "Wild"))+
  theme_bw() +
  theme(text=element_text(size=20), 
        legend.position="none") +
  labs(
    x = "Month",
    y = "Mean Number of Fires",
    fill=" Fire Type",
    tag="A")

# mean monthly area
fires_summary <- NoNA %>% # monthly means
  group_by(Year,FireMo, FireType) %>%
  summarise(mean_area = mean(Area_km2))
average_fires <- NoNA %>%
  group_by(FireMo, FireType) %>%
  summarise(avg_area_per_month = mean(Area_km2))
average_fires$FireMo <- as.factor(average_fires$FireMo)
average_fires <- average_fires[which(average_fires$FireMo != "NA"),]
mean_monthly_area <- ggplot(average_fires, aes(x=FireMo, y= avg_area_per_month, fill=FireType)) +
  geom_col(position="dodge") +
  scale_fill_manual(values = c("#efc86d","#6f9969"),
                    labels=c("Prescribed", "Wild"))+
  theme_bw() +
  theme(text=element_text(size=20), 
        legend.position="bottom") +
  labs(
    x = "Month",
    y = "Mean Burned Area (Km2)",
    fill=" Fire Type",
    tag="C")

# combined plots
plot_fireHist <- ggarrange( 
  mean_monthly_fires,
  annual_fires,
  mean_monthly_area,
  annual_area ,  
  nrow=2, ncol=2, 
  #labels= c( "A","B", "C", "D", "E","F", "G", "H"),
  common.legend = TRUE)

ggsave("SUP.FireHist.png", plot = plot_fireHist, width = 14, height = 10, dpi = 300)


# Map: Recovery Status (>T80).........................................................................................
# format data
drivers.sp$rec.status <- factor(drivers.sp$rec.status, levels = c(">80%", "<80%"))
drivers.sp$rec.status.num <- 0
drivers.sp$rec.status.num[drivers.sp$rec.status == ">80%"] <- 1
# hexagon grid
hex_aggregated <- st_join(hex_grid_sf, drivers.sp, join = st_intersects) %>%
  group_by(geometry) %>%
  summarise(
    avg_rec.stat = mean(rec.status.num, na.rm = TRUE),
    .groups = "drop")
# round hexagon means back to bianary 
hex_aggregated$avg_rec.stat[hex_aggregated$avg_rec.stat < 0.5] <- 0
hex_aggregated$avg_rec.stat[hex_aggregated$avg_rec.stat >= 0.5] <- 1
hex_aggregated$avg_rec.stat <- as.factor(hex_aggregated$avg_rec.stat)
hex_aggregated$avg_rec.stat[hex_aggregated$avg_rec.stat == "NaN"] <- NA
# map
map_rec.status <- hex_aggregated %>% 
  filter(!is.na(avg_rec.stat)) %>%
  ggplot( ) + 
  geom_sf(data = FL_clip, fill="lightgrey", color="lightgrey") + 
  geom_sf(data = EVG_bound, fill = NA, color = "black", size = 1.5) +    
  geom_sf(aes(fill = factor(avg_rec.stat)), color=NA, size = 0.3) +
  theme_classic() +
  theme(text = element_text(size = 20),
        legend.position = "right") +
  labs(fill="Recovery Status") +
  scale_fill_manual(values = c("0" = "#9bc184", "1" = "#1e5b24"),
                    labels = c("0" = "<80%", "1" = ">80%"),
                    na.value = "transparent", drop=T) +
  scale_x_continuous(
    breaks = seq(-81.5, -80, by = 0.5)) + 
  scale_y_continuous(
    breaks = seq(24, 27, by = 0.5))+ 
  # FL inlay
  annotation_custom( 
    grob = ggplotGrob(inlay_plot), 
    xmin = -81.6, xmax = -81.2, ymin = 24.8, ymax = 25.2) +
  # north arrow 
  annotation_north_arrow(
    location = "tr", 
    width = unit(1.5, "cm"), 
    height = unit(2, "cm"),
    style = north_arrow_orienteering()) +
  # scale bar
  annotation_scale(
    location = "bl",
    width_hint = 0.2, 
    height = unit(0.4, "cm"), 
    style = "ticks", 
    bar_cols = c("black", "white"),
    text_cex = 2)

ggsave("SUP.map.T80.png", plot = map_rec.status, width = 12, height = 10, dpi = 300)


# Map: Recovery Time..................................................................................
# Spatial join to assign points to hexagons
hex_aggregated <- st_join(hex_grid_sf, drivers.sp, join = st_intersects) %>%
  group_by(geometry) %>%
  summarise(
    avg_RecYrs = mean(Rec_Yrs, na.rm = TRUE),
    avg_PDSI = mean(pdsi.mean, na.rm= TRUE),
    .groups = "drop")
# map
map_RecTime<-hex_aggregated %>%
  arrange(avg_RecYrs) %>%  # Reorder data so lower values are plotted first
  ggplot() + 
  geom_sf(data = FL_clip, fill = "lightgrey", color = "lightgrey") + 
  geom_sf(data = EVG_bound, fill = NA, color = "black", size = 1.5) +    
  geom_sf(aes(fill = avg_RecYrs), color=NA) +
  theme_classic() +
  scale_fill_gradientn(
    colors = met.brewer("Ingres", n = 7, direction = 1), 
    na.value = "transparent",                             
    breaks = c(1,  5, 10, 15, 20),
    labels = c("1",  "5",  "10", "15", "20")) +
  scale_x_continuous(
    breaks = seq(-81.5, -80, by = 0.5)) + 
  scale_y_continuous(
    breaks = seq(24, 27, by = 0.5))+ 
  theme(text = element_text(size = 20),
        legend.position = "right") +
  labs(fill = "Recovery Time (years)") +
  # FL inlay
  annotation_custom( 
    grob = ggplotGrob(inlay_plot), 
    xmin = -81.6, xmax = -81.2, ymin = 24.8, ymax = 25.2) +
  # north arrow 
  annotation_north_arrow(
    location = "tr", 
    width = unit(1.5, "cm"), 
    height = unit(2, "cm"),
    style = north_arrow_orienteering()) +
  # scale bar
  annotation_scale(
    location = "bl",
    width_hint = 0.2, 
    height = unit(0.4, "cm"), 
    style = "ticks", 
    bar_cols = c("black", "white"),
    text_cex = 2)

ggsave("SUP.map.RecTime.png", plot = map_RecTime, width = 12, height = 10, dpi = 300)



# Map: PDSI......................................................................................
map.pdsi <- hex_aggregated %>%
  arrange(avg_PDSI) %>%  
  ggplot() + 
  geom_sf(data = FL_clip, fill = "lightgrey", color = "lightgrey") + 
  geom_sf(data = EVG_bound, fill = NA, color = "black", size = 1.5) +    
  geom_sf(aes(fill = avg_PDSI), color=NA) +
  theme_classic() +
  scale_fill_gradientn(
    colors = met.brewer("Ingres", n = 7, direction = 1), 
    na.value = "transparent") +
  scale_x_continuous(
    breaks = seq(-81.5, -80, by = 0.5)) + 
  scale_y_continuous(
    breaks = seq(24, 27, by = 0.5))+ 
  theme(text = element_text(size = 20),
        legend.position = "right") +
  labs(fill = "Mean Recovery Period PDSI") +
  # FL inlay
  annotation_custom( 
    grob = ggplotGrob(inlay_plot), 
    xmin = -81.6, xmax = -81.2, ymin = 24.8, ymax = 25.2) +
  # north arrow 
  annotation_north_arrow(
    location = "tr", 
    width = unit(1.5, "cm"), 
    height = unit(2, "cm"),
    style = north_arrow_orienteering()) +
  # scale bar
  annotation_scale(
    location = "bl",
    width_hint = 0.2, 
    height = unit(0.4, "cm"), 
    style = "ticks", 
    bar_cols = c("black", "white"),
    text_cex = 2)


ggsave("SUP.map.PDSI.png", plot = map.pdsi, width = 12, height = 10, dpi = 300)






