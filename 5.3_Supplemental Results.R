# Results : Supplemental 

# This script generates figures for the Everglades Post-Fire Recovery Manuscript
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(randomForest)
library(sf)
library(terra)
library(ggspatial)
library(ggmap)
library(MetBrewer)
library(grid)
library(ggpubr)

rm(list=ls())

# Baseline:
load(file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/BL_Master_df_08082025.RDATA")

load( file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/NDVI_rf_08082025.RDATA")

load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/BL_train_test_08082025.RDATA")

# Add model projections in the files and make 121 with NDVI model and enhanced NDVI model.
BL_Master_filtered.nona.pineland %>% names()

BL_Master_filtered.nona.pineland$model.NDVI <- predict(NDVI_rf, newdata=BL_Master_filtered.nona.pineland )

NDVI.ENHANCED <- function( df){
  df.enhanced <- df
  df.enhanced$Prev.Int <- train$Prev.Int[ train$TotalFires == 9 ] %>% mean
  df.enhanced$TotalFires <- 9
  df.enhanced$NIR.SWIR1 <- train$NIR.SWIR1[ train$TotalFires == 9 ] %>% mean
  df.enhanced$SWIR1.SWIR2 <- train$SWIR1.SWIR2[ train$TotalFires == 9 ] %>% mean
  df.enhanced$Pt.B4max <- train$Pt.B4max[ train$TotalFires == 9 ] %>% mean
  df.enhanced$Pt.B5max <- train$Pt.B5max[ train$TotalFires == 9 ] %>% mean
  df.enhanced$Pt.B4min <- train$Pt.B4min[ train$TotalFires == 9 ] %>% mean
 
  df.enhanced$NDVI.enhanced <- predict(NDVI_rf, newdata=df.enhanced )
  df$NDVI.enhanced <- df.enhanced$NDVI.enhanced 
  return(df)
}

predictionsNDVI <- NDVI.ENHANCED(df =  BL_Master_filtered.nona.pineland)

predictionsNDVI %>% names()

one2one.enhanced <- ggplot(data = predictionsNDVI)  + geom_point(aes( x= model.NDVI, y = NDVI.enhanced)) +
  theme_bw() +
  ylab('Enhanced NDVI') + xlab('NDVI Model') + theme(text = element_text(size = 8)) + 
  geom_abline(intercept = 0, slope = 1, col="red", linetype="dashed") + 
  ylim(0,0.5)+ xlim(0,0.5)


ggplot(data = predictionsNDVI)  + geom_density(aes( x= model.NDVI), linetype="dotted") + 
  geom_density(aes( x= NDVI.enhanced)) + theme_bw() + xlab('NDVI') + ylab('Density')





# Recovery ####

load(file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/Recov_Combo_082025.RDATA")
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/Recov_Drivers_082025.RDATA")


# Filter recovery by Driver Analysis DATA
Recov_Combo.adj <- subset(Recov_Combo, ptID %in% Recov_Drivers$ptID)

Recov_Combo.adj$ptID %>% unique %>% length


Recov_Combo.postfire <- Recov_Combo.adj %>% select( ptID, StartDate, Obs_Date, model.NDVI, NDVI, PreNDVI) %>% mutate(
  Fire.date = as.Date(StartDate, "%Y-%m-%d"),
  NDVI.date = as.Date(Obs_Date, "%Y-%m-%d"),
  Time.diff = (NDVI.date - Fire.date) %>% as.numeric) %>% filter(NDVI.date >= Fire.date)

minTimeDiff.postfire  <- Recov_Combo.postfire %>% reframe(.by=ptID, Time.diff.min= min(Time.diff) ) 

Recov_Combo_minT <- Recov_Combo.postfire %>% left_join( minTimeDiff.postfire, by='ptID') %>% filter(Time.diff == Time.diff.min) %>% filter( Time.diff < 30) %>% mutate( DeltsNDVIpost.model = NDVI - model.NDVI,
                                                                                                                                                                        DeltsNDVIpost = abs(NDVI - PreNDVI),
                                                                                                                                                                        DIFF.en = abs(model.NDVI-PreNDVI),
                                                                                                                                                                        recovery_percent.en = (DIFF.en/PreNDVI)*100)
PreNDVI_pt <- Recov_Combo.postfire %>% reframe(.by=ptID, PreNDVI= mean(PreNDVI), model.NDVI= mean(model.NDVI)  ) 


(PreNDVI_pt$PreNDVI %>% mean  - Recov_Combo_minT$NDVI %>% mean)/PreNDVI_pt$PreNDVI %>% mean * 100 # reduction in NDVI by fire

# Prefire NDVI
PreNDVI_pt$PreNDVI %>% range
PreNDVI_pt$PreNDVI %>% mean # Mean prefire NDVI
PreNDVI_pt$PreNDVI %>% sd / sqrt(length(PreNDVI_pt$PreNDVI ))  # SE of prefire NDVI

# Immediately postfire NDVI
Recov_Combo_minT$NDVI %>% mean # Average Postfire NDVI
Recov_Combo_minT$NDVI %>% range # NDVI right After the fire

# The Enhanced NDVI:
PreNDVI_pt$model.NDVI %>% range
PreNDVI_pt$model.NDVI %>% mean # Mean prefire NDVI
PreNDVI_pt$model.NDVI %>% sd / sqrt(length(PreNDVI_pt$PreNDVI ))  # SE of prefire NDVI


# Recovery Rate Comparison:
Recov_Combo_minT$DeltsNDVIpost %>% mean
Recov_Combo_minT$DeltsNDVIpost %>% sd / sqrt(length(Recov_Combo_minT$PreNDVI )) 
Recov_Combo_minT$DeltsNDVIpost.model%>% mean
Recov_Combo_minT$DeltsNDVIpost.model %>% sd / sqrt(length(Recov_Combo_minT$PreNDVI )) 

# Recovery mean 
Recov_Combo_minT$recovery_percent.en %>% mean

# Recovery Time by Threshold
load(file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/Max_Thresh_082025.RDATA")

Recov_Combo$ptID %>% unique %>% length
rec_tbl_new %>% mutate( )
rec_tbl_new$locations %>% sum

Max_Thresh <- Max_Thresh %>% mutate( PreNDVI_Recovery=  NDVI - PreNDVI)
Max_Thresh %>% names
Recovery_table <- Max_Thresh %>% mutate(
  DIFF.en = abs(model.NDVI-PreNDVI),
  recovery_percent.en = (DIFF.en/PreNDVI)*100) %>% reframe( .by = thrshold,
                                                            R_percent =  R_percent %>% as.numeric,
                                                            locations = length(ptID %>% unique),
                                                            mean_time =  mean(Rec_Yrs, na.rm=T), 
                                                            max_time =  max(Rec_Yrs, na.rm=T), 
                                                            min_time=  min(Rec_Yrs, na.rm=T),
                                                            model.NDVI = mean(model.NDVI),
                                                            PreNDVI = mean(PreNDVI),
                                                            NDVI = mean(NDVI)) %>% distinct %>%  mutate( 
                                                              locations_percent = locations/sum(locations)*100,
                                                              locations_percent = cumsum(locations_percent),
                                                              sample_size =  cumsum(locations),
                                                              rec_time_range = max_time - min_time,
                                                              PreNDVI_Recovery=  NDVI - PreNDVI)


load( "/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Spatial_files.RDATA")

load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/RF_threshold_index_082025.RDATA")

load(file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/BL_Master_df_08082025.RDATA")

# assign frequency categories
BL_Master_filtered.nona.pineland <- BL_Master_filtered.nona.pineland %>% mutate( 
  TSF =  as.factor(TSF),
  TotalFires =  as.factor(TotalFires))


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

ggsave("/Users/sm3466/YSE Dropbox/Sparkle Malone/Research/ENP_Post-Fire_Recovery/Figures/SUP.map.T80.png", plot = map_rec.status, width = 12, height = 10, dpi = 300)


# Map: Recovery Time..................................................................................
# Spatial join to assign points to hexagons
hex_aggregated <- st_join(hex_grid_sf, drivers.sp, join = st_intersects) %>%
  group_by(geometry) %>%
  summarise(
    avg_RecYrs = mean(Rec_Yrs, na.rm = TRUE),
    avg_PDSI = mean(pdsi.mean, na.rm= TRUE),
    .groups = "drop")
# map
map_RecTime <-hex_aggregated %>%
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

ggsave("/Users/sm3466/YSE Dropbox/Sparkle Malone/Research/ENP_Post-Fire_Recovery/Figures/SUP.map.RecTime.png", plot = map_RecTime, width = 12, height = 10, dpi = 300)