rm(list=ls())


load(file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/Recov_Combo_082025.RDATA")
load( "/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Spatial_files.RDATA")
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/BL_train_test_08082025.RDATA")
load( file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/NDVI_rf_08082025.RDATA")
load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/Recov_Drivers_082025.RDATA")

library(tidyverse)


# RUN TO SEE WHAT ALL THIS IS ABOUT! 
# Subset the data needed:
Recov_Combo.adj <- subset(Recov_Combo, ptID %in% Recov_Drivers$ptID) %>% distinct()

Recovery.baseline <- Recov_Combo.adj %>% select( ptID, StartDate, model.NDVI, NDVI, PreNDVI, precip, tmax,
                                                  Prev.Int, TotalFires, Obs_Date) %>% 
  mutate(
  Fire.date = as.Date(StartDate, "%Y-%m-%d"),
  NDVI.date = as.Date(Obs_Date, "%Y-%m-%d"),
  Obs_month = format(NDVI.date, '%m') %>% as.numeric(),
  Time.diff = (NDVI.date - Fire.date) %>% as.numeric())

Recovery.baseline.minimum <- Recovery.baseline %>% group_by(ptID) %>% summarise( Time.diff.min = min(Time.diff)) %>% 
  full_join(Recovery.baseline, by = 'ptID')


Recovery.baseline.start <- Recovery.baseline.minimum %>%  filter( Time.diff.min == Time.diff) %>% mutate(NDVI.Diff = model.NDVI - NDVI)

# For the MS we need to show how the baseline compares to NDVI right before the fire: 
time.diff.90.df <- Recovery.baseline.start %>% filter( Time.diff %>% as.numeric <= 90) 

mean( Recovery.baseline.start$model.NDVI) - mean(Recovery.baseline.start$NDVI)

(mean( Recov_Combo$PreNDVI ) - mean(Recovery.baseline.start$model.NDVI))/ mean( Recov_Combo$PreNDVI )*100


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
ndvi.enhanced <- function( df, train){
  
  df$TotalFires <- 9
  df$Prev.Int <- mean(train$Prev.Int[train$TotalFires == 9])
  df$Pt.B4max <- mean(train$Pt.B4max[train$TotalFires == 9])
  df$Pt.B4min <- mean(train$Pt.B4min[train$TotalFires == 9])
  df$Pt.B5max <- mean(train$Pt.B5max[train$TotalFires == 9])
  df$SWIR1.SWIR2 <- mean(train$SWIR1.SWIR2[train$TotalFires == 9])
  df$NIR.SWIR1 <-  mean(train$NIR.SWIR1[train$TotalFires == 9])
  
  library(randomForest)
  
  df$NDVI.enhanced <- predict(NDVI_rf, newdata=df )
  
  return(df)
}

time.diff.90.df <- time.diff.90.df  %>% ndvi.enhanced( train)

train.predtions <- train %>% ndvi.enhanced ( train)
train.predtions$NDVI.model <- predict(NDVI_rf, data =train.predtions )

save(train.predtions, time.diff.90.df, Recovery.baseline.start,Recovery.baseline.start.1,
     Recovery.baseline.start.6, file = "/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Figure3_Data_082025.RDATA")


Recov_NDVI <-  Recovery.baseline  %>% filter(Obs_Date > StartDate) %>% select( ptID, model.NDVI, NDVI, PreNDVI ) %>% reframe(.by = ptID, model.NDVI = max(model.NDVI) ,NDVI = min(NDVI) ,PreNDVI = mean(PreNDVI) )

recovery.density.plot <- ggplot(data =Recov_NDVI) + 
  geom_density(aes(x= NDVI, y = ..scaled..),fill="#a01202", alpha=0.2 ,linetype = "solid" , color="#a01202") + 
  geom_density(data=Recov_Combo, aes(x= model.NDVI, y = ..scaled..), fill="#1e5b24", alpha=0.4,linetype = "dashed" , color="#1e5b24") + 
  geom_density(aes( x= PreNDVI, y = ..scaled..), fill="goldenrod", alpha=0.2 ,linetype = "solid" , color="goldenrod") + 
  labs(y="Density", x="NDVI") +
  annotate( geom="text", x= 0.1, y=1.13, label="Post-fire", size=5 , fontface = "bold", color="#a01202") +
  annotate( geom="text", x= 0.25, y=1.13, label="Pre-fire", size=5 , fontface = "bold", color="goldenrod") + 
  annotate( geom="text", x= 0.37, y=0.6, label="Enhanced", size=5 , fontface = "bold", color="#1e5b24", angle=-90) +
  theme_bw() +
  theme(text = element_text(size = 20))

ggsave("baseline_plot.png", plot = baseline_plot, width = 12, height = 8, dpi = 300)
ggsave(filename="/Users/sm3466/YSE Dropbox/Sparkle Malone/Research/ENP_Post-Fire_Recovery/Figures/Figure4_NDVI-pre-post-enhanced.png", plot = recovery.density.plot, width = 6, height = 4, dpi = 300)




# Additional plots

train.NDVI.plot <- ggplot(data =train.predtions) + 
  geom_density(aes(x= NDVI, y = ..scaled..), color="black" , size=1.5) + 
  geom_density(aes(x= NDVI.model, y = ..scaled..),linetype = "dashed", color="#1e5b24", size=1.5 ) + 
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


recovery.density.plot.new <- ggplot(data = Recov_Combo) + 
  geom_density(aes(x= NDVI, y = ..scaled..),fill="#a01202", alpha=0.2 ,linetype = "solid" , color="#a01202") + 
  geom_density(aes(x= model.NDVI, y = ..scaled..), fill="#1e5b24", alpha=0.4,linetype = "dashed" , color="#1e5b24") + 
  geom_density( data= Recov_Combo, aes( x= PreNDVI, y = ..scaled..), fill="goldenrod", alpha=0.2 ,linetype = "solid" , color="goldenrod") + 
  annotate( geom="text", x= 0.1, y=1.13, label="Post-fire", size=5 , fontface = "bold", color="#a01202") +
  annotate( geom="text", x= 0.25, y=1.13, label="Pre-fire", size=5 , fontface = "bold", color="goldenrod") + 
  annotate( geom="text", x= 0.37, y=0.6, label="Enhanced", size=5 , fontface = "bold", color="#1e5b24", angle=-90) +
  theme_bw() +
  theme(text = element_text(size = 20))

library(ggpubr)
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

