# Recovery Analysis Results

# This script generates figures for the Everglades Post-Fire Recovery Manuscript
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(randomForest)
library(ggridges)
library(sf)
library(terra)
library(leaflet)
library(ggspatial)
library(ggmap)
library(MetBrewer)
library(grid)
library(ggpubr)

rm(list=ls())

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



Recov_Combo_minT <- Recov_Combo.postfire %>% 
  left_join( minTimeDiff.postfire, by='ptID') %>% 
  filter(Time.diff == Time.diff.min) %>% filter( Time.diff < 30) %>% 
  mutate( DeltsNDVIpost.model = NDVI - model.NDVI,
          DeltsNDVIpost = abs(NDVI - PreNDVI),
          DeltsNDVIpost.en = abs(NDVI - model.NDVI),
          DIFF.en = abs(model.NDVI-PreNDVI))


Recov_Combo_minT$DeltsNDVIpost %>% mean
Recov_Combo_minT$DeltsNDVIpost.en %>% mean



PreNDVI_pt <- Recov_Combo.postfire %>% reframe(.by=ptID, PreNDVI= mean(PreNDVI), model.NDVI= mean(model.NDVI)  ) 

# Reduction in NDVI
(PreNDVI_pt$PreNDVI %>% mean  - Recov_Combo_minT$NDVI %>% mean)/PreNDVI_pt$PreNDVI %>% mean * 100 # reduction in NDVI by fire

(PreNDVI_pt$PreNDVI %>% mean  - Recov_Combo_minT$NDVI %>% mean)/PreNDVI_pt$PreNDVI %>% sd/sqrt( length( PreNDVI_pt$PreNDVI))

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
Recov_Combo_minT$recovery_percent.en %>% mean

Recov_Combo_minT$DeltsNDVIpost %>% range %>% round
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



write.csv( Recovery_table, "/Users/sm3466/YSE Dropbox/Sparkle Malone/Research/ENP_Post-Fire_Recovery/Figures/Recovery_Table.csv")

ggplot(Recovery_table) +
  geom_point(aes(x=R_percent, y=PreNDVI_Recovery)) +
  geom_smooth(aes(x=R_percent, y=PreNDVI_Recovery), linetype="dotted", se = FALSE, formula=(y ~ exp(x)))


Recovery_table$PreNDVI[Recovery_table$R_percent >= 80] %>% mean
Recovery_table$model.NDVI[Recovery_table$R_percent >= 80] %>% mean

(Recovery_table$PreNDVI[Recovery_table$R_percent >= 80] %>% mean - Recovery_table$model.NDVI[Recovery_table$R_percent >= 80] %>% mean)/ Recovery_table$PreNDVI[Recovery_table$R_percent >= 80] %>% mean


TSF <- Recov_Combo %>% reframe( .by=ptID, TSF = mean(Prev.Int.og))

LUB <- Max_Thresh %>% left_join(TSF, by='ptID') %>% 
  mutate(LUB = case_when( TSF >=15 ~ "Long Unburned", TSF <15 ~ "Maintained"),
         R.threshold= thrshold %>% factor( levels= c(  "<50%","50%" , "60%" , "70%" ,"80%" ,"90%" ,"100%"))) %>% filter( R_percent > 60)

LUB.plot <- LUB %>% ggplot(aes( x=  R.threshold, y=Rec_Yrs, col=LUB) ) + 
  geom_violin( trim=TRUE) +
  stat_summary(fun.y=mean, geom="point", size=2) +
  scale_color_manual(values=c("black", "#E69F00"), name="") + theme_bw() + 
  ylab("Recovery Time (Years)") +
  xlab('Recovery Threshold')


ggsave("/Users/sm3466/YSE Dropbox/Sparkle Malone/Research/ENP_Post-Fire_Recovery/Figures/Supplemenatary_S2.png",
       LUB.plot, width = 6, height = 3, dpi = 300)


# plot recovery to x threshold over time
plot_recovery_rate <- ggplot(Recovery_table) +
  geom_smooth(aes(x=mean_time, y=R_percent),color="#555555", linetype="dotted") +
  geom_point(aes(x=mean_time, y=R_percent, size=locations_percent, color=R_percent)) +
  scale_size(range = c(3, 10)) +
  labs(y="% Recovered",
       x="Recovery Time (years)", 
       size="Locations (%)",
       color="Recovery Threshold") +
  scale_color_met_c("VanGogh3") +
  theme_bw()+
  theme( text = element_text(size = 15),
    legend.position= "top") + 
  guides(color = "none")



# RECOVERY NDVI DISTRIBUTION BY THREHOLD

plot_recovery_densityNDVI.new <- ggplot(data = thresholds) +
  stat_density( color="darkgreen", alpha=.7, aes(x = NDVI, y=..scaled.., fill= thrshold)) +
  stat_density(aes(x=PreNDVI, y=..scaled..), color="goldenrod", fill='transparent', size=1.2) + scale_fill_brewer(palette = 5) +
  theme_bw() +  ylim(0, 3) +
  theme(legend.position = "none",
        text = element_text(size = 15),
        panel.grid.major.y = element_blank(),  
        panel.grid.minor.y = element_blank()) + 
  labs(x = "NDVI", y = "Density", fill=" Recovery Threshold") +
  geom_vline(xintercept = mean(thresholds$PreNDVI), color="goldenrod", linetype="dashed", size=1.2) +
  annotate(geom="text", x= 0.15, y=2.8, label="Pre-Fire NDVI", size=6 , color="goldenrod") 


# Combined Recovery Plots

plot_recovery.3 <- ggarrange(plot_recovery_rate,
                             plot_recovery_densityNDVI.new, 
                             nrow=2, ncol=1,
                             font.label = list(size = 15))


load( "/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Spatial_files.RDATA")


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
  theme( text = element_text(size = 15))+
  scale_fill_gradientn(
    colors = met.brewer("VanGogh3", n = 7, direction = 1),     
    na.value = "transparent",                             
    breaks = c(40, 50, 60, 70, 80, 90, 100),
    labels = c("<50", "50", "60","70", "80", "90","100")) + 
  scale_x_continuous(
    breaks = seq(-81.5, -80, by = 0.5)) + 
  scale_y_continuous(
    breaks = seq(24, 27, by = 0.5))+ 
  theme(text = element_text(size = 15),
        legend.position = "bottom") +
  labs(color = "Max Threshold") +
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
    text_cex = 1.5)+
  guides(fill = "none")


# Final Figure


# save
# getwd()
ggsave("/Users/sm3466/YSE Dropbox/Sparkle Malone/Research/ENP_Post-Fire_Recovery/Figures/Recovery_thresholds_Final_1.png", map_recovery_MapMaxThresh, width = 4, height = 6, dpi = 300)

ggsave("/Users/sm3466/YSE Dropbox/Sparkle Malone/Research/ENP_Post-Fire_Recovery/Figures/Recovery_thresholds_Final_2.png", plot_recovery.3, width = 6, height = 6, dpi = 300)

# Longunburned:

long_unburned <- Recov_Combo.adj %>% filter( TotalFires.og < 2 , Prev.Int.og > 15)
long_unburned$NDVI %>% mean(na.rm=T)
(long_unburned$ptID %>% unique %>% length / Recov_Combo.adj$ptID %>% unique %>% length)*100

long_unburned.higerpreNDVI <- long_unburned %>% filter(PreNDVI > model.NDVI)

long_unburned.higerpreNDVI$ptID %>% unique %>% length / Recov_Combo.adj$ptID %>% unique %>% length *100

### FIGURE 6: Drivers of Enhancement (>T80) ######
# Recivery Threshold Index (T80_rf_index) 
# uses an 80% threshold recovery index (T80_rf_index)
# to compare locations that exceeded 80% recovered to those that did not.


# Drivers 
library(randomForest)
library(caret)

load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/RF_threshold_index_082025.RDATA")
T80_rf_index

varImpPlot(T80_rf_index)

train$T80_rf_index <- predict(T80_rf_index , train)
confusionMatrix(train$T80_rf_index, train$rec.status )

load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Manuscript/Sensitivity_data.RDATA")


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
    panel.border = element_blank(),
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
    panel.border = element_blank(),
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
    panel.border = element_blank(),
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
    panel.border = element_blank(),
    legend.position = "none") + 
  guides(color = "none")


plot.T80_rf_index.pdsi.min <-sensitivity.df %>% filter(target.var == 'pdsi.min') %>%  
  ggplot( aes(x = pdsi.min, y = level, fill = T80_rf_index), alpha=0.5) + 
  geom_density_ridges_gradient() +  
  theme_bw()+
  scale_color_manual(
    values = c("#9bc184","#1e5b24")) +
  scale_fill_manual(
    values = c("#9bc184","#1e5b24")) +
  theme(
    text = element_text(size = 20),
    legend.position = "none",
    panel.border = element_blank(),
    axis.title.x = element_text(hjust = 0.5),
    axis.title.y = element_text(hjust = 0.5)) +
  labs(x="PDSI minimum", y="", fill="Recovery Status") + 
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
    panel.border = element_blank(),
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
    panel.border = element_blank(),
    legend.position = "none",
    axis.title.x = element_text(hjust = 0.5),
    axis.title.y = element_text(hjust = 0.5)) +
  labs(x="Time Since Fire", y="Index", fill="Recovery Status") + 
  guides(color = "none")

# combine
plot_T80_drivers <- ggarrange( 
  plot.T80_rf_index.PreNDVI,
  plot.T80_rf_index.pdsi.max,
  plot.T80_rf_index.pdsi.min,  
  nrow=3, ncol=1, labels= c( "A","B", "C"),common.legend = TRUE,
  font.label = list(size = 16))


ggsave("/Users/sm3466/YSE Dropbox/Sparkle Malone/Research/ENP_Post-Fire_Recovery/Figures/Drivers_T80_pannel.png", plot = plot_T80_drivers, width = 5, height = 10, dpi = 300)

### FIGURE 7: Drivers of Recovery Time ####

load(file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Recovery/RF_Rec_Yrs_index_082025.RDATA")

load("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Manuscript/Sensitivity_data.RDATA")

#MODELS:
Rc_Yrs_rf_index 

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
    panel.border = element_blank(),
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
    panel.border = element_blank(),
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
    panel.border = element_blank(),
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
    panel.border = element_blank(),
    legend.position = "none")

plot_pdsi <- ggarrange( 
  plot.Rc_Yrs_rf_inde.pdsi.max ,
  plot.Rc_Yrs_rf_inde.pdsi.mean,
  plot.Rc_Yrs_rf_inde.pdsi.min, 
  plot.Rc_Yrs_rf_inde.pdsi.sd, 
  nrow=2, ncol=2, labels = c("A","B", "C", "D"),
  common.legend = TRUE,
  font.label = list(size = 20))

ggsave("/Users/sm3466/YSE Dropbox/Sparkle Malone/Research/ENP_Post-Fire_Recovery/Figures/Drivers_RecTime.png", plot = plot_pdsi, width = 12, height = 10, dpi = 300)



