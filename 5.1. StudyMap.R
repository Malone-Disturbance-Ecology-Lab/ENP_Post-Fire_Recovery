library(ggplot2)
library(tidyverse)
library(sf)
library(terra)
library(ggpubr)

rm(list=ls())

load( "/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Spatial_files.RDATA")

### FIGURE 1: Pinelands Map 
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
ggsave("/Users/sm3466/YSE Dropbox/Sparkle Malone/Research/ENP_Post-Fire_Recovery/Figures/Map_Pinelands_06212025.png", plot = map.pinelands_plot_new, width = 12, height = 10, dpi = 300)
