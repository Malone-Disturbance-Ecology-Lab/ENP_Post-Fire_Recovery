## --------------------------------------------- ##
#                 Fire Severity
## --------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script calculates fire severity as the change in NBR pre- and post-fire for the Recovery Sample Points.
## 1. Identifies last available pre-fire image and first available post-fire image
## 2. Calculates the time between imaging dates and fire start and end dates

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

rm(list=ls())

# Load necessary libraries
# If you don't have the "librarian" package, uncomment the next line and run it to install the package
# install.packages("librarian")
librarian::shelf(tidyverse)

## --------------------------------------------- ##
#              Calculate Severity -----
## --------------------------------------------- ##

# Load Recovery Master Dataframe
load(file.path("/", "Volumes", "MaloneLab", "Research", "ENP", "ENP Fire", "Grace_McLeod", "Recovery", "Recov_Master.RDATA"))

# Make a Severity dataframe
Sev_df <- Recov_Master %>% 
  dplyr::select(ptID, StartDate, EndDate, FireYear) %>%
  dplyr::distinct()

# Immediate pre and post fire spec ---------------------

# Create empty vectors to store info for later
max_obs_date_list <- rep(NA_character_, length(Sev_df$ptID))
max_nbr_list <- rep(NA_real_, length(Sev_df$ptID))
max_ndvi_list <- rep(NA_real_, length(Sev_df$ptID))
min_obs_date_list <- rep(NA_character_, length(Sev_df$ptID))
min_nbr_list <- rep(NA_real_, length(Sev_df$ptID))

# Start off cycling through points
# Loop through points in Sev_df
for (i in seq_along(Sev_df$ptID)){
  # Identify point ID
  pt <- Sev_df$ptID[i]
  
  message("Now calculating for ptID: ", pt, ", ", i, " out of ", length(Sev_df$ptID))
  
  # Subset obs for this point from master
  # Commenting out tidyverse filter function here because base R is way faster
  # subID <- Recov_Master %>%
  #   dplyr::filter(ptID == pt)
  subID <- Recov_Master[which(Recov_Master$ptID == pt),]
  
  # Immediate pre-fire
  pre_StartDate <- Sev_df %>%
    dplyr::filter(ptID == pt) %>%
    dplyr::pull(StartDate) %>%
    as.Date()
  
  # Pull out pre-fire obs 
  try(subPRE <- subID %>%
        dplyr::filter(Obs_Date < pre_StartDate) %>%
        dplyr::filter(Obs_Date == max(Obs_Date)))
  
  # Take max
  try(max_obs_date <- subPRE %>%
        dplyr::pull(Obs_Date) %>%
        as.character())
  
  try(max_nbr <- subPRE %>%
        dplyr::pull(NBR))
  
  try(max_ndvi <- subPRE %>%
        dplyr::pull(NDVI))
  
  # Immediate post-fire
  post_EndDate <- Sev_df %>%
    dplyr::filter(ptID == pt) %>%
    dplyr::pull(EndDate) %>%
    as.Date()
  
  # Pull out post-fire obs 
  try(subPOST <- subID %>%
        dplyr::filter(Obs_Date > post_EndDate) %>%
        dplyr::filter(Obs_Date == min(Obs_Date)))
  
  # Take min 
  try(min_obs_date <- subPOST %>%
        dplyr::pull(Obs_Date) %>%
        as.character())
  
  try(min_nbr <- subPOST %>%
        dplyr::pull(NBR))
  
  # Add stats to list
  try(max_obs_date_list[i] <- max_obs_date)
  try(max_nbr_list[i] <- max_nbr)
  try(max_ndvi_list[i] <- max_ndvi)
  try(min_obs_date_list[i] <- min_obs_date)
  try(min_nbr_list[i] <- min_nbr)
}

# Create new columns and fill them with the stats
Sev_df_v2 <- Sev_df %>%
  dplyr::mutate(PreDate = max_obs_date_list, 
                PreNBR = max_nbr_list, 
                PreNDVI = max_ndvi_list, 
                PostDate = min_obs_date_list, 
                PostNBR = min_nbr_list) %>%
  # Convert date columns
  dplyr::mutate(PreDate = as.Date(PreDate),
                StartDate = as.Date(StartDate),
                PostDate = as.Date(PostDate),
                EndDate = as.Date(EndDate))

# Check it: Are there any cases where pre date > post date? (data entry errors from the Park Service)
test <- Sev_df_v2 %>%
  dplyr::filter(StartDate > EndDate)

length(unique(test$ptID))
# If so, remove from dataframe
# Sev_df_v2 <- anti_join(Sev_df_v2, test)

# Severity Calculation ---------------------------------

Sev_df_v3 <- Sev_df_v2 %>%
  dplyr::mutate(DateDif = PostDate - PreDate,
                Severity = PreNBR - PostNBR) %>%
  # Round to just two decimal places
  dplyr::mutate(Severity = as.numeric(format(round(Severity, 2), nsmall = 2)))

# Check it out 
summary(Sev_df_v3$Severity)
summary(Sev_df_v3$PreNBR) 
summary(Sev_df_v3$PostNBR) # NAs because EndDate close to end of sample period and no clear image after OR EndDate not recorded

# Plot
ggplot(Sev_df_v3, aes(y = Severity)) +
  geom_boxplot() +
  scale_color_manual(values = c("#ff6633")) + 
  labs(y = "Fire Severity") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 20))

Sev.1 <- Sev_df_v3 %>%
  dplyr::filter(Severity >= .1)

pre <- ggplot(Sev.1, aes(y = PreNBR)) +
  geom_boxplot() +
  scale_color_manual(values = c("#ff6633"))+ 
  labs(y = "Pre-NBR") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 20))

post <- ggplot(Sev.1, aes(y = PostNBR)) +
  geom_boxplot() +
  scale_color_manual(values = c("#ff6633"))+ 
  labs(y="Post-NBR") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 20))

pre
post

## --------------------------------------------- ##
#         Pre- Post- Date Difference -----
## --------------------------------------------- ##

# Calculate range between start-date and pre-date ----------
Sev_df_v4 <- Sev_df_v3 %>%
  dplyr::mutate(PreDateDif = as.numeric(StartDate - PreDate))

summary(Sev_df_v4$PreDateDif)
#   Min.    1st Qu.    Median    Mean   3rd Qu.    Max.    
#  1.00     9.00      25.00     37.33   52.00    1028.00    

# Almost all points have a pre-fire obs within 100 days of fire, and top 3rd quartile have obs within 38 days. 
length(unique(Sev_df_v4$ptID[Sev_df_v4$PreDateDif <= 38])) #  99,577
length(unique(Sev_df_v4$ptID[Sev_df_v4$PreDateDif > 38])) # 57,746
# Reasonable to say you need pre-fire value within one month or two satalite passes! 

# PLOT to visualize distribution
ggplot(aes(x = PreDateDif, y = PreNBR), data = Sev_df_v4) +
  geom_point(color = "#ff6633") +
  geom_smooth(color = "black") +
  labs(x = "Days Pre-Fire", y = "Pre-Fire NBR") +
  ggtitle("Pineland")  + 
  theme(text = element_text(size = 20)) +
  ylim(-0.4, 0.4) +
  scale_x_reverse() +
  theme(panel.background = element_blank()) +
  theme(panel.grid.major = element_line(color = "black",
                                        linewidth = 0.05,
                                        linetype = 1)) 

# Calculate range between end-date and post-date ------------
Sev_df_v5 <- Sev_df_v4 %>%
  dplyr::mutate(PostDateDif = as.numeric(PostDate - EndDate))

summary(Sev_df_v5$PostDateDif)
#    Min.   1st Qu.    Median    Mean   3rd Qu.    Max.    NA's 
#    1.00   13.00     31.00     46.12   66.00    602.00    4266 
# Lots of NAs due to fires with no recorded EndDate. 

no_end <- Sev_df_v5 %>%
  dplyr::filter(is.na(EndDate))
# Sev_df_v5 <- Sev_df_v5 %>% drop_na(PostDateDif)
# Almost all points have a post-fire obs within 100 days of fire. 
length(unique(Sev_df_v5$ptID[Sev_df_v5$PostDateDif <= 41])) # 112,942
length(unique(Sev_df_v5$ptID[Sev_df_v5$PostDateDif > 41])) # 36,886
# Reasonable to say you need post-fire value within one month or two satalite passes! 

# SAVE updates 
write.csv(Sev_df_v5, file = "Severity_updated.csv")
save(Sev_df_v5, file = "Sev_df_updated.RDATA")



