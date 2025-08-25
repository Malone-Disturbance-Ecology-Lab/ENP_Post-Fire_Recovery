# Spatial auto-correlation:
rm(list=ls())

library(randomForest)
library(sf)
library(dplyr)
library(ggplot2)

load(file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/BL_train_test_08082025.RDATA")

# load upland sample points
BL_smpl_pts <- sf::st_read(dsn = "/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Sampling/BL_smpl_pts.shp")  %>% st_transform(crs = 4326) # sf file

# Add lat and long into the file:
BL_smpl_pts$lat <- sf::st_coordinates(BL_smpl_pts)[,1]
BL_smpl_pts$lon <- sf::st_coordinates(BL_smpl_pts)[,2]
BL_df <- as.data.frame(BL_smpl_pts)

# Load the models:
load(file="/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/NDVI_rf_08082025.RDATA")

# Residual Krigiing:
library(gstat)
library(randomForest)

train$pred1 = predict(NDVI_rf, train )

train.res <- train %>% mutate( 
  res1= NDVI - pred1) %>% filter(!is.na(res1))

train.res %>% ggplot()+geom_point( aes(x=NDVI, y = pred1))
train.res %>% ggplot()+ geom_histogram( aes(x=res1))

sf_obs <- sf::st_as_sf(train.res, coords = c("lat","lon"), crs=4326)

EVG_bound <- st_read("/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/AOI/EVG_bound.shp") %>% st_transform(crs = 4326)

ggplot() + geom_sf( data=EVG_bound ) + geom_sf( data=sf_obs, aes(col=res1 ))

# Moran's I Analysis:
library(spdep)

## Generate distance Marix:
sf_obs$lat <- sf::st_coordinates(sf_obs)[,1]
sf_obs$lon <- sf::st_coordinates(sf_obs)[,2]

# Mean residuals by the pointID:

sf_obs_summary <- sf_obs %>% as.data.frame %>%  reframe( .by=c(ptID, lat, lon), res=mean(res1))

sf_obs_summary_sp <- sf_obs_summary %>% st_as_sf(coords= c("lat", "lon"), crs=st_crs(sf_obs))

neighbors_list <- dnearneigh(st_coordinates(sf_obs_summary_sp), 
                             d1 = 0, d2 = 100)
neighbor_distances <- nbdists(neighbors_list, st_coordinates(sf_obs_summary_sp ))

# Convert the neighbor list to a spatial weights object
listw_obj <- nb2listw(neighbors_list, style = "B") # "B" for binary weights

moran.test(sf_obs_summary_sp$res, listw_obj, alternative="greater")

moran.results <- moran.mc(sf_obs_summary_sp$res, listw_obj, alternative="greater", nsim=500)

# Vizualize spatial autocorrelation in the mean residual for points:

# Create prediction grid:
grid_pts <- expand.grid(x = seq(-81.713, -80.38847, by = 0.01), y = seq( 24.85675, 26.25934, by = 0.001))

sf_grid <- st_as_sf(grid_pts, coords = c("x", "y"), crs = 4326) %>% st_intersection(EVG_bound)

ggplot() + geom_sf( data=sf_grid  ) + geom_sf( data=EVG_bound )

# 3. Create and fit variogram model

vgm_emp <- variogram(res1 ~ 1, locations = sf_obs)
vgm_mod <- fit.variogram(vgm_emp, model =vgm(psill=100, model="Sph", range=7, nugget=1))

vgm # see model options

plot(vgm_emp, vgm_mod) # check the fit of the model used

png(filename="/Users/sm3466/YSE Dropbox/Sparkle Malone/Research/ENP_Post-Fire_Recovery/Figures/Baseline_Variogram.png",
    width = 500, height=500, res=400)
plot(vgm_emp, vgm_mod)
dev.off()

# 4. Perform Kriging
kriged_output <- krige(res1 ~ 1, locations = sf_obs, newdata = sf_grid, model = vgm_mod, maxdist=1000)

# You can then visualize the kriged_output (e.g., using ggplot2 with geom_sf)
              

# Save the files for use in the ...  

save(vgm_emp,sf_obs, vgm_mod, moran.results, file='/Volumes/MaloneLab/Research/ENP/ENP Fire/Grace_McLeod/Baseline/Spatial_AutoCorrelation.RDATA'  )
