tidy_v0_EVER_BICY <- sf::st_read(file.path("/", "Volumes", "malonelab", "Research", "ENP", "ENP Fire", "FireHistory", "EVER_BICY_1978_2023_perim.shp"))

# make dataframe
RECOV.master.sp <-tidy_v0_EVER_BICY
RECOV.master.df <- as.data.frame(tidy_v0_EVER_BICY)

# save
save(RECOV.master.sp, RECOV.master.df, file = file.path("/", "Volumes", "malonelab", "Research", "ENP", "ENP Fire", "FireHistory", "RECOV.master_082025.RDATA") )
