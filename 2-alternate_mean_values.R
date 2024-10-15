#find average per cell instead
lakes <- read_sf("D:/floating_solar/Northeast_NHD_Alison")
lake_buffer <- st_buffer(lakes,5000)
lakes_vec <- vect(lake_buffer)#create SpatVector
writeVector(lakes_vec, file = "D:/floating_solar/generated/lake_5k_buffer_vec.gpkg")

#project to ebird data crs
sp <- "amerob" #species doesn't matter, using for projection
bird_data <- rast(paste0("D:/floating_solar/generated/",sp,"_max_values.tif"))
lakes_vec_pro <- project(lakes_vec, crs(bird_data))
#plot(lakes_vec_pro)#look to check it worked
#plot(bird_data, add=T)

complete <- list.files(path = "D:/floating_solar/generated/")
complete_codes <- str_extract(complete,"[^_]+")

rm(sp)
rm(bird_data)

for(a in 1:length(complete_codes)){
  
  sp <- complete_codes[a]
  
  bird_data <- rast(paste0("D:/floating_solar/generated/",sp,"_max_values.tif"))
  bird_data1 <- bird_data*10000 #transformed to make numbers nicer to deal with
  
  bird_data1[is.na(bird_data1[])] <- 0 
  
  lake_bird_data_mean <- zonal(bird_data1, z = lakes_vec_pro, fun = "mean", na.rm=T) 

  lake_bird_data_mean$species_code <- rep(sp, length(lake_bird_data$max))
  lake_bird_data_mean$Water_ID <- lakes$Water_ID
  
  save(lake_bird_data_mean, file = paste0("D:/floating_solar/data_outputs/",sp,"_lake_abd_MEAN_weight.RData"))
  
}

species <- list.files(path = "D:/floating_solar/generated/")
species_codes <- str_extract(species,"[^_]+")

lake_biodiversity <- list()

for(s in 1:length(species_codes)){
  
  sp <- species_codes[s]
  
  #the sum of abd importance for each lake 
  load(paste0("D:/floating_solar/data_outputs/",sp,"_lake_abd_MEAN_weight.RData"))
  
  lake_biodiversity[[s]] <- lake_bird_data_mean
  
}
