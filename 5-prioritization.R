library(prioritizr)
library(sf)
library(tidyverse)
library(ebirdst)
library(terra)
library(stringr)
library(viridis)
library(scales)



# 1. Create current zone layer --------------------------------------------
#create current zone layer (raster stack)

species_data <- list.files(path = "D:/floating_solar/generated/")

current_zone_layer <- rast(paste0("D:/floating_solar/generated/",species_data))




# 2. Create solar zone ----------------------------------------------------
#create solar zone layer using VI for each species


VI_data <- read.csv("data_outputs/final_analysis_data.csv")

#assumed that species with the highest VI decrease by 90% of their max values
#and species with the lowest VI decrease by 10% of their max values

VI_data$percent_reduction <- rescale(VI_data$VI, to=c(0.1,0.9))

#sort so that values align with raster stack
VI_data_arr <- VI_data %>%
  arrange(species_code)

#only regions where solar panels could possibly be implemented will be affected
#ie, reductions will only occur on lakes where there is solar
#start by getting lake buffer
lakes <- read_sf("D:/floating_solar/Northeast_NHD_Alison")
buf <- 5000
lake_buffer <- st_buffer(lakes,buf)

#filter out unsuitable lakes - these values won't change
lake_buffer_s <- lake_buffer %>%
  filter(Suitabl_FP ==1)

lakes_vec <- vect(lake_buffer_s)#create SpatVector
#project to ebird data crs
sp <- "amerob" #species doesn't matter, using for projection
bird_data <- rast(paste0("D:/floating_solar/generated/",sp,"_max_values.tif"))
lakes_vec_pro <- project(lakes_vec, crs(bird_data))

#mask the raster stack to the lake buffer
masked_raster <- mask(current_zone_layer,lakes_vec_pro)

#calculate new values based on the scaled VI
solar_reduced_values <- masked_raster*VI_data_arr$percent_reduction

#solar zone layer
#solar_zone_layer <- ifelse(is.na(masked_raster),current_zone_layer,solar_reduced_values)
#solar_zone_layer <- rast(solar_zone_layer)

solar_zone_layer <- terra::merge(solar_reduced_values,current_zone_layer, first=T)


# 3. Trim both zones to study region --------------------------------------

#map data
state_lines <- st_read("D:/maps/NA_politicalboundaries_shapefile/PoliticalBoundaries_Shapefile/NA_PoliticalDivisions/data/bound_p/boundaries_p_2021_v3.shp")

NE <- state_lines %>%
  filter(NAME_En %in% c("Maine","New Hampshire","Vermont","Massachusetts","Rhode Island","Connecticut","New York","Pennsylvania","New Jersey","Delaware","Maryland","West Virginia","Virginia","District of Columbia"))%>%
  st_geometry()

NE_vec <- vect(NE)#create SpatVector

NE_pro <- project(NE_vec, crs(bird_data))

current_zone_layer_trim <- mask(current_zone_layer,NE_pro)

solar_zone_layer_trim <- mask(solar_zone_layer,NE_pro)

#save raster stacks
writeRaster(current_zone_layer_trim, "D:/floating_solar/data_outputs/current_zone_stack.tif")
writeRaster(solar_zone_layer_trim, "D:/floating_solar/data_outputs/solar_zone_stack.tif")
