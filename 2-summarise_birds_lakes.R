library(sf)
library(tidyverse)
library(ebirdst)
library(terra)
library(stringr)

#load in lake data
lakes <- read_sf("D:/floating_solar/Northeast_NHD_Alison")
plot(st_geometry(lakes[1,]))
plot(st_geometry(lakes))

lake_ext <- st_bbox(lakes)
#deal with later, data too big with rasters


#load in eBird data
#this is the original species selection
# species_selection <- read.csv("data/species_selection.csv")
# species_selection <- species_selection[1:177,]
# species_selection1 <- species_selection %>%
#   filter(include == 1)

# ebd_data <- ebirdst_runs
# 
# species_data <- left_join(species_selection1,ebd_data, by = "common_name")
# 
# #test <- rast("D:/floating_solar/generated/osprey_max_values.tif")
# 
# #adding more species to be more comprehensive
# #do in batches so that you can keep code running while you sort through new species
# 
# ebd_with_trends <- ebirdst_runs %>%
#   filter(has_trends==T)
# ebd_trends_am <- ebd_with_trends %>%
#   filter(trends_region == "north_america")
# write.csv(ebd_trends_am, file = "data/species_selection_updated.csv")
#add species vulnerability codes manually?
###

#new species to run
updated_selection <- read.csv("data/species_selection_updated.csv")
updated_selection1 <- updated_selection %>%
  filter(northeast_america == 1)

complete <- list.files(path = "D:/floating_solar/generated/")
complete_codes <- str_extract(complete,"[^_]+")

new_codes <- setdiff(updated_selection1$species_code,complete_codes)

species_codes_round2 <- updated_selection1 %>%
  filter(species_code %in% new_codes)

species_data <- species_codes_round2

#sp <- "amewig" #loop here eventually

for(s in 1:length(species_data$species_code)){
  
#for(s in 1:3){
  
  skip_to_next <- FALSE
  
  sp <- species_data$species_code[s]

  # tryCatch(bird_data <- rast(paste0("D:/floating_solar/ebird/2022/",sp,"/weekly/",sp,"_abundance_median_3km_2022.tif")),
  #          error = function(e){skip_to_next <<- TRUE})
  
  tryCatch(bird_data <- rast(paste0("D:/big_data/eBird_FAC/2022/",sp,"/weekly/",sp,"_abundance_median_3km_2022.tif")),
           error = function(e){skip_to_next <<- TRUE})
  
  
  if(skip_to_next) {next}

  #crop to manageable size
  extent <- ext(-13000000, -4500000, 0, 9000000)
  bird_data_cr <- crop(bird_data, extent)

  #calculate proportion of total rel abundance that each cell represents
  #replace NAs with 0s
  bird_data[is.na(bird_data)] <- 0
  
  #sum of all raster cells in the extent
  sum <- global(bird_data_cr, "sum", na.rm=T)
  
  #for each week, divide each raster cell by the total amount to get the prop rel abd
  week_list <- list()
  
  for(i in 1:nlyr(bird_data_cr)){
    
    skip_to_next <- FALSE

    tryCatch(bird_prop_abd <- bird_data_cr[[i]]/sum$sum[i],
    error = function(e){skip_to_next <<- TRUE})

    if(skip_to_next) {next}
    
    week_list[[i]] <- bird_prop_abd
    
  }
  
  species_stack <- rast(week_list)
  
  #use max to get max across stacks (ie max abd index across weeks)
  max_values <- max(species_stack)
  
  writeRaster(max_values, file = paste0("D:/floating_solar/generated/",sp,"_max_values.tif"))

}


####summarising by lake####

#load in lake data
lakes <- read_sf("D:/floating_solar/Northeast_NHD_Alison")
#plot(st_geometry(lakes[1,]))
#plot(st_geometry(lakes))

#test on a couple lakes
lake1 <- lakes[4,]
lake_buffer <- st_buffer(lake1,5000)
#plot(st_geometry(lake_buffer))
#plot(st_geometry(lake1), add=TRUE)
lakes_vec <- vect(lake_buffer)
sp <- "amerob"
bird_data <- rast(paste0("D:/floating_solar/generated/",sp,"_max_values.tif"))
lakes_vec_pro <- project(lakes_vec, crs(bird_data))
plot(lakes_vec_pro)#look to check it worked
plot(bird_data, add=T) #this explodes :/
#lake_buffer <- st_buffer(lakes, 5000)
#plot(st_geometry(lake_buffer))

#actual starting point
#create buffer for all lakes and project
#lakes <- vect("D:/floating_solar/Northeast_NHD_Alison")
lakes <- read_sf("D:/floating_solar/Northeast_NHD_Alison")
lake_buffer <- st_buffer(lakes,5000)
#st_write(lake_buffer, file = "D:/floating_solar/generated/lake_5k_buffer.shp")
# lake_buffer <- load("D:/floating_solar/generated/lake_5k_buffer.RData") this didn't work, saved wrong
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
  
  lake_bird_data <- zonal(bird_data1, z = lakes_vec_pro, fun = "sum", na.rm=T)
  lake_bird_data$species_code <- rep(sp, length(lake_bird_data$max))
  lake_bird_data$Water_ID <- lakes$Water_ID

  save(lake_bird_data, file = paste0("D:/floating_solar/data_outputs/",sp,"_lake_abd_weight.RData"))
  
  }

