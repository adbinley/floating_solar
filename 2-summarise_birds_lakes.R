library(sf)
library(tidyverse)
library(ebirdst)
library(terra)

#load in lake data
lakes <- read_sf("D:/floating_solar/Northeast_NHD_Alison")
plot(st_geometry(lakes[1,]))


#load in eBird data
species_selection <- read.csv("data/species_selection.csv")
species_selection <- species_selection[1:177,]
species_selection <- species_selection %>%
  filter(include == 1)

ebd_data <- ebirdst_runs
species_data <- left_join(species_selection,ebd_data, by = "common_name")

#sp <- "amewig" #loop here eventually

#for(s in 1:length(species_data$species_code)){
  
for(s in 1:3){
  
  sp <- species_data$species_code[s]

  bird_data <- rast(paste0("D:/floating_solar/ebird/2022/",sp,"/weekly/",sp,"_abundance_median_3km_2022.tif"))

  #crop to manageable size
  extent <- ext(-10000000, -4500000, 4000000, 9000000)
  bird_data_cr <- crop(bird_data, extent)

  #calculate proportion of total rel abundance that each cell represents
  #replace NAs with 0s
  bird_data[is.na(bird_data)] <- 0
  
  #sum of all raster cells in the extent
  sum <- global(bird_data_cr, "sum", na.rm=T)
  
  #for each week, divide each raster cell by the total amount to get the prop rel abd
  week_list <- list()
  
  for(i in 1:nlyr(bird_data_cr)){

    bird_prop_abd <- bird_data_cr[[i]]/sum$sum[i]
    
    week_list[[i]] <- bird_prop_abd
    
  }
  
  species_stack <- rast(week_list)
  
  #use max to get max across stacks (ie max abd index across weeks)
  max_values <- max(species_stack)
  
  writeRaster(max_values, file = paste0("D:/floating_solar/generated/",sp,"_max_values.tif"))

}



