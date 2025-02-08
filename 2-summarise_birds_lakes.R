library(sf)
library(tidyverse)
library(ebirdst)
library(terra)
library(stringr)
library(viridis)

#load in lake data
lakes <- read_sf("D:/floating_solar/Northeast_NHD_Alison")
lakes <- read_sf("C:/Users/allis/OneDrive/Post-doc/big_data/floating_solar/Northeast_NHD_Alison")
#take a look
plot(st_geometry(lakes[1,]))
plot(st_geometry(lakes))
#extent?
lake_ext <- st_bbox(lakes)



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

#when loop fails (memory shortage) use this to figure out where to restart
complete <- list.files(path = "D:/floating_solar/generated/")
complete_codes <- str_extract(complete,"[^_]+")
#codes that still need to run
#new_codes <- setdiff(updated_selection1$species_code,complete_codes)
new_codes <- setdiff(data$species_code,complete_codes)

species_codes_round2 <- updated_selection1 %>%
  filter(species_code %in% new_codes)

species_data <- species_codes_round2

#coming back to troubleshoot some species with NAs
# na_species <- c("chiswi","chwwid","pecsan","yebcuc", "purmar", "veery",  "bkbcuc", "miskit", "amgplo", "baisan", "bubsan", "hudgod", "pursan", "sabgul", "uplsan",
# "bicthr", "boboli")

#loop through all species selected
#for(s in 1:length(species_data$species_code)){

#for(s in 1:length(na_species)){
  
#for(s in 1:3){ #test loop on a few species first to make sure its working

which(complete_codes=="wilfly")

for(s in 1:length(new_codes)){
  
  skip_to_next <- FALSE
  
  #sp <- species_data$species_code[s]
  #sp <- na_species[s]
  #sp <- complete_codes[s]
  sp <- new_codes[s]
  
  #trycatch added in case one raster fails to load - won't break loop
  #load raster into R for each species
  # tryCatch(bird_data <- rast(paste0("D:/big_data/eBird_FAC/2022/",sp,"/weekly/",sp,"_abundance_median_3km_2022.tif")),
  #          error = function(e){skip_to_next <<- TRUE})
  
  tryCatch(bird_data <- rast(paste0("D:/floating_solar/ebird/2022/",sp,"/weekly/",sp,"_abundance_median_3km_2022.tif")),
           error = function(e){skip_to_next <<- TRUE})
  
  
  if(skip_to_next) {next}

  #crop to manageable size
  #extent still captures entire western hemisphere
  #necessary so that the importance values are based on entire range at each point in the year
  extent <- ext(-13000000, -1500000, -5000000, 9000000)
  bird_data_cr <- crop(bird_data, extent)
  #plot(bird_data_cr[[20]])

  #calculate proportion of total rel abundance that each cell represents
  #replace NAs with 0s
  bird_data_cr[is.na(bird_data_cr)] <- 0
  
  #sum of all raster cells in the extent
  sum <- global(bird_data_cr, "sum", na.rm=T)
  
  #don't divide by zero...
  sum$sum[which(sum$sum ==0)] <- 0.1 #all cell values are zero anyways, so will  give zeros 
  
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
  #species_stack <- do.call(c,week_list)

  
  #use max to get max across stacks (ie max abd index across weeks)
  max_values <- max(species_stack)
  #max_raster <- app(species_stack, fun="max")

  
  writeRaster(max_values, file = paste0("D:/floating_solar/generated/",sp,"_max_values.tif"), overwrite=TRUE)
  
  rm(species_stack)
  rm(week_list)
  rm(max_values)
  rm(bird_data)
  rm(bird_data_cr)
  rm(bird_prop_abd)

}


#here - check species above worked ok, then rerun code below and all subsequent code
#for a handful of species, one or two cells were left at one point in the FAC and therefore there are values of 1 for relative importance
#because they are the only cells "in range"


####summarising by lake####

#load in lake data
# lakes <- read_sf("D:/floating_solar/Northeast_NHD_Alison")
# #plot(st_geometry(lakes[1,]))
# #plot(st_geometry(lakes))
# 
# #test on a couple lakes
# lake1 <- lakes[4,]
# lake_buffer <- st_buffer(lake1,5000)
# #plot(st_geometry(lake_buffer))
# #plot(st_geometry(lake1), add=TRUE)
# lakes_vec <- vect(lake_buffer)
# sp <- "amerob"
# bird_data <- rast(paste0("D:/floating_solar/generated/",sp,"_max_values.tif"))
# lakes_vec_pro <- project(lakes_vec, crs(bird_data))
# plot(lakes_vec_pro)#look to check it worked
# plot(bird_data, add=T) #this explodes :/
# #lake_buffer <- st_buffer(lakes, 5000)
# #plot(st_geometry(lake_buffer))

#actual starting point
#create buffer for all lakes and project
#lakes <- vect("D:/floating_solar/Northeast_NHD_Alison")
lakes <- read_sf("D:/floating_solar/Northeast_NHD_Alison")
sensitivty_analysis <- c(1000,3000,7000) #different buffer sizes, in meters
buf <- sensitivty_analysis[1]
lake_buffer <- st_buffer(lakes,buf)
#st_write(lake_buffer, file = "D:/floating_solar/generated/lake_5k_buffer.shp")
# lake_buffer <- load("D:/floating_solar/generated/lake_5k_buffer.RData") this didn't work, saved wrong
lakes_vec <- vect(lake_buffer)#create SpatVector
#writeVector(lakes_vec, file = "D:/floating_solar/generated/lake_5k_buffer_vec.gpkg")
#project to ebird data crs
sp <- "amerob" #species doesn't matter, using for projection
bird_data <- rast(paste0("D:/floating_solar/generated/",sp,"_max_values.tif"))
lakes_vec_pro <- project(lakes_vec, crs(bird_data))
#plot(lakes_vec_pro)#look to check it worked
#plot(bird_data, add=T)

complete <- list.files(path = "D:/floating_solar/generated/")
complete_codes <- str_extract(complete,"[^_]+")

#checking to see where computer crashed...
which(complete_codes=="gycthr")

#na species who have data but weird results
#na_species <- c("chiswi","chwwid","yebcuc", "purmar", "veery",  "bkbcuc", "miskit", "baisan", "pursan","uplsan",
 #               "bicthr", "boboli")

#sp <- na_species[1]

rm(sp)
rm(bird_data)

for(a in 136:length(complete_codes)){
#for(a in 1:length(na_species)){
  
  sp <- complete_codes[a]
  #sp <- na_species[a]

  bird_data <- rast(paste0("D:/floating_solar/generated/",sp,"_max_values.tif"))
  bird_data1 <- bird_data*10000 #transformed to make numbers nicer to deal with
  
  lake_bird_data <- zonal(bird_data1, z = lakes_vec_pro, fun = "sum", na.rm=T) 
  #instead of doing this, should probably change NAs to 0s
  #rerun eventually and fix, ignoring for now NOT DONE YET
  #only matters for mean values???
  lake_bird_data$species_code <- rep(sp, length(lake_bird_data$max))
  lake_bird_data$Water_ID <- lakes$Water_ID

  save(lake_bird_data, file = paste0("D:/floating_solar/data_outputs/",buf,"_",sp,"_lake_abd_weight.RData"))
  
}

species <- list.files(path = "D:/floating_solar/generated/")
species_codes <- str_extract(species,"[^_]+")

lake_biodiversity <- list()

for(s in 1:length(species_codes)){
  
  sp <- species_codes[s]

  #the sum of abd importance for each lake 
  load(paste0("D:/floating_solar/data_outputs/",sp,"_lake_abd_weight.RData")) 

  lake_biodiversity[[s]] <- lake_bird_data
  
}

lake_biodiversity_df_updated <- bind_rows(lake_biodiversity)
save(lake_biodiversity_df_updated, file = "D:/floating_solar/data_outputs/lake_ave_biodiversity_updated.RData")

load("D:/floating_solar/data_outputs/lake_ave_biodiversity_updated.RData")

lake_biodiversity_df <- lake_biodiversity_df_updated

#not by species
lake_bio_sum <- lake_biodiversity_df %>%
  group_by(Water_ID)%>%
  summarise(richness = sum(max>0,na.rm=T),
            importance = sum(max,na.rm=T))#does this last one make sense? add actual diversity metric?

rm(lake_biodiversity_df)
rm(lake_biodiversity)

lake_data <- st_drop_geometry(lakes)
rm(lakes)

all_data <- left_join(lake_data, lake_bio_sum)
only_selected_lakes <- all_data %>%
  filter(Suitabl_FP==1)
ggplot(only_selected_lakes, aes(year1_ener))+
  geom_histogram(binwidth = 10)

#how does energy production relate to avian richness and importance? 
#both with all lakes included, and with only those suitable for solar
plot(all_data$richness,all_data$year1_ener)
plot(only_selected_lakes$richness,only_selected_lakes$year1_ener)
plot(all_data$importance, all_data$year1_ener)
plot(only_selected_lakes$importance, only_selected_lakes$year1_ener)

#how do importance and richness compare between suitable and unsuitable lakes
all_data$Suitabl_FP <- as.factor(all_data$Suitabl_FP)
ggplot(all_data)+
  geom_boxplot(aes(x=Suitabl_FP, y = importance))
ggplot(all_data)+
  geom_boxplot(aes(x=Suitabl_FP, y = richness))

#importance and richness compared to sites suitable based on "biodiversity" scenario from TNC
ggplot(all_data)+
  geom_boxplot(aes(x=as.factor(Biodiversi), y = importance))
ggplot(all_data)+
  geom_boxplot(aes(x=as.factor(Biodiversi), y = richness))
ggplot(only_selected_lakes)+
  geom_boxplot(aes(x=as.factor(Biodiversi), y = importance))
ggplot(only_selected_lakes)+
  geom_boxplot(aes(x=as.factor(Biodiversi), y = richness))

#how do richness and importance relate to lake size?
ggplot(all_data)+
  geom_boxplot(aes(x=Water_Type, y = importance))
ggplot(all_data)+
  geom_boxplot(aes(x=Water_Type, y = richness))
plot(all_data$importance, all_data$Shape_Area)
plot(all_data$richness, all_data$Shape_Area)

#can we plot lakes on the map as points instead of polygons?
#size of point = energy, colour = birds
#try hex map instead?
points <- st_as_sf(all_data, coords=c('water_lon','water_lat'))
points <- st_set_crs(points, crs = 4326)
plot(st_geometry(points))
#too many still
points_select <- points %>%
  filter(Suitabl_FP == 1)
plot(st_geometry(points_select))

ggplot()+
  geom_point(data = only_selected_lakes, aes(x = water_lon, y = water_lat, color = richness, size = year1_ener))+
  scale_color_viridis()

#rank correlation between energy and birds in each state
#need to think on this a bit more - zeros dont work with rank
corr_all_data <- all_data %>%
  group_by(STATE_ABBR)%>%
  mutate(rank_corr = cor.test(x=year1_ener, y = richness))

corr_total <- cor.test(all_data$year1_ener,all_data$richness, method = 'spearman')
