library(sf)
library(tidyverse)
library(ebirdst)
library(terra)
library(stringr)
library(viridis)

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

#na species who have data but weird results
na_species <- c("chiswi","chwwid","yebcuc", "purmar", "veery",  "bkbcuc", "miskit", "baisan", "pursan","uplsan",
                "bicthr", "boboli")

rm(sp)
rm(bird_data)

#for(a in 1:length(complete_codes)){
for(a in 1:length(na_species)){
  
  #sp <- complete_codes[a]
  sp <- na_species[a]
  
  bird_data <- rast(paste0("D:/floating_solar/generated/",sp,"_max_values.tif"))
  bird_data1 <- bird_data*10000 #transformed to make numbers nicer to deal with
  
  bird_data1[is.na(bird_data1[])] <- 0 
  
  lake_bird_data_mean <- zonal(bird_data1, z = lakes_vec_pro, fun = "mean", na.rm=T) 

  lake_bird_data_mean$species_code <- rep(sp, length(lake_bird_data_mean$max))
  lake_bird_data_mean$Water_ID <- lakes$Water_ID
  
  save(lake_bird_data_mean, file = paste0("D:/floating_solar/data_outputs/",sp,"_lake_abd_MEAN_weight.RData"))
  
}

species <- list.files(path = "D:/floating_solar/generated/")
species_codes <- str_extract(species,"[^_]+")

mean_lake_biodiversity <- list()

for(s in 1:length(species_codes)){
  
  sp <- species_codes[s]
  
  #the mean of abd importance for each lake 
  load(paste0("D:/floating_solar/data_outputs/",sp,"_lake_abd_MEAN_weight.RData"))
  
  mean_lake_biodiversity[[s]] <- lake_bird_data_mean
  
}

mean_lake_biodiversity_df <- bind_rows(mean_lake_biodiversity)
save(mean_lake_biodiversity_df, file = "D:/floating_solar/data_outputs/mean_lake_ave_biodiversity_updated.RData")

#summarise by lake
mean_lake_bio_sum <- mean_lake_biodiversity_df %>%
  group_by(Water_ID)%>%
  summarise(mean_importance = mean(max))

rm(mean_lake_biodiversity_df)
rm(mean_lake_biodiversity)

#load in lake data
lakes <- read_sf("D:/floating_solar/Northeast_NHD_Alison")

lake_data <- st_drop_geometry(lakes)
rm(lakes)

all_data <- left_join(lake_data, mean_lake_bio_sum)
only_selected_lakes <- all_data %>%
  filter(Suitabl_FP==1)
# ggplot(only_selected_lakes, aes(year1_ener))+
#   geom_histogram(binwidth = 10)

#add to other importance and richness dataframe created in 2-summarise_birds_lakes
load("D:/floating_solar/data_outputs/lake_ave_biodiversity_updated.RData")

lake_biodiversity_df <- lake_biodiversity_df_updated

lake_bio_sum <- lake_biodiversity_df %>%
  group_by(Water_ID)%>%
  summarise(richness = sum(max>0,na.rm=T),
            sum_importance = sum(max,na.rm=T))#does this last one make sense? add actual diversity metric?

rm(lake_biodiversity_df)

all_data1 <- left_join(all_data, lake_bio_sum, by = "Water_ID")
save(all_data1, file = "D:/floating_solar/data_outputs/all_importance_data_updated.RData")

plot(all_data1$mean_importance,all_data1$sum_importance)
#strong relationship between the different importance metrics, but with sum outliers
#are the points that have high sum importance and low mean importance really large lakes?


#how does energy production relate to mean avian importance? 
#both with all lakes included, and with only those suitable for solar
plot(all_data$mean_importance,all_data$year1_ener)
plot(only_selected_lakes$mean_importance,only_selected_lakes$year1_ener)

#how does mean avian importance compare between suitable and unsuitable lakes
all_data$Suitabl_FP <- as.factor(all_data$Suitabl_FP)
ggplot(all_data)+
  geom_boxplot(aes(x=Suitabl_FP, y = mean_importance))


#mean importance compared to sites suitable based on "biodiversity" scenario from TNC
ggplot(all_data)+
  geom_boxplot(aes(x=as.factor(Biodiversi), y = mean_importance))
ggplot(only_selected_lakes)+
  geom_boxplot(aes(x=as.factor(Biodiversi), y = mean_importance))


#how does mean importance relate to lake size?
ggplot(all_data)+
  geom_boxplot(aes(x=Water_Type, y = mean_importance))
plot(all_data$mean_importance, all_data$Shape_Area)


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
  geom_point(data = only_selected_lakes, aes(x = water_lon, y = water_lat, color = year1_ener, size = mean_importance))+
  scale_color_viridis()

#rank correlation between energy and birds in each state
#need to think on this a bit more - zeros dont work with rank
# corr_all_data <- all_data %>%
#   group_by(STATE_ABBR)%>%
#   mutate(rank_corr = cor.test(x=year1_ener, y = richness))
# 
# corr_total <- cor.test(all_data$year1_ener,all_data$richness, method = 'spearman')