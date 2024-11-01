#scenarios
library(sf)
library(tidyverse)
library(ebirdst)
library(terra)
library(stringr)
library(viridis)

#map data
state_lines <- st_read("D:/maps/NA_politicalboundaries_shapefile/PoliticalBoundaries_Shapefile/NA_PoliticalDivisions/data/bound_p/boundaries_p_2021_v3.shp")

NE <- state_lines %>%
  filter(NAME_En %in% c("Maine","New Hampshire","Vermont","Massachusetts","Rhode Island","Connecticut","New York","Pennsylvania","New Jersey","Delaware","Maryland","West Virginia","Virginia","District of Columbia"))%>%
  st_geometry()

NE_pro <- st_transform(NE, crs = st_crs(4326))


#### richness ####
#1: overlap of species richness/diversity/importance and solar energy

load("D:/floating_solar/data_outputs/all_importance_data.RData")

only_selected_lakes <- all_data2 %>%
  filter(Suitabl_FP==1)%>%
  filter()

only_selected_lakes$bird_rank <- rank(-only_selected_lakes$richness, ties.method = "first")
only_selected_lakes$energy_scaled <- scale(only_selected_lakes$year1_ener)

plot_data <- only_selected_lakes %>%
  arrange((richness))%>%
  filter(richness!=0)

png("figures/richness_solar_overlay.png", height = 12, width = 12, units = "in",res=300)

ggplot()+
  geom_sf(data = NE_pro)+
  theme_void()+
  geom_point(data = plot_data, aes(x = water_lon, y = water_lat, color = richness, size = year1_ener))+
  scale_color_viridis()

dev.off()


ggplot(data = plot_data, aes(x = energy_scaled, y = richness)) +
  geom_point()

png("figures/mean_importance_solar_overlay.png", height = 12, width = 12, units = "in",res=300)

ggplot()+
  geom_sf(data = NE_pro)+
  theme_void()+
  geom_point(data = plot_data, aes(x = water_lon, y = water_lat, color = mean_importance, size = year1_ener))+
  scale_color_viridis()

dev.off()

png("figures/sum_importance_solar_overlay.png", height = 12, width = 12, units = "in",res=300)

ggplot()+
  geom_sf(data = NE_pro)+
  theme_void()+
  geom_point(data = plot_data, aes(x = water_lon, y = water_lat, color = sum_importance, size = year1_ener))+
  scale_color_viridis()

dev.off()

#do by state? could calculate the rank correlation between biodiversity importance and solar importance and summarise by state

#### biofouling ####

#2. look at relationship between the concentration of biofouling species and solar energy

load("D:/floating_solar/data_outputs/all_importance_data.RData")

only_selected_lakes <- all_data2 %>%
  filter(Suitabl_FP==1)%>%
  filter()

only_selected_lakes$bird_rank <- rank(-only_selected_lakes$sum_biofoul_risk, ties.method = "first")

plot_data <- only_selected_lakes %>%
  arrange((sum_biofoul_risk))%>%
  filter(sum_biofoul_risk!=0)

png("figures/sum_biofouling_solar_overlay.png", height = 12, width = 12, units = "in",res=300)

ggplot()+
  geom_sf(data = NE_pro)+
  theme_void()+
  geom_point(data = plot_data, aes(x = water_lon, y = water_lat, color = sum_biofoul_risk, size = year1_ener))+
  scale_color_viridis()

dev.off()

ggplot(data = plot_data, aes(x = (year1_ener), y = (sum_biofoul_risk)))+
  geom_point()



#### comparisons ####

ggplot(data = plot_data, aes(x = (richness), y = (sum_biofoul_risk)))+
  geom_point() #high species richness does not necessarily indicate high biofouling risk



