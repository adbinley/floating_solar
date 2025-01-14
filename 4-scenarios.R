#scenarios
library(sf)
library(tidyverse)
library(ebirdst)
library(terra)
library(stringr)
library(viridis)
library(scales)

#map data
state_lines <- st_read("D:/maps/NA_politicalboundaries_shapefile/PoliticalBoundaries_Shapefile/NA_PoliticalDivisions/data/bound_p/boundaries_p_2021_v3.shp")

NE <- state_lines %>%
  filter(NAME_En %in% c("Maine","New Hampshire","Vermont","Massachusetts","Rhode Island","Connecticut","New York","Pennsylvania","New Jersey","Delaware","Maryland","West Virginia","Virginia","District of Columbia"))%>%
  st_geometry()

NE_pro <- st_transform(NE, crs = st_crs(4326))

#### species vulnerability ####
#1.: which species have the highest vulnerability to floating solar based on our equation?

#this is done, skip to bottom now for final data
data <- read.csv("data/final_analysis_data.csv")

#first reverse the quantiles for axial length such that high values been smaller axial length and greater risk
data <- data %>%
  mutate(vis_acuity_risk = 6 - axiallength_quantile)#%>%
  #mutate(VI = ((vis_acuity_risk+wingloading_quantile)/2)*CCS.max*habitat_score)

#what about when weighted by their relative abundance in this region
# for each species, take the summed relative importance across all lakes identified as appropriate for floating solar
# add this to the risk part of the assessment: VI = (VA+FM)/2 * CCS + (HS+(exp))/2
#exposure is also multiplied by percent coverage of solar

lakes <- read_sf("D:/floating_solar/Northeast_NHD_Alison")

#suitable lakes only, and their percent coverage
lake_coverage <- lakes %>%
  select(c("Water_ID","FPV_Pct_co","Suitabl_FP"))%>%
  st_drop_geometry()%>%
  filter(Suitabl_FP ==1)

rm(lakes)

exposure_vector <- c()

for(s in 1:length(data$species_code)){

  sp <- data$species_code[s]
  
  load(paste0("D:/floating_solar/data_outputs/",sp,"_lake_abd_weight.RData")) #this is sum, mean is also available
  
  lake_bird_data1 <- lake_bird_data %>%
    filter(Water_ID %in% lake_coverage$Water_ID)
  
  lake_bird_data2 <- left_join(lake_bird_data1, lake_coverage)
  
  #calculating exposure based on max abundance at each waterbody and solar coverage
  lake_bird_data2 <- lake_bird_data2 %>%
    #multiply importance at each lake by the proportion of the lake to be covered
    mutate(exposure = max*FPV_Pct_co)
  
  sum_exposure <- sum(lake_bird_data2$exposure) #sum of exposure across study range - will eventually be scaled based on other species values
  
  exposure_vector <- c(exposure_vector, sum_exposure)
  
}

data$exposure <- exposure_vector

#scale between 1 and 5, but leave zeros as zeros
data <- data %>%
mutate(exposure_scaled = ifelse(exposure != 0,
                                scales::rescale(exposure, to = c(1,5)),
                                0))

#calculate new VI that accounts for exposure
data <- data %>%
  mutate(VI = ((vis_acuity_risk+wingloading_quantile)/2)*CCS.max*((habitat_score+exposure_scaled)/2))


save(data, file = "data_outputs/final_analysis_data.RData")


load("data_outputs/final_analysis_data.RData")
write.csv(data, file = "data_outputs/final_analysis_data.csv")

#### richness ####
#2: overlap of species richness/diversity/importance and solar energy

load("D:/floating_solar/data_outputs/all_importance_data_updated.RData") 

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



