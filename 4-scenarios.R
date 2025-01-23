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
  mutate(vis_acuity_risk = 6 - axiallength_quantile)

#VI itself is not measured using exposure, that comes in later
#VI = (VA+FM)/2 * CCS * HS

lakes <- read_sf("D:/floating_solar/Northeast_NHD_Alison")

#suitable lakes only, and their percent coverage
lake_coverage <- lakes %>%
  select(c("Water_ID","FPV_Pct_co","Suitabl_FP"))%>%
  st_drop_geometry()%>%
  filter(Suitabl_FP ==1)

rm(lakes)

#calculate exposure to FVP
exposure_vector <- c()

for(s in 1:length(data$species_code)){

  sp <- data$species_code[s]
  
  load(paste0("D:/floating_solar/data_outputs/",sp,"_lake_abd_weight.RData")) #this is sum, mean is also available
  
  lake_bird_data1 <- lake_bird_data %>%
    filter(Water_ID %in% lake_coverage$Water_ID)
  
  lake_bird_data2 <- left_join(lake_bird_data1, lake_coverage)
  
  # #calculating exposure based on max abundance at each waterbody and solar coverage
  lake_bird_data2 <- lake_bird_data2 %>%
    #multiply importance at each lake by the proportion of the lake to be covered
    mutate(exposure = max*FPV_Pct_co)
  
  sum_exposure <- sum(lake_bird_data2$exposure) #sum of exposure across study range - will eventually be scaled based on other species values

  exposure_vector <- c(exposure_vector, sum_exposure)
  
}

data$exposure <- exposure_vector

#scale between 1 and 5, but leave zeros as zeros
# data <- data %>%
# mutate(exposure_scaled = ifelse(exposure != 0,
#                                 scales::rescale(exposure, to = c(1,5)),
#                                 0))

#new - scaling exposure between 0 and 1
data <- data %>%
mutate(exposure_scaled_2 = scales::rescale(exposure, to = c(0,1)))

data$risk <- data$exposure_scaled_2*data$VI

# #calculate new VI that accounts for exposure
# data <- data %>%
#   mutate(VI = ((vis_acuity_risk+wingloading_quantile)/2)*CCS.max*((habitat_score+exposure_scaled)/2))
#VI will now be calculated without adding exposure
#exposure added later

#reverting to original calculation of VI
#does not account for exposure yet
data$VI = ((data$wingloading_quantile + data$vis_acuity_risk)/2) * data$CCS_quantile * data$habitat_score


save(data, file = "data_outputs/final_analysis_data.RData")
write.csv(data, file = "data_outputs/final_analysis_data.csv")



#### richness ####
#2: overlap of species richness/diversity/importance and solar energy

load("D:/floating_solar/data_outputs/all_importance_data_updated.RData") 

only_selected_lakes <- all_data1 %>%
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

#### VI for each lake ####

load("data_outputs/final_analysis_data.RData")

lakes <- read_sf("D:/floating_solar/Northeast_NHD_Alison")

#suitable lakes only, and their percent coverage
lake_coverage <- lakes %>%
  select(c("Water_ID","FPV_Pct_co","Suitabl_FP"))%>%
  st_drop_geometry()%>%
  filter(Suitabl_FP ==1)

rm(lakes)

lake_VI_df <- data.frame(Water_ID = lake_bird_data$Water_ID)%>%
  filter(Water_ID %in% lake_coverage$Water_ID)

#which(data$species_code == "gbbgul")

for(s in 1:length(data$species_code)){
  
  sp <- data$species_code[s]
  
  load(paste0("D:/floating_solar/data_outputs/",sp,"_lake_abd_weight.RData")) #this is sum, mean is also available
  
  lake_bird_data1 <- lake_bird_data %>%
    filter(Water_ID %in% lake_coverage$Water_ID)
  
  lake_bird_data2 <- left_join(lake_bird_data1, lake_coverage)
  
  #calculating exposure based on max abundance at each waterbody and solar coverage
  #this gives us a scaled value that represents the range of exposures for each species
  lake_bird_data2 <- lake_bird_data2 %>%
    #multiply importance at each lake by the proportion of the lake to be covered
    #note here the scaling is for each species
    mutate(exposure = max*FPV_Pct_co)%>%
    mutate(exposure_scaled_2 = scales::rescale(exposure, to = c(0,1)))
    # mutate(exposure_scaled = ifelse(exposure != 0,
    #                                 scales::rescale(exposure, to = c(1,5)),
    #                                 0))
  
  sp_data <- data %>%
    filter(species_code == sp)
  
  #calculating the VI for each species at each lake, using the specific exposure values for each lake
  #lake_bird_data2$lake_VI = ((sp_data$vis_acuity_risk+sp_data$wingloading_quantile)/2)*sp_data$CCS.max*((sp_data$habitat_score+lake_bird_data2$exposure_scaled)/2)
  lake_bird_data2$risk = sp_data$VI * lake_bird_data2$exposure_scaled_2
  
  sp_df <- lake_bird_data2 %>%
    select(c("Water_ID","risk"))
  
  colnames(sp_df) <- c("Water_ID",sp)
  
  
  lake_VI_df <- left_join(lake_VI_df,sp_df)
  
}

save(lake_VI_df, file="data_outputs/lake_VI_df.RData")
load("data_outputs/lake_VI_df.RData")

lake_VI_df$mean_risk = rowMeans(lake_VI_df[,2:ncol(lake_VI_df)])
lake_VI_df$sum_risk = rowSums(lake_VI_df[,2:ncol(lake_VI_df)])

lake_risk_df <- lake_VI_df %>%
  select(c("Water_ID","mean_risk","sum_risk"))

save(lake_risk_df, file = "data_outputs/lake_risk_df.RData")

#### start here ####
load("data_outputs/lake_risk_df.RData")

load("D:/floating_solar/data_outputs/all_importance_data_updated.RData") 

only_selected_lakes <- all_data1 %>%
  filter(Suitabl_FP==1)%>%
  filter()

all_data2 <- left_join(only_selected_lakes,lake_risk_df)

all_data2$bird_rank <- rank(-all_data2$mean_risk, ties.method = "first")
all_data2$energy_scaled <- scale(all_data2$year1_ener)

plot_data <- all_data2 %>%
  arrange((mean_risk))

png("figures/risk_solar_overlay.png", height = 12, width = 12, units = "in",res=300)

ggplot()+
  geom_sf(data = NE_pro)+
  theme_void()+
  geom_point(data = plot_data, aes(x = water_lon, y = water_lat, color = mean_risk, size = year1_ener))+
  scale_color_viridis()

dev.off()


ggplot(data = plot_data, aes(x = energy_scaled, y = mean_risk)) +
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

#### biofouling ####

#2. look at relationship between the concentration of biofouling species and solar energy

#load("D:/floating_solar/data_outputs/all_importance_data.RData")
load("data_outputs/final_biofouling_data.RData")

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




# energy scenarios --------------------------------------------------------
library(sf)
library(tidyverse)
library(prioritizr)
library(tmap)



#simple approach where lakes are selected based on two criteria:
#first, starting with the lake that has the best energy/cost ratio - plot energy per unit cost
#second, starting with the lake that minimizes overall biodiversity risk - plot energy per unit cost
load("data_outputs/lake_risk_df.RData")
load("data/suitable_lakes.RData")


energy_data <- left_join(lakes1,lake_risk_df)

energy_data <- energy_data %>%
  mutate(energy_cost_ratio = year1_ener/fpv_ha)%>%
  arrange(desc(energy_cost_ratio))

energy_data1 <- data.frame(no_lakes_selected = 1,energy_total = energy_data[1,]$year1_ener, cost = energy_data[1,]$fpv_ha)

no_lakes_selected = 1

for(i in 1:length(energy_data$Water_ID)){
  
  no_lakes_selected <- no_lakes_selected+1
  
  lake_sel <- energy_data[1:no_lakes_selected,]
  
  energy_total <- sum(lake_sel$year1_ener)
  
  cost <- sum(lake_sel$fpv_ha)
  
  df <- data.frame(no_lakes_selected = no_lakes_selected, energy_total = energy_total, cost = cost)
  
  energy_data1 <- rbind(energy_data1,df)
  
  
}

energy_data1 <- energy_data1[-16621,]
  
save(energy_data1, file = "data_outputs/selected_lakes_energy_priority.RData")

energy_data1 <- energy_data %>%
  arrange(desc(energy_cost_ratio)) #arrange with highest energy/cost ratio first

energy_data1$no_lakes_selected <- 1:length(energy_data1$Water_ID)
energy_data1$energy_total <- cumsum(energy_data1$year1_ener)
energy_data1$cost <- cumsum(energy_data1$fpv_ha)
energy_data1$cumulative_risk <- cumsum(energy_data1$sum_risk)

#this approach is a million times faster than above

bio_data <- energy_data %>%
  arrange(sum_risk)

bio_data$no_lakes_selected = 1:length(bio_data$Water_ID)
bio_data$energy_total = cumsum(bio_data$year1_ener)
bio_data$cost = cumsum(bio_data$fpv_ha)
bio_data$cumulative_risk <- cumsum(bio_data$sum_risk)

energy_data1$scenario <- rep("energy_priority",length(energy_data1$no_lakes_selected))

bio_data$scenario <- rep("avian_priority",length(bio_data$Water_ID))

scen_data <- rbind(energy_data1,bio_data)

png("figures/2_scenario_energy_total.png", height = 10, width = 10, units="in", res = 300)

ggplot(scen_data, aes(x=no_lakes_selected,y=energy_total, col=scenario))+
  geom_line()+
  theme_classic(base_size = 22)

dev.off()


difference_energy_production <- scen_data$energy_total[scen_data$scenario=="energy_priority"] - scen_data$energy_total[scen_data$scenario=="avian_priority"]
max(difference_energy_production) #23711.21


png("figures/2_scenario_cum_risk.png", height = 10, width = 10, units="in", res = 300)

ggplot(scen_data, aes(x=no_lakes_selected,y=cumulative_risk, col=scenario))+
  geom_line()+
  theme_classic(base_size = 22)

dev.off()


#most basic approach - removing priority avian biodiversity lakes from the analysis

sum(lakes1$Climate_Cr)#16620
sum(lakes1$Social_B_1)#14318, 2302 lakes removed
sum(lakes1$Biodiversi)#8777, 7843 lakes removed
sum(lakes1$Precaution)#7477, 9143 lakes removed

mean(energy_data$mean_risk)

energy_data$avian_scenario <- ifelse(energy_data$sum_risk > median(energy_data$sum_risk), 0,1)
sum(energy_data$avian_scenario)#8310 lakes removed

#which lakes are selected based on all three criteria?
#which lakes meet the other two criteria, but not the avian biodiversity?
energy_data$new_precaution <- ifelse(energy_data$Social_B_1 ==1 & energy_data$Biodiversi==1 & energy_data$avian_scenario ==1,"all",
                                     ifelse(energy_data$Social_B_1 + energy_data$Biodiversi ==1 & energy_data$avian_scenario ==1, "birds +",
                                            ifelse(energy_data$Social_B_1 + energy_data$Biodiversi >=1 & energy_data$avian_scenario ==0,"non-bird value",
                                                   ifelse(energy_data$Social_B_1 + energy_data$Biodiversi == 0 & energy_data$avian_scenario ==1, "birds_only","none"))))
#look
scenarios <- energy_data %>%
  select(c("Water_ID","Social_B_1","Biodiversi","avian_scenario","new_precaution"))%>%
  st_drop_geometry()

lakes <- read_sf("D:/floating_solar/Northeast_NHD_Alison")
lakes1 <- lakes %>%
  filter(Suitabl_FP ==1)

scenarios_sf <- left_join(lakes1,scenarios)

NE_pro1 <- st_transform(NE, st_crs(scenarios_sf))

png("figures/all_scenarios1.png", height = 10, width = 10, units="in", res = 300)

plot(st_geometry(NE_pro1))
plot(scenarios_sf[,"new_precaution"], border=NA, pal = c("blue","orange","green","yellow","black"), main=NULL, add=T)

dev.off()

# plot map of prioritization
plot(
  s_min_birds[, "map_1"], pal = c("grey90","purple"),
  main = NULL, key.pos = 1
)

#bar chart or sexier equivalent? probably shows data better

#how much do the avian and TNC scenarios overlap?





#load("C:/Users/allis/OneDrive/Post-doc/floating_solar/floating_solar/data_outputs/final_analysis_data.RData")
#lakes <- read_sf("C:/Users/allis/OneDrive/Post-doc/big_data/floating_solar/Northeast_NHD_Alison")
#lakes <- read_sf("D:/floating_solar/Northeast_NHD_Alison")
# load("data_outputs/lake_risk_df.RData")
# 
# lakes1 <- lakes %>%
#   filter(Suitabl_FP==1)
# rm(lakes)
# st_write(lakes1, dsn ="data/suitable_lakes.shp")

#load("data/suitable_lakes.RData")

#maximizing energy production with no constraints

#using original dataset, and including the hectares of FPV as the "cost"
p_max_energy <- prioritizr::problem(x=lakes1, features="year1_ener", cost_column = "fpv_ha")%>%
  add_relative_targets(1)%>% #aiming to maximize energy production
  add_min_shortfall_objective(budget = 26000)%>% #"budget" of max FVP available, approx. 1/5 of total climate crisis coverage
  add_binary_decisions()%>%
  add_default_solver(gap=0, verbose=F)

print(p_max_energy)

s_max_energy <- solve(p_max_energy, force = T)
sum(s_max_energy$solution_1) #6451 of the total 16620 suitable lakes were selected

selection_max_energy <- s_max_energy %>%
  filter(solution_1==1)
sum(selection_max_energy$year1_ener)#generating 33587.63 units energy


#maximizing energy production with avian biodiversity constraint

lakes2 <- left_join(lakes1,lake_risk_df)
lakes2$mean_risk_scaled <- scale(lakes2$mean_risk)[,1]
lakes2$energy_scaled <- scale(lakes2$year1_ener)[,1]

penalty_data <- lakes2 %>%
  select("Water_ID","Climate_Cr","Biodiversi","Social_B_1","mean_risk","year1_ener")%>%
  st_drop_geometry()


p_min_birds <- prioritizr::problem(x=lakes2, features ="year1_ener", cost_column = "fpv_ha")%>%
  add_relative_targets(1)%>% #aiming to maximize energy production
  add_min_shortfall_objective(budget = 26000)%>% #"budget" of max FVP available, approx. 1/5 of total climate crisis coverage
  add_linear_penalties(penalty = 1, data = "mean_risk")%>% #(penalty="mean_risk", data=penalty_data)%>%
  add_binary_decisions()%>%
  add_default_solver(gap=0, verbose=F)

print(p_min_birds)

s_min_birds <- solve(p_min_birds, force = T)

s_min_birds$map_1 <- case_when(
  s_min_birds$solution_1 > 0.5 ~ "priority",
  TRUE ~ "other"
)

# plot map of prioritization
plot(
  s_min_birds[, "map_1"], pal = c("grey90","purple"),
  main = NULL, key.pos = 1
)

#different way of looking at it
s_min_birds$solution_1 <- as.factor(s_min_birds$solution_1)

ggplot(s_min_birds, aes(x=mean_risk,y=log(year1_ener), col=solution_1))+
  geom_point()+
  theme_classic()
