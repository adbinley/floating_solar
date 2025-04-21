#scenarios
library(sf)
library(tidyverse)
library(ebirdst)
library(terra)
library(stringr)
library(viridis)
library(scales)
library(ggsci)

#map data
state_lines <- st_read("A:/maps/NA_politicalboundaries_shapefile/PoliticalBoundaries_Shapefile/NA_PoliticalDivisions/data/bound_p/boundaries_p_2021_v3.shp")

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

#calculate exposure to FPV
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

#manually removing species out of range and readjusting HS based on Song et al. 2024
data <- read.csv("data_outputs/final_analysis_data_n291.csv")

#recalculate VI
data$VI = ((data$wingloading_quantile + data$vis_acuity_risk)/2) * data$CCS_quantile * data$habitat_score
write.csv(data, "data_outputs/final_analysis_data_n291.csv")

####quantify uncertainty in HS####

data <- read.csv("data_outputs/final_analysis_data_n291.csv")

min_range <- function(x){
  ifelse(x==1,1,x-1)
}

max_range <- function(x){
  ifelse(x==5,5,x+1)
}

data$median_VI_un <- rep(NA,length(data$species_code))
data$hdi_l <- rep(NA,length(data$species_code))
data$hdi_u <- rep(NA,length(data$species_code))

for(s in 1:length(data$species_code)){
  
  df <- data[s,]
  
  set.seed(1)
  VI_unc <- ((runif(1000, min=min_range(df$wingloading_quantile),max=max_range(df$wingloading_quantile)) + 
              runif(1000,min = min_range(df$vis_acuity_risk),max=max_range(df$vis_acuity_risk)))/2) * 
              runif(1000, min = min_range(df$CCS_quantile),max=max_range(df$CCS_quantile)) * 
              runif(1000, min = 1, max = 5)

  data$median_VI_un[s] <- median(VI_unc)
  data$hdi_l[s] <- HPDI(VI_unc,0.9)[1]
  data$hdi_u[s] <- HPDI(VI_unc,0.9)[2]

}

#plot VIs and uncertainty pdf

data_arr <- data %>%
  arrange(desc(VI))


rows_per_batch <- 75

pdf("figures/VIs_with_uncertainty.pdf")

for(s in seq(1, 291, by = rows_per_batch)){
  
  subset_df <- data_arr[s:min(s + rows_per_batch - 1, nrow(data_arr)), ]
  #subset_df <- data_arr[1:97,]
  
  VI_plot <- ggplot(subset_df, aes(x=reorder(as.factor(common_name.x),VI),y=VI))+
    geom_point()+
    geom_hline(yintercept = mean(data_arr$VI), col="blue",linetype = "dashed")+
    geom_point(aes(x=reorder(as.factor(common_name.x),VI), y=median_VI_un), col="grey",alpha = 0.3)+
    geom_errorbar(aes(x=as.factor(common_name.x),ymin=hdi_l,ymax=hdi_u), col="grey", alpha = 0.3)+
    coord_flip()+
    #ylab("mean slope (90% CI)")+
    #ggtitle("Ecoregion Average Slope")+
    theme_classic(base_size = 8)
  
  print(VI_plot)
  
}

dev.off()

#spider/radar charts
library(fmsb)

radar_data <- data %>%
  filter(species_code %in% c("horgre","wessan","osprey","mallar3","leabit","marwre"))%>%
  select(c("common_name.x","VI","vis_acuity_risk","CCS_quantile","wingloading_quantile","habitat_score"))

min <- c("min",1,1,1,1,1)
max <- c("max",5,5,5,5,5)

radar_data <- rbind(max,min, radar_data)

colnames(radar_data) <- c("species","VI","VA","CCS","WL","HS")
radar_data[,3:6] <- sapply(radar_data[,3:6],as.numeric)


# Define colors and titles
colors <- c("#00AFBB", "#E7B800", "#FC4E07","#660000","#003300","#000066")
titles <- c("Horned Grebe (VI: 125.0)","Least Bittern (VI: 26.4)","Mallard (VI: 14.0)","Marsh Wren (VI: 5)","Osprey (VI: 7.5)","Western Sandpiper (VI: 56.5)")

# Reduce plot margin using par()
# Split the screen in 3 parts
op <- par(mar = c(1, 1, 1, 1))
par(mfrow = c(2,3))

# Create the radar chart
create_beautiful_radarchart <- function(data, color = "#00AFBB", 
                                        vlabels = colnames(data), vlcex = 1.2,
                                        caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "grey", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}

png("figures/radar_plot.png",height = 9, width = 12, units = "in",res=300)

op <- par(mar = c(1, 1, 1, 1))
par(mfrow = c(2,3))

for(i in 1:6){
  create_beautiful_radarchart(
    data = radar_data[c(1, 2, i+2), c(3:6)], caxislabels = c(1, 2, 3, 4, 5),
    color = colors[i], title = titles[i]
  )
}
par(op)

dev.off()

mean(data_arr$VI)
sd(data_arr$VI)
top10 <- data_arr[1:10,]
mean(top10$habitat_score)
sd(top10$habitat_score)/sqrt(10)
mean(top10$vis_acuity_risk)
sd(top10$vis_acuity_risk)/sqrt(10)
mean(top10$wingloading_quantile)
sd(top10$wingloading_quantile)/sqrt(10)
mean(top10$CCS_quantile)
sd(top10$CCS_quantile)/sqrt(10)

par(mfrow = c(1,1))

#### richness ####
#2: overlap of species richness/diversity/importance and solar energy

load("A:/floating_solar/data_outputs/all_importance_data_updated.RData") 

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

#load("data_outputs/final_analysis_data.RData")
data <- read.csv("data_outputs/final_analysis_data_n291.csv")


lakes <- read_sf("D:/floating_solar/Northeast_NHD_Alison")
#lakes <- read_sf("C:/Users/adb299/OneDrive/Post-doc/big_data/floating_solar/Northeast_NHD_Alison")

#suitable lakes only, and their percent coverage
lake_coverage <- lakes %>%
  select(c("Water_ID","FPV_Pct_co","Suitabl_FP"))%>%
  st_drop_geometry()%>%
  filter(Suitabl_FP ==1)

rm(lakes)

# lake_VI_df <- data.frame(Water_ID = lake_bird_data$Water_ID)%>%
#   filter(Water_ID %in% lake_coverage$Water_ID)

lake_VI_df <- data.frame(Water_ID = lake_coverage$Water_ID)
lake_exp_df <- data.frame(Water_ID = lake_coverage$Water_ID)

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
    mutate(importance = max)%>%
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
  
  sp_exp_df <- lake_bird_data2 %>%
    select(c("Water_ID","exposure"))
  
  colnames(sp_df) <- c("Water_ID",sp)
  
  colnames(sp_exp_df)<- c("Water_ID",sp)
  
  
  lake_VI_df <- left_join(lake_VI_df,sp_df)
  
  lake_exp_df <- left_join(lake_exp_df,sp_exp_df)
  
}
save(lake_exp_df, file = "data_outputs/lake_exp_df.RData")
save(lake_VI_df, file="data_outputs/lake_VI_df.RData")

load("data_outputs/lake_VI_df.RData")

lake_risk_df <- data.frame(Water_ID = lake_VI_df$Water_ID,
                      mean_risk = rowMeans(lake_VI_df[,2:ncol(lake_VI_df)]),
                      sum_risk = rowSums(lake_VI_df[,2:ncol(lake_VI_df)]))

lake_risk_df$mean_risk_scaled <- scale(risk_df$mean_risk)[,1]

save(lake_risk_df, file = "data_outputs/lake_risk_df.RData")

#alternate approach

load("data_outputs/lake_exp_df.RData")
data <- read.csv("data_outputs/final_analysis_data_n291.csv")

w_mean_risk <- c()

for(l in 1:length(lake_exp_df$Water_ID)){

  # lake_risk_df_weighted <- data.frame(Water_ID = lake_VI_df$Water_ID,
  #                                     w_mean_risk = weighted.mean(lake_VI_df[,2:ncol(lake_VI_df)], data$VI))
  
  d <- lake_exp_df[l,2:ncol(lake_exp_df)]
  w_risk <- weighted.mean(d,data$VI)
  
  w_mean_risk <- c(w_mean_risk,w_risk)
  
}


save(lake_risk_df_weighted, file = "data_outputs/lake_risk_df_weighted.RData")

#### start here ####
load("data_outputs/lake_risk_df.RData")
load("data_outputs/lake_risk_df_weighted.RData")

lake_risk_df <- lake_risk_df_weighted

data <- read.csv("data_outputs/final_analysis_data_n291.csv")

lakes <- read_sf("A:/floating_solar/Northeast_NHD_Alison")
lakes <- read_sf("D:/floating_solar/Northeast_NHD_Alison")

selected_lakes <- lakes %>%
  filter(Suitabl_FP ==1)

#load("A:/floating_solar/data_outputs/all_importance_data_updated.RData") 

# only_selected_lakes <- all_data1 %>%
#   filter(Suitabl_FP==1)%>%
#   filter()

all_data2 <- left_join(selected_lakes,lake_risk_df)

all_data2$bird_rank <- rank(-all_data2$w_mean_risk, ties.method = "first")
all_data2$energy_scaled <- scale(all_data2$year1_ener)[,1]

cor(all_data2$energy_scaled,all_data2$mean_risk_scaled)

plot_data <- all_data2 %>%
  arrange((mean_risk))
# 
# plot_data2 <- plot_data %>%
#   filter(mean_risk_scaled>0)

png("figures/risk_solar_overlay1.png", height = 9, width = 11, units = "in",res=300)

ggplot()+
  geom_sf(data = NE_pro)+
  theme_classic(base_size = 15)+
  geom_point(data = plot_data, aes(x = water_lon, y = water_lat, col = mean_risk_scaled, size = energy_scaled))+
  scale_color_viridis(option="inferno",limits = c(-3,19))+
  ylab("")+
  xlab("")+
  scale_size(guide="none")
  #scale_color_gsea(reverse = TRUE, limits = c(-19,19))

dev.off()


plot_data2 <- plot_data %>%
  arrange(-energy_scaled)

png("figures/energy_risk_corr.png", height = 6, width = 8, units = "in",res=300)

ggplot(plot_data, aes(x=energy_scaled, y=mean_risk_scaled) ) +
  geom_hex() +
  scale_fill_continuous(type = "viridis") +
  theme_bw(base_size = 20)+
  xlab("energy")+
  ylab("risk")

dev.off()


ggplot(data = plot_data, aes(x = energy_scaled, y = mean_risk_scaled)) +
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
# load("data_outputs/final_biofouling_data.RData")

# only_selected_lakes <- all_data2 %>%
#   filter(Suitabl_FP==1)%>%
#   filter()

data <- read.csv("data_outputs/final_analysis_data_n291.csv")

biofoul_data <- data %>%
  filter(group %in% c("waterbird","shorebird","waterfowl"))

#scale risk between 0 and 1, as percentage of max
biofoul_data$biofoul_risk <- rescale(biofoul_data$predicted_defecation_rate, to=c(0,1))
lakes <- read_sf("A:/floating_solar/Northeast_NHD_Alison")

#suitable lakes only, and their percent coverage
lake_coverage <- lakes %>%
  select(c("Water_ID","FPV_Pct_co","Suitabl_FP"))%>%
  st_drop_geometry()%>%
  filter(Suitabl_FP ==1)

lake_biofoul_df <- data.frame(Water_ID = lake_coverage$Water_ID)

#which(data$species_code == "gbbgul")

for(s in 1:length(biofoul_data$species_code)){
  
  sp <- biofoul_data$species_code[s]
  
  load(paste0("A:/floating_solar/data_outputs/",sp,"_lake_abd_weight.RData")) #this is sum, mean is also available
  
  lake_bird_data1 <- lake_bird_data %>%
    filter(Water_ID %in% lake_coverage$Water_ID)
  
  lake_bird_data2 <- left_join(lake_bird_data1, lake_coverage)
  
  #calculating exposure based on max abundance at each waterbody and solar coverage
  #this gives us a scaled value that represents the range of exposures for each species
  lake_bird_data2 <- lake_bird_data2 %>%
    #multiply importance at each lake by the proportion of the lake to be covered
    #note here the scaling is for each species
    mutate(importance = max)%>%
    mutate(exposure = max*FPV_Pct_co)%>%
    mutate(exposure_scaled_2 = scales::rescale(exposure, to = c(0,1)))
  # mutate(exposure_scaled = ifelse(exposure != 0,
  #                                 scales::rescale(exposure, to = c(1,5)),
  #                                 0))
  
  sp_data <- biofoul_data %>%
    filter(species_code == sp)
  
  #calculating the biofoul risk for each species at each lake, using the specific exposure values for each lake
  lake_bird_data2$biofoul_risk = sp_data$biofoul_risk * lake_bird_data2$exposure_scaled_2
  
  sp_df <- lake_bird_data2 %>%
    select(c("Water_ID","biofoul_risk"))
  
  colnames(sp_df) <- c("Water_ID",sp)
  
  
  lake_biofoul_df <- left_join(lake_biofoul_df,sp_df)
  
}

save(lake_biofoul_df, file="data_outputs/lake_biofoul_df.RData")
load("data_outputs/lake_biofoul_df.RData")


lake_biofoul_df1 <- data.frame(Water_ID = lake_biofoul_df$Water_ID,
                           mean_risk = rowMeans(lake_biofoul_df[,2:ncol(lake_biofoul_df)]),
                           sum_risk = rowSums(lake_biofoul_df[,2:ncol(lake_biofoul_df)]))

lake_biofoul_df1$mean_risk_scaled <- scale(lake_biofoul_df1$mean_risk)[,1]

colnames(lake_biofoul_df1) <- c("Water_ID","mean_biof_risk","sum_biof_risk","mean_biof_risk_scaled")

save(lake_biofoul_df1, file = "data_outputs/lake_biofoul_risk_df.RData")

####start here####

load("data_outputs/lake_biofoul_risk_df.RData")

data <- read.csv("data_outputs/final_analysis_data_n291.csv")

lakes <- read_sf("A:/floating_solar/Northeast_NHD_Alison")

selected_lakes <- lakes %>%
  filter(Suitabl_FP ==1)

all_data_biof <- left_join(selected_lakes,lake_biofoul_df1)

all_data_biof$bird_rank <- rank(-all_data_biof$mean_biof_risk, ties.method = "first")
all_data_biof$energy_scaled <- scale(all_data_biof$year1_ener)[,1]

plot_data <- all_data_biof %>%
  arrange((mean_biof_risk))

all_data_ranked <- data.frame(Water_ID = plot_data$Water_ID,
                              bird_rank = rank(all_data_biof$mea))

png("figures/biofouling_risk.png", height = 6, width = 8, units = "in",res=300)

ggplot()+
  geom_sf(data = NE_pro)+
  theme_classic(base_size = 15)+
  geom_point(data = plot_data, aes(x = water_lon, y = water_lat, col = mean_biof_risk_scaled, size = energy_scaled))+
  scale_color_viridis(option="inferno",limits = c(-5,30))+
  ylab("")+
  xlab("")+
  scale_size(guide="none")
#scale_color_gsea(reverse = TRUE, limits = c(-19,19))

dev.off()

#is biofouling risk simply high where risk to avian biodiversity is high?

biof_risk <- data.frame(Water_ID = all_data_biof$Water_ID,
                              biof_risk = all_data_biof$mean_biof_risk_scaled,
                              energy = all_data2$energy_scaled)

VI_risk <- data.frame(Water_ID = all_data2$Water_ID,
                      VI_risk = all_data2$mean_risk_scaled)

risk_comparison <- left_join(biof_risk,VI_risk)

plot(risk_comparison$VI_risk, risk_comparison$biof_risk)

risk_comparison$energy_rank <- rank(-risk_comparison$energy, ties.method = "first")

risk_comparison <- risk_comparison %>%
  arrange(energy)
png("figures/corr_biofoul_VI_risk.png",height = 9, width = 9, units = "in",res=300)
ggplot(risk_comparison, aes(VI_risk,biof_risk))+
  geom_point(aes(col = (energy)))+
  geom_abline(slope=1, col="red",linetype="dashed")+
  theme_classic(base_size = 20)+
  scale_color_viridis(option="inferno")
dev.off()

cor(risk_comparison$VI_risk,risk_comparison$biof_risk)
  

#outliers?

outlier_inds <- which(lake_biofoul_df1$mean_biof_risk %in% boxplot.stats(lake_biofoul_df1$mean_biof_risk)$out)

outlier_data <- lake_biofoul_df1[outlier_inds,]

lakes <- read_sf("C:/Users/allis/OneDrive/Post-doc/big_data/floating_solar/Northeast_NHD_Alison")
lakes_selected <- lakes %>%
  filter(Suitabl_FP ==1)

outlier_data_lakes <- left_join(outlier_data,lakes_selected)
outlier_data_lakes$energy_scaled <- scale(outlier_data_lakes$year1_ener)[,1]

outlier_data_lakes <- outlier_data_lakes %>%
  arrange((mean_biof_risk))

png("figures/biofouling_risk_outliers.png", height = 9, width = 11, units = "in",res=300)

ggplot()+
  geom_sf(data = NE_pro)+
  theme_classic(base_size = 15)+
  geom_point(data = outlier_data_lakes, aes(x = water_lon, y = water_lat, col = mean_biof_risk_scaled, size = energy_scaled), alpha = 0.5)+
  scale_color_viridis(option="inferno",limits = c(0,30))+
  ylab("")+
  xlab("")+
  scale_size(guide="none")
#scale_color_gsea(reverse = TRUE, limits = c(-19,19))

dev.off()




# plot_data <- only_selected_lakes %>%
#   arrange((sum_biofoul_risk))%>%
#   filter(sum_biofoul_risk!=0)
# 
# png("figures/sum_biofouling_solar_overlay.png", height = 12, width = 12, units = "in",res=300)
# 
# ggplot()+
#   geom_sf(data = NE_pro)+
#   theme_void()+
#   geom_point(data = plot_data, aes(x = water_lon, y = water_lat, color = sum_biofoul_risk, size = year1_ener))+
#   scale_color_viridis()
# 
# dev.off()
# 
# ggplot(data = plot_data, aes(x = (year1_ener), y = (sum_biofoul_risk)))+
#   geom_point()
# 





#### comparisons ####
library(rstan)

load("data_outputs/lake_biofoul_risk_df.RData")
load("data_outputs/lake_risk_df.RData")

data <- read.csv("data_outputs/final_analysis_data_n291.csv")

lakes <- read_sf("C:/Users/allis/OneDrive/Post-doc/big_data/floating_solar/Northeast_NHD_Alison")

selected_lakes <- lakes %>%
  filter(Suitabl_FP ==1)

all_data_biof <- left_join(selected_lakes,lake_biofoul_df1)
all_data <- left_join(all_data_biof,lake_risk_df)

all_data_ranked <- data.frame(Water_ID = all_data$Water_ID,
                              VI_value = all_data$mean_risk_scaled,
                              VI_rank = rank(all_data$mean_risk_scaled),
                              biof_rank = rank(all_data$mean_biof_risk_scaled),
                              water_quality = all_data$Biodiversi,
                              social_value = all_data$Social_B_1)

all_data_ranked <- all_data_ranked %>%
  arrange(VI_rank)

save(all_data_ranked, file = "data/all_data_model.RData")
load("data/all_data_model.RData")

scen_mod_data <- list(N = length(all_data_ranked$Water_ID),
                 n_WQ = as.integer(2),
                 n_SOC = as.integer(2),
                 VI_value = all_data_ranked$VI_value,
                 WQ_value = all_data_ranked$water_quality+1,
                 SOC_value = all_data_ranked$social_value+1)



scen_mod_fit <- stan(file = "models/scenario_comparison_model.stan",
                   data = scen_mod_data)

save(scen_mod_fit, file = "mod_outputs/scen_mod_fit.RData")
load("mod_outputs/scen_mod_fit.RData")

summary(scen_mod_fit)

library(shinystan)
launch_shinystan(scen_mod_fit)

draws <- as.data.frame(scen_mod_fit)

#lakes have both WQ and SOC value
WQ_SOC <- draws$`WQ_intercept[2]`+ draws$`SOC_intercept[2]`
WQ_noSOC <- draws$`WQ_intercept[2]`+ draws$`SOC_intercept[1]`
SOC_noWQ <- draws$`WQ_intercept[1]`+ draws$`SOC_intercept[2]`

#high probability density interval
#library(rethinking)

comp_plot_data1 <- data.frame(scenario = c("Freshwater & Social","Freshwater Only","Social Only"),
                              mean_VI = c(mean(WQ_SOC),mean(WQ_noSOC),mean(SOC_noWQ)),
                              hpdi_u = c(HPDI(WQ_SOC, prob = 0.9)[2],HPDI(WQ_noSOC, prob = 0.9)[2],HPDI(SOC_noWQ, prob = 0.9)[2]),
                              hpdi_l = c(HPDI(WQ_SOC, prob = 0.9)[1],HPDI(WQ_noSOC, prob = 0.9)[1],HPDI(SOC_noWQ, prob = 0.9)[1]))

png("figures/scenario_comparison_plot.png", height = 10, width = 10, units="in", res = 300)

ggplot(comp_plot_data1, aes(x=scenario, y = mean_VI))+
  geom_errorbar(aes(ymin=hpdi_l,ymax=hpdi_u), width = 0.1, size = 1)+
  geom_point(size = 3, shape = 21, fill = "white")+
  geom_hline(aes(yintercept=0), linetype = "dashed", col = "blue")+
  theme_classic(base_size = 18)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1))+
  scale_y_continuous(limits=c(-0.5,0.5))+
  ylab("Mean Avian Risk")+
  xlab("")
  
dev.off()

#try violin plot

comp_plot_data <- data.frame(WQ_SOC=WQ_SOC,
                             WQ_noSOC=WQ_noSOC,
                             SOC_noWQ=SOC_noWQ)

comp_plot_data1 <- pivot_longer(comp_plot_data, names_to = "water_body_value", cols=1:3, values_to="draws")

ggplot(comp_plot_data1, aes(x=water_body_value, y = draws))+
  geom_violin()


ggplot(data = plot_data, aes(x = (richness), y = (sum_biofoul_risk)))+
  geom_point() #high species richness does not necessarily indicate high biofouling risk




# energy scenarios --------------------------------------------------------
library(sf)
library(tidyverse)
library(prioritizr)
library(tmap)

corr()

#simple approach where lakes are selected based on two criteria:
#first, starting with the lake that has the best energy/cost ratio - plot energy per unit cost
#second, starting with the lake that minimizes overall biodiversity risk - plot energy per unit cost
load("data_outputs/lake_risk_df.RData")
load("data/suitable_lakes.RData")


energy_data <- left_join(lakes1,lake_risk_df)

energy_data1 <- energy_data %>%
  mutate(energy_cost_ratio = year1_ener/fpv_ha)%>%
  arrange(desc(energy_cost_ratio))

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

save(energy_data, file = "data_outputs/energy_scenario_data.RData")
load("data_outputs/energy_scenario_data.RData")

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


#bar chart or sexier equivalent? probably shows data better
#include energy produced in each of these combinations

# scenario_comp <- energy_data %>%
#   select(c("Water_ID","Social_B_1","Biodiversi","avian_scenario","new_precaution"))

plot_data <- energy_data %>%
  st_drop_geometry()%>%
  group_by(new_precaution)%>%
  summarise(total_energy = sum(year1_ener))

plot_data$new_precaution <- factor(plot_data$new_precaution, levels = c("none","non-bird value","birds +","birds_only","all"))

#total energy produced in each scenario
ggplot(plot_data, aes(x=new_precaution, y=total_energy))+
  geom_bar(stat="identity")+
  theme_classic()

#energy production tradeoff for each scenario
plot_data$energy_deficit <- plot_data$total_energy-sum(energy_data$year1_ener)
plot_data$new_precaution <- factor(plot_data$new_precaution, levels = c("all","birds_only","birds +","non-bird value","none"))

png("figures/energy_deficit_scenarios.png", height = 10, width = 10, units="in", res = 300)

ggplot(plot_data, aes(x=new_precaution, y=energy_deficit))+
  geom_bar(fill=c("#420D09","#8D021F","#B80F0A","#CD5C5C","#FA8072"), stat="identity")+
  theme_classic(base_size = 22)

dev.off()

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
