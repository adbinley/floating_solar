#tests and troubleshooting
library(sf)
library(tidyverse)

#trouble: ~300 lakes have 0 species richness and 0 biodiversity value... error?
load("D:/floating_solar/data_outputs/all_importance_data.RData")

no_birds <- all_data2 %>%
  filter(richness == 0)
unique(no_birds$TNC_Bio_Va)
range(no_birds$year1_ener)

# load gis data for making maps
state_lines <- st_read("D:/maps/NA_politicalboundaries_shapefile/PoliticalBoundaries_Shapefile/NA_PoliticalDivisions/data/bound_p/boundaries_p_2021_v3.shp")

NE <- state_lines %>%
  filter(NAME_En %in% c("Maine","New Hampshire","Vermont","Massachusetts","Rhode Island","Connecticut","New York","Pennsylvania","New Jersey","Delaware","Maryland","West Virginia","Virginia","District of Columbia"))%>%
  st_geometry()

NE_pro <- st_transform(NE, crs = st_crs(4326))

ggplot()+
  geom_sf(data = NE_pro)+
  theme_void()


ggplot()+
  geom_point(data = no_birds, aes(x = water_lon, y = water_lat, col = "red"))+
  geom_sf(data = NE_pro)+
  theme_void()
#trouble is shot - weird waterbodies are off the coast and not captured well be eBIrd, plus not relevant for floating solar
#removed from analysis - did you actually remove these though?


no_poo <- all_data2 %>%
  filter(sum_biofoul_risk == 0) #same deal


#####

sp <- "cangoo"
sp <- na_species[5]

bird_data <- rast(paste0("D:/big_data/eBird_FAC/2022/",sp,"/weekly/",sp,"_abundance_median_3km_2022.tif"))
test <- rast(paste0("D:/floating_solar/generated/",sp,"_max_values.tif"))
test
plot(test)
test1 <- test*10000

for(t in 1:52){
  
  max <- global(week_list[[t]], fun="max")
  
  print(max)
  
}

freq(week_list[[46]], value=0)

subset_test <- test[test>0 & test<1]

freq(test, value=0)

bobo <- lake_biodiversity_df %>%
  filter(species_code == "boboli")
purmar <- lake_biodiversity_df %>%
  filter(species_code == "purmar")
chiswi <- lake_biodiversity_df %>%
  filter(species_code == "chiswi")

#get rid of annoying 1s
test <- bird_data
values(test)[values(test)==1]<-0
plot(test)
freq(bird_data, value=0)

mat <- as.matrix(bird_data1)

still.na <- lake_biodiversity_df %>%
  filter(is.na(max))
unique(still.na$species_code)
load("D:/floating_solar/data_outputs/lake_ave_biodiversity.RData")


test <- rescale(data$exposure, to=c(1,5))

tas_pu <- prioritizr::get_tas_pu()


lakes_clean <- lakes2 %>%
  select(c("Water_ID","fpv_ha","year1_ener","mean_risk","Social_B_1"))

#lake overlap problem

lakes <- read_sf("D:/floating_solar/Northeast_NHD_Alison")
buf <- 1000
lake_buffer <- st_buffer(lakes,buf)
lakes_sel <- lake_buffer %>%
  filter(Suitabl_FP==1)
#lakes_vec <- vect(lake_buffer)

pdf("figures/1km_buffer_overlap.pdf")
plot(st_geometry(lakes_sel))
dev.off()
