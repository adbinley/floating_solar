#tests and troubleshooting
library(sf)
library(tidyverse)

#trouble: ~300 lakes have 0 species richness and 0 biodiversity value... error?
load("D:/floating_solar/data_outputs/all_importance_data.RData")

no_birds <- all_data1 %>%
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
#removed from analysis