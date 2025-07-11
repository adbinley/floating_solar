library(prioritizr)
library(sf)
library(tidyverse)
library(ebirdst)
library(terra)
library(stringr)
library(viridis)
library(scales)



# 1. Create current zone layer --------------------------------------------
#create current zone layer (raster stack)

species_data <- list.files(path = "D:/floating_solar/generated/")

current_zone_layer <- rast(paste0("D:/floating_solar/generated/",species_data))




# 2. Create solar zone ----------------------------------------------------
#create solar zone layer using VI for each species


VI_data <- read.csv("data_outputs/final_analysis_data.csv")

#assumed that species with the highest VI decrease by 90% of their max values
#and species with the lowest VI decrease by 10% of their max values

VI_data$reduction_factor <- 1-rescale(VI_data$VI, to=c(0.1,0.9))

#sort so that values align with raster stack
VI_data_arr <- VI_data %>%
  arrange(species_code)

#only regions where solar panels could possibly be implemented will be affected
#ie, reductions will only occur on lakes where there is solar
#start by getting lake buffer
lakes <- read_sf("D:/floating_solar/Northeast_NHD_Alison")
buf <- 5000
lake_buffer <- st_buffer(lakes,buf)

#filter out unsuitable lakes - these values won't change
lake_buffer_s <- lake_buffer %>%
  filter(Suitabl_FP ==1)

lakes_vec <- vect(lake_buffer_s)#create SpatVector
#project to ebird data crs
sp <- "amerob" #species doesn't matter, using for projection
bird_data <- rast(paste0("D:/floating_solar/generated/",sp,"_max_values.tif"))
lakes_vec_pro <- project(lakes_vec, crs(bird_data))

#mask the raster stack to the lake buffer
masked_raster <- mask(current_zone_layer,lakes_vec_pro)

#calculate new values based on the scaled VI
solar_reduced_values <- masked_raster*VI_data_arr$reduction_factor

#solar zone layer
#solar_zone_layer <- ifelse(is.na(masked_raster),current_zone_layer,solar_reduced_values)
#solar_zone_layer <- rast(solar_zone_layer)

#merge with raster values with no floating solar
#reductions will occur in regions where solar is possible, but values outside of these regions are retained
solar_zone_layer <- terra::merge(solar_reduced_values,current_zone_layer, first=T)


# 3. Trim both zones to study region --------------------------------------

#map data
state_lines <- st_read("D:/maps/NA_politicalboundaries_shapefile/PoliticalBoundaries_Shapefile/NA_PoliticalDivisions/data/bound_p/boundaries_p_2021_v3.shp")

NE <- state_lines %>%
  filter(NAME_En %in% c("Maine","New Hampshire","Vermont","Massachusetts","Rhode Island","Connecticut","New York","Pennsylvania","New Jersey","Delaware","Maryland","West Virginia","Virginia","District of Columbia"))%>%
  st_geometry()

NE_vec <- vect(NE)#create SpatVector

NE_pro <- project(NE_vec, crs(bird_data))

current_zone_layer_trim <- mask(current_zone_layer,NE_pro)

solar_zone_layer_trim <- mask(solar_zone_layer,NE_pro)

#crop to study extent
current_zone_cr <- crop(current_zone_layer_trim,NE_pro)
solar_zone_cr <- crop(solar_zone_layer_trim, NE_pro)

#make values easier to deal with
current_zone_tr <- current_zone_cr*10000
solar_zone_tr <- solar_zone_cr*10000

#save raster stacks
writeRaster(current_zone_tr, "D:/floating_solar/data_outputs/current_zone_stack.tif", overwrite=T)
writeRaster(solar_zone_tr, "D:/floating_solar/data_outputs/solar_zone_stack.tif", overwrite = T)


# 4. Get absolute target values ----------------------------------------------

#target is the total "abundance" of each species across the extent
#idea is to try to preserve as close to this maximum value as possible

#load raster stacks
current_zone <- rast("A:/floating_solar/data_outputs/current_zone_stack.tif")
solar_zone <- rast("D:/floating_solar/data_outputs/solar_zone_stack.tif")

#get sum across study extent for each species

abs_target <- global(current_zone, fun = "sum", na.rm=T)

VI_data_arr$abs_target <- abs_target$sum

#filter out species out of range
 VI_data_arr1 <- VI_data_arr %>%
   filter(abs_target > 2.47)
 
write.csv(VI_data_arr1, file = "data_outputs/final_analysis_data_n291.csv")

#get indices for species removed
keep_ind <- which(VI_data_arr$abs_target > 2.47)

#targets for lakes only (unsure which is more useful atm)
current_zone_lakes_cr <- crop(masked_raster,NE_pro)

current_zone_lakes_cr1 <- subset(current_zone_lakes_cr, keep_ind)

abs_target_lakes <- global(current_zone_lakes_cr1, fun = "sum", na.rm = T)

VI_data_arr1$abs_target_lakes <- abs_target_lakes$sum

save(VI_data_arr1, file = "data_outputs/final_analysis_data_targets.RData")

current_zone1 <- subset(current_zone, keep_ind)
solar_zone1 <- subset(solar_zone, keep_ind)

writeRaster(current_zone1, "D:/floating_solar/data_outputs/current_zone_final.tif", overwrite=T)
writeRaster(solar_zone1, "D:/floating_solar/data_outputs/solar_zone_final.tif", overwrite = T)


# 5. Deal with overlapping MUs --------------------------------------------

#new approach

#create rasters of lakes and determine degree of overlap

#lakes <- read_sf("D:/floating_solar/Northeast_NHD_Alison")
lakes <- read_sf("C:/Users/allis/OneDrive/Post-doc/big_data/floating_solar/Northeast_NHD_Alison")
lakes1 <- lakes %>%
  filter(Suitabl_FP ==1)
rm(lakes)

#using current zone raster to get raster dimensions and crs
current_zone <- rast("A:/floating_solar/data_outputs/current_zone_final.tif")

#create raster of 1s
pu_raster <- rast(
  nrows = 409, ncols = 789,
  xmin = -7470586, xmax = -5132931, ymin = 4064160, ymax = 5275949, vals = 1
)

#project lakes, rasterize 
lake_raster <- lakes1 %>%
  vect()%>%
  project(crs(current_zone))%>%
  terra::rasterize(pu_raster, field = 1, touches=T)

lake_proj <- lakes1 %>%
  vect()%>%
  project(crs(current_zone))

pdf("figures/lake_raster.pdf")
plot(lake_raster)
plot(lake_proj, add=T, lwd=0.02)
dev.off()

# create a separate layer for each sea mount
lake_stack <-
  lake_raster %>%
  terra::patches() %>%
  as("Raster") %>%
  {lapply(sort(unique(values(.))), function(i) {
    raster::Which(. == i)
  })} %>%
  raster::stack() %>%
  setNames(paste0("mu_", seq_len(raster::nlayers(.))))

#ugh

#old approach - too many dimensions, not enough computing power

#several of the buffered lakes overlap
#to avoid double counting in overlapping regions, we treat these as a potential management unit that can be selected

#lakes <- read_sf("D:/floating_solar/Northeast_NHD_Alison")
# lakes <- read_sf("C:/Users/allis/OneDrive/Post-doc/big_data/floating_solar/Northeast_NHD_Alison")
# lakes1 <- lakes %>%
#   filter(Suitabl_FP ==1)
# buf <- 1000
# lake_buffer <- st_buffer(lakes1,buf)
# 
# #get indices for which lakes overlap
# lake_overlap <- st_overlaps(lake_buffer)
# 
# #test to see if working properly
# # i <- 8
# # plot(st_geometry((lake_buffer[lake_overlap[[i]],])))
# # plot(st_geometry(lake_buffer[i,]), col="red",add=T)
# 
# #calculate the different combinations that overlap with each lake
# #note that these DO NOT necessarily overlap with one another
# #only with the original indexed lake
# powerset <- function(set) {
#   n <- length(set)
#   result <- list()
#   # Iterate over all possible subset sizes
#   for (k in 1:n) {
#     # Generate combinations of size k
#     result <- c(result, combn(set, k, simplify = FALSE))
#   }
#   return(result)
# }
# 
# #loop through all lakes and find all combinations
# 
# combos_list <- list()
# 
# for(i in 1:length(lake_overlap)){
# #for(i in 18:21){ #- testing to see how 1s and zeros are working
# 
#   #if the me does not overlap with any lakes, skip to next
#     if(length(lake_overlap[[i]])==0) {
#       next
#     }
# 
#     #combn() does something weird if there is only one mu overlapping, corrected here
#            ifelse(length(lake_overlap[[i]])==1, powerset_results <- lake_overlap[[i]],
#                   powerset_results <- powerset(lake_overlap[[i]]))
# 
#     #add original indexed lake into the vector
#     powerset_results <- lapply(powerset_results, function(x) c(x,i))
# 
#     combos_list <- c(combos_list,powerset_results)
# 
# 
# }



#if all are selected, non of the component buffers can be also individually selected
#set linear penalties to prevent this


# prioritization ----------------------------------------------------------

#load raster stacks
current_zone <- rast("D:/floating_solar/data_outputs/current_zone_final.tif")
solar_zone <- rast("D:/floating_solar/data_outputs/solar_zone_final.tif")

VI_data <- load("data_outputs/final_analysis_data_targets.RData")

lakes <- read_sf("D:/floating_solar/Northeast_NHD_Alison")
lakes1 <- lakes %>%
  filter(Suitabl_FP ==1)


p_1 <- prioritizr::problem(x=lakes1, features="year1_ener", cost_column = "fpv_ha")%>%
