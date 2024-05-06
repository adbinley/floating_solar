#download eBird data
library(tidyverse)
library(ebirdst)

species_selection <- read.csv("data/species_selection.csv")
species_selection <- species_selection[1:177,]
species_selection <- species_selection %>%
  filter(include == 1)

ebd_data <- ebirdst_runs
species_data <- left_join(species_selection,ebd_data, by = "common_name")
#eurasian wigeon, little gull, mountain plover, red phalarope,
# Red-faced Cormorant, Red-necked Phalarope, Sharp-tailed Sandpiper,
# Whooping Crane, Yellow Rail, 	Yellow-billed Loon, Yellow-crowned Night-Heron are missing

sps_code <- species_data$species_code %>%
  na.omit()

#download abundance data

for(a in 1:length(sps_code)){
  
  skip_to_next <- FALSE
  
  sps <- sps_code[a]
  
  tryCatch(ebirdst_download_status(sps,
                                   path = "A:/floating_solar/ebird",
                                   download_abundance = T,
                                   download_occurrence = F,
                                   pattern = "abundance_median_3km",
                                   #dry_run = T,
                                   force=T),
                                   error = function(e){skip_to_next <<- TRUE})
  
  if(skip_to_next) {next}
  
}

