# floating_solar

##Scripts and descriptions

### 1-download_data.R

Downloading all relative abundance data from eBird for all species in the region. Data used in this project are the 2022 version:

Fink, D., T. Auer, A. Johnston, M. Strimas-Mackey, S. Ligocki, O. Robinson, 
  W. Hochachka, L. Jaromczyk, C. Crowley, K. Dunham, A. Stillman, I. Davies, 
  A. Rodewald, V. Ruiz-Gutierrez, C. Wood. 2023.
  eBird Status and Trends, Data Version: 2022; Released: 2023. Cornell Lab of
  Ornithology, Ithaca, New York. https://doi.org/10.2173/ebirdst.2022

  Here we use the median relative abundance at a 3km by 3 km resolution.

  ### 2-summarise_birds_lakes.R

**Importance metric (sum):** Loading the weekly relative abundance data into a raster stack (52 rasters, one for each week). Relative abundance values in each cell are divided by the sum of relative abundance values across the range, yielding the proportion of the total relative abundance for the species that each cell represents. This tells us the relative importance of each cell relative to other cells. We then took the maximum of these values across all weeks.

For each lake, we created a 5km buffer around the shoreline. We extracted the importance values (calculated as per above) that fell within this buffer for each lake, and summed the value of all cells within the buffer for each lake. Thus, larger lakes/buffers will have more cells and likely greater summed values, reflective of the fact that they provide more habitat for birds.

### 2-alternate_mean_values.R

**Importance metric (mean):** Loading the maximum relative importance values for each species and cell (as calculated in *2-summarise_birds_lakes.R*), but this time calculating the **average** cell value for each lake, using the cells extracted within a 5km buffer around the shoreline. This metric does not weight larger lakes as more important than smaller lakes, but rather looks at the relative value of a the average cell near a lake compared to the average value of a cell near another lake.

### 3-vulnerability_scoring
