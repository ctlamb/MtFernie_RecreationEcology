
#####################################################################
#
# Title: Extracting and annotating covariate data to iSSF data frame
# Author: Eric Palm
#
#####################################################################

require(terra)
require(sf)
require(data.table)
require(dplyr)
require(sfheaders)
require(tidyverse)
#require(distanceto)

covs.ann <- "data/spatial/covs/annual"
covs.static <- "data/spatial/covs/static"
aoi <- st_read("data/aoi_buffered.shp")

# Importing the ssf data frame without covariates.
df_no_covs <- readRDS(file.path("data", 'derived', "ssf_2_hrs_no_covs.rds"))

# This annotates covariates to each location, using either the starts of steps ("x1_", "y1_"),
# or ends of steps ("x2_", "y2_"). Then creates a stratum variable, cos_ta (cosine of turn angle),
# and a year variable for matching annual covariate rasters to location years.
# I used 30-m covariates for the finest spatial resolution
# Each raster needs to be named (internally) whatever you want to call that variable in the model, and you can add "start" or "end" suffix.
df_covs <- 
  df_no_covs %>% 
  dplyr::bind_cols(., terra::extract(list.files(covs.static, full.names = T, pattern = ".tif") %>%
                                       #stringr::str_subset(., "d_|01") %>% 
                                       terra::rast(), cbind(.$x2_, .$y2_)) %>%
                     dplyr::rename_with(~stringr::str_c(., "_end"), dplyr::everything())) %>% 
  dplyr::bind_cols(., terra::extract(list.files(covs.static, full.names = T, pattern = ".tif") %>%
                                       #stringr::str_subset(., "elevation") %>%
                                       terra::rast(), cbind(.$x1_, .$y1_)) %>%
                     dplyr::rename_with(~stringr::str_c(., "_start"), dplyr::everything())) %>%
  dplyr::mutate(stratum = stringr::str_c(id, t2_, sep = "_"),
                cos_ta = cos(ta_),
                year = lubridate::year(t2_),
                delta_elev = abs(elevation_end - elevation_start)
                ) 

list.files(covs.static, full.names = T, pattern = ".tif")[3]%>%rast()

##remove any steps way outside AOI, where 3km buffer of elev still doesn't catch starting steps.
df_covs <- df_covs%>%
  drop_na(delta_elev)


# Selecting annual rasters whose years match the years with locations
# List all .tif files
r_files <- list.files(covs.ann, full.names = TRUE, pattern = "\\.tif$", ignore.case = TRUE)

# Extract covariate name and year from each file
r_tbl <- tibble(
  file = r_files,
  fname = basename(r_files),
  year = str_extract(fname, "\\d{4}") %>% as.integer(),
  var = str_remove(fname, "_?\\d{4}\\.tif$")  # Strip trailing _YYYY.tif or YYYY.tif
) %>%
  filter(!is.na(year)) %>%
  mutate(rast = map(file, terra::rast))

df_covs_ann <- df_covs


for (v in unique(r_tbl$var)) {
  r_subset <- filter(r_tbl, var == v)
  
  # Max year for this variable
  var_year_max <- max(r_subset$year, na.rm = TRUE)
  
  # Use variable-specific year match
  df_covs_ann <- df_covs_ann %>%
    mutate(year_match = pmin(lubridate::year(t2_), var_year_max)) %>%
    group_by(year_match) %>%
    nest() %>%
    left_join(r_subset, by = c("year_match" = "year")) %>%
   mutate(
      vals = map2(data, rast, ~ terra::extract(.y, cbind(.x$x2_, .x$y2_))[, 1]),
      data = map2(data, vals, ~ mutate(.x, !!sym(v) := .y))
    ) %>%
    select(data) %>%
    unnest(data) %>%
    ungroup() %>%
    select(-year_match)  # optional: clean up
}

# #### Distance to disturbances ####
# # I used the 'distanceto' R package for calculating distances,
# # as it is way faster than sf's st_join and st_distance options.
# # First I create an sf object of the ssf dataframe to use for distance calculations.
# # I only use the "end" of steps for these, but you could add the beginning as well if desired.
#   
# df_sf <- 
#   df_covs_static %>%
#   sf::st_as_sf(., coords = c("x2_", "y2_"), crs = 3005, remove = F)
#   
# 
# # Load in all the disturbance shapefiles and make sure they all have the same column name for year.  
# burn <- sf::read_sf(file.path(shp, "burn.shp"))
# seismic <- sf::read_sf(file.path(shp, "seismic.shp")) 
# cutblock <- sf::read_sf(file.path(shp, "cutblock.shp")) 
#   
# well <- sf::read_sf(file.path(shp, "well.shp")) %>% 
#   dplyr::select(year_sf = YEAR_ES) %>% 
#   st_cast(., "POLYGON")
#   
# mine <- sf::read_sf(file.path(shp, "mine.shp")) %>% 
#   dplyr::select(year_sf = YEAR_ES) %>% 
#   st_cast(., "POLYGON") 
#   
# pipeline <- sf::read_sf(file.path(shp, "pipeline_no_rupert.shp")) %>% 
#   dplyr::select(year_sf = YEAR_ES) %>% 
#   st_cast(., "POLYGON") 
#   
# ag <- sf::read_sf(file.path(shp, "ag.shp")) %>% 
#   dplyr::select(year_sf = YEAR_ES) %>% 
#   st_cast(., "POLYGON")
# 
# wind <- sf::read_sf(file.path(shp, "wind.shp")) %>% 
#   dplyr::select(year_sf = YEAR) %>% 
#   st_cast(., "POLYGON") 
#   
# # Here I extract distances to nearest disturbances, making sure the location year
# # is after the disturbance year, plus I make sure the location year is no more than 50 years
# # after burns. You could change this threshold though.
# # This entire step took about 6 minutes on the server I'm using.
# tictoc::tic()
# df_dist <-
#   df_sf %>% 
#   dplyr::nest_by(year) %>% 
#   dplyr::mutate(d_burn_end = list(distanceto::distance_to(data, burn[burn$year_sf < year & year - burn$year_sf < 50, ])),
#                 d_cutblock_end = list(distanceto::distance_to(data, cutblock[cutblock$year_sf < year, ])),
#                 d_seismic_end = list(distanceto::distance_to(data, seismic[seismic$year_sf < year, ])),
#                 d_well_end = list(distanceto::distance_to(data, well[well$year_sf < year, ])),
#                 d_mine_end = list(distanceto::distance_to(data, mine[mine$year_sf < year, ])),
#                 d_pipeline_end = list(distanceto::distance_to(data, pipeline[pipeline$year_sf < year, ])),
#                 d_ag_end = list(distanceto::distance_to(data, ag[ag$year_sf < year, ]))) %>% 
#   tidyr::unnest(cols = c(data:d_ag_end)) %>% 
#   dplyr::ungroup() %>% 
#   sf::st_sf()
# tictoc::toc()
#   
# # I had to do wind differently because there weren't any wind turbines until 2007
# df_dist_wind <-
#   df_dist %>% 
#   dplyr::filter(year <= 2007) %>% 
#   dplyr::mutate(d_wind_end = 130000) %>% 
#   dplyr::bind_rows(., df_dist %>% 
#                      dplyr::filter(year > 2007) %>% 
#                      dplyr::nest_by(year) %>% 
#                      dplyr::mutate(d_wind_end = list(distanceto::distance_to(data, wind[wind$year_sf < year, ]))) %>% 
#                      tidyr::unnest(cols = c(data, d_wind_end)) %>% 
#                      dplyr::ungroup() %>% 
#                      sf::st_sf()) %>% 
#   sf::st_drop_geometry()
# 
# # Vector of unique years in ssf dataframe
# years <- sort(unique(df_sf$year))
# 
# # Add in canopy cover (separate raster for each year) and start and end of steps
# df_covs_all_dist <-
#   df_dist_wind %>%
#   dplyr::nest_by(year) %>% 
#   dplyr::ungroup() %>% 
#   dplyr::mutate(rast_cc = list.files(file.path(covs90m, "cc"), full.names = T) %>%
#                   stringr::str_subset(., stringr::str_c(years, collapse = "|")) %>%
#                   purrr::map(., rast)) %>% 
#     dplyr::rowwise() %>% 
#     dplyr::mutate(cc_end = list(as.double(terra::extract(rast_cc, cbind(data$x2_, data$y2_))[,1])),
#                   cc_start = list(as.double(terra::extract(rast_cc, cbind(data$x1_, data$y1_))[,1]))) %>% 
#     dplyr::select(-c(rast_cc)) %>% 
#     tidyr::unnest(cols = c(data, cc_end:cc_start)) %>% 
#   dplyr::ungroup()



# Annotating semi-permeable barriers for towns, mines, glacier/rock, lakes/reservoirs, highways.
# Need to create line segments for successive locations.
# I use 'data.table' and 'sfheaders' packages because they are super fast, but it's
# probably not necessary...I was just looking for ways to speed up the simulations
# and ended up using that code here as well.

# Preparing data frame for data.table
dt <- df_covs_ann %>%
  dplyr::mutate(index = dplyr::row_number()) %>% 
  dplyr::select(x1_:y2_, index) %>% 
  data.table::as.data.table()

stopifnot(length(unique(dt$index)) == nrow(df_covs_ann))  # ✅ Should now be TRUE

## To use `sfheaders` the data needs to be in long form
dt1 <- dt[, .(index, lon = x1_, lat = y1_)]
dt2 <- dt[, .(index, lon = x2_, lat = y2_)]

## Add on a 'sequence' variable so we know which one comes first
dt1[, seq := 1L ]
dt2[, seq := 2L ]

# Binding rows with start locations (x1_, y1_)
# to rows with end locations (x2_, y2_)
dt <- data.table::rbindlist(list(dt1, dt2), use.names = TRUE)

# Order the data.table
data.table::setorder(dt, index, seq)

# Create sf line for each set of successive locations
line_sf <-
  sfheaders::sf_linestring(
    obj = dt
    , x = "lon"
    , y = "lat"
    , linestring_id = "index"
  ) %>% 
  sf::st_set_crs(3005) 


# Use sf::st_intersects to determine whether each line intersects landscape features (0 if no, 1 if yes).
# Create "log distance to" variables and add a "log_" prefix for those.
trails <- st_read("data/spatial/covs/linear/trails.shp")
trails_hu <- st_read("data/spatial/covs/linear/trails_hu.shp")
trails_multi <- st_read("data/spatial/covs/linear/trails_multi.shp")

rd <- st_read("data/spatial/covs/linear/rd_no_hwy.shp")
hwy <- st_read("data/spatial/covs/linear/hwy.shp")

builtup <- st_read("data/spatial/covs/linear/builtup.shp")

river <- st_read("data/spatial/covs/linear/river.shp")

df_covs_barriers <- df_covs_ann %>%
  dplyr::mutate(
    barrier_hwy = as.integer(lengths(sf::st_intersects(line_sf, hwy)) > 0),
    barrier_rd = as.integer(lengths(sf::st_intersects(line_sf, rd)) > 0),
    barrier_trail = as.integer(lengths(sf::st_intersects(line_sf, trails)) > 0),
    barrier_trail_hu = as.integer(lengths(sf::st_intersects(line_sf, trails_hu)) > 0),
    barrier_trail_multi = as.integer(lengths(sf::st_intersects(line_sf, trails_multi)) > 0),
    barrier_builtup = as.integer(lengths(sf::st_intersects(line_sf, builtup)) > 0),
    barrier_river = as.integer(lengths(sf::st_intersects(line_sf, river)) > 0),
    rd_cross_count = lengths(sf::st_intersects(line_sf, rd)),
    trail_cross_count = lengths(sf::st_intersects(line_sf, trails)),
    
    # NEW: mean riders2024 for intersecting trails
    trail_use_mean = map_dbl(
      sf::st_intersects(line_sf, trails),
      function(idxs) {
        if (length(idxs) > 0) {
          mean(trails$Riders2024[idxs], na.rm = TRUE)
        } else {
          0
        }
      }
    )
  ) %>%
  dplyr::mutate(across(starts_with("dist_to"), ~log(. + 1), .names = "log_{col}")) %>%
  tibble() %>%
  select(-geometry)

##log some others too
df_covs_barriers <- df_covs_barriers%>%
  dplyr::mutate(dplyr::across(tidyselect::starts_with(c("tri", "bldg.dens")), ~log(. + 1), .names = "log_{col}"))


# Save the annotated data frame 
saveRDS(df_covs_barriers, file.path('data', 'derived',"ssf_2_hrs_unscaled.rds"))


df_covs_barriers%>%
  filter(case_==1)%>%
  summarise_at(vars(barrier_hwy:barrier_river),  sum)

df_covs_barriers%>%
  group_by(case_==1)%>%
  summarise_at(vars(barrier_hwy:barrier_river),  sum)%>%
  pivot_longer(barrier_hwy:barrier_river)%>%
  pivot_wider(names_from = `case_ == 1`, values_from="value")%>%
  mutate(`TRUE`/(`FALSE`)*10)

            