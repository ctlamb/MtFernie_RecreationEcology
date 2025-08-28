library(sf)
library(tidyverse)
library(terra)
library(mapview)
library(here)

##load AOI
aoi <- st_read("data/aoi_buffered.shp")

##load linear features
trails <- st_read("/Users/claytonlamb/Dropbox/Documents/University/Work/RecreationMonitoring/MtFernie/data/spatial/trails/FTA_TrailNetwork_UPDATE/FTA_TrailNetwork.shp")%>%
  st_make_valid()%>%
  st_transform(3005)

hwy <- st_read("/Users/claytonlamb/Dropbox/Documents/University/Geographic_Data/BCAB_Layers/static/roads/processed/hwy.shp")%>%
  st_intersection(aoi)

rd <- st_read("/Users/claytonlamb/Dropbox/Documents/University/Geographic_Data/BCAB_Layers/static/roads/processed/all_roads.shp")%>%
  st_intersection(aoi)

river <- st_read("/Users/claytonlamb/Dropbox/Documents/University/Geographic_Data/Water/rivers.shp")%>%
  st_transform(3005)%>%
  st_intersection(aoi)

##remove hwy from rd
# Buffer highway by a small distance (e.g., 5 meters)
hwy_buf <- hwy %>%
  st_transform(3005) %>%
  st_buffer(5) %>%
  st_union() %>%         # union speeds up difference ops
  st_make_valid()

# Make sure roads are valid and in same CRS
rd_clean <- rd %>%
  st_transform(3005) %>%
  st_make_valid()

# Erase buffered highway from roads (fast version)
rd_no_hwy <- st_difference(rd_clean, hwy_buf)


##remove trails from rd
trails_buf <- trails %>%
  st_transform(3005) %>%
  st_buffer(10) %>%
  st_union() %>%         # union speeds up difference ops
  st_make_valid()

rd_no_hwy_no_trail <- st_difference(rd_no_hwy, trails_buf)

##remove town areas from rd
bldg.dens<- rast("data/spatial/covs/static/bldg.dens.1km.tif")
bldg.dens.binary <-bldg.dens>20
town_poly <- as.polygons(bldg.dens.binary, dissolve = TRUE) %>%
  st_as_sf() %>%
  filter(bldg.dens.1km > 0) # keep only TRUE areas

rd_no_hwy_no_trail_notown <- st_difference(rd_no_hwy_no_trail, town_poly)

##remove road fragments
rd_cleaned <- rd_no_hwy_no_trail_notown %>%
  mutate(length_m = st_length(geometry)) %>%
  filter(as.numeric(length_m) >= 250)

mapview(rd_cleaned) + mapview(rd_no_hwy_no_trail_notown, color="black")

##create another town area to use as builtup area
town_poly <- as.polygons(bldg.dens>60, dissolve = TRUE) %>%
  st_as_sf() %>%
  filter(bldg.dens.1km > 0) # keep only TRUE areas


##select high use trails
hist(trails$Riders2024)
trails_highuse <- trails%>%filter(Riders2024>500)
mapview(trails_highuse)

##remove hike only trails and some other low use ones around the edge
trails$d_usage%>%unique()
trails$d_season%>%unique()
trails_multi <- trails%>%filter(d_usage%in%c("Multi-use", NA,"Bike Only"),
                                #d_season!="Multi-Use Groomed",
                                !is.na(d_season),
                                !name%in%c("TCT/EVT/CDT (Porky Blue to Road crossing)", 
                                           "TCT/EVT/CDT (road crossing to Hosmer Ruins)",
                                           "Silk Trail - Scandia Loop",
                                           "Mustard Tiger",
                                           "Big Money",
                                           "Coal Discovery Trail",
                                           "Porky Pine Rim",
                                           "Porky Blue",
                                           "Silk Trail",
                                           "Manchiria"))
mapview(trails_multi, zcol="d_season")

trails_not_multi <- trails%>%filter(!d_usage%in%c("Multi-use", NA,"Bike Only")|
                                #d_season!="Multi-Use Groomed",
                                is.na(d_season)|
                                name%in%c("TCT/EVT/CDT (Porky Blue to Road crossing)", 
                                           "TCT/EVT/CDT (road crossing to Hosmer Ruins)",
                                           "Silk Trail - Scandia Loop",
                                           "Mustard Tiger",
                                           "Big Money",
                                           "Coal Discovery Trail",
                                           "Porky Pine Rim",
                                           "Porky Blue",
                                           "Silk Trail",
                                           "Manchiria"))

mapview(trails_not_multi, color="red")+mapview(trails_multi, color="black")

##export LFs
st_write(trails, "data/spatial/covs/linear/trails.shp", delete_dsn = TRUE)
st_write(trails_highuse, "data/spatial/covs/linear/trails_hu.shp", delete_dsn = TRUE)
st_write(trails_multi, "data/spatial/covs/linear/trails_multi.shp", delete_dsn = TRUE)
st_write(hwy, "data/spatial/covs/linear/hwy.shp", delete_dsn = TRUE)
st_write(rd_cleaned, "data/spatial/covs/linear/rd_no_hwy.shp", delete_dsn = TRUE)
st_write(town_poly , "data/spatial/covs/linear/builtup.shp", delete_dsn = TRUE)
st_write(river , "data/spatial/covs/linear/river.shp", delete_dsn = TRUE)


####################
###Distance to
###################
# Create blank raster using AOI extent
# Load raster template
r_template <- rast("data/rast_template.tif")

# Rasterize
trails_vect <- vect(trails)  # convert to SpatVector
trails_hu_vect <- vect(trails_highuse)  # convert to SpatVector
trails_multi_vect <- vect(trails_multi)  # convert to SpatVector
hwy_vect <- vect(hwy)  # convert to SpatVector
rd_no_hwy_vect <- vect(rd_cleaned)  # convert to SpatVector


##Dist to
r_dist_trails <- distance(r_template, trails_vect)
names(r_dist_trails) <- "dist_to_trail"

r_dist_trails_hu <- distance(r_template, trails_hu_vect)
names(r_dist_trails_hu) <- "dist_to_trail_hu"

r_dist_trails_multi <- distance(r_template, trails_multi_vect)
names(r_dist_trails_multi) <- "dist_to_trail_multi"

r_dist_hwy <- distance(r_template, hwy_vect)
names(r_dist_hwy) <- "dist_to_hwy"

r_dist_rd_no_hwy <- distance(r_template, rd_no_hwy_vect)
names(r_dist_rd_no_hwy) <- "dist_to_rd_no_hwy"

# Save output rasters
writeRaster(r_dist_trails, here("data/spatial/covs/static/dist_to_trail.tif"), overwrite = TRUE)
writeRaster(r_dist_trails_hu, here("data/spatial/covs/static/dist_to_trail_hu.tif"), overwrite = TRUE)
writeRaster(r_dist_trails_multi, here("data/spatial/covs/static/dist_to_trail_multi.tif"), overwrite = TRUE)
writeRaster(r_dist_hwy, here("data/spatial/covs/static/dist_to_hwy.tif"), overwrite = TRUE)
writeRaster(r_dist_rd_no_hwy, here("data/spatial/covs/static/dist_to_rdnohwy.tif"), overwrite = TRUE)