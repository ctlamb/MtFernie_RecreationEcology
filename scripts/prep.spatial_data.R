library(here)
library(sf)
library(terra)
library(sp)
library(spatialEco)
library(tidyverse)


# 1. List all .tif files recursively
tif_files <- list.files(path = here::here("data/spatial/raw"), pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)

# 2. Identify the land cover file
lc_file <- "/Users/claytonlamb/Dropbox/Documents/University/Work/RecreationMonitoring/MtFernie/data/spatial/raw/CA_forest_VLCE2_2022/CA_forest_VLCE2_2022.tif"

# 3. Separate the files
other_files <- setdiff(tif_files, lc_file)

# 4. Load rasters
lc_raster <- rast(lc_file)
other_rasters <- lapply(other_files, rast)

# 5. Choose reference raster (first non-landcover one)
ref <- other_rasters[[1]]

# 6. Align continuous rasters to reference (bilinear OK)
aligned_continuous <- lapply(other_rasters, function(r) {
  r_proj <- terra::project(r, ref)
})

# 7. Align land cover raster (use nearest-neighbor)
lc_raster_proj <- terra::project(lc_raster, ref, method = "near")
names(lc_raster_proj) <- "lc"

##mask to project extent

# 1. Get Fernie's coordinates in lat/long
fernie_ll <- st_sf(
  name = "Fernie",
  geometry = st_sfc(st_point(c(-115.0631, 49.5044)), crs = 4326)  # EPSG:4326 WGS84
)

# 2. Transform to EPSG:3348 (Statistics Canada Lambert)
fernie_proj <- st_transform(fernie_ll, crs = 3348)

# 3. Buffer by 30 km (30,000 meters)
fernie_buffer <- st_buffer(fernie_proj, dist = 30000)

# 4. Convert to SpatVector (terra) and make extent
fernie_vect <- vect(fernie_buffer)
fernie_extent <- ext(fernie_vect)

# 5. Use to crop rasters
aligned_continuous_clipped <- terra::crop(rast(aligned_continuous), fernie_extent)
lc_raster_proj_clipped <- terra::crop(lc_raster_proj, fernie_extent)


# 8. Combine everything into one stack
r_stack <- c(aligned_continuous_clipped, lc_raster_proj_clipped)


##export
# Loop through each layer in the stack
for (i in 1:nlyr(r_stack)) {
  layer_name <- names(r_stack)[i]
  output_file <- here::here("data/spatial/stacked", paste0(layer_name, ".tif"))
  
  writeRaster(r_stack[[i]], filename = output_file, overwrite = TRUE)
}


##spatial scales to summarize across
scale <- c(250, 500)  #meters

#loop
for(i in 1:length(names(r_stack))){
  rast.i <- r_stack[[i]]
  layer_name <- names(rast.i)
  
  for(x in 1: length(scale)){
    
    circle_window <- focalMat(rast.i, type = "circle", d = scale[x]) ##moving window
    
    if(names(rast.i)!="lc"){
    rast.x <- rast.i%>%
      terra::focal(w = circle_window, fun = mean, na.rm = TRUE)}
    else{
      rast.x <- rast.i%>%
        terra::focal(w = circle_window, fun = modal, na.rm = TRUE)
      }
    names(rast.x) <-paste0(layer_name,"_",scale[x])
   
     #export
    writeRaster(rast.x, here::here("data/spatial/stacked_focal",paste0(layer_name,"_",scale[x],"m.tif")), overwrite=TRUE)
    rm(rast.x)
  }
  rm(rast.i)
}


##do a far one for building density as an alternative
circle_window <- focalMat(r_stack[["bldg.count"]], type = "circle", d = 3000) ##moving window

rast.x <- r_stack[["bldg.count"]]%>%
    terra::focal(w = circle_window, fun = mean, na.rm = TRUE)

plot(rast.x)
names(rast.x) <- "bldg.count_3000"

writeRaster(rast.x, here::here("data/spatial/stacked_focal",paste0("bldg.count_3000","m.tif")), overwrite=TRUE)

