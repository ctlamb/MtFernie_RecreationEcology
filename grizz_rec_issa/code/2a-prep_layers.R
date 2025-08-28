library(terra)
library(sf)
library(tidyverse)

# Set base path and layer filters
lyr.base <- "/Users/claytonlamb/Dropbox/Documents/University/Geographic_Data/BCAB_Layers"
want <- c("landcover/nfls/binary/Broadleaf", 
          "landcover/nfls/binary/Coniferous",
          "landcover/nfls/binary/Exposed_Barren_Land",
          "landcover/nfls/binary/Herb",
          "landcover/nfls/binary/Mixedwood",
          "landcover/nfls/binary/Rock_Rubble",
          "landcover/nfls/binary/Shrubland",
          "landcover/nfls/binary/Wetland",
          "veg/evi", "buildings", "topo")

# Load buffered AOI 
aoi <- vect("data/aoi_buffered.shp")

# Load raster template
template <- rast("data/rast_template.tif")

# Step 1: List all GeoTIFFs in desired folders
folder_pattern <- paste(want, collapse = "|")
tif_paths <- list.files(
  path = lyr.base,
  full.names = TRUE,
  recursive = TRUE,
  pattern = "\\.tif{1,2}$",
  ignore.case = TRUE
) %>%
  str_subset(folder_pattern)%>%
  ##Remove files with years > 2015
discard(~ {
  years <- str_extract_all(.x, "\\d{4}")[[1]]
  any(as.numeric(years) < 2016)
})

# Step 2: Categorize rasters
df_paths <- tibble(path = tif_paths) %>%
  mutate(category = if_else(str_detect(path, "/annual/"), "annual", "static"))

# Helper function to resample to template
process_raster <- function(path, template, aoi) {
  r <- rast(path)
  if (crs(r) != crs(template)) {
    r <- project(r, crs(template))
  }
  r <- resample(r, template, method = "near")
  r <- mask(r, aoi)
  return(r)
}

# Step 3: Process and align static rasters
static_stack <- df_paths %>%
  filter(category == "static") %>%
  pull(path) %>%
  map(~ process_raster(.x, template, aoi)) %>%
  rast()

# Step 4: Process and align annual rasters
annual_stack <- df_paths %>%
  filter(category == "annual") %>%
  pull(path) %>%
  map(~ process_raster(.x, template, aoi)) %>%
  rast()

# Optional: name the layers based on filenames
# names(static_stack) <- df_paths %>%
#   filter(category == "static") %>%
#   pull(path) %>%
#   tools::file_path_sans_ext() %>%
#   basename()
# 
# names(annual_stack) <- df_paths %>%
#   filter(category == "annual") %>%
#   pull(path) %>%
#   tools::file_path_sans_ext() %>%
#   basename()


# Define output folders
dir.create("data/spatial/covs/annual", recursive = TRUE, showWarnings = FALSE)
dir.create("data/spatial/covs/static", recursive = TRUE, showWarnings = FALSE)

# Export annual layers
for (i in 1:nlyr(annual_stack)) {
  lyr_name <- names(annual_stack)[i]
  writeRaster(
    annual_stack[[i]],
    filename = file.path("data/spatial/covs/annual", paste0(lyr_name, ".tif")),
    overwrite = TRUE,
    wopt = list(gdal = c("COMPRESS=LZW", "TILED=YES"))
  )
}

# Export static layers
for (i in 1:nlyr(static_stack)) {
  lyr_name <- names(static_stack)[i]
  writeRaster(
    static_stack[[i]],
    filename = file.path("data/spatial/covs/static", paste0(lyr_name, ".tif")),
    overwrite = TRUE,
    wopt = list(gdal = c("COMPRESS=LZW", "TILED=YES"))
  )
}
