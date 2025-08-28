# 00_run_all_buffers.R
# --------------------
# This “master” script loops over buffer distances, calls each step in order,
# and organizes outputs into buffer‐specific folders.

# 1) Define the buffer distances (in metres) you want to try:
buffer_list <- c(5000, 10000, 15000, 17000, 23000) 
# (i.e., 5 km, 10 km, 15 km, 20 km)

# 2) Loop over each buffer size:
for (buf in buffer_list) {
  
  message("==== Starting run with buf_dist = ", buf, " metres (", buf/1000, " km) ====")
  
  # 3) Expose `buf_dist` so that sourced scripts see it:
  buf_dist <- buf
  
  # Optionally create a parent folder to collect all outputs for this buffer:
  out_base <- file.path("analysis_runs", paste0("buf_", buf/1000, "km"))
  dir.create(out_base, showWarnings = FALSE, recursive = TRUE)
  
  #
  # 4) Step 1: Prep telemetry
  #
  message("-> Running step 1: 01_prep_telemetry.R")
  # Source the telemetry‐prep script. It will write to data/aoi_buf_<X>km/…
  source("code/1-import_filter_resample.R")
  

  #
  # 5) Step 2: Prep rasters
  #
  message("-> Running step 2: 02_extract_rasters.R")
  source("code/2c-extract_covariates.R")
  
  
  #
  # 6) Step 3: Run iSSA
  #
  message("-> Running step 3: 03_run_issa.R")
  # The iSSA script should read from data/aoi_buf_<X>km/… and write model results.
  source("code/3-fit_iSSF.R")
  
  message("==== Finished run with buf_dist = ", buf, " metres ====\n\n")
}