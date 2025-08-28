

################################################################
#
# Title: Simulate UD from mixed-effects iSSF model 
# Author: Julie Turner adapted from Eric Palm 
#
###############################################################

# Load required packages
library(readr)
library(tidyr)  
library(stringr)     
require(raster)
library(sf)                  
library(terra)
library(dplyr)
require(glmmTMB)
require(tictoc)
require(data.table)
require(sfheaders)
require(foreach)
require(caret)
require(stringr)

# data locations
covs30m <- file.path('data', 'covs', '30m')
covs90m <- file.path('data', 'covs', '90m')
shp <- file.path('data', 'shp')
tmpSim90m <- file.path(covs90m, "temp_sim")

# Import function for simulating coefficients from full variance-covariance
# matrix of the mode. Note that including random slopes on  
# movement covariates and interactions with those movement covariates
# might result in negative step lengths in some instances.

source("code/functions/simulate_coefs_random.R")

# Define season for model and simulation
#m_season <- "winter"

# define name of simulation
sim.name <- "disturbances_trueStart_2yrs"

# Import model from which to simulate
model <- readRDS(file.path('data', 'derived', 'mods',"mod_disturbances.rds")) 

# Import input dataset
df_unscaled <- readRDS(file.path('data', 'derived',"ssf_5_9_hrs_unscaled.rds")) #%>% 
 # dplyr::filter(season == m_season)


# Calculate means of all variables to use for centering/scaling the raster values so the 
# raster values are scaled the same way as the variables were in the models.
# Save the values as a list.
# You have to manually include the "_start" variables because they only appear in interactions
# and not on their own. Probably a better way to do this step.
means <- df_unscaled %>%
  dplyr::select(dplyr::any_of(rownames(data.frame(glmmTMB::fixef(model)$cond))),
                cc_start, tpi_5_start, tri_start) %>% 
  dplyr::summarize(dplyr::across(dplyr::everything(), ~mean(.))) %>% 
  as.list()

# Repeat the same step as above but calculate the standard deviation for scaling the raster values.
sds <- df_unscaled %>%
  dplyr::select(dplyr::any_of(rownames(data.frame(glmmTMB::fixef(model)$cond))),
                cc_start, tpi_5_start, tri_start) %>% 
  dplyr::summarize(dplyr::across(dplyr::everything(), ~sd(.))) %>% 
  as.list()

# This saves the scaling parameters of the delta elev variable 
# from the input dataset, so we can apply these to center and scale this variable
# while we're simulating.
delta_elev_scaled <- df_unscaled %>%
  dplyr::select(delta_elev_sim = delta_elev) %>% 
  caret::preProcess(., method = c("center", "scale"))

###################### CENTER AND SCALE RASTERS ###################
# Import and scale all rasters in model and save to temporary folder.
# This way we can call these rasters using terra package directly within the 
# foreach loop below, because we can't import terra spatRaster objects into the 
# foreach loop (terra uses C++ which apparently doesn't work with parallelization)
# Probably could automate this step to avoid hard coding.

terra::rast(file.path(covs90m, "tri.tif")) %>% 
  terra::scale(., center = means$tri_end, scale = sds$tri_end) %>% 
  `names<-`("tri_end") %>% 
  terra::writeRaster(., file.path(covs90m, "temp_sim", "tri_end.tif"), datatype = "FLT4S", overwrite = T)

# terra::rast(file.path(covs90m, "tri.tif")) %>% 
#   terra::scale(., center = means$tri_start, scale = sds$tri_start) %>% 
#   `names<-`("tri_start") %>% 
#   terra::writeRaster(., file.path(covs90m, "temp_sim", "tri_start.tif"), datatype = "FLT4S", overwrite = T)

# terra::rast("covs/90m/tpi_5.tif") %>% 
#   terra::scale(., center = means$tpi_5_end, scale = sds$tpi_5_end) %>% 
#   `names<-`("tpi_5_end") %>% 
#   terra::writeRaster(., "covs/90m/temp_sim/tpi_5_end.tif", datatype = "FLT4S", overwrite = T)
# 
# terra::rast("covs/90m/tpi_5.tif") %>% 
#   terra::scale(., center = means$tpi_5_start, scale = sds$tpi_5_start) %>% 
#   `names<-`("tpi_5_start") %>% 
#   terra::writeRaster(., "covs/90m/temp_sim/tpi_5_start.tif", datatype = "FLT4S", overwrite = T)
# 
# terra::rast("covs/90m/heatload.tif") %>% 
#   terra::scale(., center = means$heatload_end, scale = sds$heatload_end) %>% 
#   `names<-`("heatload_end") %>% 
#   terra::writeRaster(., "covs/90m/temp_sim/heatload_end.tif", datatype = "FLT4S", overwrite = T)

terra::rast(file.path(covs90m, 'cc',  "cc_2024.tif")) %>% 
  terra::scale(., center = means$cc_end, scale = sds$cc_end) %>% 
  `names<-`("cc_end") %>% 
  terra::writeRaster(., file.path(covs90m, "temp_sim", "cc_end.tif"), datatype = "FLT4S", overwrite = T)

terra::rast(file.path(covs90m, 'cc',  "cc_2024.tif")) %>% 
  terra::scale(., center = means$cc_start, scale = sds$cc_start) %>% 
  `names<-`("cc_start") %>% 
  terra::writeRaster(., file.path(covs90m, "temp_sim", "cc_start.tif"), datatype = "FLT4S", overwrite = T)

# This applies the negative exponential decay for distance to before scaling.
terra::rast(file.path(covs90m, "d_burn_2024.tif")) %>% 
  terra::app(., function(x) log(x + 1)) %>% 
  terra::scale(., center = means$log_d_burn_end, scale = sds$log_d_burn_end) %>% 
  `names<-`("log_d_burn_end") %>% 
  terra::writeRaster(., file.path(covs90m, "temp_sim","log_d_burn_end.tif"), datatype = "FLT4S", overwrite = T)

terra::rast(file.path(covs90m, "d_cutblock_2024.tif")) %>% 
  terra::app(., function(x) log(x + 1)) %>% 
  terra::scale(., center = means$log_d_cutblock_end, scale = sds$log_d_cutblock_end) %>% 
  `names<-`("log_d_cutblock_end") %>% 
  terra::writeRaster(., file.path(covs90m, "temp_sim", "log_d_cutblock_end.tif"), datatype = "FLT4S", overwrite = T)

terra::rast(file.path(covs90m, "d_transmission_2024.tif")) %>% 
  terra::app(., function(x) log(x + 1)) %>% 
  terra::scale(., center = means$log_d_transmission_end, scale = sds$log_d_transmission_end) %>% 
  `names<-`("log_d_transmission_end") %>% 
  terra::writeRaster(., file.path(covs90m, "temp_sim", "log_d_transmission_end.tif"), datatype = "FLT4S", overwrite = T)

terra::rast(file.path(covs90m, "d_seismic_2024.tif")) %>% 
  terra::app(., function(x) log(x + 1)) %>% 
  terra::scale(., center = means$log_d_seismic_end, scale = sds$log_d_seismic_end) %>% 
  `names<-`("log_d_seismic_end") %>% 
  terra::writeRaster(., file.path(covs90m, "temp_sim", "log_d_seismic_end.tif"), datatype = "FLT4S", overwrite = T)
####################################################################################################

###############################
# Import semi-permeable barrier shapefiles

sf_barrier_transmission <- sf::read_sf(file.path(shp, "barrier_transmission.shp")) 
###############################

###############################
# Extract parameters for turning angle (von Mises concentration parameter, kappa) and
# step length (gamma scale and shape) from the input dataset and correct season.
# Use these to create random steps in simulation. In this case, I generated random
# turn angles in the SSF data frame using a uniform distribution because there
# were so many instances of "turning around" that a von Mises was a particularly
# bad fit. Uniform might not be much better though.

# Only use used locations for this.
df_used <-
  df_unscaled %>%
  dplyr::filter(case_ == T) 

# ta_kappa <-
#   amt::fit_distr(df_used$ta_, "vonmises") %>%
#   .$params %>% .$kappa

sl_km_scale <-  
  amt::fit_distr(df_used$sl_km, "gamma") %>% 
  .$params %>% .$scale

sl_km_shape <-  
  amt::fit_distr(df_used$sl_km, "gamma") %>% 
  .$params %>% .$shape

# find end points for potential starts
df_end <- setDT(df_used)[year>=2020, .SD[which.max(t1_)], by = id]

###########################################################

# Find the 99.9th percentile of (centered, scaled) TRI values in the real GPS data
# We don't allow simulated animals to start in pixels higher than this 
# because sometimes with the negative coefficient for log_sl:tri_end, simulated
# animals can get "stuck" in cliff areas, as the model might say that step lengths are
# nearly zero if terrain is extremely rugged.
tri_max <- df_unscaled %>% 
  dplyr::mutate(tri_end = as.numeric(scale(tri_end)[,1])) %>% 
  dplyr::filter(case_ == T) %>% 
  dplyr::select(tri_end) %>% 
  dplyr::pull(tri_end) %>% 
  quantile(., 0.999) %>% 
  as.numeric()

# Define the simulation start zones. You might want to use a 
# 2nd-order RSF to determine this, but I didn't do that here. Instead, I just
# prevented start locations from being in water or cliffs.
ssf_start <- terra::rast(file.path(covs90m, "r_study.tif")) %>%
  terra::mask(., rast(file.path(covs90m, "temp_sim", "tri_end.tif")) %>% 
                terra::classify(., cbind(tri_max, Inf, NA))) %>% 
  terra::mask(., terra::rast(file.path(covs90m, "water_01.tif")), maskvalues = 1) %>% 
  raster::raster()

# This is the simulation start extent. It is the full study area with a -1km
# buffer, so simulated animals can move outwards and still be in the study area.
# This step might not really be necessary.
start_extent <- terra::rast(file.path(covs90m, "r_study.tif")) %>% 
  terra::as.polygons() %>% 
  terra::buffer(., -1000) %>% 
  sf::st_as_sf()

###########################################
#### 2.  Simulate movement paths       ####
###########################################

# Define the number of "monster" loop' I set this to be equal to the number of cores I'm using,
# but you could change this depending on how many cores you have available.
# You'll have to play with these values depending on your goals and computing power.

gc()
parallelly::availableCores(constraints = "connections")
#32, picking fewer cores to not overload -- 16
n_monster <- 100               
#n_start_locs <- 1000     # The number of random paths within one loop  # this is what i did for random start disturbance model
# start from where ended, but have many start points to create many tracks
start_rep <- 15
n_start_locs <- uniqueN(df_end$id)*start_rep
n_rand <- 20                # Number of random locations to generate for each step  # this is what i did for random start disturbance model
#n_rand <- 100                # Number of random locations to generate for each step  
n_days <- 365*2               # Number of days to simulate
reloc_int <- 7 # We're simulating one location every seven hours for n_days days
nt <- round(n_days / (reloc_int/24), 0) # number of steps per path

# Total number of locations generated (actual number will be less due to 
# NAs when enough candidate locations fall outside study area)
n_total <- n_monster * n_start_locs * nt  

# Sequence variable to be used for the monster loops
iseq <- 1:n_monster

# Define the number of cores for the foreach loop and initiate the cluster
#noCLS <- n_monster   # this is what i did for random start disturbance model
noCLS <- 16
cl <- parallel::makeCluster(noCLS)
doParallel::registerDoParallel(cl)

# Create a log text file so you can track progress 
writeLines(c(" "), file.path("results", paste0("ssf_", sim.name, "_log.txt")))

gc()
# Use tictoc to record the overall time for the simulation.
tictoc::tic()

# Start the foreach loop and remember to include all packages here that you
# will use within the loop. Recall that you can't import previously created objects that are
# terra SpatRasters into the foreach loop, but you can create and call them inside the loop.
path_monster <- foreach::foreach(yyy=iseq, .combine = rbind,
                                 .packages = c("sf", "Rfast", "dplyr", "tidyr", "data.table", 
                                               "sfheaders", "terra", "raster", "caret", "stringr")) %dopar% 
  {
    # Create an empty data frame with the necessary columns, including path id, stratum,
    # previous and current locations.
    df_path <- dplyr::tibble(path_id = 1:n_start_locs,
                             i_strata = 1,
                             x_prev = as.numeric(NA),
                             y_prev = as.numeric(NA),
                             x = as.numeric(NA),
                             y = as.numeric(NA)) 
    
    # Create random start locations; I create too many locations first and 
    # then subsample the exact amount I want after making sure they intersect
    # with the start extent.
    
    # xy <- raster::sampleRandom(ssf_start, n_start_locs*3, method = "random", na.rm = T, xy = T)[,-3] %>%
    #   sfheaders::sf_point() %>%
    #   sf::st_set_crs(., 3005) %>%
    #   dplyr::filter(lengths(sf::st_intersects(., start_extent)) > 0) %>%
    #   dplyr::slice_sample(., n = n_start_locs)
####SWITCH ----
    #### select last end points to be potential start points for sims
    xy <- df_end[,.(x1_, y1_)] %>%
      slice(rep(1:n(), each = start_rep)) %>%
      st_as_sf( coords = c("x1_", "y1_")) %>%
      st_set_crs(3005)
    
    
    # Create an initial location for each simulated path using the start locations above
    # Create a random initial bearing
    # Then simulate a set of coefficients using the 'sim_coef_ran' function, which 
    # draws from multivariate normal distributions of each coefficient in the model,
    # and includes random slope variation (except for some variables in this example).
    # So each path (or simulated animal) will have its
    # own set of coefficients throughout the simulation.
    tmp_0 <- df_path %>% 
      mutate(i_strata = 1,
             x_prev = st_coordinates(xy)[,1], 
             y_prev = st_coordinates(xy)[,2], 
             direction_prev = round(runif(n_start_locs, -pi, pi), 3),
             good = 1) %>% 
      dplyr::bind_cols(., sim_coef_ran(model, n = nrow(.)))
    
    # start simulations by each individual step across all animals
    path_list <- vector(mode = 'list', length = nt + 1)
    
    for (i in 1:nt){
      
      # Add each step to the log file created earlier so you can monitor progress 
      cat(paste("Starting iteration: Monster loop", yyy, "smaller loop",  i, lubridate::now(), "\n"),
          file = file.path("results", paste0("ssf_", sim.name, "_log.txt")), append=TRUE)
      
      tmp <- tmp_0 %>% filter(good == 1)
      
      # Start the n_keep loop (this separate loop might not be necessary but I've left it in from Jesse W.)
      n_keep <- nrow(tmp)
      if(n_keep > 0){
        
        # Add the simulated strata number (starting at 1 and going to nt)
        tmp <- tmp %>% mutate(i_strata = i)
        
        # At all the start points, extract the covariates that are used in interactions with the 
        # movement parameters (log_sl and cos_ta), so the distribution of step lengths and turning angles
        # depends on the habitat conditions at the start of steps.
        # I'm using the gamma distribution for the step lengths, which has shape and scale parameters.
        # Also, I only used log_sl_km (modifies the gamma distribution shape parameter) in this model, but you could also use sl_km.
        # The code for modifying sl_km (scale) is a bit more complicated but can be found in amt's source code on github.
        # Most importantly, this next code generates all the proposed (candidate) steps.
        tmp_2 <- tmp %>%
          dplyr::bind_cols(., terra::extract(list.files(tmpSim90m, full.names = T, pattern = "start") %>% rast(), cbind(.$x_prev, .$y_prev))) %>% 
          # Extract elevation at the start of steps so you can calculate delta elevation later 
          dplyr::mutate(elev_start = terra::extract(terra::rast(file.path(covs90m, "elev.tif")), cbind(.$x_prev, .$y_prev))[,1]) %>% 
          dplyr::rowwise() %>% 
          # setting adjusted shape, this is hard coded dependent on model right now
          dplyr::mutate(step = list(rgamma(n=n_rand,
                                           scale = sl_km_scale,
                                           shape = sl_km_shape + `log_sl_km` +
                                             #tri_start*`log_sl_km:tri_start` +
                                             cc_start*`log_sl_km:cc_start`# +
                                             #tpi_5_start*`log_sl_km:tpi_5_start`
                                           ))) %>% 
          tidyr::unnest(step) %>% 
          dplyr::filter(!is.na(step)) %>%
          # Update the concentration parameter for a von Mises distribution
          # And if it's less than 0, flip the mean direction to pi (backwards)
          # Could use the ta_kappa object created earlier if you had used a von Mises distribution
          # to generate available turn angles for your ssf data frame
          dplyr::mutate(conc = `cos_ta` + log(step)*`cos_ta:log_sl_km`, # + ta_kappa
                        mu = dplyr::if_else(conc < 0, pi, 0)) %>% 
          dplyr::rowwise() %>% 
          # Use the absolute value of the concentration parameter so you can generate turn angles.
          # centered around going backwards when conc is negative.
          # Otherwise this function will draw from a uniform distirbution if the conc is negative.
          dplyr::mutate(angle = as.numeric(Rfast::rvonmises(1, m=mu, k = abs(conc)))) %>% 
          dplyr::ungroup() %>% 
          dplyr::mutate(angle = if_else(angle > pi, angle - (2 * pi), angle),
                        direction = direction_prev + angle,
                        direction = dplyr::if_else(direction < -2*pi, direction + 2*pi,
                                                   dplyr::if_else(direction > 2*pi, direction - 2*pi, direction)),
                        dx = step*1000 * cos(direction), # *1000 converts meters to km
                        dy = step*1000 * sin(direction), # *1000 converts meters to km
                        x = round(x_prev + dx), 
                        y = round(y_prev + dy),
                        # determines whether candidate locations are within study area
                        in_study = terra::extract(terra::rast(file.path(covs90m, "r_study.tif")), cbind(x, y))[,1],
                        delta_elev_sim = abs(terra::extract(rast(file.path(covs90m, "elev.tif")), cbind(x, y))[,1] - elev_start)) %>%  
          dplyr::group_by(path_id) %>% 
          # Calculates proportion of locs are within study area
          dplyr::mutate(prop_in_study = sum(in_study) / n_rand) %>% 
          dplyr::ungroup() %>% 
          dplyr::filter(in_study == 1) %>% # only keeps locations within the study area
          # centers/scales the delta elevation variable based on saved parameters from before (this uses the caret package)
          predict(delta_elev_scaled, .) 
        
        ############# CREATE LINE SEGMENTS
        # Use the 'data.table' and 'sfheaders' packages to quickly create line segments for semi-permeable barrier intersections
        dt <- data.table::as.data.table(tmp_2 %>% dplyr::select(x_prev:y) %>% dplyr::mutate(id = dplyr::row_number()))
        
        ## To use `sfheaders` the data needs to be in long form
        dt1 <- dt[, .(id, lon = x_prev, lat = y_prev)]
        dt2 <- dt[, .(id, lon = x, lat = y)]
        
        ## Add on a 'sequence' variable so we know which one comes first
        dt1[, seq := 1L ]
        dt2[, seq := 2L ]
        
        ## put back together
        dt <- data.table::rbindlist(list(dt1, dt2), use.names = TRUE)
        data.table::setorder(dt, id, seq)
        
        # Create an sf line object, where every individual step is a line    
        line_sf <- sfheaders::sf_linestring(
          obj = dt
          , x = "lon"
          , y = "lat"
          , linestring_id = "id"
        ) %>% 
          sf::st_set_crs(3005) 
        ########################################################## 
        
        # Get probability of each step, using the "end" covariates and the line segment intersections
        # This step could be more automated to avoid manual hard coding
        # Then calculate the exponentiated linear predictor by multiplying extracted covariate values
        # by the simulated coefficients for that covariate.
        # Then calculate the relative probability of each candidate step being chosen by dividing exp_lp
        # by the sum of all exp_lp values within that step (for that path)
        tmp_3 <- 
          tmp_2 %>%
          dplyr::bind_cols(., terra::extract(list.files(tmpSim90m, full.names = T, pattern = "end") %>% terra::rast(), cbind(.$x, .$y)) %>%
                             dplyr::rename_with(~stringr::str_c(., "_sim"), dplyr::everything())) %>%
          # Determine whether candidate steps cross (intersect with) barriers and which are in binary variables like water and avalanche terrain.
          # Set values in extremely rugged terrain to NA.
          dplyr::mutate(barrier_transmission_sim = as.numeric(lengths(sf::st_intersects(line_sf, sf_barrier_transmission)) > 0),
                        water_01_end_sim = terra::extract(terra::rast(file.path(covs90m, "water_01.tif")), cbind(tmp_2$x, tmp_2$y))[,1],
                        #av_risk_01_end_sim = terra::extract(terra::rast("covs/90m/av_risk_01.tif"), cbind(tmp_2$x, tmp_2$y))[,1],
                        # Calculate the relative probability of selection (this is hard coded by model right now)
                        exp_lp = exp(cc_end_sim*cc_end +
                                       # tpi_5_end_sim*tpi_5_end +
                                       tri_end_sim*tri_end +# tri_end_sim^2*`I(tri_end^2)` +
                                       # heatload_end_sim*heatload_end +
                                       # av_risk_01_end_sim*av_risk_01_end +
                                       water_01_end_sim*water_01_end +
                                       log_d_seismic_end_sim*log_d_seismic_end +
                                       log_d_transmission_end_sim*log_d_transmission_end +
                                       log_d_burn_end_sim*log_d_burn_end +
                                       log_d_cutblock_end_sim*log_d_cutblock_end +
                                       delta_elev_sim*delta_elev +
                                       barrier_transmission_sim*barrier_transmission))  %>%
          dplyr::filter(exp_lp > 0 & is.finite(exp_lp)) %>%
          dplyr::group_by(path_id) %>%
          dplyr::mutate(prob = exp_lp / sum(exp_lp, na.rm = TRUE)) %>%
          dplyr::ungroup() %>%
          dplyr::filter(prob > 0 & is.finite(prob)) 
        
        
        # For each path, use the slice_sample function (weighted by relative probability of selection)
        # to probabilistically select one location from the candidate set of locations. 
        # If more than 60% of the candidate locations are outside of the study area, terminate the path. 
        if (nrow(tmp_3) > 0){
          tmp_4 <- tmp_3 %>% 
            dplyr::group_by(path_id) %>% 
            dplyr::slice_sample(n = 1, weight_by = prob) %>% 
            dplyr::ungroup() %>% 
            dplyr::arrange(path_id) %>% 
            dplyr::mutate(good = dplyr::if_else(prop_in_study >= 0.6, 1, 0))   
          
          # Quick rounding to make cleaner (and smaller file sizes) output tables
          tmp_4 <- tmp_4 %>% 
            dplyr::mutate(step = round(step),
                          angle = round(angle, 3),
                          direction = round(direction, 3))
          
          # Only retain the necessary information for saving
          tmp_save <- tmp_4 %>% 
            dplyr::select(path_id, i_strata, x, y) 
          
          # Store the selected steps for this iteration in a list
          path_list[[i]] <- tmp_save
          
          # Then set these x,y locations to x_prev and y_prev for use in the next iteration
          # Do the same for the previous direction
          tmp_5 <- tmp_4 %>% 
            dplyr::mutate(x_prev = x,
                          y_prev = y,
                          x = NA,
                          y = NA,
                          direction_prev = direction) %>% 
            dplyr::select(tidyselect::all_of(names(tmp_0)))
          
          # use tmp_0 in next iteration at top of loop
          tmp_0 <- tmp_5  
          
        } 
        
      }  # END of if n_keep > 0
      
    }  # End of i (individual steps)
    
    # Take the list with all saved steps and add the i_monster sequence
    dplyr::bind_rows(path_list) %>% 
      dplyr::mutate(i_monster = yyy, .before = path_id)
  }

# END OF MONSTER LOOP
tictoc::toc()

# Stop the cluster
parallel::stopCluster(cl)

# Delete temporary rasters used in simulation
unlink(file.path(tmpSim90m, "*"))

# Save the path_monster object, which has all the actual simulated locations
saveRDS(path_monster, file.path('results', paste0('path_all_', sim.name, '.rds')))

# Rasterize the locations so pixel values are the number of simulated locations that fall within that pixel
tictoc::tic()

path_monster_sub <- setDT(path_monster)[i_strata>=(nt/2)]
ud <- 
  terra::rasterize(cbind(path_monster_sub$x, path_monster_sub$y), terra::rast(file.path(covs90m, "r_study.tif")), fun = length) %>% 
  terra::classify(., cbind(NA, NA, 0)) %>% 
  terra::mask(., terra::rast(file.path(covs90m, "r_study.tif"))) %>% 
  terra::writeRaster(., file.path('maps', paste0(sim.name, '_ud.tif')), overwrite = T, datatype = "INT2U")

tictoc::toc()

##### END ####
require(amt)

#paths <- setDT(readRDS(file.path('results', paste0('path_all_', sim.name, '.rds'))))
paths <- path_monster_sub

df_id <- paths %>%
  #read_rds(file.path('results', paste0('path_all_', sim.name, '.rds'))) %>%
  dplyr::mutate(path_id = stringr::str_c(i_monster, path_id, sep = "_")) 

# This basically binds two data frames together, one with starting coordinates for steps, and one with ending coordinates
# Then uses the 'sfheaders' package to create the sf linestrings
path_sf <-
  df_id %>%
  dplyr::mutate(seq = 1) %>%
  dplyr::select(-i_monster) %>%
  dplyr::bind_rows(., df_id %>%
                     dplyr::group_by(path_id) %>%
                     dplyr::mutate(x2 = dplyr::lead(x),
                                   y2 = dplyr::lead(y)) %>%
                     dplyr::ungroup() %>%
                     dplyr::filter(!is.na(x2)) %>%
                     dplyr::select(x = x2, y = y2, path_id, i_strata) %>%
                     dplyr::mutate(seq = 2)) %>%
  dplyr::arrange(., path_id, seq, i_strata) %>%
  sfheaders::sf_linestring(.,
                           x = "x",
                           y = "y",
                           linestring_id = "path_id",
                           keep = T) %>%
  sf::st_set_crs(3005)



#path_sf_valid <- st_make_valid(path_sf)




hwy <- st_read(file.path(shp, "barrier_hwy.shp"))

#hwy.intersct <- st_intersects(hwy, path_sf)


#vect.paths <- vect(path_sf)

####

start = as.POSIXct("2023-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%OS")
end = as.POSIXct("2023-12-31 23:59:59", format = "%Y-%m-%d %H:%M:%OS")
times <- data.table(i_strata = 1:nt, 
                    t = seq(start, end, length.out = nt))
paths.times <- merge(paths, times, all.x = T)
paths.times[,id:=paste(path_id, i_monster, sep = "_")]
#paths.simp <- paths.times[, .SD[1], by = .(id, i_strata)]
paths.times[,n:=.N, by = .(id)]

hwy.buff <- st_buffer(hwy, 5000)
coords <- paths.times[,id:= paste(path_id, i_monster, sep = "_")] %>% 
  st_as_sf( coords = c("x", "y")) %>% 
  st_set_crs(3005) 

hwy.points <- coords %>%
  st_join(hwy.buff, join = st_within) %>%
  filter(!is.na(ROAD_CLASS))
plot(hwy.points$geometry)
st_write(hwy.points, file.path('results', paste0('hwy5km_pts_', sim.name, '.shp')))

gc()
# trk <- paths.times %>%
#   group_by(id) %>%
#   amt::make_track(., .x= x, .y = y, t = t, crs = sf::st_set_crs(3005) ) %>%
#   ungroup()


# sf.paths  <- paths.times[n >= 3] %>% 
#   st_as_sf( coords = c("x", "y")) %>% 
#   group_by(id) %>% 
#   arrange(i_strata) %>%
#   summarize(do_union = F) %>%
#   st_cast("LINESTRING") %>% 
#   st_set_crs(3005) 



# xing <- st_intersection(sf.paths, hwy)
# vect.paths <- vect(sf.paths)

gc()
path_sf2 <-
  paths.times[n >= 3] %>%
  dplyr::mutate(seq = 1) %>%
  dplyr::select(-i_monster) %>%
  dplyr::bind_rows(., paths.times[n >= 3] %>%
                     dplyr::group_by(id) %>%
                     dplyr::mutate(x2 = dplyr::lead(x),
                                   y2 = dplyr::lead(y)) %>%
                     dplyr::ungroup() %>%
                     dplyr::filter(!is.na(x2)) %>%
                     dplyr::select(x = x2, y = y2, id, i_strata) %>%
                     dplyr::mutate(seq = 2)) %>%
  dplyr::arrange(., id, seq, i_strata) %>%
  sfheaders::sf_linestring(.,
                           x = "x",
                           y = "y",
                           linestring_id = "id",
                           keep = T) %>%
  sf::st_set_crs(3005)

st_write(path_sf2, file.path('maps', paste0(sim.name, '_paths.shp')), append = F)

gc()
path_vect<- vect(path_sf2) 

xing <- st_intersection(path_sf2, hwy)
writeVector(vect(xing), file.path('maps', paste0(sim.name, '_hwy_xing.shp')))




hwy.buff.rast <- rasterize(vect(hwy.buff), terra::rast(file.path(covs90m, "r_study.tif")))
#xing.rast <- rasterize(vect(xing), hwy.buff.rast, fun = 'length')

paths.ud <- rasterize(path_vect,  terra::rast(file.path(covs90m, "r_study.tif")), fun = 'count')
plot(paths.ud)
writeRaster(paths.ud, file.path('maps', paste0(sim.name, '_paths_ud.tif')))


plot(vect(file.path('maps', paste0(sim.name, '_hwy_xing.shp'))))


paths.hwy.ud  <- mask(paths.ud, hwy.buff.rast, inverse = F)
plot(paths.hwy.ud)
