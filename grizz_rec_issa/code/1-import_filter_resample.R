

#################################################################
#
# Study: Pine and Quintette LPU caribou iSSA/connectivity
# Script title: Import raw telemetry data and filter
# Author: Eric Palm
#
#################################################################

require(lubridate)
require(amt)
require(sf)
require(purrr)
require(stringr)
require(readr)
require(readxl)
require(ggplot2)
library(terra)
library(mapview)
require(dplyr)

##set buffer
buf_dist <- 10000

# sanity check: buf_dist must exist
if (!exists("buf_dist")) stop("`buf_dist` must be defined before sourcing 01_prep_telemetry.R")


# AOI
aoi <- st_point(c(-115.0631, 49.5040)) %>% ##fernie
  st_sfc(crs = 4326)%>%
  st_transform(crs = 3005)%>%  # BC Albers
  st_buffer(dist = buf_dist)%>%  # 10,000 meters = 10 km, keep smaller radius so no carcass pits
  st_as_sf

#mapview(aoi)

st_write(aoi, "data/aoi.shp", delete_dsn = TRUE)


# make raster template
aoi_vect <- aoi%>%
  st_buffer(5000)
st_write(aoi_vect, "data/aoi_buffered.shp", delete_dsn = TRUE)
aoi_vect <- vect(aoi_vect)

# Create a 30m template raster aligned to AOI
template <- rast(
  extent = ext(aoi_vect),
  resolution = 30,
  crs = crs(aoi)
)

values(template) <- NaN

writeRaster(template, "data/rast_template.tif", overwrite=TRUE)


r_elev <- rast("/Users/claytonlamb/Dropbox/Documents/University/Geographic_Data/BCAB_Layers/static/topo/bc_ab_100km_dem.tiff")
# library(mapview)
# mapview(aoi)
# Import raw telemetry data
telem_raw <-
  read_csv("data/telem/EVcollar_Relocs_raw.csv")



##make spatial
telem_spat <- telem_raw%>%
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326)%>%
  st_transform(3005)%>%
  cbind(st_coordinates(.))%>%
  mutate(
    month = month(DateTime),
    elev = terra::extract(r_elev, .)[,2]  # assuming you have a SpatRaster of elevation
  ) 


##remove dens
trk_dens <- telem_spat%>%
  filter(month %in% c(1, 2, 3),  # Dec to March
              elev > 1500)%>%
  group_by(Name, year(DateTime))%>%
  summarise(geometry = st_union(geometry)) %>%   # collapse points into one group geometry
  st_centroid()%>%
  rename(den.name=Name)%>%
  st_buffer(300)

mapview(trk_dens)


telem_spat.den <- st_join(telem_spat, trk_dens, left=TRUE)%>%
  filter(!(Name==den.name & month%in%c(10:12,1:5) & !is.na(den.name)))



# Make track
telem_trk <- telem_spat.den%>%
  amt::make_track(., .x = X, .y = Y, .t = DateTime, id = Name, month, elev, crs = 3005) 
  

# remove dups
all_trk_no_dups <-
  telem_trk %>% 
  dplyr::group_by(id, t_) %>% 
  dplyr::slice(1) %>% 
  dplyr::ungroup() %>% 
  `class<-`(class(telem_trk))






# Flag fast steps using a 10 km/hr speed threshold in a 1 hour period as the delta.
delta <- amt::calculate_sdr(speed = 10, time = lubridate::period("1 hour"), speed_unit = c("km/h"))


# I added simple seasons just as an example and to make models smaller for the example 
# not filtering by season for now
all_trk_filtered <-
  all_trk_no_dups %>% 
  dplyr::nest_by(id) %>% 
  dplyr:::mutate(fast = list(amt::flag_fast_steps(data, delta = delta))) %>% 
  tidyr::unnest(fast) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(fast_step_ == F) %>% 
  dplyr::select(-c(data, fast_step_)) %>% 
  dplyr::arrange(id, t_) %>% 
  `class<-`(class(telem_trk)) #%>% 
  #dplyr::mutate(season = dplyr::if_else(lubridate::month(t_) %in% c(1:3, 11:12), "winter", "summer"))
  

# Histogram showing common time lag for each animal, across all animals
all_trk_filtered %>% 
  dplyr::group_by(id) %>% 
  dplyr::mutate(t_lag = round(as.numeric(difftime(t_, lag(t_), units = "hours")), 0)) %>% 
  dplyr::group_by(id, t_lag) %>%
  dplyr::summarize(n=n()) %>% 
  dplyr::slice(which.max(n)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(t_lag < 30) %>% 
  ggplot(., aes(x = t_lag)) +
  geom_histogram() +
  theme_classic() +
  labs(x = "Time lag (hours)")




# For this example, I required each animal to have at least 30 unique days of locations to stay in the dataset,
# but you could change this. This chunk of codes identifies the animals that will be retained
# using that criteria. That way I can calculate the population-level step length and turning angles
# using only these retained animals in the next step. Can nest by season if desired.
to_keep <- 
  all_trk_filtered %>% 
  dplyr::nest_by(id) %>%
  dplyr::mutate(steps = list(amt::track_resample(data, rate = lubridate::hours(2), tolerance = lubridate::minutes(90)) %>%
                               amt::filter_min_n_burst(min_n = 3))) %>% 
  tidyr::unnest(steps) %>% 
  dplyr::select(-data) %>% 
  dplyr::filter(dplyr::n_distinct(lubridate::as_date(t_)) >= 30) %>% 
  dplyr::ungroup() 

# This saves the population-level (data pooled across all animals) parameters for step length
# (gamma distribution: scale, shape) and turning angle (von Mises: kappa), which I'll use when generating
# available steps. You could also generate steps using parameters fit to each animal or animal-season. Whatever you pick
# you just want to do the same thing when simulating.
move_params <-
  all_trk_filtered  %>% 
  dplyr::filter(id %in% unique(to_keep$id)) %>% 
  dplyr::nest_by(id) %>%
  dplyr::mutate(steps = list(amt::track_resample(data, rate = lubridate::hours(2), tolerance = lubridate::minutes(90)) %>%
                               amt::filter_min_n_burst(min_n = 3) %>%
                               amt::steps_by_burst())) %>% 
  dplyr::select(-data) %>% 
  tidyr::unnest(cols = steps) %>% 
  group_by()%>%
  dplyr::summarize(d_gamma = list(amt::fit_distr(sl_, "gamma")),
                   d_von_mises = list(amt::fit_distr(ta_, "vonmises"))) %>%
  dplyr::ungroup()%>%
  mutate(pop="Elk Valley")
 

hist(rgamma(100000, shape = move_params$d_gamma[[1]]$params$shape, scale = move_params$d_gamma[[1]]$params$scale))
library(circular)

# Extract parameters

# Step 1: Create theta from -π to π
theta <- seq(-pi, pi, length.out = 360)

# Step 2: Compute von Mises density
kappa <- move_params$d_von_mises[[1]]$params$kappa
mu <- move_params$d_von_mises[[1]]$params$mu  # assumed to be 0

density_vals <- dvonmises(theta, mu = mu, kappa = kappa)

# Step 3: Build tibble for ggplot
vm_df <- tibble(
  theta = theta,
  density = density_vals
)

# Step 4: Plot
ggplot(vm_df, aes(x = theta, y = density)) +
  geom_line() +
  labs(x = "Angle (radians)", y = "Density", title = "Von Mises Distribution (Centered at 0)") +
  theme_minimal() +
  ylim(0, 0.3) +
  scale_x_continuous(
    breaks = c(-pi, -pi/2, 0, pi/2, pi),
    labels = c("-π", "-π/2", "0", "π/2", "π")
  )


# Resample the tracks to one location every 2 hours and generate 15 available locations
# per used location; you can cut this down to 10 available:used after cropping data by spatial extent
# And filtering out available step lengths that are longer than the maximum used step length.
# This uses the movement parameters from above in the amt::random_steps() function

# Set seed for reproducibility 
set.seed(537)

df_resampled <- 
  all_trk_filtered %>% 
  dplyr::filter(id %in% unique(to_keep$id), !id%in%c("Leslie", "Honey")) %>% 
  dplyr::nest_by(id) %>%
  mutate(pop="Elk Valley")%>%
  dplyr::inner_join(., move_params, by="pop") %>% 
  dplyr::mutate(steps = list(amt::track_resample(data, rate = lubridate::hours(2), tolerance = lubridate::minutes(70)) %>%
                               amt::filter_min_n_burst(min_n = 3) %>% 
                               amt::steps_by_burst() %>% 
                               amt::random_steps(n_control = 500, sl_distr = d_gamma, rand_ta = random_numbers(make_unif_distr(), 100000)))) %>%
  dplyr::select(-c(data, d_gamma, d_von_mises)) %>% 
  tidyr::unnest(cols = steps) %>% 
  dplyr::ungroup() 

# Find the longest real step length in the GPS data. Use this to filter longer 
# available steps, which can happen due to long tail of step length gamma distribution.
# Again, there are still implausible steps included due to pen locations, so this isn't realistic at the moment.
sl_real_max <-
  df_resampled %>% 
  dplyr::filter(case_ == T) %>% 
  dplyr::summarize(max = max(sl_)) %>% 
  dplyr::pull(max)

# Drop available locations that have longer step lengths than the longest
# real step length. Make sure each stratum has more than one location.
# Add variables for step length (km) and log step length (km), which makes their
# magnitudes similar to centered/scaled covariates in the model and will aid
# model convergence. We will not center and scale the movement parameters.
df_adjusted <-
  df_resampled %>% 
  dplyr::filter(!sl_ > sl_real_max,
                !is.na(x2_)) %>% 
  dplyr::group_by(id, t2_) %>% 
  dplyr::filter(dplyr::n() > 1) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(sl_km = dplyr::if_else(sl_ < 1, .001, sl_/1000),
                log_sl_km = log(sl_km),
                stratum = stringr::str_c(id, t2_, sep = "_"))

# Cropping resampled data by study area extent, which might eliminate some available locations
# that fall outside of this area. This spatial filter step and removing super long step lengths
# may have produced different sized strata. 
# Here I make sure each stratum has exactly 10 available locations by randomly selecting 5.
# Usually I'd generate 15 and subsample 10, but you'll be modifying this anyway so this is a faster example.
# I'll create a "study_extent.shp" file and drop it in the folder. I haven't been using Quintette
# data so my study extent is too small for you. But I'll be expanding all my covariate layers to a bigger
# extent and you can use as many of those as you want.

set.seed(1234)

df_used <- dplyr::filter(df_adjusted, case_ == 1)%>%
  sf::st_as_sf(., coords = c("x2_", "y2_"), crs = 3005, remove = F) %>% 
  sf::st_filter(., aoi)


df_avail <- df_adjusted %>% 
  dplyr::filter(case_ == 0,
                stratum%in%df_used$stratum) %>% 
  sf::st_as_sf(., coords = c("x2_", "y2_"), crs = 3005, remove = F) %>% 
  sf::st_filter(., aoi) %>% 
  sf::st_drop_geometry() %>% 
  dplyr::group_by(id, t2_) %>% 
  dplyr::slice_sample(n = 10) %>% 
  dplyr::ungroup()

##make sure there always 10 per strata
avail.check <- df_avail%>%
  group_by(stratum)%>%
  add_count()%>%
  filter(n!=10)%>%
  ungroup%>%
  sf::st_as_sf(., coords = c("x2_", "y2_"), crs = 3005, remove = F)

#stopifnot(avail.check$stratum%>%unique%>%length()==0)

avail.check$stratum%>%unique%>%length()/df_avail$stratum%>%unique%>%length() ##0.04% of steps

strata_remove <- avail.check$stratum%>%unique

# mapview(df_adjusted%>%filter(stratum%in%strata_remove[10])%>%  sf::st_as_sf(., coords = c("x2_", "y2_"), crs = 3005, remove = F),
#         zcol="case_")+mapview(aoi)


df_ssf <-df_used%>%
  filter(!stratum%in%strata_remove)%>%
  dplyr::bind_rows(df_avail%>%
                     filter(!stratum%in%strata_remove)) 


# select months to use that correspond to high trail use
df_ssf <- df_ssf%>%
  mutate(month_num = month(t1_))

hist(df_ssf$month_num)

df_ssf_summer <-
  df_ssf %>%
  filter(month_num %in% 5:10)


##keep animals with at least 50 locations
fernie.bears <-df_ssf_summer%>%
  group_by(id)%>%
  dplyr::filter(case_ == 1)%>%
  count()%>%
  filter(n>50)%>%
  pull(id)

df_ssf_summer_fernie <- df_ssf_summer%>%
  filter(id%in%fernie.bears)

df_ssf_summer_fernie$sl_%>%quantile()

##make sure there are used locations for each strata
issues <- df_ssf_summer_fernie%>%
  tibble%>%
  group_by(stratum)%>%
  summarise(
    used = sum(case_ == 1),
    mean_case = mean(case_))%>%
  filter(used!=1 | mean_case!=(1/11))

#stopifnot(nrow(issues)==0)

df_ssf_summer_fernie <- df_ssf_summer_fernie%>%
  filter(!stratum%in%issues$stratum)

# Save the data frame without annotated habitat covariates
saveRDS(df_ssf_summer_fernie, file.path("data", 'derived', "ssf_2_hrs_no_covs.rds"))











df_adjusted %>% 
       dplyr::filter(case_ == 1)

df_ssf_summer_fernie %>% 
  dplyr::filter(case_ == 1)

df_ssf_summer_fernie$id%>%unique()%>%length

a <- df_ssf_summer_fernie%>%
  group_by(id)%>%
  dplyr::filter(case_ == 1)%>%
  count()

mean(a$n)

range(df_ssf_summer_fernie$t1_)
