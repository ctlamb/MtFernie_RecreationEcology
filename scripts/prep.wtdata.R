
library(remotes)
#remotes::install_github("ABbiodiversity/wildrtrax@development")
library(wildrtrax)
library(here)
library(lubridate)
library(tidyverse)
library(tidylog)
source("helpers/ind_det.R")
source("helpers/summarise_detections.R")
source(here::here("helpers", "gg_theme_custom.R"))

Sys.setenv(WT_USERNAME = "lamb.eco.research", WT_PASSWORD = "Wild$%1991")
wt_auth()
my_projects <- wt_get_download_summary(sensor_id = "CAM") %>% filter(organization == "LER")

## Add in regions and also control for season
## control for enviro, trail across/up, tortuosity, etc


# Load data ---------------------------------------------------------------
# load species lookup
sp.lkup <- read_csv(here::here("data/species_lookup.csv"))
# Download the project reports
proj.ids <- c(998, 1401, 1971, 2626)
image <- map_df(
  .x = proj.ids,
  .f = ~ wildrtrax::wt_download_report(
    project_id = .x,
    sensor_id = "CAM",
    report = "image_report",
    weather_cols = FALSE
  ) %>% select(-image_fire)
)

tag <- map_df(
  .x = proj.ids,
  .f = ~ wildrtrax::wt_download_report(
    project_id = .x,
    sensor_id = "CAM",
    report = "tag",
    weather_cols = FALSE
  )
)

location <- map_df(
  .x = proj.ids,
  .f = ~ wildrtrax::wt_download_report(
    project_id = .x,
    sensor_id = "CAM",
    report = "location",
    weather_cols = FALSE
  )
) %>% distinct()

# load species lookup
sp.lkup <- read_csv(here::here("data/species_lookup.csv"))





# Prep WT data ------------------------------------------------

##get trail types prepped
location$location_comments%>%unique()
location <- location%>%
  mutate(location_comments=case_when(location_comments%in%c("Road", "Single Track")~"Rec trail",
                                     location_comments%in%c("Wildlife")~"Wildlife trail"))

## join image and tag data##
df <- image %>%
  select(
    project_id, location_id, location, image_id, image_date_time,image_set_id,
    equipment_make,
    image_comments, image_fov, image_trigger_mode
  ) %>%
  left_join(
    tag %>%
      select(
        project_id, location_id, image_id,
        species_scientific_name, species_common_name, age_class, sex_class, individual_count, tag_comments
      ),
    by = c("project_id", "location_id", "image_id")
  ) %>%
  mutate(species_common_name = replace(species_common_name, species_common_name == "Elk (wapiti)", "Elk"))


## collapse to detections

df$image_trigger_mode%>%unique()
detections3min <-wt_ind_detect_custom(df%>%filter(image_trigger_mode!="Time Lapse"),
                                  threshold = 3,
                                  units = "minutes",
                                  remove_human = FALSE, 
                                  remove_domestic = FALSE
)


detections <-wt_ind_detect_custom(df%>%filter(image_trigger_mode!="Time Lapse"),
                            threshold = 10,
                            units = "minutes",
                            remove_human = FALSE, 
                            remove_domestic = FALSE
)

##prep effort data using timelapse images
eff.dat <- df%>%
  group_by(project_id, location,image_date_time) %>%
  mutate(day = date(image_date_time),
         year = as.integer(format(image_date_time, "%Y"))) %>%
  ungroup%>%
  select(project_id, location, year, day)%>%
  distinct(location, day, .keep_all=TRUE)


cam.summary <- 
  wt_summarise_cam_custom(
    detect_data = detections,
    effort_data = eff.dat,
    time_interval = "month",
    variable = "all",
    output_format = "long",
    exclude_out_of_range = TRUE
  ) %>%
  pivot_wider(names_from = "variable", values_from = "value") %>%
  left_join(
    location %>%
      select(
        location, latitude, longitude, location_comments
      ),
    by = c("location")
  ) %>%
  left_join(sp.lkup, by = "species_common_name")%>%
  mutate(CameraModel="HP2X")%>%
  ungroup()

cam.summary.week <- 
  wt_summarise_cam_custom(
    detect_data = detections,
    effort_data = eff.dat,
    time_interval = "week",
    variable = "all",
    output_format = "long",
    exclude_out_of_range = TRUE
  ) %>%
  pivot_wider(names_from = "variable", values_from = "value") %>%
  left_join(
    location %>%
      select(
        location, latitude, longitude, location_comments
      ),
    by = c("location")
  ) %>%
  left_join(sp.lkup, by = "species_common_name")%>%
  mutate(CameraModel="HP2X")%>%
  ungroup()

## species to drop
cam.summary %>%
  filter(is.na(group1)) %>%
  pull(species_common_name) %>%
  unique()

# drop and count freq of detections
cam.summary <- cam.summary %>%
  drop_na(group1) %>%
  group_by(species_common_name) %>%
  mutate(prevalence = sum(presence))%>%
  ungroup()

cam.summary.week <- cam.summary.week %>%
  drop_na(group1) %>%
  group_by(species_common_name) %>%
  mutate(prevalence = sum(presence))%>%
  ungroup()


hist(cam.summary$prevalence)


##export
write_csv(location, here::here("output/data/mtfernie_cams", "cam.locs.csv"))
write_csv(df, here::here("output/data/mtfernie_cams", "raw.data.csv"))
write_csv(detections, here::here("output/data/mtfernie_cams", "detections_10min.csv"))
write_csv(detections3min, here::here("output/data/mtfernie_cams", "detections_3min.csv"))
write_csv(cam.summary, here::here("output/data/mtfernie_cams", "cam.summary.monthly.csv"))
write_csv(cam.summary.week, here::here("output/data/mtfernie_cams", "cam.summary.week.csv"))