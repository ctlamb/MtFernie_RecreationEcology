library(readxl)
#library(wildrtrax)
library(sf)
library(tidyverse)
source("helpers/ind_det.R")
source("helpers/summarise_detections.R")
# load species lookup
sp.lkup <- read_csv(here::here("data/species_lookup.csv"))
sr16 <- read_excel("/Users/claytonlamb/Dropbox/Documents/University/Work/RecreationMonitoring/MtFernie/data/other_projects/FH_EChow/2016_Classified Image Data_Kootenays FINAL.xlsx")
sr17 <- read_excel("/Users/claytonlamb/Dropbox/Documents/University/Work/RecreationMonitoring/MtFernie/data/other_projects/FH_EChow/2017_Classified Image Data_Kootenays FINAL.xlsx")
sr18 <- read_excel("/Users/claytonlamb/Dropbox/Documents/University/Work/RecreationMonitoring/MtFernie/data/other_projects/FH_EChow/2018-19_Classified Image Data_Kootenays IN PROGRESS.xlsx")
sr.deployments <- read_excel("/Users/claytonlamb/Dropbox/Documents/University/Work/RecreationMonitoring/MtFernie/data/other_projects/FH_EChow/2021_Deployment Records_Kootenays.xlsx")
sr.trailtype <- read_excel("/Users/claytonlamb/Dropbox/Documents/University/Work/RecreationMonitoring/MtFernie/data/other_projects/FH_EChow/Kootenay Remote Camera Project_Camera Info_20Nov2019.xlsx")


sr <- rbind(sr16,
            sr17,
            sr18)

##remove WK cams
sr <- sr%>%
  filter(!str_detect(Location, "NCC WL0"))
sr.deployments <- sr.deployments%>%
  filter(!str_detect(Location, "NCC WL0"))

##get camera names sorted
clean_string <- function(x) {
  gsub("[^A-Za-z0-9]", "", x)
}

sr.deployments <- sr.deployments%>%
  mutate(location=clean_string(Location))%>%
  left_join(sr.trailtype%>%mutate(location=clean_string(Location))%>%select(location,TrailType,RubTree_Lick), 
            by="location")

sr <- sr%>%
  mutate(location=clean_string(Location))



# camera_days <- sr.deployments %>%
#   filter(CameraWorking == "Yes") %>%
#   mutate(days = as.numeric(difftime(DateLastWorking, DateStart, units = "days"))) %>%
#   group_by(location) %>%
#   summarize(total_operating_days = sum(days, na.rm = TRUE))%>%
#   left_join(sr.deployments%>%distinct(location, CameraModel)%>%group_by(location)%>%slice(1),
#             by="location")%>%
#   mutate(CameraModel = ifelse(CameraModel == "Reconyx_Unknown", "Hyperfire PC900", CameraModel))

##align species names
sr$Species%>%unique()
sr<- sr%>%
  mutate(Species = case_when(Species == "Wolf"~"Gray Wolf",
                          Species=="Biker"~"Bike",
                          Species=="Hiker"~"Human",
                          Species=="Runner"~"Human",
                          Species=="Snowshoer"~"Human",
                          Species=="Quad/Vehicle/Motorbike"~"Vehicle",
                          Species=="Snowmobile"~"All Terrain Vehicle",
                          Species=="Dog"~"Domestic Dog",
                          TRUE~Species))


library(fuzzyjoin)
library(lubridate)

# Ensure DateImage is Date class
sr <- sr %>%
  mutate(DateImage = as.Date(DateImage),
    image_date_time = ymd_hms(
      paste(
        DateImage,               # Get just the date part
        format(TimeImage, "%H:%M:%S")     # Extract just the time
      )
    )
  )

sr.deployments.keep <- sr.deployments %>%
  filter(CameraWorking == "Yes",
         location%in%sr$location) %>%
  mutate(project_id="Emcams",
         start_time = as.Date(DateStart),
         end_time = as.Date(DateLastWorking))


##extract utms
sr.locs <- sr.deployments %>%
  distinct(location, .keep_all=TRUE)%>%
  mutate(Easting = as.numeric(str_match(Location, "([0-9]+)\\s+([0-9]+)$")[,2]),
         Northing = as.numeric(str_match(Location, "([0-9]+)\\s+([0-9]+)$")[,3]))%>%
  select(-Location)%>%
  st_as_sf(coords=c("Easting", "Northing"), crs=26911)%>%
  st_transform(4326)%>%
  cbind(st_coordinates(.))
  

# Interval join: keep only images within deployment windows
sr.filtered <- fuzzy_left_join(
  sr,
  sr.deployments.keep %>% select(location, CameraModel, start_time, end_time, TrailType, RubTree_Lick),
  by = c(
    "location" = "location",
    "DateImage" = "start_time",
    "DateImage" = "end_time"
  ),
  match_fun = list(`==`, `>`, `<=`)
) %>%
  filter(!is.na(end_time))  # Only keep matched rows


# match names to wildrtrax data
sr.filtered.wt <- sr.filtered%>%
  filter(!RubTree_Lick%in%c("Rub Tree","Mineral Lick"),
         TrailType%in%c("Wildlife"))%>%
  mutate(project_id="Emcams",
         species_common_name=Species)%>%
  filter(Species!="LastImage")%>%
  select(project_id, location=location.x, species_common_name, start_time=image_date_time, max_animals=Total, CameraModel)

##custom effort data
sr.deployments.keep <- sr.deployments.keep%>%
  select(project_id, location, start_time,end_time, CameraModel)

sr.deployments.effort <- sr.deployments.keep%>%
  group_by(location) %>%
  mutate(day = map2(start_time, end_time, ~ seq.Date(.x, .y, by = "day"))) %>%
  unnest(day) %>%
  mutate(year = as.integer(format(day, "%Y"))) %>%
  ungroup() %>%
  select(location, year, day) %>%
  distinct(location, day, .keep_all = TRUE)

##get a summary
cam.summary <- 
  wt_summarise_cam_custom(
    detect_data = sr.filtered.wt,
    effort_data = sr.deployments.effort,
    time_interval = "month",
    variable = "all",
    output_format = "long"
  )%>%
  pivot_wider(names_from = "variable", values_from = "value")%>%
  left_join(sp.lkup, by = "species_common_name")%>%
  left_join(sr.deployments.keep%>%distinct(location, CameraModel), by="location")%>%
  ungroup()


cam.summary.week <- 
  wt_summarise_cam_custom(
    detect_data = sr.filtered.wt,
    effort_data =sr.deployments.effort,
    time_interval = "week",
    variable = "all",
    output_format = "long"
  )%>%
  pivot_wider(names_from = "variable", values_from = "value")%>%
  left_join(sp.lkup, by = "species_common_name")%>%
  left_join(sr.deployments.keep%>%distinct(location, CameraModel), by="location")%>%
  ungroup()


##export
write_csv(sr.locs, here::here("output/data/sr_cams", "cam.locs.csv"))
write_csv(sr.filtered.wt, here::here("output/data/sr_cams", "detections.csv"))
write_csv(cam.summary, here::here("output/data/sr_cams", "cam.summary.monthly.csv"))
write_csv(cam.summary.week, here::here("output/data/sr_cams", "cam.summary.week.csv"))


