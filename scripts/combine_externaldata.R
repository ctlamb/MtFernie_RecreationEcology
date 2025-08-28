library(tidyverse)
library(sf)

source(here::here("helpers", "gg_theme_custom.R"))
##LOAD DATA
##BEC ZONES
bec <- st_read("/Users/claytonlamb/Dropbox/Documents/University/Geographic_Data/BEC/BEC_BIOGEOCLIMATIC_POLY/BEC_POLY_polygon.shp")

###CAMS
##em cam data
#mapview(read_csv(here::here("output/data/sr_cams", "cam.locs.csv"))%>%st_as_sf(coords=c("X","Y"), crs=4326))
# sr.cams.select <- c("LC16393005481619", "Saskatoon6241335497808",
# "Southgal6276245501552", "BullSulphur6400405507741", "CDT6480825494218", "WestHosmer6461175496498",
# "NCC6425065476463", "13KmRiverRoad6448405467741", "LowerCoal6491055485051", "Morrisey6539105476756")
sr.cam.locs <- read_csv(here::here("output/data/sr_cams", "cam.locs.csv"))
sr.detections <- read_csv(here::here("output/data/sr_cams", "detections.csv"))%>%
  mutate(month=month(start_time, label=TRUE, abbr=FALSE),
         project="BCgov")
sr.cam.summary <- read_csv(here::here("output/data/sr_cams", "cam.summary.monthly.csv"))%>%
  mutate(project="BCgov")


##hwy 3
hwy3.cam.summary <- read_csv(here::here("output/data/hwy3", "monthly_summary10min.csv"))%>%
  #filter(location%in%c("DUMP-W1", "OLSENO-C2","DUMP-E2","DUMP-W1B","DUMP-W2"))%>%
  filter(!str_detect(location, "LIZARD|LOOP-C2"))%>% ##both em cams and should be in her data
  filter(type=="Control")%>%
  rename(year=project)%>%
  mutate(project="rtrbc")
hwy3.cam.keep <- hwy3.cam.summary$location%>%unique()
hwy3.detections <- read_csv(here::here("output/data/hwy3", "detections10min.csv"))%>%
  filter(location%in%hwy3.cam.keep)%>%
  mutate(month=month(start_time, label=TRUE, abbr=FALSE),
         project="rtrbc")
hwy3.cam.locs <- read_csv("/Users/claytonlamb/Dropbox/Documents/University/Work/BCGov/RTR_BC/output/F25_RTR/data/clean/cams_clean.csv")%>%
  select(location=camera_name,Y=location_latitude, X=location_longitude)%>%
  mutate(location=str_to_upper(location))%>%
  filter(location%in%hwy3.cam.keep)%>%
  distinct()

##Heikos
heikos.detections <- read_csv(here::here("output/data/heikos", "detections10min.csv"))%>%
  filter(type2=="Wildlife-far")%>%
  mutate(month=month(start_time, label=TRUE, abbr=FALSE),
         project="heikos")
heikos <- read_csv(here::here("output/data/heikos", "monthly_summary10min.csv"))%>%
  filter(type2=="Wildlife-far")%>%
  mutate(project="heikos")
  

##hosmer
hoz.detections <- read_csv(here::here("output/data/transrockies", "detections10min.csv"))%>%
  filter(location%in%c("COAL4", "COAL2"))%>%
  mutate(month=month(start_time, label=TRUE, abbr=FALSE),
         project="transrockies")
  
hoz <- read_csv(here::here("output/data/transrockies", "monthly_summary10min.csv"))%>%
  filter(location%in%c("COAL4", "COAL2"))%>%
  mutate(project="transrockies")
  
hoz.cam.locs <- tibble(location=c("COAL4", "COAL2"),
                       X=c(-114.944031, -114.934000),
                       Y=c(49.600758,  49.610346))



##mtn pass
mtnpass.detections <- read_csv(here::here("output/data/mtnpass", "detections10min.csv"))%>%
  mutate(month=month(start_time, label=TRUE, abbr=FALSE),
         project="mtnpass")

mtnpass <- read_csv(here::here("output/data/mtnpass", "monthly_summary10min.csv"))%>%
  mutate(project="mtnpass")%>%
  pivot_wider(names_from = "variable", values_from = "value")



##combine other cams
external.cam.detections <- sr.detections%>%
  left_join(sr.cam.locs%>%select(location:Y)%>%distinct(), by="location")%>%
  mutate(location_comments="Wildlife-sr")%>%
  select(project, location, CameraModel, location_comments, species_common_name, start_time, max_animals, latitude=Y, longitude=X)%>%
  rbind(hwy3.detections%>%
          left_join(hwy3.cam.locs%>%distinct(), by="location")%>%
          mutate(location_comments="Wildlife-hwy3",
                 CameraModel="HP2X")%>%
          select(project, location, CameraModel, location_comments, species_common_name, start_time, max_animals, latitude=Y, longitude=X))%>%
  rbind(heikos.detections%>%
          mutate(location_comments="Wildlife-Heikos",
                 CameraModel="HP2X")%>%
          select(project, location, CameraModel, location_comments, species_common_name, start_time, max_animals, latitude=Y, longitude=X))%>%
  rbind(hoz.detections%>%
          left_join(hoz.cam.locs%>%distinct(), by="location")%>%
          mutate(location_comments="Wildlife-Hoz",
                 CameraModel="HP2X")%>%
          select(project, location, CameraModel, location_comments, species_common_name, start_time, max_animals, latitude=Y, longitude=X))%>%
  rbind(mtnpass.detections%>%
          mutate(location_comments="Wildlife-Mtnpass",
                 CameraModel="HP2X")%>%
          select(project, location, CameraModel, location_comments, species_common_name, start_time, max_animals, latitude=Y, longitude=X))%>%
  mutate(month=month(start_time, label=TRUE, abbr=FALSE),
         species_common_name=case_when(species_common_name=="Elk (wapiti)"~"Elk", TRUE~species_common_name))

external.cam.summary <- sr.cam.summary%>%
  left_join(sr.cam.locs%>%select(location:Y)%>%distinct(), by="location")%>%
  mutate(location_comments="Wildlife-sr")%>%
  select(project, location, CameraModel, location_comments, species_common_name,counts, detections,  month, year, wild, n_days_effort, latitude=Y, longitude=X)%>%
  rbind(hwy3.cam.summary%>%
          left_join(hwy3.cam.locs%>%distinct(), by="location")%>%
          mutate(location_comments="Wildlife-hwy3",
                 CameraModel="HP2X")%>%
          select(project, location, CameraModel, location_comments, species_common_name,counts, detections, month, year, wild, n_days_effort, latitude=Y, longitude=X))%>%
  rbind(heikos%>%
          mutate(location_comments="Wildlife-Heikos",
                 CameraModel="HP2X")%>%
          select(project, location, CameraModel, location_comments, species_common_name,counts, detections, month, year, wild, n_days_effort, latitude=Y, longitude=X))%>%
  rbind(hoz%>%
          mutate(location_comments="Wildlife-Hoz",
                 CameraModel="HP2X")%>%
          select(project, location, CameraModel, location_comments, species_common_name,counts, detections, month, year, wild, n_days_effort, latitude, longitude))%>%
  rbind(mtnpass%>%
          left_join(mtnpass.detections%>%distinct(location,X,Y), by="location")%>%
          mutate(location_comments="Wildlife-Mtnpass",
                 CameraModel="HP2X")%>%
          select(project, location, CameraModel, location_comments, species_common_name,counts, detections, month, year, wild, n_days_effort, latitude=Y, longitude=X))%>%
  mutate(species_common_name=case_when(species_common_name=="Elk (wapiti)"~"Elk", TRUE~species_common_name))



##subset to same BEC Zone as Mt Fernie cams "ICH" and within 25k of Fernie

# Define Fernie coordinates (in lon/lat)
fernie_coords <- st_sfc(st_point(c(-115.058, 49.504)), crs = 4326) %>%
  st_transform(st_crs(bec))

# Create a 25 km buffer around Fernie
fernie_buffer <- st_buffer(fernie_coords, dist = 25000)  # 25 km in meters

# Filter detections to ICH/MS zones and within 25 km of Fernie
external.cam.detections.bec <- external.cam.detections %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(st_crs(bec)) %>%
  st_join(bec %>% select(ZONE), left = TRUE) %>%
  filter(ZONE %in% c("ICH", "MS")) %>%
  filter(st_intersects(geometry, fernie_buffer, sparse = FALSE))

external.cam.summary.bec <- external.cam.summary %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(st_crs(bec)) %>%
  st_join(bec %>% select(ZONE), left = TRUE) %>%
  filter(ZONE %in% c("ICH", "MS")) %>%
  filter(st_intersects(geometry, fernie_buffer, sparse = FALSE))

external.cam.detections$project%>%unique()
external.cam.detections.bec.inverse <- external.cam.detections %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(st_crs(bec)) %>%
  st_join(bec %>% select(ZONE), left = TRUE) %>%
  filter(project%in%c("BCgov", "heikos","mtnpass"))%>% ##no transrockies or rtrbc as these are valley bottom front country
  filter(!ZONE %in% c("ICH", "MS", "IDF"))

human.cams <- external.cam.detections.bec.inverse%>%
  filter(species_common_name%in%c("Human","Vehicle"))%>%
  group_by(location)%>%
  count()%>%
  filter(n>10)%>%
  pull(location)

external.cam.detections.bec.inverse <- external.cam.detections.bec.inverse%>%
  filter(!location%in%human.cams)


#check
ggplot(external.cam.summary.bec%>%
         filter(wild=="Wild") %>%
         select(location, location_comments, species_common_name,counts, month, wild,n_days_effort)%>%
         filter(!month%in%c("December", "January", "February", "March"))%>%
         group_by(location, species_common_name, location_comments,wild) %>%
         mutate(hitrate=(counts/n_days_effort)*30)%>%
         summarise(hitrate = mean(hitrate)), aes(x = location, y = hitrate, fill = species_common_name)) +
  geom_col() +
  theme.custom +
  theme(
    plot.title = element_text(size = 20),
    axis.title.y = element_text(size = 17),
    axis.title.x = element_text(size = 17),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    strip.text.x = element_text(size = 12),
    strip.text.y = element_text(size = 18),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 20),
    legend.position = "none"
  ) +
  labs(x = "Location", y = "Counts per month", title = "", colour = "Type", fill = "Type")+
  facet_wrap(vars(location_comments), scales="free_x")


#remove some cameras with issues or those that aren't legit comparisons
external.cam.summary.bec <- 
  external.cam.summary.bec%>%
  filter(!str_detect(location,"TUNNEL"),##in rock slope for sheep
         !str_detect(location,"BuffaloHead|Hartley|SulphurLick"))%>%
  st_transform(4326)%>%
  cbind(st_coordinates(.))

external.cam.detections.bec <- external.cam.detections.bec%>%
  filter(!str_detect(location,"TUNNEL"),##in rock slope for sheep
         !str_detect(location,"BuffaloHead|Hartley|SulphurLick")) %>%
  st_transform(4326)%>%
  cbind(st_coordinates(.))

external.cam.detections.bec.inverse <- external.cam.detections.bec.inverse%>%
  filter(!str_detect(location,"TUNNEL"),##in rock slope for sheep
         !str_detect(location,"BuffaloHead|Hartley|SulphurLick|NFORK|RACEHORSE2|RACEHORSE5|RACEHORSE6|DEADMANTRAIL3|DEADMANTRAIL5|RamCabin|McLatchie|Foisey|Morrisey|Morel|NorthFork")) %>%
  st_transform(4326)%>%
  cbind(st_coordinates(.))



##write
write_csv(external.cam.detections.bec, "data/other_projects/external.detections.csv")
write_csv(external.cam.detections.bec.inverse, "data/other_projects/external.detections.inverse.csv")
write_csv(external.cam.summary.bec, "data/other_projects/external.summary.csv") 


# a <- external.cam.detections.bec.inverse%>%
#   distinct(X,Y,location)%>%
#   st_as_sf(coords=c("X","Y"), crs=4326)
# 
# mapview(a)
