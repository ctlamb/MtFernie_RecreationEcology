# library(renv)
# library(remotes)
# remotes::install_github("ABbiodiversity/wildrtrax@development")
library(wildrtrax)
library(sf)
library(mapview)
library(patchwork)
library(basemaps)
library(raster)
library(terra)
library(here)
library(activity)
library(lubridate)
library(readxl)
library(overlap)
library(glmmTMB)
library(furrr)
library(MetBrewer)
library(tidyverse)
library(tidylog)
options(scipen = 999)

source(here::here("helpers", "gg_theme_custom.R"))



## check times for MST vs MDT. Using daylight in images correlated with sunset:
## Heikos:MDT
## Fernie Rec: MDT
## Hwy 3: seems to be on MDT
## Mtn pass on MDT

# Load data ---------------------------------------------------------------

# load species lookup
sp.lkup <- read_csv(here::here("data/species_lookup.csv"))


### set months to keep
months.keep <- c("May", "June", "July", "August", "September", "October")

# Create a sequence of all dates in a year
all.dates <- tibble(date = seq(ymd("2023-01-01"), ymd("2023-12-31"), by = "1 day"))

# Extract weeks for your selected months
weeks.keep <- all.dates %>%
  mutate(
    month_name = month(date, label = TRUE, abbr = FALSE),
    week_num = isoweek(date)
  ) %>%
  filter(month_name %in% months.keep) %>%
  distinct(week_num) %>%
  pull(week_num)

## fernie rec data
cam.locs <- read_csv(here::here("output/data/mtfernie_cams", "cam.locs.csv")) %>%
  filter(!location %in% c("ROAD3", "ROAD4"))
df <- read_csv(unz(here::here("output/data/mtfernie_cams", "raw.data.csv.zip"), "raw.data.csv")) %>%
  mutate(month = month(image_date_time, label = TRUE, abbr = FALSE)) %>%
  filter(month %in% months.keep) %>%
  mutate(species_common_name = str_to_lower(species_common_name)) %>%
  filter(
    !location %in% c("ROAD3", "ROAD4"),
    species_common_name != "domestic dog"
  )
detections <- read_csv(here::here("output/data/mtfernie_cams", "detections_10min.csv")) %>%
  mutate(
    month = month(start_time, label = TRUE, abbr = FALSE),
    species_common_name = str_to_lower(species_common_name)
  ) %>%
  filter(month %in% months.keep) %>%
  filter(
    !location %in% c("ROAD3", "ROAD4"),
    species_common_name != "domestic dog"
  )
detections3min <- read_csv("output/data/mtfernie_cams/detections_3min.csv") %>%
  mutate(species_common_name = str_to_lower(species_common_name)) %>%
  filter(
    !location %in% c("ROAD3", "ROAD4"),
    species_common_name != "domestic dog"
  )
cam.summary.week <- read_csv(here::here("output/data/mtfernie_cams", "cam.summary.week.csv")) %>%
  filter(week %in% weeks.keep) %>%
  mutate(project = "mt_fernie") %>%
  mutate(species_common_name = str_to_lower(species_common_name)) %>%
  filter(
    !location %in% c("ROAD3", "ROAD4"),
    species_common_name != "domestic dog"
  )
cam.summary <- read_csv(here::here("output/data/mtfernie_cams", "cam.summary.monthly.csv")) %>%
  filter(month %in% months.keep) %>%
  mutate(project = "mt_fernie") %>%
  mutate(species_common_name = str_to_lower(species_common_name)) %>%
  filter(
    !location %in% c("ROAD3", "ROAD4"),
    species_common_name != "domestic dog"
  )

## external cam data
external.detections <- read_csv(here::here("data/other_projects/external.detections.csv")) %>%
  filter(month %in% months.keep) %>%
  mutate(species_common_name = str_to_lower(species_common_name))
external.detections.inverse <- read_csv(here::here("data/other_projects/external.detections.inverse.csv")) %>%
  filter(month %in% months.keep) %>%
  mutate(species_common_name = str_to_lower(species_common_name))
external.cam.summary <- read_csv(here::here("data/other_projects/external.summary.csv")) %>%
  filter(month %in% months.keep) %>%
  mutate(species_common_name = str_to_lower(species_common_name))


##cam locs
all.cam.locs <- external.detections.inverse%>%
  select(project, location, CameraModel, location_comments, X, Y)%>%
  mutate(location_comments="Landscape-remote")%>%
  rbind(external.detections%>%
          select(project, location, CameraModel, location_comments, X, Y)%>%
          mutate(location_comments="Landscape-valley"))%>%
  rbind(cam.locs%>%
          mutate(project="MtFernie",
                 CameraModel="HP2X",
                 location_comments)%>%
          select(project, location, CameraModel, location_comments, X=longitude, Y=latitude))%>%
  ungroup%>%
  distinct(location, .keep_all=TRUE)%>%
  st_as_sf(coords=c("X", "Y"), crs=4326)

mapview(all.cam.locs, zcol="project")
mapview(all.cam.locs, zcol="location_comments")
st_write(all.cam.locs, "output/spatial/cams/analysis_cams.shp", delete_dsn = TRUE)


### summary for paper results
external.cam.summary %>%
  group_by(wild) %>%
  summarise(detections = sum(detections))

external.cam.summary %>%
  distinct(location, month, year, n_days_effort) %>%
  summarise(monitoringdays = sum(n_days_effort))



### summary for paper results
range(df$image_date_time)
nrow(df)
nrow(detections)
cam.summary %>%
  group_by(wild) %>%
  summarise(detections = sum(detections))

cam.summary %>%
  distinct(location, month, year, n_days_effort) %>%
  summarise(monitoringdays = sum(n_days_effort))


cam.table <- cam.summary %>%
  filter(!species_common_name %in% c("Bird", "Wolves, Coyotes and Allies")) %>%
  group_by(wild, species_common_name) %>%
  summarise(
    sites = n_distinct(location[presence == 1]),
    detections = sum(detections),
    trap.nights=sum(n_days_effort)
  ) %>%
  arrange(wild, -detections)

cam.hits <- cam.summary %>%
  dplyr::group_by(species_common_name) %>%
  dplyr::summarise(
    animals = sum(counts),
    days = sum(n_days_effort)
  )%>%
  ungroup()%>%
  mutate(hitrate_per100 = round((animals / days) * 100, 1))
  


cam.table.join <- cam.table %>%
  select(-trap.nights)%>%
  left_join(cam.hits%>%select(species_common_name, hitrate_per100), by = "species_common_name") %>%
  filter(
    hitrate_per100 > 0,
    species_common_name != "domestic dog"
  ) %>%
  arrange(wild, -hitrate_per100)

write_csv(cam.table.join, "output/tables/detection_summary.csv")


####################
### EXPLORATORY PLOTS
####################


cam.summary %>%
  filter(group2 %in% c("bike", "carnivore", "Human", "ungulate", "motorized vehicle", "domestic animal")) %>%
  group_by(location, species_common_name, location_comments, group2) %>%
  summarise(hitrate = mean(counts)) %>%
  ggplot(aes(x = location, y = hitrate, color = group2, fill = group2)) +
  geom_col(position = position_dodge(width = 0.5)) +
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
    legend.position = "bottom"
  ) +
  labs(x = "Location", y = "Counts per month", title = "Human vs Wildlife Use", colour = "Type", fill = "Type")


cam.summary %>%
  mutate(month = factor(month, levels = month.name)) %>%
  group_by(wild, month) %>%
  summarise(hitrate = mean(counts), .groups = "drop") %>%
  ggplot(aes(x = month, y = hitrate)) +
  geom_col() +
  facet_wrap(vars(wild), scales = "free")

### !!!! Different (smaller) values than used in spring 2024 ppt. check #s
cam.summary %>% distinct(species_common_name, prevalence)
sp.keep <- cam.summary %>%
  filter(prevalence > 30 & group1 != "bird", !species_common_name %in% c("Deer")) %>%
  dplyr::select(location, location_comments, species_common_name, counts) %>%
  distinct(species_common_name) %>%
  pull(species_common_name)






cam.plot.raw <- cam.summary %>%
  select(project, location, CameraModel, location_comments, species_common_name, counts, month, wild, n_days_effort) %>%
  rbind(external.cam.summary %>%
    # filter(location_comments!="Wildlife-sr")%>%
    select(project, location, CameraModel, location_comments, species_common_name, counts, month, wild, n_days_effort)) %>%
  filter(species_common_name %in% sp.keep) %>%
  group_by(location, CameraModel, species_common_name, location_comments, wild) %>%
  mutate(hitrate = (counts / n_days_effort) * 30) %>%
  summarise(hitrate = mean(hitrate))%>%
  ungroup()

# Reusable base plot setup
# Define shared facet layout

cam.plot.raw$location_comments %>% unique()
cam.plot <- cam.plot.raw %>%
  #filter(!location_comments %in% c("Road")) %>%
  mutate(location_comments = case_when(
    location_comments %in% c("Wildlife-hwy3", "Wildlife-Heikos", "Wildlife-Hoz", "Wildlife-sr") ~ "Landscape-valley",
    TRUE ~ location_comments
  )) %>%
  ungroup() %>%
  mutate(location_comments = fct_relevel(location_comments, "Landscape-valley", "Wildlife trail", "Rec trail"))
cam.plot$location_comments %>% unique()

unique(cam.plot$location)

cam.plot %>%
  group_by(location_comments) %>%
  summarise(n = n_distinct(location))

#### Stats
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(rstatix)
library(stringr)

# Prep
location_levels <- c("Landscape-valley", "Wildlife trail", "Rec trail")
human_species <- c("Bike", "Human") %>% str_to_lower()
wild_species.all <- c("Black Bear", "Elk", "Grizzly Bear", "Moose", "Mule Deer", "Red Fox", "White-tailed Deer", "Gray Wolf", "Cougar") %>% str_to_lower()
wild_species <- c("Black Bear", "Elk", "Grizzly Bear", "Moose", "Mule Deer", "Red Fox", "White-tailed Deer") %>% str_to_lower()


plot.dat <- cam.plot %>%
  mutate(species_common_name = case_when(
    species_common_name %in% c("all terrain Vehicle", "vehicle") ~ "motorized vehicle",
    TRUE ~ species_common_name
  )) %>%
  mutate(location_comments = factor(location_comments, levels = location_levels)) %>%
  group_by(species_common_name) %>%
  mutate(hitrate = if_else(row_number() == 1, hitrate + 0.0001, hitrate)) %>%
  ungroup()

# Kruskal + pairwise Wilcoxon
pairwise_results <- plot.dat %>%
  filter(species_common_name %in% c(human_species, wild_species.all)) %>%
  mutate(location_comments2 = ifelse(location_comments == "Landscape-valley", "Landscape-valley", "Mt Ferie Rec")) %>%
  group_by(species_common_name) %>%
  pairwise_wilcox_test(hitrate ~ location_comments2, p.adjust.method = "holm")

pairwise_results <- plot.dat %>%
  filter(species_common_name %in% c(human_species, wild_species.all)) %>%
  group_by(species_common_name) %>%
  pairwise_wilcox_test(hitrate ~ location_comments, p.adjust.method = "holm")
pairwise_results %>%
  filter(group1 == "Wildlife trail")
# y-position reference
y_pos <- plot.dat %>%
  group_by(species_common_name, location_comments) %>%
  summarise(y = mean(hitrate) + (sd(hitrate) / sqrt(n())), .groups = "drop") %>%
  group_by(species_common_name) %>%
  summarise(y = max(y))

# Add significance labels
plot_labels <- pairwise_results %>%
  mutate(
    label = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01 ~ "**",
      p.adj <= 0.05 ~ "*",
      p.adj > 0.05 ~ ""
    )
  ) %>%
  # Join for group1
  left_join(y_pos, by = c("species_common_name")) %>%
  group_by(species_common_name) %>%
  mutate(
    y_rank = row_number(),
    y_base = y,
    y_position = y_base * 1.05 + ((y_rank - 1) * 0.1 * y_base) # Adjust spacing here
  ) %>%
  ungroup()

# Final plot
hum_dat <- plot.dat %>%
  filter(species_common_name %in% human_species) %>%
  mutate(species_common_name = factor(species_common_name, levels = c("", " ", human_species)))
# Add two blank rows to force empty facets
blank_rows <- tibble(
  location = NA,
  species_common_name = factor(c("", " "), levels = c("", " ", human_species)),
  location_comments = NA,
  wild = NA,
  hitrate = NA
)
# Combine real and blank data
hum_dat <- bind_rows(hum_dat, blank_rows)

human_plot <- ggplot(hum_dat, aes(x = location_comments, y = hitrate, fill = location_comments)) +
  stat_summary(fun = mean, geom = "col") +
  stat_summary(fun.data = mean_se, geom = "linerange", linewidth = 1) +
  facet_wrap(~ fct_relevel(species_common_name, "bike", "domestic dog", "human", "", "motorized vehicle", " "), ncol = 1, scales = "free") +
  geom_segment(
    data = plot_labels %>% filter(species_common_name %in% human_species),
    aes(x = group1, xend = group2, y = y_position, yend = y_position),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = plot_labels %>% filter(species_common_name %in% human_species),
    aes(
      x = as.numeric(factor(group1, levels = location_levels)) * 0.5 +
        as.numeric(factor(group2, levels = location_levels)) * 0.5,
      y = y_position * 1.01,
      label = label
    ),
    inherit.aes = FALSE,
    size = 4
  ) +
  scale_fill_manual(values = c(
    "Rec trail" = "#E76F51",
    "Wildlife trail" = "#2A9D8F",
    "Landscape-valley" = "grey60"
  )) +
  labs(x = "Location", y = "Counts per month", title = "Human use") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Subsets
plot.dat %>% filter(species_common_name %in% wild_species.all)%>%group_by(species_common_name, location_comments)%>%summarise(mean=mean(hitrate))

wild_plot <- ggplot(plot.dat %>% filter(species_common_name %in% wild_species.all), aes(x = location_comments, y = hitrate, fill = location_comments)) +
  stat_summary(fun = mean, geom = "col") +
  stat_summary(fun.data = mean_se, geom = "linerange", linewidth = 1) +
  facet_wrap(~species_common_name, scales = "free", ncol = 3) +
  geom_segment(
    data = plot_labels %>% filter(species_common_name %in% wild_species.all),
    aes(x = group1, xend = group2, y = y_position, yend = y_position),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = plot_labels %>% filter(species_common_name %in% wild_species.all),
    aes(
      x = as.numeric(factor(group1, levels = location_levels)) * 0.5 +
        as.numeric(factor(group2, levels = location_levels)) * 0.5,
      y = y_position * 1.01,
      label = label
    ),
    inherit.aes = FALSE,
    size = 4
  ) +
  scale_fill_manual(values = c(
    "Rec trail" = "#E76F51",
    "Wildlife trail" = "#2A9D8F",
    "Landscape-valley" = "grey60"
  )) +
  labs(x = "Location", y = "", title = "Wildlife use") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )



# Combine side-by-side using patchwork
human_plot + wild_plot + plot_layout(widths = c(1, 3))

## save
ggsave("output/plots/hit_rate_compare.png", dpi = 300, width = 10, height = 9, unit = "in")


### test significance
lm_results <- cam.plot %>%
  filter(location_comments %in% c("Landscape-valley", "Rec trail", "Wildlife trail")) %>%
  group_by(species_common_name) %>%
  group_modify(~ {
    model <- lm(hitrate ~ location_comments, data = .x)
    broom::tidy(model)
  }) %>%
  ungroup() %>%
  filter(term != "(Intercept)")


lm_results2 <- cam.plot %>%
  filter(location_comments %in% c("Landscape-valley", "Rec trail", "Wildlife trail")) %>%
  mutate(location_comments = ifelse(location_comments == "Landscape-valley", "Landscape-valley", "Mt_Fernie")) %>%
  group_by(species_common_name) %>%
  group_modify(~ {
    model <- lm(hitrate ~ location_comments, data = .x)
    broom::tidy(model)
  }) %>%
  ungroup() %>%
  filter(term != "(Intercept)")




# sp.loop <- (weekly.mod.wild$species_common_name%>%unique())[-c(2:3,5,10)]
# loop.dat <- tibble()
# # + `Domestic Dog_scale`
# for(i in 1:length(sp.loop)){
# library(broom.mixed)
#   loop.dat <- glmmTMB(
#   counts ~ Bike_scale +  Human_scale  + presence.lag + elev_scale +canopy.all_scale + slope_deg_250_scale +  bldg.dist_scale + trail.width +
#     (1 | week) + (1 | year) + (1 | location),                 # Count model (for counts > 0)
#   ziformula = ~ Bike_scale +  Human_scale + presence.lag + elev_scale +canopy.all_scale + slope_deg_250_scale +  bldg.dist_scale + trail.width,  # Zero-hurdle model
#   data = weekly.mod.wild %>%
#     #filter(species_common_name == sp.loop[i], !location %in% c("ROAD3", "ROAD4")) %>%
#   filter(species_common_name == sp.loop[i])%>%mutate(trail.width=case_when(location_comments=="Road"~2, TRUE~1))%>%
#   #filter(species_common_name == sp.loop[i]),
#   family = truncated_poisson()
# )%>%
#   tidy%>%
#   mutate(sp=sp.loop[i])%>%
#     rbind(loop.dat)
#
# }
mod.dat <- cam.summary %>%
  select(project, location, CameraModel, location_comments, species_common_name, counts, month, wild, n_days_effort, year) %>%
  rbind(external.cam.summary %>%
    # filter(location_comments!="Wildlife-sr")%>%
    select(project, location, CameraModel, location_comments, species_common_name, counts, month, wild, n_days_effort, year)) %>%
  filter(species_common_name %in% c(human_species, wild_species.all)) %>%
  # mutate(species_common_name = case_when(
  #   species_common_name %in% c("all terrain Vehicle", "vehicle") ~ "motorized vehicle",
  #   TRUE ~ species_common_name
  # )) %>%
  mutate(hitrate = (counts / n_days_effort) * 30) %>%
  mutate(location_comments = case_when(
    location_comments %in% c("Wildlife-hwy3", "Wildlife-Heikos", "Wildlife-Hoz", "Wildlife-sr") ~ "Landscape-valley",
    TRUE~location_comments
  )) %>%
  filter(location_comments %in% c("Landscape-valley", "Rec trail", "Wildlife trail"))
mod.dat %>%
  group_by(location_comments) %>%
  count()
mod.dat %>%
  group_by(location_comments) %>%
  distinct(location) %>%
  count()


mod.dat$location_comments %>% unique()



## add in spatial data

# Load and stack rasters
# ppt <- rast("data/spatial/raw/MAP.tif")
# evi <- rast("data/spatial/raw/evi_summer_2022.tif")
# 
# ppt <- project(ppt, ppt)
# evi <-  project(evi, ppt)
# writeRaster(ppt, "data/spatial/stacked_larger/ppt.tif")
# writeRaster(evi, "data/spatial/stacked_larger/evi.tif")
ppt <- rast("data/spatial/stacked_larger/ppt.tif")
evi <- rast("data/spatial/stacked_larger/evi.tif")

r_stack <- c(ppt, evi)

# Convert sf to SpatVector
points_vect <- all.cam.locs %>%
  st_transform(st_crs(r_stack)) %>%
  vect()

# Extract values
vals <- terra::extract(r_stack, points_vect)

# Combine with original point IDs
cam.locs.attributed <- bind_cols(all.cam.locs, vals)

mod.dat.spat <- mod.dat%>%
  left_join(cam.locs.attributed%>%select(location, evi=evi_summer_2022, MAP), by="location")%>%
  tibble()

## standardize

mod.dat.spat <- mod.dat.spat %>%
  mutate(across(evi:MAP, ~ scale(.)[, 1], .names = "{.col}_scale"))


library(broom.mixed)






model_list <- list(
  MAP = counts ~ location_comments + MAP_scale + (1 | CameraModel) + (1 | month) + (1 | year) + (1 | location),
  evi = counts ~ location_comments + evi_scale + (1 | CameraModel) + (1 | month) + (1 | year) + (1 | location),
  Simple = counts ~ location_comments + (1 | CameraModel) + (1 | month) + (1 | year) + (1 | location)
)

zi_list  <- list(
  MAP = ~ location_comments + MAP_scale + (1 | CameraModel) + (1 | month) + (1 | year) + (1 | location),
  evi = ~ location_comments + evi_scale + (1 | CameraModel) + (1 | month) + (1 | year) + (1 | location),
  Simple = ~ location_comments + (1 | CameraModel) + (1 | month) + (1 | year) + (1 | location)
)


# Prep data
filtered_data <- mod.dat.spat %>%
  mutate(location_comments = ifelse(location_comments == "Landscape-valley", "Landscape-valley", "Mt_Fernie")) %>%
  filter(species_common_name %in% str_to_lower(c("Black Bear", "Elk", "Grizzly Bear", "Moose", "Mule Deer", "Red Fox", "White-tailed Deer")))

##check 0's to assess if hurdle model needed
filtered_data%>%
  group_by(species_common_name)%>%
  summarise(
    mean_count = mean(counts),
    median_count = median(counts),
    pct_zeros = mean(counts == 0) * 100
  )

# Expand over species x model
aic_table <- filtered_data %>%
  group_by(species_common_name) %>%
  group_map(~ {
    species <- unique(.y$species_common_name)
    
    map_dfr(names(model_list), function(mod_name) {
      # Fit the model
      fit <- tryCatch({
        glmmTMB(
          formula = model_list[[mod_name]],
          ziformula = zi_list[[mod_name]],
          offset = .x$n_days_effort,
          family = truncated_poisson(),
          data = .x
        )
      }, error = function(e) NULL)
      
      # Extract AIC if successful
      if (!is.null(fit)) {
        tibble(
          species_common_name = species,
          model_name = mod_name,
          AIC = AIC(fit)
        )
      } else {
        tibble(
          species_common_name = species,
          model_name = mod_name,
          AIC = NA_real_
        )
      }
    })
  }) %>%
  bind_rows()

ggplot(aic_table, aes(x = AIC, y = species_common_name, color = model_name)) +
  geom_point(size = 3) +
  labs(title = "Model comparison by AIC", y = "Species", x = "AIC") +
  theme_minimal()+
  facet_wrap(vars(species_common_name), scales="free")+
  scale_x_continuous(breaks = function(x) seq(floor(min(x, na.rm = TRUE)), ceiling(max(x, na.rm = TRUE)), by = 2))+
  theme(panel.grid.minor.x = element_blank())


## predict
# New data frame with the levels of location_comments for prediction
# model <- glmmTMB( counts ~ location_comments + (1|CameraModel) + (1 | month) + (1 | location),
#                   ziformula = ~location_comments,
#                   offset=n_days_effort,
#                   data = mod.dat%>%filter(species_common_name=="white-tailed deer"),
#                   family = truncated_poisson())




location_lvs <- c("Landscape-valley", "Wildlife trail", "Rec trail")


mean(mod.dat.spat$evi_scale)
mean(mod.dat.spat$MAP_scale)
cor(mod.dat.spat$evi, mod.dat.spat$MAP)

##check make sure makes sense with table etc
a <- mod.dat.spat %>%
  filter(n_days_effort > 10)%>%
  group_by(species_common_name, location_comments)%>%
  summarise(count=sum(counts),
            effort=sum(n_days_effort),
            daily.rate=count/effort,
            monthly.rate=daily.rate*30,
            per100=daily.rate*100)
print(a, n=100)
  



lm_results4 <- mod.dat.spat %>%
  filter(n_days_effort > 10) %>%
  group_by(species_common_name) %>%
  group_modify(~ {
    species <- .y$species_common_name
    if(species %in%c("bike", "cougar", "gray wolf")){
      presence <- glmmTMB(presence ~ location_comments + (1 | month) + (1 | year) + (1 | location),
                          offset = log(n_days_effort),
                          data = .x %>% mutate(presence = ifelse(counts == 0, 0, 1)),
                          family = binomial(link = "logit")
      )
      
      abundance <- glmmTMB(counts ~ location_comments + (1 | month) + (1 | year) + (1 | location),
                           offset = log(n_days_effort),
                           data = .x %>% filter(counts > 0),
                           family = truncated_poisson()
      )
    }else if(species %in%c("elk", "white-tailed deer")){
      presence <- glmmTMB(presence ~ location_comments + MAP_scale + CameraModel + (1 | month) + (1 | year) + (1 | location),
                          offset = log(n_days_effort),
                          data = .x %>% mutate(presence = ifelse(counts == 0, 0, 1)),
                          family = binomial(link = "logit")
      )
      
      abundance <- glmmTMB(counts ~ location_comments + MAP_scale + CameraModel + (1 | month) + (1 | year) + (1 | location),
                           offset = log(n_days_effort),
                           data = .x %>% filter(counts > 0),
                           family = truncated_poisson()
      )
    }else if(species %in%c("mule deer")){
      presence <- glmmTMB(presence ~ location_comments + evi_scale + CameraModel + (1 | month) + (1 | year) + (1 | location),
                          offset = log(n_days_effort),
                          data = .x %>% mutate(presence = ifelse(counts == 0, 0, 1)),
                          family = binomial(link = "logit")
      )
      
      abundance <- glmmTMB(counts ~ location_comments + evi_scale + CameraModel + (1 | month) + (1 | year) + (1 | location),
                           offset = log(n_days_effort),
                           data = .x %>% filter(counts > 0),
                           family = truncated_poisson()
      )
    }else{
      presence <- glmmTMB(presence ~ location_comments +  CameraModel + (1 | month) + (1 | year) + (1 | location),
                          offset = log(n_days_effort),
                          data = .x %>% mutate(presence = ifelse(counts == 0, 0, 1)),
                          family = binomial(link = "logit")
      )
      
      abundance <- glmmTMB(counts ~ location_comments + CameraModel + (1 | month) + (1 | year) + (1 | location),
                           offset = log(n_days_effort),
                           data = .x %>% filter(counts > 0),
                           family = truncated_poisson()
      )
    }
    
    
    newdat <- data.frame(
      location_comments = c("Landscape-valley", "Wildlife trail", "Rec trail"),
      CameraModel = factor("HP2X", levels = levels(as.factor(mod.dat.spat$CameraModel))),
      month = factor("July", levels = levels(as.factor(mod.dat.spat$month))),
      year  = factor("2024", levels = levels(as.factor(mod.dat.spat$year))),
      location = factor("PHATB2", levels = levels(as.factor(mod.dat.spat$location))),
      evi_scale=mean(mod.dat.spat$evi_scale),
      MAP_scale=mean(mod.dat.spat$MAP_scale),
      n_days_effort = 30 # Adjust effort scaling as needed
    )
    # --- Predict conditional (count) component on link scale ---
    pred_abund <- predict(abundance, newdat, type = "response", se.fit = TRUE, re.form = NULL)
    
    # Back-transform to response scale (expected count given presence)
    ### do se
    newdat$cond_fit <- pred_abund$fit
    newdat$cond_se <- pred_abund$se.fit
    newdat$cond_lower <- pred_abund$fit - pred_abund$se.fit
    newdat$cond_upper <- pred_abund$fit + pred_abund$se.fit
    
    
    # --- Predict probability of presence ---
    zi_link <- predict(presence, newdat, type = "response", se.fit = TRUE, re.form = NULL)
    newdat$prob_pres <- zi_link$fit # already on probability scale
    
    # --- Combine: expected total = Pr(presence) × conditional abundance ---
    newdat$expected_total <- newdat$cond_fit * newdat$prob_pres
    newdat$expected_lower <- newdat$cond_lower * newdat$prob_pres
    newdat$expected_upper <- newdat$cond_upper * newdat$prob_pres
    
    ##add ranefs
    # newdat$cam.ranef.pres <- ranef(presence)$cond$CameraModel
    # newdat$cam.ranef.abund <- ranef(abundance)$cond$CameraModel
    newdat
  }) %>%
  ungroup() %>%
  mutate(
    location_comments = factor(location_comments, levels = location_lvs)
    # expected_upper = case_when(
    #   expected_upper > 100000 ~ 0,
    #   TRUE ~ expected_upper
    

  )



mean.detections <- mod.dat %>%
  filter(n_days_effort > 10) %>%
  group_by(species_common_name, location_comments) %>%
  summarise(
    mean = mean(hitrate),
    prop.one = mean(hitrate > 0)
  )


compare <-
  lm_results4 %>%
  left_join(mean.detections, by = c("species_common_name", "location_comments")) %>%
  mutate(
    abs.dif = prob_pres - prop.one,
    abs.perc = (abs(abs.dif) / prop.one) * 100,
    count.dif = expected_total - mean,
    count.perc = (abs(count.dif) / mean) * 100
  ) %>%
  select(species_common_name, location_comments, cond_fit, expected_total, mean, count.dif, count.perc, prob_pres, count.perc, prop.one, abs.dif, abs.perc)


plot(log(compare$expected_total), log(compare$mean))
abline(0, 1, col = "red", lwd = 2) 
cor(log(compare$expected_total), log(compare$mean))

plot(compare$prob_pres, compare$prop.one)
abline(0, 1, col = "red", lwd = 2) 
cor(compare$prob_pres, compare$prop.one)







ggplot(
  lm_results4 %>%
    mutate(location_comments = fct_relevel(location_comments, "Landscape-valley", "Wildlife trail", "Rec trail")),
  aes(
    x = location_comments, y = expected_total,
    ymin = expected_lower, ymax = expected_upper,
    fill = location_comments
  )
) +
  geom_col() +
  geom_linerange(linewidth = 1) +
  facet_wrap(~species_common_name, scales = "free", ncol = 3) +
  scale_fill_manual(values = c(
    "Rec trail" = "#E76F51",
    "Wildlife trail" = "#2A9D8F",
    "Landscape-valley" = "grey60"
  )) +
  labs(x = "Location", y = "Counts per month", title = "Wildlife use") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )


# ─── 0. LIBRARIES ─────────────────────────────────────────────────────────────────
library(glmmTMB)
library(emmeans)
library(dplyr)
library(forcats)
library(ggplot2)
library(stringr)
library(tidyr)
library(ggtext)

# ─── 1. FACTOR LEVELS ───────────────────────────────────────────────────────────────
mod.dat <- mod.dat.spat %>%
  mutate(location_comments = factor(location_comments,
    levels = c("Landscape-valley", "Wildlife trail", "Rec trail")
  ))

# ─── 2. HELPER: TWO CUSTOM CONTRASTS PER SPECIES ────────────────

get_two_contrasts <- function(dat, spname) {
  dat <- dat %>%
    mutate(
      presence = ifelse(counts == 0, 0, 1)
    )

  # Fit presence model (binary)
  if(spname %in%c("bike", "cougar", "gray wolf")){
    presence <- glmmTMB(presence ~ location_comments + (1 | month) + (1 | year) + (1 | location),
                        offset = log(n_days_effort),
                        data = dat %>% mutate(presence = ifelse(counts == 0, 0, 1)),
                        family = binomial(link = "logit")
    )
    
    abundance <- glmmTMB(counts ~ location_comments + (1 | month) + (1 | year) + (1 | location),
                         offset = log(n_days_effort),
                         data = dat %>% filter(counts > 0),
                         family = truncated_poisson()
    )
  }else if(spname  %in%c("elk", "white-tailed deer")){
    presence <- glmmTMB(presence ~ location_comments + MAP_scale + CameraModel + (1 | month) + (1 | year) + (1 | location),
                        offset = log(n_days_effort),
                        data =  dat%>% mutate(presence = ifelse(counts == 0, 0, 1)),
                        family = binomial(link = "logit")
    )
    
    abundance <- glmmTMB(counts ~ location_comments + MAP_scale + CameraModel + (1 | month) + (1 | year) + (1 | location),
                         offset = log(n_days_effort),
                         data = dat %>% filter(counts > 0),
                         family = truncated_poisson()
    )
  }else if(spname  %in%c("mule deer")){
    presence <- glmmTMB(presence ~ location_comments + evi_scale + CameraModel + (1 | month) + (1 | year) + (1 | location),
                        offset = log(n_days_effort),
                        data = dat %>% mutate(presence = ifelse(counts == 0, 0, 1)),
                        family = binomial(link = "logit")
    )
    
    abundance <- glmmTMB(counts ~ location_comments + evi_scale + CameraModel + (1 | month) + (1 | year) + (1 | location),
                         offset = log(n_days_effort),
                         data = dat %>% filter(counts > 0),
                         family = truncated_poisson()
    )
  }else{
    presence <- glmmTMB(presence ~ location_comments +  CameraModel + (1 | month) + (1 | year) + (1 | location),
                        offset = log(n_days_effort),
                        data = dat %>% mutate(presence = ifelse(counts == 0, 0, 1)),
                        family = binomial(link = "logit")
    )
    
    abundance <- glmmTMB(counts ~ location_comments + CameraModel + (1 | month) + (1 | year) + (1 | location),
                         offset = log(n_days_effort),
                         data = dat %>% filter(counts > 0),
                         family = truncated_poisson()
    )
  }
  # Define custom contrasts in the order: Landscape, Wildlife trail, Rec trail
  my_contrasts <- list(
    Landscape_vs_MtFernie = c(1, -0.5, -0.5),
    Wild_vs_Single        = c(0, 1, -1),
    Landscape_vs_Wildlife = c(1, -1, 0)
  )

  # --- Conditional abundance component ---
  emm_c <- emmeans(
    abundance,
    ~location_comments,
    type = "response",
    offset = log(30) # keep consistent with predictions
  )
  pw_c <- contrast(emm_c, my_contrasts) %>%
    summary(infer = TRUE) %>%
    as_tibble() %>%
    mutate(component = "frequency")

  # --- Presence probability component ---
  emm_p <- emmeans(
    presence,
    ~location_comments,
    type = "response"
  )
  pw_p <- contrast(emm_p, my_contrasts) %>%
    summary(infer = TRUE) %>%
    as_tibble() %>%
    mutate(component = "presence")

  # Combine contrasts
  bind_rows(pw_c, pw_p) %>%
    mutate(
      species_common_name = spname,
      # map contrast names to clean group labels
      group1 = case_when(
        contrast == "Landscape_vs_MtFernie" ~ "Landscape-valley",
        contrast == "Wild_vs_Single" ~ "Wildlife trail",
        contrast == "Landscape_vs_Wildlife" ~ "Landscape-valley"
      ),
      group2 = case_when(
        contrast == "Landscape_vs_MtFernie" ~ "Mt Fernie",
        contrast == "Wild_vs_Single" ~ "Rec trail",
        contrast == "Landscape_vs_Wildlife" ~ "Wildlife trail"
      ),
      # star labels for quick plot labels
      label = case_when(
        p.value < 0.05 ~ "**",
        p.value < 0.1 ~ "*",
        TRUE ~ "ns"
      )
    )
}

##check for overdispersion
od.compile <- tibble()
for (i in 1:length(wild_species.all)){
od.test <- glmmTMB(
  counts ~ location_comments + 
    (1 | CameraModel) + (1 | month) + (1  | year)+ 
    (1 | location),
  offset = log(n_days_effort),
  family = truncated_poisson(),
  data = mod.dat %>% filter(species_common_name==wild_species.all[i], counts > 0))
  
  a <- DHARMa::testDispersion(od.test)
  od.compile <- rbind(od.compile,
                      tibble(sp=wild_species.all[i],
                             dispersion=a$statistic,
                             p=a$p.value))
}

simres <- DHARMa::simulateResiduals(od.test)
plot(simres)

write_csv(od.compile , here::here("output", "tables", "overdispersion.csv"))
###test for spatial autocorrelation
library(DHARMa)
library(spdep)
ac.compile <- tibble()
for (i in 1:length(wild_species.all)){
  monthly <- mod.dat %>% filter(species_common_name==wild_species.all[i], counts > 0)%>%
    st_as_sf()%>%
    st_transform(3005)%>%
    cbind(st_coordinates(.))%>%
    st_drop_geometry()
  ac.test <- glmmTMB(
    counts ~ location_comments + 
      (1 | CameraModel) + (1 | month) + (1  | year)+ 
      (1 | location),
    offset = log(n_days_effort),
    family = truncated_poisson(),
    data = monthly)
  
  res <- simulateResiduals(ac.test)
  sp_res <- residuals(res)
  coords <- cbind(monthly$X, monthly$Y)
  
  # create neighbor list: within 5 km for example
  nb <- dnearneigh(coords, 0, 20000)
  lw <- nb2listw(nb, style="W", zero.policy=TRUE)
  
  
  
  a <- moran.test(sp_res, lw, zero.policy=TRUE)
  ac.compile <- rbind(ac.compile,
                      tibble(sp=wild_species.all[i],
                             statistic=a$statistic,
                             estimate=a$estimate[[1]],
                             p=a$p.value))
  
}

ac.compile
write_csv(ac.compile, here::here("output", "tables", "spatial_autocorr.csv"))
library(ncf)
correlog_elkin <- correlog(monthly$X, monthly$Y,
                           sp_res,
                           increment = 250)
plot(correlog_elkin, main="elk spatial correlelogram")

##distances
cam.locs.m <-mod.dat%>%
  distinct(location, .keep_all=TRUE)%>%
  st_as_sf()%>%
  st_transform(3005)%>%
  cbind(st_coordinates(.))%>%
  select(location, location_comments, X, Y)%>%
  st_drop_geometry()

# 1. Pairwise distance matrix
dist_mat <- as.matrix(dist(cam.locs.m[, c("X","Y")]))

# 2. Add row/colnames for cameras
rownames(dist_mat) <- cam.locs.m$location
colnames(dist_mat) <- cam.locs.m$location

# 3. Convert to long format
dist_long <- as.data.frame(dist_mat) %>%
  mutate(cam1 = rownames(.)) %>%
  pivot_longer(-cam1, names_to = "cam2", values_to = "distance")

dist_clean <- dist_long %>%
  filter(cam1 < cam2)   # this automatically keeps only one of each pair


dist_clean <- dist_clean %>%
  left_join(cam.locs.m %>% select(location, location_comments),
            by = c("cam1" = "location")) %>%
  rename(comment1 = location_comments) %>%
  left_join(cam.locs.m %>% select(location, location_comments),
            by = c("cam2" = "location")) %>%
  rename(comment2 = location_comments)

hist(dist_clean$distance)


dist_clean_mtfernie <- dist_clean%>%
  filter(comment1%in%c("Rec trail", "Wildlife trail"),
         comment2%in%c("Rec trail", "Wildlife trail")
         )

hist(dist_clean_mtfernie$distance)


# ─── 3. RUN FOR EACH SPECIES ───────────────────────────────────────────────────────
pairwise_df_raw <- mod.dat %>%
  group_by(species_common_name) %>%
  group_map(~ get_two_contrasts(.x, .y$species_common_name)) %>%
  bind_rows()

## convert to simpler format for inputting values into paper
##flip direction
iou_stats <- pairwise_df_raw%>%
  mutate(ratio=1/ratio,
         odds.ratio=1/odds.ratio)%>%
  select(species_common_name,
         a=group2,
         b=group1,
         ratio,
         odds.ratio,
         p.value,
         label)


# ─── 4. MAX ERROR‐BAR HEIGHT PER SPECIES ──────────────────────────────────────────
max_err <- lm_results4 %>%
  group_by(species_common_name) %>%
  summarize(y_max = max(expected_upper, na.rm = TRUE), .groups = "drop")

# ─── 5. PREPARE STAGGERED ANNOTATIONS ──────────────────────────────────────────────

plot_labels <- pairwise_df_raw %>%
  # filter(label != "ns") %>%                      # drop non‐significant
  # filter(contrast!="Landscape_vs_MtFernie")%>%
  select(species_common_name, contrast, group1, group2, component, label) %>%
  pivot_wider(
    names_from  = component,
    values_from = label,
    values_fill = "ns"
  ) %>%
  left_join(max_err, by = "species_common_name") %>%
  mutate(
    # assign order: 1 = MtFernie bar, 2 = ST vs Wild bar
    idx = case_when(
      contrast == "Landscape_vs_MtFernie" ~ 3L,
      contrast == "Wild_vs_Single" ~ 2L,
      contrast == "Landscape_vs_Wildlife" ~ 1L
    ),
    # stagger heights: e.g. 1.05× and 1.10×
    y_position = (y_max * .98) * (1 + 0.15 * idx),
    # numeric x positions for segment endpoints
    g1 = as.numeric(factor(group1, levels = location_lvs)),
    g2 = case_when(
      contrast == "Landscape_vs_MtFernie" ~ 3,
      TRUE ~ as.numeric(factor(group2, levels = location_lvs))
    ),
    xmid = case_when(
      contrast == "Landscape_vs_MtFernie" ~ 1.75,
      TRUE ~ (g1 + g2) / 2
    ),
    # combined label
    comp_label = case_when(
      presence != "ns" & frequency != "ns" ~ paste0("p", presence, ", ", "f", frequency),
      presence != "ns" & frequency == "ns" ~ paste0("p", presence),
      presence == "ns" & frequency != "ns" ~ paste0("f", frequency),
      presence == "ns" & frequency == "ns" ~ ""
    )
  )


# ─── 6. FINAL PLOT ────────────────────────────────────────────────────────────────
wild_plot <- ggplot(
  lm_results4 %>%
    filter(species_common_name %in% wild_species.all) %>%
    mutate(location_comments = factor(location_comments, levels = location_lvs)),
  aes(
    x    = location_comments,
    y    = expected_total,
    fill = location_comments,
    ymin = expected_lower,
    ymax = expected_upper
  )
) +
  geom_col() +
  geom_linerange(linewidth = 1) +
  facet_wrap(~species_common_name, scales = "free", ncol = 3) +
  scale_fill_manual(values = c(
    "Landscape-valley" = "grey60",
    "Wildlife trail"  = "#2A9D8F",
    "Rec trail"    = "#E76F51"
  )) +

  # two staggered comparison bars
  geom_segment(
    data = plot_labels %>% filter(species_common_name %in% wild_species.all),
    aes(x = g1, xend = g2, y = y_position, yend = y_position),
    inherit.aes = FALSE,
    size = 0.8
  ) +
  geom_text(
    data = plot_labels %>% filter(species_common_name %in% wild_species.all),
    aes(x = xmid, y = y_position * .93, label = comp_label),
    inherit.aes = FALSE,
    fill = NA, # no box around text
    label.color = NA, # no outline
    size = 3,
    hjust = 0.5
  ) +
  labs(
    x     = "Location",
    y     = "Counts per month",
    title = "Wildlife"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )


human_plot <- ggplot(
  lm_results4 %>%
    filter(species_common_name %in% human_species),
  aes(
    x    = location_comments,
    y    = expected_total,
    fill = location_comments,
    ymin = expected_lower,
    ymax = expected_upper
  )
) +
  geom_col() +
  geom_linerange(linewidth = 1) +
  facet_wrap(~species_common_name, scales = "free", ncol = 1) +
  scale_fill_manual(values = c(
    "Landscape-valley" = "grey60",
    "Wildlife trail"  = "#2A9D8F",
    "Rec trail"    = "#E76F51"
  )) +

  # two staggered comparison bars
  geom_segment(
    data = plot_labels %>% filter(species_common_name %in% human_species),
    aes(x = g1, xend = g2, y = y_position, yend = y_position),
    inherit.aes = FALSE,
    size = 0.8
  ) +
  geom_text(
    data = plot_labels %>% filter(species_common_name %in% human_species),
    aes(x = xmid, y = y_position * .93, label = comp_label),
    inherit.aes = FALSE,
    fill = NA, # no box around text
    label.color = NA, # no outline
    size = 3,
    hjust = 0.5
  ) +
  labs(
    x     = "Location",
    y     = "Counts per month",
    title = "Human"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x     = element_text(angle = 50, hjust = 1),
    legend.position = "none"
  )


# Combine side-by-side using patchwork
human_plot + wild_plot + plot_layout(widths = c(1, 3))

## save
ggsave("output/plots/hit_rate_compare_hurdle.png", dpi = 300, width = 10, height = 9, unit = "in")




# 3rd order  ----------------------------------------------


## prep data for modelling
weekly.mod <- cam.summary.week %>%
  filter(species_common_name %in% c(wild_species, human_species)) %>%
  mutate(hitrate = (counts / n_days_effort) * 7)

weekly.mod.wild <- weekly.mod %>%
  filter(wild == "Wild")

weekly.mod.human <- weekly.mod %>%
  filter(wild != "Wild") %>%
  select(location, year, week, species_common_name, hitrate) %>%
  pivot_wider(names_from = species_common_name, values_from = hitrate)

weekly.mod.wild <- weekly.mod.wild %>%
  left_join(weekly.mod.human, by = c("location", "year", "week"))


## add in lag

weekly.mod.wild.lag <- weekly.mod.wild %>%
  select(location:week, species_common_name, presence.lag = presence) %>%
  mutate(
    year = case_when(
      week == 52 ~ year + 1,
      TRUE ~ year
    ),
    week = case_when(
      week < 52 ~ week + 1,
      week == 52 ~ 1
    )
  )


weekly.mod.wild <- weekly.mod.wild %>%
  left_join(weekly.mod.wild.lag, by = c("location", "year", "week", "species_common_name")) %>%
  drop_na(presence.lag)

## extract spatial data to points
# Define your folders
folder1 <- here::here("data/spatial/stacked")
folder2 <- here::here("data/spatial/stacked_focal")

# List all .tif files recursively
files <- list.files(path = c(folder1, folder2), pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)

# Load and stack rasters
r_stack <- rast(files) # terra::rast auto-stacks them if same extent/proj

# Convert sf to SpatVector
points_vect <- cam.locs %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(st_crs(r_stack)) %>%
  vect()

# Extract values
vals <- terra::extract(r_stack, points_vect)

# Combine with original point IDs
cam.locs.attributed <- bind_cols(cam.locs, vals) %>%
  select(location, bldg.count_250:tri)

weekly.mod.wild <- weekly.mod.wild %>%
  left_join(cam.locs.attributed, by = "location")

## standardize

weekly.mod.wild <- weekly.mod.wild %>%
  mutate(
    bike.sqrt = sqrt(bike),
    human.sqrt = sqrt(human),
    ppl.sqrt = sqrt(bike + human)
  ) %>%
  mutate(across(bike:ppl.sqrt, ~ scale(.)[, 1], .names = "{.col}_scale"))

library(corrplot)
library(dplyr)

# Assume your tibble is called `covariates`
# Step 1: Compute correlation matrix
cor_mat <- weekly.mod.wild %>%
  select(list.files(path = c(folder1), pattern = "\\.tif$", full.names = TRUE, recursive = TRUE) %>% basename() %>% str_remove(".tif")) %>%
  select(where(is.numeric)) %>%
  select(where(~ sd(., na.rm = TRUE) > 0)) %>%
  cor(use = "pairwise.complete.obs")

# Step 2: Mask correlations below 0.5 (both positive and negative)
cor_mat_filtered <- cor_mat
cor_mat_filtered[abs(cor_mat) < 0.5] <- NA

# Step 3: Plot
library(corrplot)
corrplot(cor_mat_filtered,
  method = "color", na.label = " ", type = "upper",
  tl.col = "black", tl.cex = 0.8, number.cex = 0.7, addCoef.col = "black"
)

## what sort of variation do we have to work with?
weekly.mod.wild %>%
  filter(location_comments != "Wildlife trail", !location %in% c("ROAD3", "ROAD4")) %>%
  distinct(location, year, week, bike, human) %>%
  mutate(value = bike + human) %>%
  ggplot(aes(x = value)) +
  geom_histogram()

weekly.mod.wild %>%
  filter(location_comments != "Wildlife trail", !location %in% c("ROAD3", "ROAD4")) %>%
  distinct(location, year, week, bike, human) %>%
  mutate(value = bike + human) %>%
  group_by(location) %>%
  summarise(mean = mean(value)) %>%
  ggplot(aes(x = location, y = mean)) +
  geom_col()

weekly.mod.wild %>%
  filter(location_comments != "Wildlife trail", !location %in% c("ROAD3", "ROAD4")) %>%
  distinct(location, year, week, bike, human) %>%
  group_by(location) %>%
  summarise(
    bike = mean(bike),
    human = mean(human)
  ) %>%
  ggplot(aes(x = bike, y = human)) +
  geom_point()


# cor(weekly.mod.wild$human, weekly.mod.wild$`domestic dog`)
sp.loop <- (weekly.mod.wild$species_common_name %>% unique())
loop.dat <- tibble()
# + `Domestic Dog_scale`
for (i in 1:length(sp.loop)) {
  library(broom.mixed)
  loop.dat <- glmmTMB(
    counts ~ bike_scale + human_scale + presence.lag + elev_scale + canopy.all_scale + slope_deg_250_scale + bldg.dist_scale +
      (1 | week) + (1 | year) + (1 | location), # Count model (for counts > 0)
    ziformula = ~ bike_scale + human_scale + presence.lag + elev_scale + canopy.all_scale + slope_deg_250_scale + bldg.dist_scale, # Zero-hurdle model
    data = weekly.mod.wild %>%
      filter(species_common_name == sp.loop[i], !location %in% c("ROAD3", "ROAD4")) %>%
      filter(
        species_common_name == sp.loop[i],
        location_comments != "Wildlife trail"
      ),
    # mutate(trail.width=case_when(location_comments=="Road"~2, TRUE~1)),
    # filter(species_common_name == sp.loop[i]),
    family = truncated_poisson()
  ) %>%
    tidy() %>%
    mutate(sp = sp.loop[i]) %>%
    rbind(loop.dat)
}

loop.dat %>%
  filter(!term %in% c("sd__(Intercept)", "(Intercept)")) %>%
  mutate(component = case_when(
    component == "zi" ~ "Absence",
    TRUE ~ "Frequency"
  )) %>%
  # filter(term%in%c("Bike_scale", "Human_scale"))%>%
  ggplot(aes(y = term, x = estimate, color = sp, xmin = estimate - std.error, xmax = estimate + std.error)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(position = position_dodge(width = 0.5)) +
  facet_wrap(vars(component), scales = "free") +
  geom_vline(xintercept = 0, linetype = "dashed")




loop.dat %>%
  filter(term %in% c("bike_scale", "human_scale")) %>%
  filter(p.value < .1) %>%
  arrange(component)
loop.dat %>%
  filter(p.value < .05) %>%
  filter(term != "(Intercept)") %>%
  print(n = 100)



library(broom.mixed)
#cor(weekly.mod.wild$human, weekly.mod.wild$`domestic dog`)
cor(weekly.mod.wild$human, weekly.mod.wild$bike)
sp.loop <- (weekly.mod.wild$species_common_name %>% unique())
loop.dat <- tibble()
# + `Domestic Dog_scale`
for (i in 1:length(sp.loop)) {
  # presence <- glmmTMB(presence ~ bike.sqrt_scale + human.sqrt_scale + presence.lag + canopy.all_scale + slope_deg_250_scale + bldg.dist_scale + (1|week) + (1|year) + (1 | location),
  #                  data = weekly.mod.wild%>%filter(species_common_name==sp.loop[i], !location %in% c("ROAD3", "ROAD4")),
  #                  family = binomial(link = "logit"))
  #
  # abundance  <- glmmTMB(counts ~ bike.sqrt_scale + human.sqrt_scale + canopy.all_scale + slope_deg_250_scale + bldg.dist_scale + (1|week) + (1|year) + (1 | location),
  #                       data = weekly.mod.wild%>%filter(species_common_name==sp.loop[i], !location %in% c("ROAD3", "ROAD4"))%>%filter(counts>0),
  #                       family = truncated_poisson())

  # presence <- glmmTMB(presence ~ location_comments  + presence.lag + canopy.all_scale + slope_deg_250_scale + bldg.dist_scale + (1|week) + (1|year) + (1 | location),
  #                     data = weekly.mod.wild%>%filter(species_common_name==sp.loop[i], !location%in%c("ROAD3", "ROAD4"))%>%ilter(!str_detect(location, "ROAD|WILD2")),
  #                     family = binomial(link = "logit"))
  #
  # abundance  <- glmmTMB(counts ~ location_comments  + presence.lag + canopy.all_scale + slope_deg_250_scale + bldg.dist_scale + (1|week) + (1|year) + (1 | location),
  #                       data = weekly.mod.wild%>%filter(species_common_name==sp.loop[i], !location%in%c("ROAD3", "ROAD4"))%>%ilter(!str_detect(location, "ROAD|WILD2"))%>%filter(counts>0),
  #                       family = truncated_poisson())

  presence <- glmmTMB(presence ~ ppl.sqrt_scale + canopy.all_scale + slope_deg_250_scale + bldg.dist_scale + (1 | week) + (1 | year) + (1 | location),
    data = weekly.mod.wild %>% filter(species_common_name == sp.loop[i]),
    family = binomial(link = "logit")
  )

  abundance <- glmmTMB(counts ~ ppl.sqrt_scale + canopy.all_scale + slope_deg_250_scale + bldg.dist_scale + (1 | week) + (1 | year) + (1 | location),
    data = weekly.mod.wild %>% filter(species_common_name == sp.loop[i]) %>% filter(counts > 0),
    family = truncated_poisson()
  )



  loop.dat <- presence %>%
    tidy() %>%
    mutate(
      sp = sp.loop[i],
      type = "Presence"
    ) %>%
    rbind(abundance %>%
      tidy() %>%
      mutate(
        sp = sp.loop[i],
        type = "Frequency"
      )) %>%
    rbind(loop.dat)
}

my_sp_order <- c(
  "black bear",
  "elk",
  "grizzly bear",
  "moose",
  "mule deer",
  "red fox",
  "white-tailed deer"
)



# pick your palette
library(RColorBrewer)
pal <- MRColorBrewerpal <- MetBrewer::met.brewer("Hokusai1", length(my_sp_order))
pal <- brewer.pal(7,"Dark2")

loop.dat %>%
  filter(!term %in% c("sd__(Intercept)", "(Intercept)")) %>%
  mutate(
    type = fct_relevel(type, "Presence", "Frequency"),
    term = case_when(
      term == "slope_deg_250_scale" ~ "slope",
      term == "bike.sqrt_scale" ~ "bike",
      term == "ppl.sqrt_scale" ~ "recreationists",
      term == "human.sqrt_scale" ~ "human",
      term == "canopy.all_scale" ~ "canopy cover",
      term == "bldg.dist_scale" ~ "building distance",
      TRUE ~ term
    ) %>% fct_relevel("human", "bike", "slope", "canopy cover", "building distance", "presence.lag"),
    sig = p.value < 0.1,
    sp = factor(sp, levels = rev(my_sp_order)),
    fill_col = ifelse(sig, as.character(sp), "white")
  ) %>%
  ggplot(
    aes(
      y = term,
      x = estimate,
      color = sp, # outline
      fill = fill_col, # fill directly
      xmin = estimate - std.error,
      xmax = estimate + std.error,
      group = sp
    )
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_linerange(position = position_dodge(width = 0.6)) +
  geom_point(
    shape = 21,
    position = position_dodge(width = 0.6),
    size = 3
  ) +
  facet_wrap(vars(type), scales = "free") +
  scale_color_manual(
    values = setNames(pal, my_sp_order),
    name = "Species",
    breaks = my_sp_order
  ) +
  scale_fill_manual(
    values = c(setNames(pal, my_sp_order), white = "white"),
    guide = "none"
  ) +
  labs(x = "Estimate ± SE", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90"),
    legend.position = "right"
  )

ggsave("output/plots/rec.users.impact.png", dpi = 300, width = 8, height = 6, unit = "in", bg = "white")



loop.dat %>%
  filter(term %in% c("bike.sqrt_scale", "human.sqrt_scale", "ppl.sqrt_scale")) %>%
  filter(p.value < .1) %>%
  arrange(type)



loop.dat %>%
  filter(term %in% c("ppl.sqrt_scale")) %>%
  arrange(sp, type)




loop.dat %>%
  filter(!term %in% c("sd__(Intercept)", "(Intercept)")) %>%
  group_by(term, type) %>%
  summarise(
    est = mean(estimate),
    n = sum(p.value < .1)
  ) %>%
  arrange(type)


# Daily activity  ----------------------------------------------

sp.loop <- wild_species

fernie.activity <- detections %>%
  ## add lat long
  left_join(
    cam.locs %>%
      select(
        location, latitude, longitude, location_comments
      ),
    by = c("location")
  ) %>%
  ##
  mutate(mid_time = as_datetime((as.numeric(start_time) + as.numeric(end_time)) / 2)) %>%
  ## calculate sun rise/set
  mutate(solar_time = solartime(mid_time, # the date time column
    latitude,
    longitude,
    tz = -6, # an offset in numeric hours to UTC (Alberta is 6 hours behind on MDT, summer time when the cams were deployed)
    format = "%Y-%m-%d %H:%M:%S"
  )$solar)

## other data


sr.activity <- external.detections %>%
  filter(species_common_name %in% sp.loop) %>%
  # filter(location_comments!="Wildlife-sr")%>%
  ## calculate sun rise/set
  mutate(solar_time = solartime(start_time, # the date time column
    X,
    Y,
    tz = -6, # an offset in numeric hours to UTC (Alberta is 6 hours behind on MDT, summer time when the cams were deployed)
    format = "%Y-%m-%d %H:%M:%S"
  )$solar)


sr.activity.inverse <- external.detections.inverse %>%
  filter(species_common_name %in% sp.loop) %>%
  # filter(location_comments!="Wildlife-sr")%>%
  ## calculate sun rise/set
  mutate(solar_time = solartime(start_time, # the date time column
    X,
    Y,
    tz = -6, # an offset in numeric hours to UTC (Alberta is 6 hours behind on MDT, summer time when the cams were deployed)
    format = "%Y-%m-%d %H:%M:%S"
  )$solar)

# # Fit an activity model
# bike.actv <- fernie.activity %>%
#   filter(species_common_name == "Bike") %>%
#   pull(solar_time) %>%
#   fitact(sample = "data", reps = 10)
#
# hiker.actv <- fernie.activity %>%
#   filter(species_common_name == "Human") %>%
#   pull(solar_time) %>%
#   fitact(sample = "data", reps = 10)
#
#
#
#
# plot(bike.actv,
#   yunit = "density", data = "none", las = 1, lwd = 2,
#   tline = list(lwd = 2,lty=2), # Thick line
#   cline = list(lty = 0)
# ) # Supress confidence intervals
#
# plot(hiker.actv,
#   yunit = "density", data = "none", add = TRUE,
#   tline = list(col = "red", lwd = 2),
#   cline = list(lty = 0)
# )
#
# legend("topright", c("Bike", "Human"), col = 1:2, lty = c(2,1), lwd = 2)


df.activity <- fernie.activity %>%
  select(location, latitude, longitude, species_common_name, mid_time, solar_time, location_comments) %>%
  rbind(
    sr.activity %>%
      tibble() %>%
      mutate(location_comments = "Wildlife-sr") %>%
      select(location, latitude = Y, longitude = X, species_common_name, mid_time = start_time, solar_time, location_comments)
  ) %>%
  rbind(
    sr.activity.inverse %>%
      tibble() %>%
      mutate(location_comments = "Wildlife-sr-inverse") %>%
      select(location, latitude = Y, longitude = X, species_common_name, mid_time = start_time, solar_time, location_comments)
  )


df.activity%>%
  group_by(location_comments)%>%
  summarise(n=n_distinct(location))


df.activity %>%
  filter(species_common_name == "elk") %>%
  group_by(location_comments) %>%
  count()

## do bikes outside to speed things up
human <- df.activity %>%
  filter(species_common_name %in% c("bike", "human")) %>%
  pull(solar_time)

human.actv <- human %>% fitact(sample = "data", reps = 10)

vehic <- df.activity %>%
  filter(species_common_name %in% c("vehicle")) %>%
  pull(solar_time)

vehic.actv <- vehic %>% fitact(sample = "data", reps = 10)





###LOOP
plot(human.actv, yunit = "density")
plot(vehic.actv, yunit = "density", add = TRUE)
diel.results <- tibble()
diel.dif.results <- tibble()
for (sp in sp.loop) {
  ## pull out species data
  st <- df.activity %>%
    filter(
      species_common_name == sp,
      location_comments == "Rec trail"
    ) %>%
    pull(solar_time)

  wild.near <- df.activity %>%
    filter(
      species_common_name == sp,
      location_comments == "Wildlife trail"
    ) %>%
    pull(solar_time)

  wild.far <- df.activity %>%
    filter(
      species_common_name == sp,
      location_comments == "Wildlife-sr"
    ) %>%
    pull(solar_time)

  wild.inverse <- df.activity %>%
    filter(
      species_common_name == sp,
      location_comments == "Wildlife-sr-inverse"
    ) %>%
    pull(solar_time)

  min.dat <- min(c(length(st), length(wild.near), length(wild.far), length(wild.inverse)))

  # if( min.dat>10){
  # Estimate overlap coefficients
  # overlap.st <- overlapEst(human, st)
  # overlap.wt <- overlapEst(human, wild.near)

  # Bootstrap to estimate confidence intervals
  set.seed(123)
  boot.st <- if (length(st) > 10) {
    bootstrap(human, st, 10, type = "Dhat4")
  } else {
    NA
  }
  boot.wt <- if (length(wild.near) > 10) {
    bootstrap(human, wild.near, 10, type = "Dhat4")
  } else {
    NA
  }
  boot.sr <- if (length(wild.far) > 10) {
    bootstrap(human, wild.far, 10, type = "Dhat4")
  } else {
    NA
  }
  boot.sr.inv <- if (length(wild.inverse) > 10) {
    bootstrap(human, wild.inverse, 10, type = "Dhat4")
  } else {
    NA
  }
  
  # Get confidence intervals
  ci.st <- if (length(length(st) > 10)) {
    quantile(boot.st, c(0.05, 0.5, 0.95)) %>% round(2)
  } else {
    rep(NA, 3)
  }
  ci.wt <- if (length(wild.near) > 10) {
    quantile(boot.wt, c(0.05, 0.5, 0.95)) %>% round(2)
  } else {
    rep(NA, 3)
  }
  ci.sr <- if (length(wild.far) > 10) {
    quantile(boot.sr, c(0.05, 0.5, 0.95)) %>% round(2)
  } else {
    rep(NA, 3)
  }
  ci.sr.inv <- if (length(wild.inverse) > 10) {
    quantile(boot.sr.inv, c(0.05, 0.5, 0.95)) %>% round(2)
  } else {
    rep(NA, 3)
  }

  # Bootstrap the differences between the two overlap estimates
  bootDiff <- if (length(wild.near) > 10) {
    boot.st - boot.wt
  } else {
    NA
  }
  bootDiff.sr <- if (length(wild.far) > 10) {
    boot.st - boot.sr
  } else {
    NA
  }
  bootDiff.sr.inv <- if (length(wild.inverse) > 10) {
    boot.st - boot.sr.inv
  } else {
    NA
  }

  # Confidence interval of the difference
  ci.dif <- if (length(wild.near) > 10) {
    quantile(bootDiff, c(0.05, 0.5, 0.95)) %>% round(2)
  } else {
    rep(NA, 3)
  }
  ci.dif.sr <- if (length(wild.far) > 10) {
    quantile(bootDiff.sr, c(0.05, 0.5, 0.95)) %>% round(2)
  } else {
    rep(NA, 3)
  }
  ci.dif.sr.inv <- if (length(wild.inverse) > 10) {
    quantile(bootDiff.sr.inv, c(0.05, 0.5, 0.95)) %>% round(2)
  } else {
    rep(NA, 3)
  }

  # If CI doesn't include 0 → significant difference in overlap

  # Build the annotation text
  label_text <- paste0(
    "Overlap with human period on trails:",
    "\nRec trail = ", ci.st[2],
    "\nWildlife Trail = ", ci.wt[2],
    "\nLandscape-Valley = ", ci.sr[2],
    "\nLandscape-Remote = ", ci.sr.inv[2]
    # "\nDif (st-wt) = ", ci.dif[2], " [", ci.dif[1], ", ", ci.dif[3], "]",
    # "\nDif (st-land.val) = ", ci.dif.sr[2], " [", ci.dif.sr[1], ", ", ci.dif.sr[3], "]",
    # "\nDif (st-land.rem) = ", ci.dif.sr.inv[2], " [", ci.dif.sr.inv[1], ", ", ci.dif.sr.inv[3], "]"
  )

  # save results
  diel.results <- diel.results %>%
    rbind(tibble(
      trail = c("wildlife trail", "Rec trail", "landscape-valley", "landscape-remote"),
      overlap = c(ci.wt[2], ci.st[2], ci.sr[2], ci.sr.inv[2]),
      overlap.lower = c(ci.wt[1], ci.st[1], ci.sr[1], ci.sr.inv[1]),
      overlap.upper = c(ci.wt[3], ci.st[3], ci.sr[3], ci.sr.inv[3]),
      species = sp
    ))

  diel.dif.results <- diel.dif.results %>%
    rbind(tibble(
      trail = c("Rec trail—wildlife trail", "Rec trail—landscape-valley", "Rec trail—landscape-remote"),
      group1 = c("Rec trail"),
      group2 = c("wildlife trail", "landscape-valley", "landscape-remote"),
      dif = c(ci.dif[2], ci.dif.sr[2], ci.dif.sr.inv[2]),
      dif.lower = c(ci.dif[1], ci.dif.sr[1], ci.dif.sr.inv[1]),
      dif.upper = c(ci.dif[3], ci.dif.sr[3], ci.dif.sr.inv[3]),
      species = sp
    ))




  plot.act <- function(sp = "Mule Deer", rep.num = 10, bandwidth = NULL) {
    ## estimate activity
    actv.st <- st %>%
      activity::fitact(sample = "data", reps = rep.num, bw = bandwidth)

    actv.wild.near <- wild.near %>%
      activity::fitact(sample = "data", reps = rep.num, bw = bandwidth)

    actv.wild.far <- wild.far %>%
      activity::fitact(sample = "data", reps = rep.num, bw = bandwidth)

    actv.wild.inv <- wild.inverse %>%
      activity::fitact(sample = "data", reps = rep.num, bw = bandwidth)



    plot(human.actv,
      yunit = "density", data = "none",
      tline = list(col = "black", lwd = 2, lty = "dashed"),
      cline = list(lty = 0),
      main = sp
    )

    plot(actv.st,
      yunit = "density", data = "none", las = 1, lwd = 2, add = TRUE,
      tline = list(lwd = 2), # Thick line
      cline = list(lty = 0)
    ) # Supress confidence intervals

    plot(actv.wild.near,
      yunit = "density", data = "none", add = TRUE,
      tline = list(col = "red", lwd = 2),
      cline = list(lty = 0)
    )

    plot(actv.wild.far,
      yunit = "density", data = "none", add = TRUE,
      tline = list(col = "blue", lwd = 2),
      cline = list(lty = 0)
    )

    plot(actv.wild.inv,
      yunit = "density", data = "none", add = TRUE,
      tline = list(col = "orange", lwd = 2),
      cline = list(lty = 0)
    )


    # legend("topleft", c("Humans", "Rec trail", "Wildlife Trail-near", "Wildlife Trail-far"), col = c(1, 1, 2, 4), lty = c(2, 1, 1, 1), lwd = 2)
    legend("topleft", c("Bikers and Hikers", "Rec trail", "Wildlife Trail", "Landscape-valley", "Landscape-remote"),
      col = c("black", "black", "red", "blue", "orange"),
      lty = c("dashed", "solid", "solid", "solid", "solid"),
      lwd = 2
    )
    # Add to top right of the plot
    text(
      x = par("usr")[2] - 0.2, y = par("usr")[4] - 0.002, labels = label_text,
      adj = c(1, 1), cex = 0.9, font = 2
    )
  }


  # Save the activity plot to a file
  png(
    filename = paste0("/Users/claytonlamb/Dropbox/Documents/University/Work/RecreationMonitoring/MtFernie_RecreationEcology/output/plots/diel/activity_overlap_", gsub(" ", "_", sp), ".png"),
    width = 1200, height = 900, res = 150
  )


  bw <- case_when(
    min.dat <= 20 ~ 5,
    min.dat >= 20 & min.dat <= 50 ~ 10,
    min.dat > 50 ~ 15
  )


  # Create the plot
  plot.act(sp = sp, bandwidth = bw)

  # Close the graphics device
  dev.off()

  ## remove data
  rm(
    st, wild.near, wild.far, wild.inverse,
    boot.st, boot.wt, boot.sr, boot.sr.inv,
    bootDiff, bootDiff.sr, bootDiff.sr.inv,
    ci.st, ci.wt, ci.sr, ci.sr.inv, ci.dif, ci.dif.sr, ci.dif.sr.inv,
    actv.st, actv.wild.near, actv.wild.far, actv.wild.inv,
    label_text, bw
  )

  # }
}


### plot

# 1. Mark significance (non-overlap with 0)
diel.dif.results$trail %>% unique()
diel.dif.results <- diel.dif.results %>%
  mutate(
    sig_label = ifelse(dif.lower > 0 | dif.upper < 0, "*", ""),
    trail = fct_relevel(trail, "Rec trail—wildlife trail", "Rec trail—landscape-valley", "Rec trail—landscape-remote")
  )

# 2. Plot
# Reverse species order
diel.dif.results <- diel.dif.results %>%
  mutate(species = factor(species, levels = rev(unique(species))))

# Plot
ggplot(diel.dif.results, aes(x = species, y = dif, ymin = dif.lower, ymax = dif.upper, color = trail)) +
  geom_point(position = position_dodge(width = 0.2)) +
  geom_linerange(position = position_dodge(width = 0.2)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(
    aes(label = sig_label),
    position = position_dodge(width = 0.2),
    vjust = -0.5,
    show.legend = FALSE
  ) +
  coord_flip() +
  labs(
    y = "Difference in overlap with hiker/biker period",
    x = "Species",
    color = "Comparison"
  ) +
  theme.custom +
  scale_color_manual(values = c(
    "Rec trail—landscape-valley" = "#E76F51",
    "Rec trail—landscape-remote" = "#2A9D8F"
  ))

## save
ggsave("output/plots/diel/dif.overlap.png", dpi = 300, width = 7, height = 6, unit = "in", bg = "white")

ggplot(diel.results %>% mutate(species = factor(species, levels = rev(unique(species)))), aes(x = species, y = overlap, ymin = overlap.lower, ymax = overlap.upper, color = fct_relevel(trail, "across landscape", "wildlife trail", "Rec trail"))) +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_linerange(position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(
    y = "Overlap with peak hiker/biker period",
    x = "Species",
    color = "Comparison"
  ) +
  theme.custom +
  scale_color_brewer(palette = "Dark2")

## save
ggsave("output/plots/diel/diel.overlap.png", dpi = 300, width = 7, height = 6, unit = "in", bg = "white")


###
max.upper <- max(diel.results$overlap.upper, na.rm = TRUE)
diel_plot_stats <- diel.dif.results %>%
  group_by(species) %>%
  mutate(
    y_rank = row_number(),
    y_base = max.upper,
    y_position = y_base * 1.05 + ((y_rank - 1) * 0.1 * y_base) # Adjust spacing here
  )%>%
  drop_na(dif)

trail.order <- c("Rec trail", "wildlife trail", "landscape-valley", "landscape-remote")
ggplot(
  diel.results,
  aes(x = fct_relevel(trail, trail.order), y = overlap, ymin = overlap.lower, ymax = overlap.upper, fill = fct_relevel(trail, trail.order))
) +
  geom_col() +
  geom_linerange(linewidth = 1) +
  facet_wrap(~species, , ncol = 3) +
  geom_segment(
    data = diel_plot_stats,
    aes(x = group1, xend = group2, y = y_position, yend = y_position),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = diel_plot_stats,
    aes(
      x = as.numeric(factor(group1, levels = c(trail.order))) * 0.5 +
        as.numeric(factor(group2, levels = c(trail.order))) * 0.5,
      y = y_position * 1.02,
      label = sig_label
    ),
    inherit.aes = FALSE,
    size = 4
  ) +
  scale_fill_manual(values = c(
    "Rec trail" = "#E76F51",
    "wildlife trail" = "#2A9D8F",
    "landscape-valley" = "grey80",
    "landscape-remote" = "grey50"
  )) +
  labs(x = "Location", y = "Overlap (proportion)", title = "") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave("output/plots/diel/diel.comparison.png", dpi = 300, width = 7, height = 8, unit = "in", bg = "white")

# ggplot(diel.results%>%
#          mutate(trail=fct_relevel(trail, "Rec trail","wildlife trail", "across landscape")), aes(x=species, y=overlap, ymin=overlap.lower, ymax=overlap.upper, color=trail))+
#   geom_point(position = position_dodge(width=0.2))+
#   geom_linerange(position = position_dodge(width=0.2))+
#   geom_hline(yintercept=0, linetype="dashed")







## For FTA sign
# Save the plot as PNG
png("output/plots/diel/fta_activity_plot.png", width = 3000, height = 2000, res = 300)

plot(human.actv,
  xunit = "hour",
  yunit = "density", data = "none",
  tline = list(col = "black", lwd = 2, lty = "dashed"),
  cline = list(lty = 0),
  ylab = "Activity intensity",
  yaxt = "n", # suppress y-axis labels
  main = "Human and Wildlife Detections on Mt. Fernie Trails"
)
plot(
  fernie.activity %>%
    filter(species_common_name == "black bear") %>%
    pull(solar_time) %>%
    fitact(sample = "data", reps = 10),
  xunit = "hour",
  yunit = "density", data = "none", las = 1, lwd = 2, add = TRUE,
  tline = list(lwd = 2, col = "black"), # Thick line
  cline = list(lty = 0)
) # Supress confidence intervals

plot(
  fernie.activity %>%
    filter(species_common_name == "grizzly bear") %>%
    pull(solar_time) %>%
    fitact(sample = "data", reps = 10),
  xunit = "hour",
  yunit = "density", data = "none", las = 1, lwd = 2, add = TRUE,
  tline = list(lwd = 2, col = "brown"), # Thick line
  cline = list(lty = 0)
) # Supress confidence intervals

plot(
  fernie.activity %>%
    filter(species_common_name == "moose") %>%
    pull(solar_time) %>%
    fitact(sample = "data", reps = 10),
  xunit = "hour",
  yunit = "density", data = "none", las = 1, lwd = 2, add = TRUE,
  tline = list(lwd = 2, col = "orange"), # Thick line
  cline = list(lty = 0)
) # Supress confidence intervals

plot(
  fernie.activity %>%
    filter(species_common_name == "elk") %>%
    pull(solar_time) %>%
    fitact(sample = "data", reps = 10),
  xunit = "hour",
  yunit = "density", data = "none", las = 1, lwd = 2, add = TRUE,
  tline = list(lwd = 2, col = "gold"), # Thick line
  cline = list(lty = 0),
) # Supress confidence intervals

plot(
  fernie.activity %>%
    filter(species_common_name == "mule deer") %>%
    pull(solar_time) %>%
    fitact(sample = "data", reps = 10),
  xunit = "hour",
  yunit = "density", data = "none", las = 1, lwd = 2, add = TRUE,
  tline = list(lwd = 2, col = "red"), # Thick line
  cline = list(lty = 0),
) # Supress confidence intervals

plot(
  fernie.activity %>%
    filter(species_common_name == "red fox") %>%
    pull(solar_time) %>%
    fitact(sample = "data", reps = 10),
  xunit = "hour",
  yunit = "density", data = "none", las = 1, lwd = 2, add = TRUE,
  tline = list(lwd = 2, col = "darkgoldenrod"), # Thick line
  cline = list(lty = 0),
) # Supress confidence intervals

plot(
  fernie.activity %>%
    filter(species_common_name == "white-tailed deer") %>%
    pull(solar_time) %>%
    fitact(sample = "data", reps = 10),
  xunit = "hour",
  yunit = "density", data = "none", las = 1, lwd = 2, add = TRUE,
  tline = list(lwd = 2, col = "forestgreen"), # Thick line
  cline = list(lty = 0),
) # Supress confidence intervals




## get anchor times using all the data

# # ##get sunset times for a specific day
# date <- ymd_hms("2021-07-31 18:02:20", tz="MST")
# lat <- 49.49707
# long <- -115.1126
tz <- -6
# activity::solartime(date, lat, long, tz)
#
# # suntimes <- get_suntimes(date, lat, long, tz)[, -3] * pi/12
# # apply(suntimes, 2, cmean)


## radians
# suntimes <- get_suntimes(df.activity$mid_time, df.activity$latitude, df.activity$longitude, tz)[, -3] * pi/12
# anchors <- apply(suntimes, 2, cmean)
# anchors

## hours
suntimes_hr <- get_suntimes(df.activity$mid_time, df.activity$latitude, df.activity$longitude, tz)[, -3]
anchors <- apply(suntimes_hr, 2, mean)
anchors


# Add vertical dashed lines at sunrise and sunset
## radians
# abline(v = 1.69, lty = 2, col = "orange")  # Sunrise
# abline(v = 5.45, lty = 2, col = "blue")    # Sunset
abline(v = anchors[1], lty = 2, col = "orange") # Sunrise
abline(v = anchors[2], lty = 2, col = "blue") # Sunset
axis(side = 1, at = seq(1, 24, by = 1), labels = FALSE)

# Add labels for sunrise and sunset
text(x = anchors[1], y = par("usr")[4] * 0.95, labels = "Sunrise", col = "orange", pos = 3, cex = 0.9, font = 2)
text(x = anchors[2], y = par("usr")[4] * 0.95, labels = "Sunset", col = "blue", pos = 3, cex = 0.9, font = 2)

legend("topleft", c("bikers & hikers", "black bear", "grizzly bear", "moose", "elk", "mule deer", "red fox", "white-tailed deer"), col = c("black", "black", "brown", "orange", "gold", "red", "darkgoldenrod", "forestgreen"), lty = c(2, 1, 1, 1, 1, 1, 1, 1), lwd = 2, bg="white")

# Close the device
dev.off()







# AAR ---------------------------------------------------------------------

## determine common species and rec types to focus on
det.n <- detections %>%
  count(species_common_name) %>%
  arrange(-n)
# rec.aar <- c("Bike", "Vehicle", "Human", "All Terrain Vehicle")
# wild.aar <- c("White-tailed Deer", "Mule Deer", "Red Fox", "Elk", "Moose", "Black Bear", "Grizzly Bear")
# cams <- cam.locs$location
# ## calculate times
# ## think about blackout times for cam, when not active, extract timed photo dates?
#
# aar.combined <- tibble()
# for(i in 1:length(wild.aar)){
# sp <- wild.aar[i]
# aar.bike <- tibble()
# aar.human <- tibble()
# aar.vehicle <- tibble()
# for(j in 1:length(cams)){
# loc <- cams[j]
# aar.bike.j <- detections %>%
#   filter(species_common_name%in%c(rec.aar, sp), location==loc)%>%
#   ## order by time
#   arrange(end_time) %>%
#   ## calculate time since last observation and fill down to next
#   mutate(
#     last_wildlife_time = case_when(species_common_name == sp ~ end_time, TRUE ~ NA_POSIXct_),
#     last_human_time = case_when(species_common_name == "Human" ~ end_time, TRUE ~ NA_POSIXct_),
#     last_bike_time = case_when(species_common_name == "Bike" ~ end_time, TRUE ~ NA_POSIXct_),
#     last_vehicle_time = case_when(species_common_name == "Vehicle" ~ end_time, TRUE ~ NA_POSIXct_),
#     last_atv_time = case_when(species_common_name == "All Terrain Vehicle" ~ end_time, TRUE ~ NA_POSIXct_)
#   ) %>%
#   fill(last_wildlife_time, .direction = "down") %>%
#   fill(last_human_time, .direction = "down") %>%
#   fill(last_bike_time, .direction = "down") %>%
#   fill(last_vehicle_time, .direction = "down") %>%
#   fill(last_atv_time, .direction = "down")%>%
#   ## this only works for a single cam so need to loop through all cams, unfort, need to figure out the group_by with location
#   ## Wildlife count to eventually remove subsequent detections without rec between
#   mutate(wildlife.count = if_else(species_common_name == sp, 1, NA_integer_)) %>%
#   group_by(grp = cumsum(lag(species_common_name, default = first(species_common_name)) != sp)) %>%
#   mutate(wildlife.count = if_else(species_common_name == sp, row_number(), NA_integer_)) %>%
#   ungroup() %>%
#   select(-grp) %>%
#   ## bike count to eventually remove subsequent detections without deer between
#   mutate(bike.count = if_else(species_common_name == "Bike", 1, NA_integer_)) %>%
#   group_by(grp = cumsum(lag(species_common_name, default = first(species_common_name)) != "Bike")) %>%
#   mutate(bike.count = if_else(species_common_name == "Bike", row_number(), NA_integer_)) %>%
#   ungroup() %>%
#   select(-grp) %>%
#   ## calculate HW and WH
#   mutate(
#     HW = case_when(species_common_name == sp & # species is wildlife
#       wildlife.count == 1 & # first count of wildlife in back to back deer sequences
#       !is.na(last_bike_time) & # not before the first bike is detected
#       (last_bike_time > last_human_time | is.na(last_human_time)) & # no humans between
#       (last_bike_time > last_vehicle_time | is.na(last_vehicle_time)) & # no vehicles between
#       (last_bike_time > last_atv_time | is.na(last_atv_time)) # no atvs between
#     ~ as.numeric(difftime(start_time, last_bike_time, units = "hours")), TRUE ~ NA_real_), # calculate time difference in hours.
#     WH = case_when(species_common_name == "Bike" &
#       bike.count == 1 &
#       !is.na(last_wildlife_time) &
#         (last_wildlife_time > last_human_time | is.na(last_human_time)) & # no humans between
#         (last_wildlife_time > last_vehicle_time | is.na(last_vehicle_time)) & # no vehicles between
#         (last_wildlife_time > last_atv_time | is.na(last_atv_time)) # no atvs between
#     ~ as.numeric(difftime(start_time, last_wildlife_time, units = "hours")), TRUE ~ NA_real_)
#   )%>%
#   mutate(species_of_interest=sp,
#          recreation="Bike")%>%
#   pivot_longer(HW:WH)%>%
#   drop_na(value)%>%
#   select(location, species_of_interest, recreation, species_common_name, start_time, end_time, name, value)%>%
#   arrange(location,end_time)
#
#
#
#
# aar.human.j <- detections %>%
#   filter(species_common_name%in%c(rec.aar, sp), location==loc)%>%
#   ## order by time
#   arrange(end_time) %>%
#   ## calculate time since last observation and fill down to next
#   mutate(
#     last_wildlife_time = case_when(species_common_name == sp ~ end_time, TRUE ~ NA_POSIXct_),
#     last_human_time = case_when(species_common_name == "Human" ~ end_time, TRUE ~ NA_POSIXct_),
#     last_bike_time = case_when(species_common_name == "Bike" ~ end_time, TRUE ~ NA_POSIXct_),
#     last_vehicle_time = case_when(species_common_name == "Vehicle" ~ end_time, TRUE ~ NA_POSIXct_),
#     last_atv_time = case_when(species_common_name == "All Terrain Vehicle" ~ end_time, TRUE ~ NA_POSIXct_)
#   ) %>%
#   fill(last_wildlife_time, .direction = "down") %>%
#   fill(last_human_time, .direction = "down") %>%
#   fill(last_bike_time, .direction = "down") %>%
#   fill(last_vehicle_time, .direction = "down") %>%
#   fill(last_atv_time, .direction = "down")%>%
#   ## this only works for a single cam so need to loop through all cams, unfort, need to figure out the group_by with location
#   ## Wildlife count to eventually remove subsequent detections without rec between
#   mutate(wildlife.count = if_else(species_common_name == sp, 1, NA_integer_)) %>%
#   group_by(grp = cumsum(lag(species_common_name, default = first(species_common_name)) != sp)) %>%
#   mutate(wildlife.count = if_else(species_common_name == sp, row_number(), NA_integer_)) %>%
#   ungroup() %>%
#   select(-grp) %>%
#   ## human count to eventually remove subsequent detections without deer between
#   mutate(human.count = if_else(species_common_name == "Human", 1, NA_integer_)) %>%
#   group_by(grp = cumsum(lag(species_common_name, default = first(species_common_name)) != "Human")) %>%
#   mutate(human.count = if_else(species_common_name == "Human", row_number(), NA_integer_)) %>%
#   ungroup() %>%
#   select(-grp) %>%
#   ## calculate HW and WH
#   mutate(
#     HW = case_when(species_common_name == sp & # species is wildlife
#                      wildlife.count == 1 & # first count of wildlife in back to back deer sequences
#                      !is.na(last_human_time) & # not before the first human is detected
#                      (last_human_time > last_bike_time | is.na(last_bike_time)) & # no bikes between
#                      (last_human_time > last_vehicle_time | is.na(last_vehicle_time)) & # no vehicles between
#                      (last_human_time > last_atv_time | is.na(last_atv_time)) # no atvs between
#                    ~ as.numeric(difftime(start_time, last_human_time, units = "hours")), TRUE ~ NA_real_), # calculate time difference in hours.
#     WH = case_when(species_common_name == "Human" &
#                      human.count == 1 &
#                      !is.na(last_wildlife_time) &
#                      (last_wildlife_time > last_bike_time | is.na(last_bike_time)) & # no bikes between
#                      (last_wildlife_time > last_vehicle_time | is.na(last_vehicle_time)) & # no vehicles between
#                      (last_wildlife_time > last_atv_time | is.na(last_atv_time)) # no atvs between
#                    ~ as.numeric(difftime(start_time, last_wildlife_time, units = "hours")), TRUE ~ NA_real_)
#   )%>%
#   mutate(species_of_interest=sp,
#          recreation="Human")%>%
#   pivot_longer(HW:WH)%>%
#   drop_na(value)%>%
#   select(location, species_of_interest, recreation, species_common_name, start_time, end_time, name, value)%>%
#   arrange(location,end_time)
#
#
# aar.vehicle.j <- detections %>%
#   filter(species_common_name%in%c(rec.aar, sp), location==loc)%>%
#   ## order by time
#   arrange(end_time) %>%
#   ## calculate time since last observation and fill down to next
#   mutate(
#     last_wildlife_time = case_when(species_common_name == sp ~ end_time, TRUE ~ NA_POSIXct_),
#     last_human_time = case_when(species_common_name == "Human" ~ end_time, TRUE ~ NA_POSIXct_),
#     last_bike_time = case_when(species_common_name == "Bike" ~ end_time, TRUE ~ NA_POSIXct_),
#     last_vehicle_time = case_when(species_common_name == "Vehicle" ~ end_time, TRUE ~ NA_POSIXct_),
#     last_atv_time = case_when(species_common_name == "All Terrain Vehicle" ~ end_time, TRUE ~ NA_POSIXct_)
#   ) %>%
#   fill(last_wildlife_time, .direction = "down") %>%
#   fill(last_human_time, .direction = "down") %>%
#   fill(last_bike_time, .direction = "down") %>%
#   fill(last_vehicle_time, .direction = "down") %>%
#   fill(last_atv_time, .direction = "down")%>%
#   ## this only works for a single cam so need to loop through all cams, unfort, need to figure out the group_by with location
#   ## Wildlife count to eventually remove subsequent detections without rec between
#   mutate(wildlife.count = if_else(species_common_name == sp, 1, NA_integer_)) %>%
#   group_by(grp = cumsum(lag(species_common_name, default = first(species_common_name)) != sp)) %>%
#   mutate(wildlife.count = if_else(species_common_name == sp, row_number(), NA_integer_)) %>%
#   ungroup() %>%
#   select(-grp) %>%
#   ## human count to eventually remove subsequent detections without deer between
#   mutate(vehicle.count = if_else(species_common_name == "Vehicle", 1, NA_integer_)) %>%
#   group_by(grp = cumsum(lag(species_common_name, default = first(species_common_name)) != "Vehicle")) %>%
#   mutate(vehicle.count = if_else(species_common_name == "Vehicle", row_number(), NA_integer_)) %>%
#   ungroup() %>%
#   select(-grp) %>%
#   ## calculate HW and WH
#   mutate(
#     HW = case_when(species_common_name == sp & # species is wildlife
#                      wildlife.count == 1 & # first count of wildlife in back to back deer sequences
#                      !is.na(last_vehicle_time) & # not before the first human is detected
#                      (last_vehicle_time > last_bike_time | is.na(last_bike_time)) & # no bikes between
#                      (last_vehicle_time > last_human_time | is.na(last_human_time)) & # no humans between
#                      (last_vehicle_time > last_atv_time | is.na(last_atv_time)) # no atvs between
#                    ~ as.numeric(difftime(start_time, last_human_time, units = "hours")), TRUE ~ NA_real_), # calculate time difference in hours.
#     WH = case_when(species_common_name == "Vehicle" &
#                      vehicle.count == 1 &
#                      !is.na(last_wildlife_time) &
#                      (last_wildlife_time > last_bike_time | is.na(last_bike_time)) & # no bikes between
#                      (last_wildlife_time > last_human_time | is.na(last_human_time)) & # no humans between
#                      (last_wildlife_time > last_atv_time | is.na(last_atv_time)) # no atvs between
#                    ~ as.numeric(difftime(start_time, last_wildlife_time, units = "hours")), TRUE ~ NA_real_)
#   )%>%
#   mutate(species_of_interest=sp,
#          recreation="Vehicle")%>%
#   pivot_longer(HW:WH)%>%
#   drop_na(value)%>%
#   select(location, species_of_interest, recreation, species_common_name, start_time, end_time, name, value)%>%
#   arrange(location,end_time)
#
# aar.bike <- rbind(aar.bike, aar.bike.j)
# aar.human <- rbind(aar.human, aar.human.j)
# aar.vehicle <- rbind(aar.vehicle, aar.vehicle.j)
# }
# aar.combined <- rbind(aar.combined, aar.bike, aar.human,aar.vehicle)
# }


#
# ###TEST DATA TO MAKE SURE THIS WORKS
# test.dat <- tribble(
#   ~location,  ~species_common_name, ~start_time,
#   "A", "Bike", "2021-10-21 18:02:20",
#   "A", "Bike", "2021-10-21 18:05:21",
#   "A", "Human", "2021-10-21 18:10:21",
#   "A", "Human", "2021-10-21 23:05:21",
#   "A", "Bike", "2021-10-22 06:05:21",
#   "A", "Mule Deer", "2021-10-22 06:20:21",
#   "A", "Bike", "2021-10-22 06:40:21",
#   "A", "Bike", "2021-10-22 06:50:21",
#   "A", "Motorized vehicle", "2021-10-22 06:55:21",
#   "B", "Motorized vehicle", "2021-10-22 06:55:21",
#   "B", "Mule Deer", "2021-10-22 08:30:21",
#   "B", "Mule Deer", "2021-10-22 08:40:21",
#   "B", "Motorized vehicle", "2021-10-22 08:50:21",
#   "B", "Motorized vehicle", "2021-10-22 08:55:21",
#   "B", "Mule Deer", "2021-10-22 09:30:21",
#   "B", "Human", "2021-10-22 09:31:21",
#   "B", "Mule Deer", "2021-10-22 09:50:21",
#   "B", "Grizzly Bear", "2021-10-22 09:51:21",
#   "B", "Motorized vehicle", "2021-10-22 09:55:21"
# )%>%
#   arrange(start_time)%>%
#   mutate(
#     start_time = as.POSIXct(start_time))
# ##expected
# wild.aar <- c("Moose", "Black Bear", "Grizzly Bear", "Elk", "Mule Deer", "Red Fox", "White-tailed Deer")
# rec.aar <- c("Bike", "Human", "Motorized vehicle")
# cams <- unique(test.dat$location)
# aar.test <- run_t2b_t1_analysis(test.dat, wild.aar, rec.aar, cams)%>%
#   arrange(location, wildlife_time,recreation_time)

##


## POST TAAL PAPER
library(dplyr)
library(furrr)
library(lubridate)
library(future)

# Parallel plan: safe
plan(multisession, workers = parallel::detectCores() - 1)

process_one_group <- function(df, sp, rec_type) {
  df <- df %>%
    arrange(mid_time) %>%
    mutate(
      next_species = lead(species_common_name),
      next_time = lead(mid_time),
      gap_hours = as.numeric(difftime(next_time, mid_time, units = "hours"))
    )

  df %>%
    mutate(
      interaction = case_when(
        species_common_name == sp & next_species == rec_type ~ "WH",
        species_common_name == rec_type & next_species == sp ~ "HW",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(interaction)) %>%
    mutate(
      species_of_interest = sp,
      recreation = rec_type,
      wildlife_time = if_else(interaction == "WH", mid_time, next_time),
      recreation_time = if_else(interaction == "HW", mid_time, next_time)
    ) %>%
    select(location, species_of_interest, recreation, interaction,
      value = gap_hours, wildlife_time, recreation_time
    )
}

run_analysis_vectorized <- function(detections, wild_species, rec_types, cam_locations) {
  # Pre-filter the big detections once
  detections <- detections %>%
    filter(
      location %in% cam_locations,
      species_common_name %in% c(wild_species, rec_types)
    )

  # Group by camera for splitting
  cam_groups <- detections %>%
    group_by(location) %>%
    group_split()

  # Use future_map over cameras (efficient chunking)
  results <- future_map_dfr(cam_groups, function(cam_df) {
    loc <- unique(cam_df$location)
    out <- list()
    for (sp in wild_species) {
      for (rec_type in rec_types) {
        df_sub <- cam_df %>%
          filter(species_common_name %in% c(sp, rec_type))
        if (nrow(df_sub) > 1) {
          out[[length(out) + 1]] <- process_one_group(df_sub, sp, rec_type)
        }
      }
    }
    bind_rows(out)
  })

  results
}

# Example usage
wild.aar <- c("Moose", "Black Bear", "Grizzly Bear", "Cougar", "Gray wolf") %>%
  str_to_lower()
rec.aar <- c("Bike", "Human", "Motorized vehicle") %>% str_to_lower()


detections.aar <- detections3min %>%
  mutate(
    species_common_name = case_when(
      species_common_name %in% c("All Terrain Vehicle", "Vehicle") ~ "Motorized vehicle",
      TRUE ~ species_common_name
    ),
    start_time = as.POSIXct(start_time),
    end_time = as.POSIXct(end_time),
    mid_time = start_time + (as.numeric(difftime(end_time, start_time, units = "secs")) / 2),
    dur = end_time - start_time
  )

cams <- unique(detections.aar$location)

aar.combined <- run_analysis_vectorized(detections.aar, wild.aar, rec.aar, cams)

plan(sequential)
gc()

aar.combined <- aar.combined %>%
  rename(species = species_of_interest)

write_csv(aar.combined, here::here("output/tables/aar.csv"))

aar.combined <- read_csv(here::here("output/tables/aar.csv"))

ggplot(aar.combined %>% filter(value < 24), aes(x = interaction, y = value, fill = "species")) +
  geom_boxplot() +
  facet_wrap(vars(recreation))



close.call.dat <- aar.combined %>%
  left_join(
    cam.summary %>%
      distinct(location, location_comments),
    by = "location"
  ) %>%
  filter(location_comments == "Rec trail") %>%
  mutate(wild.event.id=paste0(species,location,wildlife_time),
           value = value * 60) %>%
  mutate(close = case_when(
    value < 30 ~ 1,
    TRUE ~ 0
  )) %>%
  filter(
    recreation %in% c("bike", "human")
    #species %in% c("black bear", "grizzly bear", "moose")
  ) %>%
  mutate(interaction = case_when(
    interaction == "HW" ~ "Human then wildlife",
    interaction == "WH" ~ "Wildlife then human",
  ) %>% fct_relevel("Wildlife then human", "Human then wildlife"))

close.call.dat %>%
  group_by(species, wild.event.id )%>%
  mutate(close2=max(close))%>%
  summarise(time=min(value, na.rm=TRUE),
            close.call=mean(close2))%>%
  group_by(species) %>%
  summarise(
    total_events = n_distinct(wild.event.id),
    within30 = sum(close.call),
    percwithin30 = (within30/total_events) * 100,
    within5 = sum(time <= 5, na.rm=TRUE),
    percwithin5 = (within5/total_events) * 100
  )

ggplot(
  close.call.dat %>% filter(close == 1),
  aes(x = value, fill = recreation)
) +
  geom_histogram() +
  facet_grid(species ~ interaction, scales = "free_y") +
  theme.custom +
  labs(x = "Minutes between", y = "Number of events", fill = "Recreation") +
  scale_fill_manual(values = c(
    "bike" = "#E76F51",
    "human" = "#2A9D8F"
  ))
ggsave("output/plots/close.calls.png", dpi = 300, width = 8, height = 6, unit = "in", bg = "white")






## percent within
library(dplyr)
library(lubridate)

# Example wild + rec species
wild_species <- c("moose", "black bear", "grizzly bear", "cougar", "gray wolf")
rec_types <- c("bike", "human")

detections.aar2 <- detections.aar %>%
  mutate(
    type = case_when(
      species_common_name %in% wild_species ~ "wildlife",
      species_common_name %in% rec_types ~ "rec",
      TRUE ~ NA_character_
    )
  )

wild <- detections.aar2 %>%
  filter(type == "wildlife")

rec <- detections.aar2 %>%
  filter(type == "rec")

# Pair every rec with all wildlife detections of interest
pairs <- rec %>%
  rename(rec_time = mid_time, rec_type = species_common_name) %>%
  inner_join(
    wild %>%
      rename(wild_time = mid_time, species = species_common_name),
    by = "location"
  ) %>%
  mutate(
    time_diff_mins = as.numeric(difftime(wild_time, rec_time, units = "mins")),
    abs_time_diff_mins = abs(time_diff_mins)
  ) %>%
  # Keep only the *closest* wildlife detection of each species per rec detection
  group_by(rec_time, rec_type, location, species) %>%
  slice_min(order_by = abs_time_diff_mins, n = 1, with_ties = FALSE) %>%
  ungroup()

pairs


# Count total recs per type
total_recs <- rec %>%
  count(rec_type = species_common_name) %>%
  rename(total_recs = n)

# Count matches within 30 min per rec_type + species
matched <- pairs %>%
  filter(abs_time_diff_mins <= 30) %>%
  count(species, rec_type) %>%
  rename(matched_recs = n)

# Combine and calculate %
result <- matched %>%
  left_join(total_recs, by = "rec_type") %>%
  mutate(percent = (matched_recs / total_recs) * 100) %>%
  arrange(desc(percent))

result




close.call.dat2 <- pairs %>%
  left_join(
    cam.summary %>%
      distinct(location, location_comments),
    by = "location"
  ) %>%
  filter(location_comments == "Rec trail") %>%
  mutate(close = case_when(
    abs_time_diff_mins < 30 ~ 1,
    TRUE ~ 0
  )) %>%
  filter(
    rec_type %in% c("bike", "human"),
    species %in% wild_species
  )

close.call.dat2 %>%
  group_by(species) %>%
  summarise(
    total_events = n(),
    within30 = sum(close),
    percwithin30 = mean(close) * 100,
    within5 = sum(abs_time_diff_mins <= 5),
    percwithin5 = mean(abs_time_diff_mins <= 5) * 100
  )


library(ggplot2)
library(patchwork)

p_bb <- ggplot(close.call.dat2 %>% filter(close==1, species == "black bear"),
               aes(x = abs_time_diff_mins, fill = rec_type)) +
  geom_histogram() +
  scale_y_continuous(breaks = seq(0, 12, by = 2)) +
  labs(title = "black bear") +
  theme_minimal() +
  labs(x = "", y = "Number of events", fill = "Recreation") +
  scale_fill_manual(values = c(
    "bike" = "#E76F51",
    "human" = "#2A9D8F"
  ))+
  theme(legend.position = "none")

p_gb <- ggplot(close.call.dat2 %>% filter(close==1, species == "grizzly bear"),
               aes(x = abs_time_diff_mins, fill = rec_type)) +
  geom_histogram() +
  scale_y_continuous(breaks = seq(0, 3, by = 1)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA))+
  labs(title = "grizzly bear") +
  theme_minimal() +
  labs(x = "Minutes between", y = "", fill = "Recreation") +
  scale_fill_manual(values = c(
    "bike" = "#E76F51",
    "human" = "#2A9D8F"
  ))+
  theme(legend.position = "none")


p_moose <- ggplot(close.call.dat2 %>% filter(close==1, species == "moose"),
                  aes(x = abs_time_diff_mins, fill = rec_type)) +
  geom_histogram() +
  scale_y_continuous(breaks = seq(0, 12, by = 2)) +
  labs(title = "moose") +
  theme_minimal() +
  labs(x = "", y = "", fill = "Recreation") +
  scale_fill_manual(values = c(
    "bike" = "#E76F51",
    "human" = "#2A9D8F"
  ))+
  theme(legend.position = "none")

p_coug <- ggplot(close.call.dat2 %>% filter(close==1, species == "cougar"),
                  aes(x = abs_time_diff_mins, fill = rec_type)) +
  geom_histogram() +
  scale_y_continuous(breaks = seq(0, 2, by = 1)) +
  labs(title = "cougar") +
  theme_minimal() +
  labs(x = "", y = "", fill = "Recreation") +
  scale_fill_manual(values = c(
    "bike" = "#E76F51",
    "human" = "#2A9D8F"
  ))

combined <- p_bb + p_gb + p_moose + p_coug + plot_layout(ncol = 4) & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")

ggsave("output/plots/close.calls2.png", dpi = 300, width = 8, height = 4, unit = "in", bg = "white")








# aar.combined%>%filter(value<24)%>%
#   group_by(species,recreation, name, location, month(end_time))%>%
#   summarise(median=median(value))%>%
#   pivot_wider(values_from=median, names_from=name)%>%
#   mutate(aar=HW/WH)%>%
#   left_join(cam.locs)%>%
#   ggplot(aes(x=species, y=log(aar), fill=recreation))+
#   geom_violin(alpha = 0.5) +
#   geom_dotplot(binaxis = "y",
#                stackdir = "center",
#                dotsize = 0.5,
#                position = position_dodge())+
#   geom_hline(yintercept = 0, linetype="dashed")+
#   facet_wrap(vars(location_comments), ncol=1)+
#   theme.custom
#
# aar.combined%>%filter(value<24)%>%
#   group_by(species,recreation, name, location, month(end_time))%>%
#   summarise(median=median(value))%>%
#   pivot_wider(values_from=median, names_from=name)%>%
#   mutate(aar=HW/WH)%>%
#   left_join(cam.locs)%>%
#   filter(location_comments=="Road")%>%
#   #filter(location_comments=="Rec trail", recreation!="Vehicle")%>%
#   #filter(location_comments=="Wildlife trail")%>%
#   ggplot(aes(x=species, y=log(aar), fill=recreation))+
#   geom_violin(alpha = 0.5) +
#   geom_dotplot(binaxis = "y",
#                stackdir = "center",
#                dotsize = 0.3,
#                position = position_dodge())+
#   geom_hline(yintercept = 0, linetype="dashed")+
#   facet_wrap(vars(location_comments), ncol=1)+
#   theme.custom+
#   labs(x="Species", y="Avoidance ratio")
#
# aar.combined%>%filter(value<24)%>%
#   filter(location%in%c("ROAD3", "ROAD4"))%>%
#   group_by(species,recreation, name, location, month(end_time))%>%
#   summarise(median=median(value))%>%
#   pivot_wider(values_from=median, names_from=name)%>%
#   mutate(aar=HW/WH)%>%
#   left_join(cam.locs)%>%
#   filter(location_comments=="Road")%>%
#   #filter(location_comments=="Rec trail", recreation!="Vehicle")%>%
#   #filter(location_comments=="Wildlife trail")%>%
#   ggplot(aes(fill=species, y=log(aar), x=recreation))+
#   geom_violin(alpha = 0.5) +
#   geom_dotplot(binaxis = "y",
#                stackdir = "center",
#                dotsize = 0.3,
#                position = position_dodge())+
#   geom_hline(yintercept = 0, linetype="dashed")+
#   facet_wrap(vars(location_comments), ncol=1)+
#   theme.custom+
#   labs(x="Species", y="Avoidance ratio")
#
#
# aar.combined%>%filter(value<24)%>%
#   group_by(species,recreation, name, location, month(end_time))%>%
#   summarise(median=median(value))%>%
#   pivot_wider(values_from=median, names_from=name)%>%
#   mutate(aar=HW/WH)%>%
#   left_join(cam.locs)%>%
#   #filter(location_comments=="Road")%>%
#   filter(location_comments=="Rec trail", recreation!="Vehicle")%>%
#   #filter(location_comments=="Wildlife trail")%>%
#   ggplot(aes(fill=species, y=log(aar), x=recreation))+
#   geom_violin(alpha = 0.5) +
#   geom_dotplot(binaxis = "y",
#                stackdir = "center",
#                dotsize = 0.3,
#                position = position_dodge())+
#   geom_hline(yintercept = 0, linetype="dashed")+
#   facet_wrap(vars(location_comments), ncol=1)+
#   theme.custom+
#   labs(x="Species", y="Avoidance ratio")
#
# aar.combined%>%filter(value<24)%>%
#   group_by(species,recreation, name, location, month(end_time))%>%
#   summarise(median=median(value))%>%
#   pivot_wider(values_from=median, names_from=name)%>%
#   mutate(aar=HW/WH)%>%
#   left_join(cam.locs)%>%
#   #filter(location_comments=="Road")%>%
#   filter(location_comments=="Rec trail", recreation!="Vehicle")%>%
#   #filter(location_comments=="Wildlife trail")%>%
#   ggplot(aes(fill=species, y=log(aar), x=recreation))+
#   geom_boxplot()+
#   geom_hline(yintercept = 0, linetype="dashed")+
#   facet_wrap(vars(location_comments), ncol=1)+
#   theme.custom+
#   labs(x="Species", y="Avoidance ratio")
#
# aar.combined%>%filter(value<24& value>=0)%>%
#   group_by(species,recreation, name, location, month(end_time))%>%
#   summarise(median=median(value))%>%
#   pivot_wider(values_from=median, names_from=name)%>%
#   mutate(aar=HW/WH)%>%
#   left_join(cam.locs)%>%
#   filter(location_comments=="Road")%>%
#   #filter(location_comments=="Rec trail", recreation!="Vehicle")%>%
#   #filter(location_comments=="Wildlife trail")%>%
#   ggplot(aes(fill=species, y=log(aar), x=recreation))+
#   geom_boxplot()+
#   geom_hline(yintercept = 0, linetype="dashed")+
#   facet_wrap(vars(location_comments), ncol=1)+
#   theme.custom+
#   labs(x="Species", y="Avoidance ratio")
#
# 
# 
# 
# library(glmmTMB)
# 
# 
# aar.dat <- aar.combined %>%
#   filter(value <= 16, value > 0) %>%
#   mutate(month = month(recreation_time)) %>%
#   group_by(species, recreation, wildlife_time, interaction, location, month) %>%
#   add_count() %>% ## remove a couple spots where animals and people detected at same time, which messes things up, only a few of these, but they are interesting!
#   filter(n == 1) %>%
#   select(-recreation_time) %>%
#   pivot_wider(values_from = value, names_from = interaction) %>%
#   mutate(
#     aar = HW / WH,
#     log.aar = log(aar)
#   ) %>%
#   left_join(cam.locs) %>%
#   filter(
#     location_comments != "Wildlife trail",
#     species %in% wild.aar
#   ) %>%
#   drop_na(aar) %>%
#   group_by(species, recreation, location, month) %>%
#   mutate(mean.wh = mean(WH)) %>%
#   ungroup() %>%
#   mutate(
#     WH.scale = scale(WH)[, 1],
#     mean.WH.scale = scale(mean.wh)[, 1]
#   ) %>%
#   mutate(
#     time_decimal = hour(wildlife_time) + minute(wildlife_time) / 60 + second(wildlife_time) / 3600,
#     time_decimal2 = time_decimal^2,
#     hours_from_noon = time_decimal - 12
#   ) %>%
#   mutate(
#     time_decimal.scale = scale(time_decimal)[, 1],
#     time_decimal2.scale = scale(time_decimal2)[, 1]
#   )
# 
# aar.dat %>%
#   ggplot(aes(fill = species, y = log(aar), x = recreation)) +
#   geom_boxplot() +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   theme.custom +
#   labs(x = "Species", y = "Avoidance ratio")
# 
# # Use bayesian model, converges better
# # model<- glmmTMB(log.aar ~ recreation*species_of_interest + (1|location)+ (1|month),
# #         data = aar.dat)
# #
# # model2 <- glmmTMB(log.aar ~ recreation + (1+recreation|species_of_interest) + (1|location)+ (1|month),
# #                  data = aar.dat)
# 
# 
# library(brms)
# 
# # Define the Bayesian model
# # brm_model <- brm(
# #   formula = HW ~ recreation + WH + (1 + recreation | species) + (1 | location) + (1 | month) + time_decimal.scale + time_decimal2.scale,
# #   data = aar.dat,
# #   family = gaussian(),  # log.aar is continuous
# #   chains = 4,
# #   cores = 4,
# #   iter = 4000,
# #   warmup = 1000,
# #   control = list(adapt_delta = 0.99)
# # )
# 
# 
# 
# # Define the Bayesian model
# brm_model <- brm(
#   formula = log.aar ~ recreation + (1 + recreation | species) + (1 | location) + (1 | month) + time_decimal + time_decimal2,
#   data = aar.dat,
#   family = gaussian(), # log.aar is continuous
#   chains = 4,
#   cores = 4,
#   iter = 4000,
#   warmup = 1000,
#   control = list(adapt_delta = 0.99)
# )
# 
# 
# # aar.dat%>%
# #   filter(location%in%c("ROAD4"))%>%
# #   group_by(recreation)%>%
# #   summarize(#aar=median(aar),
# #     log.aar=median(log.aar))
# # ##double check that bike and human effects consistent with overall effects on these busy road cams to ensure pooling OK
# # aar.dat%>%
# #   #filter(!str_detect(location, "WILD"))%>%
# #   #filter(str_detect(location, "ROAD"))%>%
# #   group_by(recreation, location,location_comments)%>%
# #   summarize(#aar=median(aar),
# #     log.aar=mean(log.aar))%>%
# #   ggplot(aes(x=location, y=log.aar, color=location_comments))+
# #   geom_point()+
# #   facet_wrap(vars(recreation))+
# #   geom_hline(yintercept = 0, linetype="dashed")+
# #   theme(axis.text.x = element_text(angle=45, hjust=1))
# 
# 
# 
# # # Step 3: Predict values with SEs
# # preds <- predict(model, newdata = newdata, type = "link", se.fit = TRUE)
# #
# # # Combine predictions with newdata
# # pred_df <- newdata %>%
# #   mutate(
# #     fit = preds$fit,
# #     se = preds$se.fit,
# #     lower = fit - 1.96 * se,
# #     upper = fit + 1.96 * se
# #   )
# #
# # # Step 4: Add overall fixed effects (ignores random species deviations)
# # newdata_overall <- data.frame(
# #   recreation = recreation_types,
# #   species_of_interest = NA,  # force population-level fixed effect
# #   location=NA,
# #   month=7
# # )
# #
# #
# #
# # overall_preds <- predict(model, newdata = newdata_overall, type = "link", se.fit = TRUE, re.form = NA)
# #
# # overall_df <- data.frame(
# #   recreation = recreation_types,
# #   species = "Overall",
# #   fit = overall_preds$fit,
# #   se = overall_preds$se.fit
# # ) %>%
# #   mutate(
# #     lower = fit - 1.96 * se,
# #     upper = fit + 1.96 * se
# #   )
# #
# # # Step 5: Clean up species-level data and combine
# # plot_df <- pred_df %>%
# #   rename(species = species_of_interest)
# #
# # # Step 6: Plot
# # ggplot(plot_df, aes(x = recreation, y = fit, color = species, group = species)) +
# #   geom_point(position = position_dodge(width = 0.6), size = 2) +
# #   geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2,
# #                 position = position_dodge(width = 0.6)) +
# #
# #   geom_point(data=overall_df,
# #              aes(x = recreation, y = fit,  group = species),
# #              color="black") +
# #   geom_errorbar(data=overall_df,
# #                 color="black",
# #                 aes(x = recreation, y = fit,  group = species,ymin = lower, ymax = upper),
# #                 width = 0.1) +
# #   labs(
# #     x = "Recreation Type",
# #     y = "Predicted Effect (log scale)",
# #     color = "Species"
# #   ) +
# #   geom_hline(yintercept = 0, linetype="dashed")+
# #   theme.custom+  theme(axis.text.x = element_text(angle = 30, hjust = 1))
# 
# 
# 
# ## for BRMS
# recreation_types <- unique(aar.dat$recreation)
# species <- unique(aar.dat$species)
# 
# 
# newdata <- data.frame(
#   recreation = "bike",
#   species = "mule Deer",
#   location = NA,
#   month = 7,
#   WH.scale = median(aar.dat$WH.scale),
#   mean.WH.scale = median(aar.dat$mean.WH.scale),
#   ## 8am
#   time_decimal = aar.dat$time_decimal[1:100],
#   time_decimal2 = aar.dat$time_decimal2[1:100]
#   ## 9pm
#   # time_decimal.scale=1.246416,
#   # time_decimal2.scale=1.370596
# )
# 
# 
# preds <- fitted(
#   brm_model,
#   newdata = newdata,
#   re_formula = NULL,
#   summary = TRUE,
#   allow_new_levels = TRUE
# )
# 
# pred_df <- newdata %>%
#   mutate(
#     fit = preds[, "Estimate"],
#     lower = preds[, "Q2.5"],
#     upper = preds[, "Q97.5"]
#   )
# 
# 
# ggplot(pred_df, aes(x = time_decimal, y = fit, ymin = lower, ymax = upper)) +
#   geom_point() +
#   geom_linerange() +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   labs(title = "Modelled result: Mule deer")
# 
# 
# ggplot(aar.dat %>% filter(species == "mule deer"), aes(x = hour(wildlife_time), y = log.aar)) +
#   geom_point() +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   geom_smooth(se = FALSE, method = "loess") +
#   facet_wrap(vars(recreation)) +
#   labs(title = "Raw data: Mule deer")
# 
# ggplot(aar.dat %>% filter(species == "elk"), aes(x = hour(wildlife_time), y = log.aar)) +
#   geom_point() +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   geom_smooth(se = FALSE, method = "loess") +
#   facet_wrap(vars(recreation)) +
#   labs(title = "Raw data: Mule deer")
# 
# 
# newdata <- expand.grid(
#   recreation = recreation_types,
#   species = species,
#   location = NA,
#   month = 7,
#   WH.scale = median(aar.dat$WH.scale),
#   mean.WH.scale = median(aar.dat$mean.WH.scale),
#   ## 8am
#   time_decimal.scale = -0.8964351,
#   time_decimal2.scale = -0.9348653
#   ## 9pm
#   # time_decimal.scale=1.246416,
#   # time_decimal2.scale=1.370596
# )
# 
# 
# preds <- fitted(
#   brm_model,
#   newdata = newdata,
#   re_formula = NULL,
#   summary = TRUE,
#   allow_new_levels = TRUE
# )
# 
# pred_df <- newdata %>%
#   mutate(
#     fit = preds[, "Estimate"],
#     lower = preds[, "Q2.5"],
#     upper = preds[, "Q97.5"]
#   )
# 
# newdata_overall <- data.frame(
#   recreation = recreation_types,
#   species_of_interest = NA,
#   location = NA,
#   month = 7,
#   WH.scale = median(aar.dat$WH.scale),
#   mean.WH.scale = median(aar.dat$mean.WH.scale),
#   ## 8am
#   time_decimal.scale = -0.8964351,
#   time_decimal2.scale = -0.9348653
#   ## 9pm
#   # time_decimal.scale=1.246416,
#   # time_decimal2.scale=1.370596
# )
# 
# overall_preds <- fitted(brm_model, newdata = newdata_overall, re_formula = NA)
# 
# overall_df <- data.frame(
#   recreation = recreation_types,
#   species = "Overall",
#   fit = overall_preds[, "Estimate"],
#   lower = overall_preds[, "Q2.5"],
#   upper = overall_preds[, "Q97.5"]
# ) %>%
#   mutate(recreation = fct_relevel(recreation, "Human", "Bike", "Motorized vehicle"))
# 
# # Create a manual offset by species so they don't sit directly on top of overall
# species_offsets <- c(
#   "Moose" = -0.25,
#   "Black Bear" = -0.2,
#   "Grizzly Bear" = -0.15,
#   "Elk" = -0.1,
#   "Mule Deer" = 0.1,
#   "Red Fox" = 0.15,
#   "White-tailed Deer" = 0.20
# )
# 
# # Join offset into pred_df
# overall_df <- overall_df %>%
#   mutate(recreation = fct_relevel(recreation, "Human", "Bike", "Motorized vehicle"))
# 
# pred_df <- pred_df %>%
#   mutate(
#     x_nudge = as.numeric(as.factor(recreation)) + species_offsets[species],
#     recreation = fct_relevel(recreation, "Human", "Bike", "Motorized vehicle")
#   )
# pred_df <- pred_df %>%
#   mutate(
#     x_nudge = as.numeric(as.factor(recreation)) + species_offsets[species],
#     recreation = fct_relevel(recreation, "Human", "Bike", "Motorized vehicle")
#   )
# 
# 
# 
# library(MetBrewer)
# 
# ggplot() +
#   # Species points
#   geom_point(data = pred_df, aes(x = x_nudge, y = fit, color = species), size = 2) +
#   geom_errorbar(data = pred_df, aes(x = x_nudge, ymin = lower, ymax = upper, color = species), width = 0) +
# 
#   # Overall black points at center
#   geom_point(data = overall_df, aes(x = as.numeric(as.factor(recreation)), y = fit), color = "black", size = 3) +
#   geom_errorbar(data = overall_df, aes(x = as.numeric(as.factor(recreation)), ymin = lower, ymax = upper), width = 0.1, color = "black") +
# 
#   # Dashed baseline
#   geom_hline(yintercept = 0, linetype = "dashed") +
# 
#   # Annotations
#   annotate("text", x = 3.3, y = 0.22, label = "Less frequent after", color = "black", hjust = 0, angle = 90, size = 4, fontface = "italic") +
#   annotate("text", x = 3.3, y = -2.5, label = "More frequent after", color = "black", hjust = 0, angle = 90, , size = 4, fontface = "italic") +
# 
#   # Color palette
#   scale_color_manual(values = met.brewer("Hokusai1", length(unique(pred_df$species)))) +
# 
#   # Fix the x-axis labels
#   scale_x_continuous(
#     breaks = 1:3,
#     labels = levels(as.factor(pred_df$recreation)),
#     name = "Recreation Type"
#   ) +
#   labs(
#     y = "Predicted Effect (log AAR)",
#     color = "Species"
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(
#     axis.title = element_text(size = 16),
#     legend.title = element_text(size = 14),
#     legend.text = element_text(size = 12),
#     panel.grid.minor = element_blank(),
#     panel.grid.major.x = element_blank(),
#     legend.position = c(0.8, 0.82)
#   ) +
#   guides(color = guide_legend(ncol = 2)) +
#   labs(title = "Fine-scale temporal responses")
# 
# 
# ggsave("output/plots/aar.png", dpi = 300, width = 10, height = 6, unit = "in", bg = "white")
# 
# ggplot(aar.dat, aes(x = HW, y = log.aar)) +
#   geom_point(alpha = 0.5) +
#   facet_wrap(vars(recreation, species)) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   theme.custom
# 
# 
# aar.dat %>%
#   filter(WH < 16, HW < 16) %>%
#   group_by(recreation) %>%
#   summarise(
#     HW = mean(HW),
#     WH = mean(WH)
#   ) %>%
#   mutate(aar = HW / WH)
# 
# aar.combined %>%
#   filter(value < 10) %>%
#   mutate(month = month(recreation_time)) %>%
#   group_by(species, recreation, wildlife_time, interaction, location, month) %>%
#   add_count() %>% ## remove a couple spots where animals and people detected at same time, which messes things up, only a few of these, but they are interesting!
#   filter(n == 1) %>%
#   select(-recreation_time) %>%
#   pivot_wider(values_from = value, names_from = interaction) %>%
#   mutate(
#     aar = HW / WH,
#     log.aar = log(aar)
#   ) %>%
#   # filter(WH<16, HW<16)%>%
#   group_by(recreation) %>%
#   summarise(
#     HW = mean(HW, na.rm = TRUE),
#     WH = mean(WH, na.rm = TRUE)
#   ) %>%
#   mutate(aar = HW / WH)
# 
# # ##proof of concept
# # library(dplyr)
# # library(lubridate)
# #
# # # Example data
# # df <- tribble(
# #   ~sp,     ~time,                 ~expected_HW, ~expected_WH,
# #   "deer",  "2024-06-10 01:00:00", NA, NA,
# #   "deer",  "2024-06-10 02:00:00", NA, NA,
# #   "hiker", "2024-06-10 03:00:00", NA, 1,
# #   "deer",  "2024-06-10 03:30:00", 0.5, NA,
# #   "hiker", "2024-06-10 08:00:00", NA, 4.5,
# #   "hiker", "2024-06-10 09:00:00", NA, NA,
# #   "biker", "2024-06-10 09:30:00", NA, NA,
# #   "hiker", "2024-06-10 14:00:00", NA, NA,
# #   "atv", "2024-06-10 14:35:00", NA, NA,
# #   "deer",  "2024-06-10 20:00:00", NA, NA,
# #   "deer",  "2024-06-10 20:20:00", NA, NA,
# #   "hiker", "2024-06-10 20:30:00", NA, 0.167,
# #   "deer",  "2024-06-10 22:30:00", 2, NA,
# #   "deer",  "2024-06-10 23:30:00", NA, NA,
# #   "atv", "2024-06-10 23:35:00", NA, NA,
# #   "deer",  "2024-06-10 23:40:00", NA, NA
# # )
# #
# #
# #
# # # Convert time column to POSIXct
# # df <- df %>%
# #   mutate(time = ymd_hms(time))
# #
# # # Calculate time differences
# # df <- df %>%
# #   ##order by time
# #   arrange(time) %>%
# #   ##calculate time since last observation and fill down to next
# #   mutate(
# #     last_deer_time = case_when(sp == "deer" ~ time, TRUE ~ NA_POSIXct_),
# #     last_hiker_time = case_when(sp == "hiker" ~ time, TRUE ~ NA_POSIXct_),
# #     last_bike_time = case_when(sp == "biker" ~ time, TRUE ~ NA_POSIXct_),
# #     last_atv_time = case_when(sp == "atv" ~ time, TRUE ~ NA_POSIXct_)
# #   ) %>%
# #   fill(last_deer_time, .direction = "down") %>%
# #   fill(last_hiker_time, .direction = "down") %>%
# #   fill(last_bike_time, .direction = "down") %>%
# #   fill(last_atv_time, .direction = "down") %>%
# #   ##Deer count to eventually remove subsequent detections without rec between
# #   mutate(deer.count = if_else(sp == "deer", 1, NA_integer_)) %>%
# #   group_by(grp = cumsum(lag(sp, default = first(sp)) != "deer")) %>%
# #   mutate(deer.count = if_else(sp == "deer", row_number(), NA_integer_)) %>%
# #   ungroup() %>%
# #   select(-grp)%>%
# #   ##bike count to eventually remove subsequent detections without deer between
# #   mutate(hiker.count = if_else(sp == "hiker", 1, NA_integer_)) %>%
# #   group_by(grp = cumsum(lag(sp, default = first(sp)) != "hiker")) %>%
# #   mutate(hiker.count = if_else(sp == "hiker", row_number(), NA_integer_)) %>%
# #   ungroup() %>%
# #   select(-grp)%>%
# #   ##calculate HW and WH
# #   mutate(
# #     expected_HW2 = case_when(sp == "deer" & #species is deer
# #                                deer.count==1 & # first count of deer in back to back deer sequences
# #                                !is.na(last_hiker_time) & #not before the first hiker is detected
# #                                (last_hiker_time > last_bike_time | is.na(last_bike_time)) & # no bikes between
# #                                (last_hiker_time > last_atv_time | is.na(last_atv_time))  #no atvs between
# #                                ~ as.numeric(difftime(time, last_hiker_time, units = "hours")), TRUE ~ NA_real_), #calculate time difference in hours.
# #     expected_WH2 = case_when(sp == "hiker" &
# #                                hiker.count==1 &
# #                                !is.na(last_deer_time) &
# #                                (last_deer_time > last_bike_time | is.na(last_bike_time)) &
# #                                (last_deer_time > last_atv_time | is.na(last_atv_time))
# #                              ~ as.numeric(difftime(time, last_deer_time, units = "hours")), TRUE ~ NA_real_)
# #   )%>%
# #   select(-c(last_deer_time:last_atv_time, deer.count:hiker.count))
# #
# # # Print the resulting data frame
# # print(df)
# 
# ## map cams
# library(janitor)
# cam.locs <- read_csv("data/Camera Deployment/Camera Deployment.csv") %>%
#   clean_names() %>%
#   st_as_sf(coords = c("location_longitude", "location_latitude"), crs = 4326)
# 
# 
# cam.locs %>%
#   mutate(Name = camera_name) %>% ## get labels on kml
#   st_write(here::here("output", "spatial", "cams", "rec_cams_all.kml"), driver = "KML", delete_dsn = TRUE)
# 
# cam.locs %>%
#   filter(is.na(date_removed)) %>%
#   mutate(Name = camera_name) %>% ## get labels on kml
#   st_write(here::here("output", "spatial", "cams", "rec_cams_active.kml"), driver = "KML", delete_dsn = TRUE)
# 
# cam.locs %>%
#   filter(surface == "Road") %>%
#   mutate(Name = camera_name) %>% ## get labels on kml
#   st_write(here::here("output", "spatial", "cams", "road.kml"), driver = "KML", delete_dsn = TRUE)
# 
# cam.locs %>%
#   filter(surface == "Rec trail") %>%
#   mutate(Name = camera_name) %>% ## get labels on kml
#   st_write(here::here("output", "spatial", "cams", "singletrack.kml"), driver = "KML", delete_dsn = TRUE)
# 
# cam.locs %>%
#   filter(surface == "Wildlife trail") %>%
#   mutate(Name = camera_name) %>% ## get labels on kml
#   st_write(here::here("output", "spatial", "cams", "wildlife.kml"), driver = "KML", delete_dsn = TRUE)
# 
# 
# cam.locs %>%
#   filter(area == "FerniePP") %>%
#   st_transform(26911) %>%
#   st_bbox() %>%
#   st_as_sfc() %>%
#   st_buffer(5000) %>%
#   st_as_sf() %>%
#   # mapview()%>%
#   st_write(here::here("output", "spatial", "sa.shp"), delete_dsn = TRUE)
# 
# 
# 
# 
# 
# 
# sr.cam.locs <- sr.cam.locs %>%
#   clean_names() %>%
#   st_as_sf(coords = c("location_longitude", "location_latitude"), crs = 4326)



##check conflict between cities in EV
conf <- read_csv("/Users/claytonlamb/Dropbox/Documents/University/Work/EV_coexistence/EV_Coex_Analyses/data/conflict/encounter-390896-103283.csv")%>%
  mutate(year=year(encounter_date))
conf$species_name%>%unique()
conf$enctype_name%>%unique()
conf$encounter_locality%>%unique()

conf.ev <- conf%>%filter(encounter_locality%in%str_to_upper(c("Fernie", "Elko", "Hosmer", "Sparwood", "Elkford")))
conf.ev$enctype_name%>%unique()
conf.ev%>%
  group_by(species_name, year,encounter_locality)%>%
  count()%>%
  ggplot(aes(x=year, y=n, color=species_name))+
  geom_line()+
  facet_wrap(vars(encounter_locality), scales="free_y")

conf.ev%>%
group_by(species_name, encounter_locality)%>%
  count()%>%
  ggplot(aes(x=encounter_locality, y=n, fill=species_name))+
  geom_col(position = position_dodge(width=0.7))

conf.ev%>%
  filter(enctype_name%in%c("FOOD CONDITIONED", "AGGRESSIVE", "LIVESTOCK/PETS - KILLED/INJURED", "DAMAGE TO PROPERTY", "HUMAN INJURY/DEATH"))%>%
  group_by(species_name, encounter_locality)%>%
  count()%>%
  ggplot(aes(x=encounter_locality, y=n, fill=species_name))+
  geom_col(position = position_dodge(width=0.7))





##find large files for github rep and add to .gitignore

# List all files recursively
files <- list.files(here::here(), recursive = TRUE, full.names = TRUE)

# Get file info
info <- file.info(files)

# Add file names for clarity
info$filename <- rownames(info)

# Filter >100 MB (100*1024^2 bytes)
big_files <- subset(info, size > 100 * 1024^2)

# Show results
big_files[, c("filename", "size")]






