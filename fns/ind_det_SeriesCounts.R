wt_ind_detect_series <- function(x, threshold, units = "minutes", datetime_col = image_date_time,
                                 remove_human = TRUE, remove_domestic = TRUE) {
  if (!is.data.frame(x)) {
    stop("The first argument must supply a dataframe.")
  }
  name <- enquo(datetime_col) %>% quo_name()
  if (!inherits(x[[name]], c("POSIXct"))) {
    x <- x %>% mutate(`:=`(image_date_time, lubridate::as_datetime(image_date_time)))
    message("Your datetime_col has been converted to a Date.")
  }
  req_cols <- c(
    "project_id", "location", "species_common_name",
    "individual_count", "series_count"
  ) ## CTL added "series_count"
  if (!all(req_cols %in% colnames(x))) {
    stop("Important columns are missing from the data you have supplied. All of `project_id`, `location`, `species_common_name`, and `individual_count` are required.")
  }
  if (!units %in% c("seconds", "minutes", "hours")) {
    stop("Please use 'seconds', 'minutes', or 'hours' as your threshold unit.")
  }
  if (units == "minutes") {
    threshold <- threshold * 60
  } else if (units == "hours") {
    threshold <- threshold * 60 * 60
  } else {
    threshold
  }
  t <- c("NONE", "STAFF/SETUP", "UNKNOWN")
  if (remove_human) {
    t <- c(
      t, "Human", "Vehicle", "Unknown Vehicle", "All Terrain Vehicle",
      "Train", "Heavy Equipment"
    )
  }
  x <- filter(x, !species_common_name %in% t)
  if (remove_domestic) {
    x <- filter(x, !str_detect(species_common_name, "^Domestic"))
  }
  x1 <- x %>%
    mutate(individual_count = as.numeric(ifelse(individual_count ==
      "VNA", 1, individual_count))) %>%
    group_by(
      location,
      image_date_time, species_common_name
    ) %>%
    mutate(individual_count = case_when(!is.na(max(series_count))~max(series_count),
                                        TRUE~sum(individual_count))) %>%
    distinct(location, image_date_time, species_common_name, individual_count, .keep_all = TRUE) %>%
    ungroup() %>%
    arrange(project_id, location, image_date_time, species_common_name) %>%
    group_by(
      project_id, location,
      species_common_name
    ) %>%
    mutate(interval = int_length(image_date_time %--% lag(image_date_time))) %>%
    mutate(new_detection = ifelse(is.na(interval) |
      abs(interval) >= threshold, TRUE, FALSE)) %>%
    ungroup() %>%
    mutate(detection = c(1, cumsum(new_detection[-1]) + 1))
  x2 <- x1 %>%
    group_by(detection, project_id, location, species_common_name) %>%
    summarise(
      start_time = min(image_date_time), end_time = max(image_date_time), total_duration_seconds = int_length(start_time %--%
        end_time), n_images = n(), avg_animals_per_image = mean(individual_count),
      max_animals = max(individual_count)
    ) %>%
    ungroup()
  return(x2)
}