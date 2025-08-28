#testing
detect_data = detections
raw_data = df
time_interval = "month"
variable = "all"
output_format = "long"
exclude_out_of_range = TRUE
effort_data = NULL

wt_summarise_cam_custom <- function(detect_data, raw_data, time_interval = "day",
                             variable = "detections", output_format = "wide",
                             effort_data = NULL,
                             exclude_out_of_range = FALSE,
                             start_col = start_time, end_col = end_time,
                             detection_id_col = detection,
                             start_col_det = start_time) {
  
  
  

  # Make sure one of raw_data or effort_data is supplied
  if (rlang::is_missing(raw_data) & is.null(effort_data)) {
    stop("Please supply a value for one of `raw_data` or `effort_data`.")
  }
  
  # Check that only one is supplied
  if (!rlang::is_missing(raw_data) & !is.null(effort_data)) {
    stop("Please only supply a value for one of `raw_data` or `effort_data`.")
  }
  
  if (variable == "all") {
    variable <- c("detections", "counts", "presence")
  }
  
  # Todo - check that supplied data contains all necessary columns
  
  # Check timeframe variable
  int <- c("day", "week", "month", "full")
  if (!time_interval %in% int) {
    stop("Please select a valid time interval: 'day', 'week', 'month', or 'full'.")
  }
  
  # Check output format
  formats <- c("wide", "long")
  if (!output_format %in% formats) {
    message("Please specify `wide` or `long` in the output_format argument.")
  }
  
  # Because station_col is also a function and arguments are lazily evaluated,
  # we need to deparse/substitute to a string
  #station_col <- deparse(substitute(station_col))
  
  # Parse the raw or effort data to get time ranges for each camera deployment.
  if (!rlang::is_missing(raw_data)) {

    if (exclude_out_of_range == FALSE) {
      x <- raw_data |>
        group_by(project_id, location, image_set_id) |>
        summarise(start_time = as.Date(min(image_date_time)),
                  end_time = as.Date(max(image_date_time))) |>
        ungroup()
      
      if (any(c(is.na(x$start_time), is.na(x$end_time)))) {
        message("Parsing of image date time produced NAs, these will be dropped")
        x <- drop_na(x)
      }
      
      # Expand the time ranges into individual days of operation (smallest unit)
      x <- x |>
        group_by(project_id, location,, image_set_id) |>
        mutate(day = list(seq.Date(start_time, end_time, by = "day"))) |>
        tidyr::unnest(day) |>
        mutate(year = as.integer(format(day, "%Y"))) |>
        select(project_id, location, year, day) |>
        ungroup()
      
    } else {
      x <- raw_data |>
        group_by(project_id, location, image_set_id) |>
        arrange(image_date_time) |>
        mutate(period = rep(seq_along(rle(image_fov)$lengths), rle(image_fov)$lengths)) |>
        filter(image_fov == "WITHIN") |>
        group_by(project_id, location, period, image_set_id) |>
        summarise(
          start_time = as.Date(min(image_date_time)),
          end_time = as.Date(max(image_date_time))
        ) |>
        ungroup()
      
      if (any(c(is.na(x$start_time), is.na(x$end_time)))) {
        message("Parsing of image date time produced NAs, these will be dropped")
        x <- drop_na(x)
      }
      
      # Expand the time ranges into individual days of operation (smallest unit)
      x <- x |>
        group_by(project_id, location, period, image_set_id) |>
        mutate(day = list(seq.Date(start_time, end_time, by = "day"))) |>
        tidyr::unnest(day) |>
        mutate(year = as.integer(format(day, "%Y"))) |>
        ungroup()%>%
        select(project_id, location, year, day)%>%
        distinct(location, day, .keep_all=TRUE)
    }

  }
  
  if(!is.null(effort_data)){
    x <- effort_data
  }
  
  # Based on the desired timeframe, assess when each detection occurred
  if (time_interval == "day" | time_interval == "full") {
    y <- detect_data |>
      mutate(year = as.integer(format(start_time, "%Y")),
             day = as.Date(start_time))
    grouping_cols <- c("year", "day")
  } else if (time_interval == "week") {
    y <- detect_data |>
      mutate(year = as.integer(format(start_time, "%Y")),
             week = as.integer(format(start_time, "%V")))  # ISO week
    grouping_cols <- c("year", "week")
  } else if (time_interval == "month") {
    y <- detect_data |>
      mutate(year = as.integer(format(start_time, "%Y")),
             month = format(start_time, "%B"))  # Full month name
    grouping_cols <- c("year", "month")
  }
  
  # Summarise variable of interest
  y <- y |>
    group_by(across(all_of(grouping_cols)), location, species_common_name) |>
    summarise(detections = n(),
              counts = sum(max_animals)) |>
    ungroup() |>
    mutate(presence = ifelse(detections > 0, 1, 0))%>%
    select(all_of(grouping_cols), location, species_common_name, detections:presence)
  
  # Species present in the data
  sp <- y |> select(species_common_name) |> distinct()
  
  # Create long df object of all species x location x timeframe combos
  if (time_interval == "day") {
    z <- x |>
      select(location, year, day)|>
      mutate(n_days_effort = 1) |>
      crossing(sp) |>
      left_join(y) |>
      mutate(across(6:8, ~ tidyr::replace_na(.x, 0)))
  } else if (time_interval == "week") {
    x <- x |>
      mutate(week = as.numeric(format(day, "%V"))) |>
      group_by(location, year, week) |>
      tally(name = "n_days_effort") |>
      ungroup()
    z <- x |>
      crossing(sp) |>
      left_join(y) |>
      mutate(across(6:8, ~ tidyr::replace_na(.x, 0)))
  } else if (time_interval == "month") {
    x <- x |>
      mutate(month = format(day, "%B")) |>
      group_by(location, year, month) |>
      tally(name = "n_days_effort") |>
      ungroup()
    z <- x |>
      crossing(sp) |>
      left_join(y) |>
      mutate(across(6:8, ~ tidyr::replace_na(.x, 0)))
  } else if (time_interval == "full") {
    z <- x |>
      crossing(sp) |>
      left_join(y) |>
      mutate(across(all_vars, ~ tidyr::replace_na(.x, 0))) |>
      group_by(location, year, species_common_name) |>
      summarise(detections = sum(detections),
                counts = sum(counts),
                presence = ifelse(any(presence == 1), 1, 0)) |>
      ungroup()
  }
  
  # Make wide if desired
  if (output_format == "wide") {
    z <- z |>
      tidyr::pivot_wider(id_cols = 1:5, names_from = species_common_name, values_from = {{ variable }}, names_sep = ".")
  } else if (output_format == "long") {
    z <- z |> dplyr::select(1:6, {{ variable }}) |>
      tidyr::pivot_longer(cols = {{ variable }}, names_to = "variable", values_to = "value")
  }
  
  return(z)
}
