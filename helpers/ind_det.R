wt_ind_detect_custom <- function(x, threshold, units = "minutes", datetime_col = image_date_time, remove_human = TRUE, remove_domestic = TRUE) {
  
  # Check that x is a dataframe
  if (!is.data.frame(x)) {
    stop("The first argument must supply a dataframe.")
  }
  
  # Ensure that datetime_col is of class POSIXct; if not, try to convert.
  name <- enquo(datetime_col) |> quo_name()
  if (!inherits(x[[name]], c("POSIXct"))) {
    x <- x |> mutate({{datetime_col}} := as.POSIXct({{datetime_col}}, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))
    message("Your datetime_col has been converted to a POSIXct DateTime.")
  }
  # Check if x contains the required columns - standard output from WildTrax. Probably should make this more flexible.
  req_cols <- c("project_id", "location", "species_common_name", "individual_count")
  if (!all(req_cols %in% colnames(x))) {
    stop("Important columns are missing from the data you have supplied. All of `project_id`, `location`, `species_common_name`, and `individual_count` are required.")
  }
  
  # Check that the units argument is either seconds, minutes, or hours
  if (!units %in% c("seconds", "minutes", "hours")) {
    stop("Please use 'seconds', 'minutes', or 'hours' as your threshold unit.")
  }
  
  # Convert threshold to seconds
  if (units == "minutes") {
    threshold <- threshold  * 60
  } else if (units == "hours") {
    threshold <- threshold * 60 * 60
  } else {
    threshold
  }
  
  # Tags to discard - not sure if UNKNOWN is ever wanted?
  t <- c("NONE", "STAFF/SETUP", "UNKNOWN")
  if (remove_human) {
    # Standard WildTrax tags that refer to human(ish) objects
    t <- c(t, "Human", "Vehicle", "Unknown Vehicle", "All Terrain Vehicle", "Train", "Heavy Equipment")
  }
  x <- filter(x, !species_common_name %in% t)
  if (remove_domestic) {
    # All tags in WildTrax that refer to domestic animals begin with 'Domestic __'
    x <- filter(x, !grepl("^Domestic", species_common_name))
  }
  
  # Create ordered dataframe, and calculate time interval between images.
  x1 <- x |>
    # Sometimes VNA sneaks in here
    mutate(individual_count = as.numeric(ifelse(individual_count == "VNA", 1, individual_count))) |>
    # Amalgamate tags of same species in same image; currently broken into two separate rows
    # Use image_id because sometimes {{datetime_col}} is same for 2 images
    group_by(location, image_id, species_common_name) |>
    mutate(individual_count = sum(individual_count)) |>
    distinct(location, {{datetime_col}}, species_common_name, individual_count, .keep_all = TRUE) |>
    ungroup() |>
    # Order the dataframe
    arrange(project_id, location, {{datetime_col}}, species_common_name) |>
    group_by(project_id, location, species_common_name) |>
    # Calculate the time difference between subsequent images
    mutate(interval = as.numeric(difftime({{datetime_col}}, lag({{datetime_col}}), units = "secs"))) |>
    # Is this considered a new detection?
    mutate(new_detection = ifelse(is.na(interval) | abs(interval) >= threshold, TRUE, FALSE)) |>
    # Number independent detections
    mutate(detection = paste(species_common_name,c(1, cumsum(new_detection[-1]) + 1)))|>
    ungroup() 
  
  
  # Summarise detections
  x2 <- x1 |>
    group_by(detection, project_id, location, species_common_name) |>
    summarise(start_time = min({{datetime_col}}),
              end_time = max({{datetime_col}}),
              total_duration_seconds = as.numeric(difftime(end_time, start_time, units = "secs")),
              n_images = n(),
              avg_animals_per_image = mean(individual_count),
              max_animals = max(individual_count)) |>
    ungroup()
  
  # Return x2
  return(x2)
  
}