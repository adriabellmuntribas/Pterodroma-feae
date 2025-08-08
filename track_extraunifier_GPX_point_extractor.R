# Authors: Militão, T., González-Solís, J., Bellmunt A.,
# Tittle: Predictive nesting habitat modelling of Pterodroma feae in the Cabo Verde Archipelago
# R code for developing habitat modelling based on Pterodroma feae nesting locations in Cabo Verde archipélago
# Last updated version: 18/07/2025
# For doubts and questions about the code: adria.bellmunt.ribas@gmail.com

# Code explanation: This Rstudio code allows us to extract and unify the geometries of each of 
# the .gpx documents in order to group them into a single .shp document for later processing. 
# It is very important to highlight two points. It is necessary to take into account the 
# projection 32626 and 32627 and it is also necessary to take into account that the fields 
# are not maintained, that is, the information from the attribute tables is lost. Also, inform 
# that the final result is a multilinstring. The points of the .gpx are not saved. Important 
# information: It filters for 5km/h walking speed.

# Required packages

if (!requireNamespace("librarian", quietly = TRUE)) {
  install.packages("librarian")
}
suppressPackageStartupMessages(library(librarian))

librarian::shelf(
  dplyr,
  tidyr,
  purrr,
  geosphere,
  lubridate,
  readr,
  readxl,
  stringr,
  data.table,
  car,
  tibble,
  corrplot,
  spdep,
  tidyverse,
  sf,
  terra,
  Hmisc,
  vapour,
  ggplot2,
  gam,
  mgcv
)

# Input/output folders and CRS
input_folder <- "/Users/adriabellmuntribas/Desktop/TFM/QGIS/Gis project/mask/Search_dog/africa_santoantao"
output_folder <- "/Users/adriabellmuntribas/Desktop/TFM/QGIS/Gis project/mask/mask_santoantao"
crs_target <- 32626  # UTM zone 26N

# List of GPX files
gpx_files <- list.files(input_folder, pattern = "\\.gpx$", full.names = TRUE, recursive = TRUE)

# Max speed in m/s (5 km/h)
max_speed_ms <- 5 * 1000 / 3600

# Time parser
parse_time_column <- function(time_col) {
  time_clean <- gsub(" [A-Z]+$", "", time_col)  # Remove timezone like "CET"
  return(ymd_hms(time_clean, tz = "UTC"))
}

# Function to segment valid track points into multiple lines
process_gpx <- function(path) {
  layer <- tryCatch(st_read(path, layer = "track_points", quiet = TRUE), error = function(e) NULL)
  if (is.null(layer)) return(NULL)
  if (!"POINT" %in% st_geometry_type(layer)) return(NULL)
  
  # Parse time and sort
  layer <- layer %>%
    mutate(time = parse_time_column(time)) %>%
    arrange(time)
  
  # Coordinates and time
  coords <- st_coordinates(layer)
  times <- layer$time
  
  # Compute distances and time diffs
  distances <- distHaversine(coords[-nrow(coords), ], coords[-1, ])
  time_diffs <- as.numeric(difftime(times[-1], times[-nrow(layer)], units = "secs"))
  speeds <- distances / time_diffs
  speeds <- c(NA, speeds)  # Align
  
  layer$speed <- speeds
  
  # Define segment ID: a new segment starts after a speed > max_speed_ms
  layer$segment_id <- cumsum(is.na(layer$speed) | layer$speed > max_speed_ms)
  
  # Group by segment and create a linestring only if >=2 points
  lines_list <- layer %>%
    group_by(segment_id) %>%
    filter(n() >= 2) %>%
    summarise(do_union = FALSE) %>%
    st_cast("LINESTRING") %>%
    st_transform(crs_target)
  
  return(lines_list)
}

# Process all GPX files
valid_lines <- map(gpx_files, process_gpx)
valid_lines <- valid_lines[!sapply(valid_lines, is.null)]

# Combine all valid line segments into a single sf object
if (length(valid_lines) > 0) {
  multilines <- do.call(rbind, valid_lines)
  multilines_sf <- st_sf(geometry = st_geometry(multilines))
  
  output_path <- file.path(output_folder, "tracks_cleaned_segmented_africa_santoantao.shp")
  st_write(multilines_sf, output_path, delete_layer = TRUE)
  
  message("Final shapefile successfully created at: ", output_path)
} else {
  message("No valid track segments were generated.")
}

