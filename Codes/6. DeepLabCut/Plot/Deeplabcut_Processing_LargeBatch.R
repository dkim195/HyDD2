library(data.table)
library(dplyr)
library(ggplot2)
library(sp)
library(future)
library(imputeTS)
library(cowplot)
library(corrplot)
library(keras)

# Set environment and other settings
set.seed(1234)
plan("multicore", workers = 10)
options(future.globals.maxSize = 100000 * 1024^2)  # 100 GB RAM
setwd("/media/thomaskim/Data/")
use_python("/home/thomaskim/miniconda3/envs/deeplabcut/bin/python", required = TRUE)
source('Figures/DeepLabCut/DLCAnalyzer_Functions_final.R')

####Run####
# Calculate the approximate center
calculate_center <- function(data) {
  x_center <- mean(range(unlist(lapply(data, function(d) d$x))), na.rm = TRUE)
  y_center <- mean(range(unlist(lapply(data, function(d) d$y))), na.rm = TRUE)
  list(x = x_center, y = y_center)
}

# Calculate the width and height
calculate_dimensions <- function(data) {
  x_range <- range(unlist(lapply(data, function(d) d$x)), na.rm = TRUE)
  y_range <- range(unlist(lapply(data, function(d) d$y)), na.rm = TRUE)
  list(width = x_range[2] - x_range[1], height = y_range[2] - y_range[1])
}

center <- calculate_center(Tracking$data)
dimensions <- calculate_dimensions(Tracking$data)

# Function to create quadrants
create_quadrants <- function(center, dimensions) {
  half_width <- dimensions$width / 2
  half_height <- dimensions$height / 2
  list(
    topLeft = data.frame(x = c(center$x - half_width, center$x, center$x, center$x - half_width),
                         y = c(center$y, center$y, center$y + half_height, center$y + half_height)),
    topRight = data.frame(x = c(center$x, center$x + half_width, center$x + half_width, center$x),
                          y = c(center$y, center$y, center$y + half_height, center$y + half_height)),
    bottomLeft = data.frame(x = c(center$x - half_width, center$x, center$x, center$x - half_width),
                            y = c(center$y, center$y, center$y - half_height, center$y - half_height)),
    bottomRight = data.frame(x = c(center$x, center$x + half_width, center$x + half_width, center$x),
                             y = c(center$y, center$y, center$y - half_height, center$y - half_height))
  )
}

# Apply to Tracking data
Tracking$zones <- create_quadrants(center, dimensions)

library(ggplot2)

# Function to plot zones
plot_zones <- function(zones) {
  ggplot() +
    geom_polygon(data = zones$topLeft, aes(x = x, y = y), fill = "red", alpha = 0.5) +
    geom_polygon(data = zones$topRight, aes(x = x, y = y), fill = "blue", alpha = 0.5) +
    geom_polygon(data = zones$bottomLeft, aes(x = x, y = y), fill = "green", alpha = 0.5) +
    geom_polygon(data = zones$bottomRight, aes(x = x, y = y), fill = "yellow", alpha = 0.5) +
    theme_minimal()
}

# Plot zones
plot_zones(Tracking$zones)

PlotZoneVisits <- function(t, points, zones = NULL, invert = FALSE) {  # Added default value for invert
  if (!IsTrackingData(t)) {
    stop("Object is not of type TrackingData")
  }
  if (is.null(t$zones)) {
    warning("No zones defined")
    return(NULL)
  }
  if (is.null(zones)) {
    zones <- names(t$zones)
  }
  if (length(setdiff(zones, names(t$zones))) > 0) {
    warning("Invalid zones")
    return(NULL)
  }
  if (length(setdiff(points, names(t$data))) > 0) {
    warning("Invalid points")
    return(NULL)
  }
  
  dat <- NULL
  for (j in points) {
    for (i in zones) {
      # Assume zone visitation logic here, where IsInZone(t, j, i, invert) is defined
      dat <- rbind(dat, data.frame(seconds = t$seconds, zone = ifelse(IsInZone(t, j, i, invert), i, NA), type = "automatic", points = j))
    }
  }
  
  ggplot(data = na.omit(dat), aes(seconds, zone, color = zone)) +
    geom_point(size = 2, shape = 124) +
    facet_grid(points ~ .) +
    theme_bw()
}

IsInZone <- function(t, point, zone, invert = FALSE) {
  # Assuming point data and zone data are both data.frames with x and y columns
  point_data <- t$data[[point]]
  zone_data <- t$zones[[zone]]
  
  # Check if each point is within the polygon defined by zone_data
  library(sp)
  point_sp <- SpatialPoints(coords = matrix(c(point_data$x, point_data$y), ncol = 2))
  zone_polygon <- Polygon(matrix(c(zone_data$x, zone_data$y), ncol = 2, byrow = TRUE))
  zone_sp <- SpatialPolygons(list(Polygons(list(zone_polygon), "zone")))
  
  inside <- point.in.polygon(point_data$x, point_data$y, zone_data$x, zone_data$y) > 0
  
  if (invert) {
    return(!inside)
  } else {
    return(inside)
  }
}

PlotZoneVisits <- function(t, points, zones = NULL, invert = FALSE) {
  if (!IsTrackingData(t)) {
    stop("Object is not of type TrackingData")
  }
  if (is.null(t$zones)) {
    warning("No zones defined")
    return(NULL)
  }
  if (is.null(zones)) {
    zones <- names(t$zones)
  }
  if (length(setdiff(zones, names(t$zones))) > 0) {
    warning("Invalid zones")
    return(NULL)
  }
  if (length(setdiff(points, names(t$data))) > 0) {
    warning("Invalid points")
    return(NULL)
  }
  
  dat <- NULL
  for (point in points) {  # Adjust this line to loop over each point
    for (zone in zones) {
      in_zone <- IsInZone(t, point, zone, invert)
      if (length(in_zone) == 0) {
        warning(paste("No valid data for zone:", zone))
        next
      }
      dat <- rbind(dat, data.frame(seconds = t$data[[point]]$frame, zone = ifelse(in_zone, zone, NA), point = point))
    }
  }
  
  if (is.null(dat) || nrow(dat) == 0) {
    warning("No data available to plot after applying zone filters")
    return(NULL)
  }
  
  dat <- na.omit(dat)
  if (nrow(dat) == 0) {
    warning("All data removed by NA omission")
    return(NULL)
  }
  
  ggplot(dat, aes(x = seconds, y = zone, color = zone)) +
    geom_point(size = 2, shape = 124) +
    facet_wrap(~point) +  # Use facet_wrap to handle multiple points
    theme_bw()
}


# Load and preprocess data for each mouse
load_and_process <- function(file_path) {
  Tracking <- ReadDLCDataFromCSV(file = file_path, fps = 60)
  Tracking <-CalibrateTrackingData(Tracking, method = "ratio", ratio = 2) #in this case the px to cm (px/cm) ratio would be set to 5 and the data calibrated with this ratio
  Tracking <- CalculateMovement(Tracking, movement_cutoff = 1, integration_period = 1)
  Tracking <- OFTAnalysis(Tracking, point = "Nose", movement_cutoff = 1, integration_period = 1)
  
  # Extract and return necessary metrics from Tracking data
  return(list(
    total_distance = sum(Tracking$Report$Nose.distance.moving, na.rm = TRUE),
    average_speed = mean(Tracking$data$Nose$speed, na.rm = TRUE)
  ))
}

# List CSV files in the directory
input_folder <- "Figures/DeepLabCut/"
files <- list.files(input_folder, full.names = TRUE, pattern = "*.csv")
files

# Apply the function to each file
results <- lapply(files, load_and_process)

# Combine all results into a data frame
results_df <- do.call(rbind.data.frame, results)
row.names(results_df) <- gsub(pattern = ".*/|\\.csv", replacement = "", x = files)

# Print results
print(results_df)

#names
results_df$mouse_id <- rownames(results_df)
name_mapping <- data.frame(
  old_names = c("Control_P21_Male_Litter5_Foxd1Cre(het)_Dlx1flox(het)_Dlx2flox(het)DLC_Resnet50_DlxJul25shuffle1_snapshot_200_filtered", 
                "Control_P22_Male_Litter4_Foxd1Cre(het)_Dlx1flox(het)_Dlx2flox(het)DLC_Resnet50_DlxJul25shuffle1_snapshot_200_filtered", 
                "Control_P23_Male_Litter3_Foxd1Cre(het)_Dlx1flox(het)_Dlx2flox(het)DLC_Resnet50_DlxJul25shuffle1_snapshot_200_filtered", 
                "Control_P33_Male_Litter1_Foxd1Cre(het)_Dlx1flox(het)_Dlx2flox(het)-001DLC_Resnet50_DlxJul25shuffle1_snapshot_200_filtered", 
                "DlxKO_P21_Male_Litter5_Foxd1Cre(het)_Dlx1flox(homo)_Dlx2flox(homo)DLC_Resnet50_DlxJul25shuffle1_snapshot_200_filtered", 
                "DlxKO_P22_Male_Litter4_Foxd1Cre(het)_Dlx1flox(homo)_Dlx2flox(homo)DLC_Resnet50_DlxJul25shuffle1_snapshot_200_filtered", 
                "DlxKO_P23_Male_Litter3_Foxd1Cre(het)_Dlx1flox(homo)_Dlx2flox(homo)DLC_Resnet50_DlxJul25shuffle1_snapshot_200_filtered", 
                "DlxKO_P33_Male_Litter1_Foxd1Cre(het)_Dlx1flox(homo)_Dlx2flox(homo)-002DLC_Resnet50_DlxJul25shuffle1_snapshot_200_filtered", 
                "IMG_1530_2DLC_Resnet50_DlxJul25shuffle1_snapshot_200_filtered", 
                "IMG_1531_2DLC_Resnet50_DlxJul25shuffle1_snapshot_200_filtered"),
  new_names = c("Ctrl1","Ctrl2","Ctl3","Ctrl4","Mut1","Mut2","Mut3","Mut4","Mut5","Ctrl5")
)
name_vector <- setNames(name_mapping$new_names, name_mapping$old_names)
# Replace the mouse_id names directly using the vector
results_df$mouse_id <- name_vector[results_df$mouse_id]
# Check the results
print(results_df)

group_assignments <- data.frame(
  mouse_id = c("Ctrl1", "Ctrl2", "Ctl3", "Ctrl4", "Mut1", "Mut2", "Mut3", 
               "Mut4", "Mut5", "Ctrl5"),
  group = c("Control", "Control", "Control", "Control", "Mutant", "Mutant", "Mutant",
            "Mutant", "Mutant", "Control")
)
# Merge group information with your main data
results_df <- merge(results_df, group_assignments, by = "mouse_id")
# Check the updated dataframe
print(results_df)

# ANOVA to test differences in total distance or speed
anova_results <- aov(total_distance ~ results_df$group, data = results_df)
summary(anova_results)

# Plotting results
ggplot(results_df, aes(x = group, y = average_speed, fill = group)) +
  geom_bar(stat = "identity") +
  labs(title = "Average speed Per Group", x = "Mouse ID", y = "Average speed") +
  theme_minimal()

ggplot(results_df, aes(x = group, y = total_distance, fill = group)) +
  geom_bar(stat = "identity") +
  labs(title = "Total_distance Per Group", x = "Mouse ID", y = "Total_distance") +
  theme_minimal()

library(dplyr)
# Adding ranks based on average speed
results_df$rank <- rank(-results_df$total_distance)  # Rank in descending order of speed

# Print results with ranks
results_df <- results_df %>%
  arrange(rank) %>%
  mutate(mouse_id = row.names(.))  # Include mouse IDs as a column for easier plotting

print(results_df)

# Plotting results with ranks
ggplot(results_df, aes(x = reorder(mouse_id, -total_distance), y = total_distance, fill = group)) +
  geom_bar(stat = "identity") +
  labs(title = "Ranking of Mice by Average Speed", x = "Mouse ID", y = "Average Speed") +
  theme_minimal() +
  geom_text(aes(label = paste("Rank:", rank)), vjust = -0.5)  # Add rank labels to the plot

# Save plot if needed
ggsave("ranked_speeds.png", width = 10, height = 8, dpi = 300)