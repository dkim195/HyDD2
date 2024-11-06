library(sp)
library(imputeTS) 
library(ggplot2)
library(ggmap)
library(data.table)
library(cowplot)
library(corrplot)
library(keras)
library(future)
library(presto)
set.seed(1234)
plan("multicore", workers = 10)
options(future.globals.maxSize = 100000 * 1024^2) #100 gb ram
setwd("/media/thomaskim/Data/")
use_python("/home/thomaskim/miniconda3/envs/deeplabcut/bin/python", required = TRUE)
source('Figures/DeepLabCut/DLCAnalyzer_Functions_final.R')

####Start####
#load
Tracking <- ReadDLCDataFromCSV(file = "Figures/DeepLabCut/IMG_1531_2DLC_Resnet50_DlxJul25shuffle1_snapshot_200_filtered.csv", fps = 60) #Ctrl
Tracking <- ReadDLCDataFromCSV(file = "Figures/DeepLabCut/IMG_1530_2DLC_Resnet50_DlxJul25shuffle1_snapshot_200_filtered.csv", fps = 60) #Mut

Tracking <- ReadDLCDataFromCSV(file = "Figures/DeepLabCut/Control_P21_Male_Litter5_Foxd1Cre(het)_Dlx1flox(het)_Dlx2flox(het)DLC_Resnet50_DlxJul25shuffle1_snapshot_200_filtered.csv", fps = 60)
Tracking <- ReadDLCDataFromCSV(file = "Figures/DeepLabCut/Control_P22_Male_Litter4_Foxd1Cre(het)_Dlx1flox(het)_Dlx2flox(het)DLC_Resnet50_DlxJul25shuffle1_snapshot_200_filtered.csv", fps = 60)
Tracking <- ReadDLCDataFromCSV(file = "Figures/DeepLabCut/Control_P23_Male_Litter3_Foxd1Cre(het)_Dlx1flox(het)_Dlx2flox(het)DLC_Resnet50_DlxJul25shuffle1_snapshot_200_filtered.csv", fps = 60)
Tracking <- ReadDLCDataFromCSV(file = "Figures/DeepLabCut/Control_P33_Male_Litter1_Foxd1Cre(het)_Dlx1flox(het)_Dlx2flox(het)-001DLC_Resnet50_DlxJul25shuffle1_snapshot_200_filtered.csv", fps = 60)

Tracking <- ReadDLCDataFromCSV(file = "Figures/DeepLabCut/DlxKO_P21_Male_Litter5_Foxd1Cre(het)_Dlx1flox(homo)_Dlx2flox(homo)DLC_Resnet50_DlxJul25shuffle1_snapshot_200_filtered.csv", fps = 60)
Tracking <- ReadDLCDataFromCSV(file = "Figures/DeepLabCut/DlxKO_P22_Male_Litter4_Foxd1Cre(het)_Dlx1flox(homo)_Dlx2flox(homo)DLC_Resnet50_DlxJul25shuffle1_snapshot_200_filtered.csv", fps = 60)
Tracking <- ReadDLCDataFromCSV(file = "Figures/DeepLabCut/DlxKO_P23_Male_Litter3_Foxd1Cre(het)_Dlx1flox(homo)_Dlx2flox(homo)DLC_Resnet50_DlxJul25shuffle1_snapshot_200_filtered.csv", fps = 60)
Tracking <- ReadDLCDataFromCSV(file = "Figures/DeepLabCut/DlxKO_P33_Male_Litter1_Foxd1Cre(het)_Dlx1flox(homo)_Dlx2flox(homo)-002DLC_Resnet50_DlxJul25shuffle1_snapshot_200_filtered.csv", fps = 60)

names(Tracking$data)

PlotPointData <- function(t, points = NULL, from = NULL, to = NULL, unit = "frame", type = NULL){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  
  if(is.null(points) & is.null(type)){
    points <- names(t$data)
  }
  
  if(unit == "second"){
    if(!is.null(from)){
      from <- t$frames[which(t$seconds >= from)[1]]
    }
    if(!is.null(to)){
      to <- t$frames[which(t$seconds >= to)[1]]
    }
  }
  
  if(is.null(from)){
    from <- min(t$frames)
  }
  if(is.null(to)){
    to <- max(t$frames)
  }
  
  if(!is.null(type)){
    points <- t$point.info[t$point.info$PointType == type, "PointName"]
  }
  
  range <- from:to
  
  p <- ggdraw()
  dim <- ceiling(sqrt(length(points)))
  nplot <- 0
  
  for (i in points) {
    p <- p + draw_plot(
      ggplot(data = t$data[[i]][t$data[[i]]$frame %in% range,], aes(x = x, y = y)) +
        geom_path(color = "blue") +  # Set a fixed color for the path
        ggtitle(i) +
        xlab(paste("x /", t$distance.units, sep = " ")) +
        ylab(paste("y /", t$distance.units, sep = " ")) +
        theme_bw(),
      x = (nplot %% dim / dim),
      y = ((dim - 1) / dim) - floor(nplot / dim) / dim,
      width = 1 / dim,
      height = 1 / dim
    )
    nplot <- nplot + 1
  }
  
  return(p)
}


#PlotPointData(Tracking, points = names(Tracking$data))
PlotPointData(Tracking, points ="TailStart" )

####RUN####
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
#plot_zones(Tracking$zones)

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

OFTAnalysis <- function(t, movement_cutoff, integration_period, points){
  if (!IsTrackingData(t)) {
    stop("Object is not of type TrackingData")
  }
  
  if (is.null(t$zones)) {
    warning("No zones defined for OFT analysis. Returning simple analysis only")
  }
  
  # Calculate movement
  t <- CalculateMovement(t, movement_cutoff, integration_period)
  
  # Initialize report list
  t$Report <- list()
  if (!is.null(t$labels)) {
    t$Report <- append(t$Report, LabelReport(t, integration_period))
  }
  
  # Iterate over each point to generate reports
  for (k in points) {
    # Check if point data exists
    if (is.null(t$data[[k]])) {
      stop(paste("Point", k, "not found in Tracking data"))
    }
    
    dat <- t$data[[k]]
    
    # Basic movement data
    t$Report[[paste(k, "raw.distance", sep = ".")]] <- sum(dat[,"speed"], na.rm = TRUE)
    t$Report[[paste(k, "distance.moving", sep = ".")]] <- sum(dat[dat$is.moving, "speed"], na.rm = TRUE)
    t$Report[[paste(k, "raw.speed", sep = ".")]] <- mean(dat[,"speed"], na.rm = TRUE) * t$fps
    t$Report[[paste(k, "speed.moving", sep = ".")]] <- mean(dat[dat$is.moving, "speed"], na.rm = TRUE) * t$fps
    t$Report[[paste(k, "time.moving", sep = ".")]] <- sum(dat[,"is.moving"], na.rm = TRUE) / t$fps
    t$Report[[paste(k, "total.time", sep = ".")]] <- length(dat[,"is.moving"]) / t$fps
    t$Report[[paste(k, "time.stationary", sep = ".")]] <- t$Report[[paste(k, "total.time", sep = ".")]] - t$Report[[paste(k, "time.moving", sep = ".")]]
    t$Report[[paste(k, "percentage.moving", sep = ".")]] <- t$Report[[paste(k, "time.moving", sep = ".")]] / t$Report[[paste(k, "total.time", sep = ".")]] * 100
    
    # Zone-specific reports (if zones are defined)
    if (!is.null(t$zones)) {
      t$Report <- append(t$Report, ZoneReport(t, k, "topLeft", zone.name = paste(k, "topLeft", sep = ".")))
      t$Report <- append(t$Report, ZoneReport(t, k, "topRight", zone.name = paste(k, "topRight", sep = ".")))
      t$Report <- append(t$Report, ZoneReport(t, k, "bottomLeft", zone.name = paste(k, "bottomLeft", sep = ".")))
      t$Report <- append(t$Report, ZoneReport(t, k, "bottomRight", zone.name = paste(k, "bottomRight", sep = ".")))
    }
  }
  
  return(t)
}

####Process####
#calibrate
Tracking <-CalibrateTrackingData(Tracking, method = "ratio", ratio = 5) #in this case the px to cm (px/cm) ratio would be set to 5 and the data calibrated with this ratio
Tracking$px.to.cm

#Quantify how much/fast animals move and how much time it spends
#PlotZones(Tracking)
#PlotZoneVisits(Tracking, point = c("Nose", "LeftEar",
#                                   "RightEar", "Back1",
#                                   "Back2",
#                                   "TailStart", "TailEnd"))

Tracking <- CalculateMovement(Tracking, movement_cutoff = 5, integration_period = 5)
#head(Tracking$data$Nose)
#plots <- PlotDensityPaths(Tracking,points = c("Nose","Back1","TailEnd"))
#plots$TailEnd

#Reports #go through "RUN"
#Report <- ZoneReport(Tracking, point = "TailStart", zones = "topLeft")
#t(data.frame(Report))
#PlotZoneSelection(Tracking, point = "Back1", zones = "topLeft", invert = TRUE)

#analysis
Tracking <- OFTAnalysis(Tracking, point = "Nose", movement_cutoff = 5, integration_period = 5)
#t(data.frame(Tracking$Report))

# Check the structure and contents of the report
#print(Tracking$Report)

# If the report is a list or a similar structure, you might access it like this:
total_distance <- Tracking$Report[["Nose.raw.distance"]]
average_speed <- Tracking$Report[["Nose.raw.speed"]]
time_moving <- Tracking$Report[["Nose.time.moving"]]

# Print or use the metrics
print(paste("Total Distance Moved:", total_distance, "cm"))
print(paste("Average Speed:", average_speed, "cm/s"))
print(paste("Time Spent Moving:", time_moving, "seconds"))

# Analyze or visualize
# Example: Check if time moving is above a threshold
#if (time_moving > 300) {
#  print("Subject was moving for more than 5 minutes.")
#}

#ggplot(data = Tracking$data$Nose, aes(x = frame, y = speed)) +
#  geom_line() +
#  labs(title = "Speed Over Time", x = "Frame", y = "Speed (cm/s)") +
#  theme_minimal()


