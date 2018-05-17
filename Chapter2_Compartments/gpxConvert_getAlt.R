library(plotKML)
setwd('./DataSources/GPX')

# identify gpx files that contain waypoints
files <- list.files('./Chapter2_Compartments/DataSources/GPX')[-1]

path <- './Chapter2_Compartments/DataSources/GPX/'

# readGPX('./Chapter2_Compartments/DataSources/GPX/Waypoints_08-Aug-16.gpx', tracks=FALSE, routes=FALSE)$waypoints[, c('name', 'lon', 'lat', 'time')]
# readGPX(paste(path,files[i], sep=""), tracks=FALSE, routes=FALSE)$waypoints[, c('name', 'lon', 'lat', 'time')]

# combine all waypoints into a dataframe
allwaypoints <- list()
for (i in 1:length(files)) {
  allwaypoints[[i]] <- readGPX(paste(path,files[i], sep=""), tracks=FALSE, routes=FALSE)$waypoints[, c('name', 'lon', 'lat', 'time', 'ele')]
}

allwaypoints <- do.call('rbind', allwaypoints)
names(allwaypoints)

range(na.omit(allwaypoints$ele))

