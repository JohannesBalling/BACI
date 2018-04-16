library(raster)
library(sp)
library(rgdal)
rm(list=ls())
gc()

####################################################################################################
# This script provides:
# 1. Method to extract all detected extreme events out of a raster stack
# 2. Area of each detected extreme event
# 3. Sampling of the extreme events used for the validation
#
#
#
########
# Parameters that have to be defined:
# Timefile - timefile of the dates of the CABLAB data (provided in Github repository)
timefile = "Link to this file as string"
# Classes - Raster stack of the BACIndex detected extreme events
Classes = "Link to the Classes stack detected by the BACIndex"
# Outpath - Link to the folder in which the shapefiles of the different extreme events have to be stored
outpath = "Link to results folder as string"
#####################################################################################################



timefile = read.table(file= timefile,sep = ";")
timefile$x = as.Date(timefile$x, format("%Y-%m-%d"))
Classes = brick(Classes)
area_ls=list() 
id_ls=list()
time_ls=list()
temp_time = as.numeric(seq(from=1, to =9973))


for (i in 1:nlayers(Classes)){
  print(paste(round(((i/nlayers(Classes))*100), digits=2),"%"))
  shape <- rasterToPolygons(Classes[[i]],fun=function(x){x == 3},dissolve=TRUE)
  if (maxValue(Classes[[i]] == 3)){
    shape2 <- disaggregate(shape)
    shape2$id <- factor(seq_len(length(shape2)))
    shape2$IDSel = temp_time[1:length(shape2$id)]
    temp_time = temp_time[-c(1:length(shape2$id))]
    writeOGR(obj=shape2, dsn=outpath,layer=paste('extreme_events',timefile[i,1],sep = "_"),
             driver="ESRI Shapefile",overwrite_layer=T)
    plot(Classes[[i]], main=paste(timefile[i,1]))
    plot(shape2,add=T)
    area_temp = raster::area(shape2)
    for (x in 1:length(shape2@polygons)){
      area_ls[[(length(area_ls)+1)]]= area_temp[x]
      id_ls[[(length(id_ls)+1)]] = x
      time_ls[[(length(time_ls)+1)]] = timefile[[i,1]]
    }
  }
  else{
    print("No extreme events")
  }
}


# Creating single Dataframe out of the three lists
df_area_ls = do.call(c, area_ls)
df_time_ls = do.call(c, time_ls)
df_id_ls = do.call(c, id_ls)

df_extreme_events = as.data.frame(matrix(NA, ncol = 4, nrow = length(df_area_ls)))
colnames(df_extreme_events)<-c('ID','Date','Area','IDSel')
df_extreme_events$ID = df_id_ls
df_extreme_events$Date = df_time_ls
df_extreme_events$Area = df_area_ls
df_extreme_events$IDSel = row.names(df_extreme_events)
# Sort dataframe based on areas of shapefiles --> select 20 biggest
area_sorted <- df_extreme_events[order(df_extreme_events$Area, decreasing = T),]
largest_events_20 = area_sorted[c(1:20),]
# Select 50 randomly
area_sorted_random = area_sorted[-c(1:20),]
random_samples = sample_n(area_sorted_random, 50)

