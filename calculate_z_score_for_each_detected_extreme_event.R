#Opening packages
library(raster)
library(rgdal)
library(sp)
rm(list=ls())
gc()

########################################################################################################
# This script was developed for the H2020 BACI project in order to validate the develpoed BACIndex
# with the calculation of a z-score for all detected extreme events
#
# 
# Location of a table consiting of information on (all_events)
# 1. Time of the extreme event(time)
# 2. Id of the extreme event on one date (ID_selection)
# 3. New Id of the extreme event (ID_shape)
all_events = read.table(file = "Link to the Table as string", sep = ';',header = T)
# Inpath - Folder in which the stacks of the CABLAB data are stored
inpath = "Link to the folder as string"
# shapefile - Folder in which the shapefiles of the detected extreme events are stored
shapefile_link = "Link to the folder as string"
# Classes - Stack of the Classes detected by the BACIndex
Classes = "Link to the file as string"
# timefile - timefile provided on Github giving information of the acquisition dates of the CABLAB data
timefile = "Link to the file as string"
#######################################################################################################


# Creating lists for the results of the input parameters
z_score_GPP = list()
z_score_LE = list()
z_score_NEE = list()
z_score_SE = list()
z_score_TER = list()

for(p in 5:nrow(all_events)){
  # Parameters you have to define
  time = as.character(all_events[p,2])
  ID_shape = all_events[p,4]
  ID_selection = all_events[p,1]
  shapefile = paste(shapefile_link,"extreme_events_",time,".shp", sep = "")
  time_period = 2
  Normal_Curve = T

  ################################### Z-scores
  ########
  # Selecting Time Period out of Data
  setwd(inpath)
  shape <- readOGR(dsn= paste(shapefile))
  shape = shape[shape$IDSel==ID_shape,]
  time = c(as.Date(paste(time), format("%Y-%m-%d")))
  timefile = read.table(file= timefile,sep = ";")
  timefile$x = as.Date(timefile$x, format("%Y-%m-%d"))
  timestamp = c(which(timefile$x == time))
  temp_start = as.POSIXlt(as.Date(timefile[c(which(timefile$x == time)) - time_period,]))
  temp_end = as.POSIXlt(as.Date(timefile[c(which(timefile$x == time)) + time_period,]))
  years = unique(format(timefile,'%Y'))
  time_per = integer()
  class(time_per) = "Date"
  
  time_skip = seq(from = as.Date("2001-12-15", format("%Y-%m-%d")), 
                  to = as.Date("2002-01-15", format("%Y-%m-%d")), by = 'day')
  
  
  # Accounting for Dates close to new Years Eve
  if (format(as.Date(time, format = format("%Y-%m-%d")) , "%m-%d") %in% format(as.Date(time_skip, format = format("%Y-%m-%d")) , "%m-%d") == T){
    print("Date is close to New Years eve")
    years = years[c(-1,-11),]
    for(i in 1:length(years)){
      temp_start = as.POSIXlt(as.Date(timefile[c(which(timefile$x == time)) - time_period,]))
      temp_end = as.POSIXlt(as.Date(timefile[c(which(timefile$x == time)) + time_period,]))
      temp_start$year = c(as.numeric(years[i])  - 1900)
      temp_start = as.Date(temp_start, format("%Y-%m-%d"))
      temp_end$year = c(as.numeric(years[i])  - 1899)
      temp_end = as.Date(temp_end, format("%Y-%m-%d"))
      temp = seq(from = as.Date(temp_start), to = as.Date(temp_end), by = "day")
      time_per = c(time_per, temp)
    }
  }
  
  if (format(as.Date(time, format = format("%Y-%m-%d")) , "%m-%d") %in% format(as.Date(time_skip, format = format("%Y-%m-%d")) , "%m-%d") == F){
    print("Date is not close to New Years eve")
    for(i in 1:nrow(years)){
      temp_start = as.POSIXlt(as.Date(timefile[c(which(timefile$x == time)) - time_period,]))
      temp_end = as.POSIXlt(as.Date(timefile[c(which(timefile$x == time)) + time_period,]))
      temp_start$year = c(as.numeric(years[i,1])  - 1900)
      temp_start = as.Date(temp_start, format("%Y-%m-%d"))
      temp_end$year = c(as.numeric(years[i,1])  - 1900)
      temp_end = as.Date(temp_end, format("%Y-%m-%d"))
      temp = seq(from = as.Date(temp_start), to = as.Date(temp_end), by = "day")
      time_per = c(time_per, temp)
    }
  }
  
  
  timefile$num = row.names(timefile)
  timefile <- timefile[timefile$x %in% time_per, ]
  as.numeric(timefile$num)
  timefile = timefile[with(timefile, !((num >= (timestamp - time_period) & num <= (timestamp+ time_period)))), ]
  
  
  
  # clipping for the desired subset based on a shapefile 
  temp_list = list.files(pattern= ".tif$")
  z_score_list = list()
  mrF = list()
  mtF = list()
  yhist = c()
  xhist = c()
  temp_not_na = c()
  temp_exclude = c()
  highestDensity = c()
  floor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level)
  ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)
  print("Preparing Class Stack")
  Classes = brick(Classes)
  Classes = Classes[[c(as.numeric(timefile$num))]]
  cr <- crop(Classes, extent(shape), snap="out")                    
  fr <- rasterize(shape, cr)
  Classes <- mask(x=cr, mask=fr)
  Classes[Classes ==3] = NA
  
  for(i in 1:length(temp_list)){
    print(paste("Analyzing raster stack", i, "of", length(temp_list)))
    temp = brick(temp_list[[i]])
    cr <- crop(temp, extent(shape), snap="out")                    
    fr <- rasterize(shape, cr)   
    temp <- mask(x=cr, mask=fr)
    temp1 = temp
    # Creating plots and z-score
    median_time = cellStats(temp[[timestamp]], "mean")
    sd_time = cellStats(temp[[timestamp]], "sd")
    min_time = floor_dec(cellStats(temp[[timestamp]], "min"), level=0)
    max_time = ceiling_dec(cellStats(temp[[timestamp]], "max"), level=0)
    temp = temp[[c(as.numeric(timefile$num))]]
    # Setting Values in Raster Layers of Class 3 NA
    for (x in 1: nrow(timefile)){
      classes_x = Classes[[x]]
      temp_x    = temp[[x]]
      adjust_test <- overlay(temp_x, classes_x, fun = function(x, y) {
        x[is.na(y[])] <- NA
        return(x)
      })
      temp[[x]] = adjust_test
    }
    # Calculating z-score
    for (u in 1:nlayers(temp)){
      temp_not_na = c(temp_not_na,sum(!is.na(temp[[u]]@data@values)))
    }
    temp_not_na = max(temp_not_na)
    for (u in 1:nlayers(temp)){
      if (temp[[u]]@data@min == Inf | temp[[u]]@data@max == -Inf | 
          (sum(!is.na(temp[[u]]@data@values))/ temp_not_na) <=0.6){
        temp_exclude = c(temp_exclude,u)
      }
    }
    if (length(temp_exclude) >= 1){
      temp = temp[[-c(temp_exclude)]] 
    }
    median_raster = overlay(temp, fun = mean)
    sd = cellStats(median_raster, "sd")
    median = cellStats(median_raster, "mean")
    # test to skip errourus data
    if(median == 0 | is.nan(median) == T | is.na(median) == T)
    {
      #error handling code, maybe just skip this iteration using
      if(i == 1){
        z_score_GPP[[p]] = NA
      }
      if(i == 2){
        z_score_LE[[p]] = NA
      }
      if(i == 3){
        z_score_NEE[[p]] = NA
      }
      if(i == 4){
        z_score_SE[[p]] = NA
      }
      if(i == 5){
        z_score_TER[[p]] = NA
      }
      next
    }
    min = floor_dec(cellStats(median_raster, "min"), level=0)
    max = ceiling_dec(cellStats(median_raster, "max"), level=0)
    z_score = round(((median_time - median) / sd), digits = 3)
    if(i == 1){
      z_score_GPP[[p]] = z_score
    }
    if(i == 2){
      z_score_LE[[p]] = z_score
    }
    if(i == 3){
      z_score_NEE[[p]] = z_score
    }
    if(i == 4){
      z_score_SE[[p]] = z_score
    }
    if(i == 5){
      z_score_TER[[p]] = z_score
    }
  }
}
