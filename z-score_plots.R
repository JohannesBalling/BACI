#Opening packages
library(raster)
library(rgdal)
library(sp)

#####################################################################################################
# This script creates histograms of the normal distribution of the extremne event 
# and a multi-year normal and their smoothing functions 
# Parameters have to be defined:
# 1. inpath - path of the CABLAB data stacks
inpath = "Path as a string"
# 2. shapefile - shapefile of the extreme event
shapefile = "Shapefile link as a string"
# 3. Classes - Raster stack of the Classes detected by the BACIndex
Classes = "Raster stack link as a string"
# 4. time - time of the extreme event as a string
time = "format YYYY-mm-dd"
# 5. outpath - folder in which the plots shall be saved
outpath = ""
# 6. timefile - timefile in repository 
timefile = "D:/balli001/000_Data/time/time.csv"
# 7. time_period - time period that should be taken into account aorund the extreme date
time_period = 2
# 8. Logical value, whether a normal curve should be drawn or not
Normal_Curve = T
# 9. Units as strings for all the CABLAB input variables
units_ta = c('1','2', '3', '4', '5')
#####################################################################################################





# Selecting Time Period out of Data
setwd(inpath)
shape <- readOGR(dsn= paste(shapefile))
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
highestDensity = c()
floor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level)
ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)
print("Preparing Class Stack")
Classes = brick(Classes)
Classes = Classes[[c(as.numeric(timefile$num))]]
Classes = crop(Classes, extent(shape), snap="near")
Classes[Classes >2] = NA

png(filename=paste(outpath,time,".png",sep = ""), width = 850, height = 500)
par(mfrow=c(2,3))
for(i in 1:length(temp_list)){
  print(paste("Analyzing raster stack", i, "of", length(temp_list)))
  temp = brick(temp_list[[i]])
  temp = crop(temp, extent(shape), snap="near")
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
  median_raster = overlay(temp, fun = mean)
  sd = cellStats(median_raster, "sd")
  median = cellStats(median_raster, "mean")
  min = floor_dec(cellStats(median_raster, "min"), level=0)
  max = ceiling_dec(cellStats(median_raster, "max"), level=0)
  z_score = round(((median_time - median) / sd), digits = 3)
  x = seq(min(min, min_time), max(max, max_time), length = 1000)
  # finding highest density to adjust ylim in histogram
  yhist <- hist(median_raster,breaks=20, plot=FALSE)
  xhist <- hist(temp1[[timestamp]], breaks=20, plot=FALSE)
  highestDensity <- max(xhist$density, yhist$density)
  # plotting histograms and due to if loop maybe dnorm curve

  hist(median_raster,col=rgb(0.1,0.1,0.1,0.4),xlim=c(min(min,min_time),max(max,max_time)),ylim=c(0,highestDensity),
       probability = TRUE,breaks=seq(from=min(min, min_time), to= max(max,max_time),by= (max(max,max_time)-min(min, min_time))/20), 
       xlab='',ylab="", lty="blank", 
       include.lowest=T, cex.axis=1.4,  yaxt='l',cex.lab=1.5, main=paste(gsub(".tif","",temp_list[[i]])), cex.main=2.05)
  # title(ylab=paste(gsub(".tif","",temp_list[[i]])), line=2.25, cex.lab=1.95,font.lab=1)
  if(Normal_Curve == T){
    mrF  = fitdistr(na.omit(c(getValues(median_raster))),"normal")
    curve(dnorm(x,mean=mrF$estimate[1],sd=mrF$estimate[2]),add=TRUE,col="black",lwd=2)
    mtF  = fitdistr(na.omit(c(getValues(temp1[[timestamp]]))),"normal")
    hist(temp1[[timestamp]],col=rgb(0.8,0.1,0.1,0.4), xlim=c(min(min,min_time),max(max,max_time)) , 
         ylim=c(0,highestDensity), probability = TRUE,
         breaks=seq(from=min(min, min_time), to= max(max,max_time), by= (max(max,max_time)-min(min, min_time))/20),
         main="", xlab="",ylab="Density", lty="blank",add=TRUE, include.lowest=T)
    curve(dnorm(x,mean=mtF$estimate[1],sd=mtF$estimate[2]),add=TRUE,col="red",lwd=2)
    z_score_list[[i]] = c(z_score, temp_list[[i]])}
  else{
    hist(temp1[[timestamp]],col=rgb(0.8,0.1,0.1,0.4), xlim=c(min(min,min_time), max(max,max_time)) ,
         ylim=c(0,highestDensity), probability = TRUE,breaks=20,main="", xlab="",ylab="Density", lty="blank",add=TRUE)
    z_score_list[[i]] = c(z_score, temp_list[[i]])}
  # mtext(paste("[Z-score: ",as.character(z_score),']', sep = ""), side = 3, line = 0, adj = 0.45, cex = 1.1)
  mtext(paste("[z-score: ",as.character(z_score),']', sep = ""), side = 1, line = 2.65, adj = 0.9, cex = 1.2, font = 1)
  mtext(paste(units_ta[i]), side = 1, line = 2.65, adj = 0.01, cex = 1.1)
}
plot.new()
legend("center", c("Extreme Event","Multi-Year Normal"), 
       col = c(col=rgb(0.8,0.1,0.1,0.6),rgb(0.1,0.1,0.1,0.6)), cex = 2, pch = 15, ncol = 1, bty='n')
dev.off()
z_score = array(unlist(z_score_list), dim=c(2,length(z_score_list)))
rm(list=setdiff(ls(), c("z_score_list","z_score")))






