#Opening packages
library(raster)
library(rgdal)
library(sp)

rm(list=ls())
gc()

########################################################################################################
# This script is computing a t-test based on yearly downstream data (raster stack) 
# and detected extreme event (shapefile) --> also boxplots are created!
# Parameters that have to be defined
# 1. all_together = table containing all the information on the extreme events
# everthing under 1. should be provided within the table itself
#     1.1 time = time of the extreme event
#     1.2 ID = ID of the detected extreme event after extracting them for each day
#     1.3 ID_shape = New ID after computing the t-test
all_together = read.table(file = "Link to the table", 
                          sep = '\t',header = T)
# 2. inpath = link to the folder containing the yearly stacks of the Downstream products
inpath = "link to the folder of the raster stacks"
# 3. shapefile_1 = link to the folder containing all the BACI detected shapefiles
shapefile_1 = "Link to the folder containing all the shapefiles of the extreme events"
# 4. outpath = link to the folder in which all the results are supposed to be stored
outpath = "Link to the folder in which the results are supposed to be saved in"
########################################################################################################







for (p in 1:nrow(all_together)){
time = as.character(all_together[p,2])
ID = all_together[p,1]
ID_shape = all_together[p,3]
shapefile = paste(shapefile_1,time,".shp", sep = "")

setwd(inpath)
shape <- readOGR(dsn= paste(shapefile))
shape = shape[shape$IDSel==ID_shape,]
temp_list = list.files(pattern= ".tif$")

mean = as.data.frame(matrix(NA, nrow = nlayers(brick(temp_list[[1]])), ncol = 3))
colnames(mean) = c("Mean","Year","SD")
mean$Year = seq(from= 2001, to = 2014)
year_temp = as.numeric(seq(from= 2001, to = 2014))
cols = vector(mode="character", length=nlayers(brick(temp_list[[1]])))
cols2 = cols
for(i in 1:nlayers(brick(temp_list[[1]]))){
  cols[i] = "gray" 
}

for(i in 1:nlayers(brick(temp_list[[1]]))){
  cols2[i] = "gray20" 
}

cols[as.numeric(row.names(mean[which(mean$Year == as.numeric(format(as.Date(time, format = "%Y-%m-%d"),'%Y'))) ,]))] = "red"
cols2[as.numeric(row.names(mean[which(mean$Year == as.numeric(format(as.Date(time, format = "%Y-%m-%d"),'%Y'))) ,]))] = "red4"

ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)
floor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level)


for(x in 1:length(temp_list)){
  print(paste("Processing stack",x,"of", length(temp_list)))
  rstack = brick(temp_list[[x]])
  cr <- crop(rstack, extent(shape), snap="out")                    
  fr <- rasterize(shape, cr)   
  temp <- mask(x=cr, mask=fr)
  # rm(box_tab)
  box_tab = matrix(NA, ncol = 2 , nrow = (nlayers(temp) * length(getValues(temp[[1]]))))
  box_tab = data.frame(box_tab)
  
  for(i in 1:nlayers(rstack)){
    box_values = getValues(temp[[i]])
    box_tab[c(seq(from= (((i-1)*length(getValues(temp[[i]])))+1), 
                  to=(((i-1)*length(getValues(temp[[i]])))+1)+(-1+length(getValues(temp[[i]]))), 
                  by=1)),1] = box_values
    box_tab[c(seq(from= (((i-1)*length(getValues(temp[[i]])))+1), 
                  to=(((i-1)*length(getValues(temp[[i]])))+1)+(-1+length(getValues(temp[[i]]))), 
                  by=1)),2] = rep(year_temp[i],each=length(getValues(temp[[i]])))
    if(i == nlayers(rstack)){
      png(filename=paste(outpath,'ID_',ID,"_",time,"_",gsub(".tif","",temp_list[[x]]),".png",sep = ""), width = 650, height = 480)
      box_tab = na.omit(box_tab)
      boxplot(box_tab$X1~box_tab$X2,main=paste(gsub(".tif","",temp_list[[x]])),col= cols, 
              ylim = c(floor_dec(min(box_tab$X1)),ceiling_dec(max(box_tab$X1))), las = 2, cex.main=1.4)
      legend("bottom" , inset=-0.22, c("Extreme event", "No extreme event"), col = c("red","gray"), 
             cex = 1.2, pch = c(15,15), xpd = T, xjust = 0, yjust = 0.5, y.intersp=1, x.intersp = 0.75, ncol=2,bty = "n")
      test = t.test(as.numeric(box_tab[which(box_tab$X2 == as.numeric(format(as.Date(time, format = "%Y-%m-%d"),'%Y'))),1]), 
             as.numeric(box_tab[which(box_tab$X2 != as.numeric(format(as.Date(time, format = "%Y-%m-%d"),'%Y'))),1]))
      print(paste(test$statistic))
      print(paste(test$p.value))
      mtext(paste("Value - t-statistic:",ceiling_dec(test$statistic, level = 1)), side = 3, at = 12.95, line = 1.75, font=2, cex=1.3)
      if(test$p.value >= 0.01){
      mtext(paste("p-value:",ceiling_dec(test$p.value , level = 2)), side = 3, at = 13.8, line = 0.5, font=2, cex=1.3)}
      else{
        mtext(paste("p-value <= 0.001"), side = 3, at = 13.4, line = 0.5, font=2, cex=1.3)  
      }
      dev.off()
    }
  }
}
}

