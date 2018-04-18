#Opening packages
library(raster)
library(rgdal)
library(sp)
rm(list=ls())
gc()

library(ggplot2)
######################################################################################################
# This script creating an overview of all Input CABLAB data parameters via the pairs function,
# including:
# 1.
# 2.
# 3.
##########
# Data_boxplot - Table has to be provided giving all the z-scores or values for all extreme events 
# and all CABLAB data
Data_boxplot = read.table(file = "Link to the Table containing needed information", 
                          sep = ';',header = T)
#######################################################################################################


panel.cor <- function(x, y, digits=2, cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  test <- cor.test(x,y)
  text(0.5, 0.5, paste("R =",txt), cex = sqrt(sqrt(abs(r)))*3)
}


panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "red3", ...)
}

panel.regression <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
                              cex = 1, col.regres = "red", ...) 
{ 
  points(x, y, pch = pch, col = col, bg = bg, cex = cex) 
  ok <- is.finite(x) & is.finite(y) 
  if (any(ok)) 
    abline(stats::lm(y[ok] ~ x[ok]), col = col.regres, ...) 
} 

pairs(test,lower.panel=panel.regression, upper.panel=panel.cor,
      diag.panel = panel.hist, cex.labels =2.5)



