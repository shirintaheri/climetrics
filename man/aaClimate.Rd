\name{aaClimate}
\alias{aaClimate}
\alias{aaClimateC}
\alias{aaClimate,RasterStackBrick-method}
\alias{aaClimate,SpatRaster-method}
\alias{aaClimateC,RasterStackBrick-method}
\alias{aaClimateC,SpatRaster-method}
\alias{aaClimate,RasterStackBrickTS-method}
\alias{aaClimate,SpatRasterTS-method}


\title{Changes in the Area of Analogous Climate}

\usage{
  aaClimate(precip,tmin, tmax, tmean,t1,t2)
  aaClimateC(c1,c2)
}

\arguments{
  \item{precip}{A time series of precipitation as a Raster or Raster Time Series object}
  \item{tmin}{A time series of minimum temperature as a Raster or Raster Time Series object}
  \item{tmax}{A time series of maximum temperature as a Raster or Raster Time Series object}
  \item{tmean}{A time series of mean temperature as a Raster or Raster Time Series object; if not provided, it will be calculated from tmin and tmax}
  \item{t1}{a chanracter or a numeric vector, specifying the index of raster layers for time 1}
  \item{t2}{a chanracter or a numeric vector, specifying the index of raster layers for time 2}
  \item{c1}{A single layer Raster layer contains Climate classes (e.g., Koggen Geiger climate classification) for time 1}
  \item{c2}{A single layer Raster layer contains Climate classes (e.g., Koggen Geiger climate classification) for time 2}
  
}

\description{
To quantify the change in area of analogous climates (classes),the change in area occupied by a given class between time 1 (t1) and time 2 (t2) periods is quantified. For a given cell with a given climate class, the change in area of analogous climates represented the ratio (in percentage) of the difference between time 2 and time 1 area of that class to the time 1 area of the same class. Positive values indicated gains in area, negative values indicated losses, and null values reflect no change.

The aaClimate first uses \code{apply.months} function to generate the monthly mean of climate parameter (results a Raster object with 12 layers correspond to 12 months) over each time period (t1 and t2); then uses the \code{kgc} function to generate the Koggen Geiger climate classification for each time period. Then the changes in the area of each class is calculated between time 1 and time 2. If the climate classification (regions) are available for time 1 and time 2, then the \code{aaClimateC} can be used instead.

}
\examples{
\dontrun{
#-------
filePath <- system.file("external/", package="climetrics") # path to the dataset folder

# read the climate variables using the terra package (you can use the raster package as well):

pr <- rast(paste0(filePath,'/precip.tif'))
tmin <- rast(paste0(filePath,'/tmin.tif'))
tmax <- rast(paste0(filePath,'/tmax.tif'))
tmean <- rast(paste0(filePath,'/tmean.tif'))

pr # has 360 layers corresponds to months of the years 1991-2020

n <- readRDS(paste0(filePath,'/dates.rds')) # read corresoinding dates

class(n)
length(n)

head(n) # Dates corresponds to the layers in climate variables (pr, tmin, tmax, tmean)

####################

# use rts function in the rts package to make a raster time series:

pr.t <- rts(pr,n) 
tmin.t <- rts(tmin,n)
tmax.t <- rts(tmax,n)
tmean.t <- rts(tmean,n)
#------

pr.t # see the summary report of the raster time series object


###########################
# test of the metric:
#---------
#---------
# t1 (time1) = '1991/1990' takes all layers corrspond to years between 1991-01-01 to 2000-12-31
# t2 (time2) = '2010/2020' takes all layers corrspond to years between 2010-01-01 to 2020-12-31

aa <- aaClimate(precip=pr.t,tmin=tmin.t,tmax=tmax.t,tmean=tmean.t,t1='1991/2000',t2='2010/2020')

plot(aa,main="Changes in the area of Analogous Climates")


########
# Alternatively, if the climate in two times are presented as climate classes 
# (e.g., Koggen Geiger climate classification), changes in distance to analogous 
# climate can be quantified using the daClimateC function:

# Here, we first generate Koggen Geiger climate classification for time 1 and time 2, separately,
# we need monthly average of climate variables (over 12 months):

# take average of climate variables over each of 12 months:

p12.1 <- apply.months(pr.t[['1991/2000']],'mean')
p12.2 <- apply.months(pr.t[['2010/2020']],'mean')

p12.1

plot(p12.1)
#--
tmin12.1 <- apply.months(tmin.t[['1991/2000']],'mean')
tmin12.2 <- apply.months(tmin.t[['2010/2020']],'mean')
#--
tmax12.1 <- apply.months(tmax.t[['1991/2000']],'mean')
tmax12.2 <- apply.months(tmax.t[['2010/2020']],'mean')
#--
tmean12.1 <- apply.months(tmean.t[['1991/2000']],'mean')
tmean12.2 <- apply.months(tmean.t[['2010/2020']],'mean')
#--


##-------- now, the kgc function can be used to generate the climate clasification map:
k1 <- kgc(p=p12.1,tmin = tmin12.1,tmax=tmax12.1, tmean = tmean12.1)
k2 <- kgc(p=p12.2,tmin = tmin12.2,tmax=tmax12.2, tmean = tmean12.2)

plot(k1, main= "Koggen Geiger climate classification - 1995")
plot(k2, main= "Koggen Geiger climate classification - 2015")

aa <- aaClimateC(k1,k2)

plot(aa)


}


}
