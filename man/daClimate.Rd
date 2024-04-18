\name{daClimate}
\alias{daClimate}
\alias{daClimateC}
\alias{daClimate,RasterStackBrick-method}
\alias{daClimate,SpatRaster-method}
\alias{daClimateC,RasterStackBrick-method}
\alias{daClimateC,SpatRaster-method}
\alias{daClimate,RasterStackBrickTS-method}
\alias{daClimate,SpatRasterTS-method}


\title{Changes in distance to Analogous Climates}


\description{
Quantifies the changes in distances to analogous climates between two time periods.

}

\usage{
  daClimate(precip,tmin, tmax, tmean,t1,t2)
  daClimateC(c1,c2)
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

\details{
For a given cell, the function quantifies the distances to all cells with analogous climates in the time 1 period, i.e., belonging to the same climate class of the given cell. It also quantifies the distances to all cells that experience analogous climates in the time 2. Then, for each cell, it calculates the median of the great-circle distances (in km) below the 10th percentile of the distribution of all values, for both time 1 and time 2 periods, and maps the change over time. Negative values indicated a temporal decrease in distance, whereas positive values indicated an increase.

The daClimate first uses \code{apply.months} function to generate the monthly mean of climate parameter (results a Raster object with 12 layers correspond to 12 months) over each time period (t1 and t2); then uses the \code{kgc} function to generate the Koppen Geiger climate classification for each time period. Then the changes in the area of each class is calculated between time 1 and time 2. If the climate classification (regions) are available for time 1 and time 2, then the \code{daClimateC} can be used instead.

}
\value{
A single Raster layer (RasterLayer or SpatRaster depending on the input)
}


\author{Shirin Taheri; Babak Naimi

\email{taheri.shi@gmail.com}; \email{naimi.b@gmail.com}

}



\examples{

\donttest{
#-------
filePath <- system.file("external/", package="climetrics") # path to the dataset folder

# read the climate variables using the terra package (you can use the raster package as well):

pr <- rast(paste0(filePath,'/precip.tif'))
tmin <- rast(paste0(filePath,'/tmin.tif'))
tmax <- rast(paste0(filePath,'/tmax.tif'))
tmean <- rast(paste0(filePath,'/tmean.tif'))

pr # has 360 layers corresponds to months of the years 1991-2020

n <- readRDS(paste0(filePath,'/dates.rds')) # read corresponding dates

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
# t1 (time1) = '1991/1990' takes all layers correspond to years between 1991-01-01 to 2000-12-31
# t2 (time2) = '2010/2020' takes all layers correspond to years between 2010-01-01 to 2020-12-31

da <- daClimate(precip=pr.t,tmin=tmin.t,tmax=tmax.t,tmean=tmean.t,t1='1991/2000',t2='2010/2020')

plot(da,main="Changes in distance to Analogous Climates")


########
# Alternatively, if the climate in both times are available as climate classes
# (e.g., Koppen Geiger climate classification), changes in distance to analogous 
# climate can be quantified using the daClimateC function:

# Here, we first generate Koppen Geiger climate classification for time 1 and time 2, separately,
# To do so, we need the monthly average of climate variables (over years for each month):

# take the average of climate variables for each month over different years:

p12.1 <- apply.months(pr.t[['1991/2000']],'mean')
p12.2 <- apply.months(pr.t[['2010/2020']],'mean')

p12.1 # has 12 layers corresponding to 12 months

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

##-------- now, the kgc function can be used to generate the climate classification map:
k1 <- kgc(p=p12.1,tmin = tmin12.1,tmax=tmax12.1, tmean = tmean12.1)
k2 <- kgc(p=p12.2,tmin = tmin12.2,tmax=tmax12.2, tmean = tmean12.2)

plot(k1, main= "Koppen Geiger climate classification - 1995")
plot(k2, main= "Koppen Geiger climate classification - 2015")

# Now, given the k1 and k2 climate classes, the changes in distance to Analogous climate classes 
# between two times is calculated using the daClimateC function:

da <- daClimateC(k1,k2)

plot(da)


}


}
