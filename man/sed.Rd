\name{sed}
\alias{sed}
\alias{sed,SpatRasterTS-method}
\alias{sed,SpatRaster-method}
\alias{sed,RasterStackBrickTS-method}
\alias{sed,RasterStackBrick-method}
\alias{sed,missing-method}

\title{Standardized Local Anomalies}

\usage{
  sed(x,...,t1,t2)
}

\arguments{
  \item{x}{a Raster object or a Raster Time Series of a climate variable}
  \item{...}{additional climate variables can be included as a Raster object or a Raster Time Series}
  \item{t1}{a chanracter or a numeric vector, specifying the index of raster layers for time 1}
  \item{t2}{a chanracter or a numeric vector, specifying the index of raster layers for time 2}
}

\description{
Standardized Local Anomalies is one of the climate change metrics to quantify the change between two time periods.

}
\examples{
\dontrun{
filePath <- system.file("external/", package="climetircs") # path to the dataset folder

# read the climate variables using the terra package (you can use the raster package as well):

pr <- rast(paste0(filePath,'/precip.tif'))
tmin <- rast(paste0(filePath,'/tmin.tif'))
tmax <- rast(paste0(filePath,'/tmax.tif'))

n <- readRDS(paste0(filePath,'/dates.rds')) # read corresoinding dates

head(n) # Dates corresponds to the layers in climate variables (pr, tmin, tmax, tmean)

####################

# use rts function in the rts package to make a raster time series:

pr.t <- rts(pr,n) 
tmin.t <- rts(tmin,n)
tmax.t <- rts(tmax,n)

###########################
# test of the metric:

s <- novelClimate(pr.t,tmin.t,tmax.t,t1='1971/1990',t2='2001/2020') 
plot(s, main='Standardized Local Anomalies')



}



}
