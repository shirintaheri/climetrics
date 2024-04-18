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

\value{
A single Raster layer (RasterLayer or SpatRaster depending on the input)
}


\author{Shirin Taheri; Babak Naimi

\email{taheri.shi@gmail.com}; \email{naimi.b@gmail.com}

}

\examples{
\donttest{
filePath <- system.file("external/", package="climetrics") # path to the dataset folder

# read the climate variables using the terra package (you can use the raster package as well):

pr <- rast(paste0(filePath,'/precip.tif'))
tmin <- rast(paste0(filePath,'/tmin.tif'))
tmax <- rast(paste0(filePath,'/tmax.tif'))

n <- readRDS(paste0(filePath,'/dates.rds')) # read corresponding dates

head(n) # Dates corresponds to the layers in climate variables (pr, tmin, tmax, tmean)

####################

# use rts function in the rts package to make a raster time series:

pr.t <- rts(pr,n) 
tmin.t <- rts(tmin,n)
tmax.t <- rts(tmax,n)

###########################
# test of the metric:

s <- sed(pr.t,tmin.t,tmax.t,t1='1991/2000',t2='2010/2020') 

plot(s, main='Standardized Local Anomalies')

######
# if the input object is SpatRaster (or RasterStack or RasterBrick) object:
# t1 and t2 would be the numbers specifying which layers correspond to time1 and time2:

s <- sed(pr,tmin,tmax,t1=1:120,t2=229:360)  

plot(s)


}



}
