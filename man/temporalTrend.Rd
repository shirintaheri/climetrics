\name{temporalTrend}
\alias{temporalTrend}
\alias{temporalTrend,RasterStackBrickTS-method}
\alias{temporalTrend,RasterStackBrick-method}
\alias{temporalTrend,SpatRasterTS-method}
\alias{temporalTrend,SpatRaster-method}



\title{Temporal Trend (Slope) of a climate variable}

\usage{
  temporalTrend(x,...)
}

\arguments{
  \item{x}{a Raster object or a Raster Time Series of climate variable}
  
  \item{...}{additional arguments; not implemented}
}

\description{
The method calculates the trend of a climate variale changes over time. The function use a Raster time series as input and returns a Raster object that represent the trend at the pixel level.
}
\examples{
filePath <- system.file("external/", package="climetrics") # path to the dataset folder

# read the climate variables using the terra package (you can use the raster package as well):

pr <- rast(paste0(filePath,'/precip.tif'))

tr <- temporalTrend(pr)

# plot(tr, main='Trend (slope) of Precipitation time series')

}
