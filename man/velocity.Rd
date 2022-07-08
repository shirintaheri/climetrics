\name{velocity}
\alias{velocity}
\alias{velocity,RasterStackBrickTS,RasterStackBrickTS-method}
\alias{velocity,RasterStackBrick,RasterStackBrick-method}
\alias{velocity,SpatRaster,SpatRaster-method}
\alias{velocity,SpatRasterTS,SpatRasterTS-method}

\title{Velocity of Climate Change}

\usage{
  velocity(x1,x2,t1, t2,...)
}

\arguments{
\item{x1}{The first input climate variable that can be a RasterStack or RasterBrick, or a raster time series}
\item{x2}{The second input climate variable that can be a RasterStack or RasterBrick, or a raster time series}
  
  \item{t1}{a chanracter or a numeric vector, specifying the index of raster layers for time 1}
  \item{t2}{a chanracter or a numeric vector, specifying the index of raster layers for time 2}
  \item{...}{not implemented}
  
}

\description{
The method is developed based on the method represented in the paper Hamnan et al. (2015). Two additional functions including Distance-based (\code{dVelocity}) and Gradiant-based velocity  (\code{gVelocity}) are also implemented and available in the package.

- Hamann, A., Roberts, D.R., Barber, Q.E., Carroll, C. and Nielsen, S.E. (2015). Velocity of climate change algorithms for guiding conservation and management. Global Change Biology, 21(2), pp.997-1004.


}
\examples{

# path <- system.file("external/", package="climetrics") # path to the dataset folder



}
