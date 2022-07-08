\name{dVelocity}
\alias{dVelocity}
\alias{dVelocity,RasterStackBrickTS-method}
\alias{dVelocity,RasterStackBrick-method}
\alias{dVelocity,SpatRaster-method}
\alias{dVelocity,SpatRasterTS-method}

\title{Distant-based Velocity of Climate Change}

\usage{
  dVelocity(x,...,t1, t2, ny)
}

\arguments{
\item{x}{The first input climate variable that can be a RasterStack or RasterBrick, or a raster time series}
  \item{...}{additional input climate variables in case of using multiple climate variables that can be entered as \code{v1, v2, ...}}
  
  \item{t1}{a chanracter or a numeric vector, specifying the index of raster layers for time 1}
  \item{t2}{a chanracter or a numeric vector, specifying the index of raster layers for time 2}
  \item{ny}{integer; specifies the number of years between time 1 and time2, if the input is Raster time series, it is not needed as it will be identified by the function}
  
}

\description{
The method is developed based on the approach used in Garcia et al. (2014)
Gradiant-based velocity is also implemented in the package and is available through the \code{gVelocity} function.

Garcia, R.A., Cabeza, M., Rahbek, C. and Araújo, M.B. (2014). Multiple dimensions of climate change and their implications for biodiversity. Science, 344(6183).
}
\examples{
#-------
#path <- system.file("external/", package="climetrics") # path to the dataset folder



}
