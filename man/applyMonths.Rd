\name{apply.months}
\alias{apply.months}
\alias{apply.months,RasterStackBrickTS-method}
\alias{apply.months,RasterStackBrick-method}
\alias{apply.months,SpatRaster-method}
\alias{apply.months,SpatRasterTS-method}

\title{Apply a function on rasters of each month over years}

\usage{
  apply.months(x,FUN,...)
}

\arguments{
  \item{x}{A Raster or Raster Time Series object}
  \item{FUN}{A function to apply on rasters of each month over different years; The defauls if \code{mean}}
  \item{...}{additional argument: if x is a Raster or SpatRaster object, corresponsing dates should be specified in the \code{dates} argument}
  
}

\description{
the output contains 12 raster layers corresponding to the 12 months (Jan-Dec.). For each of the available 12 months (from January to December), the corresponding rasters are selected and the \code{FUN} is applied. For example if the \code{FUN} is \code{mean}, the first layer in the output would be the mean of the rasters in January over different years, and the last later (12th layer) would be the mean of the rasters in December over different years.

}
\examples{

filePath <- system.file("external/", package="climetrics") # path to the dataset folder


pr <- rast(paste0(filePath,'/precip.tif'))

n <- readRDS(paste0(filePath,'/dates.rds')) # read corresoinding dates

head(n) # Dates corresponds to the layers in climate variables (pr, tmin, tmax, tmean)

####################

# use rts function in the rts package to make a raster time series:

pr.t <- rts(pr,n) 


###########################

p12 <- apply.months(pr.t,'mean')


p12

# plot(p12)


}
