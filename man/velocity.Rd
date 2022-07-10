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


filePath <- system.file("external/", package="climetrics") # path to the dataset folder

# read the climate variables using the terra package (you can use the raster package as well):

pr <- rast(paste0(filePath,'/precip.tif'))
tmax <- rast(paste0(filePath,'/tmax.tif'))

pr # has 360 layers corresponds to months of the years 1991-2020

n <- readRDS(paste0(filePath,'/dates.rds')) # read corresoinding dates

head(n) # Dates corresponds to the layers in climate variables (pr, tmin, tmax, tmean)

####################

# use rts function in the rts package to make a raster time series:

pr.t <- rts(pr,n) 
tmax.t <- rts(tmax,n)

###########################
# test of the metric:

# The extreme argument corresponds to the first and second climate variables
# (i.e., x1 and x2; precipitation and temperature) that specify the percentile of the extreme 
# condition in climate variable; here, 0.05 is used for precipitation; and 0.95 for temperature

ve <- velocity(x1=pr.t,x2=tmax.t,t1='1991/2000',t2='2010/2020')

# plot(ve, main='Velocity of Climate Change')

}
