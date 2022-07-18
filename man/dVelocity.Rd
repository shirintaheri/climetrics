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

\value{
A single Raster layer (RasterLayer or SpatRaster depending on the input)
}


\author{Shirin Taheri; Babak Naimi

\email{taheri.shi@gmail.com}; \email{naimi.b@gmail.com}

}


\examples{

filePath <- system.file("external/", package="climetrics") # path to the dataset folder

pr <- rast(paste0(filePath,'/precip.tif'))
tmean <- rast(paste0(filePath,'/tmean.tif'))

n <- readRDS(paste0(filePath,'/dates.rds')) # corresoinding dates

head(n) # Dates corresponds to the layers in climate variables

####################

# use rts function in the rts package to make a raster time series:

pr.t <- rts(pr,n) 
tmean.t <- rts(tmean,n)
###########################


dv <- dVelocity(pr.t,tmean.t,t1='1991/2000',t2='2010/2020')

plot(dv)


}
