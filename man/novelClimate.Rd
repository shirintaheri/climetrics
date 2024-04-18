\name{novelClimate}
\alias{novelClimate}
\alias{novelClimate,SpatRasterTS-method}
\alias{novelClimate,SpatRaster-method}
\alias{novelClimate,RasterStackBrickTS-method}
\alias{novelClimate,RasterStackBrick-method}
\alias{novelClimate,missing-method}

\title{Novel Climates}

\usage{
  novelClimate(x,...,t1,t2)
}

\arguments{
  \item{x}{a Raster object or a Raster Time Series of a climate variable}
  \item{...}{additional climate variables can be included as a Raster object or a Raster Time Series}
  \item{t1}{a character or a numeric vector, specifying the index of raster layers for time 1; if \code{x} is missing, \code{t1} and \code{t2} can be list of climate variables in time 1 and time 2, respectively.}
  \item{t2}{a character or a numeric vector, specifying the index of raster layers for time 2;  if \code{x} is missing, \code{t1} and \code{t2} can be list of climate variables in time 1 and time 2, respectively.}
}

\description{
Novel climate is one of the climate change metrics to quantify the dissimilarities between climate parameters between two time periods. This metric uses standard euclidean distance between each cell in Time 2 and all cells in Time 1 and  retains the minimum of those distances. The inter-annual standard deviation for each parameter is used for the standardization. The larger the score, the more dissimilar the climate in Time 2 is in relation to the global pool of potential climatic analogues.

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
tmean <- rast(paste0(filePath,'/tmean.tif'))

pr # has 360 layers corresponds to months of the years 1991-2020

n <- readRDS(paste0(filePath,'/dates.rds')) # read corresponding dates

head(n) # Dates corresponds to the layers in climate variables (pr, tmin, tmax, tmean)

####################

# use rts function in the rts package to make a raster time series:

pr.t <- rts(pr,n) 
tmean.t <- rts(tmean,n)

###########################
# test of the metric:

n1 <- novelClimate(pr.t,tmean.t,t1='1991/2000',t2='2010/2020') 

plot(n1, main='Novel Climate')

######
# if the input object is SpatRaster (or RasterStack or RasterBrick) object:
# t1 and t2 would be the numbers specifying which layers correspond to time1 and time2:

n2 <- novelClimate(pr,tmean,t1=1:120,t2=229:360)  

plot(n2)


}

}
