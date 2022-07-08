\name{localExtreme}
\alias{localExtreme}
\alias{localExtreme,SpatRasterTS-method}
\alias{localExtreme,SpatRaster-method}
\alias{localExtreme,RasterStackBrickTS-method}
\alias{localExtreme,RasterStackBrick-method}
\alias{localExtreme,list,list-method}
\alias{localExtreme,missing,missing,list,list-method}
\alias{localExtreme,SpatRasterTS,ANY,ANY,ANY-method}
\alias{localExtreme,SpatRaster,ANY,ANY,ANY-method}
\alias{localExtreme,RasterStackBrickTS,ANY,ANY,ANY-method}
\alias{localExtreme,RasterStackBrick,ANY,ANY,ANY-method}
\alias{localExtreme,list,list,ANY,ANY-method}



\title{Changes in probability of local extremes}

\usage{
  localExtreme(x1,x2,t1,t2,extreme)
}

\arguments{
  \item{x1}{Time Series of the first climate variable as a Raster or Raster Time Series Object}
  \item{x2}{Time Series of the second climate variable as a Raster or Raster Time Series Object; \code{x1} and \code{x2} can be lists of time series of climate variables in time 1 and time 2, respectively}
  \item{t1}{a chanracter or a numeric vector, specifying the index of raster layers for time 1}
  \item{t2}{a chanracter or a numeric vector, specifying the index of raster layers for time 2}
  \item{extreme}{a numeric vector with a length of one or two (depends on whether one or two climate variables are used as input), specifying the percentile of extreme value in the input (first value corresponds to x1 and second corresponds to x2 climate variables) }
}

\description{
Probability of changes in local extremes is calculated at each cell based on either one or two climate variables (e.g., precipitation and temperature). For each climate variable, the extreme value (percentile of the extreme in the variable distribution) should be specified by user in the extreme parameter. For instance, if temperature and precipitation are specified as the first and second inputs (x1 and x2), extreme=c(0.95,0.05) may be used that consideres 95th and 5th percentiles of the distribution of the temperature and precipitation time series, respectively over the first time period, and then it calculates the probability that the time1 percentile will be exceeded in time2. 
If two climate variables are used (e.g., temperature and precipitation), for each cell, the two probabilities are summed from which the product of the two probabilities are subtracted to obtain a measure of the probability of occurrence of either of the two extreme events (to avoid counting probabilities twice). The  probability of time_1 extreme climates in time_2, is then subtracted from the probability of time_1 extreme climates to obtain the changes in the probability of extreme climate variables. Positive values indicated increased probability in time 2, whereas negative values indicated a decrease. 

}
\examples{

\dontrun{

filePath <- system.file("external/", package="climertics") # path to the dataset folder

# read the climate variables using the terra package (you can use the raster package as well):

pr <- rast(paste0(filePath,'/precip.tif'))
tmax <- rast(paste0(filePath,'/tmax.tif'))

pr # has 612 layers corresponds to months of the years 1970-2020

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

le <- localExtreme(x1=pr.t,x2=tmax.t,t1='1971/1990',t2='2001/2020', extreme = c(0.05, 0.95))

plot(le, main='Probability of Changes in Local Climate Extreme')


}



}
