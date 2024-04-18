\name{kgc}
\alias{kgc}
\alias{kgc,RasterStackBrick,RasterStackBrick,RasterStackBrick-method}
\alias{kgc,SpatRaster,SpatRaster,SpatRaster-method}

\title{Koggen Geiger climate classification}


\description{
Given the monthly mean (12 months) of Precipitation, minimum, maximum, and mean temperature datasets, this function identify climate zones following Koppen Geiger climate classification rules.


}

\usage{
  kgc(p,tmin, tmax, tmean)
}

\arguments{
  \item{p}{A Raster object contains 12 months precipitation}
  \item{tmin}{A Raster object contains 12 months minimum temperature}
  \item{tmax}{A Raster object contains 12 months maximum temperature}
  \item{tmean}{A Raster object contains 12 months mean temperature}
}

\details{
If the input object is not a 12 month Raster object but a time series, the function first calls the \code{apply.months} to get the 12 month raster object.

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
tmean <- rast(paste0(filePath,'/tmean.tif'))

pr # has 360 layers corresponds to months of the years 1991-2020

n <- readRDS(paste0(filePath,'/dates.rds')) # read corresponding dates

class(n)
length(n)

head(n) # Dates corresponds to the layers in climate variables (pr, tmin, tmax, tmean)

####################

# use rts function in the rts package to make a raster time series:

pr.t <- rts(pr,n) 
tmin.t <- rts(tmin,n)
tmax.t <- rts(tmax,n)
tmean.t <- rts(tmean,n)
#------

pr.t # see the summary report of the raster time series object


###########################

# To generate Koppen Geiger climate classification, we need the monthly average of 
# climate variables (over years for each month):

# you can use apply.months function to take the average of climate variables for each month
# over different years that generates 12 layers corresponds to 12 months:

p12 <- apply.months(pr.t[['1991/2000']],'mean')


p12

# plot(p12)

#--
tmin12 <- apply.months(tmin.t[['1991/2000']],'mean')

tmax12 <- apply.months(tmax.t[['1991/2000']],'mean')

tmean12 <- apply.months(tmean.t[['1991/2000']],'mean')


##-------- now, the kgc function can be used to generate the climate classification map:
k1 <- kgc(p=p12,tmin = tmin12,tmax=tmax12, tmean = tmean12)

k1

plot(k1)
}



}
