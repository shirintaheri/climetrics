\name{ccm}
\alias{ccm}
\alias{ccm,SpatRasterTS-method}
\alias{ccm,SpatRaster-method}
\alias{ccm,RasterStackBrickTS-method}
\alias{ccm,RasterStackBrick-method}

\title{Climate Change Metrics}

\description{
A function to quantify a single or multiple climate change metrics between two time periods given the climatic raster time series variables.
}

\usage{
  ccm(x,...,stat,t1,t2,extreme,longlat,ny,names,verbose=TRUE)
}

\arguments{
  \item{x}{a Raster object or a Raster Time Series of a climate variable}
  \item{...}{additional climate variables can be included as a Raster object or a Raster Time Series}
  \item{stat}{a character specifying the metric used to generate climate change metrics, can be one or multiple metrics including 'sed','novelClimate' ,'localExtreme', 'aaClimate', 'daClimate', 'velocity', 'gVelocity', 'dVelocity' (or their abbreviations)}
  \item{t1}{a chanracter or a numeric vector, specifying the index of raster layers for time 1}
  \item{t2}{a chanracter or a numeric vector, specifying the index of raster layers for time 2}
  \item{extreme}{the percentile used for stat='lce', local climate extreme event change}
  \item{longlat}{TRUE or FALSE, specifies whether the input data is geographic or projected}
  \item{ny}{integer; specifies the number of years between time 1 and time2, if the input is Raster time series, it will be identified by the function}
  \item{names}{character with a length equal to the number of input climate variables; specifies the names of the input climate variables; for some metrics the names (precipitation, tmin, tmax, tmean) is required}
  \item{verbose}{logical (default=TRUE), specifies whether the function print messages as a text output on the screen}
}

\details{
multiple climate change metrics are available to quantify multiple dimensions of climate change between two times.

The metrics include:

- sed: Standardized Local Anomalies

- localExtreme (lce): Changes in probability of local extremes

- aaClimate (aac): Changes in areas of analogous climates

- novelClimate (nc): Novel climates

- daClimate (dac): Changes in distance to analogous climates

- velocity (ve): Climate change velocity (based on the method represented in the paper Hamnan et al. (2015))

- dVelocity (dVe): Distance-based climate change velocity (based on the method used in Garcia et al. (2014))

- gVelocity (gVe): Gradiant-based climate change velocity (based on the method represented in Burrow et al. (2011))


References:

- Burrows MT, Schoeman DS, Buckley LB, Moore P, Poloczanska ES, Brander KM, Brown C, Bruno JF, Duarte CM, Halpern BS, Holding J. (2011) The pace of shifting climate in marine and terrestrial ecosystems. Science. Nov 4;334(6056):652-5.

- Garcia, R.A., Cabeza, M., Rahbek, C. and Araújo, M.B. (2014). Multiple dimensions of climate change and their implications for biodiversity. Science, 344(6183).

- Hamann, A., Roberts, D.R., Barber, Q.E., Carroll, C. and Nielsen, S.E. (2015). Velocity of climate change algorithms for guiding conservation and management. Global Change Biology, 21(2), pp.997-1004.



}


\value{
A Raster object (RasterLayer or RasterStackBrick or SpatRaster depending on the input) with single or multiple layers correspond to the selected metrics (stat)
}


\author{Shirin Taheri; Babak Naimi

\email{taheri.shi@gmail.com}; \email{naimi.b@gmail.com}

\url{https://r-gis.net/}
}


\examples{
\donttest{
#-------
filePath <- system.file("external/", package="climetrics") # path to the dataset folder

# read the climate variables using the terra package (you can use the raster package as well):

pr <- rast(paste0(filePath,'/precip.tif'))
tmin <- rast(paste0(filePath,'/tmin.tif'))
tmax <- rast(paste0(filePath,'/tmax.tif'))
tmean <- rast(paste0(filePath,'/tmean.tif'))

pr # has 360 layers corresponds to months of the years 1991-2020


# n <- names(pr)
# n <- substring(n,2,11)
# head(n)
# n <- as.Date(n,format = '%Y.%m.%d')


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
# test of the metric:
#---------
#---------
# t1 (time1) = '1991/1990' takes all layers correspond to years between 1991-01-01 to 2000-12-31
# t2 (time2) = '2010/2020' takes all layers correspond to years between 2010-01-01 to 2020-12-31


sed <- ccm(pr.t,tmin.t,tmean.t,tmax.t,t1='1991/2000',t2='2010/2020',stat='sed')

plot(sed, main='Standardized Local Anomalies')


s2 <- ccm(pr.t,tmin.t,tmean.t,tmax.t,t1='1991/2000',t2='2010/2020',stat=c('nc','sed'))

plot(s2)



#ve <- ccm(pr.t,tmin.t,tmean.t,tmax.t,t1='1991/2000',t2='2010/2020',stat=c('gve','dve'))

#plot(ve)

# following, the extreme argument is needed for localExtreme (lce) metric:
s3 <- ccm(pr.t,tmin.t,t1='1991/2000',t2='2010/2020',stat=c('localExtreme','sed')
          ,extreme = c(0.05,0.95))

plot(s3)

}




}
