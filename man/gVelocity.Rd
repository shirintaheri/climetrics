\name{gVelocity}
\alias{gVelocity}
\alias{gVelocity,RasterStackBrickTS-method}
\alias{gVelocity,RasterStackBrick-method}
\alias{gVelocity,SpatRasterTS-method}
\alias{gVelocity,SpatRaster-method}


\title{Gradient-based Velocity of Climate Change}

\usage{
  gVelocity(x,...)
}

\arguments{
  \item{x}{a Raster object or a Raster Time Series of climate variable}
  
  \item{...}{additional arguments; not implemented}
}

\description{
The method is developed based on Burrow et al. (2011) which is a gradient-based velocity of climate change. Climate velocity can inform conservation planning and action in a warming world.

Distance-based velocity is also implemented in the package and is available through the ccm function.

Burrows MT, Schoeman DS, Buckley LB, Moore P, Poloczanska ES, Brander KM, Brown C, Bruno JF, Duarte CM, Halpern BS, Holding J. The pace of shifting climate in marine and terrestrial ecosystems. Science. 2011 Nov 4;334(6056):652-5..
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

#gv <- gVelocity(pr)

#plot(gv, main='Gradiant-based Velocity for Precipitation')

}

}
