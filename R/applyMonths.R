# Authors: Shirin Taheri (taheri.shi@gmail.com); Babak Naimi (Naimi.b@gmail.com)
# Date :  March 2021
# Last update :  March 2022
# Version 1.2
# Licence GPL v3
#--------



if (!isGeneric("apply.months")) {
  setGeneric("apply.months", function(x,FUN,...)
    standardGeneric("apply.months"))
}
#------------

setMethod('apply.months', signature(x='RasterStackBrickTS'),
          function(x,FUN,...) {
            if (missing(FUN)) FUN <- mean
            
            mbase <- c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
            
            .m <- months(index(x),TRUE)
            
            w <- which(.m == mbase[1])
            xm <- calc(x@raster[[w]],FUN, na.rm=TRUE,...)
            
            for (j in mbase) {
              w <- which(.m == j)
              xm <- addLayer(xm,calc(x@raster[[w]],FUN, na.rm=TRUE,...))
            }
            names(xm) <- mbase
            xm
          }
)
#---------


setMethod('apply.months', signature(x='RasterStackBrick'),
          function(x,FUN,dates,...) {
            
            if (missing(dates) || is.null(dates)) stop('dates should be specified')
            
            
            if (missing(FUN)) FUN <- mean
            
            mbase <- c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
            
            .m <- months(dates,TRUE)
            if (!all(mbase %in% unique(.m))) stop('All the 12 months are not covered in dates!')
            
            w <- which(.m == mbase[1])
            xm <- calc(x[[w]],FUN, na.rm=TRUE,...)
            
            for (j in mbase[-1]) {
              w <- which(.m == j)
              xm <- addLayer(xm,calc(x[[w]],FUN, na.rm=TRUE,...))
            }
            names(xm) <- mbase
            xm
          }
)
#---------

setMethod('apply.months', signature(x='SpatRaster'),
          function(x,FUN,dates,...) {
            
            if (missing(dates) || is.null(dates)) stop('dates should be specified')
            
            
            if (missing(FUN)) FUN <- mean
            
            mbase <- c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
            
            .m <- months(dates,TRUE)
            if (!all(mbase %in% unique(.m))) stop('All the 12 months are not covered in dates!')
            
            
            w <- which(.m == mbase[1])
            xm <- app(x[[w]],FUN, na.rm=TRUE,...)
            
            for (j in mbase[-1]) {
              w <- which(.m == j)
              xm <- c(xm,app(x[[w]],FUN, na.rm=TRUE,...))
            }
            names(xm) <- mbase
            xm
          }
)
#=====

setMethod('apply.months', signature(x='SpatRasterTS'),
          function(x,FUN,...) {
            if (missing(FUN)) FUN <- mean
            
            mbase <- c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
            
            .m <- months(index(x),TRUE)
            
            w <- which(.m == mbase[1])
            xm <- app(x@raster[[w]],FUN, na.rm=TRUE,...)
            for (j in mbase[-1]) {
              w <- which(.m == j)
              xm <- c(xm,app(x@raster[[w]],FUN, na.rm=TRUE,...))
            }
            names(xm) <- mbase
            xm
          }
)
#---------