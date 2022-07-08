# Authors: Shirin Taheri (taheri.shi@gmail.com); Babak Naimi (naimi.b@gmail.com)
# Date :  Nov. 2020
# Last update :  July 2022
# Version 2.4
# Licence GPL v3
#--------



.analogousClimate <- function(k1,k2) {
  r <- rast(k1)
  u <- unique(k1)[,1]
  u <- u[!is.na(u)]
  for (v in u) {
    w <- which(k1[] == v)
    w1 <- length(w)
    w2 <- length(which(k2[] == v))
    
    r[w] <- ((w2-w1) / w1) * 100
  }
  r
}
#--------------
.analogousClimateR <- function(k1,k2) {
  r <- raster(k1)
  u <- unique(k1)
  u <- u[!is.na(u)]
  for (v in u) {
    w <- which(k1[] == v)
    w1 <- length(w)
    w2 <- length(which(k2[] == v))
    
    r[w] <- ((w2-w1) / w1) * 100
  }
  r
}
#--------------
.disAnalogous <- function(k1,k2,.lon=FALSE) {
  
  if (.getProj(k1) == 'longlat') .lon <- TRUE
  
  w <- which(!is.na(k1[]))
  r <- rast(k1)
  
  for (.c in w) {
    .gv <- k1[.c][1,1]
    .xy <- xyFromCell(k1,.c)
    .xy2 <- xyFromCell(k1,which(k1[] == .gv))
    .d1 <- distance(.xy,.xy2,lonlat = .lon)[1,]
    .xy2 <- xyFromCell(k2,which(k2[] == .gv))
    if (nrow(.xy2) > 0) {
      .d2 <- distance(.xy,.xy2,lonlat = .lon)[1,]
      q <- quantile(c(.d1,.d2),0.1)
      v <- median(.d2[.d2 < q],na.rm=TRUE) - median(.d1[.d1 < q],na.rm=TRUE)
      r[.c] <- v
    } 
  }
  
  r
}
#-------------
.disAnalogousR <- function(k1,k2,.lon=FALSE) {
  
  if (.getProj(k1) == 'longlat') .lon <- TRUE
  
  w <- which(!is.na(k1[]))
  r <- raster(k1)
  
  for (.c in w) {
    .gv <- k1[.c]
    .xy <- xyFromCell(k1,.c)
    .xy2 <- xyFromCell(k1,which(k1[] == .gv))
    .d1 <- pointDistance(.xy,.xy2,lonlat = .lon)
    .xy2 <- xyFromCell(k2,which(k2[] == .gv))
    if (nrow(.xy2) > 0) {
      .d2 <- pointDistance(.xy,.xy2,lonlat = .lon)
      q <- quantile(c(.d1,.d2),0.1)
      v <- median(.d2[.d2 < q],na.rm=TRUE) - median(.d1[.d1 < q],na.rm=TRUE)
      r[.c] <- v
    }
  }
  
  r
}
#-------------
###############
if (!isGeneric("aaClimate")) {
  setGeneric("aaClimate", function(precip,tmin,tmax,tmean,t1,t2)
    standardGeneric("aaClimate"))
}


setMethod('aaClimate', signature(precip='SpatRasterTS'),
          function(precip,tmin,tmax,tmean,t1,t2) {
            
            if (missing(precip) || missing(tmin) || missing(tmax)) stop('One or more climate parameters are missing...!')
            
            if (missing(t1) || missing(t2)) stop("t1 and t2 (layers' indicators corresponding to time1 and time2) are not provided!")
            
            tmin1 <- tmin[[t1]]
            tmin2 <- tmin[[t2]]
            tmax1 <- tmax[[t1]]
            tmax2 <- tmax[[t2]]
            if (missing(tmean) || is.null(tmean)) {
              tmean <- (tmin@raster + tmax@raster) / 2
              tmean <- rts(tmean,index(tmin))
            }
            
            tmean1 <- tmean[[t1]]
            tmean2 <- tmean[[t2]]
            prt1 <- precip[[t1]]
            prt2 <- precip[[t2]]
            
            k1 <- kgc(apply.months(prt1),apply.months(tmin1),apply.months(tmax1),apply.months(tmean1))
            k2 <- kgc(apply.months(prt2),apply.months(tmin2),apply.months(tmax2),apply.months(tmean2))
            
            .analogousClimate(k1,k2)
          }
)
#-----------


setMethod('aaClimate', signature(precip='SpatRaster'),
          function(precip,tmin,tmax,tmean,t1,t2) {
            
            if (missing(precip) || missing(tmin) || missing(tmax)) stop('One or more climate parameters are missing...!')
            
            if (missing(t1) || missing(t2)) stop("t1 and t2 (layers' indicators corresponding to time1 and time2) are not provided!")
            
            if (!is.numeric(t1) || !is.numeric(t2)) stop("t1 and t2 (layers' indicators corresponding to time1 and time2) should be a numeric vector")
            
            tmin1 <- tmin[[t1]]
            tmin2 <- tmin[[t2]]
            tmax1 <- tmax[[t1]]
            tmax2 <- tmax[[t2]]
            
            if (missing(tmean) || is.null(tmean)) {
              tmean <- (tmin@raster + tmax@raster) / 2
              tmean <- rts(tmean,index(tmin))
            }
            
            tmean1 <- tmean[[t1]]
            tmean2 <- tmean[[t2]]
            prt1 <- precip[[t1]]
            prt2 <- precip[[t2]]
            
            k1 <- kgc(apply.months(prt1),apply.months(tmin1),apply.months(tmax1),apply.months(tmean1))
            k2 <- kgc(apply.months(prt2),apply.months(tmin2),apply.months(tmax2),apply.months(tmean2))
            
            .analogousClimate(k1,k2)
          }
)
#--------------

setMethod('aaClimate', signature(precip='RasterStackBrickTS'),
          function(precip,tmin,tmax,tmean,t1,t2) {
            
            if (missing(precip) || missing(tmin) || missing(tmax)) stop('One or more climate parameters are missing...!')
            
            if (missing(t1) || missing(t2)) stop("t1 and t2 (layers' indicators corresponding to time1 and time2) are not provided!")
            
            tmin1 <- tmin[[t1]]
            tmin2 <- tmin[[t2]]
            tmax1 <- tmax[[t1]]
            tmax2 <- tmax[[t2]]
            if (missing(tmean) || is.null(tmean)) {
              tmean <- (tmin@raster + tmax@raster) / 2
              tmean <- rts(tmean,index(tmin))
            }
            
            tmean1 <- tmean[[t1]]
            tmean2 <- tmean[[t2]]
            prt1 <- precip[[t1]]
            prt2 <- precip[[t2]]
            
            k1 <- kgc(apply.months(prt1),apply.months(tmin1),apply.months(tmax1),apply.months(tmean1))
            k2 <- kgc(apply.months(prt2),apply.months(tmin2),apply.months(tmax2),apply.months(tmean2))
            
            .analogousClimateR(k1,k2)
          }
)
#-----------

setMethod('aaClimate', signature(precip='RasterStackBrick'),
          function(precip,tmin,tmax,tmean,t1,t2) {
            
            if (missing(precip) || missing(tmin) || missing(tmax)) stop('One or more climate parameters are missing...!')
            
            if (missing(t1) || missing(t2)) stop("t1 and t2 (layers' indicators corresponding to time1 and time2) are not provided!")
            
            if (!is.numeric(t1) || !is.numeric(t2)) stop("t1 and t2 (layers' indicators corresponding to time1 and time2) should be a numeric vector")
            
            tmin1 <- tmin[[t1]]
            tmin2 <- tmin[[t2]]
            tmax1 <- tmax[[t1]]
            tmax2 <- tmax[[t2]]
            
            if (missing(tmean) || is.null(tmean)) {
              tmean <- (tmin@raster + tmax@raster) / 2
              tmean <- rts(tmean,index(tmin))
            }
            
            tmean1 <- tmean[[t1]]
            tmean2 <- tmean[[t2]]
            prt1 <- precip[[t1]]
            prt2 <- precip[[t2]]
            
            k1 <- kgc(apply.months(prt1),apply.months(tmin1),apply.months(tmax1),apply.months(tmean1))
            k2 <- kgc(apply.months(prt2),apply.months(tmin2),apply.months(tmax2),apply.months(tmean2))
            
            .analogousClimateR(k1,k2)
          }
)
#--------------




if (!isGeneric("aaClimateC")) {
  setGeneric("aaClimateC", function(c1,c2)
    standardGeneric("aaClimateC"))
}

setMethod('aaClimateC', signature(c1='SpatRaster'),
          function(c1,c2) {
            if (nlyr(c1) != 1 || nlyr(c2) != 1) stop('The input Raster object should have a single layer of climate classification!')
            
            .analogousClimate(c1,c2)
          }
)
#--------

setMethod('aaClimateC', signature(c1='RasterStackBrick'),
          function(c1,c2) {
            if (nlyr(c1) != 1 || nlyr(c2) != 1) stop('The input Raster object should have a single layer of climate classification!')
            
            .analogousClimateR(c1,c2)
          }
)

#--------------



###############
if (!isGeneric("daClimate")) {
  setGeneric("daClimate", function(precip,tmin,tmax,tmean,t1,t2)
    standardGeneric("daClimate"))
}


setMethod('daClimate', signature(precip='SpatRasterTS'),
          function(precip,tmin,tmax,tmean,t1,t2) {
            
            if (missing(precip) || missing(tmin) || missing(tmax)) stop('One or more climate parameters are missing...!')
            
            if (missing(t1) || missing(t2)) stop("t1 and t2 (layers' indicators corresponding to time1 and time2) are not provided!")
            
            tmin1 <- tmin[[t1]]
            tmin2 <- tmin[[t2]]
            tmax1 <- tmax[[t1]]
            tmax2 <- tmax[[t2]]
            if (missing(tmean) || is.null(tmean)) {
              tmean <- (tmin@raster + tmax@raster) / 2
              tmean <- rts(tmean,index(tmin))
            }
            
            tmean1 <- tmean[[t1]]
            tmean2 <- tmean[[t2]]
            prt1 <- precip[[t1]]
            prt2 <- precip[[t2]]
            
            k1 <- kgc(apply.months(prt1),apply.months(tmin1),apply.months(tmax1),apply.months(tmean1))
            k2 <- kgc(apply.months(prt2),apply.months(tmin2),apply.months(tmax2),apply.months(tmean2))
            
            .disAnalogous(k1,k2)
          }
)
#-----------


setMethod('daClimate', signature(precip='SpatRaster'),
          function(precip,tmin,tmax,tmean,t1,t2) {
            
            if (missing(precip) || missing(tmin) || missing(tmax)) stop('One or more climate parameters are missing...!')
            
            if (missing(t1) || missing(t2)) stop("t1 and t2 (layers' indicators corresponding to time1 and time2) are not provided!")
            
            if (!is.numeric(t1) || !is.numeric(t2)) stop("t1 and t2 (layers' indicators corresponding to time1 and time2) should be a numeric vector")
            
            tmin1 <- tmin[[t1]]
            tmin2 <- tmin[[t2]]
            tmax1 <- tmax[[t1]]
            tmax2 <- tmax[[t2]]
            
            if (missing(tmean) || is.null(tmean)) {
              tmean <- (tmin@raster + tmax@raster) / 2
              tmean <- rts(tmean,index(tmin))
            }
            
            tmean1 <- tmean[[t1]]
            tmean2 <- tmean[[t2]]
            prt1 <- precip[[t1]]
            prt2 <- precip[[t2]]
            
            k1 <- kgc(apply.months(prt1),apply.months(tmin1),apply.months(tmax1),apply.months(tmean1))
            k2 <- kgc(apply.months(prt2),apply.months(tmin2),apply.months(tmax2),apply.months(tmean2))
            
            .disAnalogous(k1,k2)
          }
)
#--------------

setMethod('daClimate', signature(precip='RasterStackBrickTS'),
          function(precip,tmin,tmax,tmean,t1,t2) {
            
            if (missing(precip) || missing(tmin) || missing(tmax)) stop('One or more climate parameters are missing...!')
            
            if (missing(t1) || missing(t2)) stop("t1 and t2 (layers' indicators corresponding to time1 and time2) are not provided!")
            
            tmin1 <- tmin[[t1]]
            tmin2 <- tmin[[t2]]
            tmax1 <- tmax[[t1]]
            tmax2 <- tmax[[t2]]
            if (missing(tmean) || is.null(tmean)) {
              tmean <- (tmin@raster + tmax@raster) / 2
              tmean <- rts(tmean,index(tmin))
            }
            
            tmean1 <- tmean[[t1]]
            tmean2 <- tmean[[t2]]
            prt1 <- precip[[t1]]
            prt2 <- precip[[t2]]
            
            k1 <- kgc(apply.months(prt1),apply.months(tmin1),apply.months(tmax1),apply.months(tmean1))
            k2 <- kgc(apply.months(prt2),apply.months(tmin2),apply.months(tmax2),apply.months(tmean2))
            
            .disAnalogousR(k1,k2)
          }
)
#-----------

setMethod('daClimate', signature(precip='RasterStackBrick'),
          function(precip,tmin,tmax,tmean,t1,t2) {
            
            if (missing(precip) || missing(tmin) || missing(tmax)) stop('One or more climate parameters are missing...!')
            
            if (missing(t1) || missing(t2)) stop("t1 and t2 (layers' indicators corresponding to time1 and time2) are not provided!")
            
            if (!is.numeric(t1) || !is.numeric(t2)) stop("t1 and t2 (layers' indicators corresponding to time1 and time2) should be a numeric vector")
            
            tmin1 <- tmin[[t1]]
            tmin2 <- tmin[[t2]]
            tmax1 <- tmax[[t1]]
            tmax2 <- tmax[[t2]]
            
            if (missing(tmean) || is.null(tmean)) {
              tmean <- (tmin@raster + tmax@raster) / 2
              tmean <- rts(tmean,index(tmin))
            }
            
            tmean1 <- tmean[[t1]]
            tmean2 <- tmean[[t2]]
            prt1 <- precip[[t1]]
            prt2 <- precip[[t2]]
            
            k1 <- kgc(apply.months(prt1),apply.months(tmin1),apply.months(tmax1),apply.months(tmean1))
            k2 <- kgc(apply.months(prt2),apply.months(tmin2),apply.months(tmax2),apply.months(tmean2))
            
            .disAnalogousR(k1,k2)
          }
)
#--------------




if (!isGeneric("daClimateC")) {
  setGeneric("daClimateC", function(c1,c2)
    standardGeneric("daClimateC"))
}

setMethod('daClimateC', signature(c1='SpatRaster'),
          function(c1,c2) {
            if (nlyr(c1) != 1 || nlyr(c2) != 1) stop('The input Raster object should have a single layer of climate classification!')
            
            .disAnalogous(c1,c2)
          }
)
#--------

setMethod('daClimateC', signature(c1='RasterStackBrick'),
          function(c1,c2) {
            if (nlyr(c1) != 1 || nlyr(c2) != 1) stop('The input Raster object should have a single layer of climate classification!')
            
            .disAnalogousR(c1,c2)
          }
)

#--------------
