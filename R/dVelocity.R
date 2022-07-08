# Authors: Shirin Taheri, taheri.shi@gmail.com; Babak Naimi (naimi.b@gmail.com)
# Date :  Oct. 2021
# Last update :  July 2022
# Version 1.2
# Licence GPL v3
#--------



.spatialgrRaster <- function(x,.lonlat) {
  if (missing(.lonlat) || !is.logical(.lonlat)) {
    if (.getProj(x) == 'longlat') .lonlat <- TRUE
    else .lonlat <- FALSE
  }
  
  nc = raster::ncol(x)
  nr = raster::nrow(x) 
  resolution = raster::xres(x) 
  basm = matrix(raster::values(x), nrow = nr, ncol = nc, byrow = TRUE)
  
  ### getting the distance in m between two cells (x direction), which means a 3 x 3 moving window
  lats = raster::yFromRow(x,1:nrow(x)) 
  longDist = NULL
  if (.lonlat) {
    for(i in 1:length(lats)) {
      longDist[i] = (pointDistance(c(0,lats[i]),c(resolution,lats[i]),longlat=TRUE)) / 1000
    } # results in km (divided by 1000), since the function gives us the result in meters, conveting the lat long values
  } else {
    for(i in 1:length(lats)) {
      longDist[i] = (pointDistance(c(0,lats[i]),c(resolution,lats[i]),longlat=FALSE))
    }
  }
  
  
  ### computing the spatial gradient
  spatialg = matrix(NA, nrow=nr, ncol=nc) 
  
  if (.lonlat) {
    for(i in 2:(nr-1)) {
      for(j in 2:(nc-1)) {
        if(!is.na(basm[i,j])) {
          xDist = longDist[i]
          
          NS = (basm[i-1,j] - basm[i+1,j]) / (2*111.3195*resolution) #number of km in 1 degree of latitude, multiplied by resolution to get cell size (y direction)
          EW = (basm[i,j-1] - basm[i,j+1]) / (2*xDist) 
          
          spatialg[i,j] = sqrt((NS^2) + (EW^2))
        }
      }
    } 
  } else {
    for(i in 2:(nr-1)) {
      for(j in 2:(nc-1)) {
        if(!is.na(basm[i,j])) {
          xDist = longDist[i]
          
          NS = (basm[i-1,j] - basm[i+1,j]) / resolution 
          EW = (basm[i,j-1] - basm[i,j+1]) / (2*xDist)
          
          spatialg[i,j] = sqrt((NS^2) + (EW^2))
        }
      }
    } 
  }
  
  #----------------
  
  ### rasterizing the spatial gradient values   
  spatialgr <- raster(x)
  raster::values(spatialgr)=c(t(spatialg))
  
  zero <- raster::Which(spatialgr < 0.000005, cells=TRUE) #see Sandel et al 2011 in Science
  spatialgr[zero] <- 0.000005
  
  spatialgr
}
#----

.spatialgrTerra <- function(x,.lonlat) {
  if (missing(.lonlat) || !is.logical(.lonlat)) {
    if (.getProj(x) == 'longlat') .lonlat <- TRUE
    else .lonlat <- FALSE
  }
  
  nc = terra::ncol(x)
  nr = terra::nrow(x) 
  resolution = terra::xres(x) 
  basm = matrix(values(x), nrow = nr, ncol = nc, byrow = TRUE)
  
  ### getting the distance in m between two cells (x direction), which means a 3 x 3 moving window
  lats = terra::yFromRow(x,1:nrow(x)) 
  
  longDist = NULL
  if (.lonlat) {
    for(i in 1:length(lats)) {
      longDist[i] = distance(matrix(c(c(0,lats[i],resolution,lats[i])),nrow=2,byrow = TRUE),lonlat=.lonlat)[1] / 1000
    } # results in km (divided by 1000), since the function gives us the result in meters, conveting the lat long values
  } else {
    for(i in 1:length(lats)) {
      longDist[i] = distance(matrix(c(c(0,lats[i],resolution,lats[i])),nrow=2,byrow = TRUE),lonlat=FALSE)[1]
    } # results in km (divided by 1000), since the function gives us the result in meters, conveting the lat long values
  }
  
  
  ### computing the spatial gradient
  spatialg = matrix(NA, nrow=nr, ncol=nc) 
  
  if (.lonlat) {
    for(i in 2:(nr-1)) {
      for(j in 2:(nc-1)) {
        if(!is.na(basm[i,j])) {
          xDist = longDist[i]
          
          NS = (basm[i-1,j] - basm[i+1,j]) / (2*111.3195*resolution) #number of km in 1 degree of latitude, multiplied by resolution to get cell size (y direction)
          EW = (basm[i,j-1] - basm[i,j+1]) / (2*xDist) 
          
          spatialg[i,j] = sqrt((NS^2) + (EW^2))
        }
      }
    } 
  } else {
    for(i in 2:(nr-1)) {
      for(j in 2:(nc-1)) {
        if(!is.na(basm[i,j])) {
          xDist = longDist[i]
          
          NS = (basm[i-1,j] - basm[i+1,j]) / resolution 
          EW = (basm[i,j-1] - basm[i,j+1]) / (2*xDist)
          
          spatialg[i,j] = sqrt((NS^2) + (EW^2))
        }
      }
    } 
  }
  
  #----------------
  
  ### rasterizing the spatial gradient values   
  spatialgr <- rast(x)
  terra::values(spatialgr) <- c(t(spatialg))
  
  zero <- which(spatialgr[] < 0.000005) #see Sandel et al 2011 in Science
  spatialgr[zero] <- 0.000005
  
  spatialgr
}
#-----------
.tempgr <- function(xt1,xt2,ny) {
  (xt2 - xt1) / ny
}
#-----------
.getVelocity <- function(s, t) {
  v <- t / s
  
  if (inherits(v,'Raster')) {
    .o <- quantile(v, prob=c(0.05,0.95))
    v[v < .o[1]] <- .o[1]
    v[v > .o[2]] <- .o[2]
    v
  } else {
    .o <- global(v, fun=quantile, prob=c(0.05,0.95),na.rm=TRUE)
    v[v < .o[1,1]] <- .o[1,1]
    v[v > .o[1,2]] <- .o[1,2]
    v
  }
  
}



if (!isGeneric("dVelocity")) {
  setGeneric("dVelocity", function(x,...,t1,t2,ny)
    standardGeneric("dVelocity"))
}

setMethod('dVelocity', signature(x='SpatRasterTS'),
          function(x,...,t1,t2,ny) {
            xx <- list(x,...)
            if (missing(t1) || missing(t2)) stop("t1 and t2 (layers' indicators corresponding to time1 and time2) are not provided!")
            
            if (missing(ny)) ny <- .getNyears(x[[1]]@time,t1,t2)
            
            xt1 <- xt2 <- list()
            
            for (i in 1:length(xx)) {
              xt1[[i]] <- mean(xx[[i]][[t1]]@raster,na.rm=TRUE)
              xt2[[i]] <- mean(xx[[i]][[t2]]@raster,na.rm=TRUE)
            }
            #---------------
            
            r <- .getScaledMultiVariateIntoOne(xt1,xt2)
            s <- .spatialgrTerra(r$t1)
            t <- .tempgr(r$t1,r$t2,ny = ny)
            .getVelocity(s,t)
          }
)



setMethod('dVelocity', signature(x='RasterStackBrickTS'),
          function(x,...,t1,t2,ny) {
            xx <- list(x,...)
            if (missing(t1) || missing(t2)) stop("t1 and t2 (layers' indicators corresponding to time1 and time2) are not provided!")
            
            if (missing(ny)) ny <- .getNyears(x[[1]]@time,t1,t2)
            
            xt1 <- xt2 <- list()
            
            for (i in 1:length(xx)) {
              xt1[[i]] <- mean(xx[[i]][[t1]]@raster,na.rm=TRUE)
              xt2[[i]] <- mean(xx[[i]][[t2]]@raster,na.rm=TRUE)
            }
            #---------------
            
            r <- .getScaledMultiVariateIntoOne(xt1,xt2)
            s <- .spatialgrRaster(r$t1)
            t <- .tempgr(r$t1,r$t2,ny = ny)
            .getVelocity(s,t)
          }
)
#---------


setMethod('dVelocity', signature(x='RasterStackBrick'),
          function(x,...,t1,t2,ny) {
            xx <- list(x,...)
            
            if (missing(ny)) stop('ny (number of years between time 1 and time 2) is not provided!')
            
            .single <- FALSE
            
            if (missing(t1) || missing(t2)) {
              if (nlayers(x) == 1 & length(xx) == 2) {
                warning('It is assumed that the first and second input raster variables correspond to a single climate variable in time 1 and time 2, respectively!')
                .single <- TRUE
              } else {
                stop("t1 and t2 (layers' indicators corresponding to time1 and time2) are not provided!")
              }
            } else {
              if (!is.numeric(t1) || !is.numeric(t2)) stop("t1 and t2 (layers' indicators corresponding to time1 and time2) should be a numeric vector")
            }
            
            
            if (.single) {
              xt1 <- list(xx[[1]])
              xt2 <- list(xx[[2]])
            } else {
              xt1 <- xt2 <- list()
              for (i in 1:length(xx)) {
                xt1[[i]] <- mean(xx[[i]][[t1]],na.rm=TRUE)
                xt2[[i]] <- mean(xx[[i]][[t2]],na.rm=TRUE)
              }
            }
            
            #---------------
            
            r <- .getScaledMultiVariateIntoOne(xt1,xt2)
            s <- .spatialgrRaster(r$t1)
            t <- .tempgr(r$t1,r$t2,ny = ny)
            .getVelocity(s,t)
            
          }
)


setMethod('dVelocity', signature('SpatRaster'),
          function(x,...,t1,t2,ny) {
            xx <- list(x,...)
            
            if (missing(ny)) stop('ny (number of years between time 1 and time 2) is not provided!')
            
            .single <- FALSE
            
            if (missing(t1) || missing(t2)) {
              if (nlyr(x) == 1 & length(xx) == 2) {
                warning('It is assumed that the first and second input raster variables correspond to a single climate variable in time 1 and time 2, respectively!')
                .single <- TRUE
              } else {
                stop("t1 and t2 (layers' indicators corresponding to time1 and time2) are not provided!")
              }
            } else {
              if (!is.numeric(t1) || !is.numeric(t2)) stop("t1 and t2 (layers' indicators corresponding to time1 and time2) should be a numeric vector")
            }
            
            
            if (.single) {
              xt1 <- list(xx[[1]])
              xt2 <- list(xx[[2]])
            } else {
              xt1 <- xt2 <- list()
              for (i in 1:length(xx)) {
                xt1[[i]] <- mean(xx[[i]][[t1]],na.rm=TRUE)
                xt2[[i]] <- mean(xx[[i]][[t2]],na.rm=TRUE)
              }
            }
            
            #---------------
            
            r <- .getScaledMultiVariateIntoOne(xt1,xt2)
            s <- .spatialgrTerra(r$t1)
            t <- .tempgr(r$t1,r$t2,ny = ny)
            .getVelocity(s,t)
            
          }
)


