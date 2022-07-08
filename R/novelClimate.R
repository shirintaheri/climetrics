# Authors: Shirin Taheri (taheri.shi@gmail.com); Babak Naimi (naimi.b@gmail.com)
# Date :  Nov. 2020
# Last update :  July 2022
# Version 2.3
# Licence GPL v3
#--------

.getRasterList <- function(x1,...,t1,t2) {
  xx <- list(x1,...)
  
  xt1 <- list(xx[[1]][[t1]])
  xt2 <- list(xx[[1]][[t2]])
  
  if (length(xx) > 1) {
    for (i in 2:length(xx)) {
      xt1 <- c(xt1,xx[[i]][[t1]])
      xt2 <- c(xt2,xx[[i]][[t2]])
    }
  }
  list(T1=xt1,T2=xt2)
}

#--------------
.nc <- function(xx,t1,t2) {
  
  for (i in 1:length(xx)) {
    x <- xx[[i]]
    x1 <- x[[t1]]
    x2 <- x[[t2]]
    mn1 <- app(x1@raster,mean,na.rm=TRUE)
    mn2 <- app(x2@raster,mean,na.rm=TRUE)
    v <- app(x1@raster,var,na.rm=TRUE)
    
    w <- which(!is.na(mn2[]))
    ww <- data.frame(mn1[w],v[w])
    
    wm <- sapply(mn2[w][,1],function(x,...) {
      min(((x - ww[,1]) ^ 2) / ww[,2])
    },na.rm=TRUE)
    wm2 <- rast(mn2)
    wm2[w] <- wm
    
    if (i == 1) .sed <- wm2
    else .sed <- .sed + wm2
  }
  
  .s <- sqrt(.sed)
  .s <- ifel(is.infinite(.s),0,.s)
  .s
}
#-----------------
.ncR <- function(xx,t1,t2) {
  for (i in 1:length(xx)) {
    x <- xx[[i]]
    x1 <- x[[t1]]
    x2 <- x[[t2]]
    mn1 <- app(x1,mean,na.rm=TRUE)
    mn2 <- app(x2,mean,na.rm=TRUE)
    v <- app(x1,var,na.rm=TRUE)
    
    w <- which(!is.na(mn2[]))
    ww <- data.frame(mn1[w],v[w])
    
    wm <- sapply(mn2[w][,1],function(x,...) {
      min(((x - ww[,1]) ^ 2) / ww[,2])
    },na.rm=TRUE)
    wm2 <- rast(mn2)
    wm2[w] <- wm
    
    if (i == 1) .sed <- wm2
    else .sed <- .sed + wm2
  }
  
  .s <- sqrt(.sed)
  .s <- ifel(is.infinite(.s),0,.s)
  .s
}
#----
.ncRas <- function(xx,t1,t2) {
  for (i in 1:length(xx)) {
    x <- xx[[i]]
    x1 <- x[[t1]]
    x2 <- x[[t2]]
    if (inherits(x1,'RasterStackBrickTS')) {
      x1 <- x1@raster
      x2 <- x2@raster
    }
    mn1 <- calc(x1,mean,na.rm=TRUE)
    mn2 <- calc(x2,mean,na.rm=TRUE)
    v <- calc(x1,var,na.rm=TRUE)
    
    w <- which(!is.na(mn2[]))
    ww <- data.frame(mn1[w],v[w])
    
    wm <- sapply(mn2[w][,1],function(x,...) {
      min(((x - ww[,1]) ^ 2) / ww[,2])
    },na.rm=TRUE)
    wm2 <- raster(mn2)
    wm2[w] <- wm
    
    if (i == 1) .sed <- wm2
    else .sed <- .sed + wm2
  }
  
  .s <- sqrt(.sed)
  .s[] <- ifelse(is.infinite(.s[]),0,.s[])
  .s
}
#------------
.ncList <- function(x1,x2) {
  for (i in 1:length(x1)) {
    mn1 <- app(x1[[i]],mean,na.rm=TRUE)
    mn2 <- app(x2[[i]],mean,na.rm=TRUE)
    v <- app(x1[[i]],var,na.rm=TRUE)
    
    w <- which(!is.na(mn2[]))
    ww <- data.frame(mn1[w],v[w])
    
    wm <- sapply(mn2[w][,1],function(x,...) {
      min(((x - ww[,1]) ^ 2) / ww[,2])
    },na.rm=TRUE)
    wm2 <- rast(mn2)
    wm2[w] <- wm
    
    if (i == 1) .sed <- wm2
    else .sed <- .sed + wm2
  }
  
  .n <- sqrt(.sed)
  .n <- ifel(is.infinite(.n),0,.n)
  .q <- global(.n,quantile,probs=0.99,na.rm=TRUE)[1,1]
  .n <- ifel(.n >= .q,.q,.n)
  .n
}
#-----------------

.ncListR <- function(x1,x2) {
  for (i in 1:length(x1)) {
    mn1 <- calc(x1[[i]],mean,na.rm=TRUE)
    mn2 <- calc(x2[[i]],mean,na.rm=TRUE)
    v <- calc(x1[[i]],var,na.rm=TRUE)
    
    w <- which(!is.na(mn2[]))
    ww <- data.frame(mn1[w],v[w])
    
    wm <- sapply(mn2[w][,1],function(x,...) {
      min(((x - ww[,1]) ^ 2) / ww[,2])
    },na.rm=TRUE)
    wm2 <- rast(mn2)
    wm2[w] <- wm
    
    if (i == 1) .sed <- wm2
    else .sed <- .sed + wm2
  }
  
  .n <- sqrt(.sed)
  .n <- ifelse(is.infinite(.n[]),0,.n[])
  .q <- quantile(.n,probs=0.99,na.rm=TRUE)
  .n[] <- ifelse(.n[] >= .q,.q,.n[])
  .n
}
#-----------------


###############
if (!isGeneric("novelClimate")) {
  setGeneric("novelClimate", function(x,...,t1,t2)
    standardGeneric("novelClimate"))
}


setMethod('novelClimate', signature(x='SpatRasterTS'),
          function(x,...,t1,t2) {
            xx <- list(x,...)
            .n <- .nc(xx,t1=t1,t2=t2)
            .q <- global(.n,quantile,probs=0.99,na.rm=TRUE)[1,1]
            .n <- ifel(.n >= .q,.q,.n)
            .n
          }
)


setMethod('novelClimate', signature(x='RasterStackBrickTS'),
          function(x,...,t1,t2) {
            xx <- list(x,...)
            
            if (.require("terra")) {
              for (i in 1:length(xx)) {
                xx[[i]] <- rts(rast(xx[[i]]@raster),index(xx[[i]]@time))
              }
              .n <- .nc(xx,t1=t1,t2=t2)
              .q <- global(.n,quantile,probs=0.99,na.rm=TRUE)[1,1]
              .n <- ifel(.n >= .q,.q,.n)
              raster(.n)
            } else {
              .n <- .ncRas(xx,t1=t1,t2=t2)
              .q <- quantile(.n,probs=0.99,na.rm=TRUE)
              .n[] <- ifelse(.n[] >= .q,.q,.n[])
              .n
            }
          }
)


setMethod('novelClimate', signature(x='SpatRaster'),
          function(x,...,t1,t2) {
            xx <- list(x,...)
            if (missing(t1) || missing(t2)) stop("t1 and/or t2 (layers' indicators corresponding to time1 and time2) are not provided!")
            
            if (!is.numeric(t1) || !is.numeric(t2)) stop("t1 and t2 (layers' indicators corresponding to time1 and time2) should be a numeric vector")
            
            .n <- .ncR(xx,t1=t1,t2=t2)
            .q <- global(.n,quantile,probs=0.99,na.rm=TRUE)[1,1]
            .n <- ifel(.n >= .q,.q,.n)
            .n
          }
)


setMethod('novelClimate', signature(x='RasterStackBrick'),
          function(x,...,t1,t2) {
            xx <- list(x,...)
            if (missing(t1) || missing(t2)) stop("t1 and/or t2 (layers' indicators corresponding to time1 and time2) are not provided!")
            
            if (!is.numeric(t1) || !is.numeric(t2)) stop("t1 and t2 (layers' indicators corresponding to time1 and time2) should be a numeric vector")
            
            if (.require("terra")) {
              for (i in 1:length(xx)) {
                xx[[i]] <- rast(xx[[i]])
              }
              .n <- .ncR(xx,t1=t1,t2=t2)
              .q <- global(.n,quantile,probs=0.99,na.rm=TRUE)[1,1]
              .n <- ifel(.n >= .q,.q,.n)
              raster(.n)
            } else {
              .n <- .ncRas(xx,t1=t1,t2=t2)
              .q <- quantile(.n,probs=0.99,na.rm=TRUE)
              .n[] <- ifelse(.n[] >= .q,.q,.n[])
              .n
            }
          }
)
#===============

setMethod('novelClimate', signature(x='missing'),
          function(x,...,t1,t2) {
            if (missing(t1) && missing(t2)) stop('No Raster Object are provided!')
            if (!is.list(t1) && !is.list(t2)) stop('The input Raster objects should be provided as time series, or as separate lists in t1 and t2')
            if (length(t1) != length(t2)) stop('The lists provided in t1 and t2 has not the same length!')
            
            cls <- unique(c(sapply(t1,function(x) class(x)),sapply(t2,function(x) class(x))))
            
            
            if (length(cls) > 1) stop('The classes of objects in the input lists (t1 and t2) are not the same!')
            
            if (!cls %in% c('SpatRasterTS','SpatRaster','RasterStack','RasterBrick','RasterStackTS','RasterBrickTS')) stop('The input objects provided in t1 and t2 lists are unknown (should be either of SpatRaster, RasterStack, RasterBrick, SpatRasterTS, RasterStackTS, RasterBrickTS!')
            
            if (cls == 'SpatRaster') .ncList(x1=t1,x2=t2)
            else if (cls == 'SpatRasterTS') {
              for (i in 1:length(t1)) {
                t1[[i]] <- t1[[i]]@raster
                t2[[i]] <- t2[[i]]@raster
              }
              .ncList(x1=t1,x2=t2)
            } else {
              
              if (.require("terra")) {
                if (cls %in% c('RasterStack','RasterBrick')) {
                  t1[[i]] <- rast(t1[[i]])
                  t2[[i]] <- rast(t2[[i]])
                } else {
                  t1[[i]] <- rast(t1[[i]]@raster)
                  t2[[i]] <- rast(t2[[i]]@raster)
                }
                
                raster(.ncList(x1=t1,x2=t2))
              } else {
                if (cls %in% c('RasterStack','RasterBrick')) {
                  .ncListR(t1,t2)
                } else {
                  for (i in 1:length(t1)) {
                    t1[[i]] <- t1[[i]]@raster
                    t2[[i]] <- t2[[i]]@raster
                  }
                  .ncListR(t1,t2)
                }
              }
            }
          }
)




# .nc <- function(p=NULL,tmin=NULL,tmax=NULL,tmean=NULL,t1,t2) {
#   xx <- list()
#   if (!is.null(p)) xx <- c(p,xx)
#   if (!is.null(tmin)) xx <- c(tmin,xx)
#   if (!is.null(tmax)) xx <- c(tmax,xx)
#   if (!is.null(tmean)) xx <- c(tmean,xx)
#   
#   for (i in 1:length(xx)) {
#     x <- xx[[i]]
#     x1 <- x[[t1]]
#     x2 <- x[[t2]]
#     mn1 <- app(x1@raster,mean,na.rm=TRUE)
#     mn2 <- app(x2@raster,mean,na.rm=TRUE)
#     v <- app(x1@raster,var,na.rm=TRUE)
#     
#     w <- which(!is.na(mn2[]))
#     ww <- data.frame(mn1[w],v[w])
#     
#     wm <- sapply(mn2[w][,1],function(x,...) {
#       min(((x - ww[,1]) ^ 2) / ww[,2])
#     },na.rm=TRUE)
#     wm2 <- rast(mn2)
#     wm2[w] <- wm
#     # wm <- app(mn2,function(x) {
#     #   min(((x - ww[,1]) ^ 2) / ww[,2])
#     # })
#     
#     if (i == 1) .sed <- wm2
#     else .sed <- .sed + wm2
#   }
#   
#   sqrt(.sed)
# }
#-----------------

#--------
# .ncR <- function(p=NULL,tmin=NULL,tmax=NULL,tmean=NULL,t1,t2) {
#   xx <- list()
#   if (!is.null(p)) xx <- c(p,xx)
#   if (!is.null(tmin)) xx <- c(tmin,xx)
#   if (!is.null(tmax)) xx <- c(tmax,xx)
#   if (!is.null(tmean)) xx <- c(tmean,xx)
#   
#   for (i in 1:length(xx)) {
#     x <- xx[[i]]
#     x1 <- x[[t1]]
#     x2 <- x[[t2]]
#     mn1 <- app(rast(x1),mean,na.rm=TRUE)
#     mn2 <- app(rast(x2),mean,na.rm=TRUE)
#     v <- app(rast(x1),var,na.rm=TRUE)
#     
#     w <- which(!is.na(mn2[]))
#     ww <- data.frame(mn1[w],v[w])
#     
#     wm <- sapply(mn2[w][,1],function(x,...) {
#       min(((x - ww[,1]) ^ 2) / ww[,2])
#     },na.rm=TRUE)
#     wm2 <- rast(mn2)
#     wm2[w] <- wm
#     
#     # wm <- app(mn2,function(x) {
#     #   min(((x - ww[,1]) ^ 2) / ww[,2])
#     # })
#     # 
#     if (i == 1) .sed <- wm2
#     else .sed <- .sed + wm2
#   }
#   
#   sqrt(.sed)
# }
#-----------------