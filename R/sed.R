# Authors: Shirin Taheri (taheri.shi@gmail.com); Babak Naimi (naimi.b@gmail.com)
# Date :  Nov. 2020
# Last update :  July 2022
# Version 2.3
# Licence GPL v3
#--------


# standardized local anomalies


.sed <- function(xx,t1,t2) {
  
  for (i in 1:length(xx)) {
    x <- xx[[i]]
    x1 <- x[[t1]]
    x2 <- x[[t2]]
    mn1 <- app(x1@raster,mean,na.rm=TRUE)
    mn2 <- app(x2@raster,mean,na.rm=TRUE)
    v <- app(x1@raster,var,na.rm=TRUE)
    
    if (i == 1) .s <- (mn1 - mn2) ^ 2 / v
    else .s <- .s + ((mn1 - mn2) ^ 2 / v)
  }
  
  .s <- sqrt(.s)
  .s <- ifel(is.infinite(.s),0,.s)
  .s
}
# -----------
.sed0 <- function(xx,t1,t2) {
  
  for (i in 1:length(xx)) {
    x <- xx[[i]]
    x1 <- x[[t1]]
    x2 <- x[[t2]]
    mn1 <- app(x1,mean,na.rm=TRUE)
    mn2 <- app(x2,mean,na.rm=TRUE)
    v <- app(x1,var,na.rm=TRUE)
    
    if (i == 1) .s <- (mn1 - mn2) ^ 2 / v
    else .s <- .s + ((mn1 - mn2) ^ 2 / v)
  }
  
  .s <- sqrt(.s)
  .s <- ifel(is.infinite(.s),0,.s)
  .s
}
#------
.sedRas <- function(xx,t1,t2) {
  
  for (i in 1:length(xx)) {
    x <- xx[[i]]
    x1 <- x[[t1]]
    x2 <- x[[t2]]
    mn1 <- calc(x1,mean,na.rm=TRUE)
    mn2 <- calc(x2,mean,na.rm=TRUE)
    v <- calc(x1,var,na.rm=TRUE)
    
    if (i == 1) .s <- (mn1 - mn2) ^ 2 / v
    else .s <- .s + ((mn1 - mn2) ^ 2 / v)
  }
  
  .s <- sqrt(.s)
  .s[] <- ifelse(is.infinite(.s[]),0,.s[])
  .s
}
#------


.sedList <- function(x1,x2) {
  # x1 and x2 are list of RasterObject (SpatRaster)
  for (i in 1:length(x1)) {
    mn1 <- app(x1[[i]],mean,na.rm=TRUE)
    mn2 <- app(x2[[i]],mean,na.rm=TRUE)
    v <- app(x1[[i]],var,na.rm=TRUE)
    
    if (i == 1) .s <- (mn1 - mn2) ^ 2 / v
    else .s <- .s + ((mn1 - mn2) ^ 2 / v)
  }
  
  .s <- sqrt(.s)
  .s <- ifel(is.infinite(.s),0,.s)
  .q <- global(.s,quantile,probs=0.995,na.rm=TRUE)[1,1]
  .s <- ifel(.s >= .q,.q,.s)
  .s
}
#----------
.sedListR <- function(x1,x2) {
  # x1 and x2 are list of RasterObject (RasterStackBrick)
  for (i in 1:length(x1)) {
    mn1 <- calc(x1[[i]],mean,na.rm=TRUE)
    mn2 <- calc(x2[[i]],mean,na.rm=TRUE)
    v <- calc(x1[[i]],var,na.rm=TRUE)
    
    if (i == 1) .s <- (mn1 - mn2) ^ 2 / v
    else .s <- .s + ((mn1 - mn2) ^ 2 / v)
  }
  
  .s <- sqrt(.s)
  .s[] <- ifelse(is.infinite(.s[]),0,.s[])
  .q <- quantile(.s,probs=0.995,na.rm=TRUE)
  .s[] <- ifelse(.s[] >= .q,.q,.s[])
  .s
}
#----------
if (!isGeneric("sed")) {
  setGeneric("sed", function(x,...,t1,t2)
    standardGeneric("sed"))
}


setMethod('sed', signature(x='SpatRasterTS'),
          function(x,...,t1,t2) {
            xx <- list(x,...)
            .s <- .sed(xx,t1=t1,t2=t2)
            .q <- global(.s,quantile,probs=0.995,na.rm=TRUE)[1,1]
            .s <- ifel(.s >= .q,.q,.s)
            .s
          }
)


setMethod('sed', signature(x='RasterStackBrickTS'),
          function(x,...,t1,t2) {
            xx <- list(x,...)
            
            if (.require('terra')) {
              for (i in 1:length(xx)) {
                xx[[i]] <- rts(rast(xx[[i]]@raster),index(xx[[i]]@time))
              }
              
              .s <- .sed(xx,t1=t1,t2=t2)
              .q <- global(.s,quantile,probs=0.995,na.rm=TRUE)[1,1]
              .s <- ifel(.s >= .q,.q,.s)
              raster(.s)
            } else {
              .s <- .sedRas(xx,t1=t1,t2=t2)
              .q <- quantile(.s,probs=0.995,na.rm=TRUE)
              .s[] <- ifelse(.s[] >= .q,.q,.s[])
              .s
            }
          }
)


setMethod('sed', signature(x='SpatRaster'),
          function(x,...,t1,t2) {
            xx <- list(x,...)
            if (missing(t1) || missing(t2)) stop("t1 and/or t2 (layers' indicators corresponding to time1 and time2) are not provided!")
            
            if (!is.numeric(t1) || !is.numeric(t2)) stop("t1 and t2 (layers' indicators corresponding to time1 and time2) should be a numeric vector")
            
            .s <- .sed0(xx,t1=t1,t2=t2)
            .q <- global(.s,quantile,probs=0.995,na.rm=TRUE)[1,1]
            .s <- ifel(.s >= .q,.q,.s)
            .s
          }
)


setMethod('sed', signature(x='RasterStackBrick'),
          function(x,...,t1,t2) {
            xx <- list(x,...)
            if (missing(t1) || missing(t2)) stop("t1 and/or t2 (layers' indicators corresponding to time1 and time2) are not provided!")
            
            if (!is.numeric(t1) || !is.numeric(t2)) stop("t1 and t2 (layers' indicators corresponding to time1 and time2) should be a numeric vector")
            
            if (.require("terra")) {
              for (i in 1:length(xx)) {
                xx[[i]] <- rast(xx[[i]])
              }
              .sed0(xx,t1=t1,t2=t2)
            } else {
              .sedRas(xx,t1=t1,t2=t2)
            }
          }
)
#---------


setMethod('sed', signature(x='missing'),
          function(x,...,t1,t2) {
            if (missing(t1) && missing(t2)) stop('No Raster Objects are provided!')
            if (!is.list(t1) && !is.list(t2)) stop('The input Raster objects should be provided as time series, or as separate lists in t1 and t2')
            if (length(t1) != length(t2)) stop('The lists provided in t1 and t2 has not the same length!')
            
            cls <- unique(c(sapply(t1,function(x) class(x)),sapply(t2,function(x) class(x))))
            
            
            if (length(cls) > 1) stop('The classes of objects in the input lists (t1 and t2) are not the same!')
            
            if (!cls %in% c('SpatRasterTS','SpatRaster','RasterStack','RasterBrick','RasterStackTS','RasterBrickTS')) stop('The input objects provided in t1 and t2 lists are unknown (should be either of SpatRaster, RasterStack, RasterBrick, SpatRasterTS, RasterStackTS, RasterBrickTS!')
            
            if (cls == 'SpatRaster') .sedList(x1=t1,x2=t2)
            else if (cls == 'SpatRasterTS') {
              for (i in 1:length(t1)) {
                t1[[i]] <- t1[[i]]@raster
                t2[[i]] <- t2[[i]]@raster
              }
              .sedList(x1=t1,x2=t2)
            } else {
              if (.require('terra')) {
                if (cls %in% c('RasterStack','RasterBrick')) {
                  t1[[i]] <- rast(t1[[i]])
                  t2[[i]] <- rast(t2[[i]])
                } else {
                  t1[[i]] <- rast(t1[[i]]@raster)
                  t2[[i]] <- rast(t2[[i]]@raster)
                }
                
                raster(.sedList(x1=t1,x2=t2))
              } else {
                if (cls %in% c('RasterStack','RasterBrick')) {
                  .sedListR(t1,t2)
                } else {
                  for (i in 1:length(t1)) {
                    t1[[i]] <- t1[[i]]@raster
                    t2[[i]] <- t2[[i]]@raster
                  }
                  .sedListR(t1,t2)
                }
              }
            }
          }
)



# .sed <- function(p=NULL,tmin=NULL,tmax=NULL,tmean=NULL,t1,t2) {
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
#     if (i == 1) .s <- (mn1 - mn2) ^ 2 / v
#     else .s <- .s + ((mn1 - mn2) ^ 2 / v)
#   }
#   
#   .s <- sqrt(.s)
#   .s <- ifel(is.infinite(.s),0,.s)
#   .s
# }

