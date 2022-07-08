# Authors: Shirin Taheri (taheri.shi@gmail.com); Babak Naimi (naimi.b@gmail.com)
# Date :  Nov. 2020
# Last update :  July 2022
# Version 2.4
# Licence GPL v3
#--------

# for a single climate time series periods, available in time 1 (x1) and time 2 (x2),
# quantifies the probability of extreme value of variable over T1 (prob: quantile of ), in both T1 and T2:
.getExtremeProb <- function(x1,x2,prob) {
  .tail <- prob < 0.5
  q <- quantile(x1,prob=prob,na.rm=TRUE)
  m <- mean(x2,na.rm=TRUE)
  s <- app(x2,'sd',na.rm=TRUE)
  p2 <- app(c(q,m,s),function(x,...) {
    if (all(!is.na(x))) {
      pnorm(x[1],x[2],x[3],lower.tail = .tail)
    } else NA
  })
  #-------
  m <- mean(x1,na.rm=TRUE)
  s <- app(x1,'sd',na.rm=TRUE)
  p1 <- app(c(q,m,s),function(x,...) {
    if (all(!is.na(x))) {
      pnorm(x[1],x[2],x[3],lower.tail = .tail)
    } else NA
  })
  list(p1,p2)
}
#---------
.getExtremeProbR <- function(x1,x2,prob) {
  .tail <- prob < 0.5
  q <- calc(x1,fun=quantile,prob=prob,na.rm=TRUE)
  m <- mean(x2,na.rm=TRUE)
  s <- calc(x2,fun=sd,na.rm=TRUE)
  p2 <- calc(stack(q,m,s),function(x,...) {
    if (all(!is.na(x))) {
      pnorm(x[1],x[2],x[3],lower.tail = .tail)
    } else NA
  })
  #-------
  m <- mean(x1,na.rm=TRUE)
  s <- calc(x1,fun=sd,na.rm=TRUE)
  p1 <- calc(stack(q,m,s),function(x,...) {
    if (all(!is.na(x))) {
      pnorm(x[1],x[2],x[3],lower.tail = .tail)
    } else NA
  })
  list(p1,p2)
}
#---------

.eeChange <- function(x1,x2, extreme) {
  # x1 and x2 should be a list of climate variables correspond to time 1 and time 2
  # can be with length of 1 or 2 (one or two climate variables)
  o <- list()
  for (i in 1:length(x1)) {
    o[[i]] <- .getExtremeProb(x1[[i]],x2[[i]],prob = extreme[i])
  }
  
  if (length(o) == 2) {
    p1 <- o[[1]][[1]] + o[[2]][[1]] - (o[[1]][[1]] * o[[2]][[1]]) 
    p2 <- o[[1]][[2]] + o[[2]][[2]] - (o[[1]][[2]] * o[[2]][[2]]) 
    p2 - p1
  } else {
    o[[1]][[2]] - o[[1]][[1]]
  }
}
#-----

.eeChangeR <- function(x1,x2, extreme) {
  # x1 and x2 should be a list of climate variables correspond to time 1 and time 2
  # can be with length of 1 or 2 (one or two climate variables)
  o <- list()
  for (i in 1:length(x1)) {
    o[[i]] <- .getExtremeProbR(x1[[i]],x2[[i]],prob = extreme[i])
  }
  
  if (length(o) == 2) {
    p1 <- o[[1]][[1]] + o[[2]][[1]] - (o[[1]][[1]] * o[[2]][[1]]) 
    p2 <- o[[1]][[2]] + o[[2]][[2]] - (o[[1]][[2]] * o[[2]][[2]]) 
    p2 - p1
  } else {
    o[[1]][[2]] - o[[1]][[1]]
  }
}
#--------
# .eeChange <- function(xx, t1, t2, extreme) {
#   # local extreme event change
#   # extreme for temperature (>= 0.95), for precipitation ( <= (1-0.95)=0.05)
#   x1.t1 <- xx[[1]][[t1]]
#   x1.t2 <- xx[[1]][[t2]]
#   x2.t1 <- xx[[2]][[t1]]
#   x2.t2 <- xx[[2]][[t2]]
#   #------
#   .ext1 <- app(x1.t1@raster,function(x) quantile(x,extreme[1],na.rm=TRUE))
#   
#   .pt <- app(c(x1.t2@raster,.ext1),function(x,...) {
#     x <- x[!is.na(x)]
#     l <- length(x)
#     if (l > 4) {
#       .ex <- x[l]
#       x <- x[1:c(l-1)]
#       1 - (length(which(x >= .ex)) / length(x))
#     } else NA
#   })
#   #--------
#   .ext2 <- app(x2.t1@raster,function(x) quantile(x,extreme[2],na.rm=TRUE))
#   
#   .pp <- app(c(x2.t2@raster,.ext2),function(x,...) {
#     x <- x[!is.na(x)]
#     l <- length(x)
#     if (l > 5) {
#       .ex <- x[l]
#       x <- x[1:c(l-1)]
#       length(which(x <= .ex)) / length(x)
#     } else NA
#   })
#   
#   extreme[1] - ((.pt + .pp) - (.pt * .pp))
# }
#-------

# 
# .eeChangeR <- function(xx, t1, t2, extreme) {
#   # local extreme event change
#   # extreme for temperature (>= 0.95), for precipitation ( <= (1-0.95)=0.05)
#   x1.t1 <- xx[[1]][[t1]]
#   x1.t2 <- xx[[1]][[t2]]
#   x2.t1 <- xx[[2]][[t1]]
#   x2.t2 <- xx[[2]][[t2]]
#   #------
#   .ext1 <- app(x1.t1@raster,function(x) quantile(x,extreme[1],na.rm=TRUE))
#   
#   .pt <- app(c(x1.t2@raster,.ext1),function(x,...) {
#     x <- x[!is.na(x)]
#     l <- length(x)
#     if (l > 5) {
#       .ex <- x[l]
#       x <- x[1:c(l-1)]
#       1 - (length(which(x >= .ex)) / length(x))
#     } else NA
#   })
#   #--------
#   .ext2 <- app(x2.t1,function(x) quantile(x,extreme[2],na.rm=TRUE))
#   
#   .pp <- app(c(x2.t2,.ext2),function(x,...) {
#     x <- x[!is.na(x)]
#     l <- length(x)
#     if (l > 5) {
#       .ex <- x[l]
#       x <- x[1:c(l-1)]
#       length(which(x <= .ex)) / length(x)
#     } else NA
#   })
#   
#   extreme[1] - ((.pt + .pp) - (.pt * .pp))
# }



###############
if (!isGeneric("localExtreme")) {
  setGeneric("localExtreme", function(x1,x2,t1,t2,extreme)
    standardGeneric("localExtreme"))
}


setMethod('localExtreme', signature(x1='SpatRasterTS'),
          function(x1,x2,t1,t2,extreme) {
            
            if (missing(t1) || missing(t2)) stop("t1 and t2 (layers' indicators corresponding to time1 and time2) are not provided!")
            
            if (missing(x2) || is.null(x2)) xx <- list(x1)
            else xx <- list(x1,x2)
            
            
            if (missing(extreme)) stop('"extreme" is needed...!')
            #-
            if (length(xx) > 2) {
              xx <- xx[1:2]
              warning('localExtreme metric can be calculated using either one or two climate variables... only the first two variables are used!')
            }
            #-
            if (length(extreme) > length(xx)) {
              warning(paste0('number of extreme probabilities (',length(extreme),') should be equal to the number of climate parameters (',length(xx),'); Only the first ',length(xx),' values from extreme is used!'))
              extreme <- extreme[1:length(xx)]
            }
            #-
            if (length(extreme) < length(xx)) stop('extreme probabilities should be provided for both climate parameters')
            #-
            if (any(extreme > 1 | extreme < 0)) stop('extreme is a probability value that should be within the range of 0 and 1')  
            #-------------------------
            x1 <- x2 <- list()
            for (i in 1:length(xx)) {
              x1[[i]] <- xx[[i]][[t1]]@raster
              x2[[i]] <- xx[[i]][[t2]]@raster
            }
            
            .eeChange(x1,x2,extreme)
          }
)
#------------

setMethod('localExtreme', signature(x1='SpatRaster'),
          function(x1,x2,t1,t2,extreme) {
            
            if (missing(t1) || missing(t2)) stop("t1 and t2 (layers' indicators corresponding to time1 and time2) are not provided!")
            
            if (!is.numeric(t1) || !is.numeric(t2)) stop("t1 and t2 (layers' indicators corresponding to time1 and time2) should be a numeric vector")
            
            if (missing(x2) || is.null(x2)) xx <- list(x1)
            else xx <- list(x1,x2)
            
            
            if (missing(extreme)) stop('"extreme" is needed...!')
            #-
            if (length(xx) > 2) {
              xx <- xx[1:2]
              warning('localExtreme metric can be calculated using either one or two climate variables... only the first two variables are used!')
            }
            #-
            if (length(extreme) > length(xx)) {
              warning(paste0('number of extreme probabilities (',length(extreme),') should be equal to the number of climate parameters (',length(xx),'); Only the first ',length(xx),' values from extreme is used!'))
              extreme <- extreme[1:length(xx)]
            }
            #-
            if (length(extreme) < length(xx)) stop('extreme probabilities should be provided for both climate parameters')
            #-
            if (any(extreme > 1 | extreme < 0)) stop('extreme is a probability value that should be within the range of 0 and 1')  
            #-------------------------
            x1 <- x2 <- list()
            for (i in 1:length(xx)) {
              x1[[i]] <- xx[[i]][[t1]]
              x2[[i]] <- xx[[i]][[t2]]
            }
            
            .eeChange(x1,x2,extreme)
          }
)
#--------

setMethod('localExtreme', signature(x1='RasterStackBrickTS'),
          function(x1,x2,t1,t2,extreme) {
            
            if (missing(t1) || missing(t2)) stop("t1 and t2 (layers' indicators corresponding to time1 and time2) are not provided!")
            
            if (missing(x2) || is.null(x2)) xx <- list(x1)
            else xx <- list(x1,x2)
            
            
            if (missing(extreme)) stop('"extreme" is needed...!')
            #-
            if (length(xx) > 2) {
              xx <- xx[1:2]
              warning('localExtreme metric can be calculated using either one or two climate variables... only the first two variables are used!')
            }
            #-
            if (length(extreme) > length(xx)) {
              warning(paste0('number of extreme probabilities (',length(extreme),') should be equal to the number of climate parameters (',length(xx),'); Only the first ',length(xx),' values from extreme is used!'))
              extreme <- extreme[1:length(xx)]
            }
            #-
            if (length(extreme) < length(xx)) stop('extreme probabilities should be provided for both climate parameters')
            #-
            if (any(extreme > 1 | extreme < 0)) stop('extreme is a probability value that should be within the range of 0 and 1')  
            #-------------------------
            
            
            if (.require('terra')) {
              for (i in 1:length(xx)) {
                xx[[i]] <- rts(rast(xx[[i]]@raster),index(xx[[i]]@time))
              }
              x1 <- x2 <- list()
              for (i in 1:length(xx)) {
                x1[[i]] <- xx[[i]][[t1]]@raster
                x2[[i]] <- xx[[i]][[t2]]@raster
              }
              
              raster(.eeChange(x1,x2,extreme))
            } else {
              x1 <- x2 <- list()
              for (i in 1:length(xx)) {
                x1[[i]] <- xx[[i]][[t1]]@raster
                x2[[i]] <- xx[[i]][[t2]]@raster
              }
              .eeChangeR(x1,x2,extreme) 
            }
          }
)
#------------

setMethod('localExtreme', signature(x1='RasterStackBrick'),
          function(x1,x2,t1,t2,extreme) {
            
            if (missing(t1) || missing(t2)) stop("t1 and t2 (layers' indicators corresponding to time1 and time2) are not provided!")
            
            if (!is.numeric(t1) || !is.numeric(t2)) stop("t1 and t2 (layers' indicators corresponding to time1 and time2) should be a numeric vector")
            
            if (missing(x2) || is.null(x2)) xx <- list(x1)
            else xx <- list(x1,x2)
            
            
            if (missing(extreme)) stop('"extreme" is needed...!')
            #-
            if (length(xx) > 2) {
              xx <- xx[1:2]
              warning('localExtreme metric can be calculated using either one or two climate variables... only the first two variables are used!')
            }
            #-
            if (length(extreme) > length(xx)) {
              warning(paste0('number of extreme probabilities (',length(extreme),') should be equal to the number of climate parameters (',length(xx),'); Only the first ',length(xx),' values from extreme is used!'))
              extreme <- extreme[1:length(xx)]
            }
            #-
            if (length(extreme) < length(xx)) stop('extreme probabilities should be provided for both climate parameters')
            #-
            if (any(extreme > 1 | extreme < 0)) stop('extreme is a probability value that should be within the range of 0 and 1')  
            #-------------------------
            
            
            if (.require('terra')) {
              for (i in 1:length(xx)) {
                xx[[i]] <- rast(xx[[i]])
              }
              x1 <- x2 <- list()
              for (i in 1:length(xx)) {
                x1[[i]] <- xx[[i]][[t1]]
                x2[[i]] <- xx[[i]][[t2]]
              }
              
              raster(.eeChange(x1,x2,extreme))
            } else {
              x1 <- x2 <- list()
              for (i in 1:length(xx)) {
                x1[[i]] <- xx[[i]][[t1]]
                x2[[i]] <- xx[[i]][[t2]]
              }
              .eeChangeR(x1,x2,extreme) 
            }
          }
)
#------------

setMethod('localExtreme', signature(x1='list',x2='list'),
          function(x1,x2,t1,t2,extreme) {
            
            cls <- unique(c(sapply(x1,function(x) class(x)),sapply(x2,function(x) class(x))))
            
            if (length(cls) > 1) stop('The classes of objects in the input lists (x1 and x2) are not the same!')
            
            if (length(x1) != length(x2)) stop('x1 and x2 should have the same length (should be a list with either one or two time series of climate variables)!')
            
            if (length(x1) > 4)  stop('x1 and x2 should be a list with either one or two time series of climate variables!')
            
            if (length(x1) > 2) {
              x1 <- x1[1:2]
              x2 <- x2[1:2]
              warning('localExtreme metric can be calculated using either one or two climate variables... only the first two time series variables in the list are used!')
            }
            #######
            
            if (missing(extreme)) stop('"extreme" is needed...!')
            #-
            if (length(extreme) > length(x1)) {
              warning(paste0('number of extreme probabilities (',length(extreme),') should be equal to the number of climate parameters (',length(x1),'); Only the first ',length(x1),' values from extreme is used!'))
              extreme <- extreme[1:length(x1)]
            }
            #-
            if (length(extreme) < length(x1)) stop('extreme probabilities should be provided for both climate parameters')
            #-
            if (any(extreme > 1 | extreme < 0)) stop('extreme is a probability value that should be within the range of 0 and 1')  
            #-------------------------
            if (cls == 'SpatRasterTS') {
              for (i in 1:length(x1)) {
                x1[[i]] <- x1[[i]]@raster
                x2[[i]] <- x2[[i]]@raster
              }
              .eeChange(x1,x2,extreme) 
            } else if (cls == 'SpatRaster') {
              .eeChange(x1,x2,extreme) 
            } else if (cls == 'RasterStackBrickTS') {
              
              if (.require('terra')) {
                for (i in 1:length(x1)) {
                  x1[[i]] <- rast(x1[[i]]@raster)
                  x2[[i]] <- rast(x2[[i]]@raster)
                }
                raster(.eeChange(x1,x2,extreme))
              } else {
                for (i in 1:length(x1)) {
                  x1[[i]] <- x1[[i]]@raster
                  x2[[i]] <- x2[[i]]@raster
                }
                .eeChangeR(x1,x2,extreme)
              }
              
            } else if (cls == 'RasterStackBrick') {
              if (.require('terra')) {
                for (i in 1:length(x1)) {
                  x1[[i]] <- rast(x1[[i]])
                  x2[[i]] <- rast(x2[[i]])
                }
                raster(.eeChange(x1,x2,extreme))
              } else {
                .eeChangeR(x1,x2,extreme)
              }
            }
            
          }
)
#--------
setMethod('localExtreme', signature(x1='missing',x2='missing',t1='list',t2='list'),
          function(x1,x2,t1,t2,extreme) {
            
            cls <- unique(c(sapply(t1,function(x) class(x)),sapply(t2,function(x) class(x))))
            
            if (length(cls) > 1) stop('The classes of objects in the input lists (t1 and t2) are not the same!')
            
            if (length(t1) != length(t2)) stop('t1 and t2 should have the same length (should be a list with either one or two time series of climate variables)!')
            
            if (length(t1) > 4)  stop('x1 and x2 should be a list with either one or two time series of climate variables!')
            
            if (length(t1) > 2) {
              t1 <- t1[[1:2]]
              t2 <- t2[[1:2]]
              warning('localExtreme metric can be calculated using either one or two climate variables... only the first two time series variables in the list are used!')
            }
            #######
            
            x1 <- t1
            x2 <- t2
            
            if (missing(extreme)) stop('"extreme" is needed...!')
            #-
            if (length(extreme) > length(x1)) {
              warning(paste0('number of extreme probabilities (',length(extreme),') should be equal to the number of climate parameters (',length(x1),'); Only the first ',length(x1),' values from extreme is used!'))
              extreme <- extreme[1:length(x1)]
            }
            #-
            if (length(extreme) < length(x1)) stop('extreme probabilities should be provided for both climate parameters')
            #-
            if (any(extreme > 1 | extreme < 0)) stop('extreme is a probability value that should be within the range of 0 and 1')  
            #-------------------------
            if (cls == 'SpatRasterTS') {
              for (i in 1:length(x1)) {
                x1[[i]] <- x1[[i]]@raster
                x2[[i]] <- x2[[i]]@raster
              }
              .eeChange(x1,x2,extreme) 
            } else if (cls == 'SpatRaster') {
              .eeChange(x1,x2,extreme) 
            } else if (cls == 'RasterStackBrickTS') {
              
              if (.require('terra')) {
                for (i in 1:length(x1)) {
                  x1[[i]] <- rast(x1[[i]]@raster)
                  x2[[i]] <- rast(x2[[i]]@raster)
                }
                raster(.eeChange(x1,x2,extreme))
              } else {
                for (i in 1:length(x1)) {
                  x1[[i]] <- x1[[i]]@raster
                  x2[[i]] <- x2[[i]]@raster
                }
                .eeChangeR(x1,x2,extreme)
              }
              
            } else if (cls == 'RasterStackBrick') {
              if (.require('terra')) {
                for (i in 1:length(x1)) {
                  x1[[i]] <- rast(x1[[i]])
                  x2[[i]] <- rast(x2[[i]])
                }
                raster(.eeChange(x1,x2,extreme))
              } else {
                .eeChangeR(x1,x2,extreme)
              }
            }
            
          }
)
#--------