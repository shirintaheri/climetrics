# Authors: Shirin Taheri (taheri.shi@gmail.com); Babak Naimi (naimi.b@gmail.com)
# Date :  Sep. 2021
# Last update :  April 2022
# Version 1.3
# Licence GPL v3
#--------


# temporal gradients using time series:
# threshold of > 5 obs. are considered to get slope 
.tempgradFun <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) > 1) {
    s <- lm(x~c(1:length(x)))
    s$coefficients[2]*10
  } else NA
}


# x: RasterStackBrick obect:
.tempgradRaster <- function(x) {
  calc(x,.tempgradFun)
}
#--------------
.tempgradTerra <- function(x) {
  app(x,.tempgradFun)
}
#--------------

.spatialgradTerra <- function(rx, y_diff = 1) {
  # based on the function provided in the vocc package (https://github.com/cbrown5/vocc)
  if (.getProj(rx) == 'longlat') y_dist <- res(rx) * c(111.325, 111.325)
  else {
    y_dist <- res(rx) / 1000
    y_diff <- NA
  }
  
  nlats <- nrow(rx)
  nlons <- ncol(rx)
  y <- data.frame(adjacent(rx, cells=1:ncell(rx), directions=8,pairs=TRUE))
  y <- y[order(y$from, y$to),]
  y <- na.omit(y)
  y$sst <- rx[y$to][,1]
  y$sy <- rowFromCell(rx, y$from)-rowFromCell(rx, y$to)
  y$sx <- colFromCell(rx, y$to)-colFromCell(rx, y$from)
  y$sx[y$sx > 1] <- -1
  y$sx[y$sx < -1] <- 1
  y$code <- paste(y$sx, y$sy)
  
  y$code1 <- dplyr::recode(y$code,
                           `1 0` = "sstE",
                           `-1 0` = "sstW",
                           `-1 1` = "sstNW",
                           `-1 -1` = "sstSW",
                           `1 1` = "sstNE",
                           `1 -1` = "sstSE",
                           `0 1` = "sstN",
                           `0 -1` = "sstS")
  
  y3b <- eval(parse(text="dplyr::select(y,from, code1, sst)"),envir =environment())
  y3b <- eval(parse(text="tidyr::spread(y3b,code1, sst)"),envir =environment())
  y3b$sstFocal <- rx[y3b$from][,1]
  y3b$LAT <- yFromCell(rx, y3b$from)
  
  if(!is.na(y_diff)) {
    y3b <- eval(parse(text="dplyr::mutate(y3b,
                         latpos = cos(.rad(LAT + y_diff)),
                         latneg = cos(.rad(LAT - y_diff)),
                         latfocal = cos(.rad(LAT)))"),envir =environment())
  } else {
    y3b <- dplyr::mutate(y3b,
                         latpos = 1,
                         latneg = 1,
                         latfocal = 1)
  }
  
  y3c <- "dplyr::mutate(y3b,
                       gradWE1 = (sstN-sstNW)/
                         (latpos *  y_dist[1]),
                       gradWE2 = (sstFocal - sstW)/(latfocal * y_dist[1]),
                       gradWE3 = (sstS-sstSW)/(latneg * y_dist[1]),
                       gradWE4 = (sstNE-sstN)/(latpos * y_dist[1]),
                       gradWE5 = (sstE-sstFocal)/(latfocal * y_dist[1]),
                       gradWE6 = (sstSE-sstS)/(latneg*y_dist[1]),
                       gradNS1 = (sstNW-sstW)/y_dist[2],
                       gradNS2 = (sstN-sstFocal)/y_dist[2],
                       gradNS3 = (sstNE-sstE)/y_dist[2],
                       gradNS4 = (sstW-sstSW)/y_dist[2],
                       gradNS5 = (sstFocal-sstS)/y_dist[2],
                       gradNS6 = (sstE-sstSE)/y_dist[2])" 
  
  y3c <- eval(parse(text=y3c),envir=environment())
  
  y3c <- dplyr::rowwise(y3c)
  y3c <- eval(parse(text="dplyr::mutate(y3c,
      WEgrad = .mnwm(gradWE1, gradWE2, gradWE3, gradWE4, gradWE5, gradWE6),
      NSgrad = .mnwm(gradNS1, gradNS2, gradNS3, gradNS4, gradNS5, gradNS6),
      angle = .ang(WEgrad, NSgrad))"),envir=environment())
  
  y3c <- eval(parse(text="dplyr::select(y3c,icell = from, WE = WEgrad, NS = NSgrad, angle = angle)"),envir=environment())
  return(y3c)
}
#----------

.spatialgradRaster <- function(rx, y_diff = 1) {
  # based on the function provided in the vocc package (https://github.com/cbrown5/vocc)
  if (.getProj(rx) == 'longlat') y_dist <- res(rx) * c(111.325, 111.325)
  else {
    y_dist <- res(rx) / 1000
    y_diff <- NA
  }
  
  nlats <- nrow(rx)
  nlons <- ncol(rx)
  y <- data.frame(adjacent(rx, cells=1:ncell(rx), directions=8,pairs=TRUE))
  y <- y[order(y$from, y$to),]
  y <- na.omit(y)
  y$sst <- rx[y$to]
  y$sy <- rowFromCell(rx, y$from)-rowFromCell(rx, y$to)
  y$sx <- colFromCell(rx, y$to)-colFromCell(rx, y$from)
  y$sx[y$sx > 1] <- -1
  y$sx[y$sx < -1] <- 1
  y$code <- paste(y$sx, y$sy)
  
  y$code1 <- dplyr::recode(y$code,
                           `1 0` = "sstE",
                           `-1 0` = "sstW",
                           `-1 1` = "sstNW",
                           `-1 -1` = "sstSW",
                           `1 1` = "sstNE",
                           `1 -1` = "sstSE",
                           `0 1` = "sstN",
                           `0 -1` = "sstS")
  
  y3b <- eval(parse(text="dplyr::select(y,from, code1, sst)"),envir =environment())
  y3b <- eval(parse(text="tidyr::spread(y3b,code1, sst)"),envir =environment())
  y3b$sstFocal <- rx[y3b$from]
  y3b$LAT <- yFromCell(rx, y3b$from)
  
  if(!is.na(y_diff)) {
    y3b <- eval(parse(text="dplyr::mutate(y3b,
                         latpos = cos(.rad(LAT + y_diff)),
                         latneg = cos(.rad(LAT - y_diff)),
                         latfocal = cos(.rad(LAT)))"),envir =environment())
  } else {
    y3b <- dplyr::mutate(y3b,
                         latpos = 1,
                         latneg = 1,
                         latfocal = 1)
  }
  
  y3c <- "dplyr::mutate(y3b,
                       gradWE1 = (sstN-sstNW)/
                         (latpos *  y_dist[1]),
                       gradWE2 = (sstFocal - sstW)/(latfocal * y_dist[1]),
                       gradWE3 = (sstS-sstSW)/(latneg * y_dist[1]),
                       gradWE4 = (sstNE-sstN)/(latpos * y_dist[1]),
                       gradWE5 = (sstE-sstFocal)/(latfocal * y_dist[1]),
                       gradWE6 = (sstSE-sstS)/(latneg*y_dist[1]),
                       gradNS1 = (sstNW-sstW)/y_dist[2],
                       gradNS2 = (sstN-sstFocal)/y_dist[2],
                       gradNS3 = (sstNE-sstE)/y_dist[2],
                       gradNS4 = (sstW-sstSW)/y_dist[2],
                       gradNS5 = (sstFocal-sstS)/y_dist[2],
                       gradNS6 = (sstE-sstSE)/y_dist[2])" 
  
  y3c <- eval(parse(text=y3c),envir=environment())
  
  y3c <- dplyr::rowwise(y3c)
  y3c <- eval(parse(text="dplyr::mutate(y3c,
      WEgrad = .mnwm(gradWE1, gradWE2, gradWE3, gradWE4, gradWE5, gradWE6),
      NSgrad = .mnwm(gradNS1, gradNS2, gradNS3, gradNS4, gradNS5, gradNS6),
      angle = .ang(WEgrad, NSgrad))"),envir=environment())
  
  y3c <- eval(parse(text="dplyr::select(y3c,icell = from, WE = WEgrad, NS = NSgrad, angle = angle)"),envir=environment())
  
  
  return(y3c)
}
#----------
.calcvelocity <- function(grad, slope,.terra) {
  #slope$w <- y_dist * cos(.rad(slope$y))
  grad$NS[is.na(grad$NS)] <- 0
  grad$WE[is.na(grad$WE)] <- 0
  grad$NAsort <- ifelse((abs(grad$NS)+abs(grad$WE)) == 0, NA, 1)
  grad$Grad <- grad$NAsort * sqrt((grad$WE^2) + (grad$NS^2))
  
  if (.terra) {
    v <- rast(slope)
    
    v[grad$icell] <- slope[grad$icell] / grad$Grad
    
    .o <- as.matrix(global(v,fun=quantile,probs=c(0.05,0.95),na.rm=TRUE))[1,]
    v[v < .o[1]] <- .o[1]
    v[v > .o[2]] <- .o[2]
  } else {
    v <- raster(slope)
    
    v[grad$icell] <- slope[grad$icell] / grad$Grad
    
    
    .o <- quantile(v, probs=c(0.05,0.95),na.rm = TRUE)
    v[v < .o[1]] <- .o[1]
    v[v > .o[2]] <- .o[2]
  }
  
  return(v)
}
#--------

if (!isGeneric("temporalTrend")) {
  setGeneric("temporalTrend", function(x,...)
    standardGeneric("temporalTrend"))
}

setMethod('temporalTrend', signature(x='RasterStackBrickTS'),
          function(x,...) {
            .tempgradRaster(x@raster)
          }
)
#----
setMethod('temporalTrend', signature(x='SpatRasterTS'),
          function(x,...) {
            .tempgradTerra(x@raster)
          }
)
#----------
setMethod('temporalTrend', signature(x='RasterStackBrick'),
          function(x,...) {
            .tempgradRaster(x)
          }
)

#----------
setMethod('temporalTrend', signature(x='SpatRaster'),
          function(x,...) {
            .tempgradTerra(x)
          }
)

##################################


if (!isGeneric("gVelocity")) {
  setGeneric("gVelocity", function(x,...)
    standardGeneric("gVelocity"))
}

setMethod('gVelocity', signature(x='RasterStackBrickTS'),
          function(x,...) {
            xm <- mean(x@raster,na.rm=TRUE)
            g <- .spatialgradRaster(xm)
            s <- .tempgradRaster(x@raster)
            .calcvelocity(g,s,.terra = FALSE)
          }
)


setMethod('gVelocity', signature(x='RasterStackBrick'),
          function(x,...) {
            g <- .spatialgradRaster(x)
            s <- .tempgradRaster(x)
            .calcvelocity(g,s,.terra = FALSE)
          }
)
#-------

setMethod('gVelocity', signature(x='SpatRasterTS'),
          function(x,...) {
            g <- .spatialgradTerra(x@raster)
            s <- .tempgradTerra(x@raster)
            .calcvelocity(g,s,.terra=TRUE)
          }
)
#-------

setMethod('gVelocity', signature(x='SpatRaster'),
          function(x,...) {
            g <- .spatialgradTerra(x)
            s <- .tempgradTerra(x)
            .calcvelocity(g,s,.terra=TRUE)
          }
)
