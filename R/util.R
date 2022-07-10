# Authors: Shirin Taheri, taheri.shi@gmail.com; Babak Naimi (naimi.b@gmail.com)
# Date :  Nov. 2020
# Last update :  July 2022
# Version 2.0
# Licence GPL v3
#--------



.getProj <- function(x) {
  if (inherits(x,'Raster')) {
    if (!is.na(projection(x))) strsplit(strsplit(projection(x),'\\+proj=')[[1]][2],' ')[[1]][1]
    else {
      if (all(extent(x)[1:2] >= -180 & extent(x)[1:2] <= 180 & extent(x)[3:4] >= -90 & extent(x)[3:4] <= 90)) 'longlat'
      else 'projected'
    }
  } else {
    if (!is.na(crs(x))) strsplit(strsplit(crs(x,proj=TRUE),'\\+proj=')[[1]][2],' ')[[1]][1]
    else {
      if (all(extent(x)[1:2] >= -180 & extent(x)[1:2] <= 180 & extent(x)[3:4] >= -90 & extent(x)[3:4] <= 90)) 'longlat'
      else 'projected'
    }
  }
  
}
#=-===============
.is.projected <- function(x) {
  if (inherits(x,'SpatRaster')) {
    e <- as.vector(terra::ext(x))
  } else e <- as.vector(extent(x))
  
  !all(e >= -180 & e <= 180)
}
#-----

.rad <- function (degree) {
  (degree * pi) / 180
}
#---
.deg <-  function (radian) {
  (radian * 180) / pi
}



# following two functions are copied from the package vocc in "https://github.com/cbrown5/vocc"
.mnwm <- function(d1, d2, d3, d4, d5, d6){
  X <- sum(c(d1, d2*2, d3, d4, d5*2, d6), na.rm = T)
  w <- sum(c(1,2,1,1,2,1) * is.finite(c(d1, d2, d3, d4, d5, d6)))
  return(X/w)
}
#-----
.ang <- function(dx, dy){
  ifelse(dy < 0, 180 + .deg(atan(dx/dy)),
         ifelse(dx < 0, 360 + .deg(atan(dx /dy )), .deg(atan(dx/dy))))
}
#--------
######################



# this function stretch each variable between 0 and 1 (in T1, also in T2 but based on parameters of T1), 
# and then the multiple variable are averaged to have a single variable for each time.
# x1 and x2 is the list of Raster(s) in T1 and T2 (if multiple variables, then multiple raster in each list)"
# The raster of each variable in x1 and x2 should be a single layer (mean of Time series)
.getScaledMultiVariateIntoOne <- function(x1,x2) {
  
  if (length(x1) != length(x2)) stop('x1 and x2 should have the same number of variables!')
  
  .scaleParam <- vector(mode='list',length=length(x1))
  names(.scaleParam) <- sapply(x1,names)
  
  
  if (inherits(x1[[1]],'Raster')) {
    for (i in 1:length(x1)) {
      .scaleParam[[i]] <- cellStats(x1[[1]],'range',na.rm=TRUE)
      x1[[i]] <- (x1[[i]] - .scaleParam[[i]][1]) / (.scaleParam[[i]][2] - .scaleParam[[i]][1])
      x2[[i]] <- (x2[[i]] - .scaleParam[[i]][1]) / (.scaleParam[[i]][2] - .scaleParam[[i]][1])
    }
  } else {
    for (i in 1:length(x1)) {
      .scaleParam[[i]] <- as.matrix(global(x1[[1]],'range',na.rm=TRUE))[1,]
      x1[[i]] <- (x1[[i]] - .scaleParam[[i]][1]) / (.scaleParam[[i]][2] - .scaleParam[[i]][1])
      x2[[i]] <- (x2[[i]] - .scaleParam[[i]][1]) / (.scaleParam[[i]][2] - .scaleParam[[i]][1])
    }
  }
  #-----
  .x1 <- x1[[1]]
  .x2 <- x2[[1]]
  
  if (length(x1) > 1) {
    for (i in 2:length(x1)) {
      .x1 <- .x1 + x1[[i]]
      .x2 <- .x2 + x2[[i]]
    }
    .x1 <- .x1 / length(x1)
    .x2 <- .x2 / length(x1)
  }
  list(t1=.x1,t2=.x2)
}
#----
.getNyears <- function(x,t1,t2) {
  if (missing(t1) | missing(t2)) {
    nyears(x)
  } else {
    y1 <- mean(index(x[t1]))
    y2 <- mean(index(x[t2]))
    as.numeric(format(y2,'%Y')) - as.numeric(format(y1,'%Y'))
  }
}



# .getScaledMultiVariateIntoOne <- function(x,...,t1,t2) {
#   lst <- list(x,...)
#   xt1 <- xt2 <- list()
#   
#   if (length(t1) > 1 | length(t2) > 1) .multi <- TRUE
#   else .multi <- FALSE
#   
#   if (inherits(x,'Raster')) .raster <- TRUE
#   else .raster <- FALSE
#   
#   
#   
#   for (i in 1:length(lst)) {
#     if (.multi) {
#       if (.raster) {
#         xt1[[i]] <- mean(rast(lst[[i]])[[t1]])
#         xt2[[i]] <- mean(rast(lst[[i]])[[t2]])
#       } else {
#         xt1[[i]] <- mean(lst[[i]][[t1]])
#         xt2[[i]] <- mean(lst[[i]][[t2]])
#       }
#     } else {
#       if (.raster) {
#         xt1[[i]] <- rast(lst[[i]])[[t1]]
#         xt2[[i]] <- rast(lst[[i]])[[t2]]
#       } else {
#         xt1[[i]] <- lst[[i]][[t1]]
#         xt2[[i]] <- lst[[i]][[t2]]
#       }
#     }
#   }
#   #----
#   
#   .scaleParam <- vector(mode='list',length=length(xt1))
#   names(.scaleParam) <- sapply(xt1,names)
#   
#   
#   for (i in 1:length(xt1)) {                                                                                                                                                                                                                                                                                      
#     .scaleParam[[i]] <- minmax(xt1[[i]])[,1]
#     xt1[[i]] <- (xt1[[i]] - .scaleParam[[i]][1]) / (.scaleParam[[i]][2] - .scaleParam[[i]][1] )
#     xt2[[i]] <- (xt2[[i]] - .scaleParam[[i]][1]) / (.scaleParam[[i]][2] - .scaleParam[[i]][1])
#   }
#   #-----
#   x1 <- xt1[[1]]
#   x2 <- xt2[[1]]
#   
#   if (length(xt1) > 1) {
#     for (i in 2:length(xt1)) {
#       x1 <- x1 + xt1[[i]]
#       x2 <- x2 + xt2[[i]]
#     }
#     x1 <- x1 / length(xt1)
#     x2 <- x2 / length(xt1)
#   }
#   list(t1=x1,t2=x2)
# }

#----------
.getLongLat <- function(x,crs) {
  # x is data.frame
  if (ncol(x) > 2) {
    if (all(c('x','y') %in% colnames(x))) x <- x[,c('x','y')]
    else {
      x <- x[,1:2]
      colnames(x) <- c('x','y')
      warning('".getLongLat function:" the first two columns in the input data.frame is assumed as the coordinate columns!')
    }
  } else if (ncol(x) < 2) stop('".getLongLat function:" input x should be a data.frame with 2 columns, i.e., x and y coordinates!')
  
  
  geom(project(vect(x,geom=colnames(x),crs=crs),"epsg:4326"))[,c('x','y')]
}
#-----
.agrep <- function(n,choices, r=seq(0,0.3,0.05)) {
  # r is a range can be used for max distance
  for (i in r) {
    w <- agrep(n,choices,ignore.case = TRUE,max.distance = list(all=i))
    if (length(w) > 0) break
  }
  if (length(w) > 1) {
    d <- unlist(lapply(choices[w],function(x) .LD(n,x)))
    w2 <- which(d == min(d))
    if (length(w2) == 1) choices[w][w2]
    else NA
  } else if (length(w) == 1) choices[w]
  else NA
}
#----
.LD <- function(s,t) {
  sl <- unlist(strsplit(s,''))
  tl <- unlist(strsplit(t,''))
  if (s == t) return(0)
  else if (length(sl) == 0) return(length(tl))
  else if (length(tl) == 0) return(length(sl))
  v0 <- 0:length(tl)
  v1 <- rep(NA,length(tl)+1)
  for (i in seq_along(sl)) {
    v1[1] <- i
    for (j in seq_along(tl)) {
      if (sl[i] == tl[j]) cost <- 0
      else cost <- 1
      v1[j+1] <- min(v1[j] + 1, v0[j + 1] + 1, v0[j] + cost)
    }
    for (j in seq_along(v0)) {
      v0[j] <- v1[j]
    }
  }
  return(v1[length(tl)+1])
}
#---
.argMatch <- function(n,choices) {
  for (i in seq_along(n)) {
    if (n[i] != '') {
      if (!n[i] %in% choices) {
        u <- try(match.arg(tolower(n[i]),tolower(choices)),silent=TRUE)
        if (!inherits(u,"try-error")) {
          n[i] <- choices[which(tolower(choices) == u)]
        } else {
          n[i] <- .agrep(n[i],choices)
        }
      }
    } else n[i] <- NA
  }
  n
}
#------
.require <- function(x) {
  x <- as.character(x)
  xx <- unlist(lapply(.libPaths(), function(lib) find.package(x, lib, quiet=TRUE, verbose=FALSE)))
  if (length(xx) > 0) {
    .loaded <- eval(parse(text=paste0('require(',x,')')))
    return (.loaded)
  } else FALSE
}
#----
.is_package_installed <- function(n) {
  names(n) <- n
  sapply(n, function(x) length(unlist(lapply(.libPaths(), function(lib) find.package(x, lib, quiet=TRUE, verbose=FALSE)))) > 0)
}
