# Authors: Shirin Taheri (taheri.shi@gmail.com); Babak Naimi (naimi.b@gmail.com)
# Date :  Nov. 2020
# Last update :  April 2024
# Version 2.3
# Licence GPL v3
#--------




.velocM <- function(p1,p2,f1,f2,log=TRUE) {
  # based on the script provided in Hamnan et al. (2015) -> Licence: "cc-by-nc"
  w <- which(!is.na(p1[]))
  
  present1 <- p1[w]
  present2 <- p2[w]
  future1  <- f1[w]
  future2  <- f2[w]
  
  idxy <- data.frame(cbind(id=1:length(present1),xyFromCell(p1,w)) )  # data frame of IDs and XY coords
  b <- (max(present1,na.rm=TRUE)-min(present1,na.rm=TRUE)) / 120  # bin size for 120 PC1 bins
  
  pr1 <- round(present1/b)              # convert PC1 to 120 bins via rounding
  pr2 <- round(present2/b)              # convert PC2 to <120 bins via rounding
  fu1 <- round(future1/b)               # same for future PC1
  fu2 <- round(future2/b)               # same for future PC2
  p  <- paste(pr1,pr2)                         # PC1/PC2 combinations in present climate
  f  <- paste(fu1,fu2)                         # PC1/PC2 combinations in future climate
  u  <- unique(p)[order(unique(p))]          # list of unique PC1/PC2 combinations
  
  sid <- c()                                 # empty vector for source IDs
  tid <- c()                                 # empty vector for target IDs
  d   <- c()                                 # empty vector for distances
  
  for(i in u){                          # loop for each unique PC1/PC2 combination
    pxy <- idxy[which(p==i),]           # coordinates of i-th combination in present
    fxy <- idxy[which(f==i),]           # coordinates of i-th combination in future
    sid <- c(sid, pxy$id)               # append i-th PC1/PC2 combination to previous 
    
    if(nrow(fxy)>0){                    # kNN search unless no-analogue climate
      knn <- data.frame(yaImpute::ann(as.matrix(fxy[,-1]), as.matrix(pxy[,-1]), k=1,verbose=FALSE)$knnIndexDist)      
      tid <- c(tid, fxy[knn[,1],"id"]) # the IDs of the closest matches  
      d <- c(d, sqrt(knn[,2]))         # their corresponding geographic distances
    } else {                              # else statement for no-analogue climates
      tid <- c(tid, rep(NA,nrow(pxy))) # flag destinations as missing for no analogues
      d <- c(d, rep(Inf,nrow(pxy)))    # flag distances as infinity for no analogues
    }
  }
  
  sxy <- merge(sid, idxy, by.y="id", all.x=TRUE, all.y=FALSE, sort=FALSE)[2:3]  # source coordinates
  txy <- merge(tid, idxy, by.y="id", all.x=TRUE, all.y=FALSE, sort=FALSE)[2:3]  # target coordinates
  names(txy)=c("target_x","target_y")
  names(sxy)=c("x","y")
  outtab <- cbind(id=sid, sxy, txy, distance=d)   
  # writes out log10 velocities and distances multiplied by 100 in ESRI ASCII format
  # conversion: -200=0.01km, -100=0.1km, 0=1km, 100=10km, 200=100km etc.
  out=merge(idxy, outtab[,c(2,3,6)], by=c("y","x"), sort=F)
  out$distance[out$distance==Inf] <- 10000  # sets no analogue to 10,000km
  out$distance[out$distance==0] <- 0.5  # sets zero distance to 0.5km (1/2 cell size)
  #out$logDist=round(log10(out$distance)*100)
  #out$logSpeed=round(log10(out$distance/50)*100)
  r <- raster(p1)
  if (log) r[cellFromXY(r,out[,c(2,1)])] <- round(log10(out$distance)*100)
  else r[cellFromXY(r,out[,c(2,1)])] <- out$distance
  r
}
#--------

.velocMTerra <- function(p1,p2,f1,f2,log=TRUE) {
  # based on the script provided in Hamnan et al. (2015) -> Licence: "cc-by-nc"
  w <- which(!is.na(p1[]))
  
  present1 <- p1[w][,1]
  present2 <- p2[w][,1]
  future1  <- f1[w][,1]
  future2  <- f2[w][,1]
  
  idxy <- data.frame(cbind(id=1:length(present1),xyFromCell(p1,w)) )  # data frame of IDs and XY coords
  b <- (max(present1,na.rm=TRUE)-min(present1,na.rm=TRUE)) / 120  # bin size for 120 PC1 bins
  
  pr1 <- round(present1/b)              # convert PC1 to 120 bins via rounding
  pr2 <- round(present2/b)              # convert PC2 to <120 bins via rounding
  fu1 <- round(future1/b)               # same for future PC1
  fu2 <- round(future2/b)               # same for future PC2
  p  <- paste(pr1,pr2)                         # PC1/PC2 combinations in present climate
  f  <- paste(fu1,fu2)                         # PC1/PC2 combinations in future climate
  u  <- unique(p)[order(unique(p))]          # list of unique PC1/PC2 combinations
  
  sid <- c()                                 # empty vector for source IDs
  tid <- c()                                 # empty vector for target IDs
  d   <- c()                                 # empty vector for distances
  
  for(i in u){                          # loop for each unique PC1/PC2 combination
    pxy <- idxy[which(p==i),]           # coordinates of i-th combination in present
    fxy <- idxy[which(f==i),]           # coordinates of i-th combination in future
    sid <- c(sid, pxy$id)               # append i-th PC1/PC2 combination to previous 
    
    if(nrow(fxy)>0){                    # kNN search unless no-analogue climate
      knn <- data.frame(yaImpute::ann(as.matrix(fxy[,-1]), as.matrix(pxy[,-1]), k=1,verbose=FALSE)$knnIndexDist)      
      tid <- c(tid, fxy[knn[,1],"id"]) # the IDs of the closest matches  
      d <- c(d, sqrt(knn[,2]))         # their corresponding geographic distances
    } else {                              # else statement for no-analogue climates
      tid <- c(tid, rep(NA,nrow(pxy))) # flag destinations as missing for no analogues
      d <- c(d, rep(Inf,nrow(pxy)))    # flag distances as infinity for no analogues
    }
  }
  
  sxy <- merge(sid, idxy, by.y="id", all.x=TRUE, all.y=FALSE, sort=FALSE)[2:3]  # source coordinates
  txy <- merge(tid, idxy, by.y="id", all.x=TRUE, all.y=FALSE, sort=FALSE)[2:3]  # target coordinates
  names(txy)=c("target_x","target_y")
  names(sxy)=c("x","y")
  outtab <- cbind(id=sid, sxy, txy, distance=d)
  # writes out log10 velocities and distances multiplied by 100 in ESRI ASCII format
  # conversion: -200=0.01km, -100=0.1km, 0=1km, 100=10km, 200=100km etc.
  out=merge(idxy, outtab[,c(2,3,6)], by=c("y","x"), sort=FALSE)
  out$distance[out$distance==Inf] <- 10000  # sets no analogue to 10,000km
  out$distance[out$distance==0] <- 0.5  # sets zero distance to 0.5km (1/2 cell size)
  
  r <- rast(p1)
  
  if (log) r[cellFromXY(r,out[,c(2,1)])] <- round(log10(out$distance)*100)
  else r[cellFromXY(r,out[,c(2,1)])] <- out$distance
  r
}
#-----------

if (!isGeneric("velocity")) {
  setGeneric("velocity", function(x1,x2,t1,t2,...)
    standardGeneric("velocity"))
}

setMethod('velocity', signature(x1='SpatRasterTS',x2='SpatRasterTS'),
          function(x1,x2,t1,t2,...) {
            
            if (missing(t1) || missing(t2)) stop("t1 and t2 (layers' indicators corresponding to time1 and time2) are not provided!")
            
            
            p1 <- app(x1[[t1]]@raster, 'mean',na.rm=TRUE)
            f1 <- app(x1[[t2]]@raster, 'mean',na.rm=TRUE)
            
            p2 <- app(x2[[t1]]@raster, 'mean',na.rm=TRUE)
            f2 <- app(x2[[t2]]@raster, 'mean',na.rm=TRUE)
            
            .velocMTerra(p1,p2,f1,f2,...)
            
          }
)
#------------
setMethod('velocity', signature(x1='SpatRaster',x2='SpatRaster'),
          function(x1,x2,t1,t2,...) {
            
            if (missing(t1) || missing(t2)) stop("t1 and t2 (layers' indicators corresponding to time1 and time2) are not provided!")
            
            if (!is.numeric(t1) || !is.numeric(t2)) stop("t1 and t2 (layers' indicators corresponding to time1 and time2) should be a numeric vector")
            
            p1 <- app(x1[[t1]], 'mean',na.rm=TRUE)
            f1 <- app(x1[[t2]], 'mean',na.rm=TRUE)
            
            p2 <- app(x2[[t1]], 'mean',na.rm=TRUE)
            f2 <- app(x2[[t2]], 'mean',na.rm=TRUE)
            
            .velocMTerra(p1,p2,f1,f2,...)
            
          }
)
#--------

setMethod('velocity', signature(x1='RasterStackBrickTS',x2='RasterStackBrickTS'),
          function(x1,x2,t1,t2,...) {
            
            if (missing(t1) || missing(t2)) stop("t1 and t2 (layers' indicators corresponding to time1 and time2) are not provided!")
            
            
            p1 <- mean(x1[[t1]]@raster, na.rm=TRUE)
            f1 <- mean(x1[[t2]]@raster, na.rm=TRUE)
            
            p2 <- mean(x2[[t1]]@raster,na.rm=TRUE)
            f2 <- mean(x2[[t2]]@raster,na.rm=TRUE)
            
            .velocM(p1,p2,f1,f2,...)
            
          }
)
#------------
setMethod('velocity', signature(x1='RasterStackBrick',x2='RasterStackBrick'),
          function(x1,x2,t1,t2,...) {
            
            if (missing(t1) || missing(t2)) stop("t1 and t2 (layers' indicators corresponding to time1 and time2) are not provided!")
            
            if (!is.numeric(t1) || !is.numeric(t2)) stop("t1 and t2 (layers' indicators corresponding to time1 and time2) should be a numeric vector")
            
            p1 <- calc(x1[[t1]], mean,na.rm=TRUE)
            f1 <- calc(x1[[t2]], mean,na.rm=TRUE)
            
            p2 <- calc(x2[[t1]], mean,na.rm=TRUE)
            f2 <- calc(x2[[t2]], mean,na.rm=TRUE)
            
            .velocM(p1,p2,f1,f2,...)
            
          }
)
