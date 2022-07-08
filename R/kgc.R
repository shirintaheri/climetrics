# Authors: Shirin Taheri (taheri.shi@gmail.com); Babak Naimi (Naimi.b@gmail.com)
# Date :  March 2021
# Last update :  July 2022
# Version 1.2
# Licence GPL v3
#--------

#-------------- Koppen Geiger climate classification




# # whether 70% of precipitation fall in summer:
# .p70 <- function(pm) {
#   if (inherits(pm,'Raster')) {
#     if (nlayers(pm) != 12) stop('monthly precipitation does not have 12 layers!')
#     calc(pm,function(x) {
#       .s <- sum(x,na.rm = TRUE) # total precipitation
#       (sum(x[4:9],na.rm = TRUE) / .s) >= 0.7 # whether 70% of precipitation is in summer
#     })
#   } else if (inherits(pm,'SpatRaster')) {
#     if (nlyr(pm) != 12) stop('monthly precipitation oes not have 12 layers!')
#     app(pm,function(x) {
#       .s <- sum(x,na.rm = TRUE) # total precipitation
#       (sum(x[4:9],na.rm = TRUE) / .s) >= 0.7 # whether 70% of precipitation is in summer
#     })
#   } else stop('The input object is not a Raster!')
# }
#---------
# if 70% of precipitation fall in summer, return 1;
# if 70% of precipitation fall in winter, return 2,
# otherwise, returns 3:
.p70 <- function(pm) {
  if (inherits(pm,'Raster')) {
    if (nlayers(pm) != 12) stop('monthly precipitation should have 12 layers!')
    r <- raster(pm)
  } else if (inherits(pm,'SpatRaster')) {
    if (nlyr(pm) != 12) stop('monthly precipitation should have 12 layers!')
    r <- rast(pm[[1]])
  } else stop('The input object is not a Raster!')
  #------
  
  pdf <- as.data.frame(pm,xy=TRUE,na.rm=TRUE)
  xy <- pdf[,1:2]
  pdf <- pdf[,-c(1:2)]
  
  if (.getProj(pm) == 'projected') xy <- .getLongLat(xy,projection(pm))
  
  if (!all(xy$y >= 0)) {
    if (!all(xy$y < 0)) {
      w <- which(xy$y >= 0)
      v <- rep(NA,nrow(pdf))
      
      
      v[w] <- apply(pdf[w,],1,FUN = function(x,...) {
        .s <- sum(x,na.rm = TRUE) # total precipitation
        
        if (.s == 0) return(3)
        
        s1 <- (sum(x[4:9],na.rm = TRUE) / .s) >= 0.7 # whether 70% of precipitation is in summer
        if (s1) return(1)
        else {
          s1 <- (sum(x[-c(4:9)],na.rm = TRUE) / .s) >= 0.7 # whether 70% of precipitation is in summer
          if (s1) return(2)
          else 3
        }
      })
      
      v[-w] <- apply(pdf[-w,],1,FUN = function(x) {
        .s <- sum(x,na.rm = TRUE) # total precipitation
        
        if (.s == 0) return(3)
        
        s1 <- (sum(x[-c(4:9)],na.rm = TRUE) / .s) >= 0.7 # whether 70% of precipitation is in summer
        if (s1) return(1)
        else {
          s1 <- (sum(x[c(4:9)],na.rm = TRUE) / .s) >= 0.7 # whether 70% of precipitation is in summer
          if (s1) return(2)
          else 3
        }
      })
    } else {
      v <- apply(pdf,1,FUN = function(x) {
        .s <- sum(x,na.rm = TRUE) # total precipitation
        
        if (.s == 0) return(3)
        
        s1 <- (sum(x[-c(4:9)],na.rm = TRUE) / .s) >= 0.7 # whether 70% of precipitation is in summer
        if (s1) return(1)
        else {
          s1 <- (sum(x[c(4:9)],na.rm = TRUE) / .s) >= 0.7 # whether 70% of precipitation is in summer
          if (s1) return(2)
          else 3
        }
      })
    }
  } else {
    v <- apply(pdf,1,FUN = function(x) {
      .s <- sum(x,na.rm = TRUE) # total precipitation
      
      if (.s == 0) return(3)
      
      s1 <- (sum(x[4:9],na.rm = TRUE) / .s) >= 0.7 # whether 70% of precipitation is in summer
      if (s1) return(1)
      else {
        s1 <- (sum(x[-c(4:9)],na.rm = TRUE) / .s) >= 0.7 # whether 70% of precipitation is in summer
        if (s1) return(2)
        else 3
      }
    })
  }
  #----
  
  r[cellFromXY(r,as.matrix(xy))] <- v
  r
}
#---------
####################################

# tmean: monthly mean temperature
.MAT <- function(tmean) {
  if (inherits(tmean,'Raster')) {
    calc(tmean,mean,na.rm=TRUE)
  } else if (inherits(tmean,'SpatRaster')) {
    app(tmean,mean,na.rm=TRUE)
  } else stop('The input object is not a Raster!')
  
}
#-----
# pm: monthly precipitation
.MAP <- function(pm) {
  if (inherits(pm,'Raster')) {
    calc(pm,mean,na.rm=TRUE)
  } else if (inherits(pm,'SpatRaster')) {
    app(pm,mean,na.rm=TRUE)
  } else stop('The input object is not a Raster!')
}
#------


# Pthreshold: check .p70, if true, it would be 2*MAT+28, otherwise 2*MAT+14
# pm: monthly precipitation
# mat: mean annual temperature
.Pth <- function(pm,mat) {
  
  p7 <- .p70(pm)
  
  if (inherits(pm,'Raster')) {
    r <- raster(p7)
    
  } else if (inherits(pm,'SpatRaster')) {
    r <- rast(p7)
  } else stop('The input object is not a Raster!')
  
  w <- which(p7[] == 1)
  if (length(w) > 0) r[w] <- (2*mat[w]) + 28
  
  w <- which(p7[] == 2)
  if (length(w) > 0) r[w] <- (2*mat[w])
  
  w <- which(p7[] == 3)
  if (length(w) > 0) r[w] <- (2*mat[w]) + 14
  
  r
}

#---------

#precipitation in driest month
.Pdry <- function(pm) {
  
  if (inherits(pm,'Raster')) {
    
    if (nlayers(pm) != 12) stop('monthly precipitation should have 12 layers!')
    calc(pm,function(x,...) {
      x <- x[!is.na(x)]
      if (length(x) > 0) x[which.min(x)]
      else NA
    })
  } else if (inherits(pm,'SpatRaster')) {
    
    if (nlyr(pm) != 12) stop('monthly precipitation should have 12 layers!')
    app(pm,function(x,...) {
      x <- x[!is.na(x)]
      if (length(x) > 0) x[which.min(x)]
      else NA
    })
  } else stop('The input object is not a Raster!')
  
}
#####
# precipitation in driest month in summer
.Psdry <- function(pm) {
  if (inherits(pm,'Raster')) {
    if (nlayers(pm) != 12) stop('monthly precipitation should have 12 layers!')
    r <- raster(pm)
  } else if (inherits(pm,'SpatRaster')) {
    if (nlyr(pm) != 12) stop('monthly precipitation should have 12 layers!')
    r <- rast(pm[[1]])
  } else stop('The input object is not a Raster!')
  
  #--------------
  pdf <- as.data.frame(pm,xy=TRUE,na.rm=TRUE)
  xy <- pdf[,1:2]
  pdf <- pdf[,-c(1:2)]
  
  if (.getProj(pm) == 'projected') xy <- .getLongLat(xy,projection(pm))
  
  if (!all(xy$y >= 0)) {
    if (!all(xy$y < 0)) {
      w <- which(xy$y >= 0)
      v <- rep(NA,nrow(pdf))
      v[w] <- apply(pdf[w,4:9],1,FUN = function(x) {
        x <- x[!is.na(x)]
        if (length(x) > 0) x[which.min(x)]
        else NA
      })
      
      v[-w] <- apply(pdf[-w,-c(4:9)],1,FUN = function(x) {
        x <- x[!is.na(x)]
        if (length(x) > 0) x[which.min(x)]
        else NA
      })
    } else {
      v <- apply(pdf[,-c(4:9)],1,FUN = function(x) {
        x <- x[!is.na(x)]
        if (length(x) > 0) x[which.min(x)]
        else NA
      })
    }
  } else {
    v <- apply(pdf[,c(4:9)],1,FUN = function(x) {
      x <- x[!is.na(x)]
      if (length(x) > 0) x[which.min(x)]
      else NA
    })
  }
  #----
  
  r[cellFromXY(r,as.matrix(xy))] <- v
  r
  
}
####

# precipitation in driest month in winter
.Pwdry <- function(pm) {
  if (inherits(pm,'Raster')) {
    if (nlayers(pm) != 12) stop('monthly precipitation should have 12 layers!')
    r <- raster(pm)
  } else if (inherits(pm,'SpatRaster')) {
    if (nlyr(pm) != 12) stop('monthly precipitation should have 12 layers!')
    r <- rast(pm[[1]])
  } else stop('The input object is not a Raster!')
  #-----------------
  
  pdf <- as.data.frame(pm,xy=TRUE,na.rm=TRUE)
  xy <- pdf[,1:2]
  pdf <- pdf[,-c(1:2)]
  
  if (.getProj(pm) == 'projected') xy <- .getLongLat(xy,projection(pm))
  
  if (!all(xy$y >= 0)) {
    if (!all(xy$y < 0)) {
      w <- which(xy$y >= 0)
      v <- rep(NA,nrow(pdf))
      v[w] <- apply(pdf[w,-c(4:9)],1,FUN = function(x) {
        x <- x[!is.na(x)]
        if (length(x) > 0) x[which.min(x)]
        else NA
      })
      
      v[-w] <- apply(pdf[-w,c(4:9)],1,FUN = function(x) {
        x <- x[!is.na(x)]
        if (length(x) > 0) x[which.min(x)]
        else NA
      })
    } else {
      v <- apply(pdf[,c(4:9)],1,FUN = function(x) {
        x <- x[!is.na(x)]
        if (length(x) > 0) x[which.min(x)]
        else NA
      })
    }
  } else {
    v <- apply(pdf[,-c(4:9)],1,FUN = function(x) {
      x <- x[!is.na(x)]
      if (length(x) > 0) x[which.min(x)]
      else NA
    })
  }
  #----
  
  r[cellFromXY(r,as.matrix(xy))] <- v
  r
}
#------
# precipitation in wettest month in summer
.Pswet <- function(pm) {
  if (inherits(pm,'Raster')) {
    if (nlayers(pm) != 12) stop('monthly precipitation should have 12 layers!')
    r <- raster(pm)
  } else if (inherits(pm,'SpatRaster')) {
    if (nlyr(pm) != 12) stop('monthly precipitation should have 12 layers!')
    r <- rast(pm[[1]])
  } else stop('The input object is not a Raster!')
  #-----------------
  
  pdf <- as.data.frame(pm,xy=TRUE,na.rm=TRUE)
  xy <- pdf[,1:2]
  pdf <- pdf[,-c(1:2)]
  
  if (.getProj(pm) == 'projected') xy <- .getLongLat(xy,projection(pm))
  
  if (!all(xy$y >= 0)) {
    if (!all(xy$y < 0)) {
      w <- which(xy$y >= 0)
      v <- rep(NA,nrow(pdf))
      v[w] <- apply(pdf[w,4:9],1,FUN = function(x) {
        x <- x[!is.na(x)]
        if (length(x) > 0) x[which.max(x)]
        else NA
      })
      
      v[-w] <- apply(pdf[-w,-c(4:9)],1,FUN = function(x) {
        x <- x[!is.na(x)]
        if (length(x) > 0) x[which.max(x)]
        else NA
      })
    } else {
      v <- apply(pdf[,-c(4:9)],1,FUN = function(x) {
        x <- x[!is.na(x)]
        if (length(x) > 0) x[which.max(x)]
        else NA
      })
    }
  } else {
    v <- apply(pdf[,c(4:9)],1,FUN = function(x) {
      x <- x[!is.na(x)]
      if (length(x) > 0) x[which.max(x)]
      else NA
    })
  }
  #----
  
  r[cellFromXY(r,as.matrix(xy))] <- v
  r
  
}
####
# precipitation in wettest month in winter
.Pwwet <- function(pm) {
  if (inherits(pm,'Raster')) {
    if (nlayers(pm) != 12) stop('monthly precipitation should have 12 layers!')
    r <- raster(pm)
  } else if (inherits(pm,'SpatRaster')) {
    if (nlyr(pm) != 12) stop('monthly precipitation should have 12 layers!')
    r <- rast(pm[[1]])
  } else stop('The input object is not a Raster!')
  #-----------------
  
  pdf <- as.data.frame(pm,xy=TRUE,na.rm=TRUE)
  xy <- pdf[,1:2]
  pdf <- pdf[,-c(1:2)]
  
  if (.getProj(pm) == 'projected') xy <- .getLongLat(xy,projection(pm))
  
  if (!all(xy$y >= 0)) {
    if (!all(xy$y < 0)) {
      w <- which(xy$y >= 0)
      v <- rep(NA,nrow(pdf))
      v[w] <- apply(pdf[w,-c(4:9)],1,FUN = function(x) {
        x <- x[!is.na(x)]
        if (length(x) > 0) x[which.max(x)]
        else NA
      })
      
      v[-w] <- apply(pdf[-w,c(4:9)],1,FUN = function(x) {
        x <- x[!is.na(x)]
        if (length(x) > 0) x[which.max(x)]
        else NA
      })
    } else {
      v <- apply(pdf[,c(4:9)],1,FUN = function(x) {
        x <- x[!is.na(x)]
        if (length(x) > 0) x[which.max(x)]
        else NA
      })
    }
  } else {
    v <- apply(pdf[,-c(4:9)],1,FUN = function(x) {
      x <- x[!is.na(x)]
      if (length(x) > 0) x[which.max(x)]
      else NA
    })
  }
  #----
  
  r[cellFromXY(r,as.matrix(xy))] <- v
  r
}
#temperature in coldest month
.Tcold <- function(tmin) {
  if (inherits(tmin,'Raster')) {
    if (nlayers(tmin) != 12) stop('monthly temperature should have 12 layers!')
    calc(tmin,function(x,...) {
      x <- x[!is.na(x)]
      if (length(x) > 0) x[which.min(x)]
      else NA
    })
  } else if (inherits(tmin,'SpatRaster')) {
    if (nlyr(tmin) != 12) stop('monthly temperature should have 12 layers!')
    app(tmin,function(x,...) {
      x <- x[!is.na(x)]
      if (length(x) > 0) x[which.min(x)]
      else NA
    })
  } else stop('The input object is not a Raster!')
}
#####
#temperature in coldest month
.Thot <- function(tmax) {
  if (inherits(tmax,'Raster')) {
    if (nlayers(tmax) != 12) stop('monthly temperature should have 12 layers!')
    calc(tmax,function(x,...) {
      x <- x[!is.na(x)]
      if (length(x) > 0) x[which.max(x)]
      else NA
    })
  } else if (inherits(tmax,'SpatRaster')) {
    if (nlyr(tmax) != 12) stop('monthly temperature should have 12 layers!')
    app(tmax,function(x,...) {
      x <- x[!is.na(x)]
      if (length(x) > 0) x[which.max(x)]
      else NA
    })
  } else stop('The input object is not a Raster!')
}
#####
#number of month with air temperature > 10
.Tm10 <- function(tmean) {
  if (inherits(tmean,'Raster')) {
    if (nlayers(tmean) != 12) stop('monthly temperature should have 12 layers!')
    calc(tmean,function(x,...) {
      x <- x[!is.na(x)]
      if (length(x) > 0) length(which(x > 10))
      else NA
    })
  } else if (inherits(tmean,'SpatRaster')) {
    if (nlyr(tmean) != 12) stop('monthly temperature should have 12 layers!')
    app(tmean,function(x,...) {
      x <- x[!is.na(x)]
      if (length(x) > 0) length(which(x > 10))
      else NA
    })
  } else stop('The input object is not a Raster!')
}
#####

.KGC_classes <- list(Af=c(1,1),Am=c(2,1),Aw=c(3,1),BWh=c(4,2),BWk=c(5,2),BSh=c(6,2),BSk=c(7,2),
                    Csa=c(8,3),Csb=c(9,3),Csc=c(10,3),Cwa=c(11,3),Cwb=c(12,3),Cwc=c(13,3),Cfa=c(14,3),
                    Cfb=c(15,3),Cfc=c(16,3),Dsa=c(17,4),Dsb=c(18,4),Dsc=c(19,4),Dsd=c(20,4),Dwa=c(21,4),
                    Dwb=c(22,4),Dwc=c(23,4),Dwd=c(24,4),Dfa=c(25,4),Dfb=c(26,4),Dfc=c(27,4),Dfd=c(28,4),
                    ET=c(29,5),EF=c(30,5))
#------------
#----------------
# .kgcR <- function(pm,tmean,tmin,tmax) {
#   k <- data.frame(t(as.data.frame(.KGC_classes)),names(.KGC_classes))
#   names(k) <- c('code','class','names')
#   #---------
#   .mat <- .MAT(tmean)
#   .pth <- .Pth(pm, .mat)
#   .map <- .MAP(pm)
#   .tc <- .Tcold(tmin)
#   .th <- .Thot(tmax)
#   .pdry <- .Pdry(pm)
#   .psdry <- .Psdry(pm)
#   .pww <- .Pwwet(pm)
#   .tm10 <- .Tm10(tmean)
#   .pwd <- .Pwdry(pm)
#   .psw <- .Pswet(pm)
#   #------
#   # empty raster
#   if (inherits(pm,'Raster')) r <- raster(.mat)
#   else r <- rast(.mat)
#   
#   #-----
#   .c <- which(!is.na(.map[])) # non-NA cells
#   #------------
#   
#   #### BWh, BWk (4:7):
#   w <- .map[.c] < (10 * .pth[.c])
#   w2 <- .map[.c] < (5 * .pth[.c])
#   w3 <- .mat[.c] >= 18
#   
#   .w <- .c[w & w2 & w3]
#   if (length(.w) > 0) r[.w] <- 4
#   .w <- .c[w & w2 & !w3]
#   if (length(.w) > 0) r[.w] <- 5
#   
#   .w <- .c[w & !w2 & w3]
#   if (length(.w) > 0) r[.w] <- 6
#   .w <- .c[w & !w2 & !w3]
#   if (length(.w) > 0) r[.w] <- 7
#   
#   #=======================
#   #### Af, Am, Aw (1,2,3):
#   w <- (!r[.c] %in% c(4:7)) & (.tc[.c] >= 18) # Not(B) & Tcold >= 18
#   w2 <- .pdry[.c] >= 60
#   .w <- .c[w & w2]
#   if (length(.w) > 0) r[.w] <- 1
#   
#   w2 <- .pdry[.c] >= (100 - (.map[.c] / 25))
#   w3 <- r[.c] != 1
#   .w <- .c[w & w2 & w3]
#   if (length(.w) > 0) r[.w] <- 2
#   
#   .w <- .c[w & !w2 & w3]
#   if (length(.w) > 0) r[.w] <- 3
#   #========================
#   #### Csa.... (8:16):
#   w <- (!r[.c] %in% c(4:7)) & (.th[.c] > 10) & (.tc[.c] > 0) & (.tc[.c] < 18) # Not(B) & Thot>10 & 0<Tcold <10
#   #Csa,Csb,Csc:
#   w1 <- (.psdry[.c] < 40) & (.psdry[.c] < (.pww[.c]/3))
#   w2 <- .th[.c] >= 22
#   .w <- .c[w & w1 & w2]
#   if (length(.w) > 0) r[.w] <- 8 # CSa
#   
#   w3 <- .tm10[.c] >= 4
#   .w <- .c[w & w1 & !w2 & w3]
#   if (length(.w) > 0) r[.w] <- 9 # CSb
#   w4 <- .tm10[.c] >= 1
#   .w <- .c[w & w1 & !w2 & !w3 & w4]
#   if (length(.w) > 0) r[.w] <- 10 # CSc
#   #--
#   #Cwa,Cwb,C.w:
#   w.1 <- (.pwd[.c] < (.psw[.c]/10))
#   .w <- .c[w & w.1 & w2]
#   if (length(.w) > 0) r[.w] <- 11 # Cwa
#   
#   .w <- .c[w & w.1 & !w2 & w3]
#   if (length(.w) > 0) r[.w] <- 12 # Cwb
#   w4 <- .tm10[.c] >= 1
#   .w <- .c[w & w.1 & !w2 & !w3 & w4]
#   if (length(.w) > 0) r[.w] <- 13 # C.w
#   #--
#   #Cfa,Cfb,Cfc:
#   w1 <- !(w1 | w.1)
#   .w <- .c[w & w1 & w2]
#   if (length(.w) > 0) r[.w] <- 14 # Cfa
#   
#   .w <- .c[w & w1 & !w2 & w3]
#   if (length(.w) > 0) r[.w] <- 15 # Cfb
#   w4 <- .tm10[.c] >= 1
#   .w <- .c[w & w1 & !w2 & !w3 & w4]
#   if (length(.w) > 0) r[.w] <- 16 # Cfc
#   #=============
#   #### D (17:28):
#   w <- (!r[.c] %in% c(4:7)) & (.th[.c] > 10) & (.tc[.c] <= 0) # Not(B) & Thot>10 & Tcold <=0
#   #Dsa,Dsb,Dsc,Dsd:
#   w1 <- (.psdry[.c] < 40) & (.psdry[.c] < (.pww[.c]/3))
#   wa <- .th[.c] >= 22
#   .w <- .c[w & w1 & wa]
#   if (length(.w) > 0) r[.w] <- 17 # Dsa
#   
#   wb <- !wa & .tm10[.c] >= 4 # not(a) & Tmon10 >=4
#   .w <- .c[w & w1 & wb]
#   if (length(.w) > 0) r[.w] <- 18 # Dsb
#   
#   wd <- !(wa | wb) & (.tc[.c] < -18)
#   .w <- .c[w & w1 & wd]
#   if (length(.w) > 0) r[.w] <- 20 # Dsd
#   
#   wc <- !(wa | wb | wd)
#   .w <- .c[w & w1 & wc]
#   if (length(.w) > 0) r[.w] <- 19 # Dsc
#   
#   #--
#   #Dwa,Dwb,Dwc,Dwd:
#   w.1 <- (.pwd[.c] < (.psw[.c]/10))
#   .w <- .c[w & w.1 & wa]
#   if (length(.w) > 0) r[.w] <- 21 # Dwa
#   
#   .w <- .c[w & w.1 & wb]
#   if (length(.w) > 0) r[.w] <- 22 # Dsb
#   
#   .w <- .c[w & w.1 & wd]
#   if (length(.w) > 0) r[.w] <- 24 # Dsd
#   
#   .w <- .c[w & w.1 & wc]
#   if (length(.w) > 0) r[.w] <- 23 # Dsc
#   
#   #--
#   #Dfa,Dfb,Dfc,Dfd:
#   w1 <- !(w1 | w.1)
#   .w <- .c[w & w1 & wa]
#   if (length(.w) > 0) r[.w] <- 25 # Dfa
#   
#   .w <- .c[w & w1 & wb]
#   if (length(.w) > 0) r[.w] <- 26 # Dfb
#   
#   .w <- .c[w & w1 & wd]
#   if (length(.w) > 0) r[.w] <- 28 # Dfd
#   
#   .w <- .c[w & w1 & wc]
#   if (length(.w) > 0) r[.w] <- 27 # Dfc
#   ################
#   # ET, EF (29:30):
#   w <- (!r[.c] %in% c(4:7)) & (.th[.c] <= 18) # Not(B) & Thot <= 10
#   wT <- .th[.c] > 0
#   .w <- .c[w & wT]
#   if (length(.w) > 0) r[.w] <- 29 #ET
#   
#   wF <- .th[.c] <= 0
#   .w <- .c[w & wF]
#   if (length(.w) > 0) r[.w] <- 30 #EF
#   #========================
#   r
# }

.kgcR <- function(pm,tmean,tmin,tmax) {
  k <- data.frame(t(as.data.frame(.KGC_classes)),names(.KGC_classes))
  names(k) <- c('code','class','names')
  #---------
  .mat <- .MAT(tmean)
  .pth <- .Pth(pm, .mat)
  .map <- .MAP(pm)
  .tc <- .Tcold(tmin)
  .th <- .Thot(tmax)
  .pdry <- .Pdry(pm)
  .psdry <- .Psdry(pm)
  .pww <- .Pwwet(pm)
  .tm10 <- .Tm10(tmean)
  .pwd <- .Pwdry(pm)
  .psw <- .Pswet(pm)
  #------
  # empty raster
  if (inherits(pm,'Raster')) r <- raster(.mat)
  else r <- rast(.mat)
  
  #-----
  .c <- which(!is.na(.map[])) # non-NA cells
  #------------
  
  #=======================
  #### Af, Am, Aw (1,2,3):
  w1 <- .tc[.c] >= 18 # Tcold >= 18
  w2 <- .pdry[.c] >= 60
  .wn <- .c[w1 & w2]
  if (length(.wn) > 0) r[.wn] <- 1 # Af
  
  w3 <- .pdry[.c] >= (100 - (.map[.c] / 25))
  
  .wn <- .c[w1 & !c(w1 & w2) & w3]
  if (length(.wn) > 0) r[.wn] <- 2 # Am
  
  w3 <- .pdry[.c] < (100 - (.map[.c] / 25))
  .wn <- .c[w1 & !c(w1 & w2) & w3]
  if (length(.wn) > 0) r[.wn] <- 3 # Aw
  #========================
  
  w <- is.na(r[.c])
  #### BWh, BWk (4:7):
  w1 <- .map[.c] < (10 * .pth[.c])
  w2 <- .map[.c] < (5 * .pth[.c])
  w3 <- .map[.c] >= (5 * .pth[.c])
  w4 <- .mat[.c] >= 18
  w5 <- .mat[.c] < 18
  
  .wn <- .c[w & w1 & w2 & w4]
  if (length(.wn) > 0) r[.wn] <- 4 #BW h
  .wn <- .c[w & w1 & w2 & w5]
  if (length(.wn) > 0) r[.wn] <- 5 #BW k
  
  .wn <- .c[w & w1 & w3 & w4]
  if (length(.wn) > 0) r[.wn] <- 6
  .wn <- .c[w & w1 & w3 & w5]
  if (length(.wn) > 0) r[.wn] <- 7
  
  #### Csa.... (8:16):
  # w <- c(!r[.c] %in% c(1:7))
  # w1 <- (.th[.c] > 10) & (.tc[.c] > 0) & (.tc[.c] < 18) # C
  # w2 <- (.psdry[.c] < 40) & (.psdry[.c] < (.pww[.c] / 3)) #..s
  # w3 <- (.pwdry[.c] < (.psw[.c] / 10)) # ..w
  # 
  # w4 <- .th[.c] >= 22 # a
  # .wn <- .c[w & w1 & w2 & w4]
  # if (length(.wn) > 0) r[.wn] <- 8 # CSa
  # 
  # w5 <- .tm10[.c] >= 4 # b
  # 
  # .wn <- .c[w & w1 & w2 & !w4 & w5]
  # if (length(.wn) > 0) r[.wn] <- 9 # CSb
  # 
  # w6 <- (.tm10[.c] >= 1) & (.tm10[.c] < 4) #... c                                                                                                         c
  # 
  # .wn <- .c[w & w1 & w2 & !w4 & !w5 & w6]
  # 
  # if (length(.wn) > 0) r[.wn] <- 10 # CSc
  # #----
  # 
  # .wn <- .c[w & w1 & w3 & w4]
  # if (length(.wn) > 0) r[.wn] <- 11 # CWa
  # 
  # .wn <- .c[w & w1 & w3 & !w4 & w5]
  # if (length(.wn) > 0) r[.wn] <- 12 # CWb
  # 
  # .wn <- .c[w & w1 & w3 & !w4 & !w5 & w6]
  # if (length(.wn) > 0) r[.wn] <- 13 # CWc
  # w23 <- !c((w1 & w2) | (w1 & w3)) # ..f
  
  #.wn <- .c[w &  !w4 & !w5 & w6]
  
  w <- (!r[.c] %in% c(1:7)) & (.th[.c] > 10) & (.tc[.c] > 0) & (.tc[.c] < 18) # Not(B) & Thot>10 & 0<Tcold <10
  #Csa,Csb,Csc:
  w1 <- (.psdry[.c] < 40) & (.psdry[.c] < (.pww[.c]/3))
  w2 <- .th[.c] >= 22
  .w <- .c[w & w1 & w2]
  if (length(.w) > 0) r[.w] <- 8 # CSa
  
  w3 <- .tm10[.c] >= 4
  .w <- .c[w & w1 & !w2 & w3]
  if (length(.w) > 0) r[.w] <- 9 # CSb
  w4 <- .tm10[.c] >= 1
  .w <- .c[w & w1 & !w2 & !w3 & w4]
  if (length(.w) > 0) r[.w] <- 10 # CSc
  #--
  #Cwa,Cwb,C.w:
  w.1 <- (.pwd[.c] < (.psw[.c]/10))
  .w <- .c[w & w.1 & w2]
  if (length(.w) > 0) r[.w] <- 11 # Cwa
  
  .w <- .c[w & w.1 & !w2 & w3]
  if (length(.w) > 0) r[.w] <- 12 # Cwb
  w4 <- .tm10[.c] >= 1
  .w <- .c[w & w.1 & !w2 & !w3 & w4]
  if (length(.w) > 0) r[.w] <- 13 # C.w
  #--
  #Cfa,Cfb,Cfc:
  w1 <- !(w1 | w.1)
  .w <- .c[w & w1 & w2]
  if (length(.w) > 0) r[.w] <- 14 # Cfa
  
  .w <- .c[w & w1 & !w2 & w3]
  if (length(.w) > 0) r[.w] <- 15 # Cfb
  w4 <- .tm10[.c] >= 1
  .w <- .c[w & w1 & !w2 & !w3 & w4]
  if (length(.w) > 0) r[.w] <- 16 # Cfc
  #=============
  #### D (17:28):
  
  w <- (!r[.c] %in% c(1:16)) & (.th[.c] > 10) & (.tc[.c] <= 0) # Not(B) & Thot>10 & Tcold <=0
  #Dsa,Dsb,Dsc,Dsd:
  w1 <- (.psdry[.c] < 40) & (.psdry[.c] < (.pww[.c]/3))
  wa <- .th[.c] >= 22
  .w <- .c[w & w1 & wa]
  if (length(.w) > 0) r[.w] <- 17 # Dsa
  
  wb <- !wa & .tm10[.c] >= 4 # not(a) & Tmon10 >=4
  .w <- .c[w & w1 & wb]
  if (length(.w) > 0) r[.w] <- 18 # Dsb
  
  wd <- !(wa | wb) & (.tc[.c] < -18)
  .w <- .c[w & w1 & wd]
  if (length(.w) > 0) r[.w] <- 20 # Dsd
  
  wc <- !(wa | wb | wd)
  .w <- .c[w & w1 & wc]
  if (length(.w) > 0) r[.w] <- 19 # Dsc
  
  #--
  #Dwa,Dwb,Dwc,Dwd:
  w.1 <- (.pwd[.c] < (.psw[.c]/10))
  .w <- .c[w & w.1 & wa]
  if (length(.w) > 0) r[.w] <- 21 # Dwa
  
  .w <- .c[w & w.1 & wb]
  if (length(.w) > 0) r[.w] <- 22 # Dsb
  
  .w <- .c[w & w.1 & wd]
  if (length(.w) > 0) r[.w] <- 24 # Dsd
  
  .w <- .c[w & w.1 & wc]
  if (length(.w) > 0) r[.w] <- 23 # Dsc
  
  #--
  #Dfa,Dfb,Dfc,Dfd:
  w1 <- !(w1 | w.1)
  .w <- .c[w & w1 & wa]
  if (length(.w) > 0) r[.w] <- 25 # Dfa
  
  .w <- .c[w & w1 & wb]
  if (length(.w) > 0) r[.w] <- 26 # Dfb
  
  .w <- .c[w & w1 & wd]
  if (length(.w) > 0) r[.w] <- 28 # Dfd
  
  .w <- .c[w & w1 & wc]
  if (length(.w) > 0) r[.w] <- 27 # Dfc
  ################
  # ET, EF (29:30):
  w <- (!r[.c] %in% c(1:27)) & (.th[.c] <= 18) # Not(B) & Thot <= 10
  wT <- .th[.c] > 0
  .w <- .c[w & wT]
  if (length(.w) > 0) r[.w] <- 29 #ET
  
  wF <- .th[.c] <= 0
  .w <- .c[w & wF]
  if (length(.w) > 0) r[.w] <- 30 #EF
  #========================
  r
}

#######################


if (!isGeneric("kgc")) {
  setGeneric("kgc", function(p,tmin,tmax,tmean)
    standardGeneric("kgc"))
}
#------------

setMethod('kgc', signature(p='RasterStackBrick','RasterStackBrick','RasterStackBrick'),
          function(p,tmin,tmax,tmean) {
            if (missing(tmean)) {
              tmean <- (tmin + tmax) / 2
            }
            
            .kgcR(p,tmean,tmin,tmax)
          }
)
#---------
setMethod('kgc', signature(p='SpatRaster','SpatRaster','SpatRaster'),
          function(p,tmin,tmax,tmean) {
            if (missing(tmean)) {
              tmean <- (tmin + tmax) / 2
            }
            
            .kgcR(p,tmean,tmin,tmax)
          }
)



