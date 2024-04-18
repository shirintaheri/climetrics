# Authors: Shirin Taheri (taheri.shi@gmail.com); Babak Naimi (naimi.b@gmail.com)
# Date :  Nov. 2020
# Last update :  April 2024
# Version 2.7
# Licence GPL v3
#--------

if (!isGeneric("ccm")) {
  setGeneric("ccm", function(x,...,stat,t1,t2,extreme,longlat,ny,names,verbose=TRUE)
    standardGeneric("ccm"))
}


setMethod('ccm', signature(x='SpatRasterTS'),
          function(x,...,stat,t1,t2,extreme,longlat,ny,names,verbose=TRUE) {
            xx <- list(x=x,...)
            
            if (missing(t1) || missing(t2)) stop("t1 and t2 (Indicators of layers corresponding to time1 and time2) are not provided!")
            
            if (missing(names)) names <- NULL
            
            if (missing(verbose)) verbose <- TRUE
            
            if (missing(longlat)) {
              longlat <- .is.projected(xx[[1]]@raster)
            }
            
            if (missing(stat)) stop('stat should be specified!')
            
            stat <- tolower(stat)
            
            for (i in 1:length(stat)) {
              if (any(c('sed','sla','localanomalies','anomaly') %in% stat[i])) stat[i] <- 'sed'
              else if (any(c('le','lce','lextreme','localextreme','localclimateextreme','lech','lee','extreme') %in% stat[i])) stat[i] <- 'localExtreme'
              else if (any(c('nc','novel','novelc','novelclimate') %in% stat[i])) stat[i] <- 'novelClimate'
              else if (any(c('aac','analogous','aaclimate','analog','aanalog','ac') %in% stat[i])) stat[i] <- 'aaClimate'
              else if (any(c('dac','daclimate','disanalog','danalog') %in% stat[i])) stat[i] <- 'daClimate'
              else if (any(c('ve','velocity','vel') %in% stat[i])) stat[i] <- 'velocity'
              else if (any(c('dve','dvelocity','dvel','distancevelocity') %in% stat[i])) stat[i] <- 'dVelocity'
              else if (any(c('gve','gvelocity','gvel','gradiantvelocity','gradv','gradve') %in% stat[i])) stat[i] <- 'gVelocity'
              else stat[i] <- NA
            }
            
            if (all(is.na(stat))) stop('the metrics specified in stat are unknown!')
            
            if (any(is.na(stat))) {
              warning(paste0('Some of the metrics specified in stat (',length(which(is.na(stat))),' metrics) are unknown that are discarded!'))
              stat <- stat[!is.na(stat)]
            }
            
            o <- vector('list',length = length(stat))
            names(o) <- stat
            
            if ('sed' %in% stat) {
              o[['sed']] <- .sed(xx,t1=t1,t2=t2)
            }
            #--------
            if ('localExtreme' %in% stat) {
              
              if (missing(extreme)) stop('"extreme" is needed...!')
              #-
              if (length(xx) > 2) {
                .xx <- xx[1:2]
                
                if (verbose) cat('\nlocalExtreme metric can be calculated using either one or two climate variables... only the first two variables are used!')
                else warning('localExtreme metric can be calculated using either one or two climate variables... only the first two variables are used!')
                
              } else .xx <- xx
              #-
              if (length(extreme) > length(.xx)) {
                warning(paste0('number of extreme probabilities (',length(extreme),') should be equal to the number of climate parameters (',length(xx),'); Only the first ',length(.xx),' values from extreme is used!'))
                extreme <- extreme[1:length(.xx)]
              }
              #-
              if (length(extreme) < length(.xx)) stop('extreme probabilities should be provided for both climate parameters')
              #-
              if (any(extreme > 1 | extreme < 0)) stop('extreme is a probability value that should be within the range of 0 and 1')  
              #-------------------------
              x1 <- x2 <- list()
              for (i in 1:length(.xx)) {
                x1[[i]] <- .xx[[i]][[t1]]@raster
                x2[[i]] <- .xx[[i]][[t2]]@raster
              }
              
              o[['localExtreme']] <- .eeChange(x1,x2,extreme)
              rm(.xx);gc()
            }
            #-------
            if ('novelClimate' %in% stat) {
              .n <- .nc(xx,t1=t1,t2=t2)
              .q <- global(.n,quantile,probs=0.99,na.rm=TRUE)[1,1]
              .n <- ifel(.n >= .q,.q,.n)
              o[['novelClimate']] <- .n
            }
            #----
            if (any(c('daClimate','aaClimate') %in% stat)) {
              nn <- names(xx)
              nn <- tolower(nn[nn != ""])
              if (length(xx) > 4) stop('The aaClimate metric needs precipitation, minimum and maximum temperature (optional also mean temperature). More than 4 input variables are provided')
              nstat <- c('precip','tmin','tmax','tmean')
              if (is.null(names)) {
                if (length(nn) < (length(xx))) {
                  stop('The "names" argument is missing! names of the input variables should be provided (or the input variables should be provided as the named arguments)!')
                } else {
                  for (i in 2:length(nn)) {
                    nn[i] <- .argMatch(nn[i],nstat)  
                  }
                  if (any(is.na(nn))) stop('The aaClimate metric needs precipitation, minimum and maximum temperature (optional also mean temperature). The input variables are not identified! Specify their names in the names argument!')
                  
                  w <- which(!nstat %in% nn)
                  if (length(w) > 1) {
                    if (!'tmean' %in% nn) {
                      nn[1] <- nstat[which(!nstat[-4] %in% nn)]
                    } else stop('The aaClimate metric needs precipitation, minimum and maximum temperature (optional also mean temperature). The input variables are not identified! Specify their names in the names argument!')
                  } else {
                    nn[1] <- nstat[w]
                  }
                }
                names <- nn
              } else {
                names <- tolower(names)
                if (length(names) != length(xx)) stop('The length of the provided "names" is not equal to the number of the input variables!')
              }
              nstat <- c('precip','tmin','tmax','tmean')
              for (i in 2:length(names)) {
                names[i] <- .argMatch(names[i],nstat)  
              }
              
              #------
              if (!'tmean' %in% names) {
                names(xx) <- names
                .tmean <- (xx$tmin@raster + xx$tmax@raster) / 2
                .tmean <- rts(.tmean,index(xx$tmin))
                xx$tmean <- .tmean
                names <- c(names,'tmean')
              }
              #----------
              
              names(xx) <- names
              
              tmin1 <- xx$tmin[[t1]]
              tmin2 <- xx$tmin[[t2]]
              tmax1 <- xx$tmax[[t1]]
              tmax2 <- xx$tmax[[t2]]
              tmean1 <- xx$tmean[[t1]]
              tmean2 <- xx$tmean[[t2]]
              prt1 <- xx$precip[[t1]]
              prt2 <- xx$precip[[t2]]
              
              k1 <- kgc(apply.months(prt1),apply.months(tmin1),apply.months(tmax1),apply.months(tmean1))
              k2 <- kgc(apply.months(prt2),apply.months(tmin2),apply.months(tmax2),apply.months(tmean2))
              
              if ('daClimate' %in% stat) o[['daCliamte']] <-  .disAnalogous(k1,k2)
              if ('aaClimate' %in% stat) o[['aaCliamte']] <-  .analogousClimate(k1,k2)
              
            }
            
            if ('velocity' %in% stat) {
              
              if (length(xx) > 1) {
                if (length(xx) > 2) {
                  if (verbose) cat('\nTo claculate the stat "velocity" (ve), the first two climate variables are used as the metric is implemented based on the approach in Hamnan et al. (2014), works based on two climate variable (For multiple variables, you may either use a PCA transformation and take the first two components)!')
                  else warning('To claculate the stat "velocity" (ve), the first two climate variables are used as the metric is implemented based on the approach in Hamnan et al. (2014), works based on two climate variable (For multiple variables, you may either use a PCA transformation and take the first two components)!')
                }
              } else {
                if (length(stat) > 1) {
                  if (verbose) cat('\nThe stat "velocity" (ve) is ignored as it is calculated based on two climate variables (Hamnan et al., 2014)')
                  else warning('The stat "velocity" (ve) is ignored as it is calculated based on two climate variables (Hamnan et al., 2014)')
                } else stop('This method of velocity is implemented based on the approach in Hamnan et al. (2014) that works based on two climate variable; For single variable, you may use dVe (dVelocity) stat!')
              }
              
              
              #----
              p1 <- app(xx[[1]][[t1]]@raster, 'mean',na.rm=TRUE)
              f1 <- app(xx[[1]][[t2]]@raster, 'mean',na.rm=TRUE)
              
              p2 <- app(xx[[2]][[t1]]@raster, 'mean',na.rm=TRUE)
              f2 <- app(xx[[2]][[t2]]@raster, 'mean',na.rm=TRUE)
              
              o[['velocity']] <-  .velocMTerra(p1,p2,f1,f2)
              
            } 
            
            if ('dVelocity' %in% stat) {
              if (missing(ny)) ny <- .getNyears(xx[[1]]@time,t1,t2)
              
              xt1 <- xt2 <- list()
              
              for (i in 1:length(xx)) {
                xt1[[i]] <- mean(xx[[i]][[t1]]@raster,na.rm=TRUE)
                xt2[[i]] <- mean(xx[[i]][[t2]]@raster,na.rm=TRUE)
              }
              #---------------
              
              r <- .getScaledMultiVariateIntoOne(xt1,xt2)
              s <- .spatialgrTerra(r$t1)
              t <- .tempgr(r$t1,r$t2,ny = ny)
              o[['dVelocity']] <- .getVelocity(s,t)
            }
            
            if ('gVelocity' %in% stat) {
              if (length(xx) > 1) {
                
                if (verbose) cat('\ngVelocity is calculated for a single climate variable. Since multiple variables are provided, average velocity is returned! ')
                else warning('gVelocity is calculated for a single climate variable. Since multiple variables are provided, average velocity is returned! ')
                
              }
              
              r <- list()
              
              for (i in 1:length(xx)) {
                r[[i]] <- gVelocity(xx[[i]])
              }
              #---------------
              if (length(xx) > 1) {
                o[['gVelocity']] <- mean(rast(r),na.rm=TRUE)
              } else {
                o[['gVelocity']] <- r[[1]]
              }
            }
              
             if (length(o) > 1) {
               o <- rast(o)
               names(o) <- stat
             } else {
               o <- o[[1]]
             }
            o
          }
)



