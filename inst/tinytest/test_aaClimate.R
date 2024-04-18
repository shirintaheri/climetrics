#-------
filePath <- system.file("external/", package="climetrics") # path to the dataset folder

# read the climate variables using the terra package (you can use the raster package as well):

pr <- rast(paste0(filePath,'/precip.tif'))
tmin <- rast(paste0(filePath,'/tmin.tif'))
tmax <- rast(paste0(filePath,'/tmax.tif'))
tmean <- rast(paste0(filePath,'/tmean.tif'))

n <- readRDS(paste0(filePath,'/dates.rds')) # read corresponding dates
####################

# use rts function in the rts package to make a raster time series:

pr.t <- rts(pr,n) 
tmin.t <- rts(tmin,n)
tmax.t <- rts(tmax,n)
tmean.t <- rts(tmean,n)
#------
###########################
# test of the metric:
#---------
#---------

aa <- aaClimate(precip=pr.t,tmin=tmin.t,tmax=tmax.t,tmean=tmean.t,t1='1991/2000',t2='2010/2020')

expect_equal(round(sum(aa[],na.rm=T)),-400)
