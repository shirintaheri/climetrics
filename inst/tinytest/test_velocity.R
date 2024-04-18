filePath <- system.file("external/", package="climetrics") # path to the dataset folder

pr <- rast(paste0(filePath,'/precip.tif'))
tmean <- rast(paste0(filePath,'/tmean.tif'))

n <- readRDS(paste0(filePath,'/dates.rds')) # corresoinding dates
####################

# use rts function in the rts package to make a raster time series:

pr.t <- rts(pr,n) 
tmean.t <- rts(tmean,n)
###########################


dv <- dVelocity(pr.t,tmean.t,t1='1991/2000',t2='2010/2020')
names(dv) <- 'dVelocity'

ve <- velocity(x1=pr.t,x2=tmean.t,t1='1991/2000',t2='2010/2020')
names(ve) <- 'velocity'

v <- ccm(pr.t,tmean.t,t1='1991/2000',t2='2010/2020',stat=c('dVelocity','velocity'))


expect_equal(round(sum(dv[],na.rm=T)),3136)

expect_equal(round(sum(ve[],na.rm=T)),46577)

expect_equivalent(dv,v[[1]])

expect_equivalent(ve,v[[2]])
