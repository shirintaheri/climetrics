filePath <- system.file("external/", package="climetrics")

pr <- rast(paste0(filePath,'/precip.tif'))
tmax <- rast(paste0(filePath,'/tmax.tif'))
n <- readRDS(paste0(filePath,'/dates.rds')) # read corresponding dates

####################

# use rts function in the rts package to make a raster time series:

pr.t <- rts(pr,n) 
tmax.t <- rts(tmax,n)

###########################
# test of the metric:

# The extreme argument corresponds to the first and second climate variables
# (i.e., x1 and x2; precipitation and temperature) that specify the percentile of the extreme 
# condition in climate variable; here, 0.05 is used for precipitation; and 0.95 for temperature

se <- sed(pr.t,tmax.t,t1='1991/2000',t2='2010/2020')

expect_equal(round(sum(se[],na.rm=T)),135)

se2 <- ccm(pr.t,tmax.t,t1='1991/2000',t2='2010/2020',stat=c('sed'))

expect_equivalent(se,se2)
