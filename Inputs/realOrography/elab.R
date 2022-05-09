


rm(list=ls())

library(raster) # To handle raster data.
library(imager) # image processing

#args<-commandArgs(trailingOnly = TRUE)
# args[1]   # number of Conditional Simulations

folder = c('../45013.hpc.mate.polimi.it/Outputs/0','../45014.hpc.mate.polimi.it/Outputs/0', 
           '../45015.hpc.mate.polimi.it/Outputs/0', '../45016.hpc.mate.polimi.it/Outputs/0')
number_sim_days = 365
# Load the Digital Elevation Model (DEM). The resolution, extent and coordinate system will be used as reference.
dem=raster('DEM.tif')
mask=raster('Mask_bin.tif')

resolution=res(raster(paste0(folder[1],'/basin_mask.asc')))[[1]]

crs_ref=crs(dem)

mask_grad_x=mask
mask_grad_y=mask

mask_cimg=as.cimg(mask,maxpixels=1e+07)
mask_grad_cimg=abs(imgradient(mask_cimg,axes="x",scheme=3L))
values(mask_grad_x) = abs(as.vector(mask_grad_cimg)) #as.vector(mask_cimg) #as.vector(mask_grad_cimg[[1]])

mask_grad_cimg=abs(imgradient(mask_cimg,axes="y",scheme=3L))
values(mask_grad_y) = abs(as.vector(mask_grad_cimg)) 

values(mask_grad_x)[is_less_than(values(mask_grad_x),.02)] = NA
values(mask_grad_y)[is_less_than(values(mask_grad_y),.02)] = NA




xx_x=c(530850.173)
xx_y=c(5077721.741)

d <- data.frame(x=c(xx_x), y=c(xx_y))
coordinates(d) <- c("x", "y")
proj4string(d) <- dem@crs # WGS 84

vasca_bassa_x=sum(c(531638, 531627, 531640, 531654))/4
vasca_bassa_y=sum(c(5078678,5078681,5078770,5078764))/4

vasca_alta_x=sum(c(532755, 532764, 532777, 532784))/4
vasca_alta_y=sum(c(5079635,5079630,5079681,5079680))/4

d1 <- data.frame(x=c(vasca_bassa_x,vasca_alta_x), y=c(vasca_bassa_y,vasca_alta_y))
coordinates(d1) <- c("x", "y")
proj4string(d1) <- dem@crs # WGS 84

quartz()
plot(dem,main="orography: b")
plot(d,add=T,col='red')
plot(d1,add=T,col='blue')
plot(mask_grad_x*10000,add=TRUE,col='black',legend=FALSE)
plot(mask_grad_y*10000,add=TRUE,col='black',legend=FALSE)
quartz.save(paste0('dem.png'),type="png")





## bassa
vasca_bassa_x=c(531638, 531627, 531640, 531654)
vasca_bassa_y=c(5078678,5078681,5078770,5078764)

d_bassa <- data.frame(X=vasca_bassa_x, Y=vasca_bassa_y) 
coordinates(d_bassa) <- c("X", "Y")
proj4string(d_bassa) <- crs_ref 
h_bassa=2.20
dem_bassa=dem
values(dem_bassa) = 0.
dem_bassa[extent(d_bassa)] = -h_bassa
#dem_points<-extract(dem_bassa,extent(d_bassa))



## alta
vasca_alta_x=c(532755, 532764, 532777, 532784)
vasca_alta_y=c(5079635,5079630,5079681,5079680)
d_alta <- data.frame(X=vasca_alta_x, Y=vasca_alta_y) 
coordinates(d_alta) <- c("X", "Y")
proj4string(d_alta) <- crs_ref
h_alta=4.8
dem_alta=dem
values(dem_alta) = 0.
dem_alta[extent(d_alta)] = -h_alta
#dem_points<-extract(dem_bassa,extent(d_bassa))

dem_tanks=dem+dem_bassa+dem_alta


diff_tanks=aggregate(dem,resolution/res(dem)[[1]])-aggregate(dem_tanks,resolution/res(dem)[[1]])
dem_tanks_up=aggregate(dem_tanks,resolution/res(dem)[[1]])

quartz()
plot(diff_tanks)


slope=aggregate(terrain(dem_tanks, opt="slope", unit="tangent", neighbors=4),resolution/res(dem)[[1]])
quartz()
plot(slope)

slope[extent(d_alta)]
slope[extent(d_bassa)]


j = 1
sed_alta_1=c()
sed_bassa_1=c()
for(i in 0:number_sim_days) {
  i
  h_sd = raster(paste0(folder[j],'/hsd_',i,'.asc'))
  max(h_sd)
  #u = raster(paste0(folder,'/u_',i,'.asc'))
  #sediment_tanks = h_sd*diff_tanks h_sd[extent(d_alta)]
  #sed_alta  = c(sed_alta,sum(dem_tanks_up[extent(d_alta)]  - sediment_tanks[extent(d_alta)] *res(diff_tanks)[[1]]*res(diff_tanks)[[2]]))
  #sed_bassa = c(sed_bassa,sum(dem_tanks_up[extent(d_bassa)] - sediment_tanks[extent(d_bassa)]*res(diff_tanks)[[1]]*res(diff_tanks)[[2]]))
  
  sed_alta_1  = c(sed_alta_1, sum(h_sd[extent(d_alta)] ) / sum(diff_tanks[extent(d_alta )])*100)
  sed_bassa_1 = c(sed_bassa_1,sum(h_sd[extent(d_bassa)]) / sum(diff_tanks[extent(d_bassa)])*100)
  
  
}


j = 2
sed_alta_2=c()
sed_bassa_2=c()
for(i in 0:number_sim_days) {
  i
  h_sd = raster(paste0(folder[j],'/hsd_',i,'.asc'))
  max(h_sd)
  #u = raster(paste0(folder,'/u_',i,'.asc'))
  #sediment_tanks = h_sd*diff_tanks h_sd[extent(d_alta)]
  #sed_alta  = c(sed_alta,sum(dem_tanks_up[extent(d_alta)]  - sediment_tanks[extent(d_alta)] *res(diff_tanks)[[1]]*res(diff_tanks)[[2]]))
  #sed_bassa = c(sed_bassa,sum(dem_tanks_up[extent(d_bassa)] - sediment_tanks[extent(d_bassa)]*res(diff_tanks)[[1]]*res(diff_tanks)[[2]]))
  
  sed_alta_2  = c(sed_alta_2, sum(h_sd[extent(d_alta)] ) / sum(diff_tanks[extent(d_alta )])*100)
  sed_bassa_2 = c(sed_bassa_2,sum(h_sd[extent(d_bassa)]) / sum(diff_tanks[extent(d_bassa)])*100)
  
  
}



j = 3
sed_alta_3=c()
sed_bassa_3=c()
for(i in 0:number_sim_days) {
  i
  h_sd = raster(paste0(folder[j],'/hsd_',i,'.asc'))
  max(h_sd)
  #u = raster(paste0(folder,'/u_',i,'.asc'))
  #sediment_tanks = h_sd*diff_tanks h_sd[extent(d_alta)]
  #sed_alta  = c(sed_alta,sum(dem_tanks_up[extent(d_alta)]  - sediment_tanks[extent(d_alta)] *res(diff_tanks)[[1]]*res(diff_tanks)[[2]]))
  #sed_bassa = c(sed_bassa,sum(dem_tanks_up[extent(d_bassa)] - sediment_tanks[extent(d_bassa)]*res(diff_tanks)[[1]]*res(diff_tanks)[[2]]))
  
  sed_alta_3  = c(sed_alta_3, sum(h_sd[extent(d_alta)] ) / sum(diff_tanks[extent(d_alta )])*100)
  sed_bassa_3 = c(sed_bassa_3,sum(h_sd[extent(d_bassa)]) / sum(diff_tanks[extent(d_bassa)])*100)
  
  
}


j = 4
sed_alta_4=c()
sed_bassa_4=c()
for(i in 0:number_sim_days) {
  i
  h_sd = raster(paste0(folder[j],'/hsd_',i,'.asc'))
  max(h_sd)
  #u = raster(paste0(folder,'/u_',i,'.asc'))
  #sediment_tanks = h_sd*diff_tanks h_sd[extent(d_alta)]
  #sed_alta  = c(sed_alta,sum(dem_tanks_up[extent(d_alta)]  - sediment_tanks[extent(d_alta)] *res(diff_tanks)[[1]]*res(diff_tanks)[[2]]))
  #sed_bassa = c(sed_bassa,sum(dem_tanks_up[extent(d_bassa)] - sediment_tanks[extent(d_bassa)]*res(diff_tanks)[[1]]*res(diff_tanks)[[2]]))
  
  sed_alta_4  = c(sed_alta_4, sum(h_sd[extent(d_alta)] ) / sum(diff_tanks[extent(d_alta )])*100)
  sed_bassa_4 = c(sed_bassa_4,sum(h_sd[extent(d_bassa)]) / sum(diff_tanks[extent(d_bassa)])*100)
  
  
}


quartz()
plot(sed_alta_2,col="blue",pch=19,type = "b", xlab = "days", ylab = "tank filling percentage", main = "high tank",ylim=c(0,.9))
lines(sed_alta_1,col="blue",pch=19,type = "b", xlab = "days", ylab = "tank filling percentage", main = "high tank")
lines(sed_alta_3,col="red",pch=19,type = "b", xlab = "days", ylab = "tank filling percentage", main = "high tank")
lines(sed_alta_4,col="red",pch=19,type = "b", xlab = "days", ylab = "tank filling percentage", main = "high tank")
legend("topleft", legend=c("35 (m)","50 (m)"),
       col=c("blue", "red"), lty=c(1,1), pch=c(19,19), cex=1.2,bty="n")
quartz.save(paste0('sed_alta.png'),type="png")


quartz()
plot(sed_bassa_2,col="blue",pch=19,type = "b", xlab = "days", ylab = "tank filling percentage", main = "low tank",ylim=c(0,60))
lines(sed_bassa_1,col="blue",pch=19,type = "b", xlab = "days", ylab = "tank filling percentage", main = "low tank")
lines(sed_bassa_3,col="red",pch=19,type = "b", xlab = "days", ylab = "tank filling percentage", main = "low tank")
lines(sed_bassa_4,col="red",pch=19,type = "b", xlab = "days", ylab = "tank filling percentage", main = "low tank")
legend("topleft", legend=c("35 (m)","50 (m)"),
       col=c("blue", "red"), lty=c(1,1), pch=c(19,19), cex=1.2,bty="n")
quartz.save(paste0('sed_bassa.png'),type="png")







quartz()
plot(sed_alta,type = "b", xlab = "days", ylab = "tank filling percentage", main = "high tank")
quartz.save(paste0('sed_alta.png'),type="png")

quartz()
plot(sed_bassa,type = "b", xlab = "days", ylab = "tank filling percentage", main = "low tank")
quartz.save(paste0('sed_bassa.png'),type="png")



cord = SpatialPoints(cbind(9.397794, 45.851888 ), proj4string = CRS("+proj=longlat"))
cord.UTM <- spTransform(cord, crs(dem))

quartz()
plot(dem)
plot(cord.UTM,add=T)
plot(mask_grad_x*10000,add=TRUE,col='black',legend=FALSE)
plot(mask_grad_y*10000,add=TRUE,col='black',legend=FALSE)




