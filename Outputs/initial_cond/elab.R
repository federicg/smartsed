


rm(list=ls())

library(raster) # To handle raster data.
library(imager) # image processing
library(viridis)
library(ggplot2)
library(ggquiver)

#args<-commandArgs(trailingOnly = TRUE)
# args[1]   # number of Conditional Simulations

raster2quiver <- function(dem, uu, vv, colours = terrain.colors(6))
{
  names(dem) <- "z"
  dem_df <- as.data.frame(dem, xy=T)
  quiv <- aggregate(dem, 1)
  quiv$u <- values(uu) 
  quiv$v <- values(vv)
  quiv_df <- as.data.frame(quiv, xy = TRUE)
  
  print(ggplot(mapping = aes(x = x, y = y, fill = z)) + 
          geom_raster(data = dem_df, na.rm = TRUE) + 
          geom_quiver(data = quiv_df, aes(u = u, v = v), vecsize = 1.5) +
          scale_fill_gradientn(colours = colours, na.value = "transparent") +
          theme_bw())
  
  return(quiv_df)
}



# Load the Digital Elevation Model (DEM). The resolution, extent and coordinate system will be used as reference.
dem=raster('../visualize/DEM.tif')
mask=raster('../visualize/Mask_bin.tif')

clay = raster(paste0('clay_0.asc'))
h = raster(paste0('H_ev.asc'))
u = raster(paste0('u_ev.asc'))
v = raster(paste0('v_ev.asc'))

h = resample(h,clay,method="bilinear")
h[values(mask)==0]=NA

quartz()
plot(h)

dem_ = aggregate(dem,7)
#vv=resample(v, dem_, method="bilinear")
#uu=resample(u, dem_, method="bilinear")

uu = resample(u,dem_,method="bilinear")
vv = resample(-v,dem_,method="bilinear")

h_res = resample(h, dem, method="bilinear") 
u_res = resample(u, dem, method="bilinear") 
v_res = resample(v, dem, method="bilinear") 

writeRaster(h_res,file=paste0('h_ev_5.asc'),overwrite=T)
writeRaster(u_res,file=paste0('u_ev_5.asc'),overwrite=T)
writeRaster(v_res,file=paste0('v_ev_5.asc'),overwrite=T)


pal <- c("#B2182B", "#E68469", "#D9E9F1", "#ACD2E5", "#539DC8", "#3C8ABE", "#2E78B5")

#quartz()
raster2quiver(dem_, uu, vv, colours = pal)


u = resample(u,dem,method="bilinear")
v = resample(v,dem,method="bilinear")



quartz()
par(cex.lab=1.5,cex.main=1.5,cex.axis=1.5)
plot(h,main="water height (m)",xlab='x (m)', ylab='y (m)',col=topo.colors(20))
quartz.save(paste0('H.png'),type="png")

vel = sqrt(u^2+v^2)
vel[values(mask)==0]=NA
quartz()
par(cex.lab=1.5,cex.main=1.5,cex.axis=1.5)
plot(vel,main="velocity module (m/sec.)",xlab='x (m)', ylab='y (m)',col=topo.colors(20))
quartz.save(paste0('vel_abs.png'),type="png")




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

# via carlo porta
cord = SpatialPoints(cbind(9.397794, 45.851888 ), proj4string = CRS("+proj=longlat"))
cord.UTM <- spTransform(cord, crs(dem))

quartz()
plot(dem)
plot(cord.UTM,add=T)
plot(mask_grad_x*10000,add=TRUE,col='black',legend=FALSE)
plot(mask_grad_y*10000,add=TRUE,col='black',legend=FALSE)



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






