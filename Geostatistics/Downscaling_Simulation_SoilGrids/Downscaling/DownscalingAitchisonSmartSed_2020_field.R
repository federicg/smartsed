
########################################################################################################
########  Geostatistical downscaling (prediction and simulation) with SoilGrids data  ##################
########################################################################################################

# this script only performs Downscaling in the Aitchison geometry


##### 1. Load and process SoilGrids maps (more info at https://www.isric.org/explore/soilgrids) #####
##### Maps can be downloaded from the site https://soilgrids.org/ ###################################
##### or at the link: https://files.isric.org/soilgrids/data/recent/ ################################

### In the working directory there must be:
### - A directory 'SoilGrids' containing SoilGrids' maps of the variable in exam covering the area of interest
### - A raster (.tif or .tiff) file of the Digital Elevation Model for the area in exam.


library(raster) # To handle raster data.
library(sp)
library(rgdal)

args<-commandArgs(trailingOnly = TRUE)
# args[1]   # number of Conditional Simulations

#print("a")
# Load the Digital Elevation Model (DEM). The resolution, extent and coordinate system will be used as reference.
dem=raster('../Inputs/realOrography/DEM.tif')
## Load the slope
slp=raster('../Inputs/SoilData/Slope.tif')
## Load Curvature Files (2 of them)
#Curv1 = raster('../Outputs/SoilData/Curvatura_PF.tif')
#Curv2 = raster('../Outputs/SoilData/Curvatura_PL.tif')
#Curv3 = raster('../Outputs/SoilData/Curvatura_TG.tif')
#Curv4 = raster('../Outputs/SoilData/Curvatura_GEN.tif')

#print(str(dem))

#print(str(slp))
#print(str(Curv1))



#print("b")

# Plot
#plot(dem)

# Coordinates Reference System (CRS) of the dem used as reference.
crs_ref=crs(dem) 
# Reference extent.
extent_ref=extent(dem)  


# Load particle-size fractions (psf) into a stack (raster with multiple data associated to each pixel) 
#psf = stack( c('../Inputs/SoilGrids/CLYPPT_M_sl1_250m.tiff',    # clay
#               '../Inputs/SoilGrids/SLTPPT_M_sl1_250m.tiff',    # silt
#               '../Inputs/SoilGrids/SNDPPT_M_sl1_250m.tiff' ) ) # sand

psf = stack( c('../Inputs/SoilGrids/SoilGrids2020/clay.tif',    # clay
               '../Inputs/SoilGrids/SoilGrids2020/silt.tif',    # silt
               '../Inputs/SoilGrids/SoilGrids2020/sand.tif' ) ) # sand


### The .tiff files are searched in the directory 'SoilGrids'.
### The '1' in 'sl1' in the names of the raster files indicates that the values considered afre those of the topsoil
### (0 cm from the top). Number 2 correspond to 15 cm, and so on up to number 7 following the standard depths 
### (0-15-30-60-100-200cm).

#To load psf at different depths
# psf5cm = stack( c('SoilGrids/CLYPPT_M_sl2_250m.tiff',  # clay
#                'SoilGrids/SLTPPT_M_sl2_250m.tiff',     # silt
#                'SoilGrids/SNDPPT_M_sl2_250m.tiff' ) )  # sand

# Change the names of the layers to CLAY, SILT and SAND.
names(psf) = c('CLAY','SILT','SAND') 

# Change the coordinate system of the stack of SoilGrids maps to match the reference one.
psf=projectRaster(psf, crs = crs_ref) 

# Crop the maps to match the extent of the DEM (SoilGrids maps are originally broader)
psf=crop(psf, extent_ref) #, filename='psf.tiff')  ( The processed data can be saved )
#print(psf)

# Create a Data Frame (df) containing the values of each coarse 'pixel'.
psf_df = as.data.frame(psf) 

# Each row must sum to 100 (%), but there is a residual (res) error due to rounding
res=rep(1000,nrow(psf_df))-apply(psf_df,1,sum,na.rm=TRUE)

# Apply (uniform) correction to the data frame and the rasters. 
for(i in 1:ncol(psf_df)) {
  psf_df[,i] = psf_df[,i] + res/ncol(psf_df)
  values(psf[[i]]) = psf_df[,i] 
}

psf=psf/10
psf_df=psf_df/10

#### Add reading of OGR files for field data 
library(raster)
#library(imager)

shape <- readOGR(dsn = "../Inputs/SoilData", layer = "SoilTexture_complete")
ID_points<-shape@data$ID[1:73]
field_coords <-shape@coords[1:73,]
#print(ID_points)
psf_field=as.data.frame(ID_points)
psf_field['coarse']=shape@data$X.Coarse[1:73]
psf_field['sand']=shape@data$X.Sand[1:73]
psf_field['fines']=shape@data$X.Fines[1:73]

tot<-psf_field$sand+psf_field$fines#+psf_field$coarse
#tot
#psf_field
psf_field$coarse=psf_field$coarse/tot
psf_field$sand=psf_field$sand/tot
psf_field$fines=psf_field$fines/tot

#psf_field
psf_field$ID_points=as.data.frame(ID_points)
#print("psf_field")
#print(typeof(psf_field))

#psf_field

coords=shape@coords[1:73,]
#print(coords)
#print(psf_field)
##### Plots #####

# viridis colors (for pretty plots)
#library(viridis)

### Functions ###

# Plot the histogram (with optional fitted density) of raster values.
plot_raster_hist = function(rstr, length = 40, colors = terrain.colors(40), title = "", 
                            density_col = "blue", xlim = NA, fit.density = TRUE) {
  
  ### Arguments:
  # rstr = raster map
  # length = number of intervals 
  # density_col = color of the denisty line
  # fit.density = boolean, if true fits a density and plots it
  
  vals = values(rstr)
  min_val <-min(vals, na.rm = TRUE)
  max_val <-max(vals, na.rm = TRUE)
  if(is.na(xlim)){xlim = c(min_val,max_val)}
  
  hist(vals, breaks=seq(min_val,max_val,length = length), 
       col=colors, xlab="",
       main=title,
       xlim=xlim, 
       freq=FALSE) 
  if(fit.density) { lines(density(vals, na.rm=TRUE), col = density_col, lwd = 2) }
}

# Plot a map and the histogram of its values on the side. 
#quartz(width = 10, height = 5) ### For windows it is not 'quartz()' but 'x11()' !!! (consider using 'replace all').
#par(mfrow=c(1,2),mai = c(1,1,1,1)) # For multiple plots in the same figure. mai = margins
#plot(psf$CLAY, main = "Map of clay %")
#plot_raster_hist(psf$CLAY, title = "Histogram of clay %") # Function defined above

# Plot two rasters side by side.
#quartz(width = 10, height = 5)
#par(mfrow=c(1,2),mai = c(1,1,1,1))
#plot(CDB, main = "Map of CDB (cm)", col = viridis(20))
#plot(PRH, main = "Map of PRH (%)", col = viridis(20))


##### 2. Analysis of soil texture  ###############################################################
##### https://cran.r-project.org/web/packages/soiltexture/vignettes/soiltexture_vignette.pdf #####

#library(soiltexture)

# Plot psf on the soil texture triangle with the different soil types of the USDA classification.
#geo = TT.geo.get()
#quartz(width = 6, height = 6)
#TT.plot(class.sys = 'USDA.TT')
#TT.plot(class.sys = 'USDA.TT', tri.data = psf_df[-which(is.na(psf_df$CLAY)),], geo = geo, #grid.show = FALSE, 
#        col = "blue", cex=.5, lwd=.5) #, add = FALSE)

# Get classes (one-hot encoding) for each observation.
#soil_type_one_hot= TT.points.in.classes(class.sys = 'USDA.TT', tri.data = psf_df[-which(is.na(psf_df$CLAY)),])
#classes = colnames(soil_type_one_hot)
# Get numeric encoding.
#soil_type = apply(soil_type_one_hot,1,FUN = function(x) as.vector(x) %*% 1:ncol(soil_type_one_hot) )

# To get the name of the class of the pixel in position i
#classes[soil_type[i]]

# Map of soil types (numerical code)
#soil_type_map = psf$CLAY
#values(soil_type_map)[-which(is.na(psf_df$CLAY))] = soil_type
#plot(soil_type_map)

# Get modal soil type over the area.
#Mode = function(x) {
#  ux = unique(x)
#  ux[which.max(tabulate(match(x, ux)))]
#}
#classes[Mode(soil_type)]


##### 3. Downscaling and simulation of psf through Isometric Log-Ratio Area-to-Point Regression Kriging (ILR_ATPRK) ###

### Functions ###################
compile_polygons = function(zField,xField,yField){
  for (i in 1:length(zField)){
  xoff = 10
  yoff = 10
#  print(i)
    mydf <-data.frame(x = c(xField[i]-xoff,xField[i]-xoff,xField[i]+xoff,xField[i]+xoff),
                  y = c(yField[i]-yoff,yField[i]+yoff,yField[i]-yoff,yField[i]+yoff),
                  z = c(zField[i],zField[i],zField[i],zField[i]))
    #print(str(mydf))
    myspg <-mydf
    coordinates(myspg)<-~x+y
    gridded(myspg) <-TRUE
    rasterDF <-raster(myspg)
    #rasterDF 
    single_poly = rasterToPolygons(rasterDF,fun=NULL,n=16, na.rm=TRUE, digits=12, dissolve=FALSE)
    #single_poly
    if (i ==1){
      tot_poly<-single_poly
    }else{
      tot_poly<-(rbind(tot_poly,single_poly,makeUniqueIDs = TRUE))
    }

  }
  #print("tot_poly")
  return(tot_poly)
}

variogram_deconvolution = function(coarse_raster, vg_type = "Sph", nblocks = 4, maxiter = 100, cutoff = "default",
                                   tol1 = 1e-2, tol2 = 1e-6) { 
#variogram_deconvolution = function(coarse_raster, vg_type = "Sph", nblocks = 4, maxiter = 100, cutoff = "default",
  
  ##### Variogram deconvolution procedure (Goovaerts, 2008) #########
  ##### This is a special version for regular grids (raster data) ###
  
  require(gstat)
  require(fields)
  
  ########### Arguments: ###################################################################################
  
  # coarse_raster = coarse resolution raster. 
  # vg_type = type of variogram to be fitted.
  # nblocks = the regularized variogram is computed using R function vgmArea, but only at short distances,
  #           nblocks is the number of adjacent blocks that are considered near enough to require vgmArea.
  # maxiter = maximum number of iterations.
  # cutoff = cutoff. If = "default" it is set equal to half the raster extent.
  # tol1 = tolerance for Di/D0, if the value is lower the iterations stop.
  # tol2 = tolerance for abs(D_i-D_opt)/D_opt.
  
  ##########################################################################################################
  
  # Borders and extent of the raster map.
  xmin = coarse_raster@extent@xmin
  xmax = coarse_raster@extent@xmax
  ymin = coarse_raster@extent@ymin
  ymax = coarse_raster@extent@ymax
  coarse_res = res(coarse_raster)
  #print("res")
  #print(coarse_res)
  x_extent = xmax - xmin
  y_extent = ymax - ymin
  
  if (cutoff == "default") {
    cutoff = min(c(x_extent,y_extent))/2 # cutoff equal to half the raster extent by default.
  }
  
  # Change the name to a generic "z".
  names(coarse_raster) = 'z' 
  
  # Convert the raster to a list of square polygons corresponding to the coarse pixels (=blocks).
  poly_coarse = rasterToPolygons(coarse_raster, fun=NULL, n=16, na.rm=TRUE, digits=12, dissolve=FALSE)
  #print(str(poly_coarse))
  ##### Lettura di un raster aggiuntivo raster_field
  #poly_coarse += rasterToPolygons(fine_raster_field, fun=NULL, n=16, na.rm=TRUE, digits=12, dissolve=FALSE)
  #print("Variogram")
  #print(str(coarse_raster))
  
   ##QUI 
  #poly_coarse<-rbind(poly_coarse,point_poly)
  #print(str(newpoly_coarse,max.level=2))
  
  # 1. Compute empirical variogram on areal data and fit a model.
  gamma_v_hat = variogram(z ~ 1, data = poly_coarse, cutoff=cutoff, width = min(coarse_res)/2)
  #print("HERE")
  #print(str(gamma_v_hat))
  gamma_v_exp = fit.variogram(gamma_v_hat, vgm(mean(tail(gamma_v_hat$gamma),10), vg_type, 
                                               cutoff/2, min(gamma_v_hat$gamma)))
  
  #print("HERE2")
  #print(str(gamma_v_exp))
  # 2. Initial variogram.
  gamma_0 = gamma_v_exp
  
  # 3. Variogram regularization using 'vgmArea' for low distances, and Journel approximation for high distances.
  # This procedure is specific for regular grids and allows to significantly speed up the algorithm,
  # avoiding the repeated computation of block covariances, or the computation of block covariances at great
  # distances, which can be approximated using the Journel formula (Journel, 1978).
  # help(vgmArea)
  #new_xmax = xmin + nblocks*coarse_res[1]
  #new_ymin = ymax - nblocks*coarse_res[2]
  new_xmax = xmin + nblocks*coarse_res[1]
  new_ymin = ymax - nblocks*coarse_res[2]
  #print("here3")
  inrange = crop(coarse_raster, extent(c(xmin,new_xmax,new_ymin,ymax)))
  #print("here4")
  #print(inrange)
  poly_ref = polygons(rasterToPolygons(inrange, fun=NULL, n=16, na.rm=TRUE, digits=12, dissolve=FALSE))
  #print("here5")
  #print(str(poly_ref))
  # Now poly_ref contains a square of blocks, given the regularity of the problem we can drop the upper part
  # of the square since covariance only depends on the distance of the blocks.
  #print("here6")
  polygon_indexes = upper.tri(matrix(1:(length(poly_ref)),nrow = nblocks, ncol = nblocks), diag = TRUE)
  #print("here7")
  polygon_indexes = as.vector(polygon_indexes)
  #print("here8")
  #print(polygon_indexes)
  poly_ref = poly_ref[polygon_indexes]
  #print("here9")
  # Compute gamma^(v,v_h) for small lags. 
  gamma_A_0 = vgmArea(x = poly_ref[1],y = poly_ref, vgm = gamma_0, covariance = FALSE)
  #print("here10")
  gamma_vv_0 = gamma_A_0[1] # gamma(v,v), unique since the grid is regular.
  coords_ref = as.matrix(coordinates(poly_ref))
  short_dist =  as.vector( rdist(t(coords_ref[1,]),coords_ref[2:nrow(coords_ref),]) ) # Distances of neighboring blocks.
  drop_duplicates = !duplicated(short_dist)
  # Drop distance duplicates (the block covariance only depend on distance).
  short_dist = short_dist[drop_duplicates]
  ordered = order(short_dist)
  short_dist = short_dist[ordered]
  gamma_v_0 = as.vector( gamma_A_0[2:length(gamma_A_0)] )
  gamma_v_0 = gamma_v_0[drop_duplicates]
  gamma_v_0 = gamma_v_0[ordered] # Values of regularized variogram for small distances
  # Add lags (greater distances).
  dist_tail = seq(short_dist[length(short_dist)]+max(coarse_res),cutoff,max(coarse_res))
  ref_dist = c(short_dist, dist_tail)
  gamma_v_0_tail = variogramLine(gamma_0,dist_vector=dist_tail)$gamma
  gamma_v_0 = c(gamma_v_0,gamma_v_0_tail)
  gamma_v_0 = gamma_v_0 - gamma_vv_0 # Regularize
  
  # 4. Quantify deviation.
  D_0 = ( 1/length(ref_dist) ) * 
    sum( abs(gamma_v_0-variogramLine(gamma_v_exp,dist_vector=ref_dist)$gamma)/
           variogramLine(gamma_v_exp,dist_vector=ref_dist)$gamma )
  
  # 5. Define initial optimal variograms.
  
  gamma_opt = gamma_0
  gamma_v_opt = gamma_v_0
  D_opt = D_0
  rescaling_flag = 0
  
  # Start loop.
  for (i in 1:maxiter){
    print(paste0("iter: ", i))
    
    # 6. Compute experimental values for the new point support semivariogram through a rescaling of the optimal 
    #    point support model. 
    if (!rescaling_flag){
      w = 1 + ( 1/(gamma_v_exp$psill[2]* (i^(1/2)) ) )*
        (variogramLine(gamma_v_exp, dist_vector=ref_dist)$gamma - gamma_v_opt)
    }
    else {
      w = 1+(w-1)/2
      rescaling_flag = 0
    }
    
    # Empirical values to which a variogram model must be fitted.
    gamma_hat_i_values = variogramLine(gamma_opt,dist_vector=ref_dist)$gamma * w
    gamma_hat_i = gamma_v_hat[1:length(ref_dist),]
    gamma_hat_i$np=rep(1,length(ref_dist))
    gamma_hat_i$dist = ref_dist 
    gamma_hat_i$gamma = gamma_hat_i_values
    gamma_hat_i$dir.hor[is.na(gamma_hat_i$dir.hor)] = 0
    gamma_hat_i$dir.ver[is.na(gamma_hat_i$dir.ver)] = 0
    gamma_hat_i$id[is.na(gamma_hat_i$id)] = gamma_hat_i$id[1]
    
    # 7. Fit a new model.
    gamma_i = fit.variogram(object = gamma_hat_i, vgm(mean(tail(gamma_hat_i$gamma),2), vg_type, 
                                                      cutoff, min(gamma_hat_i$gamma)))
    
    # 8. Regularize gamma_i.
    gamma_A_i = vgmArea(x = poly_ref[1],y = poly_ref, vgm = gamma_i, covariance = FALSE)
    gamma_vv_i = gamma_A_i[1]
    gamma_v_i = as.vector( gamma_A_i[2:length(gamma_A_i)] )
    gamma_v_i = gamma_v_i[drop_duplicates]
    gamma_v_i = gamma_v_i[ordered]
    
    gamma_v_i_tail = variogramLine(gamma_i,dist_vector=dist_tail)$gamma
    gamma_v_i = c(gamma_v_i,gamma_v_i_tail)
    gamma_v_i = gamma_v_i - gamma_vv_i
    
    # 9. Compute D_i.
    D_i= ( 1/length(ref_dist) ) * 
      sum( abs(gamma_v_i-variogramLine(gamma_v_exp,dist_vector=ref_dist)$gamma)/
             variogramLine(gamma_v_exp,dist_vector=ref_dist)$gamma )
    
    # 10. Stopping criteria. 
    if ( (D_i/D_0 < tol1) || abs(D_i-D_opt)/D_opt < tol2) break
    if (D_i < D_opt) {
      gamma_opt = gamma_i
      gamma_v_opt = gamma_v_i
      D_opt = D_i
    }
    else {
      rescaling_flag = 1
    }
    
  }
  return(gamma_opt)
}

variogram_deconvolution_field = function(coarse_raster, vg_type = "Sph", nblocks = 4, maxiter = 100, cutoff = "default",
                                   tol1 = 1e-2, tol2 = 1e-6) { 
#variogram_deconvolution = function(coarse_raster, vg_type = "Sph", nblocks = 4, maxiter = 100, cutoff = "default",
  
  ##### Variogram deconvolution procedure (Goovaerts, 2008) #########
  ##### This is a special version for regular grids (raster data) ###
  
  require(gstat)
  require(fields)
  
  ########### Arguments: ###################################################################################
  
  # coarse_raster = coarse resolution raster. 
  # vg_type = type of variogram to be fitted.
  # nblocks = the regularized variogram is computed using R function vgmArea, but only at short distances,
  #           nblocks is the number of adjacent blocks that are considered near enough to require vgmArea.
  # maxiter = maximum number of iterations.
  # cutoff = cutoff. If = "default" it is set equal to half the raster extent.
  # tol1 = tolerance for Di/D0, if the value is lower the iterations stop.
  # tol2 = tolerance for abs(D_i-D_opt)/D_opt.
  
  ##########################################################################################################
  
  # Borders and extent of the raster map.
  xmin = coarse_raster@extent@xmin
  xmax = coarse_raster@extent@xmax
  ymin = coarse_raster@extent@ymin
  ymax = coarse_raster@extent@ymax
  coarse_res = res(coarse_raster)
  #print("res")
  #print(coarse_res)
  x_extent = xmax - xmin
  y_extent = ymax - ymin
  
  if (cutoff == "default") {
    cutoff = min(c(x_extent,y_extent))#/2 # cutoff equal to half the raster extent by default. RED equal to all, since we are working with few isolated points
  }
  
  # Change the name to a generic "z".
  names(coarse_raster) = 'z' 
  
  # Convert the raster to a list of square polygons corresponding to the coarse pixels (=blocks).
  poly_coarse = rasterToPolygons(coarse_raster, fun=NULL, n=16, na.rm=TRUE, digits=12, dissolve=FALSE)
  #print(str(poly_coarse))
  ##### Lettura di un raster aggiuntivo raster_field
  #poly_coarse += rasterToPolygons(fine_raster_field, fun=NULL, n=16, na.rm=TRUE, digits=12, dissolve=FALSE)
  #print("Variogram")
  #print(str(coarse_raster))
  
   ##QUI 
  #poly_coarse<-rbind(poly_coarse,point_poly)
  #print(str(newpoly_coarse,max.level=2))
  
  # 1. Compute empirical variogram on areal data and fit a model.
  gamma_v_hat = variogram(z ~ 1, data = poly_coarse, cutoff=cutoff, width = min(coarse_res)/2)
  #print("HERE")
  #print(str(gamma_v_hat))
  gamma_v_exp = fit.variogram(gamma_v_hat, vgm(mean(tail(gamma_v_hat$gamma),10), vg_type, 
                                               cutoff/2, min(gamma_v_hat$gamma)))
  
  #print("HERE2")
  #print(str(gamma_v_exp))
  # 2. Initial variogram.
  gamma_0 = gamma_v_exp
  
  # 3. Variogram regularization using 'vgmArea' for low distances, and Journel approximation for high distances.
  # This procedure is specific for regular grids and allows to significantly speed up the algorithm,
  # avoiding the repeated computation of block covariances, or the computation of block covariances at great
  # distances, which can be approximated using the Journel formula (Journel, 1978).
  # help(vgmArea)
  #new_xmax = xmin + nblocks*coarse_res[1]
  #new_ymin = ymax - nblocks*coarse_res[2]
  new_xmax = xmin + nblocks*coarse_res[1]
  new_ymin = ymax - nblocks*coarse_res[2]
  #print("here3")
  inrange = crop(coarse_raster, extent(c(xmin,new_xmax,new_ymin,ymax)))
  #print("here4")
  #print(inrange)
  poly_ref = polygons(rasterToPolygons(inrange, fun=NULL, n=16, na.rm=TRUE, digits=12, dissolve=FALSE))
  #print("here5")
  #print(str(poly_ref))
  # Now poly_ref contains a square of blocks, given the regularity of the problem we can drop the upper part
  # of the square since covariance only depends on the distance of the blocks.
  #print("here6")
  #polygon_indexes = upper.tri(matrix(1:(length(poly_ref)),nrow = nblocks, ncol = nblocks), diag = TRUE)
  ##### Parziale modifica sezione
  ##polygon_indexes = matrix(1:(length(poly_ref)))
  ##print("here7")
  ##polygon_indexes = as.vector(polygon_indexes)
  ##print("here8")
  ##print(polygon_indexes)
  ### Prendendoli tutti posso anche non selezionarlo?
  ##poly_ref = poly_ref[polygon_indexes]
  #### Fine sezione
  ##print("here9")
  # Compute gamma^(v,v_h) for small lags. 
  gamma_A_0 = vgmArea(x = poly_ref[1],y = poly_ref, vgm = gamma_0, covariance = FALSE)
  #print("here10")
  gamma_vv_0 = gamma_A_0[1] # gamma(v,v), unique since the grid is regular.
  coords_ref = as.matrix(coordinates(poly_ref))
  short_dist =  as.vector( rdist(t(coords_ref[1,]),coords_ref[2:nrow(coords_ref),]) ) # Distances of neighboring blocks.

  #print("here11")
  drop_duplicates = !duplicated(short_dist)
  # Drop distance duplicates (the block covariance only depend on distance).
  short_dist = short_dist[drop_duplicates]
  ordered = order(short_dist)
  short_dist = short_dist[ordered]
  gamma_v_0 = as.vector( gamma_A_0[2:length(gamma_A_0)] )
  gamma_v_0 = gamma_v_0[drop_duplicates]
  gamma_v_0 = gamma_v_0[ordered] # Values of regularized variogram for small distances
  # Add lags (greater distances).
  #print("here12")
  #print(short_dist[length(short_dist)]+max(coarse_res))
  #print(cutoff)
  #print(max(coarse_res))
  dist_tail = seq(short_dist[length(short_dist)]+max(coarse_res),cutoff,max(coarse_res))
  ref_dist = c(short_dist, dist_tail)

  gamma_v_0_tail = variogramLine(gamma_0,dist_vector=dist_tail)$gamma
  gamma_v_0 = c(gamma_v_0,gamma_v_0_tail)
  gamma_v_0 = gamma_v_0 - gamma_vv_0 # Regularize
  
  # 4. Quantify deviation.
  D_0 = ( 1/length(ref_dist) ) * 
    sum( abs(gamma_v_0-variogramLine(gamma_v_exp,dist_vector=ref_dist)$gamma)/
           variogramLine(gamma_v_exp,dist_vector=ref_dist)$gamma )
  
  # 5. Define initial optimal variograms.
  
  gamma_opt = gamma_0
  gamma_v_opt = gamma_v_0
  D_opt = D_0
  rescaling_flag = 0
  
  # Start loop.
  for (i in 1:maxiter){
    print(paste0("iter: ", i))
    
    # 6. Compute experimental values for the new point support semivariogram through a rescaling of the optimal 
    #    point support model. 
    if (!rescaling_flag){
      w = 1 + ( 1/(gamma_v_exp$psill[2]* (i^(1/2)) ) )*
        (variogramLine(gamma_v_exp, dist_vector=ref_dist)$gamma - gamma_v_opt)
    }
    else {
      w = 1+(w-1)/2
      rescaling_flag = 0
    }
    
    # Empirical values to which a variogram model must be fitted.
    gamma_hat_i_values = variogramLine(gamma_opt,dist_vector=ref_dist)$gamma * w
    gamma_hat_i = gamma_v_hat[1:length(ref_dist),]
    gamma_hat_i$np=rep(1,length(ref_dist))
    gamma_hat_i$dist = ref_dist 
    gamma_hat_i$gamma = gamma_hat_i_values
    gamma_hat_i$dir.hor[is.na(gamma_hat_i$dir.hor)] = 0
    gamma_hat_i$dir.ver[is.na(gamma_hat_i$dir.ver)] = 0
    gamma_hat_i$id[is.na(gamma_hat_i$id)] = gamma_hat_i$id[1]
    
    # 7. Fit a new model.
    gamma_i = fit.variogram(object = gamma_hat_i, vgm(mean(tail(gamma_hat_i$gamma),2), vg_type, 
                                                      cutoff, min(gamma_hat_i$gamma)))
    
    # 8. Regularize gamma_i.
    gamma_A_i = vgmArea(x = poly_ref[1],y = poly_ref, vgm = gamma_i, covariance = FALSE)
    gamma_vv_i = gamma_A_i[1]
    gamma_v_i = as.vector( gamma_A_i[2:length(gamma_A_i)] )
    gamma_v_i = gamma_v_i[drop_duplicates]
    gamma_v_i = gamma_v_i[ordered]
    
    gamma_v_i_tail = variogramLine(gamma_i,dist_vector=dist_tail)$gamma
    gamma_v_i = c(gamma_v_i,gamma_v_i_tail)
    gamma_v_i = gamma_v_i - gamma_vv_i
    
    # 9. Compute D_i.
    D_i= ( 1/length(ref_dist) ) * 
      sum( abs(gamma_v_i-variogramLine(gamma_v_exp,dist_vector=ref_dist)$gamma)/
             variogramLine(gamma_v_exp,dist_vector=ref_dist)$gamma )
    
    # 10. Stopping criteria. 
    if ( (D_i/D_0 < tol1) || abs(D_i-D_opt)/D_opt < tol2) break
    if (D_i < D_opt) {
      gamma_opt = gamma_i
      gamma_v_opt = gamma_v_i
      D_opt = D_i
    }
    else {
      rescaling_flag = 1
    }
    
  }
  return(gamma_opt)
}

ATPlm = function(coarse_raster, covariates_stack){
  
  ### Area-to-Point Linear Regression.
  ### Remark: the target coarse resolution raster and the covariate maps must have the same coordinate system
  ### and must overlap.
  ### The function returns an object containing a fine resolution map (matching the resolution of the covariates)
  ### of the fitted target variable and an lm-type object containig the regression results.
  
  df = as.data.frame(coarse_raster[[1]])
  xnames = names(covariates_stack) # Covariate names.
  yname = names(coarse_raster)[1] # Name of target variable to be used in the formula.
  nx = length(xnames)
  #print("covariates")
  #print(str(covariat2s_stack))
  #print("coarse_raster")
  #print(str(coarse_raster))
  # Upscale covariate maps for the regression.
  for (i in 1:nx){
    df[xnames[i]] = values( resample(covariates_stack[[i]],coarse_raster, method = "bilinear") )
  }
  f <- as.formula( paste(yname, "~ .") )
  #print(nx,str(covariates_stack))
  LM = lm(formula = f, data = df) # Linear regression
  #print("LM")
  #print(str(LM))
  predicted = predict(LM, as.data.frame(covariates_stack))
  map = covariates_stack[[1]]
  values(map) = predicted
  out = list(map = map, LM = LM)

  return(out)
}

ATPlm_field = function(coarse_raster,point_val,point_pos, covariates_stack){
  
  ### Area-to-Point Linear Regression.
  ### Remark: the target coarse resolution raster and the covariate maps must have the same coordinate system
  ### and must overlap.
  ### The function returns an object containing a fine resolution map (matching the resolution of the covariates)
  ### of the fitted target variable and an lm-type object containig the regression results.
  
  #df = as.data.frame(coarse_raster[[1]])
  #df2 = as.data.frame(point_raster[[1]])
  #print("df!!!")
  #print(str(df))
  xnames = names(covariates_stack) # Covariate names.
  #print("xnames")
  #print(xnames)
  #print("ynames")
  #print(names(coarse_raster))
  yname = names(coarse_raster)[1] # Name of target variable to be used in the formula.
  df <- data.frame(CLAY= double(),DEM.1=double(),DEM.2=double())
  nx = length(xnames)
  #for (i in 1:nx){
  #df <- data.frame(xnames[i]= double())
  #}
  #print("df!!!")
  #print(df)
  #print("covariates")
  #print(str(covariates_stack))
  #print("coarse_raster")
  #print(str(coarse_raster))
  # Upscale covariate maps for the regression.
  dataList <- cellFromXY(covariates_stack, point_pos)
  valval <-values(covariates_stack)
  #print(str(valval))
  #print(dataList)
  mydem<-valval[dataList,]
  #print(mydem)
  #df<-rbind(df, c(mydem[,xnames[i]]))
  #print(point_val)
  df<- cbind(point_val,mydem)
  #print("RES!!!")
  colnames(df)<-c("CLAY","DEM.1","DEM.2")
  df = as.data.frame(df)
  #print(df)
  #for (i in 1:nx){
  #  print(xnames[i])
  #  print(mydem[,xnames[i]])
    #df[xnames[i]] = c(mydem[,xnames[i]])
    #df[xnames[i]] = values( resample(covariates_stack[[i]],coarse_raster, method = "bilinear") )
    #df2[xnames[i]] =values( resample(covariates_stack[[i]],point_raster, method = "bilinear") )
  #}
  #print("df2!!")
  #print(df)
  f <- as.formula( paste(yname, "~ .") )
  #print(nx,str(covariates_stack))
  LM = lm(formula = f, data = df) # Linear regression
  #print("LM")
  #print(str(LM))
  predicted = predict(LM, as.data.frame(covariates_stack))
  map = covariates_stack[[1]]
  values(map) = predicted
  out = list(map = map, LM = LM)

  return(out)
}
## Per trattare separatamente i dati di campo
#ATPK_field = function(coarse_raster,point_list, fine_raster,ukField1,ukField2, npoints = 8, nsim = 0, beta = NA, nmax = 10,
#                frac = 1.0, noise_sd = 0.0,vg_type = "Sph"){
ATPK_field = function(coarse_raster,point_list, fine_raster, npoints = 8, nsim = 0, beta = NA, nmax = 10,
                frac = 1.0, noise_sd = 0.0,vg_type = "Sph"){
  
  ### Area-to-Point kriging for the downscaling of coarse raster data.
  ### Converts the coarse pixels into spatial polygons and uses the ATPK method of Kyriakydis (2004).
  ### The function also performs Block Sequential Gaussian Simulation (BSGS).
  ### Returns the downscaled map or a stack of the simulations.
  
  require(gstat)
  
  ### Arguments: 
  # coarse_raster = coarse resolution raster to be downscaled.
  # fine_raster = fine resolution raster covering the extent of the coarse raster. It will be used to define 
  #               the target resolution for the downscaling. The values in this raster are not used.
  # frac = fraction of coarse data to use when performing block sequential simulation (increases the variability)
  # npoints = number of points used to approximate the blocks (coarse pixels). Can be 4,8 or 16.
  # For the other parameters consult the documentation of 'krige' function - help(krige) 
  
  # If simulations are required, consider only the specified fraction of areal data.
  #if (nsim > 0){
  #  drop = 1.0-frac
    #print(length(point_list))
    #print(drop)
 #   ndrops = floor(drop*length(point_list))
    #print("ndrops")
 #   print(ndrops)
    #print(sample(1:length(point_list),ndrops,replace=FALSE))
 #   point_list@data$Z[sample(1:length(point_list),ndrops,replace = FALSE)] = NA

  #}
  #print("jajaja")
  #print(point_list@data)
  
  #BlockMap = rasterToPolygons(coarse_raster, fun=NULL, n=npoints, na.rm=TRUE, digits=12, dissolve=FALSE)
  
#### PArte DAniele

  #coarse_coord<- coordinates(coarse_raster)
  #names(BlockMap) = "Z"
  #names(point_poly) <-names(BlockMap)
  #centpoint<-gCentroid(point_poly,byid=TRUE)

  #point_poly@data<-cbind(point_poly@data,x=centpoint@coords[,1])
  #point_poly@data<-cbind(point_poly@data,y=centpoint@coords[,2])
  #print("point_poly")
  #print(str(point_poly@data))
  #print(names(point_poly))
  #print(names(BlockMap))
  #print(proj4string(point_poly))
  #print(identicalCRS(poly_coarse,point_poly))

  #library(rgeos)
  #centroid<- gCentroid(BlockMap,byid=TRUE)




  #xcol = centroid@coords[,1]
  #print(str(xcol))
  #BlockMap@data<-cbind(BlockMap@data, x=centroid@coords[,1])
  #BlockMap@data<-cbind(BlockMap@data, y=centroid@coords[,2])
  #### The following lines add the polygons of field data ###########
  #BlockMap<- rbind(BlockMap,point_poly,makeUniqueIDs = TRUE)
  ###################################################################
  #print(str(BlockMap@data,max.level =2))
  
  # Add optional noise to the block data
  if(noise_sd){
    print("Adding an uniform noise to the block data")
    noise = rnorm(length(point_list),mean = 0, sd = noise_sd)
    print(noise)
    point_list@data$Z = point_list@data$Z + noise
  }
  #print(str(point_list))
  #print(str(na.omit(point_list)))
  #point_list = na.omit(point_list)
  #### Get the resampled values on the 
  #fine_ukf1 = resample(ukField1, fine_raster, method="bilinear")
  #fine_ukf2 = resample(ukField2, fine_raster, method="bilinear")
  #fine_val_ukf1 = values(fine_ukf1)
  #fine_val_ukf2 = values(fine_ukf2)
  
  d_new <-coordinates(fine_raster)[!is.na(values(fine_raster)),]#, proj4string = crs(BlockMap)
  ddf_new = as.data.frame(d_new)

  coordinates(ddf_new) <- c('x','y')
  ##ddf_new$c1 <- fine_val_ukf1
  ##ddf_new$c2 <- fine_val_ukf2
  #print("ddf_new")
  #d_new <-SpatialPoints(coordinates(fine_raster))
  #print(str(ddf_new))
  #write.csv(ddf_new,"grid.csv")
  #print(str(ddf_new))
  mydata <- point_list
  #xycoords <- mydata@coords
  #xyid <- cellFromXY(ukField1,xycoords)
  #f1Vals = values(ukField1)
  #f2Vals = values(ukField2)
  #print(max(xyid))
  #print(min(xyid))
  #xyF1Vals <-f1Vals[xyid]
  #xyF2Vals <-f2Vals[xyid]
  #print(xyF1Vals)
  #print(xyF2Vals)
  #mydata@data<-cbind(mydata@data,c1=xyF1Vals)
  #mydata@data<-cbind(mydata@data,c2=xyF2Vals)
  #print(str(mydata))
  #d_new <- gstat(formula = Z~1,data= point_list)

  #print(str(d_new))
  #print("d2_new")
  #print(str(d2_new))
  #crs(d_new) <- crs(BlockMap)
  ## Ordinary kriging
  variogram <- variogram(Z~1, data = point_list)
  ## Universal kriging
  #variogram <- variogram(Z~c1+c2, data = mydata)
  #png(filename=paste('var',vg_type,'.png',sep=""))
  #fig1 <- plot(variogram)
  #print(fig1)
  #dev.off()

  varCut = min(max(d_new[,1])-min(d_new[,1]),(max(d_new[,2]-min(d_new[,2])))) #(d_new$data$var1$data@bbox[2,]-d_new$data$var1$data@bbox[1,])
  #varCut = min(extent[1],extent[2])

  #print("variogram")
  #print(variogram)
  #print(str(variogram))
  #print("Cut!!!")
  #print(varCut)

  ### Empirical factors on varSill, varNugg
  ### 
  fVarSill =  1
  fVarNugg = 0.1
  ###
  #########
  varSill <- mean(tail(variogram$gamma,2))*fVarSill
  varNugg <- min(variogram$gamma)*fVarNugg
  vgmMod <- vgm(psill=varSill,model = vg_type,range=0.1*varCut,nugget=varNugg)
  varMod <- fit.variogram(variogram,vgmMod)

  #png('var2.png')

  #png(filename=paste('var2',vg_type,'.png',sep=""))
  #fig2 <- plot(variogram,varMod)
  #print(fig2)
  #dev.off()
  #print("HERE!!!")

  #print(str(vgmMod))
  #png(paste("variogram.png"))
  #plot(variogram,model=varMod)
  #dev.off()
  #prinp(str(BlockMap,max.level=3))
  #browseURL(paste("variogram",vg_type,".png"))
  #print(proj4string(BlockMap))
  #print("Names")
  #mydata2 = SpatialPoints(coordinates(mydata), proj4string = crs(BlockMap))
  #crs(mydata) = crs(BlockMap)
  #print(names(coordinates(BlockMap)))
  #print(attr(coordinates(BlockMap),"dimnames"))
  #print(str(coordinates(BlockMap)))
  #print(proj4string(BlockMap))
  #print(attributes(d_new))
  #print(str(mydata))
  if(is.na(beta)) {
    ### Ordinary
    downscale = krige(Z~1, locations= mydata, newdata = ddf_new, model = varMod, nsim=nsim)
    ## Universal 
    #downscale = krige(Z~c1+c2, mydata, newdata = ddf_new, model = varMod, nsim=nsim)
    #name= paste('Pred-krige',vg_type,'.png',sep="")
    #png(name)
    #fig=spplot(downscale)
    #print(fig)
    #dev.off()
  } else {
    ## Ordinary
    downscale = krige(Z~1, locations = mydata, newdata = ddf_new, model = varMod, nmax = nmax, nsim=nsim, beta = beta)
    ### Universal
    #downscale = krige(Z~c1+c2, mydata, newdata = mydata, model = varMod, nmax = nmax, nsim=nsim, beta = beta)
  }
  if (nsim == 0) {
    downscaled_map = fine_raster
    values(downscaled_map)[!is.na(values(fine_raster))] = downscale$var1.pred
    return(downscaled_map)
  } else {
    sim_data = downscale@data 
    sims = stack(replicate(nsim,fine_raster))
    for (i in 1:nsim) values(sims[[i]])[!is.na(values(fine_raster))] = sim_data[,i]
    return(sims)
  }
}  

#################################


#ATPK = function(coarse_raster, fine_raster, variogram, npoints = 8, nsim = 0, beta = NA, nmax = 10,
#                frac = 1.0, noise_sd = 0.0){
ATPK = function(coarse_raster, fine_raster, variogram, npoints = 8, nsim = 0, beta = NA, nmax = 10,
                frac = 1.0, noise_sd = 0.0){
  
  ### Area-to-Point kriging for the downscaling of coarse raster data.
  ### Converts the coarse pixels into spatial polygons and uses the ATPK method of Kyriakydis (2004).
  ### The function also performs Block Sequential Gaussian Simulation (BSGS).
  ### Returns the downscaled map or a stack of the simulations.
  
  require(gstat)
  
  ### Arguments: 
  # coarse_raster = coarse resolution raster to be downscaled.
  # fine_raster = fine resolution raster covering the extent of the coarse raster. It will be used to define 
  #               the target resolution for the downscaling. The values in this raster are not used.
  # frac = fraction of coarse data to use when performing block sequential simulation (increases the variability)
  # npoints = number of oints used to approximate the blocks (coarse pixels). Can be 4,8 or 16.
  # For the other parameters consult the documentation of 'krige' function - help(krige) 
  
  # If simulations are required, consider only the specified fraction of areal data.
  if (nsim > 0){
    drop = 1.0-frac
    ndrops = floor(drop*length(coarse_raster))
    values(coarse_raster)[sample(1:length(coarse_raster),ndrops,replace = FALSE)] = NA
  }
  
  BlockMap = rasterToPolygons(coarse_raster, fun=NULL, n=npoints, na.rm=TRUE, digits=12, dissolve=FALSE)
  
#### PArte DAniele
  #coarse_coord<- coordinates(coarse_raster)
  #rast2<-coarse_raster
  #rast2@data@values<-cbind(rast2@data@values,coarse_coord)
  #colnames(rast2@data@values)<-c("a","x","y")
  #print(proj4string(poly_coarse))
  #print(str(coarse_raster))
  #proj4string(point_poly)<-proj4string(BlockMap)
  #print("COORDINATES")
  #print(str(coarse_coord))
  #val2 = values(coarse_raster)
  #names(val2)="a"
  #print("pre")
  #print(str(val2))
  #val2 <-cbind(val2,coarse_coord)
  #colnames(val2) = c("a","x","y")
  #print("post")
  #print(str(val2))
  #rast2@data@names = "t1"
  #print(str(rast2))
  #library(rgeos)
  #centroid<- gCentroid(BlockMap,byid=TRUE)
  #print("centroid")
  #print(str(centroid))
  #print("end centroid")

  #BlockMap2 = rasterToPolygons(rast2, fun=NULL, n=npoints, na.rm=TRUE, digits=12, dissolve=FALSE)
  #print("Names")
  names(BlockMap) = "Z"
  #names(point_poly) <-names(BlockMap)
  #centpoint<-gCentroid(point_poly,byid=TRUE)
  #print(str(centpoint@coords))

  #point_poly@data<-cbind(point_poly@data,x=centpoint@coords[,1])
  #point_poly@data<-cbind(point_poly@data,y=centpoint@coords[,2])
  #print("point_poly")
  #print(str(point_poly@data))
  #print(names(point_poly))
  #print(names(BlockMap))
  #print(proj4string(point_poly))
  #print(identicalCRS(poly_coarse,point_poly))




  #xcol = centroid@coords[,1]
  #print(str(xcol))
 ##BlockMap@data<-cbind(BlockMap@data, x=centroid@coords[,1])
 ##BlockMap@data<-cbind(BlockMap@data, y=centroid@coords[,2])
  #### The following lines add the polygons of field data ###########
  #BlockMap<- rbind(BlockMap,point_poly,makeUniqueIDs = TRUE)
  ###################################################################
  #print(str(BlockMap@data,max.level =2))
  
  # Add optional noise to the block data
  if(noise_sd){
    print("Adding an uniform noise to the block data")
    noise = rnorm(nrow(BlockMap),mean = 0, sd = noise_sd)
    BlockMap@data = BlockMap@data + noise
  }
  
  d_new <- SpatialPoints(coordinates(fine_raster)[!is.na(values(fine_raster)),], proj4string = crs(BlockMap))
  #print("d_new altro")
  #print(str(d_new))
  #print(str(BlockMap,max.level=3))
  #print(proj4string(BlockMap))
  #print(names(coordinates(BlockMap)))
  #print(attr(coordinates(BlockMap),"dimnames"))
  #print(str(coordinates(BlockMap)))
  #print(proj4string(BlockMap))
  #print(attributes(d_new))
  if(is.na(beta)) {
    downscale = krige(Z~1, BlockMap, newdata = d_new, model = variogram, nmax = nmax, nsim=nsim)
  } else {
    downscale = krige(Z~1, BlockMap, newdata = d_new, model = variogram, nmax = nmax, nsim=nsim, beta = beta)
  }
  if (nsim == 0) {
    downscaled_map = fine_raster
    values(downscaled_map)[!is.na(values(fine_raster))] = downscale$var1.pred
    return(downscaled_map)
  } else {
    sim_data = downscale@data 
    sims = stack(replicate(nsim,fine_raster))
    for (i in 1:nsim) values(sims[[i]])[!is.na(values(fine_raster))] = sim_data[,i]
    return(sims)
  }
}  

#################################

### Compositional data analysis (Aitchison, 1986).

library(compositions) # For compositional data analysis.

### ILR transformation (Egozcue et al., 2003).
#help(ilr)

# Apply closure operation to the psf dataset (acomp turns a df into a compositional object).
#print("psf_df")
#print(str(psf_df))
#print(str(psf_field))
psf_fieldbis<-data.frame(psf_field$fines*0.5*100,psf_field$sand*100,psf_field$fines*0.5*100)
psf_fieldbis = psf_fieldbis
#print("psf_fieldbis")
colnames(psf_fieldbis)=c("CLAY","SAND","SILT")
#print(psf_fieldbis)
#print(str(psf_fieldbis))
#print(str(psf_df))
#print(str(psf_fieldbis))
psf_comp = acomp(psf_df[,c(1,3,2)]) 
psf_fieldbiscomp =acomp(psf_fieldbis)
#print(class(psf_comp))
#print(class(psf_fieldbiscomp))

#print("psf_df")
#print(psf_df[1,c(1,2,3)])
#print(typeof(psf_df))

#print(psf_comp)
### Remark: silt and sand were switched to match the order of the thesis!
#### My part
#psf_field_comp = acomp(psf_field[,c(1,2,3)])

# Apply Isometric Log-Ratio to the data (with default basis), and put the results in a data frame.
ilr_df = data.frame( ilr(x = psf_comp) )
#print(ilr_df)
colnames(ilr_df) = c('ILR1','ILR2')

ilr_field = data.frame(ilr(x=psf_fieldbiscomp))
colnames(ilr_field) = c('ILR1','ILR2')
#print(str(ilr_df))
#print(str(ilr_field))
#point_raster1<-rasterize(coords,dem,ilr_field$ILR1,background=NA) ## ILR1 EQUIVALENT
#point_raster2<-rasterize(coords,dem,ilr_field$ILR2,background=NA) ## ILR2 EQUIVALENT
#print("point_raster!!!")
#print(str(point_raster1))

#print("ilr_df")
#print(ilr_df)

# When applying ilr, NAs are turned into 0s, set to NA the original NAs.  
ilr_df[is.na(psf_df$CLAY),] = NA 

# Create ilr rasters.
ILR1 = psf[[1]] # Copy clay map 
ILR2 = ILR1
values(ILR1) = ilr_df$ILR1 # Change values 
values(ILR2) = ilr_df$ILR2

#library(latex2exp) # For formulas using LaTeX syntax in plot titles.
#quartz(width=10,height = 5)
#par(mfrow=c(1,2))
#plot(ILR1, main = TeX("$\\bar{ilr}_1$")) # The overline indicates the coarse resolution.
#plot(ILR2, main = TeX("$\\bar{ilr}_2$"))

### ATPRK (Area-to-Point Regression Kriging).

# Regression step using Digital Elevation. 
# THIS STEP IS NOT NECESSARY. 
# Upscale covariate map using function 'resample' with method bilinear (upscales by averaging).
coarse_dem = resample(dem, ILR1, method="bilinear")

#print("coord position in raster")
#fpoints <- cellFromXY(coarse_dem,coords)
#print(fpoints)
#plot(coarse_dem)

# Add column 'DEM' to the data frame with the ilrs.
ilr_df['DEM'] = values(coarse_dem)

#library(psych) # For pretty scatterplots.
#quartz(width=6, height=6)
#pairs.panels(ilr_df[c('DEM','ILR1','ILR2')], 
#             method = "pearson", # correlation method
#             hist.col = "#00AFBB",
#             density = TRUE,  # show density plots
#             ellipses = FALSE, # show correlation ellipses
#             labels= c(TeX("$\\bar{ilr}_1$"),TeX("$\\bar{ilr}_2$"),TeX("$\\bar{DEM}$"))
#)

#ATPlm_field(ILR1,ilr_field$ILR1,coords,stack(c(dem,dem^2)) )
lm_dem1_field = ATPlm_field(ILR1,ilr_field$ILR1,coords,stack(c(dem,dem^2)) )
lm_dem2_field = ATPlm_field(ILR2,ilr_field$ILR2,coords,stack(c(dem,dem^2)) )
# Perform Area-to_Point Linear Regression using dem and dem^2 as covariates.
lm_dem1 = ATPlm(ILR1,stack(c(dem,dem^2)))
lm_dem2 = ATPlm(ILR2,stack(c(dem,dem^2)))
#print("map")
#print(lm_dem1$map)

#plot(lm_dem2$map)
#plot(values(ILR2)[!is.na(values(ILR2))], lm_dem2$LM$fitted.values, col = "blue", cex=.5,
#     xlab = "observed", ylab = "fitted", 
#     main = TeX("$\\mathbf{E}\\[ILR_2\\] = \\beta_0^{(2)}+ \\beta_1^{(2)} \\cdot DEM + \\beta_2^{(2)} \\cdot DEM^2$") )
#abline(a=0,b=1, col = "red", lwd = 2)

### Area to point kriging and simulation.

library(gstat) # Kriging, ATPK, variogram estimation, (Block) Sequential Gaussian Simulation.

# Compute residuals.
#print(length(lm_dem1$map@data@values))
#print(str(dem))
xseq<-seq(from=dem@extent@xmin,to=dem@extent@xmax,length.out=dem@ncols)
yseq<-seq(from=dem@extent@ymin,to=dem@extent@ymax,length.out=dem@nrows)
vseq1<-matrix(lm_dem1$map@data@values,nrow=dem@nrows,ncol = dem@ncols)
vseq2<-matrix(lm_dem2$map@data@values,nrow=dem@nrows,ncol = dem@ncols)
#print(length(xseq)*length(yseq))
#print(str(vseq1))
#print(str(psf_field))
#print(str(shape@coords))
x_field<-field_coords[,1]
y_field<-field_coords[,2]
#print(str(x_field))
#print(str(y_field))
#print(is.matrix(vseq1))
#print(typeof(vseq1))
#library(pracma)
#print(vseq1)
#print(dem@extent@xmin)
#print(x_field[1])
#print(dem@extent@xmax)
#print(dem@extent@ymin)
#print(y_field[1])
#print(dem@extent@ymax)
#zero_field1<- interp2(xseq,yseq,vseq1,x_field,y_field,method="nearest")
#zero_field2<- interp2(xseq,yseq,vseq2,x_field,y_field,method="nearest")
#print(zero_field1)
#print(zero_field2)
#print(lm_dem1)
field_index <-cellFromXY(lm_dem1_field$map,coords)
tmpval =values(lm_dem1_field$map)
zero_field1 <- tmpval[field_index]
field_index <-cellFromXY(lm_dem2_field$map,coords)
tmpval =values(lm_dem2_field$map)
zero_field2 <- tmpval[field_index]
#print("zero field!!!")
#print(zero_field1)
#print(zero_field2)

RES1 = ILR1 - resample(lm_dem1$map, ILR1, method="bilinear")
RES2 = ILR2 - resample(lm_dem2$map, ILR2, method="bilinear")
## Maybe add the following
#tmp1 = point_raster1 - resample(lm_dem1$map,point_raster1,method ="bilinear")
#tmp2 = point_raster2 - resample(lm_dem1$map,point_raster2,method ="bilinear")
#tmp1 <-as.data.frame(tmp1)
#tmp1 <-na.omit(tmp1)
#val1 <-tmp1$layer

RES1field = ilr_field$ILR1- zero_field1
RES2field = ilr_field$ILR2 - zero_field2
#print("RES1Field")
#print(RES1field)

tmptestdf1 = data.frame(x = coords[,1],y=coords[,2], Z= RES1field)
coordinates(tmptestdf1)=c('x','y')
#print(str(tmptestdf1))
tmptestdf2 = data.frame(x = coords[,1],y=coords[,2], Z= RES2field)
coordinates(tmptestdf2)=c('x','y')
#print(str(tmptestdf2))
#print("AFTER")
#print(val1)
#print(length(RES1field))
#library(rgeos)
#RES1FieldPoly = compile_polygons(RES1field,x_field,y_field)
#RES2FieldPoly = compile_polygons(RES2field,x_field,y_field)
#print("QUI")
#print(str(RES1FieldPoly,max.level =2))
### Here, must make a raster for the residuals

ndata = length(RES1field)
zerodata = replicate(ndata,0)
field_res1_rast<-rasterize(coords,RES1,zerodata,background=NA) ## ILR1 EQUIVALENT
field_res2_rast<-rasterize(coords,RES2,zerodata,background=NA) ## ILR2 EQUIVALENT
#print(str(field_res1_rast))
ind_res1 <- cellFromXY(field_res1_rast,coords)
ind_res2 <- cellFromXY(field_res2_rast,coords)
tab_res1 <- table(ind_res1)
tab_res2 <- table(ind_res2)
#print(as.data.frame(tab_res2))
fieldDF1 = as.data.frame(tab_res1)
nms = names(tab_res1)
lnttab = length(rownames(fieldDF1))
#print(nms)
#print('lnttab')
#print(lnttab)
univ_dataField1 = replicate(lnttab,0)
univ_dataField2 = replicate(lnttab,0)
univ_coord = matrix(0,lnttab,2)
#print(str(univ_dataField2))
#print(univ_coord)
#print(str(univ_coord))
#print(univ_coord[16,])
#print(aggregate(rep(1, nrow(d.col)), by=d.col, FUN=length))
for (i in 1:ndata){

  row = subset(fieldDF1,fieldDF1[,1]==ind_res2[i])
  orig_id = as.numeric(rownames(row))
  weight = row[2]
  ### Create the first row of a new coords dataframe
  #print('tobefound')
  #print(ind_res1[i])
  #print('row')
  #print(row)
  #print('index')
  #print(i)
  #print(orig_id)
  #a = as.numeric(univ_coord[orig_id,1])+ as.numeric(coords[i,1]/weight)
  #b = as.numeric(univ_coord[orig_id,2])+ as.numeric(coords[i,2]/weight)
  #print("AAAA")
  #print(as.numeric(a[1]))
  #print(as.numeric(b[1]))
  #print("BBBB")
  univ_coord[orig_id,1] = as.numeric(univ_coord[orig_id,1])+ as.numeric(coords[i,1]/weight)
  univ_coord[orig_id,2] = as.numeric(univ_coord[orig_id,2])+ as.numeric(coords[i,2]/weight)
  #univ_coord[orig_id,1]  <- univ_coord[orig_id,1]+ coords[i,1]/weight
  #univ_coord[orig_id,2]  <- univ_coord[orig_id,2]+ coords[i,2]/weight

  univ_dataField1[orig_id] = as.numeric(univ_dataField1[orig_id]) +as.numeric( RES1field[i]/weight)
  univ_dataField2[orig_id] = as.numeric(univ_dataField2[orig_id]) + as.numeric(RES2field[i]/weight)
}
#print("END CYCLE")
#print(ind_res1)
#print(fieldDF1)
#print(univ_dataField1)
#print(univ_dataField2)
#print(univ_coord)
#print(str(univ_dataField1))
#print(str(zerodata))
#print(str(univ_coord))
#print(str(coords))

field_res1_rast<-rasterize(univ_coord,RES1,univ_dataField1,background=NA) ## ILR1 EQUIVALENT
field_res2_rast<-rasterize(univ_coord,RES2,univ_dataField2,background=NA) ## ILR2 EQUIVALENT
#print(str(field_res1_rast))

# Compute empirical variograms (as a reference to choose the variogram model).
emp_vg1 = variogram(layer ~ 1, data = as(RES1, 'SpatialPointsDataFrame')) 
emp_vg2 = variogram(layer ~ 1, data = as(RES2, 'SpatialPointsDataFrame')) 

# v1 = fit.variogram(emp_vg1, vgm(model = "Sph"))
# plot(emp_vg1, v1, xlab=TeX("$h$"),ylab=TeX("$\\gamma(h)$")) # Spherical
# v2 = fit.variogram(emp_vg2, vgm(model = "Exp"))
# plot(emp_vg2, v2) # Exponential

# Perform variogram deconvolution. 
#vd1 = variogram_deconvolution(RES1, vg_type = "Sph", nblocks = 4, maxiter = 100) 
#vd2 = variogram_deconvolution(RES2, vg_type = "Exp", nblocks = 4, maxiter = 100) 
#vd1_field = variogram_deconvolution_field(field_res1_rast, vg_type = "Sph", nblocks = 40, maxiter = 100) 
#vd2_field = variogram_deconvolution_field(field_res2_rast, vg_type = "Exp", nblocks = 40, maxiter = 100) 
vd1 = variogram_deconvolution(RES1, vg_type = "Sph", nblocks = 4, maxiter = 100) 
vd2 = variogram_deconvolution(RES2, vg_type = "Exp", nblocks = 4, maxiter = 100) 
#print("vd1_field")
#print(str(vd1_field))


########################## Plot the variograms. ################################################
#quartz(width=8, height=6)
#plot(emp_vg1$dist, emp_vg1$gamma, ylim = c(0,0.012), col = 'blue', lwd=2,
#     main = TeX("Variogram deconvolution of $\\gamma_1$"), ylab = "semivariance", xlab = "distance")
#lines(seq(0,4500,45),variogramLine(v1,dist_vector=seq(0,4500,45))$gamma, type='l', col = 'green', lwd = 2)
#lines(seq(0,4500,45),variogramLine(vd1,dist_vector=seq(0,4500,45))$gamma, type='l', col = 'red', lwd = 2)
#grid(col = 'cornsilk2')
#legend(x = 3000,0.006, legend=c('Empirical variogram','Initial fit', 'Deconvoluted variogram'), 
#       col = c('blue','green','red'), lty = c(2,1,1), cex = .8)

#quartz(width=8, height=6)
#plot(emp_vg2$dist, emp_vg2$gamma, ylim = c(0,0.003), col = 'blue', lwd=2,
#     main = TeX("Variogram deconvolution of $\\gamma_2$"), ylab = "semivariance", xlab = "distance")
#lines(seq(0,4500,45),variogramLine(v2,dist_vector=seq(0,4500,45))$gamma, type='l', col = 'green', lwd = 2)
#lines(seq(0,4500,45),variogramLine(vd2,dist_vector=seq(0,4500,45))$gamma, type='l', col = 'red', lwd = 2)
#grid(col = 'cornsilk2')
#legend(x = 3000,0.001, legend=c('Empirical variogram','Initial fit', 'Deconvoluted variogram'), 
#       col = c('blue','green','red'), lty = c(2,1,1), cex = .8)
################################################################################################

### ATPK of the residuals.

reduce_factor = as.numeric(args[2])
dem2 = aggregate(dem, fact=reduce_factor) # Reduce the resolution
# Sequential simulation can take several minutes when simulating very high resolution maps, so the resolution
# is reduced by at least a 2 factor (this is done by fuction aggregate).

# Create fine resolution mask with null values over the lake.
mask = resample(psf[[1]], dem2, method='ngb')
#print("mask str!!")
#print(str(mask))
#values(mask)[!is.na(values(mask))] = 1
values(mask) = 1

### COMPUTE SLOPE
#coarse_slope <- terrain(dem2, opt="slope",unit="tangent",neighbors=8)
coarse_slope <-aggregate(slp, fact=reduce_factor)
#print(values(coarse_slope))
coarse_slope[is.na(coarse_slope[])]<- 2.0
measure_whereabouts<-rasterize(coords,coarse_slope,RES1field,background=NA)
out_pt = stack()
test_pt = stack(out_pt,measure_whereabouts)
aux_test_pt = test_pt[[1]]
#writeRaster( aux_test_pt, file=paste0('../Outputs/',0,'/test_pt_',0,'.tif'), overwrite=TRUE )
#print("SLOPE!!!")
#print(projection(coarse_dem,asText =TRUE))
blev = 20.0 #0.35
tlev = 30.0 #0.6
### Build values in a way that interpolates linearly between blev and tlev
coarse_slope@data@values <- (coarse_slope@data@values-blev)/(tlev-blev)

### Create the two masks
field_mask <-coarse_slope
field_mask[(field_mask>=1.0)]<-2.0
field_mask[(field_mask<1.0)]<-1.0
field_mask[(field_mask>=1.9)]<-0.0

sat_mask<-coarse_slope
sat_mask[sat_mask <= 0.0] <-0.0
sat_mask[sat_mask >0.0] <-1
#print(values(field_mask))
#print(str(coarse_slope))
#print(str(field_mask))
#print(table(is.na(field_mask@data@values)))
#print(table(is.na(sat_mask@data@values)))
intersec<-sat_mask*field_mask
#intersec@data@values = field_mask@data@values*sat_mask@data@values
#print(table((intersec@data@values)==1))
#print(table((intersec@data@values)==0))
trans_weight<- coarse_slope*intersec
#print(table((trans_weight@data@values)>1))
#print(table((trans_weight@data@values)==0))
### Now build exclusive masks
#mask_field <- mask(field_mask)
### Create the two masks
field_mask <-coarse_slope
field_mask[(field_mask>=0.0)]<-2.0
field_mask[(field_mask<0.0)]<-1.0
field_mask[(field_mask>=1.9)]<-0.0

sat_mask<-coarse_slope
sat_mask[sat_mask <= 1.0] <-0.0
sat_mask[sat_mask >1.0] <-1.0



#### Write down values of the data

#print("HEEEREEE!!!")
#mydata <- tmptestdf1
#write.csv(tmptestdf1,"data1.csv")
#write.csv(tmptestdf2,"data2.csv")
#xycoords <- mydata@coords
#xyid <- cellFromXY(Curv1,xycoords)
#tmpVals = values(Curv1) ### PF
#xyF1Vals <-tmpVals[xyid]
#write.csv(xyF1Vals,"curv_PF.csv")
#tmpVals = values(Curv2) ### PL
#xyF2Vals <-tmpVals[xyid] 
#write.csv(xyF2Vals,"curv_PL.csv")
#tmpVals = values(Curv3) ### TG
#xyF3Vals <-tmpVals[xyid]
#write.csv(xyF3Vals,"curv_TG.csv")
#tmpVals = values(Curv4) ### Gen 
#xyF4Vals <-tmpVals[xyid]
#write.csv(xyF4Vals,"curv_GEN.csv")
#tmpVals = values(slp)   ### SLP
#xySLVals <-tmpVals[xyid]
#write.csv(xySLVals,"slope.csv")

# Simulation.
if ( as.numeric(args[1]) == 0) {
  # Downscaling.
  #down1 = ATPK(RES1, mask, variogram = vd1, nmax = 8, beta = 0)
  #down2 = ATPK(RES2, mask, variogram = vd2, nmax = 8, beta = 0)
#  print(str(mask))
#  write.csv(tmptestdf1,"data_ilr1.csv")
#  write.csv(tmptestdf2,"data_ilr2.csv")
  down1_field   = ATPK_field(RES1, tmptestdf1,fine_raster = mask, nmax = 8, vg_type = "Sph")
  down2_field   = ATPK_field(RES2, tmptestdf2,fine_raster = mask, nmax = 8, vg_type = "Exp")
  #down1_field = ATPK(field_res1_rast, RES1FieldPoly, mask, variogram = vd1_field, nmax = 8)
  #down2_field = ATPK(field_res2_rast, RES2FieldPoly, mask, variogram = vd2_field, nmax = 8)
  down1       = ATPK(RES1, mask, variogram = vd1, nmax = 8)
  down2       = ATPK(RES2, mask, variogram = vd2, nmax = 8)

  ##### Weighting the rasters for satellite and field data using sat_mask, field_mask and trans_weight
  #print(str(testdown1))
  #print(str(down1_field))
  #down1 =  down1_field#*field_mask
  #down2 =  down2_field#*field_mask
  ### Write out the mask 
  test_mask1 = stack()
  test_mask1 = stack(test_mask1,sat_mask)
  test_mask2 = stack()
  test_mask2 = stack(test_mask2,trans_weight)
  test_mask3 = stack()
  test_mask3 = stack(test_mask3,field_mask)
  #test_mask4 = stack()
  #test_mask4 = stack(test_mask4,-field_mask+sat_mask+2*(trans_weight-0.5))
  

  aux_mask1 = test_mask1[[1]]
  aux_mask2 = test_mask2[[1]]
  aux_mask3 = test_mask3[[1]]
  #aux_mask4 = test_mask4[[1]]
  #down1 = down1_field
  #down2 = down2_field

  ## Commenting this as they are not used here
  #sim1 = ATPK(RES1, RES1FieldPoly, mask, variogram = vd1, nmax = 8, nsim = as.numeric(args[1]), frac = .5, noise_sd = 0.1)
  #sim2 = ATPK(RES2, RES2FieldPoly, mask, variogram = vd2, nmax = 8, nsim = as.numeric(args[1]), frac = .5, noise_sd = 0.1)
  
  # Add the baseline (regression results).
  baseline1 = aggregate(lm_dem1$map, fact = reduce_factor)
  baseline2 = aggregate(lm_dem2$map, fact = reduce_factor)
  down1 = down1 + baseline1
  down2 = down2 + baseline2

  # Add the baseline also to the field data results
  baseline1 = aggregate(lm_dem1_field$map, fact = reduce_factor)
  baseline2 = aggregate(lm_dem2_field$map, fact = reduce_factor)
  down1_field = down1_field +baseline1
  down2_field = down2_field +baseline2

  ## Combine the two masks
  down1 = down1*sat_mask + down1_field*field_mask + down1*trans_weight +down1_field*(intersec-trans_weight)
  down2 = down2*sat_mask + down2_field*field_mask + down2*trans_weight +down2_field*(intersec-trans_weight)
  
  
  # Back transform the results and create stack of rasters
  # This time particle size fractions are not stacked, but simulations of each are
  clay_sim_stack = stack()
  sand_sim_stack = stack()
  silt_sim_stack = stack()
  
  fine_ilr = as.data.frame(stack(c(down1,down2)))
  psf_sim_df = ilrInv(fine_ilr)
  # colnames(psf_sim_df) = c("CLAY", "SAND", "SILT")
  clay_sim_i = mask
  sand_sim_i = mask
  silt_sim_i = mask
  values(clay_sim_i) = psf_sim_df[,1]
  values(sand_sim_i) = psf_sim_df[,2]
  values(silt_sim_i) = psf_sim_df[,3]
  clay_sim_stack = stack(clay_sim_stack,clay_sim_i)
  sand_sim_stack = stack(sand_sim_stack,sand_sim_i)
  silt_sim_stack = stack(silt_sim_stack,silt_sim_i)
  
  
  
  # Save the solutions.
  aux_clay = clay_sim_stack[[1]]
  aux_sand = sand_sim_stack[[1]]
  aux_silt = silt_sim_stack[[1]]
  dir.create(paste0('../Outputs/',0))
  writeRaster( aux_clay, file=paste0('../Outputs/',0,'/clay_sim_',0,'.tif'), overwrite=TRUE )
  writeRaster( aux_sand, file=paste0('../Outputs/',0,'/sand_sim_',0,'.tif'), overwrite=TRUE )
  writeRaster( aux_silt, file=paste0('../Outputs/',0,'/silt_sim_',0,'.tif'), overwrite=TRUE )
  #writeRaster( aux_mask1, file=paste0('../Outputs/',0,'/mask1_sim_',0,'.tif'), overwrite=TRUE )
  #writeRaster( aux_mask2, file=paste0('../Outputs/',0,'/mask2_sim_',0,'.tif'), overwrite=TRUE )
  #writeRaster( aux_mask3, file=paste0('../Outputs/',0,'/mask3_sim_',0,'.tif'), overwrite=TRUE )
  
} else {
  #sim1 = ATPK(RES1, mask, variogram = vd1, nmax = 8, beta = 0, nsim = as.numeric(args[1]), frac = .5, noise_sd = 0.1)
  #sim2 = ATPK(RES2, mask, variogram = vd2, nmax = 8, beta = 0, nsim = as.numeric(args[1]), frac = .5, noise_sd = 0.1)

  #print("????")
  #sim1_field   = ATPK_field(RES1, tmptestdf1,fine_raster = mask, nmax = 8, vg_type = "Sph", nsim = as.numeric(args[1]))
  #sim2_field   = ATPK_field(RES2, tmptestdf2,fine_raster = mask, nmax = 8, vg_type = "Exp", nsim = as.numeric(args[1]))

  sim1 = ATPK(RES1, mask, variogram = vd1, nmax = 8, nsim = as.numeric(args[1]), frac = .5, noise_sd = 0.1)
  sim2 = ATPK(RES2, mask, variogram = vd2, nmax = 8, nsim = as.numeric(args[1]), frac = .5, noise_sd = 0.1)
  
  print("Limort")
  # Add the baseline (regression results).
  baseline1 = aggregate(lm_dem1$map, fact = reduce_factor)
  baseline2 = aggregate(lm_dem2$map, fact = reduce_factor)
  for (i in 1:length(sim1@layers)) {
    sim1[[i]] = sim1[[i]] + baseline1
    sim2[[i]] = sim2[[i]] + baseline2
  }

  # Add the baseline also to the field data results
  #baseline1 = aggregate(lm_dem1_field$map, fact = reduce_factor)
  #baseline2 = aggregate(lm_dem2_field$map, fact = reduce_factor)

  #for (i in 1:length(sim1@layers)) {
  #  sim1_field[[i]] = sim1_field[[i]] + baseline1
  #  sim2_field[[i]] = sim2_field[[i]] + baseline2
  #}

  ## Combine the two masks

  #for (i in 1:length(sim1@layers)) {
  #  sim1[[i]] = sim1[[i]]*sat_mask + sim1_field[[i]]*field_mask + sim1[[i]]*trans_weight + sim1_field[[i]]*(intersec-trans_weight)
  #  sim2[[i]] = sim2[[i]]*sat_mask + sim2_field[[i]]*field_mask + sim2[[i]]*trans_weight + sim2_field[[i]]*(intersec-trans_weight)
  #}
  
  # Back transform the results and create stack of rasters
  # This time particle size fractions are not stacked, but simulations of each are
  clay_sim_stack = stack()
  sand_sim_stack = stack()
  silt_sim_stack = stack()
  for (i in 1:length(sim1@layers)) {
    fine_ilr = as.data.frame(stack(c(sim1[[i]],sim2[[i]])))
    psf_sim_df = ilrInv(fine_ilr)
    # colnames(psf_sim_df) = c("CLAY", "SAND", "SILT")
    clay_sim_i = mask
    sand_sim_i = mask
    silt_sim_i = mask
    values(clay_sim_i) = psf_sim_df[,1]
    values(sand_sim_i) = psf_sim_df[,2]
    values(silt_sim_i) = psf_sim_df[,3]
    clay_sim_stack = stack(clay_sim_stack,clay_sim_i)
    sand_sim_stack = stack(sand_sim_stack,sand_sim_i)
    silt_sim_stack = stack(silt_sim_stack,silt_sim_i)
  }
  
  
  # Save the simulations.
  for (i in 1:length(sim1@layers)){
    aux_clay = clay_sim_stack[[i]]
    aux_sand = sand_sim_stack[[i]]
    aux_silt = silt_sim_stack[[i]]
    dir.create(paste0('../Outputs/',i))
    writeRaster( aux_clay, file=paste0('../Outputs/',i,'/clay_sim_',i,'.tif'), overwrite=TRUE )
    writeRaster( aux_sand, file=paste0('../Outputs/',i,'/sand_sim_',i,'.tif'), overwrite=TRUE )
    writeRaster( aux_silt, file=paste0('../Outputs/',i,'/silt_sim_',i,'.tif'), overwrite=TRUE )
  }
  
}




# Save the simulations.
#for (i in 1:length(sim1@layers)){
#  aux_clay = clay_sim_stack[[i]]
#  aux_sand = sand_sim_stack[[i]]
#  aux_silt = silt_sim_stack[[i]]
#  writeRaster( aux_clay, file=paste0('../Inputs/Geostatistics/clay_sim_',i,'.tif'), overwrite=TRUE )
#  writeRaster( aux_sand, file=paste0('../Inputs/Geostatistics/sand_sim_',i,'.tif'), overwrite=TRUE )
#  writeRaster( aux_silt, file=paste0('../Inputs/Geostatistics/silt_sim_',i,'.tif'), overwrite=TRUE )
#}





