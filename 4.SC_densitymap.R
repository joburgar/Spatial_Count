######################################################################
# 4.SC_densitymap.R
# SCR book Ch 18 code;
# created by Joanna Burgar, 31-May-2017,
# modified for Koala sites, 23-April-2020
# modified using Joshua Goldberg's code for density map 10-June-2020
# updated by Joanna Burgar, 04-Aug-2021 for Leroy Gonzales (streamlined)
# script for creating realised density spatial maps
# using output from JAGS spatial count models 
######################################################################
###--- Functions
### SCRdensity function
###--- spatial density object function

spatial_density_obj.function <- function(out = out){
   # summary(window(out[,c("D","lam0","sigma","psi")], start = 41001))
   # gelman.diag(window(out[,c("D","lam0","sigma","psi")], start = 1001), multivariate = F)
   # plot(window(out[,c("D","lam0","sigma","psi")], start = 1001))
   
   out1 <-out[[1]] # separate into individual chains - each 50000 iterations
   out2 <-out[[2]]
   out3 <-out[[3]]
   
   Sxout1 <- out1[40001:50000,5:504]
   Syout1 <- out1[40001:50000,505:1004]
   z1 <- out1[40001:50000, 1006:1505]
   
   Sxout2 <- out2[40001:50000,5:504]
   Syout2 <- out2[40001:50000,505:1004]
   z2 <- out2[40001:50000, 1006:1505]
   
   Sxout3 <- out3[40001:50000,5:504]
   Syout3 <- out3[40001:50000,505:1004]
   z3 <- out3[40001:50000, 1006:1505]
   
   obj1 <- list(Sx=rbind(Sxout1, Sxout2, Sxout3), Sy=rbind(Syout1,Syout2,Syout3), z=rbind(z1,z2,z3))
   
   return(obj1)
}


SCRdensity <- function (obj, nx = 30, ny = 30, Xl = NULL, Xu = NULL, 
                        Yl = NULL, Yu = NULL, scalein = 100, scaleout = 100, 
                        col = "gray", ncolors = 10, whichguy = NULL){
   #obj <- SC_obj 
 
   Sxout <- obj$Sx
   Syout <- obj$Sy
   z <- obj$z
   niter <- nrow(z)
   if (is.null(Xl)) {
      Xl <- min(Sxout) * 0.999
      Xu <- max(Sxout) * 1.001
      Yl <- min(Syout) * 0.999
      Yu <- max(Syout) * 1.001
   }
   xg <- seq(Xl, Xu, , nx)
   yg <- seq(Yl, Yu, , ny)
   guy <- col(Sxout)
   Sxout <- cut(Sxout[z == 1], breaks = xg)
   Syout <- cut(Syout[z == 1], breaks = yg)
   if (is.null(whichguy)) {
      Dn <- table(Sxout, Syout)/niter
      area <- (yg[2] - yg[1]) * (xg[2] - xg[1]) * scalein
      Dn <- (Dn/area) * scaleout
   }
   else {
      Dn <- table(Sxout[guy == whichguy], Syout[guy == whichguy])/niter
   }
   cat("mean: ", mean(Dn), fill = TRUE)
   par(mar = c(3, 3, 3, 6))
   if (col == "gray") {
      cc <- seq(3, 17, , 10)/20
      cc <- rev(gray(cc))
   }
   else cc <- terrain.colors(ncolors)
   
   image(xg, yg, Dn)
   #image.scale(Dn, col = cc)
   box()
   
   return(list(grid = cbind(xg, yg), Dn = Dn))
}


###--- density map to raster function
# Dn object from SCRdensity function
# traplocs can be either matrix or dataframe of x and y locations (utm)
# buffer in m units

SCRraster <- function (Dn = Dn, traplocs = traplocs, buffer = buffer, crs = crs){
   r <- raster(ncol=ncol(Dn), nrow=nrow(Dn))
   values(r) <- Dn
   as.raster(r)
   t_r <- t(r) # for some reason map is flipped, need to transpose
   
   # add correct extent and coordinate system
   minx <- min(traplocs[,1])
   maxx <- max(traplocs[,1])
   miny <- min(traplocs[,2])
   maxy <- max(traplocs[,2])
   
   bb <- extent(minx-buffer, maxx+buffer,miny-buffer, maxy+buffer) 
   extent(r) <- bb
   r <- setExtent(r, bb, keepres=TRUE)
   crs(r) <- crs
   
   bb <- extent(minx-buffer, maxx+buffer,miny-buffer, maxy+buffer) 
   extent(t_r) <- bb
   t_r <- setExtent(t_r, bb, keepres=TRUE)
   crs(t_r) <- crs

   return(list(Dn_Traster = t_r, Dn_raster = r))
}

######################################
#Load Packages
list.of.packages <- c("raster","tidyverse", "Cairo", "coda")

# Check you have them and load them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages, require, character.only = TRUE)


###--- Load data
info <- list.files(recursive = T, pattern="*.RData") # pulls in all mcmc files, from subfolders
info_map <- info[grepl("map", info)] # retains the density map (large) file
for(i in 1:length(info_map)) load(info_map[[i]]) 

info_dat <- info[grepl("dat", info)] # retains the input dat file
for(i in 1:length(info_dat)) load(info_dat[[i]]) 

###---
# all MCMC output is in 100 m units
# if input coordinate system is in "100 m" then input scalein=1
# so have as 1 if want in 100 m units

###--- 
info_temp <- info[grepl("map", info)]
out_temp <- substr(info_temp,5,str_locate(info_temp, "_map.RData")[1]-1)
# out_model <- as.character(str_replace(out_temp, pattern = "2020_", replacement = "_"))

SC_obj <- spatial_density_obj.function(out = Kalateenee_KSI1)
str(SC_obj)

traplocs <- dat$sitedf[c("utm_x", "utm_y")]
colnames (traplocs) <- c("x","y")
traplocs
buffer <- 750 # 750 m buffer
sa_x <- max(traplocs$x+buffer) - min(traplocs$x-buffer) 
sa_y <- max(traplocs$y+buffer) - min(traplocs$y-buffer) 

SC_map <- SCRdensity(SC_obj, nx=round(sa_x/100)+1, ny=round(sa_y/100)+1)
str(SC_map)
SC_raster <- SCRraster(Dn = SC_map$Dn, traplocs = traplocs, buffer = 750, crs = c("+init=epsg:28356"))
plot(flip(SC_raster$Dn_Traster, direction="y"))
writeRaster(flip(SC_raster$Dn_Traster, direction="y"), paste("out/",out_temp,"_raster_diry.tif",sep=""), overwrite=TRUE)

Cairo(file=paste("out/",out_temp,"_SpatialDensity.PNG",sep=""), 
      type="png",
      width=2000, 
      height=1600, 
      pointsize=12,
      bg="white",
      dpi=300)
plot(flip(SC_raster$Dn_Traster, direction="y"))
points(traplocs, pch=20, cex=1)
mtext(paste(out_temp," \nRealised Density Map",sep=""), side = 3, line = 1, cex=1.25)
dev.off()