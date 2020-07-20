## plot global pattern of soil heterotrophic respiration data
## 2020-07-19

# load required packages
library(ncdf4)
library(raster)
library(fields)
library(RColorBrewer)

setwd('E:/extdata')

# import data
nc_file <- nc_open('mimics_pool_flux_2000-2010_mean.nc')
shr_MIMICS <- ncvar_get(nc_file, "cHresp")
nc_close(nc_file)

nc_file <- nc_open('casaclm_pool_flux_2000-2010_mean.nc')
shr_CASACLM <- ncvar_get(nc_file, "cresp")
nc_close(nc_file)

nc_file <- nc_open('corpse_pool_flux_2000-2010_mean.nc')
shr_CORPSE <- ncvar_get(nc_file, "Soil_CO2")
nc_close(nc_file)


shr_SRDB <- raster('heterotrophic_resp_Bond_Lamberty_2004.tif')
#res(shr_SRDB)

# upscale SRDB data (area-averaged) to the resolution of MIMICS, CASACLM and CORPSE
lat_dim_source <- 21600
lon_dim_source <- 43200
lat_dim_target <- 96
lon_dim_target <- 144

shr_SRDB <- aggregate(shr_SRDB, fact=c(lon_dim_source/lon_dim_target, lat_dim_source/lat_dim_target),
                         fun=mean, na.rm = TRUE)

# change the lon from 0-360, the unit from y-1 to day-1 of SRDB data
shr_SRDB2 <- t(as.matrix(shr_SRDB))
shr_SRDB2 <- shr_SRDB2[,ncol(shr_SRDB2):1]
shr_SRDB2 <- shr_SRDB2[c((nrow(shr_SRDB2)/2+1):nrow(shr_SRDB2),1:(nrow(shr_SRDB2)/2)),]
shr_SRDB2 <- shr_SRDB2/365

# plot spatial pattern of SRH data
dev.new(width=500, height=500, unit="px",noRStudioGD = TRUE)
par(mfrow = c(2,2))
image.plot(shr_MIMICS, zlim = c(0, 5),
           main="MIMICS",
           col = brewer.pal(n = 8, name = "Oranges"))
image.plot(shr_CASACLM, zlim = c(0, 5),
           main="CASACLM",
           col = brewer.pal(n = 8, name = "Oranges"))
image.plot(shr_CORPSE, zlim = c(0, 5),
           main="CORPSE",
           col = brewer.pal(n = 8, name = "Oranges"))
image.plot(shr_SRDB2, zlim = c(0, 5),
           main="SRDB",
           col = brewer.pal(n = 8, name = "Oranges"))

# plot difference among different SRH data
dev.new(width=550, height=200, unit="px",noRStudioGD = TRUE)
par(mfrow = c(1,3))
image.plot(shr_MIMICS - shr_SRDB2, zlim = c(-2, 2),
           main="MIMICS-SRDB",
           col = rev(brewer.pal(n = 10, name = "RdBu")))
image.plot(shr_CASACLM - shr_SRDB2, zlim = c(-2, 2),
           main="CASACLM-SRDB",
           col = rev(brewer.pal(n = 10, name = "RdBu")))
image.plot(shr_CORPSE - shr_SRDB2, zlim = c(-2, 2),
           main="CORPSE-SRDB",
           col = rev(brewer.pal(n = 10, name = "RdBu")))
