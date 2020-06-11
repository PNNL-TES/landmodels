
#****************************************************************************************************
# http://geog.uoregon.edu/GeogR/topics/netCDF-in-R.html
#****************************************************************************************************
library(ncdf4)
library(raster)
library(tiff)
library(rgdal)

ncname <- "~/Data/Hashimoto/RH_yr_Hashimoto2015.nc"
# ncname <- "~/Data/Tang_Rh/RH.RF.720.360.1980.2016.Yearly.nc" 
ncname <- "~/Data/Wieder/casaclm_pool_flux_2000-2010_daily.nc"
dname <- "co2"  # note: tas means air temperature 

# open a netCDF file
ncin <- nc_open(ncname, readunlim=TRUE, auto_GMT=TRUE, suppress_dimvals=TRUE)
print(ncin)
summary(ncin)

x <- ncvar_get(ncin, "co2", start = c(1, 1, 1, 1), count = c(-1, -1, -1, 1))
image(x)

# Get the longtiudes and latitudes as before, using now the ncvar_get() function in ncdf4.

lon <- ncvar_get(ncin,"lon")
nlon <- dim(lon)
head(lon)


lat <- ncvar_get(ncin,"lat",verbose=F)
nlat <- dim(lat)
head(lat)

# time variable
time <- ncvar_get(ncin,"time")
time

head(time)
tail(time)


tunits <- ncatt_get(ncin,"time","units")
nt <- dim(time)
nt

tunits

print(c(nlon,nlat,nt))

# Get the input variable (tmp) and its attributes, and verify the size of the array.

tmp_array <- ncvar_get(ncin,dname)
dlname <- ncatt_get(ncin,dname,"long_name")
dunits <- ncatt_get(ncin,dname,"units")
fillvalue <- ncatt_get(ncin,dname,"_FillValue")
dim(tmp_array)


# Get the global attributes.

title <- ncatt_get(ncin,0,"title")
ncatt_get(ncin,0,"source")

# Close the netCDF file using the nc_close() function.
#nc_close(ncin)
#ls()


#****************************************************************************************************
# Reshape data and plot
#****************************************************************************************************
library(chron)
library(lattice)
library(RColorBrewer)

# split the time units string into fields
tustr <- strsplit(tunits$value, " ")
tdstr <- strsplit(unlist(tustr)[3], "-")
tmonth=as.integer(unlist(tdstr)[2])
tday=as.integer(unlist(tdstr)[3])
tyear=as.integer(unlist(tdstr)[1])
chron(time,origin=c(tmonth, tday, tyear))


#Replace netCDF fillvalues with R NAs
tmp_array[tmp_array==fillvalue$value] <- NA

length(na.omit(as.vector(tmp_array[,,1])))

# create map
m <- 1
tmp_slice <- tmp_array[,,m]
image(lon,lat,tmp_slice, col=rev(brewer.pal(10,"RdBu")))


# a better map
grid <- expand.grid(lon=lon, lat=lat)
cutpts <- seq(0,10,1)
levelplot(tmp_slice ~ lon * lat, data=grid, at=cutpts, cuts=11, pretty=T, 
          col.regions=(rev(brewer.pal(10,"RdBu"))))


#**************************************************************************
# 2. Get time information 1
#**************************************************************************
# matrix (nlon*nlon rows by 2 cols) of lons and lats
lonlat_time <- as.matrix(expand.grid(lon,lat,time))
dim(lonlat_time)

tmp_vec <- as.vector(tmp_slice)
length(tmp_vec)

#The data.frame() and cbind() functions are used to assemble the columns of the data frame, 
#which are assigned appropriate names using the names() function (on the left-hand side of 
#assignment operator). The head() function, applied on top of the na.omit() function lists 
#the first rows of values without NAs:

tmp_df01 <- data.frame(cbind(lonlat_time,tmp_vec))
names(tmp_df01) <- c("lon","lat","time",paste(dname,as.character(m), sep="_"))
head(na.omit(tmp_df01), 10)

#****************************************************************************************************
# lean raster package
#****************************************************************************************************
warner_rh <- raster("~/Data/Vargas_Warner_Rs/Rh_BondLamberty2004.tif")
# warner_rh_dat <- rasterToPoints(warner_rh)

# cellStats(warner_rh, min, na.rm = TRUE)
# view Coordinate Reference System (CRS)
warner_rh@crs
warner_rh@extent
ex <- raster::extent(c(129, 129.5, 62, 62.5))
raster::extract(warner_rh, ex, fun = mean, na.rm = TRUE)

# the distribution of values in the raster
hist(warner_rh, main="Distribution of elevation values", 
     col= "purple")

plot(warner_rh, 
     main="Rh from Warner")

image(warner_rh)
custom_bins <- seq(0, 1200, 200)

DSM_HARV_df <- DSM_HARV_df %>%
  mutate(fct_elevation_2 = cut(HARV_dsmCrop, breaks = custom_bins))

unique(DSM_HARV_df$fct_elevation_2)


ggplot() +
  geom_raster(data = warner_rh , aes(x = x, y = y)) + 
  coord_quickmap()

# plot tiff data
warner_tif <- raster("~/Data/Vargas_Warner_Rs/Rh_BondLamberty2004.tif")
plot(warner_tif)
