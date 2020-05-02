
#****************************************************************************************************
# http://geog.uoregon.edu/GeogR/topics/netCDF-in-R.html
#****************************************************************************************************
library(ncdf4)

# ncname <- "~/Data/Hashimoto/RH_yr_Hashimoto2015.cn"  
ncname <- "~/Data/Tang_Rh/RH.RF.720.360.1980.2016.Yearly.nc" 
dname <- "cresp"  # note: tas means air temperature 

# open a netCDF file
ncin <- nc_open(ncname, readunlim=TRUE,auto_GMT=TRUE,suppress_dimvals=TRUE)
print(ncin)
summary(ncin)

# Get the longtiudes and latitudes as before, using now the ncvar_get() function in ncdf4.

lon <- ncvar_get(ncin,"longitude")
nlon <- dim(lon)
head(lon)


lat <- ncvar_get(ncin,"latitude",verbose=F)
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


