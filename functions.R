
#*****************************************************************************************************************
# prepare data 
#*****************************************************************************************************************
# read file
read_file <- function(x) read.csv(file.path(DATA_DIR, x), comment.char = "#", stringsAsFactors = FALSE)

# clean MGRhD for get Rh from mimics, casa, and corpes models
clean_mgrhd <- function(dat){
  dat %>%
    mutate(Meas_Month2 = ifelse(!is.na(Meas_Month2), Meas_Month2, Meas_Month)) %>%
    mutate(Meas_DOY = ifelse(!is.na(Meas_DOY), Meas_DOY, Meas_Month * 30.5-15)) %>%
    mutate(Meas_DOY = floor(Meas_DOY)) %>%
    filter(!is.na(Meas_DOY)) %>%
    filter(!is.na(Latitude) & !is.na(Longitude)) ->
    MGRhD
}

# function to get Rh data from wieder according to latitude, longitude, year, and day
get_warner_rh <- function(rh_dat){
  warner_rh <- raster("~/Data/Vargas_Warner_Rs/Rh_BondLamberty2004.tif")
  for (i in 1:nrow(rh_dat)){
    lat <- rh_dat$Latitude[i]
    lon <- rh_dat$Longitude[i]
    lat_low <- lat - 0.1
    lat_high <- lat + 0.1
    lon_low <- lon - 0.1
    lon_high <- lon + 0.1
    
    ex_cell <- raster::extent(c(lon_low, lon_high, lat_low, lat_high))
    rh_cell <- raster::extract(warner_rh, ex_cell, fun = mean, na.rm = TRUE)
    
    rh_dat[i, "warner_rh"] <- rh_cell
    print(paste0("*****", i))
  }
  return (rh_dat)
}

get_tang <- function(rh_dat) {
  tang_rh <- nc_open("~/Data/Tang_Rh/RH.RF.720.360.1980.2016.Yearly.nc")
  for (i in 1:nrow(rh_dat)) {
    target_lat <- rh_dat$Latitude[i]
    target_lon <- rh_dat$Longitude[i]
    target_year <- rh_dat$Study_midyear[i]
    
    # find the locate of data in tang_rh for ith Rh_annual
    lat <- ncvar_get(tang_rh, "latitude")
    ilat <- which.min(abs(lat - target_lat))
    
    lon <- ncvar_get(tang_rh, "longitude")
    ilon <- which.min(abs(lon - target_lon))
    
    time <- ncvar_get(tang_rh, "z")
    time2 <- time + 1979
    min_time <- min(abs(time2 - target_year))
    itime <- which.min(abs(time2 - target_year))
    
    Rh_cell <- ncvar_get(tang_rh, "variable", start = c(ilon, ilat, itime), count = c(1, 1, 1))
    Rh_avg <- mean(ncvar_get(tang_rh, "variable", start = c(ilon, ilat, 1), count = c(1, 1, -1))) 
    
    # using Rh_avg if Rh_cell is na (because Rh_tang only have data for 1980-2017)
    rh_dat[i, "tang_rh"] <- ifelse(min_time > 1, Rh_avg, Rh_cell)
    print(paste0("*****", i))
  }
  return (rh_dat)
}

get_hashimoto <- function(rh_dat) {
  hashi_rh <- nc_open("~/Data/Hashimoto/RH_yr_Hashimoto2015.nc")
  for (i in 1:nrow(rh_dat)) {
    target_lat <- rh_dat$Latitude[i]
    target_lon <- rh_dat$Longitude[i]
    target_year <- rh_dat$Study_midyear[i]
    
    # find the locate of data in tang_rh for ith Rh_annual
    lat <- ncvar_get(hashi_rh, "lat")
    ilat <- which.min(abs(lat - target_lat))
    
    lon_orig <- ncvar_get(hashi_rh, "lon") 
    lon <- ifelse(lon_orig <= 180, lon_orig, lon_orig - 360)
    ilon <- which.min(abs(lon - target_lon))
    
    time <- ncvar_get(hashi_rh, "time")
    time2 <- time + 1901
    itime <- which.min(abs(time2 - target_year))
    
    Rh_cell <- ncvar_get(hashi_rh, "co2", start = c(ilon, ilat, 1, itime), count = c(1, 1, 1, 1))
    Rh_avg <- mean(ncvar_get(hashi_rh, "co2", start = c(ilon, ilat, 1, 1), count = c(1, 1, 1, -1)), na.rm=T) 
    
    # using Rh_avg if Rh_cell is na (because Rh_tang only have data for 1980-2017)
    rh_dat[i, "hashi_rh"] <- ifelse(is.na(Rh_cell), Rh_avg, Rh_cell)
    print(paste0("*****", i))
  }
  return (rh_dat)
}


# get Rh from land models
get_landm <- function(f, mgrhd){
  landm <- nc_open(f)
  var_rh <- c("casaclm" = "cresp", "corpse" = "Soil_CO2", "mimics" = "cHresp")
  modelname <- strsplit(basename(f), "_")[[1]][1]
  
  for (i in 1:nrow(mgrhd)) {
    target_lat <- mgrhd$Latitude[i]
    target_lon <- mgrhd$Longitude[i]
    target_year <- mgrhd$Meas_Year[i]
    target_doy <- mgrhd$Meas_DOY[i]
    target_time <- target_year + target_doy/365
    
    # find the locate of data in tang_rh for ith Rh_annual
    lat <- ncvar_get(landm, "lat")
    ilat <- which.min(abs(lat - target_lat))
    
    lon_orig <- ncvar_get(landm, "lon") 
    lon <- ifelse(lon_orig <= 180, lon_orig, lon_orig - 360)
    ilon <- which.min(abs(lon - target_lon))
    
    time <- as.numeric(ncvar_get(landm, "time"))
    itime <- which.min(abs(time - target_time))
    min_time <- min(abs(time - target_time))
    doy <- time %% 1 * 364.25 + 1
    idoy <- which(abs(doy - target_doy) < 10)
    
    Rh_cell <- ncvar_get(landm, var_rh[modelname], start = c(ilon, ilat, itime), count = c(1, 1, 1))
    Rh_avg <- mean(ncvar_get(landm, var_rh[modelname], start = c(ilon, ilat, 1), count = c(1, 1, -1))[idoy], na.rm=T)
    
    mgrhd[i, modelname] <- ifelse(min_time > 1, Rh_avg, Rh_cell)
    
    print(paste0("*****", i))
  }
  return(mgrhd)
}
