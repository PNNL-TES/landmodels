# Visualize models' global HR
library(tibble)
library(lubridate)
library(kableExtra)
library(ggplot2)
library(raster)
library(ncdf4)
library(tidyr)
library(dplyr)
library(drake)
theme_set(theme_bw())
library(magrittr)
library(data.table)
library(RColorBrewer)
library(fields)

DATA_DIR <- 'extdata'
OUT_DIR <- 'outputs'

source("plots.R")
source("functions.R")

read_warner_rh <- function(f) {
  dw <- raster(f)
  
  lat_bands <- seq(-90, 90, by = 5)
  lat_bands <- sort(lat_bands)  # make sure sorted for extent()
  
  dat <- tibble(model = "warner", 
                type = "benchmark",
                lat = (head(lat_bands, n = -1) + lat_bands[-1]) / 2,
                time = 2010,
                hr_gC_m2_yr = NA_real_, 
                hr_PgC_yr = NA_real_)
  
  for(i in head(seq_along(lat_bands), n = -1)) {
    message(basename(f), " ", lat_bands[i])
    ex <- raster::extent(c(-180, 180, lat_bands[i], lat_bands[i + 1]))
    dat$hr_gC_m2_yr[i] <- raster::extract(dw, ex, fun = mean, na.rm = TRUE)  # gC/m2/yr
    dat$hr_PgC_yr[i] <- NA
  }
  dat
}

read_hashimoto_areas <- function(f) {
  nc <- nc_open(f)
  on.exit(nc_close(nc))
  ncvar_get(nc, "area_cellarea")
}

read_hashimoto_rh <- function(f) {
  nc <- nc_open(f)
  time <- ncvar_get(nc, "time")
  latitude <- ncvar_get(nc, "lat")
  area_cellarea <- ncvar_get(nc, "area_cellarea")
  
  # Get each year in turn, sum across latitude
  dat <- expand.grid(model = "hashimoto", lat = -latitude, time = time,
                     hr_gC_m2_yr = NA_real_, hr_PgC_yr = NA_real_, stringsAsFactors = FALSE)
  for(t in tail(seq_along(time))) {
    message(basename(f), " ", time[t])
    for(lat in seq_along(latitude)) {
      hr <- ncvar_get(nc, "co2", start = c(1, lat, 1, t), count = c(-1, 1, 1, 1)) # gC/m2/yr
      i <- which(dat$lat == latitude[lat] & dat$time == time[t])
      stopifnot(length(i) == 1)
      dat$hr_gC_m2_yr[i] <- mean(hr, na.rm = TRUE)
      dat$hr_PgC_yr[i] <- sum(hr * 1e6 * area_cellarea[,lat], na.rm = TRUE) / 1e15  # PgC/yr
    }
  }
  nc_close(nc)
  as_tibble(dat) %>% mutate(time = 1901 + time, type = "benchmark")
}

read_tang_rh <- function(f, area_cellarea) {
  nc <- nc_open(f)
  time <- ncvar_get(nc, "z")
  latitude <- ncvar_get(nc, "latitude")
  
  # Get each year in turn, sum across latitude
  # Note it looks like Tang et al. flipped their hemispheres, which is why
  # there's a minus sign below
  dat <- expand.grid(model = "tang", lat = -latitude, time = time, hr = NA_real_, stringsAsFactors = FALSE)
  for(t in seq_along(time)) {
    message(basename(f), " ", time[t])
    for(lat in seq_along(latitude)) {
      hr <- ncvar_get(nc, "variable", start = c(1, lat, t), count = c(-1, 1, 1)) # gC/m2/yr?
      i <- which(dat$lat == latitude[lat] & dat$time == time[t])
      stopifnot(length(i) == 1)
      dat$hr_gC_m2_yr[i] <- mean(hr, na.rm = TRUE)
      dat$hr_PgC_yr[i] <- sum(hr * 1e6 * area_cellarea[,lat], na.rm = TRUE) / 1e15  # PgC/yr
    }
  }
  nc_close(nc)
  
  as_tibble(dat) %>% mutate(type = "benchmark", time = time + 1979)
}

read_landmodels <- function(lm_files) {
  
  # Open up each file in turn and extract land area, missing flag, and time
  resultslist <- list()
  for(f in lm_files) {
    nc <- nc_open(f)
    landarea <- ncvar_get(nc, "landarea")  # km2
    cellMissing <- !ncvar_get(nc, "cellMissing")
    # time <- ncvar_get(nc, "time")[1:365] # only 1 year
    time <- ncvar_get(nc, "time") # all data
    latitude <- ncvar_get(nc, "lat")
    
    modelname <- strsplit(basename(f), "_")[[1]][1]
    # Each model has a different output name for HR
    hr_var_name <- c("casaclm" = "cresp", "corpse" = "Soil_CO2", "mimics" = "cHresp") #g C m-2 day-1
    
    # Get each year in turn, sum across latitude
    dat <- expand.grid(model = modelname, lat = latitude, time = time, hr = NA_real_, stringsAsFactors = FALSE)
    for(t in seq_along(time)) {
      cat(basename(f), time[t], "\n")
      for(lat in seq_along(latitude)) {
        hr <- ncvar_get(nc, hr_var_name[modelname], start = c(1, lat, t), count = c(-1, 1, 1)) # gC/m2/day
        i <- which(dat$lat == latitude[lat] & dat$time == time[t])
        dat$hr_gC_m2_yr[i] <- mean(hr, na.rm = TRUE) * 365
        dat$hr_PgC_yr[i] <- sum(hr * 1e6 * cellMissing[,lat] * landarea[,lat], na.rm = TRUE) / 1e15
      }
    }
    nc_close(nc)
    
    resultslist[[f]] <- dat
  }
  bind_rows(resultslist) %>% mutate(type = "land model")
}

# get casa-cnp GPP and NPP data
read_casa_gppnpp <- function(lm_files) {
  # Open up each file in turn and extract land area, missing flag, and time
  resultslist <- list()
  for(f in lm_files) {
    nc <- nc_open(f)
    landarea <- ncvar_get(nc, "landarea")  # km2
    cellMissing <- !ncvar_get(nc, "cellMissing")
    # time <- ncvar_get(nc, "time")[1:365] # only 1 year
    time <- ncvar_get(nc, "time") # all data
    latitude <- ncvar_get(nc, "lat")
    
    modelname <- strsplit(basename(f), "_")[[1]][1]
    
    # Get each year in turn, sum across latitude
    dat <- expand.grid(model = modelname, lat = latitude, time = time, hr = NA_real_, stringsAsFactors = FALSE)
    for(t in seq_along(time)) {
      cat(basename(f), time[t], "\n")
      for(lat in seq_along(latitude)) {
        # hr
        hr <- ncvar_get(nc, "cresp", start = c(1, lat, t), count = c(-1, 1, 1)) # gC/m2/day
        i <- which(dat$lat == latitude[lat] & dat$time == time[t])
        dat$hr_gC_m2_yr[i] <- mean(hr, na.rm = TRUE) * 365
        dat$hr_PgC_yr[i] <- sum(hr * 1e6 * cellMissing[,lat] * landarea[,lat], na.rm = TRUE) / 1e15
        # gpp
        gpp <- ncvar_get(nc, "cgpp", start = c(1, lat, t), count = c(-1, 1, 1)) # gC/m2/day
        dat$gpp_gC_m2_yr[i] <- mean(gpp, na.rm = TRUE) * 365
        dat$gpp_PgC_yr[i] <- sum(gpp * 1e6 * cellMissing[,lat] * landarea[,lat], na.rm = TRUE) / 1e15
        # npp
        npp <- ncvar_get(nc, "cnpp", start = c(1, lat, t), count = c(-1, 1, 1)) # gC/m2/day
        dat$npp_gC_m2_yr[i] <- mean(npp, na.rm = TRUE) * 365
        dat$npp_PgC_yr[i] <- sum(npp * 1e6 * cellMissing[,lat] * landarea[,lat], na.rm = TRUE) / 1e15
      }
    }
    nc_close(nc)
    resultslist[[f]] <- dat
  }
  bind_rows(resultslist) %>% mutate(type = "land model")
}


# get rh from warner, hashimoto and tang for srdb-v5 data
get_srdb_rh <- function(dat) {
  dat %>% 
    dplyr::select(Latitude, Longitude, Study_midyear, Rh_annual, Manipulation) %>% 
    na.omit() %>% 
    mutate(Study_midyear = floor(Study_midyear)) %>% 
    arrange(Rh_annual) ->
    srdb_rh
  
  srdb_rh <- get_warner_rh(srdb_rh)
  srdb_rh <- get_tang(srdb_rh)
  srdb_rh <- get_hashimoto(srdb_rh)
}

# get rh from landmodels for mgrhd data
get_mgrhd_rh <- function() {
  # MGRhD_raw <- read.csv(file_in(paste0(DATA_DIR,'/MGRhD.csv')))
  MGRhD_raw <- read_file('MGRsD_matched_MODIS.csv') # need update to MGRsD_SRDB with NPP and GPP from MODIS
  MGRhD <- clean_mgrhd(MGRhD_raw)
  # MGRhD <- MGRhD [1:100,] # for testing purposeß
  mgrhd_landmodel <- get_landm("/Users/jian107/Data/Wieder//casaclm_pool_flux_1960-1969_daily.nc", MGRhD, 1960, 1969)
  mgrhd_landmodel <- get_landm("/Users/jian107/Data/Wieder//casaclm_pool_flux_1970-1979_daily.nc", mgrhd_landmodel, 1970, 1979)
  mgrhd_landmodel <- get_landm("/Users/jian107/Data/Wieder//casaclm_pool_flux_1980-1989_daily.nc", mgrhd_landmodel, 1980, 1989)
  mgrhd_landmodel <- get_landm("/Users/jian107/Data/Wieder//casaclm_pool_flux_1990-1999_daily.nc", mgrhd_landmodel, 1990, 1999)
  mgrhd_landmodel <- get_landm("/Users/jian107/Data/Wieder//casaclm_pool_flux_2000-2010_daily.nc", mgrhd_landmodel, 2000, 2010)

  mgrhd_landmodel <- get_landm("/Users/jian107/Data/Wieder//corpse_pool_flux_1960-1969_daily.nc", mgrhd_landmodel, 1960, 1969)
  mgrhd_landmodel <- get_landm("/Users/jian107/Data/Wieder//corpse_pool_flux_1970-1979_daily.nc", mgrhd_landmodel, 1970, 1979)
  mgrhd_landmodel <- get_landm("/Users/jian107/Data/Wieder//corpse_pool_flux_1980-1989_daily.nc", mgrhd_landmodel, 1980, 1989)
  mgrhd_landmodel <- get_landm("/Users/jian107/Data/Wieder//corpse_pool_flux_1990-1999_daily.nc", mgrhd_landmodel, 1990, 1999)
  mgrhd_landmodel <- get_landm("/Users/jian107/Data/Wieder//corpse_pool_flux_2000-2010_daily.nc", mgrhd_landmodel, 2000, 2010)

  mgrhd_landmodel <- get_landm("/Users/jian107/Data/Wieder//mimics_pool_flux_1960-1969_daily.nc", mgrhd_landmodel, 1960, 1969)
  mgrhd_landmodel <- get_landm("/Users/jian107/Data/Wieder//mimics_pool_flux_1970-1979_daily.nc", mgrhd_landmodel, 1970, 1979)
  mgrhd_landmodel <- get_landm("/Users/jian107/Data/Wieder//mimics_pool_flux_1980-1989_daily.nc", mgrhd_landmodel, 1980, 1989)
  mgrhd_landmodel <- get_landm("/Users/jian107/Data/Wieder//mimics_pool_flux_1990-1999_daily.nc", mgrhd_landmodel, 1990, 1999)
  mgrhd_landmodel <- get_landm("/Users/jian107/Data/Wieder//mimics_pool_flux_2000-2010_daily.nc", mgrhd_landmodel, 2000, 2010)
  
  # get gpp and npp from casaclm
  mgrhd_landmodel <- get_casa_gppnpp("/Users/jian107/Data/Wieder//casaclm_pool_flux_1960-1969_daily.nc", mgrhd_landmodel, 1960, 1969)
  mgrhd_landmodel <- get_casa_gppnpp("/Users/jian107/Data/Wieder//casaclm_pool_flux_1970-1979_daily.nc", mgrhd_landmodel, 1970, 1979)
  mgrhd_landmodel <- get_casa_gppnpp("/Users/jian107/Data/Wieder//casaclm_pool_flux_1980-1989_daily.nc", mgrhd_landmodel, 1980, 1989)
  mgrhd_landmodel <- get_casa_gppnpp("/Users/jian107/Data/Wieder//casaclm_pool_flux_1990-1999_daily.nc", mgrhd_landmodel, 1990, 1999)
  mgrhd_landmodel <- get_casa_gppnpp("/Users/jian107/Data/Wieder//casaclm_pool_flux_2000-2010_daily.nc", mgrhd_landmodel, 2000, 2010)
}

plan <- drake::drake_plan(
  lm_files = list.files("~/Data/Wieder/", pattern = "*.nc", full.names = TRUE), # need update nc data
  # landmodel_dat = read_landmodels(lm_files), # need update nc data, this takes too long, now do casa, corpse, and mimics seperately
  
  lm_files_casa = lm_files[3:5],
  landmodel_dat_casa = read_landmodels(lm_files_casa),
  casa_gpp_npp_rh = read_casa_gppnpp(lm_files_casa), # get casa gpp and npp
  modis_gpp_npp = read_file("lat_mean.csv"),

  lm_files_corpse = lm_files[8:10],
  landmodel_dat_corpse = read_landmodels(lm_files_corpse),

  lm_files_mimics = lm_files[13:15],
  landmodel_dat_mimics = read_landmodels(lm_files_mimics),
  
  hashimoto_dat = read_hashimoto_rh("~/Data/Hashimoto/RH_yr_Hashimoto2015.nc"),
  
  hda = read_hashimoto_areas("~/Data/Hashimoto/RH_yr_Hashimoto2015.nc"),
  tang_dat = read_tang_rh("~/Data/Tang_Rh/RH.RF.720.360.1980.2016.Yearly.nc", hda),
  
  warner_dat = read_warner_rh("~/Data/Vargas_Warner_Rs/Rh_BondLamberty2004.tif"),

  # this plot moved to the markdown file
  # latplot = plot_latitude(bind_rows(landmodel_dat, tang_dat, hashimoto_dat, warner_dat)),
  
  # load the global monthly heterotrophic respiration data
  MGRhD_landmodel = get_mgrhd_rh(),
  srdb = read.csv("~/Documents/PNNL/srdb/srdb-data.csv"),
  srdb_warner = get_warner_rs(srdb),
  srdb_rh = get_srdb_rh(srdb),
  
  # load CMIP6 data
  # Raw data and the code of processing the raw data for CMIP6 can be found at: GPP_NPP_Rh_LatLon.Rmd
  CMIP6_perc = read.csv(file_in(!!file.path("CMIP6",'rh_latlon_percent_global_1980-2015.csv'))),
  CMIP6_perc_GPP = read.csv(file_in(!!file.path("CMIP6",'gpp_latlon_percent_global_1980-2015.csv'))),
  CMIP6_perc_NPP = read.csv(file_in(!!file.path("CMIP6",'npp_latlon_percent_global_1980-2015.csv'))),
  lat_area = read.csv(file_in(!!file.path("CMIP6",'area_latlon_model.csv'))),
  cmip6_global = read.csv(file_in(!!file.path("CMIP6",'rh_global_1980-2015.csv'))),
  # cmip6_raw = clean_cmip6_raw(),
  # cmip6_annual_mean = get_annual_mean(cmip6_raw),
  # cmip6_global_mean = get_global_mean(cmip6_annual_mean),
  
  # spatial pattern map for CASA-CNP, corpse, mimics, and warner_2020 global Rh map
  nc_mimics = nc_open('extdata/mimics_pool_flux_2000-2010_mean.nc'),
  shr_MIMICS = ncvar_get(nc_mimics, "cHresp"),
  
  nc_casa = nc_open('extdata/casaclm_pool_flux_2000-2010_mean.nc'),
  shr_CASACLM = ncvar_get(nc_casa, "cresp"),

  nc_corpse = nc_open('extdata/corpse_pool_flux_2000-2010_mean.nc'),
  shr_CORPSE = ncvar_get(nc_corpse, "Soil_CO2"),
  
  shr_warner2020 = prepare_warner_rh()
)

make(plan)
