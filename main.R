# Visualize models' global HR

library(raster)
library(ncdf4)
library(tibble)
library(tidyr)
library(drake)
library(dplyr)
library(ggplot2)
theme_set(theme_bw())

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
    time <- ncvar_get(nc, "time")[1:365]
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


plan <- drake::drake_plan(
  lm_files = list.files("~/Data/Wieder/", pattern = "*.nc", full.names = TRUE),
  landmodel_dat = read_landmodels(lm_files),
  
  hashimoto_dat = read_hashimoto_rh("~/Data/Hashimoto/RH_yr_Hashimoto2015.nc"),
  
  hda = read_hashimoto_areas("~/Data/Hashimoto/RH_yr_Hashimoto2015.nc"),
  tang_dat = read_tang_rh("~/Data/Tang_Rh/RH.RF.720.360.1980.2016.Yearly.nc", hda),
  
  warner_dat = read_warner_rh("~/Data/Vargas_Warner_Rs/Rh_BondLamberty2004.tif"),
  
  # this plot moved to the markdown file
  # latplot = plot_latitude(bind_rows(landmodel_dat, tang_dat, hashimoto_dat, warner_dat)),
  
  # load the global monthly heterotrophic respiration data
  MGRhD = read_file('MGRhD.csv')
)
