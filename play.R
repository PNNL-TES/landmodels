# Visualize models' global HR

library(ncdf4)
library(tibble)
library(tidyr)
library(drake)
library(dplyr)
library(ggplot2)
theme_set(theme_bw())

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
                     hr_gC_m2_yr = NA_real_, hr_PgC_yr = NA_real_,stringsAsFactors = FALSE)
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
  browser()
  nc <- nc_open(f)
  time <- ncvar_get(nc, "z")
  latitude <- ncvar_get(nc, "latitude")
  
  # Get each year in turn, sum across latitude
  # Note it looks like Tang et al. flipped their hemispheres, which is why
  # there's a minus sign below
  dat <- expand.grid(model = "tang", lat = -latitude, time = time, hr = NA_real_, stringsAsFactors = FALSE)
  for(t in seq_along(time)) {
    message(basename(f), " ", time[t], " ", lat)
    for(lat in seq_along(latitude)) {
      hr <- ncvar_get(nc, "variable", start = c(1, lat, t), count = c(-1, 1, 1)) # gC/m2/yr?
      i <- which(dat$lat == latitude[lat] & dat$time == time[t])
      stopifnot(length(i) == 1)
      dat$hr_gC_m2_yr[i] <- mean(hr, na.rm = TRUE)
      dat$hr_PgC_yr[i] <- sum(hr * 1e6 * area_cellarea[,lat], na.rm = TRUE) / 1e15  # PgC/yr
    }
  }
  nc_close(nc)
  
  as_tibble(dat) %>% mutate(type = "benchmark")
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

plot_latitude <- function(dat) {
  
  standard_latitudes <- seq(-90, 90, by = 10)
  
  dat %>% 
    # for each model, mean value over time by latitude
    group_by(model, type, lat) %>% 
    summarise(hr_gC_m2_yr = mean(hr_gC_m2_yr, na.rm = TRUE)) %>% 
    # fill NAs
    replace_na(list(hr_gC_m2_yr = 0)) %>% 
    # put on common latitude bands
    complete(model, type, lat = c(standard_latitudes, .$lat)) %>% 
    arrange(model, type, lat) %>% 
    # interpolate
    group_by(model, type) %>% 
    mutate(hr_gC_m2_yr = approx(lat, hr_gC_m2_yr, xout = lat, rule = 2)$y) %>% 
    filter(lat %in% standard_latitudes) %>% 
    mutate(hr_gC_m2_yr_percent = hr_gC_m2_yr / sum(hr_gC_m2_yr, na.rm = TRUE)) ->
    plotting_data
  
  plotting_data %>% 
    ggplot(aes(lat, hr_gC_m2_yr_percent, color = model, linetype = type, size = type)) +
    geom_line() +
    scale_size_manual(values = c(1.25, 0.75)) +
    scale_y_continuous(labels = scales::percent) +
    xlab("Latitude") +
    ylab(expression(R[H]~by~latitude)) ->
    p3
  print(p3)
  ggsave("outputs/p3.pdf", width = 6, height = 4)
}

plot_landmodels_time <- function(lm_dat) {
  p1 <- ggplot(lm_dat, aes(time, hr, color = model)) + geom_line()
  print(p1)
  ggsave("outputs/p1.pdf")
  
  lm_dat %>% 
    mutate(Year = floor(time)) %>% 
    group_by(model, Year) %>% 
    summarise(n = n(), hr_PgC = sum(hr, na.rm = TRUE)) %>% 
    filter(n == 365) %>% 
    ggplot(aes(Year, hr_PgC, color = model)) + geom_line() ->
    p2
  print(p2)
  ggsave("outputs/p2.pdf")
}


plan <- drake::drake_plan(
  lm_files = list.files("~/Data/Wieder/", pattern = "*.nc", full.names = TRUE),
  landmodel_dat = read_landmodels(lm_files),
  
  hashimoto_dat = read_hashimoto_rh("~/Data/Hashimoto/RH_yr_Hashimoto2015.nc"),
  
  hda = read_hashimoto_areas("~/Data/Hashimoto/RH_yr_Hashimoto2015.nc"),
  tang_dat = read_tang_rh("~/Data/Tang_Rh/RH.RF.720.360.1980.2016.Yearly.nc", hda),
  
  latplot = plot_latitude(bind_rows(landmodel_dat, tang_dat, hashimoto_dat))
)

