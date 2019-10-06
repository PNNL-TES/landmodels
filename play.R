# Visualize models' global HR

library(ncdf4)
library(tibble)
library(tidyr)
library(drake)
library(dplyr)
library(ggplot2)
theme_set(theme_bw())

read_tang_rh <- function(f) {
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
      hr <- ncvar_get(nc, "variable", start = c(1, lat, t), count = c(-1, 1, 1)) # gC/m2/day
      dat$hr[lat + (t - 1) * length(lat)] <- sum(hr * 1e6, na.rm = TRUE) / 1e15 # PgC/day * cellMissing[,lat] * landarea[,lat]
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
        dat$hr[lat + (t - 1) * length(lat)] <- sum(hr * 1e6 * cellMissing[,lat] * landarea[,lat], na.rm = TRUE) / 1e15 # PgC/day
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
    # for each model, mean value by latitude
    group_by(model, type, lat) %>% 
    summarise(hr_PgC = mean(hr, na.rm = TRUE)) %>% 
    # put on common latitude bands
    complete(model, type, lat = c(standard_latitudes, .$lat)) %>% 
    arrange(model, type, lat) %>% 
    # interpolate
    group_by(model, type) %>% 
    mutate(hr_PgC = approx(lat, hr_PgC, xout = lat, rule = 2)$y) %>% 
    filter(lat %in% standard_latitudes) %>% 
    mutate(hr_PgC_percent = hr_PgC / sum(hr_PgC, na.rm = TRUE)) ->
    plotting_data
  
  plotting_data %>% 
    ggplot(aes(lat, hr_PgC_percent, color = model, linetype = type, size = type)) +
    geom_line() +
    scale_size_manual(values = c(1, 0.5)) +
    scale_y_continuous(labels = scales::percent) ->
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
  
  tang_dat = read_tang_rh("~/Data/Tang_Rh/RH.RF.720.360.1980.2016.Yearly.nc"),
  
  latplot = plot_latitude(bind_rows(landmodel_dat, tang_dat))
)

