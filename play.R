# Visualize models' global HR

library(ncdf4)
library(tibble)

files <- list.files("~/Data/Wieder/", pattern = "*.nc", full.names = TRUE)

# Each model has a different output name for HR
hr_var_name <- c("casaclm" = "cresp", "corpse" = "Soil_CO2", "mimics" = "cHresp")

# Open up each file in turn and extract land area, missing flag, and time
resultslist <- list()
for(f in files) {
  nc <- nc_open(f)
  landarea <- ncvar_get(nc, "landarea")  # km2
  cellMissing <- !ncvar_get(nc, "cellMissing")
  time <- ncvar_get(nc, "time")
  
  modelname <- strsplit(basename(f), "_")[[1]][1]

  # Get each year in turn, sum across globe
  dat <- tibble(model = modelname, time = time, hr = NA_real_)
  for(t in seq_along(time)) {
    cat(basename(f), time[t], "\n")
    hr <- ncvar_get(nc, hr_var_name[modelname], start = c(1, 1, t), count = c(-1, -1, 1)) # gC/m2/day
    dat$hr[t] <- sum(hr * 1e6 * cellMissing * landarea, na.rm = TRUE) / 1e15 # PgC/day
  }
  nc_close(nc)
  
  resultslist[[f]] <- dat
}

# Combine results...
library(dplyr)
results <- bind_rows(resultslist)

# ...and plot
library(ggplot2)
theme_set(theme_bw())

p1 <- ggplot(results, aes(time, hr, color = model)) + geom_line()
print(p1)

results %>% 
  mutate(Year = floor(time)) %>% 
  group_by(model, Year) %>% 
  summarise(n = n(), hr_PgC = sum(hr)) %>% 
  filter(n == 365) %>% 
  ggplot(aes(Year, hr_PgC, color = model)) + geom_line() ->
  p2
print(p2)

