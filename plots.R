# Plots

LAT_BANDS <- seq(-90, 90, by = 10)

plot_latitude <- function(dat, lat_bands = LAT_BANDS) {
  
  dat %>% 
    # for each model, mean value over time by latitude
    group_by(model, type, lat) %>% 
    summarise(hr_gC_m2_yr = mean(hr_gC_m2_yr, na.rm = TRUE)) %>% 
    # fill NAs
    replace_na(list(hr_gC_m2_yr = 0)) %>% 
    # put on common latitude bands
    complete(model, type, lat = c(lat_bands, .$lat)) %>% 
    arrange(model, type, lat) %>% 
    # interpolate
    group_by(model, type) %>% 
    mutate(hr_gC_m2_yr = approx(lat, hr_gC_m2_yr, xout = lat, rule = 2)$y) %>% 
    filter(lat %in% lat_bands) %>% 
    mutate(hr_gC_m2_yr_percent = hr_gC_m2_yr / sum(hr_gC_m2_yr, na.rm = TRUE)) ->
    plotting_data
  
  plotting_data %>% 
    ggplot(aes(lat, hr_gC_m2_yr_percent, color = model, linetype = type, size = type)) +
    geom_line() +
    scale_size_manual(values = c(1.25, 0.75)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
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
