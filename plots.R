# Plots

LAT_BANDS <- seq(-90, 90, by = 0.5)

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
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.01)) +
    xlab("Latitude") +
    ylab(expression(R[H]~by~latitude~("% of global total"))) ->
    p3
  print(p3)
  ggsave("outputs/p3.pdf", width = 6, height = 4)
}

plot_landmodels_time <- function(lm_dat) {
  # aggregate by model and time
  lm_dat %>% 
    group_by(model, time) %>% 
    summarise(hr_gC_m2_yr = mean(hr_gC_m2_yr, na.rm = TRUE)) %>% 
    # fill NAs
    replace_na(list(hr_gC_m2_yr = 0)) ->
    sdata
    
  p1 <- ggplot(sdata, aes(time, hr_gC_m2_yr, color = model)) + 
    geom_line() +
    labs (y = expression(R[H]~(g~C~m^{-2}~yr^{-1})))
  print(p1)
  # ggsave("outputs/p1.pdf")
}

global_rh_compar <- function(land_dat, t_dat, h_dat) {
  bind_rows(
    land_dat %>% 
      mutate(Year = floor(time)) %>% 
      group_by(model, Year) %>% 
      summarise(n = n(), hr_PgC = sum(hr_PgC_yr, na.rm = TRUE)) %>% 
      filter(Year == 2000) %>% 
      group_by(model) %>% 
      summarise(hr_PgC = mean(hr_PgC)),
    
    t_dat %>% 
      select(-hr) %>% 
      mutate(Year = floor(time)) %>% 
      na.omit() %>%
      group_by(model, Year) %>% 
      summarise(n = n(), hr_PgC = sum(hr_PgC_yr, na.rm = TRUE)) %>% 
      group_by(model) %>% 
      summarise(hr_PgC = mean(hr_PgC)),
    
    h_dat %>% 
      mutate(Year = floor(time)) %>% 
      na.omit() %>% 
      group_by(model, Year) %>% 
      summarise(n = n(), hr_PgC = sum(hr_PgC_yr, na.rm = TRUE)) %>% 
      group_by(model) %>% 
      summarise(hr_PgC = mean(hr_PgC))
    ) %>% 
    ggplot(aes(x = model, y=hr_PgC)) + 
    geom_bar(stat = "sum"
             , width = 0.5
             , size = 1
             , fill = "gray"
             , col = "black") +
    theme(legend.position = "none") +
    labs(x = "Model", y = expression(Global~annual~R[H]~(Pg~C~yr^{-1}))) ->
    p2
  print(p2)
  # ggsave("outputs/p2.pdf")
}
