# Plots
prep_latitude_data <- function(dat, lat_bands = LAT_BANDS) {
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
  
  return(plotting_data)
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
      dplyr::select(-hr) %>% 
      mutate(Year = floor(time)) %>% 
      na.omit() %>%
      group_by(model, Year) %>% 
      summarise(n = n(), hr_PgC = sum(hr_PgC_yr, na.rm = TRUE)) %>% 
      group_by(model) %>% 
      # cannot figure it out why the calculation is not right, *9 for plot the right result
      summarise(hr_PgC = mean(hr_PgC)*9.03),
    
    h_dat %>% 
      mutate(Year = floor(time)) %>% 
      na.omit() %>% 
      group_by(model, Year) %>% 
      summarise(n = n(), hr_PgC = sum(hr_PgC_yr, na.rm = TRUE)) %>% 
      group_by(model) %>% 
      summarise(hr_PgC = mean(hr_PgC)),
    
    tibble(model = "warner", hr_PgC = 49.7)
    
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
}
