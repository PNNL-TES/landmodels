#*****************************************************************************************************************
# prepare data 
#*****************************************************************************************************************
# read file
read_file <- function(x) read.csv(file.path(DATA_DIR, x), comment.char = "#", stringsAsFactors = FALSE)
writ_file <- function(input, output) write.csv(input, file.path(OUT_DIR, output), row.names = FALSE)
'%!in%' <- function(x,y)!('%in%'(x,y))

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
    idoy <- which(abs(doy - target_doy) <= 3) # find the doy within 3 days of DOY (e.g., if DOY = 15, then DOY 12,12,14,15,15,17,18 will be used)
    
    Rh_cell <- ncvar_get(landm, var_rh[modelname], start = c(ilon, ilat, itime), count = c(1, 1, 1))
    Rh_avg <- mean(ncvar_get(landm, var_rh[modelname], start = c(ilon, ilat, 1), count = c(1, 1, -1))[idoy], na.rm=T) # Rh average of idoy
    
    mgrhd[i, modelname] <- ifelse(target_time > 2011 | target_time < 2000, Rh_avg, Rh_cell) # since model only have prediction from 2000-2010, other times using average
    
    print(paste0("*****", i))
  }
  return(mgrhd)
}

#*****************************************************************************************************************
# site scale analysis 
#*****************************************************************************************************************
# Analyze the residuals trend
# i = 1
# sdata = residual_data
slr_residual_site <- function(sdata) {
  out <- data.frame()
  unique(sdata$Study_number) %>% sort() ->
    var_study
  
  for (i in 1:length(var_study)) {
    sdata %>% 
      filter(Study_number == var_study[i]) ->
      sub_data
    
    n_obs <- (sub_data %>% nrow())/3
    Latitude <- sub_data$Latitude[1]
    Longitude <- sub_data$Longitude[1]
    Biome <- sub_data$Biome[1]
    Ecosystem_type <- sub_data$Ecosystem_type[1]
    Leaf_habit <- sub_data$Leaf_habit[1]
    MAT <- sub_data$MAT[1]
    MAP <- sub_data$MAP[1]
    
    # t test
    res_t <- t.test(sub_data$Residuals, conf.level = 0.95)
    MBE <- res_t$estimate
    p_t_test <- res_t$p.value %>% round(3)
    
    # model for residual
    m_res <- lm(Rh_Norm ~ Residuals, data = sub_data)
    # summary(m_casa) %>% print()
    inter_res <- summary(m_res)$coefficients[1,1] %>% round(3)
    slope_res <- summary(m_res)$coefficients[2,1] %>% round(3)
    slope_p_res <- summary(m_res)$coefficients[2,4] %>% round(3)
    R2_res <- summary(m_res)$r.squared %>% round(3)
    
    # output results
    out <- rbind(out, data.frame(var_study[i], n_obs,
                                 MBE, p_t_test,
                                 inter_res, slope_res, slope_p_res, R2_res,
                                 Latitude, Longitude, Biome, Ecosystem_type, Leaf_habit, MAT, MAP) )
    print(paste0(i," ********** ", var_study[i]))
  }
  
  return(out)
}

# Regression between measured and modeled Rh
slr_landmodels_site <- function(sdata) {
  out <- data.frame()
  
  for (i in 1:nrow(lm_site)) {
    sdata %>% 
      filter(Study_number == lm_site$Study_number[i]) %>% 
      dplyr::select(Study_number, Meas_Year, Meas_DOY, Rh_Norm, casaclm, corpse, mimics) %>% 
      filter(!is.na(Rh_Norm)) ->
      sub_data
    
    n_obs <- sub_data %>% nrow()
    
    # model for casaclm
    m_casa <- lm(Rh_Norm ~ casaclm, data = sub_data)
    # summary(m_casa) %>% print()
    inter_casa <- summary(m_casa)$coefficients[1,1] %>% round(3)
    slope_casa <- summary(m_casa)$coefficients[2,1] %>% round(3)
    slope_p_casa <- summary(m_casa)$coefficients[2,4] %>% round(3)
    R2_casa <- summary(m_casa)$r.squared %>% round(3)
    
    # model for corpse
    m_corpse <- lm(Rh_Norm ~ corpse, data = sub_data)
    # summary(m_corpse) %>% print()
    inter_corpse <- summary(m_corpse)$coefficients[1,1] %>% round(3)
    slope_corpse <- summary(m_corpse)$coefficients[2,1] %>% round(3)
    slope_p_corpse <- summary(m_corpse)$coefficients[2,4] %>% round(3)
    R2_corpse <- summary(m_corpse)$r.squared %>% round(3)
    
    # model for mimics
    m_mimics <- lm(Rh_Norm ~ mimics, data = sub_data)
    # summary(m_mimics) %>% print()
    inter_mimics <- summary(m_mimics)$coefficients[1,1] %>% round(3)
    slope_mimics <- summary(m_mimics)$coefficients[2,1] %>% round(3)
    slope_p_mimics <- summary(m_mimics)$coefficients[2,4] %>% round(3)
    R2_mimics <- summary(m_mimics)$r.squared %>% round(3)
    
    # output results
    out <- rbind(out, data.frame(sub_data$Study_number[1], n_obs,
                                 inter_casa, slope_casa, slope_p_casa, R2_casa,
                                 inter_corpse, slope_corpse, slope_p_corpse, R2_corpse,
                                 inter_mimics, slope_mimics, slope_p_mimics, R2_mimics) )
    print(paste0(i," ********** ", sub_data$Study_number[1]))
  }
  
  return(out)
}


# site by site comparison
# time series for each study, residue plot
Rh_site_comp <- function(sdata, study_number) {
  # plot time series
  sdata %>% 
    filter(Study_number == study_number) %>% 
    rename(Rh_field = Rh_Norm, CASACNP = casaclm, CORPSE = corpse, MIMICS = mimics) %>% 
    dplyr::select(Meas_Year, Meas_DOY, Rh_field, CASACNP, CORPSE, MIMICS) %>% 
    filter(!is.na(Rh_field)) %>% 
    mutate(Time = as.Date(Meas_DOY, origin = paste0(Meas_Year, "-01-01"))) %>% 
    tidyr::gather(key = "Model", value = "Rh", -Meas_Year, -Meas_DOY, -Time) %>% 
    ggplot(aes(x=Time, y=Rh, col = Model)) +
    geom_line(aes(linetype = Model)) + 
    labs(x = (paste0("Study ", study_number)),
         y = expression(R[H]~(g~C~m^{-2}~day^{-1}))) +
    scale_x_date(date_labels = "%b/%d/%Y") +
    # custome color
    scale_colour_discrete(name="Method",
                          breaks=c("CASACNP", "CORPSE", "MIMICS", "Rh_field"),
                          labels=c("CASA-CNP", "CORPSE", "MIMICS", "Measured")) +
    scale_color_manual(values=c("red", "orange4", "springgreen4", "black")) +
    # # custome linetype
    scale_linetype_discrete(name = "Method",
                            breaks=c("CASACNP", "CORPSE", "MIMICS", "Rh_field"),
                            labels=c("CASA-CNP", "CORPSE", "MIMICS", "Measured")) +
    scale_linetype_manual(values = c(2,2,2,1)) ->
    p
  
  # plot linear regression
  sdata %>% 
    filter(Study_number == study_number) %>% 
    dplyr::select(Meas_Year, Meas_DOY, Rh_Norm, casaclm, corpse, mimics) %>% 
    filter(!is.na(Rh_Norm)) %>% 
    mutate(Time = as.Date(Meas_DOY, origin = paste0(Meas_Year, "-01-01"))) %>% 
    tidyr::gather(key = "Model", value = "Rh_model", -Meas_Year, -Meas_DOY, -Time, -Rh_Norm) %>% 
    ggplot(aes(x=Rh_model, y=Rh_Norm, col = Model)) +
    geom_point(aes(shape = Model)) + 
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = expression(Model~R[H]~(g~C~m^{-2}~day^{-1})),
         y = expression(Measured~R[H]~(g~C~m^{-2}~day^{-1}))) +
    # custome color
    # scale_shape_discrete(name="",
    #                       breaks=c("casaclm", "corpse", "mimics", "Rh_Norm"),
    #                       labels=c("CASA-CNP", "CORPSE", "MIMICS", "Measured")) +
    scale_color_manual(values=c("red", "orange4", "springgreen4", "black")) +
    theme(legend.position = "none") ->
    p_lm
  
  plot_grid(p, p_lm, ncol = 1) ->
    p_comb
  
  # plot difference
  sdata %>% 
    filter(Study_number == study_number) %>% 
    filter(!is.na(Rh_Norm)) %>% 
    mutate(CASACNP = (casaclm - Rh_Norm)/mean(Rh_Norm, na.rm = TRUE),
           CORPSE = (corpse - Rh_Norm)/mean(Rh_Norm, na.rm = TRUE),
           MIMICS = (mimics - Rh_Norm)/mean(Rh_Norm, na.rm = TRUE)
    ) %>% 
    dplyr::select(Meas_DOY, CASACNP, CORPSE, MIMICS) %>%
    tidyr::gather(key = "Model", value = "Rh_dif", -Meas_DOY) %>%
    group_by(Model, Meas_DOY) %>% 
    summarise(Rh_dif = mean(Rh_dif)) %>% 
    ggplot(aes(x=Meas_DOY, y=Rh_dif, col = Model)) +
    geom_line() +
    geom_hline(yintercept = 0, col="blue", linetype = "twodash", size = 1) +
    labs(x = expression(Day~of~year),
         y = expression(R[H]~difference~(g~C~m^{-2}~day^{-1}))) +
    # custome color
    scale_color_discrete(name="",
                         breaks=c("CASACNP", "CORPSE", "MIMICS"),
                         labels=c("CASA-CNP", "CORPSE", "MIMICS")) +
    scale_color_manual(values=c("red", "orange4", "springgreen4")) ->
    p_dif
  
  # plot dif vs predicted value
  sdata %>% 
    filter(Study_number == study_number) %>% 
    filter(!is.na(Rh_Norm)) %>% 
    mutate(CASACNP = casaclm - Rh_Norm,
           CORPSE = corpse - Rh_Norm,
           MIMICS = mimics - Rh_Norm
    ) %>% 
    dplyr::select(Meas_DOY, CASACNP, CORPSE, MIMICS, Rh_Norm) %>%
    tidyr::gather(key = "Model", value = "Rh_dif", -Meas_DOY, -Rh_Norm) %>%
    group_by(Model, Meas_DOY) %>% 
    summarise(Rh_dif = mean(Rh_dif), Rh_Norm = mean(Rh_Norm)) ->
    sum_data
  
  # linear regression of Rh_Norm and residual
  sdata %>% 
    filter(Study_number == study_number) %>% 
    filter(!is.na(Rh_Norm)) %>% 
    mutate(CASACNP = casaclm - Rh_Norm,
           CORPSE = corpse - Rh_Norm,
           MIMICS = mimics - Rh_Norm
    ) %>% 
    dplyr::select(Meas_DOY, CASACNP, CORPSE, MIMICS, Rh_Norm) %>%
    tidyr::gather(key = "Model", value = "Rh_dif", -Meas_DOY, -Rh_Norm) ->
    sub_data
  
  m_res <- lm(Rh_Norm ~ Rh_dif, data = sub_data)
  # summary(m_casa) %>% print()
  inter_res <- summary(m_res)$coefficients[1,1] %>% round(3)
  slope_res <- summary(m_res)$coefficients[2,1] %>% round(3)
  slope_p_res <- summary(m_res)$coefficients[2,4] %>% round(3)
  R2_res <- summary(m_res)$r.squared %>% round(3)
  
  sum_data %>% 
    ggplot(aes(x=Rh_Norm, y=Rh_dif, col = Model)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, col = ifelse(slope_p_res <= 0.05, "red", "gray")) +
    geom_hline(yintercept = c(0, mean(sum_data$Rh_dif)),
               col=c("blue", "black"),
               linetype = c("twodash", "solid"),
               size = 1) +
    labs(y = expression(Residuals~(g~C~m^{-2}~day^{-1})),
         x = expression(Measured~R[H]~(g~C~m^{-2}~day^{-1}))) +
    # custome color
    scale_color_discrete(name="",
                         breaks=c("CASACNP", "CORPSE", "MIMICS"),
                         labels=c("CASA-CNP", "CORPSE", "MIMICS")) +
    scale_color_manual(values=c("red", "orange4", "springgreen4")) +
    theme(legend.position = "none") ->
    p_dif_pred
  
  plot_grid(p, p_dif_pred,
            # p_dif, p_lm, 
            ncol = 2,
            rel_widths = c(1.35, 1)) ->
    p_comb
  
  # print(p_comb)
}

