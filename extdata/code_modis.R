# this is the code to get GPP and NPP from MODIS data for MGRsD
library(tidyverse)

# setwd('C:/Users/jian107/Documents/PNNL/landmodels/extdata') # need update accordingly

SoilData <- read_csv("extdata/MGRsD_SRDBV5_20200819.csv") %>%
  mutate(GPP = NA, NPP = NA, QC_npp = NA)

nrows <- nrow(SoilData) 


for (row_i in 1:nrows) {
 site_code <- SoilData$Study_number[row_i]
 meas_year <-  SoilData$Meas_Year[row_i]
 meas_doy <-  SoilData$Meas_DOY[row_i]
 meas_doy <- floor((meas_doy - 1)/8) * 8 + 1
 if (!file.exists(paste0("MODIS_GPP_NPP/",site_code,".csv"))){next}
 if (meas_year >= 2000) {
 npps <- read_csv(paste0("MODIS_GPP_NPP/",site_code,".csv")) %>%
   mutate(Year =  as.numeric(format(date,'%Y')), DOY =  as.numeric(format(date,'%j'))) %>%
   filter(Year == meas_year, DOY == meas_doy)
 SoilData$GPP[row_i] <- mean(npps$Gpp, na.rm = TRUE)
 SoilData$NPP[row_i] <- mean(npps$PsnNet, na.rm = TRUE)
 SoilData$QC_npp[row_i] <- mean(npps$Psn_QC, na.rm = TRUE)
 }
 else
 {
   npps <- read_csv(paste0("MODIS_GPP_NPP/",site_code,".csv")) %>%
     mutate(Year =  as.numeric(format(date,'%Y')), DOY =  as.numeric(format(date,'%j'))) %>%
     filter(DOY == meas_doy)
   SoilData$GPP[row_i] <- mean(npps$Gpp, na.rm = TRUE)
   SoilData$NPP[row_i] <- mean(npps$PsnNet, na.rm = TRUE)
 }
}
