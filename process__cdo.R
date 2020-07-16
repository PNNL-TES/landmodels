## Script to open and process heterotrophic respiration data from MIMICS, CORPSE, and CASA-CNP
## Created July 13, 2020 | Stephanie Pennington 

# We need to calculate an average Rh for each model output, this should give us
# One number for each year in the dataset, i.e. 11 numbers for each grid cell
cat("Running cdo...")

system("cdo -yearmean data/Wieder/mimics_pool_flux_2000-2010_daily.nc data/Wieder/mimics_pool_flux_2000-2010_mean.nc")
system("cdo -yearmean data/Wieder/corpse_pool_flux_2000-2010_daily.nc data/Wieder/corpse_pool_flux_2000-2010_mean.nc")
system("cdo -yearmean data/Wieder/casaclm_pool_flux_2000-2010_daily.nc data/Wieder/casaclm_pool_flux_2000-2010_mean.nc")

