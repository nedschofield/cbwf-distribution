# 1.2_CBWF_cov_proccess_temperature.R
# Author: Ned Ryan-Schofield
# Date: 6-2-2026
#
# Overview:
# Temperature was suggested as a driver of CBWF declines in the birds action plan.
# This script imports AGCD daily max temperature from the NCI THREDDS server and derives two
# rasters: (1) slope of change in annual extreme heat days (days > 40 degrees C)
# over 1990-2024, and (2) slope of change in mean annual maximum temperature.
# Absolute change is also computed for each (slope x years).

# ---- Packages and setup --------------------------------------------------

library(tidyverse)
library(sf)
library(terra)

source("./scripts/functions/nci_thredds_download.R")
source("./scripts/functions/compute_change_rast_timeseries.R")

study.region <- vect("./outputs/vector/study.region.gpkg")
crs.loc      <- "EPSG:28353"

# ---- Download daily max temperature --------------------------------------

tmax_years <- 1990:2024

download_urls <- paste0("https://thredds.nci.org.au/thredds/fileServer/zv2/agcd/v1-0-3/tmax/mean/r005/01day/agcd_v1_tmax_mean_r005_daily_", tmax_years, ".nc")

path_local_tmax_dir <- "./data/raster/agcd/max_temp"

out_paths <- nci_thredds_download(download_urls = download_urls, out_dir = path_local_tmax_dir, overwrite = FALSE)

# ---- Extreme heat days per cell ------------------------------------------

temp.thresh <- 40

ex.heat.days.list <- list()
for (j in seq_along(tmax_years)) {
 
  y1 <- rast(out_paths[j])
  message("counting days over 40 for year index: ", j)
  ex.heat.days.list[[j]] <- app(y1, fun = function(x) sum(x > temp.thresh, na.rm = TRUE))
  
}

ex.heat.days.stack <- do.call(c, ex.heat.days.list)

# Reproject study region for initial cropping
study.region.ex.heat.reproject <- study.region %>%
  project(crs(ex.heat.days.stack)) %>%
  makeValid() # correct invalid geometry from reprojection

ex.heat.days.stack <- ex.heat.days.stack %>%
  crop(study.region.ex.heat.reproject %>% buffer(10000)) %>%
  project(crs.loc) %>% 
  crop(study.region %>% buffer(5000))

# Compute slope of change and absolute change over the study period
ex.heat.days     <- compute_change_rast_timeseries(ex.heat.days.stack)
abs.ex.heat.days <- app(ex.heat.days, fun = function(x) x * length(tmax_years))

writeRaster(ex.heat.days,     file.path("./outputs/raster/agcd/extreme_heat_days_proccessed.tif"),          overwrite = TRUE)
writeRaster(abs.ex.heat.days, file.path("./outputs/raster/agcd/absolute_extreme_heat_days_proccessed.tif"), overwrite = TRUE)

# ---- Change in mean maximum temperature ----------------------------------

# Reproject study region for initial cropping
study.region.temp.reproject <- study.region %>%
  project(crs(rast(out_paths[1]))) %>%
  makeValid()

max.temp.list <- list()
for (j in seq_along(tmax_years)) {
  
  message("computing annual mean max temp for year index: ", j)
  max.temp.list[[j]] <- rast(out_paths[j]) %>%
    crop(study.region.temp.reproject %>% buffer(10000)) %>%
    project(crs.loc) %>% 
    crop(study.region %>% buffer(5000)) %>%
    app(fun = function(x) mean(x, na.rm = TRUE))  # annual mean of daily max temps
  
}

max.temp.stack <- do.call(c, max.temp.list)

max.temp     <- compute_change_rast_timeseries(max.temp.stack)
max.temp.abs <- app(max.temp, fun = function(x) x * length(tmax_years))

writeRaster(max.temp,     file.path("./outputs/raster/agcd/max_temp_proccessed.tif"),          overwrite = TRUE)
writeRaster(max.temp.abs, file.path("./outputs/raster/agcd/absolute_max_temp_proccessed.tif"), overwrite = TRUE)
