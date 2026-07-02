# 1.3_CBWF_cov_proccess_rainfall.R
# Author: Ned Ryan-Schofield
# Date: 6-2-2026
#
# Overview:
# Change in productivity due to reduced annual rainfall (as a climate change proxy)
# could also affect whiteface populations.
# This script imports AGCD monthly rainfall data from the NCI THREDDS server and derives
# change in annual rainfall from 1990-2024 for each pixel.
# It also computes 6-month rainfall totals (Jun-Sep) for 2024, matched to survey months

# ---- Packages and setup --------------------------------------------------

library(tidyverse)
library(sf)
library(terra)

source("./scripts/functions/nci_thredds_download.R")
source("./scripts/functions/compute_change_rast_timeseries.R")

study.region <- vect("./outputs/vector/study.region.gpkg")
crs.loc      <- "EPSG:28353"

# ---- Download rainfall data ----------------------------------------------

rainfall_years <- 1990:2024

download_urls <- paste0("https://thredds.nci.org.au/thredds/fileServer/zv2/agcd/v2-0-3/precip/total/r001/01month/agcd_v2_precip_total_r001_monthly_", rainfall_years, ".nc")

path_local_rainfall_dir <- "./data/raster/agcd/rainfall"

out_paths <- nci_thredds_download(download_urls = download_urls, out_dir = path_local_rainfall_dir, overwrite = FALSE)

# ---- Change in annual rainfall per cell ----------------------------------

# Reproject study region for initial cropping
study.region.rain.reproject <- study.region %>%
  project(crs(rast(out_paths[1]))) %>%
  makeValid()

annual.rain.list <- list()
for (j in seq_along(rainfall_years)) {
  
  y1 <- rast(out_paths[j]) %>%
    crop(study.region.rain.reproject %>% buffer(10000)) %>%
    project("EPSG:28353") %>%
    crop(study.region %>% buffer(5000))
  
  message("summing rainfall totals for year index: ", j)
  annual.rain.list[[j]] <- app(y1, fun = function(x) sum(x, na.rm = TRUE))
  
}

annual.rain.stack <- do.call(c, annual.rain.list)

annual.rain <- compute_change_rast_timeseries(annual.rain.stack)

writeRaster(annual.rain, file.path("./outputs/raster/agcd/annual_rainfall_proccessed.tif"), overwrite = TRUE)

# ---- Recent 6-month rainfall totals --------------------------------------
# CBWF may be highly mobile and move to areas with recent rainfall.
# Compute rolling 6-month totals for each survey month (June-September 2024).

rain2024 <- rast(file.path("./data/raster/agcd/rainfall/agcd_v2_precip_total_r001_monthly_2024.nc")) %>%
  project(crs(study.region)) %>%
  crop(buffer(study.region, 10000))

target_months <- 6:9  # June to September

# Get indices of the 6 months prior to (and including) each target month
get_6mo_indices <- function(month) {
  ((month - 5):month) %% 12 %>% replace(. == 0, 12)
}

six_month_sums <- lapply(target_months, function(m) {
  indices <- get_6mo_indices(m)
  sum(rain2024[[indices]])
})

names(six_month_sums) <- paste0("rain6mo_total_", month.name[target_months])

rain6mo_totals <- rast(six_month_sums)

writeRaster(rain6mo_totals, filename = "./outputs/raster/agcd/rain6mo_totals.tif", overwrite = TRUE)
