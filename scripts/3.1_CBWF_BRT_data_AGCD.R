# 3.1_CBWF_BRT_data_AGCD.R
# Author: Ned Ryan-Schofield
#
# Overview:
# Derive climate covariates from historical AGCD data for species distribution
# modelling across three time periods (pre-1990/1960s, 1990 surveys, 2024 surveys).
#
# For each time period:
# - Mean number of extreme heat days (days > 40 degrees C) in a 10-year window
# - Mean annual rainfall in a 10-year window
#
# These are used as predictors in the BRT species distribution model (script 3.3).

# ---- Packages ------------------------------------------------------------

library(tidyverse)
library(sf)
library(terra)

# ---- CBWF records and study region ---------------------------------------

records <- read_csv("./data/table/CBWF_records_sdm.csv")
  
records.vec <- vect(records, geom = c("long", "lat"), crs = "EPSG:4326") %>%
  project("EPSG:28353")

writeVector(records.vec, "./outputs/vector/cbwf_records_sdm.gpkg", overwrite = TRUE)

# Bounding box of records with 50km buffer as SDM region
e          <- ext(records.vec)
buf        <- 50 * 1000
ebuf       <- ext(xmin(e) - buf, xmax(e) + buf, ymin(e) - buf, ymax(e) + buf)
sdm.region <- as.polygons(ebuf, crs = crs(records.vec))

writeVector(sdm.region, "./outputs/vector/sdm_region.gpkg", overwrite = TRUE)

# ---- Download daily max temperature -------------------------------------

year <- 1950:2024

file_urls <- paste0("https://thredds.nci.org.au/thredds/fileServer/zv2/agcd/v1-0-3/tmax/mean/r005/01day/agcd_v1_tmax_mean_r005_daily_", year, ".nc")

path_local_daily_max_dir <- "./data/agcd/max_temp"
if (!file.exists(path_local_daily_max_dir)) dir.create(path_local_daily_max_dir, recursive = TRUE)

out_paths <- file_urls %>% 
  str_split("/", simplify = T) %>%
  .[, 13] %>% 
  file.path(path_local_daily_max_dir, .)

options(timeout = 180)
for (i in seq_along(out_paths)) {
  if (!file.exists(out_paths[i])) download.file(url = file_urls[i], destfile = out_paths[i], mode = "wb")
}

# ---- Extreme heat days per cell ------------------------------------------

ex.heat.days.list <- list()
for (j in seq_along(year)) {
 
  y1 <- rast(out_paths[j])
  message("counting days over 40 for year index: ", j)
  ex.heat.days.list[[j]] <- app(y1, fun = function(x) sum(x > 40, na.rm = TRUE))
  
}

ex.heat.days.stack <- do.call(c, ex.heat.days.list)

sdm.region.ex.heat.reproject <- sdm.region %>%
  project(crs(ex.heat.days.stack))
  
ex.heat.days.stack <- ex.heat.days.stack %>%
  crop(sdm.region.ex.heat.reproject %>% buffer(10000)) %>%
  project("EPSG:28353") %>%
  crop(sdm.region)

# ---- Mean extreme heat days by study period ------------------------------

periods <- list(
  "heat_1960" = 1950:1960,
  "heat_1990" = 1980:1990,
  "heat_2024" = 2014:2024
)

mean_list <- lapply(names(periods), function(pname) {
  yrs <- periods[[pname]]
  idx <- which(year %in% yrs)
  r   <- app(ex.heat.days.stack[[idx]], mean, na.rm = TRUE)
  names(r) <- pname
  r
})

ex.heat.mean.periods <- do.call(c, mean_list)
plot(ex.heat.mean.periods)

writeRaster(ex.heat.mean.periods, file.path("./outputs/raster/agcd/ex_heat_mean_periods_proccessed.tif"), overwrite = TRUE)

# ---- Download rainfall data ----------------------------------------------

year <- 1950:2024

file_urls <- paste0("https://thredds.nci.org.au/thredds/fileServer/zv2/agcd/v2-0-3/precip/total/r001/01month/agcd_v2_precip_total_r001_monthly_", year, ".nc")

path_local_rainfall_dir <- "./data/agcd/rainfall"
if (!file.exists(path_local_rainfall_dir)) dir.create(path_local_rainfall_dir, recursive = TRUE)

out_paths <- file_urls %>% 
  str_split("/", simplify = T) %>%
  .[, 13] %>% 
  file.path(path_local_rainfall_dir, .)

options(timeout = 180)
for (i in seq_along(out_paths)) {
  if (!file.exists(out_paths[i])) download.file(url = file_urls[i], destfile = out_paths[i], mode = "wb")
}

# ---- Annual rainfall per cell --------------------------------------------

SDM.region.rain.reproject <- sdm.region %>%
  project(crs(rast(out_paths[1])))

annual.rain.list <- list()
for (j in seq_along(year)) {
  
  y1 <- rast(out_paths[j]) %>%
    crop(SDM.region.rain.reproject %>% buffer(10000)) %>%
    project("EPSG:28353") %>%
    crop(sdm.region)
  
  message("summing rainfall totals for year index: ", j)
  annual.rain.list[[j]] <- app(y1, fun = function(x) sum(x, na.rm = TRUE))
  
}

annual.rain.stack <- do.call(c, annual.rain.list)

# ---- Mean annual rainfall by study period --------------------------------

periods <- list(
  "rain_1960" = 1950:1960,
  "rain_1990" = 1980:1990,
  "rain_2024" = 2014:2024
)

mean_list <- lapply(names(periods), function(pname) {
  yrs <- periods[[pname]]
  idx <- which(year %in% yrs)
  r   <- app(annual.rain.stack[[idx]], mean, na.rm = TRUE)
  names(r) <- pname
  r
})

annual.rain.mean.periods <- do.call(c, mean_list)
plot(annual.rain.mean.periods)

writeRaster(annual.rain.mean.periods, file.path("./outputs/raster/agcd/annual_rainfall_mean_periods_proccessed.tif"), overwrite = TRUE)
