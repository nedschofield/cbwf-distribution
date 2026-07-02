# 1.1_CBWF_cov_proccess_MODIS_FC.R
# Author: Ned Ryan-Schofield
# Date: 27-5-2025
#
# Overview:
# Overgrazing was suggested as a driver of CBWF decline in the Birds action plan.
# To investigate whether habitat quality has changed over the study area across
# the last ~20 years (MODIS stretches back to 2001), MODIS derived fractional cover
# is used as a proxy for habitat quality.

# ---- Packages and setup --------------------------------------------------

library(tidyverse)
library(sf)
library(terra)

source("./scripts/functions/nci_thredds_download.R")
source("./scripts/functions/modis_fc_thredds_urls.R")
source("./scripts/functions/compute_change_rast_timeseries.R")

# ---- Create study region boundary ----------------------------------------
# Select IBRA subregions where whiteface are found
# (Everard Block, Piderka, Gawler Volcanics, Northern Flinders, Kingoonya,
# Torrens, Roxby, Commonwealth Hill, Warriner, Breakaways, Oodnadatta,
# Murnpeowie, Peake-Dennison Inlier, Macumba, Baltana)

ibra_regions <- c('CER03', 'FIN04', 'FLB05', 'GAW02', 'GAW05', 'GAW06', 'GAW07', 'GAW08', 'SSD04', 'STP01', 'STP02', 'STP03', 'STP04', 'STP05', 'STP07')

crs.loc <- "EPSG:28353"

study.region <- read_sf("./data/vector/ibra_subregions.gpkg") |>
  st_transform(crs.loc) |>
  dplyr::filter(SUB_CODE_7 %in% ibra_regions) |>
  st_geometry()

write_sf(study.region, "./outputs/vector/study.region.gpkg", overwrite = TRUE)

# ---- Download MODIS fractional cover -------------------------------------

mod_h <- 29:30 # Modis horizontal cell
mod_v <- 11:12 # Modis vertical cell/s
year  <- 2001:2024
download_urls <- modis_fc_thredds_urls(mod_h = mod_h, mod_v = mod_v, years = year)

path_local_FC_dir <- "./data/raster/Guerschman_FC_Monthly"

out_paths <- nci_thredds_download(download_urls = download_urls, out_dir = path_local_FC_dir, overwrite = FALSE)

# ---- Read and merge MODIS rasters ----------------------------------------

dataPath <- out_paths

# There are four different grid cells, make a raster stack for each
dataPath_1 <- dataPath[which(str_detect(dataPath, paste0(mod_h[1], "v", mod_v[1])))]
dataPath_2 <- dataPath[which(str_detect(dataPath, paste0(mod_h[1], "v", mod_v[2])))]
dataPath_3 <- dataPath[which(str_detect(dataPath, paste0(mod_h[2], "v", mod_v[1])))]
dataPath_4 <- dataPath[which(str_detect(dataPath, paste0(mod_h[2], "v", mod_v[2])))]

fc_list <- list()
for (j in seq_along(year)) {
 
  r1 <- rast(dataPath_1[j])
  r2 <- rast(dataPath_2[j])
  r3 <- rast(dataPath_3[j])
  r4 <- rast(dataPath_4[j])
  
  message("Merging rasters for year index: ", j)
  fc_list[[j]] <- merge(r1, r2, r3, r4)
}

fc <- do.call(c, fc_list)

# ---- Separate fractional cover components --------------------------------

pv   <- fc %>%  .[[which(str_detect(names(.), "^phot"))]]
npv  <- fc %>%  .[[which(str_detect(names(.), "^nphot"))]]
bare <- fc %>%  .[[which(str_detect(names(.), "^bare"))]]

# ---- Reproject and clip to study region ----------------------------------

# Reproject study area to MODIS sinusoidal
study.region.modis <- study.region %>% 
  st_transform("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs")

# Crop to 10km buffered extent, reproject to GDA94 MGA zone 53, trim edges
reproject_crop <- function(x) {
  x |> 
    crop(study.region.modis |>  st_buffer(10000)) |> 
    project(crs.loc) |> 
    crop(study.region |>  st_buffer(5000)) # 5km buffer removes edge distortion
}

pv   <- reproject_crop(pv)
npv  <- reproject_crop(npv)
bare <- reproject_crop(bare)

# ---- Rename fractional cover rasters -------------------------------------

rast_names <- function(x) {
  month      <- lubridate::month(time(x))
  year       <- lubridate::year(time(x))
  cover_type <- deparse(substitute(x)) %>% stringr::str_to_upper()
  paste0(cover_type, "_Monthly_Medoid.v310.MCD43A4.", year, ".", month, ".tif")
}

names(pv)   <- rast_names(pv)
names(npv)  <- rast_names(npv)
names(bare) <- rast_names(bare)

# ---- Write processed FC to disc ------------------------------------------

start_date <- time(pv)[1] %>% str_replace_all("-", "")
end_date   <- time(pv)[length(time(pv))] %>% str_replace_all("-", "")

outputPath <- file.path("./outputs/raster/guerschman_fc_monthly")

writeCDF(pv,   file.path(outputPath, paste0("fc_monthly_pv_",   start_date, "_", end_date, ".nc")), varname = "guerschman_monthly_pv",   compression = 6, overwrite = T)
writeCDF(npv,  file.path(outputPath, paste0("fc_monthly_npv_",  start_date, "_", end_date, ".nc")), varname = "guerschman_monthly_npv",  compression = 6, overwrite = T)
writeCDF(bare, file.path(outputPath, paste0("fc_monthly_bare_", start_date, "_", end_date, ".nc")), varname = "guerschman_monthly_bare", compression = 6, overwrite = T)

# ---- Calculate change in fractional cover over time ----------------------

pv   <- "./outputs/raster/guerschman_fc_monthly/fc_monthly_pv_20010101_20241201.nc"
npv  <- "./outputs/raster/guerschman_fc_monthly/fc_monthly_npv_20010101_20241201.nc"
bare <- "./outputs/raster/guerschman_fc_monthly/fc_monthly_bare_20010101_20241201.nc"

pv_change   <- compute_change_rast_timeseries(pv)
npv_change  <- compute_change_rast_timeseries(npv)
bare_change <- compute_change_rast_timeseries(bare)

writeRaster(pv_change,   file.path("./outputs/raster/guerschman_fc_monthly/pv_change_proccessed.tif"),   overwrite = TRUE)
writeRaster(npv_change,  file.path("./outputs/raster/guerschman_fc_monthly/npv_change_proccessed.tif"),  overwrite = TRUE)
writeRaster(bare_change, file.path("./outputs/raster/guerschman_fc_monthly/bare_change_proccessed.tif"), overwrite = TRUE)
