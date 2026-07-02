# 1.4_CBWF_cov_proccess_final_outputs.R
# Author: Ned Ryan-Schofield
# Date: 6-2-2026
#
# Overview:
# Process additional data files and combine all spatial covariates into a
# single aligned raster stack for downstream analysis. Also creates masking
# layers (water polygons, NVIS habitat types) and a spatial vector of CBWF
# survey sites.

# ---- Packages ------------------------------------------------------------

library(terra)
library(sf)
library(tidyverse)

# ---- Water polygon mask --------------------------------------------------

study.region <- read_sf("./outputs/vector/study.region.gpkg")

water <- read_sf("./data/vector/water_courses.gpkg") %>% 
  mutate(water = 1) %>%
  st_crop(study.region %>% st_buffer(10000))

write_sf(water, "./outputs/vector/water_courses_cropped.gpkg", overwrite = TRUE)

# ---- NVIS habitat mask ---------------------------------------------------

study.region <- vect(study.region)

nvis <- rast("./data/raster/NVIS/NVIS6_0_AUST_EXT_MVG_ALB.tif")

study.region.nvis <- project(study.region, crs(nvis))

nvis_cropped <- nvis %>%
  crop(buffer(study.region.nvis, 10000)) %>%
  project(crs(study.region))

writeRaster(nvis_cropped, file.path("./outputs/raster/nvis_cropped.tif"), overwrite = TRUE)

# ---- CBWF survey sites vector --------------------------------------------

sites <- vect(read_csv("./data/table/CBWF_site_locations.csv"), geom = c("long", "lat"), crs = "EPSG:4326") %>%
  project(crs(study.region))

writeVector(sites, "./outputs/vector/CBWF_site_locations.gpkg", layer = "CBWF_survey_sites", filetype = "GPKG", overwrite = TRUE)

# ---- Load all spatial data -----------------------------------------------

# Vectors
study.region <- vect("./outputs/vector/study.region.gpkg")
water        <- vect("./outputs/vector/water_courses_cropped.gpkg")

# Rasters
nvis           <- rast(file.path("./outputs/raster/nvis_cropped.tif"))
annual.rain    <- rast(file.path("./outputs/raster/agcd/annual_rainfall_proccessed.tif"))
ex.heat        <- rast(file.path("./outputs/raster/agcd/extreme_heat_days_proccessed.tif"))
abs.ex.heat    <- rast(file.path("./outputs/raster/agcd/absolute_extreme_heat_days_proccessed.tif"))
max.temp       <- rast(file.path("./outputs/raster/agcd/max_temp_proccessed.tif"))
abs.max.temp   <- rast(file.path("./outputs/raster/agcd/absolute_max_temp_proccessed.tif"))
rain6mo_totals <- rast(file.path("./outputs/raster/agcd/rain6mo_totals.tif"))
pv             <- rast(file.path("./outputs/raster/guerschman_fc_monthly/pv_change_proccessed.tif"))
npv            <- rast(file.path("./outputs/raster/guerschman_fc_monthly/npv_change_proccessed.tif"))
bare           <- rast(file.path("./outputs/raster/guerschman_fc_monthly/bare_change_proccessed.tif"))

# ---- Combine and align covariates ----------------------------------------

covars_list <- list(nvis, annual.rain, abs.ex.heat, max.temp, pv, npv, bare, rain6mo_totals)

# Resample all layers to match the first (nvis)
covars_aligned <- lapply(covars_list, function(x) resample(x, covars_list[[1]]))

covars_stack <- rast(covars_aligned)
names(covars_stack) <- c("nvis", "annual.rain", "abs.ex.heat", "max.temp", "pv", "npv", "bare",
                         "rain6mo_Jun", "rain6mo_Jul", "rain6mo_Aug", "rain6mo_Sep")

# Drop nvis — categorical data handled separately for masking
covars_stack <- covars_stack[[-1]]

# ---- Mask by habitat type and study region -------------------------------
# Restrict prediction area to NVIS classes occupied by whiteface sites

hab_sites  <- terra::extract(nvis, sites, bind = T) %>% as.data.frame()
hab_labels <- unique(as.character(hab_sites$MVG_NAME))

nvis_levels <- levels(nvis)[[1]]
hab_ids     <- nvis_levels$Value[nvis_levels[["MVG_NAME"]] %in% hab_labels]

hab_mask <- nvis
hab_mask[!(hab_mask %in% hab_ids)] <- NA
hab_mask[hab_mask %in% hab_ids]    <- 1

covars_sa <- mask(covars_stack, hab_mask)
covars_sa <- mask(covars_sa, study.region)

writeRaster(covars_sa, file.path("./outputs/raster/covars_sa.tif"), overwrite = TRUE)

covars_sa <- rast(file.path("./outputs/raster/covars_sa.tif"))
