# 3.3_CBWF_BRT_narclim_ensemble_averaging.R
# Author: Ned Ryan-Schofield
#
# Overview:
# NARCliM2.0 climate variable averaging
# Produces per-model (60) and ensemble mean (6) NetCDFs for:
#   - Mean annual days with Tmax > 40 degrees C
#   - Mean annual rainfall (mm)
# Outputs written to outputs/raster/narclim/ and outputs/raster/narclim/ensemble/

# ---- Packages and setup --------------------------------------------------

library(ncdf4)
library(dplyr)

source("scripts/functions/narclim_functions.R")

load(file = "./outputs/table/narclim_filepaths.Rdata")

# ---- File path filtering -------------------------------------------------

temp_paths <- out_paths[grepl("tasmax|tmax", out_paths, ignore.case = TRUE)]
rain_paths <- out_paths[grepl("pradjust",    out_paths, ignore.case = TRUE)]

# ---- Build manifests -----------------------------------------------------

parse_narclim_paths <- function(paths, var_folder) {
  tibble(path = paths) |>
    mutate(
      gcm    = sub("./data/raster/narclim/([^/]+)/.*", "\\1", path),
      ssp    = sub(paste0(".*/([^/]+)/R[0-9]+/", var_folder, "/.*"), "\\1", path),
      rcm    = sub(paste0(".*/", var_folder, "/.*NARCliM2-0-WRF412(R[0-9]+)_.*"), "\\1", path),
      year   = as.integer(sub(".*_day_(\\d{4})\\d{4}-.*", "\\1", path)),
      period = case_when(
        year %in% 2040:2059 ~ "2040-2059",
        year %in% 2080:2099 ~ "2080-2099"
      )
    )
}

tmax_manifest <- parse_narclim_paths(temp_paths, "tmax")
rain_manifest <- parse_narclim_paths(rain_paths, "rain")

# ---- Step 1a: per-model mean days over 40 (60 files) --------------------

tmax_manifest |>
  group_by(gcm, ssp, rcm, period) |>
  group_walk(function(group_df, keys) {
    out_path <- file.path(
      "outputs/raster/narclim",
      paste0(keys$gcm, "_", keys$ssp, "_", keys$rcm, "_", keys$period, ".nc")
    )
    message("Processing: ", basename(out_path))
    narclim_mean_days_over_threshold(group_df$path, out_path, threshold_c = 40)
  })

# ---- Step 1b: per-model mean annual rainfall (60 files) -----------------

rain_manifest |>
  group_by(gcm, ssp, rcm, period) |>
  group_walk(function(group_df, keys) {
    out_path <- file.path(
      "outputs/raster/narclim/rainfall",
      paste0(keys$gcm, "_", keys$ssp, "_", keys$rcm, "_", keys$period, ".nc")
    )
    message("Processing: ", basename(out_path))
    narclim_mean_annual_rainfall(group_df$path, out_path)
  })

# ---- Step 2a: ensemble mean days over 40 (6 files) ----------------------

list.files("outputs/raster/narclim", pattern = "\\.nc$", full.names = TRUE) |>
  (\(x) tibble(path = x))() |>
  mutate(
    ssp    = sub(".*_(ssp[0-9]+)_.*",       "\\1", basename(path)),
    period = sub(".*(\\d{4}-\\d{4})\\.nc",  "\\1", basename(path))
  ) |>
  group_by(ssp, period) |>
  group_walk(function(group_df, keys) {
    out_path <- file.path(
      "outputs/raster/narclim/ensemble",
      paste0(keys$ssp, "_", keys$period, "_ensemble_mean.nc")
    )
    narclim_ensemble_mean(group_df$path, out_path)
  })

# ---- Step 2b: ensemble mean rainfall (6 files) --------------------------

list.files("outputs/raster/narclim/rainfall", pattern = "\\.nc$", full.names = TRUE) |>
  (\(x) tibble(path = x))() |>
  mutate(
    ssp    = sub(".*_(ssp[0-9]+)_.*",       "\\1", basename(path)),
    period = sub(".*(\\d{4}-\\d{4})\\.nc",  "\\1", basename(path))
  ) |>
  group_by(ssp, period) |>
  group_walk(function(group_df, keys) {
    out_path <- file.path(
      "outputs/raster/narclim/ensemble",
      paste0(keys$ssp, "_", keys$period, "_rainfall_ensemble_mean.nc")
    )
    narclim_ensemble_mean_rainfall(group_df$path, out_path)
  })
