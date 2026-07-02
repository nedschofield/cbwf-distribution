# 3.2_CBWF_BRT_data_narclim.R
# Author: Ned Ryan-Schofield
#
# Overview:
# Download NARCliM2.0 future climate data to align with the historical AGCD
# covariates used in the BRT model. Builds a manifest of download URLs across
# GCMs, SSPs, RCMs, variables, and years, checks which links exist, filters to
# the most current version, and downloads all files.
#
# Future periods: 2040-2059 (mid-century) and 2080-2099 (late-century)
# Variables: daily max temperature (tasmaxAdjust), precipitation (prAdjust)
# SSPs: ssp126, ssp245, ssp370
#
# Note: output file paths must be kept under 260 characters (Windows limit).
# Folder structure is kept short to avoid this issue.

# ---- Packages ------------------------------------------------------------

library(dplyr)
library(myutils) # github.com/nedschofield/myutils

# ---- Define download parameters ------------------------------------------

base_url <- "https://thredds.nci.org.au/thredds/fileServer/zz63/NARCliM2-0-derived/output-CMIP6/bias-adjusted-output/AUS-18/NSW-Government"

path_local_future_dir <- "./data/narclim"

years <- c(2040:2059, 2080:2099)

gcms <- c(
  "ACCESS-ESM1-5",
  "EC-Earth3-Veg",
  "MPI-ESM1-2-HR",
  "NorESM2-MM",
  "UKESM1-0-LL"
)

ssps <- c("ssp126", "ssp245", "ssp370")

rcms <- c(
  "NARCliM2-0-WRF412R3",
  "NARCliM2-0-WRF412R5"
)

variables <- c(
  "tasmaxAdjust",
  "prAdjust"
)

frequency    <- "day"
bias_adj_tag <- "v1-r1-NSWGovernment-CDF-AGCDv1-1990-2009"
versions     <- c("v20241219", "v20241122")
end_date     <- c("1231", "1230")

gcm_members <- tibble::tribble(
  ~gcm,              ~member,
  "ACCESS-ESM1-5",   "r6i1p1f1",
  "EC-Earth3-Veg",   "r1i1p1f1",
  "MPI-ESM1-2-HR",   "r1i1p1f1",
  "NorESM2-MM",      "r1i1p1f1",
  "UKESM1-0-LL",     "r1i1p1f2"
)

# ---- Build download manifest ---------------------------------------------

download_manifest <- tidyr::expand_grid(
    gcm      = gcms,
    ssp      = ssps,
    rcm      = rcms,
    variable = variables,
    year     = years,
    version  = versions,
    end_date = end_date
  ) |>
  left_join(gcm_members, by = "gcm") |>
  mutate(
    rcm_short = dplyr::case_when(
      grepl("R3$", rcm) ~ "R3",
      grepl("R5$", rcm) ~ "R5",
      TRUE ~ rcm
    ),
    
    var_short = dplyr::recode(
      variable,
      tasmaxAdjust = "tmax",
      prAdjust     = "rain",
      .default     = variable
    ),
    
    folder_url = paste(
      base_url, gcm, ssp, member, rcm, bias_adj_tag, frequency, variable, version,
      sep = "/"
    ),
    
    filename = paste0(
      variable, "_AUS-18_",
      gcm, "_",
      ssp, "_",
      member, "_NSW-Government_",
      rcm, "_",
      bias_adj_tag, "_",
      frequency, "_",
      year, "0101-",
      year, end_date, ".nc"
    ),
    
    download_url = paste0(folder_url, "/", filename),
    
    out_dir = file.path(path_local_future_dir, gcm, ssp, rcm_short, var_short)
  )

# ---- Check and filter manifest -------------------------------------------

# Check which links exist — takes a while
download_manifest$exists <- RCurl::url.exists(download_manifest$download_url)

copy_dlman <- download_manifest

# Keep only existing URLs, then take the most current version per file
download_manifest <- download_manifest[download_manifest$exists == TRUE, ] |>
  mutate(
    version_date = as.integer(stringr::str_remove(version, "^v"))
  ) |>
  group_by(gcm, ssp, rcm, variable, year) |>
  slice_max(order_by = version_date, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(-version_date) |>
  arrange(gcm, ssp, rcm, variable, version, year)

# ---- Download files ------------------------------------------------------

manifest_split <- split(download_manifest, download_manifest$out_dir)

out_paths <- lapply(manifest_split, function(x) {
  myutils::nci_thredds_download(
    download_urls = x$download_url,
    out_dir       = unique(x$out_dir),
    overwrite     = FALSE)
})

out_paths <- unlist(out_paths, use.names = FALSE)
save(out_paths, file = "./outputs/table/narclim_filepaths.Rdata")
