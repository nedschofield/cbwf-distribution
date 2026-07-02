#' author: Ned Ryan-Schofield
#' date: 5/02/2026

#' Helper function for getting urls of MODIS fractional cover dataset from the thredds database
#' Catalogue url to browse the server is https://thredds.nci.org.au/thredds/catalog/catalog.html
#' To download, have to use https://thredds.nci.org.au/thredds/fileServer/... instead of catalog
#' Product of interest is Rangeland and Pasture Productivity RaPP = tc43
#' some tiles are not .006.nc, they are .061.nc. MODIS v6 discontinued feb2023, and have switched to v6.1. So use v6.1 when available (since 2020)
#' Function automatically replaces with v6.1 for years >= 2020
#' 
#' @param mod_h numeric vector of modis horizontal cells of interest
#' @param mod_v numeric vector of modis vertical cells of interest
#' @param years numeric vector of years of interest. i.e. 2001:2024, or c(2020, 2024) etc

modis_fc_thredds_urls <- function(mod_h, mod_v, years) {
  base <- "https://thredds.nci.org.au/thredds/fileServer/tc43/modis-fc/v310/tiles/monthly/cover/FC_Monthly_Medoid.v310.MCD43A4.h"
  
  combo <- expand.grid(mod_h = mod_h, mod_v = mod_v, year = years)
  
  urls <- paste0(
    base,
    combo$mod_h,
    "v",
    combo$mod_v,
    ".",
    combo$year,
    ".006.nc"
  )
  
  urls[combo$year >= 2020] <- stringr::str_replace(
    urls[combo$year >= 2020],
    "\\.006\\.nc$",
    ".061.nc"
  )
  
  urls
}