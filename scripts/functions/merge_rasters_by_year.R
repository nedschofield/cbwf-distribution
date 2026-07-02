#' author: Ned Ryan-Schofield
#' date: 5/02/2026

#' Function for merging raster tiles together
#' developed for MODIS imagery
#' output file names must have 4 digits denoting year in filename
#' adjust year_regex as required to denote differentyears for merging
#' output is a list of rasters to then be stacked
#' 
#' @param raster_paths character string of output files. Must include a year in the filename
#' @param year_regex regular expression to search for. Default is "\\d{4}" which means 4 digits. First backslash is R escape character, followed by \d (regex digit 0-9) {4} regex quantifier = exactly 4 

merge_rasters_by_year <- function(
    raster_paths,
    year_regex = "\\d{4}",
    quiet = FALSE
) {
  stopifnot(length(raster_paths) > 0)
  
  years <- stringr::str_extract(basename(raster_paths), year_regex)
  
  if (any(is.na(years))) {
    stop("Could not extract year from some filenames")
  }
  
  paths_by_year <- split(raster_paths, years)
  
  out_list <- vector("list", length(paths_by_year))
  names(out_list) <- names(paths_by_year)
  
  for (yr in names(paths_by_year)) {
    if (!quiet) message("Merging rasters for year: ", yr)
    
    rasters <- lapply(paths_by_year[[yr]], terra::rast)
    out_list[[yr]] <- do.call(terra::merge, rasters)
  }
  
  terra::rast(out_list)
}
