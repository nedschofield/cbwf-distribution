### NARCliM2.0 processing functions

library(ncdf4)

# ---- Per-model processing ------------------------------------------------

#' Compute mean annual days over a temperature threshold across a set of NARCliM
#' daily tmax NetCDF files (one per year), and write the result as a new NetCDF.
#'
#' @param file_paths  Character vector of daily tmax NetCDF file paths (one per year).
#' @param out_path    Output NetCDF file path.
#' @param threshold_c Temperature threshold in degrees Celsius. Default 40.
#' @param varname     Name of the temperature variable in the NetCDF. Default "tasmaxAdjust".
narclim_mean_days_over_threshold <- function(
    file_paths,
    out_path,
    threshold_c = 40,
    varname     = "tasmaxAdjust"
) {
  threshold_k <- threshold_c + 273.15
  n_files     <- length(file_paths)

  # Get dimension metadata and geographic coords from first file
  nc0        <- nc_open(file_paths[1])
  rlon_vals  <- ncvar_get(nc0, "rlon")
  rlat_vals  <- ncvar_get(nc0, "rlat")
  lat2d      <- ncvar_get(nc0, "lat")
  lon2d      <- ncvar_get(nc0, "lon")
  rlon_units <- nc0$dim$rlon$units
  rlat_units <- nc0$dim$rlat$units
  nc_close(nc0)

  # Accumulate days over threshold across all files
  accum <- matrix(0, nrow = length(rlon_vals), ncol = length(rlat_vals))

  for (i in seq_len(n_files)) {
    nc_i  <- nc_open(file_paths[i])
    tmax  <- ncvar_get(nc_i, varname)   # [rlon, rlat, time]
    nc_close(nc_i)
    accum <- accum + apply(tmax > threshold_k, c(1, 2), sum)
    message("  Processed file ", i, " of ", n_files, ": ", basename(file_paths[i]))
  }

  mean_days <- accum / n_files

  # Write output NetCDF
  dim_rlon <- ncdim_def("rlon", rlon_units, rlon_vals)
  dim_rlat <- ncdim_def("rlat", rlat_units, rlat_vals)

  var_days <- ncvar_def(
    name     = "days_over_threshold",
    units    = "days",
    dim      = list(dim_rlon, dim_rlat),
    missval  = NA,
    longname = paste0("Mean annual days with Tmax > ", threshold_c, " degrees C")
  )
  var_lat <- ncvar_def("lat", "degrees_north", list(dim_rlon, dim_rlat),
                        longname = "geographic latitude")
  var_lon <- ncvar_def("lon", "degrees_east",  list(dim_rlon, dim_rlat),
                        longname = "geographic longitude")

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  nc_out <- nc_create(out_path, vars = list(var_days, var_lat, var_lon))
  ncvar_put(nc_out, var_days, mean_days)
  ncvar_put(nc_out, var_lat,  lat2d)
  ncvar_put(nc_out, var_lon,  lon2d)
  nc_close(nc_out)

  message("Written: ", out_path)
  invisible(out_path)
}


#' Compute mean annual rainfall (mm) across a set of NARCliM daily precipitation
#' NetCDF files (one per year), and write the result as a new NetCDF.
#' prAdjust units are kg m-2 s-1; multiplied by 86400 to convert to mm/day.
#'
#' @param file_paths Character vector of daily precipitation NetCDF file paths (one per year).
#' @param out_path   Output NetCDF file path.
narclim_mean_annual_rainfall <- function(file_paths, out_path) {
  n_files <- length(file_paths)
  accum   <- NULL

  for (i in seq_len(n_files)) {
    nc_i <- nc_open(file_paths[i])
    pr   <- ncvar_get(nc_i, "prAdjust")   # [rlon, rlat, time] in kg m-2 s-1

    if (is.null(accum)) {
      lat2d      <- ncvar_get(nc_i, "lat")
      lon2d      <- ncvar_get(nc_i, "lon")
      rlon_vals  <- nc_i$dim$rlon$vals
      rlat_vals  <- nc_i$dim$rlat$vals
      rlon_units <- nc_i$dim$rlon$units
      rlat_units <- nc_i$dim$rlat$units
      accum      <- matrix(0, nrow = nrow(pr), ncol = ncol(pr))
    }

    # Convert to mm/day, sum over time dimension -> annual total mm
    annual_mm <- apply(pr * 86400, c(1, 2), sum, na.rm = FALSE)
    accum     <- accum + replace(annual_mm, is.na(annual_mm), 0)
    nc_close(nc_i)
    message("  Processed file ", i, " of ", n_files, ": ", basename(file_paths[i]))
  }

  mean_rainfall <- accum / n_files

  # Restore NAs from first file
  nc_ref <- nc_open(file_paths[1])
  ref    <- ncvar_get(nc_ref, "prAdjust", start = c(1,1,1), count = c(-1,-1,1))
  nc_close(nc_ref)
  mean_rainfall[is.na(ref)] <- NA

  # Write output NetCDF
  dim_rlon <- ncdim_def("rlon", rlon_units, rlon_vals)
  dim_rlat <- ncdim_def("rlat", rlat_units, rlat_vals)

  var_rain <- ncvar_def("mean_annual_rainfall", "mm", list(dim_rlon, dim_rlat),
                         missval  = NA,
                         longname = "Mean annual rainfall")
  var_lat  <- ncvar_def("lat", "degrees_north", list(dim_rlon, dim_rlat),
                         longname = "geographic latitude")
  var_lon  <- ncvar_def("lon", "degrees_east",  list(dim_rlon, dim_rlat),
                         longname = "geographic longitude")

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  nc_out <- nc_create(out_path, vars = list(var_rain, var_lat, var_lon))
  ncvar_put(nc_out, var_rain, mean_rainfall)
  ncvar_put(nc_out, var_lat,  lat2d)
  ncvar_put(nc_out, var_lon,  lon2d)
  nc_close(nc_out)

  message("Written: ", out_path)
  invisible(out_path)
}


# ---- Ensemble averaging --------------------------------------------------

#' Average days_over_threshold across a set of per-model NetCDF files
#' (i.e. ensemble mean across GCMs/RCMs for one SSP/period combination).
#'
#' @param file_paths Character vector of per-model NetCDF paths to average.
#' @param out_path   Output NetCDF file path.
narclim_ensemble_mean <- function(file_paths, out_path) {
  accum <- NULL

  for (path in file_paths) {
    nc_i <- nc_open(path)
    vals <- ncvar_get(nc_i, "days_over_threshold")

    if (is.null(accum)) {
      lat2d      <- ncvar_get(nc_i, "lat")
      lon2d      <- ncvar_get(nc_i, "lon")
      rlon_vals  <- nc_i$dim$rlon$vals
      rlat_vals  <- nc_i$dim$rlat$vals
      rlon_units <- nc_i$dim$rlon$units
      rlat_units <- nc_i$dim$rlat$units
      accum      <- matrix(0, nrow = nrow(vals), ncol = ncol(vals))
    }

    accum <- accum + replace(vals, is.na(vals), 0)
    nc_close(nc_i)
  }

  ensemble_mean <- accum / length(file_paths)

  # Restore NAs from first file
  nc_ref   <- nc_open(file_paths[1])
  ref_vals <- ncvar_get(nc_ref, "days_over_threshold")
  nc_close(nc_ref)
  ensemble_mean[is.na(ref_vals)] <- NA

  # Write output NetCDF
  dim_rlon <- ncdim_def("rlon", rlon_units, rlon_vals)
  dim_rlat <- ncdim_def("rlat", rlat_units, rlat_vals)

  var_days <- ncvar_def("days_over_threshold", "days", list(dim_rlon, dim_rlat),
                         missval  = NA,
                         longname = "Ensemble mean annual days with Tmax > 40 degrees C")
  var_lat  <- ncvar_def("lat", "degrees_north", list(dim_rlon, dim_rlat),
                         longname = "geographic latitude")
  var_lon  <- ncvar_def("lon", "degrees_east",  list(dim_rlon, dim_rlat),
                         longname = "geographic longitude")

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  nc_out <- nc_create(out_path, vars = list(var_days, var_lat, var_lon))
  ncvar_put(nc_out, var_days, ensemble_mean)
  ncvar_put(nc_out, var_lat,  lat2d)
  ncvar_put(nc_out, var_lon,  lon2d)
  nc_close(nc_out)

  message("Written: ", out_path)
  invisible(out_path)
}


#' Average mean_annual_rainfall across a set of per-model NetCDF files
#' (i.e. ensemble mean across GCMs/RCMs for one SSP/period combination).
#'
#' @param file_paths Character vector of per-model NetCDF paths to average.
#' @param out_path   Output NetCDF file path.
narclim_ensemble_mean_rainfall <- function(file_paths, out_path) {
  accum <- NULL

  for (path in file_paths) {
    nc_i <- nc_open(path)
    vals <- ncvar_get(nc_i, "mean_annual_rainfall")

    if (is.null(accum)) {
      lat2d      <- ncvar_get(nc_i, "lat")
      lon2d      <- ncvar_get(nc_i, "lon")
      rlon_vals  <- nc_i$dim$rlon$vals
      rlat_vals  <- nc_i$dim$rlat$vals
      rlon_units <- nc_i$dim$rlon$units
      rlat_units <- nc_i$dim$rlat$units
      accum      <- matrix(0, nrow = nrow(vals), ncol = ncol(vals))
    }

    accum <- accum + replace(vals, is.na(vals), 0)
    nc_close(nc_i)
  }

  ens_mean <- accum / length(file_paths)

  # Restore NAs from first file
  nc_ref <- nc_open(file_paths[1])
  ref    <- ncvar_get(nc_ref, "mean_annual_rainfall")
  nc_close(nc_ref)
  ens_mean[is.na(ref)] <- NA

  # Write output NetCDF
  dim_rlon <- ncdim_def("rlon", rlon_units, rlon_vals)
  dim_rlat <- ncdim_def("rlat", rlat_units, rlat_vals)

  var_rain <- ncvar_def("mean_annual_rainfall", "mm", list(dim_rlon, dim_rlat),
                         missval  = NA,
                         longname = "Ensemble mean annual rainfall")
  var_lat  <- ncvar_def("lat", "degrees_north", list(dim_rlon, dim_rlat),
                         longname = "geographic latitude")
  var_lon  <- ncvar_def("lon", "degrees_east",  list(dim_rlon, dim_rlat),
                         longname = "geographic longitude")

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  nc_out <- nc_create(out_path, vars = list(var_rain, var_lat, var_lon))
  ncvar_put(nc_out, var_rain, ens_mean)
  ncvar_put(nc_out, var_lat,  lat2d)
  ncvar_put(nc_out, var_lon,  lon2d)
  nc_close(nc_out)

  message("Written: ", out_path)
  invisible(out_path)
}
