#' author: Ned Ryan-Schofield
#' date: 5/02/2026
#' 
#' Function for computing slope of linear model for each pixel in a time series raster stach
#' output is a single raster with each pixel the value of the slope of the lm
#' 
#' @param x character vector of file location, or terra SpatRaster ordered in a time series

compute_change_rast_timeseries <- function(x) {
  
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required")
  }
  
  # Normalise input to a SpatRaster
  if (inherits(x, "SpatRaster")) {
    r_stack <- x
    
  } else if (is.character(x)) {
    r_stack <- terra::rast(x)
    
  } else {
    stop(
      "`x` must be either a file path (character) ",
      "or a terra::SpatRaster"
    )
  }
  
  time_vec <- seq_len(terra::nlyr(r_stack))
  
  terra::app(r_stack, function(v) {
    if (sum(!is.na(v)) < 2) return(NA_real_)
    coef(lm(v ~ time_vec))[2]
  })
}