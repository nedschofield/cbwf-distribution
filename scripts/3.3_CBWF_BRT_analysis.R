# 3.3_CBWF_BRT_analysis.R
# Author: Ned Ryan-Schofield
#
# Overview:
# Fits a boosted regression tree (BRT) model of CBWF occupancy using historical
# AGCD climate data (extreme heat days, annual rainfall), predicts to three
# historical periods (1960, 1990, 2024), then predicts to NARCliM2.0 future
# climate scenarios (SSP1-2.6, SSP2-4.5, SSP3-7.0 x mid-century, late-century).

# ---- Packages ------------------------------------------------------------

library(terra)
library(sf)
library(tidyverse)
library(ncdf4)
library(dismo)
library(gbm)
library(patchwork)

# ---- Load data -----------------------------------------------------------

records     <- read_csv("./data/table/CBWF_records_sdm.csv")
sdm.region  <- vect("./outputs/vector/sdm_region.gpkg")
records.vec <- vect("./outputs/vector/cbwf_records_sdm.gpkg")

heat <- rast("./outputs/raster/agcd/ex_heat_mean_periods_proccessed.tif")
rain <- rast("./outputs/raster/agcd/annual_rainfall_mean_periods_proccessed.tif")

# ---- Prepare climate stack -----------------------------------------------

annual.rain.mean.aligned <- project(rain, heat, method = "bilinear")
climate_stack <- c(heat, annual.rain.mean.aligned)

# ---- Extract climate values for survey sites -----------------------------

sdm.dat <- terra::extract(climate_stack, records.vec, bind = T, xy = T) |>
  as.data.frame()

# Pick the period-matched heat and rain value for each record
period_chr <- as.character(sdm.dat$period)
heat_cols  <- paste0("heat_", period_chr)
rain_cols  <- paste0("rain_", period_chr)
heat_idx   <- match(heat_cols, names(sdm.dat))
rain_idx   <- match(rain_cols, names(sdm.dat))
rows       <- seq_len(nrow(sdm.dat))

sdm.dat$heat <- as.numeric(sdm.dat[cbind(rows, heat_idx)])
sdm.dat$rain <- as.numeric(sdm.dat[cbind(rows, rain_idx)])

# ---- Fit BRT -------------------------------------------------------------

cbwf.brt <- gbm.step(
  data            = sdm.dat,
  gbm.x           = 15:16,
  gbm.y           = 3,
  family          = "bernoulli",
  tree.complexity = 3,
  learning.rate   = 0.001,
  bag.fraction    = 0.5
)

# ---- Validate BRT --------------------------------------------------------

cbwf.brt$cv.statistics
summary(cbwf.brt)
gbm.plot(cbwf.brt, n.plots = 2, write.title = F)
gbm.perspec(cbwf.brt, 1, 2)

# ---- Predict to historical periods (1960, 1990, 2024) -------------------

historical_periods <- c("1960", "1990", "2024")

pred_historical <- lapply(historical_periods, function(period) {
  clim    <- climate_stack[[c(paste0("heat_", period), paste0("rain_", period))]]
  clim_df <- as.data.frame(clim, na.rm = FALSE)
  names(clim_df) <- c("heat", "rain")

  preds       <- predict.gbm(cbwf.brt, clim_df,
                              n.trees = cbwf.brt$gbm.call$best.trees,
                              type    = "response")
  pred_rast        <- clim[[1]]
  values(pred_rast) <- preds
  names(pred_rast)  <- paste0("pred_", period)

  plot(pred_rast, main = paste("BRT prediction:", period))
  plot(records.vec, add = T)

  pred_rast
})

brt.out <- do.call(c, pred_historical)

# 2024 presences and all absences for overlaying on plots
pres_2024 <- records.vec[records.vec$cbwf_pres == 1 & records.vec$period == 2024, ]
abs_all   <- records.vec[records.vec$cbwf_pres == 0, ]

plot(brt.out, nc = 3, nr = 1,
     range = c(0,0.85),
     fun = function() {
       plot(abs_all,   add = TRUE, pch = 1, cex = 0.6, col = "black")
       plot(pres_2024, add = TRUE, pch = 1, cex = 0.6, col = "red")
     })
terra_hist_plot <- recordPlot()

# ---- Predict to NARCliM future scenarios --------------------------------
# Ensemble NetCDFs: one for days_over_threshold (heat), one for rainfall (rain)
# per SSP x period combination. Predictions written to outputs/narclim/predictions/

ssps    <- c("ssp126", "ssp245", "ssp370")
periods <- c("2040-2059", "2080-2099")

for (ssp in ssps) {
  for (period in periods) {

    # Load ensemble mean heat and rainfall
    nc_tmax <- nc_open(file.path("outputs/narclim/ensemble",
                                  paste0(ssp, "_", period, "_ensemble_mean.nc")))
    nc_rain <- nc_open(file.path("outputs/narclim/ensemble",
                                  paste0(ssp, "_", period, "_rainfall_ensemble_mean.nc")))

    heat_vals  <- as.vector(ncvar_get(nc_tmax, "days_over_threshold"))
    rain_vals  <- as.vector(ncvar_get(nc_rain, "mean_annual_rainfall"))
    lat2d      <- ncvar_get(nc_tmax, "lat")
    lon2d      <- ncvar_get(nc_tmax, "lon")
    rlon_vals  <- nc_tmax$dim$rlon$vals
    rlat_vals  <- nc_tmax$dim$rlat$vals
    rlon_units <- nc_tmax$dim$rlon$units
    rlat_units <- nc_tmax$dim$rlat$units
    nc_close(nc_tmax)
    nc_close(nc_rain)

    # Build prediction data frame; exclude NAs from prediction
    pred_df <- data.frame(heat = heat_vals, rain = rain_vals)
    valid   <- complete.cases(pred_df)

    preds        <- rep(NA_real_, nrow(pred_df))
    preds[valid] <- predict.gbm(cbwf.brt, pred_df[valid, ],
                                 n.trees = cbwf.brt$gbm.call$best.trees,
                                 type    = "response")

    # Reshape to [rlon, rlat] matrix and write NetCDF
    pred_mat <- matrix(preds, nrow = length(rlon_vals), ncol = length(rlat_vals))

    dim_rlon <- ncdim_def("rlon", rlon_units, rlon_vals)
    dim_rlat <- ncdim_def("rlat", rlat_units, rlat_vals)
    var_pred <- ncvar_def("occupancy_prob", "probability", list(dim_rlon, dim_rlat),
                           missval  = NA,
                           longname = "CBWF predicted occupancy probability")
    var_lat  <- ncvar_def("lat", "degrees_north", list(dim_rlon, dim_rlat),
                           longname = "geographic latitude")
    var_lon  <- ncvar_def("lon", "degrees_east",  list(dim_rlon, dim_rlat),
                           longname = "geographic longitude")

    out_path <- file.path("outputs/narclim/predictions",
                           paste0(ssp, "_", period, "_cbwf_occupancy.nc"))
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

    nc_out <- nc_create(out_path, vars = list(var_pred, var_lat, var_lon))
    ncvar_put(nc_out, var_pred, pred_mat)
    ncvar_put(nc_out, var_lat,  lat2d)
    ncvar_put(nc_out, var_lon,  lon2d)
    nc_close(nc_out)

    message("Written: ", out_path)
  }
}

# ---- Reproject NARCliM predictions to brt.out grid ----------------------
# Reads the 6 prediction NetCDFs, uses the embedded 2D lat/lon arrays to
# create true geographic points, rasterizes to a 0.25 degree WGS84 intermediate
# grid (coarser than NARCliM spacing to avoid diagonal gaps), then projects
# and resamples to match brt.out exactly.

# WGS84 intermediate template at 0.25 degree resolution
template_wgs84 <- rast(
  xmin = 110, xmax = 160,
  ymin = -45, ymax = -10,
  resolution = 0.25,
  crs = "EPSG:4326"
)

narclim_layers <- lapply(ssps, function(ssp) {
  lapply(periods, function(period) {

    nc <- nc_open(file.path("outputs/narclim/predictions",
                             paste0(ssp, "_", period, "_cbwf_occupancy.nc")))
    preds <- as.vector(ncvar_get(nc, "occupancy_prob"))
    lon   <- as.vector(ncvar_get(nc, "lon"))
    lat   <- as.vector(ncvar_get(nc, "lat"))
    nc_close(nc)

    # Create geographic points from true lat/lon, drop NAs
    df  <- data.frame(lon = lon, lat = lat, pred = preds)
    df  <- df[!is.na(df$pred), ]
    pts <- vect(df, geom = c("lon", "lat"), crs = "EPSG:4326")

    # Rasterize to WGS84 intermediate grid, then align to brt.out
    pred_wgs84   <- rasterize(pts, template_wgs84, field = "pred", fun = mean)
    pred_aligned <- project(pred_wgs84, brt.out[[1]], method = "bilinear")
    names(pred_aligned) <- paste0(ssp, "_", period)

    pred_aligned
  })
}) |> unlist(recursive = FALSE)

# Stack all 6 scenarios into one SpatRaster
# Loop order is ssp-outer, period-inner, giving:
#   [1] ssp126_2040-2059  [2] ssp126_2080-2099
#   [3] ssp245_2040-2059  [4] ssp245_2080-2099
#   [5] ssp370_2040-2059  [6] ssp370_2080-2099
# Reorder so mid-century (2040-2059) is top row, late-century (2080-2099) bottom row:
#   [1] ssp126_2040  [2] ssp245_2040  [3] ssp370_2040
#   [4] ssp126_2080  [5] ssp245_2080  [6] ssp370_2080
narclim.out <- do.call(c, narclim_layers)[[c(1, 3, 5, 2, 4, 6)]]
narclim.out

plot(narclim.out,
     range = c(0,0.85),
     fun = function() {
       plot(abs_all,   add = TRUE, pch = 1, cex = 0.6, col = "black")
       plot(pres_2024, add = TRUE, pch = 1, cex = 0.6, col = "red")
     })
terra_fut_plot <- recordPlot()

# ---- Satellite basemap with translucent occupancy overlay ---------------

library(maptiles)
library(tidyterra)
library(ggplot2)

# Download Esri satellite tiles for the study region extent
tiles <- get_tiles(brt.out, provider = "Esri.WorldImagery", zoom = 7, crop = TRUE)

# Coober Pedy for geographic reference
coober_pedy <- project(
  vect(cbind(134.7204, -29.0135), crs = "EPSG:4326"),
  "EPSG:28353"
)
cp_df <- data.frame(
  x     = crds(coober_pedy)[, 1],
  y     = crds(coober_pedy)[, 2],
  label = "Coober Pedy"
)

# Reusable single-panel plot function
plot_occupancy <- function(layer, title) {
  ggplot() +
    geom_spatraster_rgb(data = tiles) +
    geom_spatraster(data = layer, alpha = 0.65) +
    scale_fill_whitebox_c(palette = "viridi", na.value = NA,
                          name = "Occ. prob.", limits = c(0, 1)) +
    geom_spatvector(data = abs_all,
                    aes(colour = "Survey sites (absent)"), shape = 1, size = 0.8) +
    geom_spatvector(data = pres_2024,
                    aes(colour = "Detections 2024"), shape = 1, size = 0.8) +
    scale_colour_manual(
      name   = NULL,
      values = c("Survey sites (absent)" = "black", "Detections 2024" = "red"),
      guide  = guide_legend(override.aes = list(shape = 1, size = 2))
    ) +
    geom_point(data = cp_df, aes(x = x, y = y),
               shape = 3, size = 1.5, colour = "white", stroke = 1) +
    geom_text(data = cp_df, aes(x = x, y = y, label = label),
              hjust = -0.15, colour = "white", fontface = "bold", size = 2.2) +
    labs(title = title, x = NULL, y = NULL) +
    coord_sf(datum = sf::st_crs("EPSG:28353"), expand = FALSE) +
    theme_bw() +
    theme(
      panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      panel.grid.major = element_line(colour = "grey80", linetype = "dotted"),
      panel.grid.minor = element_blank(),
      axis.ticks       = element_line(colour = "black"),
      axis.text        = element_text(colour = "black", size = 5),
      plot.title       = element_text(size = 7, hjust = 0.5),
      legend.text      = element_text(size = 7),
      legend.title     = element_text(size = 7)
    )
}

# Historical predictions: 1 row x 3 cols
hist_plots <- lapply(names(brt.out), function(nm) {
  plot_occupancy(brt.out[[nm]], sub("pred_", "", nm))
})

hist_clim_brt_preds <- wrap_plots(hist_plots, nrow = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")
hist_clim_brt_preds

# NARCliM predictions: 2 rows x 3 cols (mid-century top, late-century bottom)
# Titles removed — replaced with row/column strip labels
future_plots <- lapply(names(narclim.out), function(nm) {
  plot_occupancy(narclim.out[[nm]], title = "") + theme(plot.title = element_blank())
})

# Helper for strip label plots
strip_label <- function(txt, angle = 0) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = txt,
             fontface = "bold", size = 3.2, angle = angle) +
    theme_void()
}

# SSP display names (columns), period names (rows)
ssp_display <- c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0")

header_row <- (plot_spacer() | strip_label(ssp_display[1]) | strip_label(ssp_display[2]) | strip_label(ssp_display[3])) +
  plot_layout(widths = c(0.2, 1, 1, 1))

row1 <- (strip_label("2040\u20132059", 90) | future_plots[[1]] | future_plots[[2]] | future_plots[[3]]) +
  plot_layout(widths = c(0.2, 1, 1, 1))

row2 <- (strip_label("2080\u20132099", 90) | future_plots[[4]] | future_plots[[5]] | future_plots[[6]]) +
  plot_layout(widths = c(0.2, 1, 1, 1))

fut_clim_brt_preds <- (header_row / row1 / row2) +
  plot_layout(heights = c(0.05, 1, 1), guides = "collect") &
  theme(legend.position = "right")

fut_clim_brt_preds

# ---- Save all plots ------------------------------------------------------

# Terra plots use base R graphics — save with png() + replayPlot() + dev.off()
png("figures/brt_historical_occupancy_terra.png",
    width = 32, height = 10, units = "cm", res = 300)
replayPlot(terra_hist_plot)
dev.off()

png("figures/brt_narclim_occupancy_terra.png",
    width = 34, height = 22, units = "cm", res = 300)
replayPlot(terra_fut_plot)
dev.off()

# Satellite ggplot panels use ggsave()
ggsave("figures/brt_historical_occupancy_sat.png", hist_clim_brt_preds,
       width = 32, height = 10, units = "cm", dpi = 300)

ggsave("figures/brt_narclim_occupancy_sat.png", fut_clim_brt_preds,
       width = 34, height = 22, units = "cm", dpi = 300)
