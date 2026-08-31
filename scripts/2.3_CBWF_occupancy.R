# 2.3_CBWF_occupancy.R
# Author: Ned Ryan-Schofield
# Date: 11-5-2026
#
# Overview:
# Bayesian single-season occupancy model fitted in JAGS. Occupancy probability (psi)
# is modelled as a function of annual rainfall, 6-month pre-survey rainfall, extreme
# heat days, photosynthetic vegetation (PV) fractional cover, and non-photosynthetic
# vegetation (NPV) fractional cover. Detection probability is held constant across
# all sites and visits.

# ---- Packages ------------------------------------------------------------

library(terra) 
library(sf)
library(tidyverse)
library(lubridate)
library(jagsUI)
library(corrplot)

# ---- Read in data --------------------------------------------------------

covars_sa <- rast("./outputs/raster/covars_sa.tif")
sites     <- vect("./outputs/vector/CBWF_site_locations.gpkg")
obs       <- read.csv("./data/table/CBWF_surveys.csv")

# ---- Covariate extraction ------------------------------------------------

covars <- terra::extract(covars_sa, sites, bind = TRUE) |>
  as.data.frame()

coords <- terra::crds(sites, df = TRUE)
covars <- cbind(covars, coords)

# Attach survey date (one per site)
dates  <- obs |> distinct(site, date)
covars <- merge(covars, dates, by = "site")

# Select 6-month rainfall layer matching the survey month at each site
rain_selected <- covars |>
  mutate(month = month(dmy(date), label = TRUE, abbr = FALSE)) |>
  rowwise() |>
  mutate(rain6mo = case_when(
    month == "June"      ~ rain6mo_Jun,
    month == "July"      ~ rain6mo_Jul,
    month == "August"    ~ rain6mo_Aug,
    month == "September" ~ rain6mo_Sep,
    TRUE                 ~ NA_real_
  )) |>
  ungroup() |>
  select(site, rain6mo)

covars <- merge(covars, rain_selected, by = "site") |>
  select(-CBWF, -date) |>
  arrange(site)  # fix row order before aligning with y

# ---- Detection history ---------------------------------------------------
# All 10 survey cells per site are used as repeat visits

y <- subset(obs, select = c(site, period, CBWF)) %>%
  reshape(idvar = "site", timevar = "period", direction = "wide") %>%
  select(-1)
y[y > 0] <- 1
y <- as.matrix(y)

nSites  <- nrow(y)
nVisits <- ncol(y)

cat(sprintf("%d sites, %d visits per site\n", nSites, nVisits))
cat(sprintf("Sites with any detection: %d / %d\n",
            sum(apply(y, 1, max, na.rm = TRUE) > 0), nSites))

# ---- Collinearity check --------------------------------------------------

cov_check <- covars |>
  select(annual.rain, rain6mo, abs.ex.heat, pv, npv, bare)

corrplot(cor(cov_check, use = "complete.obs"),
         method = "number", type = "lower", diag = FALSE,
         tl.cex = 0.8, number.cex = 0.7)

# PV and NPV are strongly correlated with bare; bare is excluded. PV and NPV are retained
# as separate moderately correlated predictors — they represent distinct ecological processes
# (primary productivity vs. standing dry biomass driving the brown food web). Annual
# rainfall and 6-month rainfall have some correlation; both are retained because they
# represent different processes (long-term habitat suitability vs. recent resource
# availability).

# ---- Scale covariates ----------------------------------------------------

num_cols   <- c("annual.rain", "rain6mo", "abs.ex.heat", "pv", "npv")
scaled_mat <- scale(covars[, num_cols])

# Save centering/scaling params for later application to rasters
site_center <- attr(scaled_mat, "scaled:center")
site_scale  <- attr(scaled_mat, "scaled:scale")

covars <- covars |>
  mutate(
    annual.rain_s = scaled_mat[, "annual.rain"],
    rain6mo_s     = scaled_mat[, "rain6mo"],
    abs.ex.heat_s = scaled_mat[, "abs.ex.heat"],
    pv_s          = scaled_mat[, "pv"],
    npv_s         = scaled_mat[, "npv"]
  )

# ---- JAGS data and initialisations ---------------------------------------

jags_data_global <- list(
  y       = y,
  nSites  = nSites,
  nVisits = nVisits,
  ann     = covars$annual.rain_s,
  rain6   = covars$rain6mo_s,
  heat    = covars$abs.ex.heat_s,
  pv      = covars$pv_s,
  npv     = covars$npv_s
)

jags_data_null <- jags_data_global[c("y", "nSites", "nVisits")]

# Detected sites: init z = 1 (must be occupied); non-detected: init z = 0
det_max <- apply(y, 1, max, na.rm = TRUE)
z_init  <- ifelse(det_max > 0, 1L, 0L)

inits_null <- function() list(
  z      = z_init,
  alpha0 = 0,
  beta0  = 0
)

inits_global <- function() list(
  z          = z_init,
  alpha0     = 0,
  beta0      = 0,
  beta_ann   = 0,
  beta_rain6 = 0,
  beta_heat  = 0,
  beta_pv    = 0,
  beta_npv   = 0
)

# ---- Fit null model ------------------------------------------------------

params_null <- c("alpha0", "p", "beta0", "nocc", "mean_psi")

occu_null <- jags(
  data               = jags_data_null,
  inits              = inits_null,
  parameters.to.save = params_null,
  model.file         = "./scripts/2.1_occu_null.jags",
  n.chains  = 3,
  n.adapt   = 1000,
  n.iter    = 20000,
  n.burnin  = 5000,
  n.thin    = 5,
  parallel  = TRUE
)

# ---- Fit global model ----------------------------------------------------

params_global <- c(
  "alpha0", "p",
  "beta0", "beta_ann", "beta_rain6", "beta_heat", "beta_pv", "beta_npv",
  "psi", "nocc", "mean_psi",
  "bpv", "T_obs", "T_rep"
)

occu_global <- jags(
  data               = jags_data_global,
  inits              = inits_global,
  parameters.to.save = params_global,
  model.file         = "./scripts/2.2_occu_global.jags",
  n.chains  = 3,
  n.adapt   = 1000,
  n.iter    = 20000,
  n.burnin  = 5000,
  n.thin    = 5,
  parallel  = TRUE
)

# ---- Model comparison (WAIC) ---------------------------------------------
# WAIC is preferred over DIC for hierarchical Bayesian models with latent states
# because it is fully Bayesian and has a direct connection to leave-one-out CV.
# The latent occupancy state z is marginalised out analytically for each MCMC
# sample. Lower WAIC indicates better predictive performance; p_WAIC is the
# effective number of parameters.

# P(y_i | params) = psi_i * P(y_i | z=1, p) + (1-psi_i) * I(all y_i == 0)
compute_waic <- function(post_mat, lp_mat, y_mat) {
  n_vis   <- apply(y_mat, 1, function(x) sum(!is.na(x)))
  n_det   <- apply(y_mat, 1, sum, na.rm = TRUE)
  any_det <- n_det > 0

  psi_mat <- plogis(lp_mat)
  p_vec   <- plogis(post_mat[, "alpha0"])

  ll_z1 <- outer(n_det,          log(p_vec)) +
            outer(n_vis - n_det, log(1 - p_vec))

  ll_z0_vec <- ifelse(any_det, -Inf, 0)

  log_lik_mat <- log(
    psi_mat * exp(ll_z1) + (1 - psi_mat) * exp(ll_z0_vec)
  )

  lse    <- function(x) { m <- max(x); m + log(mean(exp(x - m))) }
  lppd   <- sum(apply(log_lik_mat, 1, lse))
  p_waic <- sum(apply(log_lik_mat, 1, var))
  list(waic = -2 * (lppd - p_waic), lppd = lppd, p_waic = p_waic)
}

post_null   <- as.matrix(occu_null$samples)
post_global <- as.matrix(occu_global$samples)

lp_null <- matrix(rep(post_null[, "beta0"], each = nSites), nrow = nSites)

lp_global <- outer(rep(1, nSites),      post_global[, "beta0"])      +
             outer(covars$annual.rain_s, post_global[, "beta_ann"])   +
             outer(covars$rain6mo_s,     post_global[, "beta_rain6"]) +
             outer(covars$abs.ex.heat_s, post_global[, "beta_heat"])  +
             outer(covars$pv_s,          post_global[, "beta_pv"])    +
             outer(covars$npv_s,         post_global[, "beta_npv"])

waic_null   <- compute_waic(post_null,   lp_null,   y)
waic_global <- compute_waic(post_global, lp_global, y)

cat(sprintf("Null model   WAIC: %.1f  (p_WAIC = %.1f)\n", waic_null$waic,   waic_null$p_waic))
cat(sprintf("Global model WAIC: %.1f  (p_WAIC = %.1f)\n", waic_global$waic, waic_global$p_waic))
cat(sprintf("DWAIC (null - global): %.1f\n", waic_null$waic - waic_global$waic))
cat("(Positive DWAIC favours the global model)\n")

# ---- Diagnostics ---------------------------------------------------------

# Convergence: all key parameters should have Rhat < 1.1
key_params <- c("alpha0", "p", "beta0", "beta_ann", "beta_rain6",
                "beta_heat", "beta_pv", "beta_npv", "mean_psi", "nocc")

rhat_df <- tibble(
  parameter = key_params,
  Rhat      = unlist(occu_global$Rhat[key_params])
)

rhat_df |>
  mutate(converged = Rhat < 1.1) |>
  arrange(desc(Rhat))

# Traceplots
traceplot(occu_global,
          parameters = c("beta0", "beta_ann", "beta_rain6",
                         "beta_heat", "beta_pv", "beta_npv", "alpha0"))

# Posterior predictive check (Freeman-Tukey statistic)
# Values near 0.5 indicate good fit; near 0 or 1 suggest systematic misfit
bpv_val <- occu_global$mean$bpv
cat(sprintf("Bayesian p-value (Freeman-Tukey): %.3f\n", bpv_val))
cat(sprintf("  T_obs mean: %.2f\n", occu_global$mean$T_obs))
cat(sprintf("  T_rep mean: %.2f\n", occu_global$mean$T_rep))

# ---- Results -------------------------------------------------------------

# Detection probability
cat(sprintf(
  "Detection probability p = %.3f  (95%% CrI: %.3f - %.3f)\n",
  occu_global$mean$p,
  occu_global$q2.5$p,
  occu_global$q97.5$p
))

# Occupancy covariate effects
coef_params <- c("beta_ann", "beta_rain6", "beta_heat", "beta_pv", "beta_npv")
coef_labels <- c("Annual rainfall", "6-mo rainfall",
                 "Extreme heat days", "PV fractional cover",
                 "NPV fractional cover")

coef_df <- tibble(
  label = factor(coef_labels, levels = rev(coef_labels)),
  mean  = unlist(occu_global$mean[coef_params]),
  lo    = unlist(occu_global$q2.5[coef_params]),
  hi    = unlist(occu_global$q97.5[coef_params])
)

coef_plot <- ggplot(coef_df, aes(x = mean, y = label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(xmin = lo, xmax = hi)) +
  labs(x = "Coefficient (logit scale, 95% CrI)", y = NULL,
       title = "Occupancy covariate effects — global model")
coef_plot

# Estimated number of occupied sites
cat(sprintf(
  "Occupied sites: %.1f  (95%% CrI: %.0f - %.0f)  of %d surveyed\n",
  occu_global$mean$nocc,
  occu_global$q2.5$nocc,
  occu_global$q97.5$nocc,
  nSites
))

# Site-level occupancy probability
psi_df <- tibble(
  site     = covars$site,
  mean     = occu_global$mean$psi,
  lo       = occu_global$q2.5$psi,
  hi       = occu_global$q97.5$psi,
  detected = as.logical(apply(y, 1, max, na.rm = TRUE))
) |>
  arrange(mean) |>
  mutate(site = factor(site, levels = site))

ggplot(psi_df, aes(y = site, x = mean, xmin = lo, xmax = hi,
                   colour = detected)) +
  geom_pointrange(size = 0.35) +
  scale_colour_manual(
    values = c("FALSE" = "#CC6633", "TRUE" = "#336699"),
    labels = c("Not detected", "Detected")
  ) +
  labs(x = "Occupancy probability psi (95% CrI)", y = "Site",
       colour = NULL,
       title = "Site-level posterior occupancy probability") +
  theme(axis.text.y = element_text(size = 7))

# ---- Spatial prediction --------------------------------------------------
# Rasters are aggregated to 1 km2 cells (10 x 10 from 1 ha) to approximate the
# effective survey footprint. The September 6-month rainfall layer is used for
# spatial prediction, as the majority of surveys were conducted in September.
# Uncertainty is propagated by drawing 500 thinned samples from the joint
# posterior and computing per-cell quantiles.

covars_pred <- covars_sa[[c("annual.rain", "rain6mo_Sep", "abs.ex.heat", "pv", "npv")]]
covars_agg  <- aggregate(covars_pred, fact = 10, fun = mean, na.rm = TRUE)

# Scale raster values using survey-data centering/scaling params
agg_df <- as.data.frame(covars_agg, na.rm = FALSE) |>
  mutate(
    annual.rain_s = (annual.rain - site_center["annual.rain"])  / site_scale["annual.rain"],
    rain6mo_s     = (rain6mo_Sep - site_center["rain6mo"])      / site_scale["rain6mo"],
    abs.ex.heat_s = (abs.ex.heat - site_center["abs.ex.heat"])  / site_scale["abs.ex.heat"],
    pv_s          = (pv          - site_center["pv"])           / site_scale["pv"],
    npv_s         = (npv         - site_center["npv"])          / site_scale["npv"]
  )

# Posterior mean prediction
b0      <- occu_global$mean$beta0
b_ann   <- occu_global$mean$beta_ann
b_rain6 <- occu_global$mean$beta_rain6
b_heat  <- occu_global$mean$beta_heat
b_pv    <- occu_global$mean$beta_pv
b_npv   <- occu_global$mean$beta_npv

lp_mean <- b0 +
  b_ann   * agg_df$annual.rain_s +
  b_rain6 * agg_df$rain6mo_s     +
  b_heat  * agg_df$abs.ex.heat_s +
  b_pv    * agg_df$pv_s          +
  b_npv   * agg_df$npv_s

psi_mean_vec <- 1 / (1 + exp(-lp_mean))

psi_r <- covars_agg[[1]]
values(psi_r) <- psi_mean_vec

# Uncertainty: draw 500 thinned samples from the joint posterior
post_mat  <- as.matrix(occu_global$samples)
thin_idx  <- round(seq(1, nrow(post_mat), length.out = 500))
post_thin <- post_mat[thin_idx,
                      c("beta0", "beta_ann", "beta_rain6", "beta_heat", "beta_pv", "beta_npv")]

X_pred <- cbind(
  1,
  agg_df$annual.rain_s,
  agg_df$rain6mo_s,
  agg_df$abs.ex.heat_s,
  agg_df$pv_s,
  agg_df$npv_s
)

LP_mat  <- X_pred %*% t(post_thin)
psi_mat <- 1 / (1 + exp(-LP_mat))

psi_lo_r <- psi_hi_r <- covars_agg[[1]]
values(psi_lo_r) <- apply(psi_mat, 1, quantile, 0.025, na.rm = TRUE)
values(psi_hi_r) <- apply(psi_mat, 1, quantile, 0.975, na.rm = TRUE)

rm(LP_mat, psi_mat)

# Maps: posterior mean and 95% CrI
# Anchor labels to each raster's own extent (rather than par("usr")), since
# terra pads the plot region to preserve aspect ratio under mfrow, which
# otherwise pushes par("usr")-based labels above the visible panel.
add_panel_label <- function(r, label) {
  e <- ext(r)
  text(
    x = e$xmax - 0.06 * (e$xmax - e$xmin),
    y = e$ymax - 0.06 * (e$ymax - e$ymin),
    labels = label, font = 2, cex = 1.3
  )
}

par(mfrow = c(1, 3))
plot(psi_r,    main = expression(paste("Posterior mean ", psi)), range = c(0, 1))
points(sites, pch = 20, cex = 0.6)
add_panel_label(psi_r, "a")
plot(psi_lo_r, main = expression(paste("Lower 95% CI ", psi)), range = c(0, 1))
points(sites, pch = 20, cex = 0.6)
add_panel_label(psi_lo_r, "b")
plot(psi_hi_r, main = expression(paste("Upper 95% CI ", psi)), range = c(0, 1))
points(sites, pch = 20, cex = 0.6)
add_panel_label(psi_hi_r, "c")
par(mfrow = c(1, 1))
occu_psi_plot <- recordPlot()

# Estimated area with psi >= 0.5 (threshold-based occupied area)
cell_area_km <- cellSize(psi_r, unit = "km")
thresh       <- 0.5

area_mean <- as.numeric(global(ifel(psi_r    >= thresh, 1, NA) * cell_area_km, "sum", na.rm = TRUE))
area_lo   <- as.numeric(global(ifel(psi_lo_r >= thresh, 1, NA) * cell_area_km, "sum", na.rm = TRUE))
area_hi   <- as.numeric(global(ifel(psi_hi_r >= thresh, 1, NA) * cell_area_km, "sum", na.rm = TRUE))

cat(sprintf(
  "Estimated area with psi >= 0.5:  %.0f km2  (95%% CrI: %.0f - %.0f km2)\n",
  area_mean, area_lo, area_hi
))

# ---- Save outputs --------------------------------------------------------

writeRaster(psi_r,    "./outputs/raster/psi_mean.tif",  overwrite = TRUE)
writeRaster(psi_lo_r, "./outputs/raster/psi_lo95.tif",  overwrite = TRUE)
writeRaster(psi_hi_r, "./outputs/raster/psi_hi95.tif",  overwrite = TRUE)

saveRDS(occu_global, "./outputs/occu_global.rds")
saveRDS(occu_null,   "./outputs/occu_null.rds")

saveRDS(list(waic_null = waic_null, waic_global = waic_global),
        "./outputs/waic_comparison.rds")

saveRDS(list(center = site_center, scale = site_scale),
        "./outputs/covariate_scaling_params.rds")

ggsave("figures/occu_global_covar_coeffs.png", coef_plot,
       width = 22, height = 18, units = "cm", dpi = 300)

png("figures/occupancy_psi_maps.png",
    width = 24, height = 10, units = "cm", res = 300)
replayPlot(occu_psi_plot)
dev.off()
