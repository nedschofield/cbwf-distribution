library(maptiles)
library(tidyterra)
library(ggplot2)
library(dplyr)
library(terra)
library(sf)
library(osmdata)
library(ozmaps)
library(cowplot)

#' Plot site counts on satellite imagery
#'
#' @param sites   SpatVector with a numeric column to visualise
#' @param col     Column name (string) to map to colour and size
#' @param zoom    Tile zoom level passed to maptiles::get_tiles()
#' @param provider Tile provider string (default Esri.WorldImagery)
#' @param alpha   Point transparency (0–1)
plot_sites_satellite <- function(sites, title,
                                 col      = "CBWF",
                                 zoom     = 7,
                                 roads    = NULL,   # optional sf lines for roads/tracks
                                 places   = NULL,   # optional sf points for place labels
                                 provider = "Esri.WorldImagery",
                                 alpha    = 0.85) {

  sites_buf <- buffer(sites, width = 50000)
  tiles     <- get_tiles(sites_buf, provider = provider, zoom = zoom, crop = TRUE)

  pts   <- as.data.frame(sites, geom = "XY")
  label <- col

  # Crop roads to plot extent (avoid drawing far outside the map)
  roads_layer <- if (!is.null(roads)) {
    ext_sf <- st_as_sf(as.polygons(ext(sites_buf), crs = crs(sites_buf)))
    st_crop(st_transform(roads, st_crs(ext_sf$geometry)), ext_sf)
  }

  p <- ggplot() +
    geom_spatraster_rgb(data = tiles)

  if (!is.null(roads)) {
    p <- p + geom_sf(
      data        = roads_layer,
      colour      = "white",
      linewidth   = 0.5,
      inherit.aes = FALSE
    )
  }

  if (!is.null(places)) {
    places_layer <- st_transform(places, st_crs(tiles))
    coords <- st_coordinates(places_layer)
    places_coords <- places_layer |>
      st_drop_geometry() |>
      mutate(
        x       = coords[, 1],
        y       = coords[, 2],
        # Compute label positions explicitly to avoid geom_text nudge_x/y conflict
        label_x = x + nudge_x,
        label_y = y + nudge_y
      )

    p <- p +
      geom_point(
        data        = places_coords,
        aes(x = x, y = y),
        colour      = "white",
        size        = 1.5,
        inherit.aes = FALSE
      ) +
      geom_text(
        data        = places_coords,
        aes(x = label_x, y = label_y, label = place_name, hjust = hjust),
        colour      = "white",
        size        = 2.5,
        fontface    = "bold",
        inherit.aes = FALSE
      )
  }

  p +
    geom_point(
      data  = pts,
      aes(x = x, y = y, colour = .data[[col]], size = .data[[col]]),
      alpha = alpha
    ) +
    scale_colour_viridis_c(
      name   = label,
      option = "plasma",
      breaks = scales::pretty_breaks(4)
    ) +
    scale_size_continuous(
      name   = label,
      range  = c(2, 5),
      breaks = scales::pretty_breaks(4)
    ) +
    guides(
      colour = guide_colorbar(title = label),
      size   = "none"
    ) +
    labs(title = title, x = NULL, y = NULL) +
    coord_sf(expand = FALSE) +
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

# --- One-time: download SA road/track vectors from OSM and cache locally ------
roads_path <- "./outputs/vector/sa_key_roads.gpkg"

if (!file.exists(roads_path)) {

  # Bounding box covering SA arid zone (these roads span roughly this extent)
  bb <- c(xmin = 128, ymin = -37, xmax = 142, ymax = -25)

  fetch_road <- function(name) {
    result <- opq(bb) |>
      add_osm_feature(key = "name", value = name) |>
      osmdata::osmdata_sf()
    lines <- result$osm_lines
    if (is.null(lines) || nrow(lines) == 0) {
      warning("No OSM lines returned for: ", name)
      return(NULL)
    }
    mutate(lines, road_name = name)
  }

  roads_sf <- bind_rows(
    fetch_road("Stuart Highway"),
    fetch_road("The Outback Highway"),  # Lyndhurst–Marree section of B83
    fetch_road("Strzelecki Track"),
    fetch_road("Oodnadatta Track"),
    fetch_road("Birdsville Track"),
    fetch_road("Olympic Dam Road"),     # Pimba to Roxby Downs
    fetch_road("Borefield Road")        # Roxby Downs to Oodnadatta Track
  ) |>
    select(road_name, geometry) |>
    # Unify the two B83 segments under one label
    mutate(road_name = if_else(road_name == "The Outback Highway", "B83", road_name))

  # B83 (Port Augusta to Marree) — query by road ref rather than name
  b83_raw <- opq(bb) |>
    add_osm_feature(key = "ref", value = "B83") |>
    osmdata::osmdata_sf()
  if (!is.null(b83_raw$osm_lines) && nrow(b83_raw$osm_lines) > 0) {
    roads_sf <- bind_rows(
      roads_sf,
      mutate(b83_raw$osm_lines, road_name = "B83") |> select(road_name, geometry)
    )
  } else {
    warning("No OSM lines returned for ref = B83")
  }

  st_write(roads_sf, roads_path, delete_dsn = TRUE)
  message("Roads saved to ", roads_path)
}

roads_sf <- st_read(roads_path, quiet = TRUE)

# --- One-time: download placename points from OSM and cache locally -----------
places_path <- "./outputs/vector/sa_key_places.gpkg"

# Bounding box covering SA arid zone — used by both roads and places downloads
bb <- c(xmin = 128, ymin = -37, xmax = 142, ymax = -25)

if (!file.exists(places_path)) {

  place_names <- c(
    "Bon Bon", "Glendambo", "Coober Pedy", "Oodnadatta",
    "Marla", "Roxby Downs", "Marree", "Lyndhurst", "Nonning"
  )

  fetch_place <- function(name) {
    result <- opq(bb) |>
      add_osm_feature(key = "name", value = name) |>
      osmdata::osmdata_sf()
    pts <- result$osm_points
    if (!is.null(pts) && "place" %in% names(pts)) {
      pts <- dplyr::filter(pts, !is.na(place))
    }
    if (is.null(pts) || nrow(pts) == 0) {
      warning("No OSM place point returned for: ", name)
      return(NULL)
    }
    # Take the first match and keep only the label and geometry
    pts[1, ] |> mutate(place_name = name) |> select(place_name, geometry)
  }

  places_sf <- bind_rows(lapply(place_names, fetch_place))

  st_write(places_sf, places_path, delete_dsn = TRUE)
  message("Places saved to ", places_path)
}

places_sf <- st_read(places_path, quiet = TRUE)

# Label displacements — edit these vectors to reposition labels.
# Units are metres (~30 km = 30000 m). Order matches places_sf$place_name:
# "Bon Bon", "Glendambo", "Coober Pedy", "Oodnadatta", "Marla", "Roxby Downs", "Marree", "Lyndhurst", "Nonning"
label_nudge <- data.frame(
  place_name = places_sf$place_name,
  nudge_x    = c(-45000,  10000,  10000,  10000,  10000, -50000,  20000,  15000,  10000),
  nudge_y    = c( 18000,  15000,  10000,  10000,  15000, -10000,      0, -10000,  10000),
  hjust      = c(     0,      0,      0,      0,      0,      0,      0,      0,      0)
)

places_sf <- left_join(places_sf, label_nudge, by = "place_name")

# --- CBWF records -------------------------------------------------------------
sites <- vect("./outputs/vector/CBWF_site_locations.gpkg")

# plot
p <- plot_sites_satellite(sites, title = "", col = "CBWF",
                          zoom = 7, roads = roads_sf, places = places_sf)

# --- Inset: Australia outline with study region marker -----------------------
aus <- ozmaps::ozmap_country

# Study region centre from site locations
sites_wgs84  <- project(sites, "EPSG:4326")
study_centre <- crds(sites_wgs84) |>
  colMeans() |>
  (\(m) st_sfc(st_point(m), crs = 4326))()

inset <- ggplot() +
  geom_sf(data = aus, fill = "white", colour = "black", linewidth = 0.4) +
  geom_sf(data = study_centre, colour = "red", size = 2.5) +
  theme_void() +
  panel_border(remove = TRUE)

# Combine: inset placed in top-right corner (adjust x, y, width, height as needed)
p_final <- ggdraw(p) +
  draw_plot(inset, x = 0.6, y = 0.7, width = 0.3, height = 0.27)

p_final

# Save
ggsave("./figures/cbwf_survey_outcomes.png", p_final,
       width = 18, height = 18, units = "cm", dpi = 300)

