# ---------------------------------------------------------------------------
# s2_watersheds_v5a.R
# Clip NDVI / NBR to the study watersheds and plot the burned-vs-reference
# pairs (CA vs AR, CA_TO vs AR_UP), with a DIFFERENT colour scale per index.
#
# Run AFTER s2_rgb_ndvi_nbr_v5a.R in the same session - it reuses load_scene()
# and scene_indices(). (The guard below rebuilds `idx` if it's missing.)
# ---------------------------------------------------------------------------

# install once: install.packages("patchwork")
library(patchwork)    # arranges several independently-scaled plots together

# Rebuild the 2016 NDVI+NBR stack if you started a fresh session.
if (!exists("idx")) {
  yw    <- "2016_summer"
  bands <- load_scene(yw)
  idx   <- scene_indices(bands)
}

# --- Watershed shapefiles (paths straight from your Python config) ---------
# Note the space in "WILDFIRE_DATA_from Pete" - R handles it fine as a string.
ws_root <- "/Users/benthosyy/Desktop/DHSVM-Pete/WILDFIRE_DATA_from Pete/GIS/WS_boundaries"
ws_files <- list(
  CA    = file.path(ws_root, "CA",    "cabr_watershed.shp"),  # Camp Branch  - burned (21% high-sev)
  AR    = file.path(ws_root, "AR",    "arwd_watershed.shp"),  # Arrowwood    - reference for CA
  CA_TO = file.path(ws_root, "CA_TO", "tobr_watershed.shp"),  # Tower Branch - burned (65% high-sev), nested in CA
  AR_UP = file.path(ws_root, "AR_UP", "upar_watershed.shp")   # Upper Arrowwood - reference for CA_TO, nested in AR
)

# Read each polygon with vect() and reproject to the raster's CRS (UTM 17N).
# lapply() runs a function over a list - like a Python list comprehension.
ws <- lapply(ws_files, function(f) project(vect(f), crs(idx)))

# --- Clip the index stack to one watershed ---------------------------------
# crop() trims to the bounding box; mask() then sets pixels outside the
# polygon to NA. So clip_to() returns just the watershed, on transparent.
clip_to <- function(r, poly) mask(crop(r, poly), poly)

# --- One map of a single index, with its OWN colour scale ------------------
# Plotting NDVI and NBR as separate ggplots (not facets) is what lets each
# one carry a different palette/colourbar.
map_index <- function(r, layer, poly, fill_scale, title) {
  ggplot() +
    geom_spatraster(data = r[[layer]]) +                      # the index raster
    geom_spatvector(data = poly, fill = NA,                   # crisp boundary
                    colour = "grey15", linewidth = 0.4) +
    fill_scale +
    labs(title = title, fill = layer) +
    theme_minimal() +
    theme(axis.text = element_text(size = 6),
          plot.title = element_text(size = 10))
}

# --- The two colour scales (THIS is the knob you asked about) --------------
# Swap palette names to taste. RColorBrewer palettes work with distiller;
# viridis options are "viridis"/"magma"/"inferno"/"plasma"/"cividis".
ndvi_scale <- scale_fill_distiller(palette = "RdYlGn", direction = 1,
                                   limits = c(-0.2, 1),
                                   na.value = "transparent")   # red->green = veg
nbr_scale  <- scale_fill_viridis_c(option = "magma",
                                   limits = c(-0.5, 1),
                                   na.value = "transparent")   # clearly != NDVI

# --- Plot a burned / reference pair: NDVI on top, NBR below -----------------
plot_pair <- function(burned, reference) {
  b <- clip_to(idx, ws[[burned]])
  r <- clip_to(idx, ws[[reference]])

  ndvi_row <- map_index(b, "NDVI", ws[[burned]],    ndvi_scale,
                        paste(burned, "- NDVI (burned)")) |
              map_index(r, "NDVI", ws[[reference]], ndvi_scale,
                        paste(reference, "- NDVI (reference)"))

  nbr_row  <- map_index(b, "NBR", ws[[burned]],    nbr_scale,
                        paste(burned, "- NBR (burned)")) |
              map_index(r, "NBR", ws[[reference]], nbr_scale,
                        paste(reference, "- NBR (reference)"))

  (ndvi_row / nbr_row) +
    plot_annotation(
      title = paste0(burned, " vs ", reference, "   -   ", scenes[[yw]]))
}

# CA (burned) vs AR (reference)
plot_pair("CA", "AR")

# CA_TO (burned) vs AR_UP (reference)
plot_pair("CA_TO", "AR_UP")
