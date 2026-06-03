# ---------------------------------------------------------------------------
# s2_paired_rgb_v5a.R
# Reproduce the "paired watersheds x 3 years" RGB grid:
#   columns = (CA + AR together) | CA_TO | AR_UP
#   rows    = 2016 (pre-fire) | 2017 (yr 1) | 2018 (yr 2)
# Burned watershed outlines in RED, unburned controls in WHITE. Scale bar +
# north arrow via ggspatial. Coordinate ticks are removed (that's the fix for
# the overlapping longitude labels, and it's what your target figure does).
#
# Run after s2_rgb_ndvi_nbr_v5a.R (reuses scenes / data_root / stretch_band).
# The burn-severity top row needs the RdNBR raster - see note at the bottom.
# ---------------------------------------------------------------------------
# install once: install.packages(c("patchwork", "ggspatial"))
library(terra); library(tidyterra); library(ggplot2)
library(patchwork); library(ggspatial)

# --- Watershed polygons (rebuild if not already in the session) ------------
if (!exists("ws")) {
  ws_root <- "/Users/benthosyy/Desktop/DHSVM-Pete/WILDFIRE_DATA_from Pete/GIS/WS_boundaries"
  ws_files <- list(
    CA    = file.path(ws_root, "CA",    "cabr_watershed.shp"),
    AR    = file.path(ws_root, "AR",    "arwd_watershed.shp"),
    CA_TO = file.path(ws_root, "CA_TO", "tobr_watershed.shp"),
    AR_UP = file.path(ws_root, "AR_UP", "upar_watershed.shp"))
  ws <- lapply(ws_files, function(f) project(vect(f), "EPSG:32617"))
}

burned_set <- c("CA", "CA_TO")            # these watersheds get the red outline

# Clean panel theme: no coordinate text/ticks (kills the label overlap) but
# keeps the y-axis title so we can use it as a row label on the left column.
panel_theme <- theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(angle = 90, face = "bold", size = 10),
        plot.title   = element_text(hjust = 0.5, face = "bold", size = 11))

# Pad a SpatVector's extent so the watersheds sit in a little context.
pad_ext <- function(v, frac = 0.10) {
  e  <- as.vector(ext(v))                 # named: xmin xmax ymin ymax
  dx <- (e["xmax"] - e["xmin"]) * frac
  dy <- (e["ymax"] - e["ymin"]) * frac
  ext(e["xmin"] - dx, e["xmax"] + dx, e["ymin"] - dy, e["ymax"] + dy)
}

# Stretched RGB for one year, cropped to a column's extent. Default = true
# colour (B04/B03/B02). For a SWIR burn composite (vivid burns), pass
# c("B12","B08","B04") - almost certainly what your magenta 2017/2018 use.
column_rgb <- function(year_window, e, rgb_bands = c("B04","B03","B02")) {
  date <- scenes[[year_window]]
  rd <- function(b) crop(rast(file.path(data_root, year_window,
                                        paste0(date, "_", b, ".tif"))), e)
  r <- rd(rgb_bands[1]); g <- rd(rgb_bands[2]); b <- rd(rgb_bands[3])
  # SWIR bands are 20 m; line them up with the 10 m band if we're mixing.
  if (!all(dim(g)[1:2] == dim(r)[1:2])) g <- resample(g, r)
  if (!all(dim(b)[1:2] == dim(r)[1:2])) b <- resample(b, r)
  s <- c(stretch_band(r), stretch_band(g), stretch_band(b))
  names(s) <- c("r", "g", "b"); s
}

# One panel: RGB + watershed outlines, optional column title / row label /
# scale bar / north arrow.
panel <- function(year_window, polys, e, rgb_bands,
                  col_title = NULL, row_label = NULL,
                  scalebar = FALSE, north = FALSE) {
  s <- column_rgb(year_window, e, rgb_bands)
  g <- ggplot() + geom_spatraster_rgb(data = s, r = 1, g = 2, b = 3,
                                      max_col_value = 1)
  for (nm in names(polys)) {
    outline <- if (nm %in% burned_set) "red" else "white"
    g <- g + geom_spatvector(data = polys[[nm]], fill = NA,
                             colour = outline, linewidth = 0.6)
  }
  g <- g + coord_sf(expand = FALSE) + panel_theme
  if (!is.null(col_title)) g <- g + labs(title = col_title)
  if (!is.null(row_label)) g <- g + ylab(row_label)
  if (scalebar) g <- g + annotation_scale(location = "bl",
                                          height = unit(0.15, "cm"))
  if (north)    g <- g + annotation_north_arrow(location = "tr",
                                                height = unit(0.8, "cm"),
                                                width  = unit(0.8, "cm"))
  g
}

# --- Columns, their combined extents, rows, band combo ---------------------
cols <- list(
  "Camp Branch (CA) + Arrowwood (AR)" = list(CA = ws$CA, AR = ws$AR),
  "Tower Branch (CA_TO)"              = list(CA_TO = ws$CA_TO),
  "Upper Arrowwood (AR_UP)"           = list(AR_UP = ws$AR_UP))
col_ext <- lapply(cols, function(polys) pad_ext(do.call(rbind, unname(polys))))

years     <- names(scenes)                          # 2016 / 2017 / 2018
row_labs  <- c("2016 (pre-fire)", "2017 (yr 1 post-fire)",
               "2018 (yr 2 post-fire)")
rgb_bands <- c("B04", "B03", "B02")    # <- try c("B12","B08","B04") for SWIR

# --- Build every panel; lay out rows = years, columns = groups -------------
plots <- list()
for (yi in seq_along(years)) {
  for (ci in seq_along(cols)) {
    plots[[length(plots) + 1]] <- panel(
      year_window = years[yi],
      polys       = cols[[ci]],
      e           = col_ext[[ci]],
      rgb_bands   = rgb_bands,
      col_title   = if (yi == 1) names(cols)[ci] else NULL,
      row_label   = if (ci == 1) row_labs[yi]    else NULL,
      scalebar    = (yi == length(years)),        # scale bar on the bottom row
      north       = (yi == 1 && ci == 1))         # one north arrow, top-left
  }
}

fig <- wrap_plots(plots, ncol = length(cols)) +
  plot_annotation(
    title    = "Sentinel-2 RGB across paired watersheds",
    subtitle = paste0("v5a: ", paste(unlist(scenes), collapse = " / "),
                      "   |   burned = red outline, control = white"))
fig

# Save (uncomment):
# ggsave("outputs/paired_rgb_v5a.png", fig, width = 12, height = 12, dpi = 200)

# --- To add the burn-severity TOP row (your Python Figure A row 1) ---------
# Send me the RdNBR raster path (your Python auto-discovers
# "RdNBR_20160609_20170726.tif"). I'll add a row that reprojects it onto the
# S2 grid, classifies it with the Caldwell cuts (62 / 181 / 541) into the four
# severity classes, and renders it with the green/yellow/orange/red palette.

