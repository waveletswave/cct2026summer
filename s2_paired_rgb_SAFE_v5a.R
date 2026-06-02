# ---------------------------------------------------------------------------
# s2_paired_rgb_SAFE_v5a.R
# A crash-resistant version of the paired RGB figure.
#  - Uses terra's own plotRGB (battle-tested, low memory) instead of the
#    ggplot RGB geom that was likely crashing your session.
#  - SELF-CONTAINED: defines everything it needs, so it does NOT depend on
#    other scripts or on what's already in memory. Just run this one file,
#    top to bottom, in a fresh session.
#
# Only one package needed.
# ---------------------------------------------------------------------------
library(terra)

# --- Everything this script needs (no dependence on other scripts) ---------
data_root <- "/Users/benthosyy/Desktop/DHSVM-Pete/dhsvm_rs_test/s2_data_v5a"
scenes <- list("2016_summer" = "2016-06-29",
               "2017_summer" = "2017-07-26",
               "2018_summer" = "2018-07-21")
CLOUD_DN   <- 5000
burned_set <- c("CA", "CA_TO")          # red outline; the rest get white

# Per-band 2-98% stretch to [0,1], over land pixels only. Has fallbacks so it
# can't error on an empty/clouded clip.
stretch_band <- function(b, lo = 2, hi = 98) {
  v    <- as.vector(values(b))
  pool <- v[is.finite(v) & v > 0 & v < CLOUD_DN]
  if (length(pool) < 10) pool <- v[is.finite(v)]   # fall back to all valid
  if (length(pool) < 2)  return(b * 0)             # nothing usable -> black
  q <- quantile(pool, c(lo, hi) / 100, names = FALSE)
  lo_v <- q[1]; hi_v <- q[2]; if (hi_v <= lo_v) hi_v <- lo_v + 1
  clamp((b - lo_v) / (hi_v - lo_v), 0, 1)
}

# Watershed polygons, rebuilt every run so session state can't bite you.
ws_root <- "/Users/benthosyy/Desktop/DHSVM-Pete/WILDFIRE_DATA_from Pete/GIS/WS_boundaries"
ws <- list(
  CA    = project(vect(file.path(ws_root, "CA",    "cabr_watershed.shp")), "EPSG:32617"),
  AR    = project(vect(file.path(ws_root, "AR",    "arwd_watershed.shp")), "EPSG:32617"),
  CA_TO = project(vect(file.path(ws_root, "CA_TO", "tobr_watershed.shp")), "EPSG:32617"),
  AR_UP = project(vect(file.path(ws_root, "AR_UP", "upar_watershed.shp")), "EPSG:32617"))

pad_ext <- function(v, frac = 0.10) {
  e  <- as.vector(ext(v))
  dx <- (e["xmax"] - e["xmin"]) * frac; dy <- (e["ymax"] - e["ymin"]) * frac
  ext(e["xmin"] - dx, e["xmax"] + dx, e["ymin"] - dy, e["ymax"] + dy)
}

# Read 3 bands, crop, DOWNSAMPLE (safety cap on pixel count), stretch -> RGB.
column_rgb <- function(year_window, e, rgb_bands = c("B04","B03","B02"),
                       max_px = 600) {
  date <- scenes[[year_window]]
  rd <- function(b) crop(rast(file.path(data_root, year_window,
                                        paste0(date, "_", b, ".tif"))), e)
  r <- rd(rgb_bands[1]); g <- rd(rgb_bands[2]); b <- rd(rgb_bands[3])
  if (!all(dim(g)[1:2] == dim(r)[1:2])) g <- resample(g, r)
  if (!all(dim(b)[1:2] == dim(r)[1:2])) b <- resample(b, r)
  s   <- c(r, g, b)
  fac <- max(1, floor(max(dim(s)[1:2]) / max_px))   # keep longest side <= max_px
  if (fac > 1) s <- aggregate(s, fact = fac, fun = "mean", na.rm = TRUE)
  s <- c(stretch_band(s[[1]]), stretch_band(s[[2]]), stretch_band(s[[3]]))
  names(s) <- c("r", "g", "b"); s
}

# The three columns and their combined extents.
cols <- list(
  "Camp Branch (CA) + Arrowwood (AR)" = list(CA = ws$CA, AR = ws$AR),
  "Tower Branch (CA_TO)"              = list(CA_TO = ws$CA_TO),
  "Upper Arrowwood (AR_UP)"           = list(AR_UP = ws$AR_UP))
col_ext <- lapply(cols, function(polys) pad_ext(do.call(rbind, unname(polys))))

# ===========================================================================
# STEP 1  ->  RUN JUST THIS BLOCK FIRST.  One small panel, base graphics.
#             If an image appears in the Plots pane, we're good.
# ===========================================================================
s1 <- column_rgb("2016_summer", col_ext[[1]])
plotRGB(s1, r = 1, g = 2, b = 3, scale = 1)
plot(ws$CA, add = TRUE, border = "red",   lwd = 2)
plot(ws$AR, add = TRUE, border = "white", lwd = 2)

# ===========================================================================
# STEP 2  ->  ONLY if STEP 1 drew an image: the full 3 x 3 grid.
#             Select this whole block and run it.
# ===========================================================================
years    <- names(scenes)
row_labs <- c("2016 (pre-fire)", "2017 (yr 1)", "2018 (yr 2)")

op <- par(mfrow = c(3, 3), oma = c(0, 3, 3, 0), mar = c(0.5, 0.5, 2.5, 0.5))
for (yi in seq_along(years)) {
  for (ci in seq_along(cols)) {
    s <- column_rgb(years[yi], col_ext[[ci]])
    plotRGB(s, r = 1, g = 2, b = 3, scale = 1, axes = FALSE,
            main = if (yi == 1) names(cols)[ci] else "")
    if (ci == 1) mtext(row_labs[yi], side = 2, line = 1, cex = 0.8, font = 2)
    for (nm in names(cols[[ci]])) {
      plot(cols[[ci]][[nm]], add = TRUE, lwd = 2,
           border = if (nm %in% burned_set) "red" else "white")
    }
  }
}
par(op)   # restore normal single-plot layout

# To switch to the SWIR burn composite (vivid magenta burns), change the band
# argument in column_rgb() calls to c("B12","B08","B04").
