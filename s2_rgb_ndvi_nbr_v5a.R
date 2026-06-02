# ---------------------------------------------------------------------------
# s2_rgb_ndvi_nbr_v5a.R
# R port of the Sentinel-2 RGB + NDVI + NBR step from the v5a pipeline.
# CCT Computational Fellowship 2026 - Yiyun Song
#
# Faithful to the Python (04D / 06): same v5a folder layout, same scene
# dates, same band roles, same index formulas, the same 20 m -> 10 m
# resample for B12, and the same cloud-aware percentile stretch for the RGB.
#
# Scope: one growing-season scene -> true-colour RGB, NDVI, NBR, and a clean
# comparison figure. (RdNBR burn classes, watershed masks, and the NDVI->LAI
# Beer-Lambert step from the Python pipeline are left for later.)
# ---------------------------------------------------------------------------

# install once:
# install.packages(c("terra", "tidyterra", "ggplot2"))
library(terra)       # the raster engine
library(tidyterra)   # lets ggplot2 draw rasters directly
library(ggplot2)

# --- CONFIG ----------------------------------------------------------------
# Your data root (from the path you gave me):
data_root <- "/Users/benthosyy/Desktop/DHSVM-Pete/dhsvm_rs_test/s2_data_v5a"

# v5a scenes, exactly as in the Python config:
scenes <- list(
  "2016_summer" = "2016-06-29",   # pre-fire
  "2017_summer" = "2017-07-26",   # post-fire year 1
  "2018_summer" = "2018-07-21"    # post-fire year 2
)

# Your pipeline's nodata sentinels and the cloud cut-off for the RGB stretch:
NODATA   <- c(0, -9999, -32768)
CLOUD_DN <- 5000

# --- small helpers ---------------------------------------------------------
# In R, functions are just values: no `def`, and the last expression is
# returned automatically (no `return` needed).

# Set the pipeline's nodata values to NA. terra::subst swaps listed values.
drop_nodata <- function(r) subst(r, NODATA, NA)

# Per-band 2-98% stretch to [0, 1], using only land pixels (0 < DN < CLOUD_DN)
# so clouds elsewhere in the tile don't bias the percentiles and darken the
# scene. Mirrors stretch_band() / make_rgb() in the Python.
stretch_band <- function(b, lo_pct = 2, hi_pct = 98) {
  v    <- as.vector(values(b))
  pool <- v[is.finite(v) & v > 0 & v < CLOUD_DN]
  qs   <- quantile(pool, probs = c(lo_pct, hi_pct) / 100, names = FALSE)
  lo   <- qs[1]; hi <- qs[2]
  if (hi <= lo) hi <- lo + 1
  clamp((b - lo) / (hi - lo), 0, 1)        # terra::clamp keeps values in [0,1]
}

# Load one scene's bands; put B12 (20 m) onto the 10 m grid.
load_scene <- function(year_window) {
  date <- scenes[[year_window]]
  p <- function(band)
    file.path(data_root, year_window, paste0(date, "_", band, ".tif"))

  b02 <- rast(p("B02"))    # blue  10 m
  b03 <- rast(p("B03"))    # green 10 m
  b04 <- rast(p("B04"))    # red   10 m
  b08 <- rast(p("B08"))    # NIR   10 m
  b12 <- rast(p("B12"))    # SWIR2 20 m
  b12 <- resample(b12, b08, method = "bilinear")  # 20 m -> 10 m, like the Python

  list(blue = b02, green = b03, red = b04, nir = b08, swir2 = b12)
}

# NDVI and NBR for a loaded scene (same formulas as compute_index()).
scene_indices <- function(bands) {
  red   <- drop_nodata(bands$red)
  nir   <- drop_nodata(bands$nir)
  swir2 <- drop_nodata(bands$swir2)

  ndvi <- (nir - red)   / (nir + red)      # (B08 - B04) / (B08 + B04)
  nbr  <- (nir - swir2) / (nir + swir2)    # (B08 - B12) / (B08 + B12)
  names(ndvi) <- "NDVI"
  names(nbr)  <- "NBR"
  c(ndvi, nbr)                             # a 2-layer SpatRaster
}

# ===========================================================================
# 1. Pick a scene (start with pre-fire 2016)
# ===========================================================================
yw    <- "2016_summer"
bands <- load_scene(yw)
print(bands$red)            # sanity check: extent, CRS, resolution, min/max

# ===========================================================================
# 2. True-colour RGB (B04 / B03 / B02)
# ===========================================================================
rgb <- c(stretch_band(bands$red),
         stretch_band(bands$green),
         stretch_band(bands$blue))
names(rgb) <- c("red", "green", "blue")

plotRGB(rgb, r = 1, g = 2, b = 3, scale = 1,
        main = paste("Sentinel-2 true colour -", scenes[[yw]]))
# Quick alternative (terra's own stretch, no cloud masking):
#   stack <- c(bands$red, bands$green, bands$blue)
#   plotRGB(stack, stretch = "lin")

# ===========================================================================
# 3. NDVI and NBR
# ===========================================================================
idx <- scene_indices(bands)

# A clean side-by-side map: tidyterra plugs the raster straight into ggplot,
# facet_wrap draws both indices at once, and they share one colour scale.
ggplot() +
  geom_spatraster(data = idx) +
  facet_wrap(~lyr) +
  scale_fill_viridis_c(na.value = "transparent", limits = c(-1, 1)) +
  labs(title = paste("NDVI and NBR -", scenes[[yw]]), fill = "index") +
  theme_minimal()

# ===========================================================================
# 4. (optional) Save outputs
# ===========================================================================
# dir.create("outputs", showWarnings = FALSE)
# writeRaster(idx, file.path("outputs", paste0("ndvi_nbr_", yw, ".tif")),
#             overwrite = TRUE)
# ggsave(file.path("outputs", paste0("ndvi_nbr_", yw, ".png")),
#        width = 9, height = 4, dpi = 200)

# ===========================================================================
# 5. (optional - the part that connects to your burn work)
#    A first dNBR: NBR before the fire minus NBR one year after. This is the
#    raw cousin of the RdNBR your Python pipeline classifies into severity.
# ===========================================================================
# pre  <- scene_indices(load_scene("2016_summer"))[["NBR"]]
# post <- scene_indices(load_scene("2017_summer"))[["NBR"]]
# dnbr <- pre - post                    # high where canopy/char changed most
# names(dnbr) <- "dNBR_2016_2017"
#
# ggplot() +
#   geom_spatraster(data = dnbr) +
#   scale_fill_viridis_c(option = "magma", na.value = "transparent") +
#   labs(title = "dNBR 2016 -> 2017 (raw burn signal)", fill = "dNBR") +
#   theme_minimal()
