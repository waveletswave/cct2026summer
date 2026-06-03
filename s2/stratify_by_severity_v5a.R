# ---------------------------------------------------------------------------
# stratify_by_severity_v5a.R
# Option 1: bring RdNBR into R, classify burn severity (Caldwell cuts
# 62 / 181 / 541), and put the PER-CLASS trajectory next to the whole-watershed
# (pooled) trajectory - so you can see exactly why they tell different stories.
#
# Self-contained: Session > Restart R first, then run this top to bottom.
# (terra objects don't survive a restart, so we rebuild everything here.)
# ---------------------------------------------------------------------------
# install once: install.packages(c("terra","dplyr","tidyr","ggplot2"))
library(terra); library(dplyr); library(tidyr); library(ggplot2)

# ---- paths & constants ----
data_root  <- "/Users/benthosyy/Desktop/DHSVM-Pete/dhsvm_rs_test/s2_data_v5a"
rdnbr_path <- "/Users/benthosyy/Desktop/DHSVM-Pete/WILDFIRE_DATA_from Pete/GIS/burn_severity/RdNBR_20160609_20170726.tif"
ws_root    <- "/Users/benthosyy/Desktop/DHSVM-Pete/WILDFIRE_DATA_from Pete/GIS/WS_boundaries"

scenes     <- list("2016_summer" = "2016-06-29",
                   "2017_summer" = "2017-07-26",
                   "2018_summer" = "2018-07-21")
NODATA     <- c(0, -9999, -32768)
RDNBR_CUTS <- c(62, 181, 541)               # severity breakpoints (your Python)
class_labs <- c("U-M", "M", "M-H", "H")     # classes 1..4

# ---- NDVI / NBR for one year ----
drop_nodata <- function(r) subst(r, NODATA, NA)
scene_indices <- function(year_window) {
  date <- scenes[[year_window]]
  p <- function(b) file.path(data_root, year_window, paste0(date, "_", b, ".tif"))
  red   <- drop_nodata(rast(p("B04")))
  nir   <- drop_nodata(rast(p("B08")))
  swir2 <- drop_nodata(rast(p("B12")))
  swir2 <- resample(swir2, nir, method = "bilinear")     # 20 m -> 10 m
  ndvi <- (nir - red)   / (nir + red)
  nbr  <- (nir - swir2) / (nir + swir2)
  names(ndvi) <- "NDVI"; names(nbr) <- "NBR"
  c(ndvi, nbr)
}

# ---- watershed polygons ----
ws <- list(
  CA    = project(vect(file.path(ws_root, "CA",    "cabr_watershed.shp")), "EPSG:32617"),
  AR    = project(vect(file.path(ws_root, "AR",    "arwd_watershed.shp")), "EPSG:32617"),
  CA_TO = project(vect(file.path(ws_root, "CA_TO", "tobr_watershed.shp")), "EPSG:32617"),
  AR_UP = project(vect(file.path(ws_root, "AR_UP", "upar_watershed.shp")), "EPSG:32617"))

# ---- NDVI/NBR stacks for each year ----
idx_by_year <- lapply(names(scenes), scene_indices)
names(idx_by_year) <- names(scenes)

# ---- RdNBR -> 4 severity classes ----
rdnbr <- rast(rdnbr_path)
rdnbr <- subst(rdnbr, c(-9999, -32768), NA)        # guard against sentinel values
# left-closed intervals [a, b) to match the Python classifier (right = FALSE)
sev <- classify(rdnbr,
                rcl = matrix(c(-Inf, RDNBR_CUTS[1], 1,
                               RDNBR_CUTS[1], RDNBR_CUTS[2], 2,
                               RDNBR_CUTS[2], RDNBR_CUTS[3], 3,
                               RDNBR_CUTS[3], Inf, 4), ncol = 3, byrow = TRUE),
                right = FALSE)
names(sev) <- "class"

# ---- per-pixel table: NDVI, NBR, severity class, per watershed & year ----
build_tbl <- function(ws_name, year_window) {
  idx   <- idx_by_year[[year_window]]
  sev_a <- project(sev, idx, method = "near")      # align class raster to this grid
  stk   <- c(idx, sev_a)                           # NDVI, NBR, class on one grid
  df    <- terra::extract(stk, ws[[ws_name]], ID = FALSE)
  df$watershed <- ws_name
  df$year      <- sub("_summer", "", year_window)
  df$group     <- if (ws_name %in% c("CA", "CA_TO")) "burned" else "reference"
  df
}

combos <- expand.grid(ws = c("CA", "AR"), yr = names(scenes), stringsAsFactors = FALSE)
tbl <- bind_rows(Map(build_tbl, combos$ws, combos$yr)) |>
  mutate(class = factor(class, levels = 1:4, labels = class_labs)) |>
  filter(!is.na(class), is.finite(NDVI), is.finite(NBR))

# Sanity check: pixel counts & area fractions per class (compare to Python table)
cat("\n--- pixels per class per watershed ---\n")
print(tbl |> count(watershed, class))
cat("\n--- CA area fractions (compare to your Python delta table) ---\n")
print(tbl |> filter(group == "burned") |> count(class) |>
        mutate(frac = round(n / sum(n), 2)))

# =====================================================================
# THE PUNCHLINE: per-class trajectory vs whole-watershed pooled trajectory
# =====================================================================
long_ca <- tbl |> filter(group == "burned") |>
  pivot_longer(c(NDVI, NBR), names_to = "index", values_to = "value")

per_class <- long_ca |>
  group_by(index, class, year) |>
  summarise(median = median(value), .groups = "drop")

pooled <- long_ca |>                                   # ignore class -> the R "whole-watershed" view
  group_by(index, year) |>
  summarise(median = median(value), .groups = "drop") |>
  mutate(class = "Pooled (all pixels)")

ggplot(per_class, aes(year, median, colour = class, group = class)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  geom_line(data = pooled, aes(year, median, group = class),
            colour = "black", linewidth = 1.6, linetype = "dashed") +
  geom_point(data = pooled, aes(year, median), colour = "black", size = 2) +
  facet_wrap(~index, scales = "free_y") +
  scale_colour_manual(values = c("U-M" = "#2ECC71", "M" = "#F1C40F",
                                 "M-H" = "#E67E22", "H" = "#C0392B")) +
  labs(title = "Camp Branch (CA): per-severity-class vs whole-watershed trajectory",
       subtitle = "Dashed black = pooled median of ALL pixels; colours = by burn-severity class",
       x = NULL, y = "median index value", colour = "Severity class") +
  theme_minimal()

# =====================================================================
# Delta table to cross-check against my Python (medians; raw, not vs control)
# =====================================================================
cat("\n--- CA per-class deltas (medians) ---\n")
per_class |>
  pivot_wider(names_from = year, values_from = median) |>
  mutate(d1_2017_2016 = round(`2017` - `2016`, 3),
         d2_2018_2016 = round(`2018` - `2016`, 3)) |>
  arrange(index, class) |>
  print(n = 20)
