# ===========================================================================
# stratify_by_severity_v5a_v5.R   ===  MULTI-DATASET (v5a + v5b_A + v5b_B)  ===
#
# SCL cloud-masked NDVI/NBR + RdNBR severity classes, run for three image
# selections that differ only in acquisition DATE (date auto-read from filenames).
#
# Figures (saved to my Desktop):
#   FIGURE 1 - combined per-class trajectory lines, faceted index x dataset
#              (dates listed in the caption).
#   FIGURE 2 - per-dataset violin + box + median line; the three S2 acquisition
#              dates are shown ONCE, in the subtitle (not on the x-axis).
#
# v4 changes vs v3: (a) FIGURE 2 row-selection uses base-R subsetting (avoids a
# dplyr multi-condition filter() quirk); (b) dates moved from x-axis to subtitle.
#
# Session > Restart R first, then: open this file, Cmd+A (select all), Source.
# ===========================================================================
# install once: install.packages(c("terra","dplyr","tidyr","ggplot2"))
library(terra); library(dplyr); library(tidyr); library(ggplot2)

# ---- config ----
base_dir   <- "/Users/benthosyy/Desktop/DHSVM-Pete/dhsvm_rs_test"
datasets   <- c("v5a", "v5b_A", "v5b_B")
rdnbr_path <- "/Users/benthosyy/Desktop/DHSVM-Pete/WILDFIRE_DATA_from Pete/GIS/burn_severity/RdNBR_20160609_20170726.tif"
ws_root    <- "/Users/benthosyy/Desktop/DHSVM-Pete/WILDFIRE_DATA_from Pete/GIS/WS_boundaries"
out_dir    <- file.path(path.expand("~"), "Desktop")

year_windows <- c("2016_summer" = "2016", "2017_summer" = "2017", "2018_summer" = "2018")
NODATA     <- c(0, -9999, -32768)
RDNBR_CUTS <- c(62, 181, 541)
class_labs <- c("U-M", "M", "M-H", "H")
class_cols <- c("U-M" = "#2ECC71", "M" = "#F1C40F", "M-H" = "#E67E22", "H" = "#C0392B")
SCL_KEEP   <- c(4, 5, 7)

drop_nodata <- function(r) subst(r, NODATA, NA)

scene_date <- function(ds, yw) {
  fld <- file.path(base_dir, paste0("s2_data_", ds), yw)
  f   <- list.files(fld, pattern = "_B04\\.tif$")
  if (length(f) == 0) stop("No *_B04.tif found in: ", fld)
  sub("_B04\\.tif$", "", f[1])
}

scene_indices <- function(ds, yw) {
  date <- scene_date(ds, yw)
  fld  <- file.path(base_dir, paste0("s2_data_", ds), yw)
  p <- function(b) file.path(fld, paste0(date, "_", b, ".tif"))
  red   <- drop_nodata(rast(p("B04")))
  nir   <- drop_nodata(rast(p("B08")))
  swir2 <- drop_nodata(rast(p("B12")))
  swir2 <- resample(swir2, nir, method = "bilinear")
  good  <- ifel(resample(rast(p("SCL")), nir, method = "near") %in% SCL_KEEP, 1, NA)
  ndvi <- (nir - red)   / (nir + red)
  nbr  <- (nir - swir2) / (nir + swir2)
  ndvi <- ndvi * good; nbr <- nbr * good
  names(ndvi) <- "NDVI"; names(nbr) <- "NBR"
  c(ndvi, nbr)
}

ws <- list(
  CA    = project(vect(file.path(ws_root, "CA",    "cabr_watershed.shp")), "EPSG:32617"),
  AR    = project(vect(file.path(ws_root, "AR",    "arwd_watershed.shp")), "EPSG:32617"),
  CA_TO = project(vect(file.path(ws_root, "CA_TO", "tobr_watershed.shp")), "EPSG:32617"),
  AR_UP = project(vect(file.path(ws_root, "AR_UP", "upar_watershed.shp")), "EPSG:32617"))

rdnbr <- subst(rast(rdnbr_path), c(-9999, -32768), NA)
sev <- classify(rdnbr,
                rcl = matrix(c(-Inf, RDNBR_CUTS[1], 1, RDNBR_CUTS[1], RDNBR_CUTS[2], 2,
                               RDNBR_CUTS[2], RDNBR_CUTS[3], 3, RDNBR_CUTS[3], Inf, 4),
                             ncol = 3, byrow = TRUE), right = FALSE)
names(sev) <- "class"

build_tbl_ds <- function(ds) {
  idx_by_year <- lapply(names(year_windows), function(yw) scene_indices(ds, yw))
  names(idx_by_year) <- names(year_windows)
  one <- function(ws_name, yw) {
    idx   <- idx_by_year[[yw]]
    sev_a <- project(sev, idx, method = "near")
    df    <- terra::extract(c(idx, sev_a), ws[[ws_name]], ID = FALSE)
    df$dataset <- ds; df$watershed <- ws_name; df$year <- year_windows[[yw]]
    df$group <- if (ws_name %in% c("CA", "CA_TO")) "burned" else "reference"
    df
  }
  combos <- expand.grid(ws = c("CA", "AR"), yw = names(year_windows), stringsAsFactors = FALSE)
  bind_rows(Map(one, combos$ws, combos$yw))
}

tbl_all <- bind_rows(lapply(datasets, build_tbl_ds)) |>
  mutate(class = factor(class, levels = 1:4, labels = class_labs),
         year  = factor(year,  levels = c("2016", "2017", "2018"))) |>
  filter(!is.na(class), is.finite(NDVI), is.finite(NBR))

dates_df <- expand.grid(dataset = datasets, yw = names(year_windows), stringsAsFactors = FALSE)
dates_df$year <- year_windows[dates_df$yw]
dates_df$date <- mapply(scene_date, dates_df$dataset, dates_df$yw)

cat("\n--- S2 acquisition date used per dataset/year ---\n")
print(dates_df[order(dates_df$dataset, dates_df$year), c("dataset", "year", "date")], row.names = FALSE)
cat("\n--- CA U-M low-NDVI fraction by dataset & year ---\n")
print(tbl_all |> filter(watershed == "CA", class == "U-M") |>
        group_by(dataset, year) |>
        summarise(n = n(), pct_low = round(100 * mean(NDVI < 0.5), 1),
                  med_ndvi = round(median(NDVI), 3), .groups = "drop"))

# ===========================================================================
# FIGURE 1 - combined per-class vs pooled trajectory, faceted index x dataset
# ===========================================================================
long_all <- tbl_all |> filter(group == "burned") |>
  pivot_longer(c(NDVI, NBR), names_to = "index", values_to = "value")
per_class_all <- long_all |> group_by(dataset, index, class, year) |>
  summarise(median = median(value), .groups = "drop")
pooled_all <- long_all |> group_by(dataset, index, year) |>
  summarise(median = median(value), .groups = "drop") |> mutate(class = "Pooled (all pixels)")

date_caption <- paste0("S2 image dates  —  ",
  paste(vapply(datasets, function(ds) {
    d <- dates_df[dates_df$dataset == ds, ]; d <- d[order(d$year), ]
    paste0(ds, ": ", paste(d$date, collapse = " / "))
  }, character(1)), collapse = "      "))

p1 <- ggplot(per_class_all, aes(year, median, colour = class, group = class)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.8) +
  geom_line(data = pooled_all, aes(year, median, group = class),
            colour = "black", linewidth = 1.3, linetype = "dashed") +
  geom_point(data = pooled_all, aes(year, median), colour = "black", size = 1.8) +
  facet_grid(index ~ dataset, scales = "free_y") +
  scale_colour_manual(values = class_cols) +
  labs(title = "Camp Branch (CA): per-severity-class vs whole-watershed, by image date",
       subtitle = "Cloud-masked (SCL). Dashed black = pooled median of ALL pixels.",
       caption = date_caption, x = NULL, y = "median index value", colour = "Severity class") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0, size = 9, colour = "grey30"))
print(p1)
f1 <- file.path(out_dir, "compare_imagedates_trajectory.png")
ggsave(f1, p1, width = 13, height = 6.5, dpi = 200)

cat("\n--- CA per-class deltas (medians) by dataset ---\n")
per_class_all |> pivot_wider(names_from = year, values_from = median) |>
  mutate(d1 = round(`2017` - `2016`, 3), d2 = round(`2018` - `2016`, 3)) |>
  arrange(dataset, index, class) |> print(n = 30)

# ===========================================================================
# FIGURE 2 - per-dataset violin; the three S2 dates shown ONCE in the subtitle.
#            (Row-selection uses base-R subsetting to avoid the dplyr quirk.)
# ===========================================================================
bg_nbr  <- data.frame(index = "NBR")
bg_ndvi <- data.frame(index = "NDVI")

saved <- c(f1)
for (ds in datasets) {
  sub_tbl <- tbl_all[tbl_all$dataset == ds & tbl_all$watershed == "CA", ]
  long_df <- pivot_longer(sub_tbl, c(NDVI, NBR), names_to = "index", values_to = "value")
  long_df <- long_df[is.finite(long_df$value), ]

  d <- dates_df[dates_df$dataset == ds, ]; d <- d[order(d$year), ]
  date_line <- paste0("S2 image dates    ", paste0(d$year, ": ", d$date, collapse = "      "))

  pV <- ggplot(long_df, aes(x = year, y = value)) +
    geom_rect(data = bg_nbr,  fill = "#E67E22", alpha = 0.08,
              xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, inherit.aes = FALSE) +
    geom_rect(data = bg_ndvi, fill = "#27AE60", alpha = 0.08,
              xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, inherit.aes = FALSE) +
    geom_violin(aes(fill = class), colour = NA, alpha = 0.55, scale = "width", width = 0.9) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white",
                 colour = "grey30", linewidth = 0.3) +
    stat_summary(aes(group = 1), fun = median, geom = "line", colour = "grey15", linewidth = 0.5) +
    stat_summary(fun = median, geom = "point", colour = "grey15", size = 1.4) +
    facet_grid(index ~ class, scales = "free_y") +
    scale_fill_manual(values = class_cols, guide = "none") +
    labs(title = paste0("Camp Branch (CA) — ", ds, ":  NDVI / NBR distribution"),
         subtitle = paste0(date_line,
                           "\nViolin = distribution; box = quartiles; line = median trajectory"),
         x = NULL, y = "Index value") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
          plot.subtitle = element_text(size = 10, colour = "grey25"))
  print(pV)
  fV <- file.path(out_dir, paste0("CA_", ds, "_distributions.png"))
  ggsave(fV, pV, width = 9, height = 9, dpi = 200)
  saved <- c(saved, fV)
}
message("Saved figures:"); for (s in saved) message("  -> ", s)
