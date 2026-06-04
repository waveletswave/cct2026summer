# ===========================================================================
# stratify_by_severity_v5a_v6.R  ===  MULTI-DATASET + SELECTABLE WATERSHED  ===
#
# Self-contained, all-in-one post-fire severity analysis. SCL cloud-masked
# NDVI/NBR + RdNBR burn-severity classes, run across the three Sentinel-2 image
# dates (v5a / v5b_A / v5b_B; dates auto-read from filenames). The watershed is
# a single parameter `focus_ws` (CA, CA_TO, AR, or AR_UP) - change one line.
# Output filenames include the watershed, so CA_TO figures do NOT overwrite CA.
#
# Figures (all printed and saved to ~/Desktop):
#   FIGURE 1 - per-severity-class vs whole-watershed (pooled) trajectory,
#              faceted index x dataset.
#   FIGURE 2 - per-dataset violin + box + median line (YEAR on x, faceted by
#              class); the three S2 dates shown once in the subtitle.
#   FIGURE 3 - boxplot with burn-severity CLASS on x and the three YEARS as
#              colours, side by side (the complementary view to FIGURE 2),
#              faceted index x dataset.
#
# Suggested use: run once with focus_ws = "CA", once with focus_ws = "CA_TO",
# then compare (pooled trajectory should differ; per-class should match).
#
# Session > Restart R first, then: open this file, Cmd+A (select all), Source.
#
# ---------------------------------------------------------------------------
# NOTES - relationship to the other scripts in s2/
#   * This file SUPERSEDES stratify_by_severity_v5a.R, _v1, _v2, _v3, and
#     severity_violin_boxplot_v5a.R: it reproduces all of their figures and is
#     fully self-contained. It builds its own per-pixel table internally and
#     does NOT depend on any object created by another script. (The old
#     severity_violin_boxplot_v5a.R needed a `tbl` object from the now-deleted
#     stratify_by_severity_v5a.R, so it could not run on its own.)
#   * stratify_by_severity_v5a_v5.R is kept on purpose: it is the only script
#     that extracts the BURNED watershed AND its unburned REFERENCE (CA + AR)
#     in the SAME run, tagging each pixel burned/reference. That is the starter
#     scaffolding for the planned "stratified + reference-differenced" analysis.
#     To get the reference here instead, set focus_ws = "AR" (or "AR_UP") and
#     run separately, or extend build_tbl_ds() to loop over a (burned,
#     reference) pair. Once that analysis is folded into this file, v5 can go.
# ---------------------------------------------------------------------------
# ===========================================================================
# install once: install.packages(c("terra","dplyr","tidyr","ggplot2"))
library(terra); library(dplyr); library(tidyr); library(ggplot2)

# ---- config ----
base_dir   <- "/Users/benthosyy/Desktop/DHSVM-Pete/dhsvm_rs_test"
datasets   <- c("v5a", "v5b_A", "v5b_B")
rdnbr_path <- "/Users/benthosyy/Desktop/DHSVM-Pete/WILDFIRE_DATA_from Pete/GIS/burn_severity/RdNBR_20160609_20170726.tif"
ws_root    <- "/Users/benthosyy/Desktop/DHSVM-Pete/WILDFIRE_DATA_from Pete/GIS/WS_boundaries"
out_dir    <- file.path(path.expand("~"), "Desktop")

# >>> the watershed to analyse: "CA", "CA_TO", "AR", or "AR_UP" <<<
focus_ws   <- "CA_TO"
ws_names   <- c(CA = "Camp Branch", CA_TO = "Tower Branch",
                AR = "Arrowwood", AR_UP = "Upper Arrowwood")
focus_name <- ws_names[[focus_ws]]

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

# ---- watershed polygons (all defined; only focus_ws is used) ---------------
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

# ---- per-pixel table for ONE dataset, restricted to focus_ws ---------------
build_tbl_ds <- function(ds) {
  idx_by_year <- lapply(names(year_windows), function(yw) scene_indices(ds, yw))
  names(idx_by_year) <- names(year_windows)
  one <- function(yw) {
    idx   <- idx_by_year[[yw]]
    sev_a <- project(sev, idx, method = "near")
    df    <- terra::extract(c(idx, sev_a), ws[[focus_ws]], ID = FALSE)
    df$dataset <- ds; df$watershed <- focus_ws; df$year <- year_windows[[yw]]
    df
  }
  bind_rows(lapply(names(year_windows), one))
}

tbl_all <- bind_rows(lapply(datasets, build_tbl_ds)) |>
  mutate(class = factor(class, levels = 1:4, labels = class_labs),
         year  = factor(year,  levels = c("2016", "2017", "2018"))) |>
  filter(!is.na(class), is.finite(NDVI), is.finite(NBR))

dates_df <- expand.grid(dataset = datasets, yw = names(year_windows), stringsAsFactors = FALSE)
dates_df$year <- year_windows[dates_df$yw]
dates_df$date <- mapply(scene_date, dates_df$dataset, dates_df$yw)

cat("\n--- watershed:", focus_name, "(", focus_ws, ") ---\n")
cat("\n--- pixels per class per dataset (after cloud mask) ---\n")
print(tbl_all |> count(dataset, class))
cat("\n--- low-NDVI fraction by dataset & year (U-M class; cloud check) ---\n")
print(tbl_all |> filter(class == "U-M") |>
        group_by(dataset, year) |>
        summarise(n = n(), pct_low = round(100 * mean(NDVI < 0.5), 1),
                  med_ndvi = round(median(NDVI), 3), .groups = "drop"))

# ===========================================================================
# FIGURE 1 - per-class vs pooled trajectory, faceted index x dataset
# ===========================================================================
long_all <- tbl_all |> pivot_longer(c(NDVI, NBR), names_to = "index", values_to = "value")
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
  labs(title = paste0(focus_name, " (", focus_ws,
                      "): per-severity-class vs whole-watershed, by image date"),
       subtitle = "Cloud-masked (SCL). Dashed black = pooled median of ALL pixels.",
       caption = date_caption, x = NULL, y = "median index value", colour = "Severity class") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0, size = 9, colour = "grey30"))
print(p1)
f1 <- file.path(out_dir, paste0("compare_imagedates_trajectory_", focus_ws, ".png"))
ggsave(f1, p1, width = 13, height = 6.5, dpi = 200)

cat("\n--- per-class deltas (medians) by dataset ---\n")
per_class_all |> pivot_wider(names_from = year, values_from = median) |>
  mutate(d1 = round(`2017` - `2016`, 3), d2 = round(`2018` - `2016`, 3)) |>
  arrange(dataset, index, class) |> print(n = 30)

# ===========================================================================
# FIGURE 2 - per-dataset violin; the three S2 dates shown ONCE in the subtitle.
# ===========================================================================
bg_nbr  <- data.frame(index = "NBR")
bg_ndvi <- data.frame(index = "NDVI")

saved <- c(f1)
for (ds in datasets) {
  sub_tbl <- tbl_all[tbl_all$dataset == ds & tbl_all$watershed == focus_ws, ]
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
    labs(title = paste0(focus_name, " (", focus_ws, ") — ", ds, ":  NDVI / NBR distribution"),
         subtitle = paste0(date_line,
                           "\nViolin = distribution; box = quartiles; line = median trajectory"),
         x = NULL, y = "Index value") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
          plot.subtitle = element_text(size = 10, colour = "grey25"))
  print(pV)
  fV <- file.path(out_dir, paste0(focus_ws, "_", ds, "_distributions.png"))
  ggsave(fV, pV, width = 9, height = 9, dpi = 200)
  saved <- c(saved, fV)
}
# ===========================================================================
# FIGURE 3 - boxplot: burn-severity CLASS on x, the three YEARS as colours,
#            side by side. Complementary view to FIGURE 2 (which puts year on
#            x). One panel per index x dataset; median trend line links the
#            medians within each class.
# ===========================================================================
yr_off  <- c("2016" = -0.26, "2017" = 0, "2018" = 0.26)
box_tbl <- long_all |>
  mutate(class_i = as.integer(class),
         xpos    = class_i + yr_off[as.character(year)])
med_box <- box_tbl |>
  group_by(dataset, index, class, class_i, year) |>
  summarise(median = median(value), xpos = first(xpos), .groups = "drop")

p3 <- ggplot(box_tbl, aes(xpos, value, group = interaction(class, year), fill = year)) +
  geom_boxplot(width = 0.24, outlier.shape = NA, colour = "grey30", linewidth = 0.3) +
  geom_line(data = med_box, aes(xpos, median, group = class),
            inherit.aes = FALSE, colour = "grey20", linewidth = 0.5) +
  geom_point(data = med_box, aes(xpos, median, fill = year),
             inherit.aes = FALSE, shape = 21, size = 2, colour = "grey20") +
  facet_grid(index ~ dataset, scales = "free_y") +
  scale_x_continuous(breaks = 1:4, labels = class_labs) +
  scale_fill_viridis_d(option = "D", end = 0.92) +
  labs(title = paste0(focus_name, " (", focus_ws,
                      "): NDVI / NBR by burn-severity class and year"),
       subtitle = "Boxes = per-pixel distribution; line = median trajectory within each class.",
       caption = date_caption,
       x = "Burn-severity class", y = "Index value", fill = "Year") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "top",
        strip.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0, size = 9, colour = "grey30"))
print(p3)
f3 <- file.path(out_dir, paste0("severity_boxplot_byclass_", focus_ws, ".png"))
ggsave(f3, p3, width = 13, height = 6.5, dpi = 200)
saved <- c(saved, f3)

message("Saved figures for ", focus_name, " (", focus_ws, "):")
for (s in saved) message("  -> ", s)
