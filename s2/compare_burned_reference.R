# ---------------------------------------------------------------------------
# compare_burned_reference.R
# Where R earns its keep for YOUR project. Pull the index values out of the
# rasters, tidy them, and compare burned vs reference across years - in a
# handful of lines of dplyr + ggplot. (Your Python 06 script does the same
# kind of thing in ~400 lines of numpy masking + matplotlib loops.)
#
# Run after s2_rgb_ndvi_nbr_v5a.R (reuses load_scene / scene_indices / scenes)
# and s2_watersheds_v5a.R (reuses ws).
# ---------------------------------------------------------------------------
# install once: install.packages(c("dplyr", "tidyr"))
library(terra); library(dplyr); library(tidyr); library(ggplot2)

# Compute NDVI+NBR once per year (instead of reloading inside the loop).
years       <- names(scenes)
idx_by_year <- lapply(years, function(y) scene_indices(load_scene(y)))
names(idx_by_year) <- years

# For one watershed + year: every pixel's NDVI and NBR as a tidy table.
# terra::extract() pulls the raster values under a polygon into a data.frame -
# one row per pixel, one column per index. This single call replaces the
# manual geometry_mask + boolean-index dance in the Python.
extract_ws_year <- function(ws_name, year_window) {
  # terra:: is required here because tidyr also defines extract() and, being
  # loaded later, masks terra's. The package::function form bypasses that.
  vals <- terra::extract(idx_by_year[[year_window]], ws[[ws_name]], ID = FALSE)
  vals$watershed <- ws_name
  vals$year      <- sub("_summer", "", year_window)
  vals$group     <- if (ws_name %in% c("CA", "CA_TO")) "burned" else "reference"
  vals
}

# Build ONE tidy table for the CA / AR pair across all three years.
combos <- expand.grid(ws = c("CA", "AR"), yr = years, stringsAsFactors = FALSE)
tbl <- bind_rows(Map(extract_ws_year, combos$ws, combos$yr))

# Long form, then plot. This is the part that stays short no matter how many
# watersheds / years / indices you add - the grammar does the bookkeeping.
tbl_long <- tbl |>
  pivot_longer(c(NDVI, NBR), names_to = "index", values_to = "value") |>
  filter(is.finite(value))

# Burned vs reference, by year, for both indices, in one declarative figure.
ggplot(tbl_long, aes(x = year, y = value, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85) +
  facet_wrap(~index, scales = "free_y") +
  scale_fill_manual(values = c(burned = "#C0392B", reference = "#2ECC71")) +
  labs(title = "Camp Branch (burned) vs Arrowwood (reference)",
       subtitle = "Per-pixel NDVI and NBR distributions, by year",
       x = NULL, y = "index value", fill = NULL) +
  theme_minimal()

# A paper-ready numeric summary in four lines: median per index/group/year,
# and the burned-minus-reference gap (the quantity your project is chasing).
tbl_long |>
  group_by(index, group, year) |>
  summarise(median = median(value), .groups = "drop") |>
  pivot_wider(names_from = group, values_from = median) |>
  mutate(burned_minus_reference = burned - reference)
