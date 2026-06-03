# ---------------------------------------------------------------------------
# boxplot_by_severity_v5a.R
# Python-style boxplots: per-pixel NDVI/NBR distribution for each burn-severity
# class and year, WITH a median trend line overlaid (boxes = distribution,
# line = the 2016->2018 median trajectory).
#
# Run AFTER stratify_by_severity_v5a.R in the SAME session (it uses `tbl`).
# ---------------------------------------------------------------------------
# install once: install.packages(c("dplyr","tidyr","ggplot2"))
library(dplyr); library(tidyr); library(ggplot2)

if (!exists("tbl"))
  stop("Object `tbl` not found - run stratify_by_severity_v5a.R first, same session.")

# --- long form for CA (burned). Change "burned"->"reference" / "CA"->"AR" -----
plot_ws    <- "CA"                       # which watershed to plot
plot_group <- if (plot_ws %in% c("CA","CA_TO")) "burned" else "reference"

long_df <- tbl |>
  filter(watershed == plot_ws) |>
  pivot_longer(c(NDVI, NBR), names_to = "index", values_to = "value") |>
  filter(is.finite(value)) |>
  mutate(year = factor(year, levels = c("2016", "2017", "2018")))

# --- explicit numeric x so the trend line connects cleanly -------------------
# Each box sits at (class position) + (small year offset); the line then links
# the three year-medians WITHIN each class (a short trajectory per class).
yr_off <- c("2016" = -0.25, "2017" = 0, "2018" = 0.25)
long_df <- long_df |>
  mutate(class_i = as.integer(class),
         xpos    = class_i + yr_off[as.character(year)])

med <- long_df |>                                  # per-class-per-year medians
  group_by(index, class, class_i, year) |>
  summarise(median = median(value), xpos = first(xpos), .groups = "drop")

class_cols <- c("U-M" = "#2ECC71", "M" = "#F1C40F", "M-H" = "#E67E22", "H" = "#C0392B")

# =====================================================================
# Boxplots (colour = class, transparency = year) + median trend line
# =====================================================================
ggplot(long_df, aes(x = xpos, y = value,
                    group = interaction(class, year),
                    fill = class, alpha = year)) +
  geom_boxplot(width = 0.22, outlier.shape = NA, colour = "grey25", linewidth = 0.3) +
  geom_line(data = med, aes(x = xpos, y = median, group = class),
            inherit.aes = FALSE, colour = "grey15", linewidth = 0.5) +
  geom_point(data = med, aes(x = xpos, y = median),
             inherit.aes = FALSE, colour = "grey15", size = 1.3) +
  facet_wrap(~ index, scales = "free_y") +
  scale_x_continuous(breaks = 1:4, labels = levels(long_df$class)) +
  scale_fill_manual(values = class_cols, guide = "none") +
  scale_alpha_manual(values = c("2016" = 0.4, "2017" = 0.7, "2018" = 1.0)) +
  labs(title = paste0(plot_ws, ": NDVI / NBR by burn-severity class and year"),
       subtitle = "Boxes = per-pixel distribution; line = median trajectory (2016 -> 2018)",
       x = "Burn-severity class", y = "Index value", alpha = "Year") +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank())

# =====================================================================
# ALTERNATIVE LAYOUT (uncomment): one panel per class, x = year.
# Often the cleanest way to read each class's trajectory.
# =====================================================================
# long_all <- tbl |>
#   filter(group == "burned") |>
#   pivot_longer(c(NDVI, NBR), names_to = "index", values_to = "value") |>
#   filter(is.finite(value)) |>
#   mutate(year = factor(year, levels = c("2016","2017","2018")))
#
# ggplot(long_all, aes(x = year, y = value)) +
#   geom_boxplot(aes(fill = class), outlier.shape = NA, width = 0.6,
#                colour = "grey25", linewidth = 0.3) +
#   stat_summary(aes(group = 1), fun = median, geom = "line", colour = "grey15") +
#   stat_summary(fun = median, geom = "point", colour = "grey15", size = 1.3) +
#   facet_grid(index ~ class, scales = "free_y") +
#   scale_fill_manual(values = class_cols, guide = "none") +
#   labs(x = NULL, y = "Index value",
#        title = "Camp Branch (CA): per-class trajectory (one panel per class)") +
#   theme_minimal()

# =====================================================================
# To show BURNED + CONTROL together (like your Python multi-watershed figure),
# build long_df from BOTH watersheds and add `watershed` to the facets, e.g.
# facet_grid(index ~ watershed). AR will populate mainly the U-M class.
# =====================================================================
