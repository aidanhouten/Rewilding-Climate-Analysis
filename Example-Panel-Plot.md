Example-Panel-Plot
================
2026-08-06

``` r
library(terra)
library(dplyr)
library(ggplot2)
library(sf)
library(ggspatial)
library(tidyterra)
library(ggrepel)
library(patchwork)
library(stringr)
library(viridis)
library(ggmagnify)
library(raster)
library(data.table)

theme_rewild2 <- theme(
    # text elements
    plot.caption.position = "plot",
    plot.caption = element_text(hjust = 0),
    plot.title = element_text(),
    plot.subtitle = element_text(),
    axis.title = element_text(size = 11),
    axis.title.y = element_text(vjust = 3),
    axis.title.x = element_text(vjust = -0.75),
    axis.text = element_text(size = 11, colour = "black"),
    legend.text = element_text(),
    legend.title = element_text(),
    strip.text = element_text(),
    # background
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),,
    panel.background = element_rect(fill = "white"),
    #axis & border
    panel.border = element_rect(linewidth = 0.75, fill = NA),
    # axis.line = element_line(size = 0.75),
    axis.ticks = element_line(linewidth = 0.75),
    # # plot margin
    plot.margin = unit(c(5.5, 0, 0, 0), "points")
    )

# creating a theme to be applied to all plots
theme_velo2 <- theme(
    # text elements
    plot.caption.position = "plot",
    plot.caption = element_text(hjust = 0),
    plot.title = element_text(),
    plot.subtitle = element_text(),
    axis.title = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 7.5, colour = "black"),
    axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
    legend.text = element_text(),
    legend.title = element_text(),
    strip.text = element_text(),
    # background
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    #axis & border
    panel.border = element_rect(linewidth = 0.75, fill = NA),
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = 0.75),
    # # plot margin
    plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "points")
    )

# creating a theme to be applied to all plots
theme_map2 <- theme(
    # text elements
    plot.caption.position = "plot",
    plot.caption = element_text(hjust = 0),
    plot.title = element_text(),
    plot.subtitle = element_text(),
    axis.title = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_blank(),
    legend.text = element_text(),
    legend.title = element_text(),
    strip.text = element_text(),
    # background
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    # axis & border
    panel.border = element_rect(color = "black", fill = NA),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    # # plot margin
    plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "points")
    )

output_dir <- "~/UCL 2024/BIOS0034/R/Output"
min_temp_dir <- "~/UCL 2024/BIOS0034/R/Output/Temp Climate Stacks/Temperature/Min"
ca_rast <- rast("~/UCL 2024/BIOS0034/R/Data/2026 Rasters/GEE_MK_250m/MK_Stacked_250m_Central_Apennines.tif")
ca_vect <- vect("~/UCL 2024/BIOS0034/Sites/Central Apennines/CentrApennineShp.shp")
world <- read_sf("~/UCL 2024/BIOS0034/Sites/MAP OF SITES/Data/ne_10m_admin_0_countries.shp")
```

``` r
sites_df_panel <- read.csv(file.path(output_dir, "sites_df.csv"))

coords_y <- c(
  57.85, 50.98, 41.96, 51.47, 56.95, 52.51, 57.21, 40.63, 51.22,
  50.98, 53.89, 52.45, 47.51, 45.25, 46.65, 46.25, 47.16, 44.43, 54.51
  )

coords_x <- c(
  -4.73, 5.75, 13.75, 29.92, -4.54, 13.05, -4.82, -7.05, 5.62,
  -0.36, 26.43, 5.35, 21.09, 22.54, 10.16, 29.44, 26.19, 15.32, -3.34
  )

site_map_label <- c("Alladale", "Border Meuse", "Central Apennines", "Chornobyl Exclusion Zone", "Creag Meagaidh", "Döberitzer Heide", "Dundreggan", "Greater Côa Valley", "Kempen~Broek", "Knepp", "Naliboki Forest", "Oostvaardersplassen", "Pentezug", "Southern Carpathians", "Swiss National Park", "Tarutino Steppe", "Vânători-Neamț Nature Park", "Velebit Mountains", "Wild Ennerdale")

coord_points <- data.frame(coords_x, coords_y, site_map_label)

eu <- world %>% 
  filter(
    CONTINENT == "Europe"
    ) %>% 
  subset(select = geometry)

restof <- world %>% 
  filter(
    CONTINENT != "Europe"
    ) %>% 
  subset(select = geometry)
```

``` r
ca_INDVI_yearly_stack <- raster::subset(ca_rast, grep("INDVI_2", names(ca_rast), value = T))

ca_INDVI_yearly_df <- ca_INDVI_yearly_stack %>% 
  as.data.frame(xy = TRUE, na.rm = TRUE) %>% 
  rename_with(~ str_remove(., "^INDVI_")) %>% 
  pivot_longer(
    cols = -c(x, y),
    names_to = "year",
    values_to = "value"
  ) %>% 
  mutate(
    year = as.integer(year),
    pixel_id = paste(x, y, sep = "_")
    )

ca_INDVI_yearly_mean <-ca_INDVI_yearly_df %>% 
  group_by(year) %>% 
  summarise(value = mean(value)) %>% 
  ungroup() %>% 
  mutate(
    pixel_id = "mean"
  )

ca_INDVI_time_plot <-
  ggplot(ca_INDVI_yearly_df, aes(x = year, y = value, group = pixel_id)) +
  geom_line(alpha = 0.05) +
  geom_line(data = ca_INDVI_yearly_mean, aes(x = year, y = value, group = pixel_id), col = "gold", size = 1) +
  labs(x = "Year", y = "\n\nINDVI") +
  theme_rewild2
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once per session.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
ca_start_year <- sites_df_panel$start_year[3]
ca_site_name  <- sites_df_panel$site[3]
ca_years <- ca_start_year:2022

all_rasters <- lapply(ca_years, function(year) {
  
  r <- rast(paste0(
    min_temp_dir, "/", ca_site_name, "_Min_", year, ".tif"
  ))
  
  as(r, "Raster")
})

ca_temp_raster <- stack(all_rasters)

ca_temp_yearly_df <- ca_temp_raster %>% 
  as.data.frame(xy = TRUE, na.rm = TRUE)

old_ca_names <- names(ca_temp_yearly_df[,3:25])
new_ca_names <- as.character(ca_years)
ca_temp_yearly_df <- setnames(ca_temp_yearly_df, old = old_ca_names, new = new_ca_names)

ca_temp_yearly_df <- ca_temp_yearly_df %>% 
  pivot_longer(
    cols = -c(x, y),
    names_to = "year",
    values_to = "value"
  ) %>% 
  mutate(
    year = as.integer(year),
    pixel_id = paste(x, y, sep = "_")
    )

ca_temp_yearly_mean <-ca_temp_yearly_df %>% 
  group_by(year) %>% 
  summarise(value = mean(value)) %>% 
  ungroup() %>% 
  mutate(
    pixel_id = "mean"
  )

ca_temp_time_plot <- ggplot(ca_temp_yearly_df, aes(x = year, y = value, group = pixel_id)) +
  geom_line(alpha = 0.05) +
  geom_line(data = ca_temp_yearly_mean, aes(x = year, y = value, group = pixel_id), col = "gold", size = 1) +
  labs(x = "Year", y = "Mean daily\nminimum\ntemperature (°C)") +
  theme_rewild2
```

``` r
ca_INDVI_mk_rast <- ca_rast[[1]]
ca_INDVI_p_rast <- ca_rast[[2]]
ca_INDVI_sen_rast <- ca_rast[[3]]

ca_velo_rast_fpath <- sites_df_panel$raster_path[[3]]
ca_velo_rast <- rast(ca_velo_rast_fpath)

# P value classification

classify_trends <- function(slope, pval) {
  # Combine into a stack
  s <- c(slope, pval)
  
  # Apply pixel-wise classification
  app(s, fun = function(vals) {
    slope <- vals[1]
    pval <- vals[2]
    
    if (is.na(slope) | is.na(pval)) {
      return(NA)
    } else if (slope > 0 & pval < 0.05) {
      return(1)
    } else if (slope > 0 & pval >= 0.05) {
      return(2)
    } else if (slope < 0 & pval < 0.05) {
      return(3)
    } else if (slope < 0 & pval >= 0.05) {
      return(4)
    } else {
      return(NA)
    }
  })
}

trend_labels <- c(
  "Significant increase",
  "Non-significant increase",
  "Non-significant decrease",
  "Significant decrease"
)

trend_labels <- as.factor(trend_labels)
trend_labels <- factor(trend_labels, levels = trend_labels)

ca_INDVI_trend <- classify_trends(ca_INDVI_sen_rast, ca_INDVI_p_rast)
ca_INDVI_trend <- subst(ca_INDVI_trend, 1:4, trend_labels)

ca_INDVI_trend_plot <- ggplot() +
  geom_spatraster(data = ca_INDVI_trend) +
  scale_fill_viridis_d(
    na.translate = FALSE,
    name = "NDVI trend",
    option = "viridis",
    direction = -1
  ) +
  annotation_scale() +
  # annotate("text", x = 13.05, y = 42.3, label = "a)") +
  theme_velo2

# Mann-Kendall INDVI
ca_INDVI_mk_plot <- ggplot() +
  geom_spatraster(data = ca_INDVI_mk_rast) +
  # scale_fill_manual() +
  scale_fill_viridis(na.value = NA) +
  labs(fill = "Mann-Kendall\nS statistic") +
  annotation_scale() +
  # annotate("text", x = 13, y = 42.35, label = "a)") +
  theme_velo2

# Sens slope INDVI
ca_INDVI_sen_plot <- ggplot() +
  geom_spatraster(data = ca_INDVI_sen_rast) +
  # scale_fill_manual() +
  scale_fill_viridis(na.value = NA) +
  labs(fill = "Average yearly\nchange in INDVI") +
  annotation_scale() +
  # annotate("text", x = 13, y = 42.35, label = "a)") +
  theme_velo2

# P value INDVI
ca_INDVI_p_plot <- ggplot() +
  geom_spatraster(data = ca_INDVI_p_rast) +
  # scale_fill_manual() +
  scale_fill_viridis(na.value = NA) +
  labs(fill = "P value") +
  annotation_scale() +
  # annotate("text", x = 13, y = 42.35, label = "a)") +
  theme_velo2

# Temperature Velocity
ca_velo_plot <- ggplot() +
  geom_spatraster(data = ca_velo_rast) +
  scale_fill_viridis_c(na.value = NA, option = "magma") +
  labs(fill = "Minimum\ntemperature\nvelocity\n(Km/year)") +
  annotation_scale() +
  # annotate("text", x = -7.3, y = 41.17, label = "c)") +
  theme_velo2
```

``` r
map_of_sites <-
  ggplot(data = world) +
  geom_spatvector(fill = "grey95") +
  geom_point(data = coord_points, aes(x = coords_x, y = coords_y), size = 1.5, shape = 23, fill = "#F68712", color = "black") +
  geom_spatvector(data = ca_vect, fill = "grey50") +
  # geom_spatvector(data = all_sites_shp, fill = "yellow") +
  xlim(-11, 32) +
  ylim(28, 60) +
  geom_magnify(
    from = c(xmin = 13, xmax = 14.4, ymin = 41.6, ymax = 42.4),
    to = c(xmin = 0, xmax = 20, ymin = 28.5715, ymax = 40)
  ) +
  geom_label_repel(
    data = coord_points,
    aes(x = coords_x, y = coords_y, label = site_map_label),
    label.size = NA,
    alpha = 0.85,
    ylim = c(41, 60),
    seed = 1234
  ) +
  geom_label_repel(
    data = coord_points,
    aes(x = coords_x, y = coords_y, label = site_map_label),
    label.size = NA,
    fill = NA,
    ylim = c(41, 60),
    seed = 1234
  ) +
  theme_map2
```

``` r
ca_inset_multi_panel <- ca_INDVI_mk_plot + ca_INDVI_sen_plot + ca_INDVI_p_plot + ca_INDVI_trend_plot

ca_time_multi_panel <- ca_INDVI_time_plot / ca_temp_time_plot + plot_layout(axis_titles = "collect_x")

full_ca_panel_plot <- map_of_sites + ca_inset_multi_panel + ca_time_multi_panel + ca_velo_plot +
  plot_annotation(tag_levels = "a", tag_suffix = ")") & 
  theme(plot.tag = element_text(size = 14))

# full_ca_panel_plot <- (map_of_sites3 | ca_inset_multi_panel) / ((ca_INDVI_time_plot / ca_temp_time_plot + plot_layout(axis_titles = "collect_x")) | ca_velo_plot)

ggsave(filename = file.path(output_dir, "centr_ap_example_panel_plot.png"), plot = full_ca_panel_plot, device = "png", dpi = 600, width = 20, height = 15)
```
