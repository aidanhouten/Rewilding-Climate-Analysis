Ridge-Plot
================
2026-08-06

``` r
library(dplyr)
library(terra)
library(ggplot2)
library(ggrepel)
library(ggridges)
library(forcats)
library(viridis)
library(patchwork)

theme_rewild <- theme(
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
    plot.margin = margin(12,12,12,12, "pt")
    )
```

``` r
sites_df <- read.csv("~/UCL 2024/BIOS0034/R/Output/sites_df.csv")

highres_paths <- list.files(path = "~/UCL 2024/BIOS0034/R/Data/2026 Rasters/GEE_MK_250m", full.names = TRUE, pattern = "\\.tif$")

# Function
sens_slope_raster_function <- function(file_path) {
  
  stack <- rast(file_path)
  
  # Identifies the site name from the imported file path
  site <- gsub("MK_Stacked_.*?_", "", tools::file_path_sans_ext(basename(file_path)))
  
  sen_raster_idx <- list(
    INDVI = 3,
    minNDVI = 6,
    maxNDVI = 9
  )
  
  # extracting values from each band as a named vector
  vals <- lapply(sen_raster_idx, function(idx) {
    na.omit(as.numeric(values(stack[[idx]])))
  })
  
  # ensuring all bands have the same length
  n <- min(lengths(vals))
  
  data.frame(
    site = site,
    INDVI = vals$INDVI[1:n],
    minNDVI = vals$minNDVI[1:n],
    maxNDVI = vals$maxNDVI[1:n]
  )
  
}

sens_slope_df <- bind_rows(lapply(highres_paths, sens_slope_raster_function)) %>% 
  left_join(sites_df %>% select(site, min_temp_mean), by = "site")

sens_slope_df$site[sens_slope_df$site == "Border_Meuse"] <- "Border Meuse"
sens_slope_df$site[sens_slope_df$site == "Central_Apennines"] <- "Central Apennines"
sens_slope_df$site[sens_slope_df$site == "Chernobyl_Exclusion_Zone"] <- "Chornobyl Exclusion Zone"
sens_slope_df$site[sens_slope_df$site == "Creag_Meagaidh"] <- "Creag Meagaidh"
sens_slope_df$site[sens_slope_df$site == "Doberitzer_Heide"] <- "Döberitzer Heide"
sens_slope_df$site[sens_slope_df$site == "Greater_Coa_Valley"] <- "Greater Côa Valley"
sens_slope_df$site[sens_slope_df$site == "Kempen_Broek"] <- "Kempen~Broek"
sens_slope_df$site[sens_slope_df$site == "Naliboki_Forest"] <- "Naliboki Forest"
sens_slope_df$site[sens_slope_df$site == "Southern_Carpathians"] <- "Southern Carpathians"
sens_slope_df$site[sens_slope_df$site == "Swiss_National_Park"] <- "Swiss National Park"
sens_slope_df$site[sens_slope_df$site == "Tarutino_Steppe"] <- "Tarutino Steppe"
sens_slope_df$site[sens_slope_df$site == "Vanatori_Neamt_Nature_Park"] <- "Vânători-Neamț Nature Park"
sens_slope_df$site[sens_slope_df$site == "Velebit_Mountains"] <- "Velebit Mountains"
sens_slope_df$site[sens_slope_df$site == "Wild_Ennerdale"] <- "Wild Ennerdale"

sens_slope_df <- sens_slope_df %>% 
  mutate(
    site = fct_reorder(site, desc(min_temp_mean))
  )

sens_slope_sd <- sens_slope_df %>% 
  summarise(
    twosd_above_mean_INDVI = mean(INDVI) + sd(INDVI)*2,
    twosd_below_mean_INDVI = mean(INDVI) - sd(INDVI)*2,
    twosd_above_mean_minNDVI = mean(minNDVI) + sd(minNDVI)*2,
    twosd_below_mean_minNDVI = mean(minNDVI) - sd(minNDVI)*2,
    twosd_above_mean_maxNDVI = mean(maxNDVI) + sd(maxNDVI)*2,
    twosd_below_mean_maxNDVI = mean(maxNDVI) - sd(maxNDVI)*2,
  )
```

``` r
INDVI_ridge_plot <- ggplot(sens_slope_df, aes(x = INDVI, y = site, fill = min_temp_mean)) +
  xlim(sens_slope_sd$twosd_below_mean_INDVI, sens_slope_sd$twosd_above_mean_INDVI) +
  stat_density_ridges(rel_min_height = 0.01, bandwidth = 0.00075, alpha = 0.9) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x = NULL, y = "Rewilding site") +
  scale_fill_viridis(name = "Mean\nminimum\ntemperature\nvelocity") +
  # annotate("text", x = 9.5, y = 0.99, label = "a) INDVI") +
  theme_rewild

min_ridge_plot <- ggplot(sens_slope_df, aes(x = minNDVI, y = site, fill = min_temp_mean)) +
  xlim(sens_slope_sd$twosd_below_mean_minNDVI, sens_slope_sd$twosd_above_mean_minNDVI) +
  stat_density_ridges(rel_min_height = 0.01, bandwidth = 0.00075, alpha = 0.9) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x = "Average change in NDVI", y = "Rewilding site") +
  scale_fill_viridis(name = "Mean\nminimum\ntemperature\nvelocity", guide = "none") +
  # annotate("text", x = -0.006, y = 0.99, label = "b) minNDVI") +
  theme_rewild

max_ridge_plot <- ggplot(sens_slope_df, aes(x = maxNDVI, y = site, fill = min_temp_mean)) +
  xlim(sens_slope_sd$twosd_below_mean_maxNDVI, sens_slope_sd$twosd_above_mean_maxNDVI) +
  stat_density_ridges(rel_min_height = 0.01, bandwidth = 0.00075, alpha = 0.9) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x = "Average change in NDVI", y = "Rewilding site") +
  scale_fill_viridis(name = "Mean\nminimum\ntemperature\nvelocity", guide = "none") +
  # annotate("text", x = 9.9, y = 0.99, label = "c) maxNDVI") +
  theme_rewild

min_max_ridge_plot <- (min_ridge_plot + max_ridge_plot) + plot_layout(axes = "collect_y", axis_titles = "collect")

whole_ridge_plot <- INDVI_ridge_plot / min_max_ridge_plot + plot_layout(axis_titles = "collect_x")

ggsave(filename = "ridge_plot.png", plot = whole_ridge_plot, device = "png", dpi = 600, width = 10, height = 8)
```
