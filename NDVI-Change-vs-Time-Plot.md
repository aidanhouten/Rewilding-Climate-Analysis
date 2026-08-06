NDVI-Change-vs-Time-Plot
================
2026-08-06

``` r
library(dplyr)
library(ggplot2)
library(ggrepel)
library(betareg)
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

sites_df$site[sites_df$site == "Border_Meuse"] <- "Border Meuse"
sites_df$site[sites_df$site == "Central_Apennines"] <- "Central Apennines"
sites_df$site[sites_df$site == "Chernobyl_Exclusion_Zone"] <- "Chornobyl Exclusion Zone"
sites_df$site[sites_df$site == "Creag_Meagaidh"] <- "Creag Meagaidh"
sites_df$site[sites_df$site == "Doberitzer_Heide"] <- "Döberitzer Heide"
sites_df$site[sites_df$site == "Greater_Coa_Valley"] <- "Greater Côa Valley"
sites_df$site[sites_df$site == "Kempen_Broek"] <- "Kempen~Broek"
sites_df$site[sites_df$site == "Naliboki_Forest"] <- "Naliboki Forest"
sites_df$site[sites_df$site == "Southern_Carpathians"] <- "Southern Carpathians"
sites_df$site[sites_df$site == "Swiss_National_Park"] <- "Swiss National Park"
sites_df$site[sites_df$site == "Tarutino_Steppe"] <- "Tarutino Steppe"
sites_df$site[sites_df$site == "Vanatori_Neamt_Nature_Park"] <- "Vânători-Neamț Nature Park"
sites_df$site[sites_df$site == "Velebit_Mountains"] <- "Velebit Mountains"
sites_df$site[sites_df$site == "Wild_Ennerdale"] <- "Wild Ennerdale"
```

``` r
pred_dist_INDVI_increase <- data.frame(
  years_rewilding = seq(min(sites_df$years_rewilding), max(sites_df$years_rewilding), length.out = 100)
)

# Create an empty data frame to store bootstrap predictions
boot_preds_INDVI_increase <- data.frame(
  years_rewilding = rep(pred_dist_INDVI_increase$years_rewilding, 1000),
  pred_INDVI.increase_250m = numeric(length(pred_dist_INDVI_increase$years_rewilding) * 1000)
)

# Set up the bootstrapping loop
n_bootstrap <- 1000 # Number of bootstrap iterations
for (i in 1:n_bootstrap) {
  # Resample the original data with replacement
  boot_sample <- sites_df %>%
    slice_sample(n = nrow(sites_df), replace = TRUE)

  # Fit the new betareg model (m1a) on the new sample
  boot_model <- betareg(
    INDVI.increase_250m ~ years_rewilding,
    data = boot_sample,
    na.action = "na.fail"
  )

  # Predict with the new model
  boot_pred <- predict(boot_model, newdata = pred_dist_INDVI_increase, type = "response")

  # Store the predictions in the results data frame
  boot_preds_INDVI_increase[((i - 1) * nrow(pred_dist_INDVI_increase) + 1):(i * nrow(pred_dist_INDVI_increase)), "pred_INDVI.increase_250m"] <- boot_pred
}

# Group the predictions by the predictor variable and calculate the CIs
ci_preds_INDVI_increase <- boot_preds_INDVI_increase %>%
  group_by(years_rewilding) %>%
  summarise(
    lower_bound = quantile(pred_INDVI.increase_250m, 0.025),
    upper_bound = quantile(pred_INDVI.increase_250m, 0.975),
    INDVI.increase_250m_pred = mean(pred_INDVI.increase_250m), # Use the mean as the central prediction
    .groups = 'drop'
  ) %>%
  ungroup()

INDVI_increase_plot <- ggplot(sites_df) +
  geom_ribbon(
    data = ci_preds_INDVI_increase,
    aes(x = years_rewilding, ymin = lower_bound, ymax = upper_bound),
    fill = "#107869",
    alpha = 0.15
    ) +
  geom_line(
    data = ci_preds_INDVI_increase,
    aes(x = years_rewilding, y = INDVI.increase_250m_pred),
    color = "#107869",
    linewidth = 1
    ) +
  geom_label_repel(
    aes(x = years_rewilding, y = INDVI.increase_250m, label = site),
    size = 2.25,
    seed = 1234,
    label.size = NA,
    alpha = 0.8,
    box.padding = 0.2,
    min.segment.length = 0
    ) +
  geom_point(
    aes(x = years_rewilding, y = INDVI.increase_250m),
    size = 3,
    shape = 18,
    color = "#107869") +
  geom_label_repel(
    aes(x = years_rewilding, y = INDVI.increase_250m, label = site),
    size = 2.25,
    seed = 1234,
    label.size = NA,
    fill = NA,
    box.padding = 0.2,
    segment.color = 'transparent'
    ) +
  scale_x_continuous(limits = c(9, 23), breaks = c(10, 15, 20)) +
  scale_y_continuous(limits = c(-0.02, 1), breaks = c(0, 0.25, 0.50, 0.75, 1.00)) +
  labs(
    x = "Number of years of MODIS-observed rewilding",
    y = "Proportion of all pixels\nthat are significantly increasing"
  ) +
  annotate("text", x = 22, y = 0.99, label = expression(R^2 == 0.3409)) +
  annotate("text", x = 9.5, y = 0.99, label = "a) INDVI") +
  theme_rewild

INDVI_increase_plot
```

![](NDVI-Change-vs-Time-Plot_files/figure-gfm/Plotting%20INDVI%20Increase%20Model-1.png)<!-- -->

``` r
# Prediction grid: vary years_rewilding, hold min_temp_mean at its mean
pred_dist_minNDVI_increase <- data.frame(
  years_rewilding = seq(min(sites_df$years_rewilding), max(sites_df$years_rewilding), length.out = 100),
  min_temp_mean = mean(sites_df$min_temp_mean, na.rm = TRUE)  # <-- held constant
)

# Bootstrap
boot_preds_minNDVI_increase <- data.frame(
  years_rewilding = rep(pred_dist_minNDVI_increase$years_rewilding, 1000),
  pred_minNDVI.increase_250m = numeric(nrow(pred_dist_minNDVI_increase) * 1000)
)

n_bootstrap <- 1000
for (i in 1:n_bootstrap) {
  boot_sample <- sites_df %>%
    slice_sample(n = nrow(sites_df), replace = TRUE)
  
  boot_model <- betareg(
    minNDVI.increase_250m ~ years_rewilding + min_temp_mean,
    data = boot_sample,
    na.action = "na.fail"
  )
  
  boot_pred <- predict(boot_model, newdata = pred_dist_minNDVI_increase, type = "response")
  
  boot_preds_minNDVI_increase[((i - 1) * nrow(pred_dist_minNDVI_increase) + 1):(i * nrow(pred_dist_minNDVI_increase)), "pred_minNDVI.increase_250m"] <- boot_pred
}

# Summarise CIs
ci_preds_minNDVI_increase <- boot_preds_minNDVI_increase %>%
  group_by(years_rewilding) %>%
  summarise(
    lower_bound = quantile(pred_minNDVI.increase_250m, 0.025),
    upper_bound = quantile(pred_minNDVI.increase_250m, 0.975),
    minNDVI.increase_250m_pred = mean(pred_minNDVI.increase_250m),
    .groups = 'drop'
  ) %>%
  ungroup()

# Plot — note the y aesthetic uses the RAW data, not model-adjusted values
minNDVI_increase_plot <- ggplot(sites_df) +
  geom_ribbon(
    data = ci_preds_minNDVI_increase,
    aes(x = years_rewilding, ymin = lower_bound, ymax = upper_bound),
    fill = "#107869",
    alpha = 0.15
    ) +
  geom_line(
    data = ci_preds_minNDVI_increase,
    aes(x = years_rewilding, y = minNDVI.increase_250m_pred),
    color = "#107869",
    linewidth = 1
    ) +
  geom_label_repel(
    aes(x = years_rewilding, y = minNDVI.increase_250m, label = site),
    size = 3,
    seed = 1234,
    label.size = NA,
    alpha = 0.8,
    box.padding = 0.2,
    min.segment.length = 0
    ) +
  geom_point(
    aes(x = years_rewilding, y = minNDVI.increase_250m, col = min_temp_mean),
    size = 3,
    shape = 18
    ) +
  geom_label_repel(
    aes(x = years_rewilding, y = minNDVI.increase_250m, label = site),
    size = 3,
    seed = 1234,
    label.size = NA,
    fill = NA,
    box.padding = 0.2,
    segment.color = 'transparent'
    ) +
  scale_x_continuous(limits = c(9, 23), breaks = c(10, 15, 20)) +
  scale_y_continuous(limits = c(-0.02, 1), breaks = c(0, 0.25, 0.50, 0.75, 1.00)) +
  scale_color_viridis_c() +
  labs(
    x = "Number of years of MODIS-observed rewilding",
    y = "Proportion of all pixels\nthat are significantly increasing",
    col = "Temperature\nVelocity"
  ) +
  annotate("text", x = 22, y = 0.99, label = expression(R^2 == "0.5099")) +
  theme_rewild

minNDVI_increase_plot
```

![](NDVI-Change-vs-Time-Plot_files/figure-gfm/Plotting%20minNDVI%20Increase%20Model-1.png)<!-- -->

``` r
pred_dist_maxNDVI_increase <- data.frame(
  years_rewilding = seq(min(sites_df$years_rewilding), max(sites_df$years_rewilding), length.out = 100)
)

# Create an empty data frame to store bootstrap predictions
boot_preds_maxNDVI_increase <- data.frame(
  years_rewilding = rep(pred_dist_maxNDVI_increase$years_rewilding, 1000),
  pred_maxNDVI.increase_250m = numeric(length(pred_dist_maxNDVI_increase$years_rewilding) * 1000)
)

# Set up the bootstrapping loop
n_bootstrap <- 1000 # Number of bootstrap iterations
for (i in 1:n_bootstrap) {
  # Resample the original data with replacement
  boot_sample <- sites_df %>%
    slice_sample(n = nrow(sites_df), replace = TRUE)

  # Fit the new betareg model (m1a) on the new sample
  boot_model <- betareg(
    maxNDVI.increase_250m ~ years_rewilding,
    data = boot_sample,
    na.action = "na.fail"
  )

  # Predict with the new model
  boot_pred <- predict(boot_model, newdata = pred_dist_maxNDVI_increase, type = "response")

  # Store the predictions in the results data frame
  boot_preds_maxNDVI_increase[((i - 1) * nrow(pred_dist_maxNDVI_increase) + 1):(i * nrow(pred_dist_maxNDVI_increase)), "pred_maxNDVI.increase_250m"] <- boot_pred
}

# Group the predictions by the predictor variable and calculate the CIs
ci_preds_maxNDVI_increase <- boot_preds_maxNDVI_increase %>%
  group_by(years_rewilding) %>%
  summarise(
    lower_bound = quantile(pred_maxNDVI.increase_250m, 0.025),
    upper_bound = quantile(pred_maxNDVI.increase_250m, 0.975),
    maxNDVI.increase_250m_pred = mean(pred_maxNDVI.increase_250m), # Use the mean as the central prediction
    .groups = 'drop'
  ) %>%
  ungroup()

maxNDVI_increase_plot <- ggplot(sites_df) +
  geom_ribbon(
    data = ci_preds_maxNDVI_increase,
    aes(x = years_rewilding, ymin = lower_bound, ymax = upper_bound),
    fill = "#107869",
    alpha = 0.15
    ) +
  geom_line(
    data = ci_preds_maxNDVI_increase,
    aes(x = years_rewilding, y = maxNDVI.increase_250m_pred),
    color = "#107869",
    linewidth = 1
    ) +
  geom_label_repel(
    aes(x = years_rewilding, y = maxNDVI.increase_250m, label = site),
    size = 2.25,
    seed = 1234,
    label.size = NA,
    alpha = 0.8,
    box.padding = 0.2,
    min.segment.length = 0
    ) +
  geom_point(
    aes(x = years_rewilding, y = maxNDVI.increase_250m),
    size = 3,
    shape = 18,
    color = "#107869"
    ) +
  geom_label_repel(
    aes(x = years_rewilding, y = maxNDVI.increase_250m, label = site),
    size = 2.25,
    seed = 1234,
    label.size = NA,
    fill = NA,
    box.padding = 0.2,
    segment.color = 'transparent'
    ) +
  scale_x_continuous(limits = c(9, 23), breaks = c(10, 15, 20)) +
  scale_y_continuous(limits = c(-0.02, 1), breaks = c(0, 0.25, 0.50, 0.75, 1.00)) +
  labs(
    x = "Number of years of MODIS-observed rewilding",
    y = "Proportion of all pixels\nthat are significantly increasing"
  ) +
  annotate("text", x = 22, y = 0.99, label = expression(R^2 == 0.2358)) +
  annotate("text", x = 9.9, y = 0.99, label = "b) maxNDVI") +
  theme_rewild

maxNDVI_increase_plot
```

![](NDVI-Change-vs-Time-Plot_files/figure-gfm/Plotting%20maxNDVI%20Increase%20Model-1.png)<!-- -->

``` r
INDVI_maxNDVI_plot <- INDVI_increase_plot / maxNDVI_increase_plot + plot_layout(axes = "collect")

ggsave(filename = "INDVI_maxNDVI_plot.png", plot = INDVI_maxNDVI_plot, device = "png", dpi = 600, width = 6, height = 8)
ggsave(filename = "minNDVI_plot.png", plot = minNDVI_increase_plot, device = "png", dpi = 600, width = 8, height = 5)
```
