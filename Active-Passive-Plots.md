Active-Passive-Plots
================
2026-08-06

``` r
library(dplyr)
library(ggplot2)
library(ggrepel)

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

output_dir <- "~/UCL 2024/BIOS0034/R/Output"
```

``` r
sites_df <- read.csv(file.path(output_dir, "sites_df.csv"))

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
sites_df$actual_start <- c(1995, 2000, 1960, 1986, 1986, 2010, 2008, 1960, 2010, 2001, 1994, 1990, 1997, 2011, 1914, 2004, 2012, 2012, 2003)

sites_df <- sites_df %>% 
  mutate(actual_years = (2022 - actual_start))

# Exploratory - Years rewilding vs area rewilding
area_time_plot <-
  ggplot(sites_df, aes(x = actual_years, y = area_km2, fill = rewild_cat)) +
  geom_vline(xintercept = 22, linetype = 2, linewidth = 0.8, col = "gray50") +
  geom_label_repel(
    aes(label = site),
    size = 4.5,
    force = 3,
    show.legend = FALSE,
    seed = 1234,
    label.size = NA,
    fill = NA,
    min.segment.length = 0
    ) +
  geom_point(size = 3, shape = 23, colour = "black") +
  geom_label_repel(
    aes(label = site),
    size = 4.5,
    force = 3,
    show.legend = FALSE,
    seed = 1234,
    label.size = NA,
    fill = "white",
    alpha = 0.75) +
  scale_fill_manual(
    values = c("#F68712","#1281F6"),
    guide = guide_legend(override.aes = list(
      size = 2
    ))
  ) +
  scale_y_continuous(limits = c(8, 10000), trans = "log10", breaks = c(10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000)) +
  # scale_x_continuous(limits = c(9, 25), breaks = c(10, 15, 20, 25)) +
  labs(y = bquote("Project area (Km"^2~")"), x = "Number of years since rewilding start year", fill = "Rewilding\nstrategy") +
  theme_rewild

ggsave(filename = file.path(output_dir, "area_time_plot.png"), plot = area_time_plot, device = "png", dpi = 600, width = 10, height = 8)
```
