Rewilding-Climate-Analysis
================
2026-08-06

``` r
knitr::opts_chunk$set(echo = TRUE, eval=FALSE)
```

``` r
library(terra)
library(easyclimate)
library(forcats)
library(purrr)
library(devtools)
library(VoCC)
library(raster)
library(tidyr)
library(dplyr)
library(betareg)
library(MuMIn)
library(statmod)
library(tidyterra)
library(NoSleepR)
library(marginaleffects)
library(stringr)
library(data.table)
```

``` r
# SET UP - dataframe with site, rewilding start year and rewilding technique created

site <- c("Alladale", "Border_Meuse", "Central_Apennines", "Chernobyl_Exclusion_Zone", "Creag_Meagaidh", "Doberitzer_Heide", "Dundreggan", "Greater_Coa_Valley", "Kempen_Broek", "Knepp", "Naliboki_Forest", "Oostvaardersplassen", "Pentezug", "Southern_Carpathians", "Swiss_National_Park", "Tarutino_Steppe", "Vanatori_Neamt_Nature_Park", "Velebit_Mountains", "Wild_Ennerdale")

start_year <- c(2000, 2000, 2000, 2000, 2000, 2010, 2008, 2000, 2010, 2001, 2000, 2000, 2000, 2011, 2000, 2004, 2012, 2012, 2003)

rewild_cat <- c("Active", "Active", "Passive", "Passive", "Active", "Active",
                "Active", "Passive", "Active", "Active", "Active", "Active", "Active",
                "Active", "Passive", "Passive", "Active", "Active", "Active")

area_km2 <- c(95, 29, 5279, 4784, 39, 18, 40, 3081, 61, 11, 881, 55, 24, 3188, 170, 51, 306, 2031, 45)

sites_df <- data.frame(site, start_year, rewild_cat, area_km2)

sites_df$years_rewilding <- 2022 - sites_df$start_year
```

``` r
highres_paths <- list.files(path = "~/UCL 2024/BIOS0034/R/Data/2026 Rasters/GEE_MK_250m", full.names = TRUE, pattern = "\\.tif$")
lowres_paths <- list.files(path = "~/UCL 2024/BIOS0034/R/Data/2026 Rasters/GEE_MK_1000m", full.names = TRUE, pattern = "\\.tif$")

# below function generates a raster stack for each site and calculates the proportion of pixels showing a statistically significant positive NDVI trend
# for each of INDVI, minNDVI and max NDVI

ndvi_stats_function <- function(file_path) {
  
  stack <- rast(file_path)
  
  # Sens slope (s) used to determine direction of trend and p value for significance 
  ndvi_types <- list(
    INDVI = c(s = 1, p = 2),
    minNDVI = c(s = 4, p = 5),
    maxNDVI = c(s = 7, p = 8)
  )
  
  # Identifies the site name from the imported file path
  site <- gsub("MK_Stacked_.*?_", "", tools::file_path_sans_ext(basename(file_path)))
  
  # function to calculate the proportion of significantly greening pixels
  percentages <- lapply(ndvi_types, function(ndvi) {
    
    mk_s <- stack[[ndvi["s"]]]
    p_value <- stack[[ndvi["p"]]]
    
    total_pix <- global(!is.na(mk_s), sum, na.rm = TRUE)[,1]
    # only interested in POSITIVE significant increases
    no_sig_increase <- global((mk_s > 0) & (p_value < 0.05), sum, na.rm = TRUE)[,1]
    perc_sig_increase <- no_sig_increase/total_pix
    
    list(increase = perc_sig_increase)
    
  })
  
  return(data.frame(site = site, as.list(unlist(percentages))))
  
}

# 250-metre results
results_250m <- lapply(highres_paths, ndvi_stats_function)
table_250m <- do.call(rbind, results_250m)
colnames(table_250m)  <- paste0(colnames(table_250m), "_250m")
colnames(table_250m)[1] <- "site"

# 1000-metre results
results_1000m <- lapply(lowres_paths, ndvi_stats_function)
table_1000m <- do.call(rbind, results_1000m)
colnames(table_1000m)  <- paste0(colnames(table_1000m), "_1000m")
colnames(table_1000m)[1] <- "site"

final_table <- full_join(table_250m, table_1000m, join_by(site))

sites_df <- full_join(final_table, sites_df, join_by(site))
```

``` r
# Creating polygons from imported MK rasters to be used for climate data creation
# Using existing rasters ensures areas masked by DynamicWorld remain masked when generating climate data

sites_df$vectors <- lapply(highres_paths, function(raster) {
  
  stack <- rast(raster)
  mask_layer <- ifel(!is.na(stack[[1]]), 1, NA)
  vectors <- as.polygons(mask_layer)
  
})
```

``` r
# Using the EasyClimate package to access daily min temperatures at a 1km resolution from 1950 until 2022
# Yearly rasters of mean daily min temperature generated for each site since project start and stored on disk
# SMALLER DATASETS MAY BE ABLE TO PERFORM THIS ACROSS THE ENTIRE PROJECT AREA AS OPPOSED TO ON A SITE-BY-SITE BASIS
# TO IMPROVE PERFORMANCE - however, this project was prohibited from doing so due to EasyClimate's area limit

with_nosleep({
  
  terraOptions(memfrac = 0.6)
  
  setwd("~/UCL 2024/BIOS0034/R/Output/Temp Climate Stacks/Temperature/Min")
  
  lapply(seq_len(nrow(sites_df)), function(i) {
    
    start_year <- sites_df$start_year[[i]]
    years <- start_year:2022
    
    out_files <- lapply(seq_along(years), function(j) {
      
      year <- years[j]
      
      # get_daily_climate generates a series of rasters containing max temperature values for each day of a given period
      stack <- get_daily_climate(
        coords = sites_df$vectors[[i]],
        climatic_var = "Tmin",
        period = year,
        output = "raster"
      )
      
      # e.g. Alladale_Min_2000.tif
      filename <- paste0(sites_df$site[[i]], "_Min_", year, ".tif")
      
      # Collapses raster stack to a single-band raster of mean annual maximum temperature
      # Files saved to disk
      app(stack, mean, na.rm = TRUE, filename = filename, overwrite = TRUE)
      
      rm(stack)
      gc()
      return(filename)
      
    })
    
    gc()
    
  })
  
})
```

``` r
# VoCC package used to calculated velocity of climate change for each site
# Stored as one VoCC raster for each site

setwd("~/UCL 2024/BIOS0034/R/Output/Temp Climate Stacks/Temperature/Min")

sites_df$min_temp_velo_stack <- lapply(seq_len(nrow(sites_df)), function(i) {
  
  start_year <- sites_df$start_year[[i]]
  years <- start_year:2022
  
  all_rasters <- lapply(seq_along(years), function(j) {
    
    year <- years[j]
    
    # VoCC requires RasterStack objects from the package raster
    raster <- rast(paste0(sites_df$site[[i]], "_Min_", year, ".tif"))
    
    return(as(raster, "Raster"))
    
  })
  
  raster_stack <- stack(all_rasters)
  
  # Produces a raster of local warming/cooling rates (C/year)
  ttrend <- tempTrend(
    r = raster_stack,
    # Pixels with fewer than 25% non-NA years are masked out
    th = 0.25*nlayers(raster_stack)
  )
  
  # Spatial temperature gradient: how steeply temperature changes across the landscape (C/Km)
  sgrad <- spatGrad(
    r = raster_stack,
    projected = FALSE
  )
  
  # Adding 0.1 to prevent infinite velocities caused by spatial gradients of 0
  sgrad[[1]] <- sgrad[[1]] + 0.1
  
  climate_raster <- gVoCC(tempTrend = ttrend, spatGrad = sgrad)
  
  rm(raster_stack, ttrend, sgrad)
  
  return(climate_raster)
  
})
```

``` r
# Converting climate velocity magnitude into a single mean velocity for each site

sites_df$min_temp_mean <- sapply(seq_len(nrow(sites_df)), function(i) {
  
  temp_mag <- subset(sites_df$min_temp_velo_stack[[i]], 1)
  cellStats(temp_mag, mean, na.rm = TRUE)
  
})
```

``` r
# small conservative transformation recommended by Smithson & Verkuilen (2006) to push values away from boundaries
# REFERENCE A Better Lemon Squeezer? Maximum-Likelihood Regression With Beta-Distributed Dependent Variables

sites_df <- sites_df %>% 
  mutate(across(contains("NDVI"), ~ (. * (nrow(sites_df) - 1) + 0.5) / nrow(sites_df)))
```

``` r
sites_df_export <- sites_df %>% 
  subset(select = -c(min_temp_velo_stack))

write.csv(sites_df_export, "~/UCL 2024/BIOS0034/R/Output/sites_df.csv", row.names=FALSE)
```

=================================================================

      MODELS INVESTIGATING SIGNIFICANTLY INCREASING NDVI
                                         ----------

=================================================================

``` r
INDVI_increase_model_250 <- betareg(
  INDVI.increase_250m ~ years_rewilding * min_temp_mean + rewild_cat,
  data = sites_df,
  na.action = "na.fail"
)

INDVI_increase_dredge_250 <- dredge(INDVI_increase_model_250)
# only one model with AICc delta < 2
print(INDVI_increase_dredge_250)

# exporting dredge as CSV to use Akaike weight statistics
INDVI_dredge_250_df <- as.data.frame(INDVI_increase_dredge_250)
write.csv(INDVI_dredge_250_df, file = "INDVI_250_dredge.csv")

best_INDVI_increase_250 <- betareg(
  INDVI.increase_250m ~ years_rewilding,
  data = sites_df,
  na.action = "na.fail"
)
# number of years spent rewilding has a significant effect on recovery
# no suggestion for the inclusion of climate velocity
summary(best_INDVI_increase_250)
confint(best_INDVI_increase_250)
# each spent rewilding associated with an average 1.96 percentage point
# increase in the proportional area experiencing a significant increase in primary productivity
avg_slopes(best_INDVI_increase_250)
```

``` r
minNDVI_increase_model_250 <- betareg(
  minNDVI.increase_250m ~ years_rewilding * min_temp_mean + rewild_cat,
  data = sites_df,
  na.action = "na.fail"
)

minNDVI_increase_dredge_250 <- dredge(minNDVI_increase_model_250)
# only one model with substantial support
print(minNDVI_increase_dredge_250)

# exporting dredge as CSV to use Akaike weight statistics
minNDVI_dredge_250_df <- as.data.frame(minNDVI_increase_dredge_250)
write.csv(minNDVI_dredge_250_df, file = "minNDVI_250_dredge.csv")

best_minNDVI_increase_250 <- betareg(
  minNDVI.increase_250m ~ years_rewilding + min_temp_mean,
  data = sites_df,
  na.action = "na.fail"
)

# both time spent rewilding and climate velocity have an effect on recovery outcomes
summary(best_minNDVI_increase_250)
confint(best_minNDVI_increase_250)
# each year spent rewilding associated with an average 2.35 percentage point
# increase in the proportional area experiencing a significant increase in primary productivity
# 
avg_slopes(best_minNDVI_increase_250)
```

``` r
maxNDVI_increase_model_250 <- betareg(
  maxNDVI.increase_250m ~ years_rewilding * min_temp_mean + rewild_cat,
  data = sites_df,
  na.action = "na.fail"
)

maxNDVI_increase_dredge_250 <- dredge(maxNDVI_increase_model_250)
print(maxNDVI_increase_dredge_250)
avg_maxNDVI_increase_250 <- model.avg(maxNDVI_increase_dredge_250[1:3,])
# no significant effect of years_rewilding or climate velocity on recovery
summary(avg_maxNDVI_increase_250)

# exporting dredge as CSV to use Akaike weight statistics
maxNDVI_dredge_250_df <- as.data.frame(maxNDVI_increase_dredge_250)
write.csv(maxNDVI_dredge_250_df, file = "maxNDVI_250_dredge.csv")

best_maxNDVI_increase_250 <- betareg(
  maxNDVI.increase_250m ~ years_rewilding,
  data = sites_df,
  na.action = "na.fail"
)
# Number of years spent rewilding has a significant positive effect on recovery
summary(best_maxNDVI_increase_250)
confint(best_maxNDVI_increase_250)
# each spent rewilding associated with an average 1.24 percentage point
# increase in the proportional area experiencing a significant increase in primary productivity
avg_slopes(best_maxNDVI_increase_250)
```

=================================

      1KM RESOLUTION MODELS
      

=================================

``` r
INDVI_increase_model_1000 <- betareg(
  INDVI.increase_1000m ~ years_rewilding * min_temp_mean + rewild_cat,
  data = sites_df,
  na.action = "na.fail"
)

INDVI_increase_dredge_1000 <- dredge(INDVI_increase_model_1000)
print(INDVI_increase_dredge_1000)
avg_INDVI_increase_1000 <- model.avg(INDVI_increase_dredge_1000[1:3,])
# Model averaging approach suggests no association between variables and recovery
summary(avg_INDVI_increase_1000)

best_INDVI_increase_1000 <- betareg(
  INDVI.increase_1000m ~ years_rewilding,
  data = sites_df,
  na.action = "na.fail"
)

# Best model also suggests that years_rewilding has no effect on recovery (p = 0.08485)
summary(best_INDVI_increase_1000)
```

``` r
minNDVI_increase_model_1000 <- betareg(
  minNDVI.increase_1000m ~ years_rewilding * min_temp_mean + rewild_cat,
  data = sites_df,
  na.action = "na.fail"
)

minNDVI_increase_dredge_1000 <- dredge(minNDVI_increase_model_1000)
print(minNDVI_increase_dredge_1000)
avg_minNDVI_increase_1000 <- model.avg(minNDVI_increase_dredge_1000[1:2,])
# model averaging suggests significant effect of climate velocity on recovery
# no support for years rewilding in avg modelling
summary(avg_minNDVI_increase_1000)

best_minNDVI_increase_1000 <- betareg(
  minNDVI.increase_1000m ~ years_rewilding + min_temp_mean,
  data = sites_df,
  na.action = "na.fail"
)

# best model suggests both time spent rewilding and temp velocity have
# a significant positive effect on rewilding outcomes
summary(best_minNDVI_increase_1000)
```

``` r
maxNDVI_increase_model_1000 <- betareg(
  maxNDVI.increase_1000m ~ years_rewilding * min_temp_mean + rewild_cat,
  data = sites_df,
  na.action = "na.fail"
)

maxNDVI_increase_dredge_1000 <- dredge(maxNDVI_increase_model_1000)
print(maxNDVI_increase_dredge_1000)
avg_maxNDVI_increase_1000 <- model.avg(maxNDVI_increase_dredge_1000[1:5,])
# model averaging suggests variables have no sig effect on recovery
summary(avg_maxNDVI_increase_1000)

best_maxNDVI_increase_1000 <- betareg(
  maxNDVI.increase_1000m ~ years_rewilding,
  data = sites_df,
  na.action = "na.fail"
)
# time spent rewilding has a mild, positive effect on recovery outcomes
summary(best_maxNDVI_increase_1000)
```
