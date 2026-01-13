library(tidyverse)
library(tidymodels)
library(spatialsample)
library(sf)
library(patchwork) 

# 1. SETUP DATA
# ---------------------------------------------------------
forested_sf <- st_as_sf(
  forested_wa, 
  coords = c("lon", "lat"), 
  crs = 4326, 
  remove = FALSE 
)

set.seed(1)
random_folds <- vfold_cv(forested_sf, v = 5)
spatial_folds <- spatial_block_cv(forested_sf, v = 5)

# 2. HELPER FUNCTION TO PLOT ANY SPLIT
# ---------------------------------------------------------
# This function manually pulls out the Train/Test data so we can map it
plot_split_map <- function(split, title) {
  # Extract the two sets
  train_data <- analysis(split) %>% mutate(type = "Train")
  test_data  <- assessment(split) %>% mutate(type = "Test")
  
  # Combine them
  combined <- bind_rows(train_data, test_data) %>%
    st_as_sf() # Ensure it stays an sf object
  
  # Plot
  ggplot() +
    geom_sf(data = combined, aes(color = type), size = 0.5, alpha = 0.6) +
    scale_color_manual(values = c("Train" = "gray80", "Test" = "red")) +
    ggtitle(title) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

# 3. GENERATE THE PLOTS
# ---------------------------------------------------------

# Plot Random (The Confetti)
# We take the first fold: random_folds$splits[[1]]
p_random <- plot_split_map(random_folds$splits[[1]], "Random CV (Confetti)")

# Plot Spatial (The Checkerboard)
# We take the first fold: spatial_folds$splits[[1]]
p_spatial <- plot_split_map(spatial_folds$splits[[1]], "Spatial CV (Checkerboard)")

# 4. COMPARE SIDE-BY-SIDE
# ---------------------------------------------------------
p_random + p_spatial


set.seed(42)
folds_hex <- spatial_block_cv(
  forested_sf,
  method = "hex", # <--- The switch from square to hex
  v = 5
)

autoplot(folds_hex)
library(spatialsample)

set.seed(42)
folds_env <- clustering_cv(
  forested_sf,
  # We cluster based on the environmental predictors, NOT lat/lon
  vars = c(elevation, precip_annual, temp_annual_mean, roughness),
  v = 5
)

# Visualize the clusters (These are your "Environmental Blocks")
autoplot(folds_env)