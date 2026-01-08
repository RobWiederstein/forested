# R/functions.R
library(magrittr)
library(dplyr)
library(psych)
library(sf)

# draw national map ----
library(tidyverse)
library(sf)
library(tigris)
library(rmapshaper)

# 1. Get Census State boundaries and simplify
# We use 'cb = TRUE' for cleaner cartographic lines
us_map <- states(cb = TRUE, resolution = "20m", year = 2024) %>%
  ms_simplify(keep = 0.05) %>% # Reduce vertices to 5% for performance
  shift_geometry()             # Moves AK and HI for a compact view

# 2. Convert forest_us to an 'sf' object
forest_us_sf <- forest_us %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(st_crs(us_map)) %>% # Match the map's projection
  shift_geometry()                # Shift points to match the shifted map

# 3. Create the Plot
ggplot() +
  # Draw the US backdrop
  geom_sf(data = us_map, fill = "grey98", color = "grey85", size = 0.25) +
  # Plot the forest data points colored by state
  geom_sf(data = forest_us_sf, aes(color = state), alpha = 0.4, size = 0.6) +
  # Add a specific theme for a clean atlas look
  theme_void() +
  scale_color_manual(values = c("Washington" = "#0072B2", "Georgia" = "#D55E00")) +
  labs(title = "The 'Two Islands' of the Forested Dataset",
       subtitle = "Visualizing the geographic gap between Washington and Georgia",
       color = "Source State") +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
# end ----
# draw forest coverage by state ----
library(tidyverse)
library(sf)
library(tigris)
library(rmapshaper)

# 1. Get just WA and GA boundaries and simplify
# We use 'cb = TRUE' for a clean cartographic look
wa_ga_map <- states(cb = TRUE, resolution = "20m") %>%
  filter(NAME %in% c("Washington", "Georgia")) %>%
  ms_simplify(keep = 0.1)

# 2. Convert forest_us to sf and filter for these states
forest_us_sf <- forest_us %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(st_crs(wa_ga_map))

# 3. Create the Side-by-Side Comparison Plot
library(patchwork)

# 1. Create the Washington Plot
p_wa <- ggplot() +
  geom_sf(data = filter(wa_ga_map, NAME == "Washington"), fill = "grey98") +
  geom_sf(data = filter(forest_us_sf, state == "Washington"), aes(color = forested), alpha = 0.6) +
  scale_color_manual(values = c("Yes" = "#228B22", "No" = "#D2691E")) +
  theme_minimal() +
  labs(title = "Washington") +
  theme(legend.position = "none")

# 2. Create the Georgia Plot
p_ga <- ggplot() +
  geom_sf(data = filter(wa_ga_map, NAME == "Georgia"), fill = "grey98") +
  geom_sf(data = filter(forest_us_sf, state == "Georgia"), aes(color = forested), alpha = 0.6) +
  scale_color_manual(values = c("Yes" = "#228B22", "No" = "#D2691E")) +
  theme_minimal() +
  labs(title = "Georgia")

# 3. Stitch them together
p_wa + p_ga + 
  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom')
# end ----


library(tidyverse)

# 1. Prepare data: Select numeric predictors and the outcome
forested_long <- forested %>%
  select(forested, elevation, eastness, northness, roughness, 
         dew_temp, precip_annual, starts_with("temp_"), starts_with("vapor_")) %>%
  pivot_longer(cols = -forested, names_to = "variable", values_to = "value")

# 2. Generate the density grid
ggplot(forested_long, aes(x = value, fill = forested)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ variable, scales = "free") +
  theme_minimal() +
  labs(title = "Density Distributions of Numeric Predictors",
       subtitle = "Separated by Forested Status",
       x = "Value",
       y = "Density") +
  scale_fill_manual(values = c("Yes" = "#228B22", "No" = "#D2691E"))


library(tidyverse)

# 1. Prepare and scale the data
forested_scaled <- forested %>%
  select(forested, where(is.numeric)) %>%
  pivot_longer(cols = -forested, names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  mutate(
    # Scale each variable (Z-score)
    scaled_value = (value - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE),
    # Identify extreme outliers (> 3 standard deviations)
    is_extreme = abs(scaled_value) > 3
  ) %>%
  ungroup()

# 2. Create the plot
ggplot(forested_scaled, aes(x = forested, y = scaled_value, fill = forested)) +
  # Draw boxplots without default outliers (we'll add our own colored ones)
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  # Add points for outliers, colored by whether they are > 3 SD
  geom_jitter(data = filter(forested_scaled, is_extreme == TRUE), 
              aes(color = "Extreme (>3 SD)"), width = 0.2, size = 1) +
  geom_jitter(data = filter(forested_scaled, is_extreme == FALSE & 
                              (scaled_value > (quantile(scaled_value, 0.75) + 1.5 * IQR(scaled_value)) | 
                                 scaled_value < (quantile(scaled_value, 0.25) - 1.5 * IQR(scaled_value)))),
              aes(color = "Standard Outlier"), width = 0.2, alpha = 0.3, size = 0.8) +
  facet_wrap(~variable, scales = "free_y") +
  scale_fill_manual(values = c("Yes" = "#228B22", "No" = "#D2691E")) +
  scale_color_manual(values = c("Extreme (>3 SD)" = "red", "Standard Outlier" = "black")) +
  labs(
    title = "Scaled Predictor Distributions & Extreme Outliers",
    subtitle = "Y-axis represents Standard Deviations from the Mean (Z-scores)",
    y = "Standard Deviations (σ)",
    color = "Outlier Type"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# 1. Define the recipe for the entire dataset
pca_full_recipe <- recipe(forested ~ ., data = forested) %>%
  # Keep metadata variables out of the calculation
  update_role(year, tree_no_tree, land_type, county, lon, lat, new_role = "ID") %>%
  step_naomit(all_predictors()) %>%
  step_zv(all_numeric_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_pca(all_numeric_predictors(), num_comp = 2)

# 2. Prep (train the PCA on the whole set) and Bake (apply it)
pca_estimates <- prep(pca_full_recipe)
pca_full_data <- bake(pca_estimates, new_data = NULL)

# 3. Plotting the 'Truth'
ggplot(pca_full_data, aes(x = PC1, y = PC2, color = forested)) +
  geom_point(alpha = 0.3, size = 1) + 
  scale_color_manual(values = c("Yes" = "#228B22", "No" = "#D2691E")) + # Forest green and Soil brown
  labs(
    title = "Full Dataset PCA: Class Separability",
    subtitle = "7,107 observations reduced to 2 Dimensions",
    x = paste0("PC1 (", round(pca_estimates$steps[[4]]$res$sdev[1]^2 / sum(pca_estimates$steps[[4]]$res$sdev^2)*100, 1), "%)"),
    y = paste0("PC2 (", round(pca_estimates$steps[[4]]$res$sdev[2]^2 / sum(pca_estimates$steps[[4]]$res$sdev^2)*100, 1), "%)")
  ) +
  theme_minimal()
