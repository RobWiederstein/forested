library(waywiser)
library(dplyr)
library(sf)
library(ggplot2)
library(tibble)

# 1. Define predictors
predictors <- c("elevation", "precip_annual", "temp_annual_mean", "roughness")

# 2. Prepare Data (Select only predictor columns for the calculation)
train_env <- forested_wa %>% select(all_of(predictors))
test_env <- forested_ga %>% select(all_of(predictors))

# 3. Create Dummy Importance Table
# This tells waywiser: "Treat all 4 variables as equally important"
equal_importance <- tibble::tibble(term = predictors, estimate = 1)

# 4. Fit the AOA Model on Washington Data
aoa_model <- ww_area_of_applicability(
  x = train_env,       
  importance = equal_importance
)

# 5. Predict on Georgia Data
# Returns a tibble with 'di' (Dissimilarity Index) and 'aoa' (Logical)
aoa_results <- predict(aoa_model, test_env)

# 6. Bind Results and Make Spatial (The Fix)
# We bind the new columns to the original data, THEN convert to sf using 'lon'/'lat'
ga_aoa_sf <- forested_ga %>%
  dplyr::bind_cols(aoa_results) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# 7. Visualization
ggplot(ga_aoa_sf) +
  geom_sf(aes(color = di), size = 0.5) +
  scale_color_viridis_c(
    option = "magma", 
    direction = -1, 
    name = "Dissimilarity\n(Risk Score)"
  ) +
  ggtitle("Transferability Risk: Washington -> Georgia") +
  theme_void()