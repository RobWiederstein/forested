library(targets)
library(tidyverse)
library(tidymodels)

# 1. Load objects from cache
tar_load(c(final_fit_results, test_data, boundary_wa_sf))

# 2. Extract predictions and add a unique ID
preds <- final_fit_results %>% 
  collect_predictions() %>%
  mutate(unique_id = row_number())

# 3. Add the SAME unique ID to test_data
# This ensures row 1 in preds matches row 1 in test_data
test_data_indexed <- test_data %>% 
  mutate(unique_id = row_number())

# 4. Join on unique_id and filter for the 126 errors
misclassified_df <- preds %>%
  # 1. Keep only the 126 mismatches
  filter(forested != .pred_class) %>%
  # 2. Join with original data
  inner_join(test_data_indexed, by = "unique_id") %>%
  # 3. Use the suffixed column names (.x or .y)
  mutate(
    error_type = case_when(
      forested.y == "Yes" & .pred_class == "No"  ~ "False Neg.",
      forested.y == "No"  & .pred_class == "Yes" ~ "False Pos."
    )
  )

# 5. Plot the Map
ggplot() +
  geom_sf(data = boundary_wa_sf, fill = "gray90", color = "gray10") +
  geom_point(data = misclassified_df, 
             aes(x = lon, y = lat, color = error_type), 
             alpha = 0.8, size = 2) +
  scale_color_manual(values = c(
    "False Neg." = "#E16A86",
    "False Pos." = "#50A315"
  )) +
  labs(color = "Error Type") +
  theme_bw() +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )
