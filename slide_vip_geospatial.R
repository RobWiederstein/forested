library(targets)
library(tidyverse)
library(tidymodels)
library(vip)
library(colorspace)

# 1. LOAD DATA FROM TARGETS
tar_load(c(results_random, results_block, results_cluster))

# 2. EXTRACTION FUNCTION
get_best_fold_vip <- function(wflow_obj, strategy_name) {
  resample_results <- wflow_obj %>% extract_workflow_set_result("spatial_rf")
  
  best_fold_id <- resample_results %>%
    collect_metrics(summarize = FALSE) %>%
    filter(.metric == "roc_auc") %>%
    slice_max(.estimate, n = 1) %>%
    slice(1) %>% 
    pull(id)
  
  best_split <- resample_results$splits[[which(resample_results$id == best_fold_id)]]
  
  wflow_obj %>%
    extract_workflow("spatial_rf") %>%
    finalize_workflow(select_best(resample_results, metric = "roc_auc")) %>%
    fit(data = analysis(best_split)) %>% 
    extract_fit_parsnip() %>%
    vi() %>%
    mutate(strategy = strategy_name)
}

# 3. CONSOLIDATE AND ORDER DATA
spatial_vip_final <- bind_rows(
  get_best_fold_vip(results_random,  "Random CV"),
  get_best_fold_vip(results_block,   "Block CV"),
  get_best_fold_vip(results_cluster, "Cluster CV")
) %>%
  filter(Variable %in% c("lat", "lon", "lat_x_lon")) %>%
  mutate(
    # Clean up names for horizontal facet labels
    Variable = case_when(
      Variable == "lon" ~ "Long",
      Variable == "lat" ~ "Lat",
      Variable == "lat_x_lon" ~ "Lat/Lon",
      TRUE ~ Variable
    ),
    # Top-to-bottom groups
    Variable = factor(Variable, levels = c("Long", "Lat", "Lat/Lon")),
    
    # Within each group: Random (Top) -> Block (Mid) -> Cluster (Bottom)
    strategy = factor(strategy, levels = c("Cluster CV", "Block CV", "Random CV"))
  )

# 4. THE FINAL "DARK 3" VERTICAL ISLAND PLOT
ggplot(spatial_vip_final, aes(x = Importance, y = strategy)) +
  # The Stem: Gray50
  geom_segment(aes(x = 0, xend = Importance, y = strategy, yend = strategy), 
               color = "gray20", size = .65) +
  
  # The Lollipop Head: Colored by strategy using Dark 3
  geom_point(aes(color = strategy), size = 5, show.legend = FALSE) +
  
  # Grouping coordinates vertically with labels on the left
  facet_grid(Variable ~ ., scales = "free_y", space = "free_y", switch = "y") +
  
  # Replaced Okabe-Ito with Dark 3
  scale_color_discrete_qualitative(palette = "Dark 3") +
  
  # Constant X-axis for absolute comparison
  scale_x_continuous(limits = c(0, 105), expand = c(0, 0)) +
  
  labs(
    title = "Spatial Variable Importance: Geographic Decay",
    subtitle = "Relative coordinate importance in the 'Champion' fold by strategy.",
    x = "Importance Score (Impurity)",
    y = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    # Rotate variable labels to be horizontal (0 degrees)
    #strip.text.y.left = element_text(angle = 0, face = "bold", size = 12, hjust = 1),
    #strip.background = element_blank(),
    #strip.placement = "outside",
    
    # Remove horizontal grid lines
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    
    # Breathing room between coordinate islands
    panel.spacing = unit(2, "lines"),
    
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(color = "black")
  )
