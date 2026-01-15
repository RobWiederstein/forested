library(targets)
library(tidyverse)
library(tidymodels)
library(colorspace)

# 1. Load results from the targets cache
tar_load(c(results_random, results_block, results_cluster))

# 2. Function to extract per-fold metrics for the spatial_rf workflow
extract_fold_metrics <- function(wflow_obj, strategy_name) {
  wflow_obj %>%
    extract_workflow_set_result("spatial_rf") %>%
    collect_metrics(summarize = FALSE) %>% # 'summarize = FALSE' is the key here
    filter(.metric == "roc_auc") %>%
    mutate(strategy = strategy_name)
}

# 3. Combine metrics from all three strategies
all_fold_metrics <- bind_rows(
  extract_fold_metrics(results_random,  "Random CV"),
  extract_fold_metrics(results_block,   "Block CV"),
  extract_fold_metrics(results_cluster, "Cluster CV")
) %>%
  mutate(strategy = factor(strategy, levels = c("Random CV", "Block CV", "Cluster CV")))

# 4. Create the Stability Boxplot
ggplot(all_fold_metrics, aes(x = strategy, y = .estimate, color = strategy)) +
  # Use a light boxplot to show the distribution
  geom_violin(alpha = 0.4, outlier.shape = NA, width = 0.5) +
  # Add jittered points to see the individual 'folds' (regions)
  geom_jitter(width = 0.1, size = 3, alpha = 0.8) +
  
  # Professional palette
  scale_fill_discrete_qualitative(palette = "Dark 3") +
  scale_color_discrete_qualitative(palette = "Dark 3") +
  
  # Scientific notation/Formatting
  scale_y_continuous(limits = c(0.80, 1.0), breaks = seq(0.80, 1.0, 0.05)) +
  
  labs(
    title = "Model Stability across Geographic Folds",
    subtitle = "Random CV is falsely confident; Cluster CV reveals regional performance variance.",
    x = NULL,
    y = "ROC AUC (per fold)",
    caption = "Each point represents one held-out geographic cluster or random slice."
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )
