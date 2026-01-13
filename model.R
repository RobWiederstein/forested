library(tidyverse)
library(tidymodels)
library(spatialsample)
library(sf)
library(xgboost) 
library(ranger)  
library(glmnet)  
library(vip)     

# 1. Data Setup ----
## 1.1 Convert to SF ----
forested_sf <- st_as_sf(
  forested_wa, 
  coords = c("lon", "lat"), 
  crs = 4326, 
  remove = FALSE 
)

# Split Data
set.seed(123)
split <- initial_split(forested_sf, strata = forested, prop = 0.80)
train_sf <- training(split)
test_sf  <- testing(split)

# Create Spatial Folds
# Using spatial blocks prevents the model from "cheating" by memorizing locations
set.seed(234)
spatial_folds <- spatial_block_cv(train_sf, v = 5)

# ---------------------------------------------------------
# 2. DEFINE RECIPES
# ---------------------------------------------------------

# Recipe A: For GLMNET (Linear)
# Needs: Normalization, Dummy Encoding, No Coords
rec_linear <- recipe(forested ~ ., data = train_sf) %>%
  update_role(geometry, new_role = "id") %>%
  step_rm(matches("^(lat|lon)$")) %>% 
  step_dummy(all_nominal_predictors()) %>%
  step_zv(all_predictors()) %>%
  step_normalize(all_numeric_predictors()) 

# Recipe B: For Trees (XGBoost / Random Forest)
# Needs: Dummy Encoding, No Normalization needed
rec_tree <- recipe(forested ~ ., data = train_sf) %>%
  update_role(geometry, new_role = "id") %>%
  step_rm(matches("^(lat|lon)$")) %>% 
  step_dummy(all_nominal_predictors()) %>% 
  step_zv(all_predictors())

# ---------------------------------------------------------
# 3. DEFINE MODEL SPECS
# ---------------------------------------------------------

# Model 1: Lasso Regression
spec_glmnet <- logistic_reg(penalty = tune(), mixture = 1) %>%
  set_engine("glmnet") %>%
  set_mode("classification")

# Model 2: Random Forest
spec_ranger <- rand_forest(mtry = tune(), min_n = tune(), trees = 1000) %>%
  set_engine("ranger", importance = "impurity") %>%
  set_mode("classification")

# Model 3: XGBoost
spec_xgboost <- boost_tree(
  trees = 1000, 
  tree_depth = tune(), 
  min_n = tune(), 
  loss_reduction = tune(), 
  learn_rate = tune()
) %>%
  set_engine("xgboost") %>%
  set_mode("classification")

# ---------------------------------------------------------
# 4. TOURNAMENT (Workflow Set)
# ---------------------------------------------------------

all_workflows <- workflow_set(
  preproc = list(linear = rec_linear, tree = rec_tree),
  models = list(glmnet = spec_glmnet, ranger = spec_ranger, xgboost = spec_xgboost),
  cross = TRUE 
) %>%
  # Filter to keep only logical pairings
  filter(wflow_id %in% c("linear_glmnet", "tree_ranger", "tree_xgboost"))

# Tune the models
# This step handles the heavy computation
grid_ctrl <- control_grid(save_pred = TRUE, save_workflow = TRUE)

tune_results <- all_workflows %>%
  workflow_map(
    fn = "tune_grid",
    resamples = spatial_folds,
    grid = 10,
    metrics = metric_set(roc_auc, accuracy),
    control = grid_ctrl,
    verbose = TRUE
  )

# ---------------------------------------------------------
# 5. SELECTION & VALIDATION
# ---------------------------------------------------------

# 1. Identify the Winner
cat("--- Model Leaderboard ---\n")
rank_results(tune_results, rank_metric = "roc_auc") %>%
  filter(.metric == "roc_auc") %>%
  select(wflow_id, mean, rank) %>%
  print()

best_model_id <- rank_results(tune_results, rank_metric = "roc_auc") %>%
  slice(1) %>%
  pull(wflow_id)

cat("\nThe winning model is:", best_model_id, "\n")

# 2. Extract Best Parameters
best_params <- tune_results %>%
  extract_workflow_set_result(best_model_id) %>%
  select_best(metric = "roc_auc")

# 3. Finalize the "Blueprint" (Workflow)
final_wflow <- tune_results %>%
  extract_workflow(best_model_id) %>%
  finalize_workflow(best_params)

# 4. Final Fit (Train on ALL Train Data, Evaluate on Held-out Test Data)
final_test_eval <- last_fit(final_wflow, split)

# 5. Review Metrics
collect_metrics(final_test_eval)

# ---------------------------------------------------------
# 6. VISUALIZATION
# ---------------------------------------------------------

# Plot 1: Confusion Matrix
collect_predictions(final_test_eval) %>%
  conf_mat(truth = forested, estimate = .pred_class) %>%
  autoplot(type = "heatmap") +
  labs(title = "Confusion Matrix: Washington Test Set")

# Plot 2: ROC Curve
collect_predictions(final_test_eval) %>%
  roc_curve(truth = forested, .pred_Yes) %>% 
  autoplot() +
  labs(title = "ROC Curve: Washington Test Set")

# Plot 3: Variable Importance
final_test_eval %>%
  extract_fit_parsnip() %>%
  vip(geom = "col", num_features = 12, aesthetics = list(fill = "darkred")) +
  theme_minimal() +
  labs(title = "Top Drivers of Forest Prediction (WA)")

# Plot 4: Benchmark (Boxplots of Folds)
best_configs <- rank_results(tune_results, rank_metric = "roc_auc") %>%
  filter(rank == 1) %>%
  select(wflow_id, .config)

plot_data <- collect_metrics(tune_results, summarize = FALSE) %>%
  inner_join(best_configs, by = c("wflow_id", ".config"))

ggplot(plot_data, aes(x = wflow_id, y = .estimate, fill = wflow_id)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA) + 
  geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
  facet_wrap(~.metric, scales = "free_y") +
  scale_fill_viridis_d(begin = 0.3, end = 0.8) +
  theme_minimal() +
  labs(title = "Model Stability Across Spatial Folds")

# ---------------------------------------------------------
# 7. GENERALIZATION: GEORGIA DATA (CORRECTED)
# ---------------------------------------------------------

# 1. Prepare the "Lite" Training Data
# We manually remove the columns that Georgia doesn't have.
# This prevents the recipe from ever "seeing" them.
wa_lite <- forested_sf %>%
  select(-matches("northness|eastness|aspect"))

# 2. Define a FRESH Recipe
# We cannot re-use the old recipe because it remembers the old columns.
rec_lite <- recipe(forested ~ ., data = wa_lite) %>%
  update_role(geometry, new_role = "id") %>%
  step_rm(matches("^(lat|lon)$")) %>% 
  step_dummy(all_nominal_predictors()) %>% 
  step_zv(all_predictors())

# 3. Create the Lite Workflow
# We combine our fresh recipe with the winning model (Ranger) and its best params.
wflow_lite <- workflow() %>%
  add_recipe(rec_lite) %>%
  add_model(spec_ranger) %>% 
  finalize_workflow(best_params)

# 4. Train on the Lite Data
cat("Training Lite Model (compatible with Georgia)...\n")
wa_model_lite <- fit(wflow_lite, data = wa_lite)

# 5. Predict on Georgia
# Now this will work because the model strictly requires only the columns in 'wa_lite'
ga_predictions <- augment(wa_model_lite, new_data = ga_sf)

# 6. Score Performance
ga_auc <- ga_predictions %>%
  roc_auc(truth = forested, .pred_Yes)

cat("\n--- Georgia Performance (Generalization) ---\n")
print(ga_auc)

# 7. Visualization
ga_predictions %>%
  conf_mat(truth = forested, estimate = .pred_class) %>%
  autoplot(type = "heatmap") +
  labs(
    title = "Georgia Predictions (Transfer Learning)", 
    subtitle = paste0("AUC: ", round(ga_auc$.estimate, 3))
  )