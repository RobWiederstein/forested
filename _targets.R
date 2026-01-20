library(targets)
library(tarchetypes)

# 1. Options ----
tar_option_set(
  packages = c(
    "colorspace",
    "elevatr", 
    "forested", 
    "ggcorrplot", 
    "ggrepel",
    "ggspatial",
    "gt", 
    "magrittr",
    "patchwork", 
    "psych", 
    "quarto",
    "ranger", 
    "rmapshaper",
    "sf",
    "spatialsample",
    "stringr",
    "terra", 
    "tidyterra",
    "tidymodels",
    "tidyverse",
    "tigris",
    "xgboost",   # For XGBoost
    "earth"
  ),
  format = "rds"
)

tar_source("R/functions.R")

# 3. The Pipeline ----
list(
  # constants 
  tar_target(n_folds, 10),
  # Data Ingestion
  tar_target(forested_wa, forested::forested_wa),
  tar_target(forested_ga, forested::forested_ga),
  tar_target(
    wa_sf,
    forested_wa %>% 
      sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  ),
  tar_target(
    name = eco_url,
    command = "https://dmap-prod-oms-edc.s3.us-east-1.amazonaws.com/ORD/Ecoregions/us/us_eco_l3.zip",
    format = "url" 
  ),
  tar_target(
    name = data_dir,
    command = "data/epa",
    format = "file" # Tracks the directory
  ),
  # Download data
  tar_target(
    name = eco_zip_file,
    command = get_epa_ecoregions(url = eco_url, dest_dir = data_dir),
    format = "file"
  ),
  # Data Processing
  tar_target(forested_us, combine_forest(wa_data = forested_wa, ga_data = forested_ga)),
  tar_target(boundary_wa_sf, fetch_state_boundary(state = "Washington")),
  tar_target(boundary_ga_sf, fetch_state_boundary(state = "Georgia")),
  tar_target(
    name = eco_data,
    command = process_ecoregions(
      zip_path = eco_zip_file, 
      target_states = c("Washington", "Georgia"),
      simplify_tol = 0.05
    )
  ),
  # Raster File Target
  tar_target(wa_elev_file, 
             create_elevation_raster(boundary_wa_sf, "data/wa_elevation.tif"), 
             format = "file"),
  
  # Maps
  tar_target(fig_us_map, draw_us_map()),
  tar_target(wa_ga_map, fetch_study_area(c("Washington", "Georgia"))),
  tar_target(map_wa_ga, plot_regional_comparison(forested_us, wa_ga_map)),
  tar_target(
    name = map_ecoregion_comparison,
    command = plot_ecoregion_comparison(eco_data)
  ),
  tar_target(
    name = save_plot,
    command = ggsave("figs/ecoregion_comparison.png", 
                     plot = map_ecoregion_comparison,
                     width = 12, 
                     height = 7, 
                     dpi = 300
                     ),
    format = "file"
  ),
  tar_target(plot_cv_comparison, plot_cv_strategies(forested_wa)),
  # fold mechanics
  tar_target(
    fig_fold_mechanics,
    plot_fold_mechanics(wa_sf, boundary_wa_sf)
  ),
  # fold diagram
  tar_target(fig_classic_cv, plot_classic_kfold_diagram()),
  # Analysis
  tar_target(data_summary, create_stats_summary(forested_wa)),
  tar_target(tbl_wa_summary, style_audit_table(data_summary, title = "WA Summary")),
  
  tar_target(variable_distrib, compare_univariate_distribution(forested_wa)), 
  tar_target(plt_outliers, identify_outliers(forested_wa)),
  
  tar_target(map_wa_outliers, 
             save_outlier_map_png(forested_wa, boundary_wa_sf, wa_elev_file, "figs/wa_outliers.png"),
             format = "file"),
  tar_target(plt_wa_pca, plot_wa_pca(forested_wa)),
  
  # correlogram
  tar_target(plt_correlogram, plot_correlations(forested_wa)),
  
  # vip plot
  tar_target(plt_vip, plot_rf_importance(forested_wa)),
  # umap plot
  tar_target(umap_plot, plot_umap_forested(forested_wa)),
  
  # 1. Data Splitting -------------------------------------------------
  
  # Define the split (80% Train, 20% Test)
  tar_target(splits, initial_split(wa_sf, prop = 0.80, strata = forested)),
  
  # Extract the Training Set (Used for Resampling/Modeling)
  tar_target(train_data, training(splits)),
  
  # Extract the Test Set (Locked away until the very end)
  tar_target(test_data, testing(splits)),
  # 2. Recipes ----
  ## A: Base (Includes Lat/Lon) ----
  tar_target(
    recipe_base,
    recipe(forested ~ ., data = train_data) %>%
      update_role(geometry, new_role = "id") %>%
      step_novel(all_nominal_predictors()) %>%
      step_dummy(all_nominal_predictors()) %>%
      step_zv(all_predictors()) %>%
      step_normalize(all_numeric_predictors())
  ),
  
  ## B: Non-Spatial (Bio Only) ----
  tar_target(
    recipe_non_spatial,
    recipe(forested ~ ., data = train_data) %>%
      update_role(geometry, lat, lon, new_role = "id") %>%
      step_novel(all_nominal_predictors()) %>%
      step_dummy(all_nominal_predictors()) %>%
      step_zv(all_predictors()) %>%
      step_normalize(all_numeric_predictors())
  ),
  
  ## C: Extensible (Feature Engineered) ----
  tar_target(
    recipe_extensible,
    recipe(forested ~ ., data = train_data) %>%
      update_role(geometry, lat, lon, new_role = "id") %>%
      step_rm(northness, county, year) %>%
      step_ratio(precip_annual, denom = denom_vars(temp_annual_max)) %>%
      step_mutate(
        temp_range = temp_annual_max - temp_annual_min,
        vpd_range = vapor_max - vapor_min
      ) %>%
      step_YeoJohnson(elevation) %>%
      step_novel(all_nominal_predictors()) %>%
      step_dummy(all_nominal_predictors()) %>%
      step_zv(all_predictors()) %>%
      step_normalize(all_numeric_predictors())
  ),
  tar_target(
    plot_yeo,
    plot_yeo_johnson(forested_wa)
  ),
  # 3. Engines ----
  ## Logistic Regression ----
  tar_target(
    spec_logistic,
    logistic_reg() %>% 
      set_engine("glm") %>% 
      set_mode("classification")
  ),
  
  ## MARS ----
  tar_target(
    spec_mars,
    mars(num_terms = 10, prod_degree = 2) %>% 
      set_engine("earth", nfold = 1) %>%  # nfold=1 prevents internal CV (speed)
      set_mode("classification")
  ),
  
  ## Random Forest ----
  tar_target(
    spec_rf,
    rand_forest(trees = 1000, min_n = 10) %>% 
      set_engine("ranger", 
                 importance = "impurity", # Calculate variable importance
                 num.threads = 1) %>%     # <--- Server Safety Lock
      set_mode("classification")
  ),
  
  ## XGBoost ----
  tar_target(
    spec_xgb,
    boost_tree(trees = 1000, tree_depth = 6, learn_rate = 0.01) %>% 
      set_engine("xgboost", 
                 nthread = 1) %>%         # <--- Server Safety Lock
      set_mode("classification")
  ),
  # 4. The Workflow Set ----
  # Crosses every recipe with every model (2 x 4 = 8 workflows)
  tar_target(
    model_set,
    workflow_set(
      preproc = list(base = recipe_base, 
                     non_spatial = recipe_non_spatial,
                     extensible = recipe_extensible),
      models = list(
        log = spec_logistic, 
        rf = spec_rf, 
        xgb = spec_xgb, 
        mars = spec_mars
      ),
      cross = TRUE
    )
  ),
  # 5. Resampling Strategies -----
  
  ## A. Random Folds ----
  tar_target(
    folds_random,
    vfold_cv(train_data, v = n_folds, strata = forested)
  ),
  
  ## B. Spatial Blocks ----
  tar_target(
    folds_block,
    spatial_block_cv(train_data, v = n_folds) 
  ),
  
  ## C. Spatial Clustering ----
  tar_target(
    folds_cluster,
    spatial_clustering_cv(train_data, v = n_folds) 
  ),
  # 6. Fit Models -----
  
  ## Branch 1: Random CV ----
  tar_target(
    results_random,
    workflow_map(
      model_set, 
      "fit_resamples", 
      resamples = folds_random,
      metrics = metric_set(roc_auc, accuracy, pr_auc),
      verbose = TRUE
    )
  ),
  
  ## Branch 2: Block CV ----
  tar_target(
    results_block,
    workflow_map(
      model_set, 
      "fit_resamples", 
      resamples = folds_block,
      metrics = metric_set(roc_auc, accuracy, pr_auc),
      verbose = TRUE
    )
  ),
  
  ## Branch 3: Cluster CV ----
  tar_target(
    results_cluster,
    workflow_map(
      model_set, 
      "fit_resamples", 
      resamples = folds_cluster,
      metrics = metric_set(roc_auc, accuracy, pr_auc),
      verbose = TRUE
    )
  ),
  # 7. Results ----
  tar_target(
    fig_cv_comparison,
    plot_spatial_cv_comparison(results_random, results_block, results_cluster)
  ),
  tar_target(
    fig_model_stability,
    plot_model_stability(results_random, results_block, results_cluster, best_model_id)
  ),
  # 8. Select and Tune the Best Model ----
  tar_target(
    best_model_id,
    results_cluster %>% 
      rank_results(rank_metric = "roc_auc", select_best = TRUE) %>% 
      slice(1) %>% 
      pull(wflow_id)
  ),
  tar_target(
    tbl_model_performance,
    results_cluster %>% 
      rank_results(rank_metric = "roc_auc", select_best = TRUE) %>% 
      filter(.metric == "roc_auc")
  ),
  
  # 9. Final Fit ----
  tar_target(
    final_fit_results,
    last_fit(
      extract_workflow(model_set, best_model_id),
      split = splits, # Your original 80/20 split
      metrics = metric_set(roc_auc, accuracy)
    )
  ),
  
  # 10. Test Set Performance Plot ----
  tar_target(
    fig_final_performance,
    plot_final_test_results(final_fit_results) # Use the specific plotting function
  ),
  tar_target(
    tbl_performance,
    create_performance_table(results_cluster, final_fit_results)
  ),
  # 11. Confusion Matrix ----
  tar_target(
    fig_confusion_matrix,
    plot_final_confusion_matrix(final_fit_results)
  ),
  tar_target(
    test_predictions,
    collect_predictions(final_fit_results) %>%
      dplyr::bind_cols(
        rsample::testing(splits) %>% 
          dplyr::select(lat, lon)
      )
  ),
  tar_target(
    map_wa_errors,
    save_error_map_png(
      data = test_predictions,  # <--- Use the extracted data here
      boundary_sf = boundary_wa_sf,
      raster_path = wa_elev_file,
      output_path = "figs/wa_errors.png"
    ),
    format = "file"
  ),
  # Georgia ----
  tar_target(
    model_predictors,
    c("elevation", "precip_annual", "temp_annual_mean", "roughness")
  ),
  tar_target(
    plot_aoa_ga,
    plot_georgia_aoa(
      train_data = forested_wa,
      test_data = forested_ga,
      predictors = model_predictors
    )
  ),
  # 3. Predict on Georgia using the Washington Model
  tar_target(
    ga_predictions,
    predict_external_region(
      final_fit = final_fit_results, 
      new_data = forested_ga         
    )
  ),
  
  # 4. Map the Predictions
  tar_target(
    map_ga_probs,
    plot_ga_comparison_map(
      pred_data = ga_predictions,
      boundaries = boundary_ga_sf # <--- CHECK THIS NAME
    )
  ),
  # 5. Confusion Matrix for Georgia
  tar_target(
    ga_conf_mat,
    plot_ga_confusion_matrix(ga_predictions)
  ),
  # 6. Map of Errors (False Positives + False Negatives)
  tar_target(
    map_ga_errors,
    plot_ga_error_map(
      pred_data = ga_predictions,
      boundaries = boundary_ga_sf
    )
  ),
  # Report----
  tar_quarto(report, "index.qmd", quiet = FALSE)
)
