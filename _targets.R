library(targets)
library(tarchetypes)
# library(forested) # Make sure this is installed

# --- 1. Options ---
tar_option_set(
  packages = c(
    "colorspace",
    "elevatr", 
    "forested", 
    "ggcorrplot", 
    "ggrepel", 
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
    "earth",
    "vip"
  ),
  format = "rds"
)

# --- 2. Functions ---

combine_forest <- function(wa_data, ga_data){
  dplyr::bind_rows(
    list(WA = wa_data, GA = ga_data),
    .id = "state"
  )
}

create_stats_summary <- function(data) {
  data %>%
    sf::st_drop_geometry() %>%
    dplyr::select(where(is.numeric)) %>%
    # Keeping exclusion here as requested in original raw code
    dplyr::select(-any_of(c("year", "lon", "lat"))) %>% 
    psych::describe() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "variable") %>%
    dplyr::select(-vars, -range, -se, -mad) %>%
    dplyr::arrange(desc(abs(kurtosis)))
}

style_audit_table <- function(data, title = NULL, subtitle = NULL) {
  data %>%
    gt::gt() %>%
    gt::tab_header(title = title, subtitle = subtitle) %>%
    gt::fmt_number(columns = where(is.numeric), decimals = 2) %>%
    gt::opt_row_striping() %>%
    gt::tab_options(
      table.font.size = gt::px(14),
      column_labels.font.weight = "bold"
    )
}

draw_us_map <- function(data){
  us_map <- tigris::states(cb = TRUE, resolution = "20m", year = 2022, progress_bar = FALSE) %>%
    rmapshaper::ms_simplify(keep = 0.1) %>%
    tigris::shift_geometry()
  
  forest_us_sf <- data %>%
    sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
    sf::st_transform(sf::st_crs(us_map)) %>%
    tigris::shift_geometry()
  
  ggplot2::ggplot() +
    ggplot2::geom_sf(data = us_map, fill = "grey98", color = "grey50", size = 0.25) +
    ggplot2::geom_sf(data = forest_us_sf, aes(color = state), alpha = 0.4, size = 0.6) +
    ggplot2::theme_void() +
    ggplot2::scale_color_manual(values = c("WA" = "#0072B2", "GA" = "#D55E00")) +
    ggplot2::labs(title = "The 'Two Islands'", color = "Source State") +
    ggplot2::theme(legend.position = "bottom")
}

plot_regional_comparison <- function(data, boundaries) {
  forest_sf <- data %>%
    sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
    sf::st_transform(sf::st_crs(boundaries))
  
  plot_state <- function(state_name) {
    ggplot2::ggplot() +
      ggplot2::geom_sf(data = dplyr::filter(boundaries, NAME == state_name), fill = "grey98") +
      ggplot2::geom_sf(data = dplyr::filter(forest_sf, state == state_name), 
                       aes(color = forested), alpha = 0.6) +
      ggplot2::scale_color_manual(values = c("Yes" = "#228B22", "No" = "#D2691E")) +
      ggplot2::theme_minimal() + 
      ggplot2::labs(title = state_name) +
      ggplot2::theme(panel.grid = element_blank())
  }
  
  p_wa <- plot_state("WA") + ggplot2::theme(legend.position = "none")
  p_ga <- plot_state("GA")
  
  p_wa + p_ga + patchwork::plot_layout(guides = "collect") & ggplot2::theme(legend.position = 'bottom')
}

get_epa_ecoregions <- function(url, dest_dir = "data/epa") {
  # Ensure directory exists
  if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE)
  
  # Define destination
  dest_file <- file.path(dest_dir, "us_eco_l3.zip")
  
  # Download only if missing
  if (!file.exists(dest_file)) {
    message("Downloading EPA Ecoregion Data...")
    download.file(url, dest_file, mode = "wb")
  }
  
  # Return the file path for the next target to track
  return(dest_file)
}

process_ecoregions <- function(zip_path, target_states = c("Washington", "Georgia"), simplify_tol = 0.05) {
  require(sf)
  require(tidyverse)
  require(tigris)
  require(rmapshaper)
  
  message("Processing shapefiles for: ", paste(target_states, collapse = ", "))
  
  # A. Fetch State Boundaries (The "Cookie Cutter")
  boundaries <- tigris::states(cb = TRUE, resolution = "20m", progress_bar = FALSE) %>%
    filter(NAME %in% target_states) %>%
    # Simplify boundaries to match the ecoregion resolution
    rmapshaper::ms_simplify(keep = simplify_tol, keep_shapes = TRUE)
  
  # B. Unzip and Read EPA Data
  clean_dir <- tempfile()
  dir.create(clean_dir)
  unzip(zip_path, exdir = clean_dir)
  
  shp_file <- list.files(clean_dir, pattern = "\\.shp$", full.names = TRUE, recursive = TRUE)[1]
  eco_raw  <- sf::read_sf(shp_file)
  
  # C. Standardize Columns (Handle variations in column names)
  if (!"US_L3NAME" %in% names(eco_raw)) {
    if ("NA_L3NAME" %in% names(eco_raw)) {
      eco_raw <- rename(eco_raw, US_L3NAME = NA_L3NAME)
    } else if ("LEVEL3_NAM" %in% names(eco_raw)) {
      eco_raw <- rename(eco_raw, US_L3NAME = LEVEL3_NAM)
    }
  }
  
  # D. Clip, Clean, and Simplify
  suppressMessages(sf::sf_use_s2(FALSE)) # Turn off spherical geometry for robust clipping
  
  eco_clean <- eco_raw %>%
    st_transform(st_crs(boundaries)) %>%
    st_make_valid() %>%
    st_intersection(boundaries) %>%
    select(US_L3NAME, STATE_NAME = NAME) %>%
    rmapshaper::ms_filter_islands(min_area = 200000000) %>% 
    rmapshaper::ms_simplify(keep = simplify_tol, keep_shapes = TRUE)
  
  suppressMessages(sf::sf_use_s2(TRUE))
  
  return(eco_clean)
}

plot_ecoregion_comparison <- function(eco_data){
  require(ggplot2)
  require(patchwork)
  require(ggrepel)
  require(sf)
  
  # Helper function for individual state plots
  plot_state <- function(state_name, title_abbr) {
    
    # Filter data
    state_e <- eco_data %>% 
      filter(STATE_NAME == state_name) %>% 
      st_transform(4326)
    
    # Calculate Label Centroids
    state_labels <- state_e %>%
      group_by(US_L3NAME) %>%
      summarize(geometry = st_union(geometry)) %>% 
      mutate(
        lon = st_coordinates(st_point_on_surface(geometry))[, 1],
        lat = st_coordinates(st_point_on_surface(geometry))[, 2],
        clean_label = stringr::str_wrap(US_L3NAME, width = 15)
      ) %>%
      st_drop_geometry()
    
    # Plot
    ggplot() +
      # Geometries
      geom_sf(data = state_e, aes(fill = US_L3NAME), 
              alpha = 1, color = "white", lwd = 0.2) +
      
      # Labels (Optimized for placement)
      geom_label_repel(
        data = state_labels,
        aes(x = lon, y = lat, label = clean_label),
        inherit.aes = FALSE, 
        size = 5,
        lineheight = 0.9,
        min.segment.length = 0,
        box.padding = 0.8,
        force = 25,
        force_pull = 0.5,
        max.overlaps = Inf,
        alpha = 0.95,
        segment.color = "grey30",
        segment.size = 0.3,
        seed = 42
      ) +
      
      # Styling
      #scale_fill_viridis_d(option = "turbo") +
      scale_fill_discrete_qualitative(palette = "Dark 3") +
      theme_light() +
      labs(title = title_abbr, x = NULL, y = NULL) +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0, face = "bold", size = 16),
        axis.text = element_text(color = "grey60", size = 8),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.2)
      )
  }
  
  # Generate Subplots
  p_wa <- plot_state("Washington", "WA")
  p_ga <- plot_state("Georgia", "GA")
  
  # Combine with Patchwork
  final_plot <- p_wa + p_ga
  
  return(final_plot)
}

compare_univariate_distribution <- function(data){
  data %>%
    dplyr::select(forested, elevation, eastness, northness, roughness, 
                  dew_temp, precip_annual, starts_with("temp_"), starts_with("vapor_")) %>%
    tidyr::pivot_longer(cols = -forested, names_to = "variable", values_to = "value") %>%
    ggplot2::ggplot(aes(x = value, fill = forested)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::facet_wrap(~ variable, scales = "free") +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_manual(values = c("Yes" = "#228B22", "No" = "#D2691E"))
}

identify_outliers <- function(data){
  # 1. Prepare the data (calculate Z-scores first)
  plot_data <- data %>%
    sf::st_drop_geometry() %>%
    dplyr::select(forested, where(is.numeric)) %>%
    dplyr::select(-any_of(c("lon", "lat", "year"))) %>% 
    tidyr::pivot_longer(cols = -forested, names_to = "variable", values_to = "value") %>%
    dplyr::group_by(variable) %>%
    # Wrap scale in as.numeric to fix the matrix warning
    dplyr::mutate(scaled_value = as.numeric(scale(value))) %>% 
    dplyr::ungroup()
  
  # 2. Plot
  ggplot2::ggplot(plot_data, aes(x = forested, y = scaled_value, fill = forested)) +
    # Draw the main distribution (suppress default outliers)
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.4) +
    
    # CRITICAL FIX: Only pass the extreme values to the points layer
    ggplot2::geom_jitter(
      data = dplyr::filter(plot_data, abs(scaled_value) > 3), 
      width = 0.2, 
      alpha = 0.6,
      size = 1.5 # Make the outliers pop
    ) +
    
    ggplot2::facet_wrap(~variable, scales = "free_y") +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_manual(values = c("Yes" = "#228B22", "No" = "#D2691E")) +
    ggplot2::labs(
      title = "Standardized Distribution & Extremes",
      subtitle = "Points indicate observations > 3 Standard Deviations from mean",
      y = "Z-Score (σ)"
    )
}

fetch_study_area <- function(target_states) {
  tigris::states(cb = TRUE, resolution = "20m", progress_bar = FALSE) %>%
    dplyr::filter(NAME %in% target_states) %>%
    rmapshaper::ms_simplify(keep = 0.1, keep_shapes = TRUE)
}

fetch_state_boundary <- function(state){
  tigris::states(cb = TRUE, resolution = "20m", progress_bar = FALSE) %>%
    dplyr::filter(NAME == state) %>%
    rmapshaper::ms_simplify(keep = 0.1, keep_shapes = TRUE)
}

create_elevation_raster <- function(boundary_sf, path) {
  elev_raster <- elevatr::get_elev_raster(locations = boundary_sf, z = 7, clip = "locations")
  elev_terra <- terra::rast(elev_raster) %>% terra::mask(terra::vect(boundary_sf))
  terra::writeRaster(elev_terra, path, overwrite = TRUE)
  return(path)
}

save_outlier_map_png <- function(data, boundary_sf, raster_path, output_path) {
  elev_terra <- terra::rast(raster_path)
  wa_boundary_sf <- sf::st_transform(boundary_sf, 4326)
  
  data_sf <- data %>%
    sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
    sf::st_filter(wa_boundary_sf)
  
  numeric_cols <- data_sf %>% sf::st_drop_geometry() %>% dplyr::select(where(is.numeric)) %>% names()
  outliers_sf <- data_sf %>%
    dplyr::filter(dplyr::if_any(all_of(numeric_cols), ~ abs(scale(.)) > 3))
  
  # 4. Restore the Peaks Data (CRITICAL MISSING PIECE)
  wa_peaks_sf <- tibble::tibble(
    peak = c("Mt. Rainier", "Mt. Adams", "Mt. Baker", "Glacier Peak", "Mt. St. Helens", "Mt. Olympus"),
    lat  = c(46.8523, 46.2024, 48.7767, 48.1125, 46.1914, 47.8013),
    lon  = c(-121.7603, -121.4909, -121.8132, -121.1139, -122.1956, -123.7111)
  ) %>%
    sf::st_as_sf(coords = c("lon", "lat"), crs = 4326)
  
  plt <- ggplot2::ggplot() +
    tidyterra::geom_spatraster(data = elev_terra) + 
    tidyterra::scale_fill_hypso_c(palette = "usgs-gswa2") +
    ggplot2::geom_sf(data = wa_boundary_sf, fill = NA) +
    ggplot2::geom_sf(data = outliers_sf, aes(color = forested), size = 2) +
    ggplot2::scale_color_manual(values = c("Yes" = "#2D5A27", "No" = "#A0522D")) +
    # Restored Labels
    ggrepel::geom_label_repel(
      data = wa_peaks_sf,
      ggplot2::aes(label = peak, geometry = geometry),
      stat = "sf_coordinates",
      size = 5,
      fontface = "bold",
      box.padding = 0.5
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "", y = "") +
    ggplot2::theme(text = ggplot2::element_text(size = 14))
  
  ggplot2::ggsave(output_path, plot = plt, width = 10, height = 6)
  return(output_path)
}

plot_wa_pca <- function(data) {
  # REVISED: Keep Lat/Lon in the numeric matrix
  numeric_mat <- data %>% 
    sf::st_drop_geometry() %>%
    dplyr::select(where(is.numeric)) %>%
    dplyr::select(-any_of(c("year"))) 
  
  pca_res <- stats::prcomp(numeric_mat, scale. = TRUE)
  
  # Calculate variance
  var_exp <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
  
  # Calculate the sum for the reader
  total_var <- var_exp[1] + var_exp[2]
  
  as.data.frame(pca_res$x) %>%
    dplyr::mutate(forested = data$forested) %>%
    ggplot2::ggplot(aes(x = PC1, y = PC2, color = forested)) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::scale_color_manual(values = c("Yes" = "#2D5A27", "No" = "#A0522D")) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(
      title = "PCA: Biophysical + Geography",
      # REVISED SUBTITLE: Sums it up cleanly
      subtitle = paste0("First two components capture ", total_var, "% of total variance"),
      x = paste0("PC1 (", var_exp[1], "%)"),
      y = paste0("PC2 (", var_exp[2], "%)")
    ) +
    ggplot2::scale_x_continuous(limits = c(-10, 5))
}

plot_correlations <- function(data) {
  # 1. Data Prep (Same as before)
  cor_data <- data %>%
    sf::st_drop_geometry() %>%
    dplyr::mutate(Forested_Binary = dplyr::if_else(forested == "Yes", 1, 0)) %>%
    dplyr::select(Forested_Binary, dplyr::where(is.numeric)) %>% 
    dplyr::select(-any_of(c("year"))) 
  
  cor_df <- cor(cor_data, use = "pairwise.complete.obs") %>%
    as.data.frame() %>%
    dplyr::select(correlation = Forested_Binary) %>%
    tibble::rownames_to_column("variable") %>%
    dplyr::filter(variable != "Forested_Binary") %>% 
    dplyr::mutate(
      sign = dplyr::if_else(correlation > 0, "Positive", "Negative"),
      is_spatial = dplyr::if_else(variable %in% c("lat", "lon"), "Spatial", "Biophysical")
    )
  
  # 2. Plot with SPACE
  ggplot2::ggplot(cor_df, aes(x = reorder(variable, correlation), y = correlation)) +
    ggplot2::geom_segment(aes(xend = variable, yend = 0), color = "gray50") +
    ggplot2::geom_point(aes(color = sign, shape = is_spatial), size = 5) +
    
    ggplot2::coord_flip() +
    ggplot2::theme_minimal(base_size = 18) +
    
    ggplot2::scale_color_manual(values = c("Positive" = "#228B22", "Negative" = "#D2691E")) +
  
  ggplot2::labs(
    title = "What drives the 'Forested' classification?",
    subtitle = "Correlations with Forested (Yes=1). Note the impact of Latitude.",
    x = NULL,
    y = "Correlation Coefficient",
    color = "Direction",
    shape = "Variable Type"
  )
}

plot_rf_importance <- function(data) {
  # 1. Clean Data & Fit Model
  model_data <- data %>%
    sf::st_drop_geometry() %>%
    dplyr::select(forested, dplyr::where(is.numeric)) %>%
    dplyr::select(-any_of(c("year"))) %>%
    dplyr::mutate(forested = as.factor(forested))
  
  rf_model <- ranger::ranger(
    formula = forested ~ ., 
    data = model_data, 
    importance = "permutation",
    seed = 123
  )
  
  # 2. Extract importance & Add Spatial Flag
  imp_df <- vip::vi(rf_model) %>%
    dplyr::mutate(
      Type = dplyr::if_else(Variable %in% c("lat", "lon"), "Spatial", "Biophysical")
    )
  
  # 3. Tufte-style Lollipop with Shapes
  ggplot2::ggplot(imp_df, aes(x = Importance, y = reorder(Variable, Importance))) +
    # Stems
    ggplot2::geom_segment(aes(xend = 0, yend = Variable), color = "gray60", linewidth = 0.8) +
    
    # Dots (Mapped to Shape)
    ggplot2::geom_point(
      aes(shape = Type, color = Type), 
      size = 5
    ) +
    
    # --- VISUAL MAPPING ---
    # Shapes: Circle (16) for Biophysical, Triangle (17) for Spatial
    ggplot2::scale_shape_manual(values = c("Biophysical" = 16, "Spatial" = 17)) +
    # Colors: Forest Green for legit vars, Dark Gray for spatial leakage (optional contrast)
    ggplot2::scale_color_manual(values = c("Biophysical" = "#228B22", "Spatial" = "black")) +
  
  ggplot2::theme_minimal(base_size = 16) +
    ggplot2::labs(
      title = "Random Forest Variable Importance",
      subtitle = "Permutation Importance. Note the high rank of spatial variables.",
      y = NULL,
      x = "Importance (Loss in Accuracy if scrambled)",
      shape = "Variable Type",
      color = "Variable Type"
    ) 
}

plot_umap_forested <- function(data) {
  
  data_clean <- data %>%
    drop_na()
  
  numeric_data <- data_clean %>%
    select(where(is.numeric)) %>%
    scale()
  
  set.seed(42)
  umap_fit <- umap::umap(numeric_data)
  
  # Format for plotting
  umap_df <- umap_fit$layout %>%
    as.data.frame() %>%
    dplyr::rename(UMAP1 = V1, UMAP2 = V2) %>%
    dplyr::mutate(forested = as.factor(data_clean$forested))
  
  # Generate Plot
  ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = forested)) +
    geom_point(alpha = 0.6, size = 1.5) +
    #scale_color_viridis_d(name = "Forested Status", begin = 0.2, end = 0.8) +
    scale_color_manual(
      name = "Forested?",
      # Assumes factor order is (No, Yes). Swap colors if order differs.
      values = c("#228B22", "#D2691E")
    ) +
    theme_minimal() +
    labs(
      title = "UMAP Projection of Forested Areas",
      subtitle = "Clustering based on numeric landscape features",
      x = "UMAP Dimension 1",
      y = "UMAP Dimension 2"
    )
}

plot_cv_strategies <- function(data) {
  
  # 1. SETUP DATA & ADD ID
  # We MUST add an ID column to track rows across splits
  forested_sf <- st_as_sf(
    data, 
    coords = c("lon", "lat"), 
    crs = 4326, 
    remove = FALSE 
  ) %>%
    mutate(id = row_number()) # <--- Critical: Adds the ID we need later
  
  # 2. DEFINE CUSTOM ECOLOGICAL DISTANCE
  dist_env <- function(x) {
    x %>%
      st_drop_geometry() %>%
      select(elevation, precip_annual, temp_annual_mean, roughness) %>%
      scale() %>% 
      dist()
  }
  
  # 3. CREATE THE SPLITS
  set.seed(42)
  
  # A. Random
  random_folds <- vfold_cv(forested_sf, v = 5)
  
  # B. Spatial Block
  spatial_folds <- spatial_block_cv(forested_sf, v = 5)
  
  # C. Environmental Clustering
  cluster_folds <- spatial_clustering_cv(
    forested_sf,
    v = 5,
    cluster_function = "hclust", 
    distance_function = dist_env
  )
  
  # 4. HELPER: EXTRACT FOLD IDs
  # This assigns a Fold Number (1-5) to every single point
  get_fold_ids <- function(splits, original_data) {
    results <- original_data %>% 
      mutate(fold = NA_character_)
    
    for(i in seq_along(splits$splits)) {
      test_ids <- assessment(splits$splits[[i]]) %>% pull(id)
      results <- results %>%
        mutate(fold = ifelse(id %in% test_ids, as.character(i), fold))
    }
    return(results)
  }
  
  # 5. PREPARE PLOTTING DATA
  df_random  <- get_fold_ids(random_folds, forested_sf)
  df_spatial <- get_fold_ids(spatial_folds, forested_sf)
  df_cluster <- get_fold_ids(cluster_folds, forested_sf)
  
  # 6. COMMON PLOT FUNCTION
  plot_folds <- function(df, title) {
    ggplot(df) +
      geom_sf(aes(color = fold), size = 0.5, alpha = 0.6) +
      scale_color_brewer(palette = "Set1", name = "Fold") + # <--- The Rainbow Colors
      #scale_color_discrete_qualitative(palette = "Harmonic", name = "Fold") +
      #scale_color_discrete_qualitative(palette = "Dark 3", name = "Fold") +
      #scale_color_discrete_qualitative(palette = "Dynamic", name = "Fold") +
      ggtitle(title) +
      theme_void() + 
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        legend.position = "none" 
      )
  }
  
  # 7. GENERATE & COMBINE
  p1 <- plot_folds(df_random, "Random\n(Confetti)")
  p2 <- plot_folds(df_spatial, "Spatial Block\n(Checkerboard)")
  p3 <- plot_folds(df_cluster, "Env. Clustering\n(Ecoregions)")
  
  p1 + p2 + p3 + plot_layout(ncol = 3)
}

plot_spatial_cv_comparison <- function(res_random, res_block, res_cluster) {
  # 1. Extract & Clean Data
  all_metrics <- bind_rows(
    res_random  %>% collect_metrics() %>% mutate(strategy = "Random CV"),
    res_block   %>% collect_metrics() %>% mutate(strategy = "Block CV"),
    res_cluster %>% collect_metrics() %>% mutate(strategy = "Cluster CV")
  ) %>%
    filter(.metric == "roc_auc") %>%
    separate(wflow_id, into = c("recipe", "model"), sep = "_") %>%
    mutate(
      recipe_label = case_when(
        recipe == "base" ~ "No Coords",
        recipe == "spatial" ~ "With Coords"
      ),
      strategy = factor(strategy, levels = c("Random CV", "Block CV", "Cluster CV"))
    )
  
  # 2. Generate Plot
  ggplot(all_metrics, aes(x = model, y = mean, color = model)) +
    geom_hline(yintercept = 0.96, linetype = "dashed", color = "gray50", alpha = 0.5) +
    geom_point(size = 4) +
    geom_errorbar(aes(ymin = mean - std_err, ymax = mean + std_err), width = 0.2) +
    scale_color_discrete_qualitative(palette = "Dark 3") +
    facet_grid(recipe_label ~ strategy) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.01),
                       limits = c(.65, 1)) +
    labs(
      #title = "Model Performance: The Effect of Spatial Validation",
      #subtitle = "Comparing ROC AUC across Models, Recipes, and Validation Strategies",
      y = "ROC AUC Score",
      x = "Model Type",
      color = "Model"
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", size = 12)
    )
}

plot_spatial_vip_comparison <- function(res_random, res_block, res_cluster) {
  
  # Internal helper to extract VIP from the best fold of the spatial_rf model
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
  
  # Consolidate and clean
  vip_data <- bind_rows(
    get_best_fold_vip(res_random,  "Random CV"),
    get_best_fold_vip(res_block,   "Block CV"),
    get_best_fold_vip(res_cluster, "Cluster CV")
  ) %>%
    filter(Variable %in% c("lat", "lon", "lat_x_lon")) %>%
    mutate(
      Variable = case_when(
        Variable == "lon" ~ "Long",
        Variable == "lat" ~ "Lat",
        Variable == "lat_x_lon" ~ "Lat/Lon",
        TRUE ~ Variable
      ),
      Variable = factor(Variable, levels = c("Long", "Lat", "Lat/Lon")),
      strategy = factor(strategy, levels = c("Cluster CV", "Block CV", "Random CV"))
    )
  
  # Generate Plot
  ggplot(vip_data, aes(x = Importance, y = strategy)) +
    geom_segment(aes(x = 0, xend = Importance, y = strategy, yend = strategy), 
                 color = "gray20", linewidth = .65) +
    geom_point(aes(color = strategy), size = 5, show.legend = FALSE) +
    facet_grid(Variable ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_color_discrete_qualitative(palette = "Dark 3") +
    scale_x_continuous(limits = c(0, 105), expand = c(0, 0)) +
    labs(
      x = "Importance Score (Impurity)",
      y = NULL
    ) +
    theme_classic(base_size = 14) +
    theme(
      #panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
      #panel.spacing = unit(2, "lines"),
      #axis.text.y = element_text(color = "black"),
      #strip.text.y.left = element_text(angle = 0, face = "bold", size = 12)
    )
}
# --- 3. The Pipeline ---
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
  tar_target(map_us, draw_us_map(forested_us)),
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
  # 2. Recipes --------------------------------------------------------
  
  ## Recipe A: The "Base"
  tar_target(
    recipe_base,
    recipe(forested ~ ., data = train_data) %>% 
      update_role(geometry, lat, lon, new_role = "id") %>% 
      step_zv(all_predictors()) %>% 
      step_normalize(all_numeric_predictors()) %>% 
      step_dummy(all_nominal_predictors())
  ),
  
  # Recipe B: Spatial
  tar_target(
    recipe_spatial,
    recipe(forested ~ ., data = train_data) %>% 
      update_role(geometry, new_role = "id") %>% 
      step_zv(all_predictors()) %>% 
      step_normalize(all_numeric_predictors()) %>% 
      step_dummy(all_nominal_predictors()) %>% 
      step_interact(terms = ~ lat:lon)
  ),
  # --- 3. Engines (Models) -----------------------------------------------
  # CRITICAL: All set to single-threaded to allow parallel folding.
  
  # 1. Logistic Regression (The Baseline)
  tar_target(
    spec_logistic,
    logistic_reg() %>% 
      set_engine("glm") %>% 
      set_mode("classification")
  ),
  
  # 2. MARS (Multivariate Adaptive Regression Splines)
  # A smooth, fast alternative that handles interactions automatically.
  tar_target(
    spec_mars,
    mars(num_terms = 10, prod_degree = 2) %>% 
      set_engine("earth", nfold = 1) %>%  # nfold=1 prevents internal CV (speed)
      set_mode("classification")
  ),
  
  # 3. Random Forest (The Standard)
  # using 'ranger' which is faster than the default 'randomForest'
  tar_target(
    spec_rf,
    rand_forest(trees = 1000, min_n = 10) %>% 
      set_engine("ranger", 
                 importance = "impurity", # Calculate variable importance
                 num.threads = 1) %>%     # <--- Server Safety Lock
      set_mode("classification")
  ),
  
  # 4. XGBoost (The Powerhouse)
  tar_target(
    spec_xgb,
    boost_tree(trees = 1000, tree_depth = 6, learn_rate = 0.01) %>% 
      set_engine("xgboost", 
                 nthread = 1) %>%         # <--- Server Safety Lock
      set_mode("classification")
  ),
  # --- 4. The Workflow Set -----------------------------------------------
  # Crosses every recipe with every model (2 x 4 = 8 workflows)
  tar_target(
    model_set,
    workflow_set(
      preproc = list(base = recipe_base, spatial = recipe_spatial),
      models = list(
        log = spec_logistic, 
        rf = spec_rf, 
        xgb = spec_xgb, 
        mars = spec_mars
      ),
      cross = TRUE
    )
  ),
  # --- 5. Resampling Strategies ------------------------------------------
  
  # A. Random Folds
  tar_target(
    folds_random,
    vfold_cv(train_data, v = n_folds, strata = forested)
  ),
  
  # B. Spatial Blocks
  tar_target(
    folds_block,
    spatial_block_cv(train_data, v = n_folds) 
  ),
  
  # C. Spatial Clustering
  tar_target(
    folds_cluster,
    spatial_clustering_cv(train_data, v = n_folds) 
  ),
  # --- 6. The Execution (Fit Models) -------------------------------------
  
  # Branch 1: Random CV
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
  
  # Branch 2: Block CV
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
  
  # Branch 3: Cluster CV
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
  tar_target(
    fig_cv_comparison,
    plot_spatial_cv_comparison(results_random, results_block, results_cluster)
  ),
  tar_target(
    fig_spatial_vip,
    plot_spatial_vip_comparison(results_random, results_block, results_cluster)
  ),
  # Report
  tar_quarto(report, "index.qmd", quiet = FALSE)
)
