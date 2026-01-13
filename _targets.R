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
    "ranger", 
    "rmapshaper",
    "sf",
    "spatialsample",
    "stringr",
    "terra", 
    "tidyterra",
    "tidymodels",
    "tidyverse",
    "tigris"
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
    # ----------------------
  
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

# --- 3. The Pipeline ---
list(
  # Data Ingestion
  tar_target(forested_wa, forested::forested_wa),
  tar_target(forested_ga, forested::forested_ga),
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
  
  # PCA Target
  tar_target(plt_wa_pca, plot_wa_pca(forested_wa)),
  
  # correlogram
  tar_target(plt_correlogram, plot_correlations(forested_wa)),
  
  # vip plot
  tar_target(plt_vip, plot_rf_importance(forested_wa)),
  # umap plot
  tar_target(umap_plot, plot_umap_forested(forested_wa)),
  
  # Report
  tar_quarto(report, "index.qmd", quiet = FALSE)
)