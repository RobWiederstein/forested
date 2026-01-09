library(targets)
library(tarchetypes)
# library(forested) # Make sure this is installed

# --- 1. Options ---
tar_option_set(
  packages = c(
    "elevatr", "forested", "ggcorrplot", "ggrepel", "gt", "magrittr",
    "tidyverse", "patchwork", "psych", "ranger", "rmapshaper",
    "sf", "terra", "tidyterra", "tigris"
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

# --- 3. The Pipeline ---
list(
  # Data Ingestion
  tar_target(forested_wa, forested::forested_wa),
  tar_target(forested_ga, forested::forested_ga),
  
  # Data Processing
  tar_target(forested_us, combine_forest(wa_data = forested_wa, ga_data = forested_ga)),
  tar_target(boundary_wa_sf, fetch_state_boundary(state = "Washington")),
  tar_target(boundary_ga_sf, fetch_state_boundary(state = "Georgia")),
  
  # Raster File Target
  tar_target(wa_elev_file, 
             create_elevation_raster(boundary_wa_sf, "data/wa_elevation.tif"), 
             format = "file"),
  
  # Maps
  tar_target(map_us, draw_us_map(forested_us)),
  tar_target(wa_ga_map, fetch_study_area(c("Washington", "Georgia"))),
  tar_target(map_wa_ga, plot_regional_comparison(forested_us, wa_ga_map)),
  
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