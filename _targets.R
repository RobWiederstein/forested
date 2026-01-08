library(targets)
library(tarchetypes)

#source("./R/functions.R")

tar_option_set(
  packages = c(
    "elevatr",
    "forested", 
    "ggrepel",
    "gt",
    "magrittr",
    "tidyverse", 
    "patchwork",
    "psych", 
    "rmapshaper",
    "sf", 
    "terra",
    "tidyterra",
    "tigris"
    ),
  # This ensures the pipeline runs in a clean environment
  format = "rds" 
)

combine_forest <- function(forest_1, forest_2){
  bind_rows(
    forest_1 %>% mutate(state = "GA"),
    forest_2 %>% mutate(state = "WA")
  )
}
create_stats_summary <- function(data) {
  data %>%
    st_drop_geometry() %>%
    select(where(is.numeric)) %>%
    select(-year, -lon, -lat) %>% 
    describe() %>%
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "variable") %>% 
    select(-vars, -range, -se, -mad) %>%
    arrange(desc(abs(kurtosis)))
}
style_audit_table <- function(data, title = NULL, subtitle = NULL) {
  data %>%
    gt() %>%
    tab_header(title = title, subtitle = subtitle) %>%
    # Standard Scientific Formatting
    fmt_number(columns = where(is.numeric), decimals = 2) %>%
    # Consistent styling
    opt_row_striping() %>%
    tab_options(
      table.font.size = px(14),
      column_labels.font.weight = "bold"
    )
}
draw_us_map <- function(data){
  us_map <- states(cb = TRUE, resolution = "20m", year = 2024) %>%
    ms_simplify(keep = 0.1) %>% 
    shift_geometry()
  
  # 2. Convert forest_us to an 'sf' object
  forest_us_sf <- data %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
    st_transform(st_crs(us_map)) %>%
    shift_geometry()
  
  # 3. Create the Plot
  ggplot() +
    # Draw the US backdrop
    geom_sf(data = us_map, fill = "grey98", color = "grey50", size = 0.25) +
    # Plot the forest data points colored by state
    geom_sf(data = forest_us_sf, aes(color = state), alpha = 0.4, size = 0.6) +
    # Add a specific theme for a clean atlas look
    theme_void() +
    scale_color_manual(values = c("WA" = "#0072B2", "GA" = "#D55E00")) +
    labs(title = "The 'Two Islands' of the Forested Dataset",
         subtitle = "Visualizing the geographic gap between Washington and Georgia",
         color = "Source State") +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
}
plot_regional_comparison <- function(data, boundaries) {
  # Prep data
  forest_sf <- data %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
    st_transform(st_crs(boundaries))
  
  # WA Plot
  p_wa <- ggplot() +
    geom_sf(data = filter(boundaries, NAME == "WA"), fill = "grey98") +
    geom_sf(data = filter(forest_sf, state == "WA"), 
            aes(color = forested), alpha = 0.6) +
    scale_color_manual(values = c("Yes" = "#228B22", "No" = "#D2691E")) +
    theme_minimal() + labs(title = "Washington") +
    theme(legend.position = "none", panel.grid = element_blank())
  
  # GA Plot
  p_ga <- ggplot() +
    geom_sf(data = filter(boundaries, NAME == "GA"), fill = "grey98") +
    geom_sf(data = filter(forest_sf, state == "GA"), 
            aes(color = forested), alpha = 0.6) +
    scale_color_manual(values = c("Yes" = "#228B22", "No" = "#D2691E")) +
    theme_minimal() + labs(title = "Georgia") +
    theme(panel.grid = element_blank())
  
  # Stitch
  p_wa + p_ga + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
}
compare_univariate_distribution <- function(data){
  # 1. Prepare data: Select numeric predictors and the outcome
  forested_long <- data %>%
    select(forested, elevation, eastness, northness, roughness, 
           dew_temp, precip_annual, starts_with("temp_"), starts_with("vapor_")) %>%
    pivot_longer(cols = -forested, names_to = "variable", values_to = "value")
  
  # 2. Generate the density grid
  ggplot(forested_long, aes(x = value, fill = forested)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~ variable, scales = "free") +
    theme_minimal(base_size = 14) +
    labs(title = "Density Distributions of Numeric Predictors",
         subtitle = "Separated by Forested Status",
         x = "Value",
         y = "Density") +
    scale_fill_manual(values = c("Yes" = "#228B22", "No" = "#D2691E"))
}
identify_outliers <- function(data){
  # 1. Prepare and scale the data
  forested_scaled <- data %>%
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
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom")
}
map_biophysical_outliers_sf <- function(data) {
  # 1. Pull Washington boundary and simplify vertices to 0.1
  # 'cb = TRUE' pulls a lower-resolution file initially to save bandwidth
  wa_boundary <- tigris::states(cb = TRUE, progress_bar = FALSE) %>%
    filter(NAME == "Washington") %>%
    rmapshaper::ms_simplify(keep = 0.1, keep_shapes = TRUE)
  
  # 2. Convert data to sf object
  data_sf <- st_as_sf(data, coords = c("lon", "lat"), crs = 4326)
  
  # 3. Identify 3-Sigma Outliers
  numeric_cols <- data %>% 
    select(where(is.numeric), -lon, -lat) %>% 
    names()
  
  outliers_sf <- data_sf %>%
    filter(if_any(all_of(numeric_cols), ~ abs(. - mean(.)) > 3 * sd(.)))
  
  # 4. Plot using geom_sf()
  ggplot() +
    geom_sf(data = wa_boundary, fill = "gray98", color = "gray80") +
    geom_sf(data = data_sf, color = "gray90", size = 0.5, alpha = 0.2) +
    geom_sf(data = outliers_sf, aes(color = forested), size = 1.2, alpha = 0.8) +
    scale_color_manual(values = c("Yes" = "#228B22", "No" = "#D2691E")) +
    theme_minimal() +
    labs(
      title = "Washington Biophysical Extremes",
      subtitle = "Simplified Geometry (10% Vertices) | ≥1 Var > 3 SD",
      color = "Forested?"
    )
}
fetch_state_boundary <- function(state){
  tigris::states(cb = TRUE, resolution = "20m", progress_bar = FALSE) %>%
    filter(NAME == state) %>%
    rmapshaper::ms_simplify(keep = 0.1, keep_shapes = TRUE)
}
create_elevation_raster <- function(boundary_sf, path = "wa_elevation.tif") {
  # 1. Fetch
  elev_raster <- elevatr::get_elev_raster(locations = boundary_sf, z = 7, clip = "locations")
  
  # 2. Process & Save to Disk
  elev_terra <- terra::rast(elev_raster) %>% 
    terra::mask(terra::vect(boundary_sf))
  
  terra::writeRaster(elev_terra, path, overwrite = TRUE)
  
  return(path) # Return the file path for targets to track
}
save_outlier_map_png <- function(data, boundary_sf, raster_path, output_path) {
  # 1. Re-load the raster from the file target
  elev_terra <- terra::rast(raster_path)
  
  # 2. Re-project and Filter (Internal to ensure fresh pointers)
  wa_boundary_sf <- sf::st_transform(boundary_sf, 4326)
  data_sf <- data %>%
    sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
    sf::st_filter(wa_boundary_sf)
  
  numeric_cols <- data_sf %>% 
    dplyr::as_tibble() %>%
    dplyr::select(dplyr::where(is.numeric), -geometry) %>% 
    names()
  
  outliers_sf <- data_sf %>%
    dplyr::filter(dplyr::if_any(dplyr::all_of(numeric_cols), ~ abs(. - mean(.)) > 3 * sd(.)))
  
  # 3. Add 5 highest peaks
  wa_peaks_sf <- tibble::tibble(
    peak = c("Mt. Rainier", "Mt. Adams", "Mt. Baker", "Glacier Peak", "Mt. St. Helens", "Mt. Olympus"),
    lat  = c(46.8523, 46.2024, 48.7767, 48.1125, 46.1914, 47.8013),
    lon  = c(-121.7603, -121.4909, -121.8132, -121.1139, -122.1956, -123.7111)
  ) %>%
    sf::st_as_sf(coords = c("lon", "lat"), crs = 4326)
  
  # 4. Create the Plot
  plt <- ggplot2::ggplot() +
    tidyterra::geom_spatraster(data = elev_terra) + 
    tidyterra::scale_fill_hypso_c(palette = "usgs-gswa2", name = "Elevation (m)") +
    ggplot2::geom_sf(data = wa_boundary_sf, fill = NA, color = "black", linewidth = 0.4) +
    ggplot2::geom_sf(data = outliers_sf, ggplot2::aes(color = forested), size = 1.8) +
    ggplot2::scale_color_manual(values = c("Yes" = "#2D5A27", "No" = "#A0522D")) +
    # New Peaks Layer
    ggrepel::geom_label_repel(
      data = wa_peaks_sf,
      ggplot2::aes(label = peak, geometry = geometry),
      stat = "sf_coordinates",
      size = 5,
      fontface = "bold",
      box.padding = 0.5
    ) +
    ggplot2::theme_minimal() +
    labs(x = "", y = "") +
    ggplot2::theme(text = ggplot2::element_text(size = 14)
    )
  
  # 5. Physical Render (Baking the pixels)
  ggplot2::ggsave(output_path, plot = plt, width = 14, height = 7, dpi = 300)
  
  return(output_path)
}
plot_wa_pca <- function(data, x_limits = c(-10, 5)) {
  
  # 1. Focus on biophysical variables only
  pca_input <- data
  
  # 2. Separate metadata from math
  labels_vec <- data$forested 
  numeric_mat <- pca_input %>% dplyr::select(dplyr::where(is.numeric))
  
  # 3. Compute PCA and variance stats
  pca_res <- stats::prcomp(numeric_mat, scale. = TRUE)
  var_exp <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
  
  # 4. Reconstruct for plotting
  pca_df <- as.data.frame(pca_res$x) %>%
    dplyr::mutate(forested = labels_vec)
  
  # 5. Return the ggplot object
  ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, color = as.factor(forested))) +
    ggplot2::geom_point(alpha = 0.5, size = 0.75) +
    ggplot2::scale_color_manual(
      values = c("No" = "#A0522D", "Yes" = "#2D5A27"),
      name = "Forested?"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::scale_x_continuous(limits = x_limits) +
    ggplot2::labs(
      title = "Washington: Coordinates Included",
      subtitle = paste0("Explained Variance: PC1 (", var_exp[1], "%), PC2 (", var_exp[2], "%)"),
      x = paste0("PC1 (", var_exp[1], "%)"),
      y = paste0("PC2 (", var_exp[2], "%)")
    )
}
list(
  # objects
  tar_target(boundary_wa_sf, fetch_state_boundary(state = "Washington")),
  tar_target(
    wa_elev_file, 
    create_elevation_raster(boundary_wa_sf, "data/wa_elevation.tif"), 
    format = "file" 
  ),
  tar_target(forest_us, combine_forest(forested_ga, forested_wa)),
  tar_target(data_summary, create_stats_summary(forest_us)),
  tar_target(tbl_forest_us_summary, style_audit_table(data_summary, title = "Summary Stats")),
  tar_target(map_us, draw_us_map(forest_us)),
  tar_target(
    wa_ga_map,
    states(cb = TRUE, resolution = "20m") %>%
      filter(NAME %in% c("Washington", "Georgia")) %>%
      ms_simplify(keep = 0.1)
  ),
  tar_target(map_wa_ga, plot_regional_comparison(forest_us, wa_ga_map)),
  tar_target(variable_distrib, compare_univariate_distribution(forested::forested)),
  tar_target(plt_outliers, identify_outliers(forested)),
  tar_target(map_wa_outliers, 
             save_outlier_map_png(
               forested::forested, 
               boundary_wa_sf, 
               wa_elev_file, 
               "figs/wa_outliers.png"),
             format = "file"),
  targets::tar_target(plt_wa_pca, plot_wa_pca(forested::forested_wa)),
  tar_quarto(report, "index.qmd", quiet = FALSE)
)
