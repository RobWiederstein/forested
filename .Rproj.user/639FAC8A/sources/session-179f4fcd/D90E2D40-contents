# --- 0. Setup & Libraries ---
library(sf)
library(tidyverse)
library(tigris)
library(rmapshaper)
library(patchwork)
library(ggrepel)
library(stringr)

# --- Define Constants ---
DATA_DIR <- "data/epa"
ECO_URL  <- "https://dmap-prod-oms-edc.s3.us-east-1.amazonaws.com/ORD/Ecoregions/us/us_eco_l3.zip"

# --- 1. Download Function ---
download_epa_data <- function(url, dest_dir) {
  if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE)
  dest_file <- file.path(dest_dir, "us_eco_l3.zip")
  if (!file.exists(dest_file)) download.file(url, dest_file, mode = "wb")
  return(dest_file)
}

# --- 2. Processing Function (0.05 Simplification) ---
process_ecoregions_refined <- function(zip_path) {
  
  message("Processing shapefiles...")
  
  # Boundaries
  boundaries <- tigris::states(cb = TRUE, resolution = "20m", progress_bar = FALSE) %>%
    filter(NAME %in% c("Washington", "Georgia")) %>%
    rmapshaper::ms_simplify(keep = 0.05, keep_shapes = TRUE)
  
  # Load EPA
  clean_dir <- tempfile()
  dir.create(clean_dir)
  unzip(zip_path, exdir = clean_dir)
  shp_file <- list.files(clean_dir, pattern = "\\.shp$", full.names = TRUE, recursive = TRUE)[1]
  eco_raw  <- sf::read_sf(shp_file)
  
  # Fix Columns
  if (!"US_L3NAME" %in% names(eco_raw)) {
    if ("NA_L3NAME" %in% names(eco_raw)) eco_raw <- rename(eco_raw, US_L3NAME = NA_L3NAME)
    else if ("LEVEL3_NAM" %in% names(eco_raw)) eco_raw <- rename(eco_raw, US_L3NAME = LEVEL3_NAM)
  }
  
  # Clip & Filter
  suppressMessages(sf::sf_use_s2(FALSE))
  
  eco_clean <- eco_raw %>%
    st_transform(st_crs(boundaries)) %>%
    st_make_valid() %>%
    st_intersection(boundaries) %>%
    select(US_L3NAME, STATE_NAME = NAME) %>%
    rmapshaper::ms_filter_islands(min_area = 200000000) %>% 
    rmapshaper::ms_simplify(keep = 0.05, keep_shapes = TRUE)
  
  suppressMessages(sf::sf_use_s2(TRUE))
  
  return(eco_clean)
}

# --- 3. Plotting Function (Longer Stems / Less Overlap) ---
plot_comparison_clean <- function(ecoregions) {
  
  plot_state_labeled <- function(state_name, title_label) {
    
    # 1. Filter & Project
    state_e <- ecoregions %>% filter(STATE_NAME == state_name) %>% st_transform(4326)
    
    # 2. Prepare Labels
    state_labels <- state_e %>%
      group_by(US_L3NAME) %>%
      summarize(geometry = st_union(geometry)) %>% 
      mutate(
        lon = st_coordinates(st_point_on_surface(geometry))[, 1],
        lat = st_coordinates(st_point_on_surface(geometry))[, 2],
        clean_label = stringr::str_wrap(US_L3NAME, width = 15)
      ) %>%
      st_drop_geometry()
    
    ggplot() +
      # LAYER 1: Ecoregions Only
      geom_sf(data = state_e, aes(fill = US_L3NAME), 
              alpha = 1, color = "white", lwd = 0.2) +
      
      # LAYER 2: Labels (Tuned for spacing)
      geom_label_repel(
        data = state_labels,
        aes(x = lon, y = lat, label = clean_label),
        inherit.aes = FALSE, 
        size = 3.2,
        lineheight = 0.8,
        min.segment.length = 0, # Always draw the stem
        
        # --- THE FIXES ---
        box.padding = 0.7,      # More personal space around labels (was 0.5)
        force = 20,             # Stronger repulsion pushing labels apart (was 15)
        force_pull = 0.5,       # WEAKER attraction to center (was 1). Allows longer stems.
        
        max.overlaps = Inf,
        alpha = 0.95,
        segment.color = "grey30",
        segment.size = 0.3,
        seed = 42
      ) +
      
      scale_fill_viridis_d(option = "turbo") +
      theme_light() +
      labs(title = title_label, x = NULL, y = NULL) +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0, face = "bold", size = 16),
        axis.text = element_text(color = "grey60", size = 8),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.2)
      )
  }
  
  p_wa <- plot_state_labeled("Washington", "WA")
  p_ga <- plot_state_labeled("Georgia", "GA")
  
  p_wa + p_ga
}

# --- 4. Execution ---
zip_file <- download_epa_data(ECO_URL, DATA_DIR)
eco_data <- process_ecoregions_refined(zip_file)
final_plot <- plot_comparison_clean(eco_data)

print(final_plot)