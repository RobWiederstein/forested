# --- 0. Setup & Libraries ---
library(sf)
library(leaflet)
library(tidyverse)
library(rmapshaper)
library(tigris)

# --- Define Constants ---
DATA_DIR <- "data/epa"
ECO_URL  <- "https://dmap-prod-oms-edc.s3.us-east-1.amazonaws.com/ORD/Ecoregions/us/us_eco_l3.zip"

# --- 1. Download & Processing Function ---
get_wa_ecoregions_fixed <- function(url, dest_dir) {
  
  if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE)
  dest_file <- file.path(dest_dir, "us_eco_l3.zip")
  if (!file.exists(dest_file)) download.file(url, dest_file, mode = "wb")
  
  # Get WA Border (The Cookie Cutter)
  wa_border <- tigris::states(cb = TRUE, resolution = "20m", progress_bar = FALSE) %>%
    filter(NAME == "Washington") %>%
    st_transform(4326)
  
  # Read EPA Data
  clean_dir <- tempfile()
  unzip(dest_file, exdir = clean_dir)
  shp_file <- list.files(clean_dir, pattern = "\\.shp$", full.names = TRUE, recursive = TRUE)[1]
  eco_raw <- sf::read_sf(shp_file)
  
  # Clean Columns
  if (!"US_L3NAME" %in% names(eco_raw)) {
    if ("NA_L3NAME" %in% names(eco_raw)) eco_raw <- rename(eco_raw, US_L3NAME = NA_L3NAME)
    else if ("LEVEL3_NAM" %in% names(eco_raw)) eco_raw <- rename(eco_raw, US_L3NAME = LEVEL3_NAM)
  }
  
  # --- STEP 1: Clip to Washington ---
  suppressMessages(sf::sf_use_s2(FALSE))
  
  eco_clipped <- eco_raw %>%
    st_transform(st_crs(wa_border)) %>%
    st_make_valid() %>%
    st_intersection(wa_border) %>% 
    select(US_L3NAME) %>%
    mutate(US_L3NAME = str_trim(US_L3NAME)) %>%
    rmapshaper::ms_filter_islands(min_area = 200000000) %>% 
    rmapshaper::ms_simplify(keep = 0.05, keep_shapes = TRUE)
  
  # --- STEP 2: The Geographic Surgery (Split Region 77) ---
  # We identify "North Cascades" polygons that are actually on the Olympic Peninsula (West of Seattle)
  # and rename them to "Olympic High Peaks".
  
  eco_final <- eco_clipped %>%
    mutate(
      centroid_lon = st_coordinates(st_centroid(geometry))[,1],
      US_L3NAME = case_when(
        US_L3NAME == "North Cascades" & centroid_lon < -122.5 ~ "Olympic High Peaks",
        TRUE ~ US_L3NAME
      )
    ) %>%
    select(-centroid_lon) # Cleanup
  
  suppressMessages(sf::sf_use_s2(TRUE))
  
  return(eco_final)
}

# --- 2. The Visual Dictionary (Updated for Split) ---
image_lookup <- tribble(
  ~US_L3NAME, ~img_url, ~desc,
  
  "Olympic High Peaks",
  "https://commons.wikimedia.org/wiki/Special:FilePath/Mount_Olympus_Washington.jpg?width=500",
  "<b>The Olympic Core.</b><br>The rugged, glaciated interior of the Olympic Peninsula, centered on Mt. Olympus.",
  
  "North Cascades", 
  "https://commons.wikimedia.org/wiki/Special:FilePath/Mount_Shuksan_reflection_in_Picture_Lake_during_autumn.jpg?width=500", 
  "<b>The American Alps.</b><br>Jagged, uplifted granite peaks distinct from the volcanoes to the south.",
  
  "Coast Range", 
  "https://commons.wikimedia.org/wiki/Special:FilePath/Forks_WA_Hoh_National_Forest_Trail.JPG?width=500", 
  "<b>The Rain Forest.</b><br>Dense coniferous forests, high rainfall, and steep slopes defined by the Pacific Ocean.",
  
  "Puget Lowland", 
  "https://commons.wikimedia.org/wiki/Special:FilePath/Mazama_Pocket_gopher_habitat.jpg?width=500", 
  "<b>Glacial Troughs & Prairies.</b><br>Home to rare glacial outwash prairies and the majority of the state's population.",
  
  "Cascades", 
  "https://commons.wikimedia.org/wiki/Special:FilePath/Mount_Rainier_from_west.jpg?width=500", 
  "<b>Volcanic Peaks.</b><br>Dominated by active and dormant volcanoes like Mt. Rainier and St. Helens.",
  
  "Columbia Plateau", 
  "https://commons.wikimedia.org/wiki/Special:FilePath/Alkali_Lake_Washington.jpg?width=500", 
  "<b>Scablands & Sage.</b><br>Arid shrub-steppe and dramatic basalt coulees carved by ancient floods.",
  
  "Blue Mountains", 
  "https://commons.wikimedia.org/wiki/Special:FilePath/John_Day_Fossil_Beds.jpg?width=500", 
  "<b>Island Ranges.</b><br>Complex mountains and volcanic basins, home to the colorful John Day Fossil Beds.",
  
  "Northern Rockies", 
  "https://commons.wikimedia.org/wiki/Special:FilePath/Salmo-Priest_Wilderness_Area.jpg?width=500", 
  "<b>Inland Highlands.</b><br>High, steep mountains with a continental climate, distinct from the coastal ranges.",
  
  "Eastern Cascades Slopes and Foothills",
  "https://commons.wikimedia.org/wiki/Special:FilePath/Sagebrush_with_shattered_trunk.jpg?width=500",
  "<b>The Rain Shadow.</b><br>A stark transition zone where mountain forests give way to dry sagebrush slopes.",
  
  "Willamette Valley",
  "https://commons.wikimedia.org/wiki/Special:FilePath/Winery_in_the_Willamette_Valley.jpg?width=500",
  "<b>The Valley.</b><br>Lush alluvial plains. Only a tiny sliver extends into southern WA."
)

# --- 3. Execution ---
wa_data <- get_wa_ecoregions_fixed(ECO_URL, DATA_DIR)

# Join Images
wa_interactive <- wa_data %>%
  left_join(image_lookup, by = "US_L3NAME") %>%
  mutate(
    popup_html = case_when(
      !is.na(img_url) ~ paste0(
        "<div style='font-family: sans-serif; width: 220px;'>",
        "<h4 style='margin: 0 0 8px 0; border-bottom: 2px solid #ccc;'>", US_L3NAME, "</h4>",
        "<img src='", img_url, "' width='100%' style='border-radius: 5px; margin-bottom: 8px;'>",
        "<p style='font-size: 13px; margin: 0; line-height: 1.4;'>", desc, "</p>",
        "</div>"
      ),
      TRUE ~ paste0("<b>", US_L3NAME, "</b><br><i>(No image available)</i>")
    )
  )

# --- 4. Render Final Map ---
# Use "Paired" palette to ensure "Olympic High Peaks" and "North Cascades" get different colors
pal <- colorFactor(palette = "Paired", domain = wa_interactive$US_L3NAME)

leaflet(wa_interactive) %>%
  addProviderTiles(providers$CartoDB.Positron) %>%
  addPolygons(
    fillColor = ~pal(US_L3NAME),
    weight = 1,
    color = "white",
    opacity = 1,
    fillOpacity = 0.7,
    highlightOptions = highlightOptions(
      weight = 3, color = "#666", fillOpacity = 0.9, bringToFront = TRUE
    ),
    popup = ~popup_html,
    label = ~US_L3NAME,
    labelOptions = labelOptions(
      style = list("font-weight" = "bold", "padding" = "3px 8px", "font-size" = "15px"),
      direction = "auto"
    )
  )