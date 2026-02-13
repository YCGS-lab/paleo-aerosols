# =============================================================================
# Map of North American FPD Entities by Fire Proxy Type + Ice Core Data
# =============================================================================


# Clear the environment
rm(list = ls(all.names = TRUE))
cat("\014")

library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(ggrepel)

# setwd("~/Documents/GitHub/paleo-aerosols")
# data_dir <- "~/Documents/GitHub/paleo-aerosols/data"
# output_dir <- "~/Documents/GitHub/paleo-aerosols/output"

base_dir <- "~/Google_Drive/RESEARCH_PROJECTS/Paleo_aerosols/paleo-aerosols"
setwd(base_dir)

data_dir <- file.path(base_dir, "data")
output_dir <- file.path(base_dir, "output")

# Load North America data (pre-filtered using fpdv1_GEOGRAPHY_GLOBAL.csv)

na_entities <- read_csv(file.path(data_dir, "fpdv1_ENTITY_NORTH_AMERICA.csv"),
                        show_col_types = FALSE)

na_data <- read_csv(file.path(data_dir, "fpdv1_DATA_NORTH_AMERICA.csv"),
                    show_col_types = FALSE)

# Load Zhang ice core data
zhang_ice <- read_csv(file.path(data_dir, "zhang_ice_data.csv"),
                      show_col_types = FALSE)

# Load NASEM fire scar data
nafss <- read_csv(file.path(data_dir, "NAFSS_Master_Metadata_V1.1.9000 - NAFSS_Meta.csv"),
                  show_col_types = FALSE)

# Get the primary (most common) fire_proxy for each NA entity
entity_fire_proxy <- na_data %>%
  filter(!is.na(fire_proxy)) %>%
  count(entity_id, fire_proxy) %>%
  group_by(entity_id) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(entity_id, fire_proxy)

# Join with entity locations (include Mexico, filter to >= 14°N for Central America)
na_map_data <- na_entities %>%
  left_join(entity_fire_proxy, by = "entity_id") %>%
  filter(!is.na(lat_dd_north) & !is.na(lon_dd_east) & !is.na(fire_proxy)) %>%
  filter(lat_dd_north >= 14)


# Define Albers Equal Area projection for North America
albers_na <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

# Get North America map (include Greenland region)
north_america <- ne_countries(scale = "medium", returnclass = "sf",
                               continent = "North America") %>%
  st_transform(albers_na)

# Get lakes for context
lakes <- ne_download(scale = "medium", type = "lakes", category = "physical",
                     returnclass = "sf") %>%
  st_transform(albers_na)

# Create GFED region polygons 

# Function to split a polygon at a given longitude
# Uses reasonable longitude bounds to avoid dateline wrapping issues
split_polygon_at_longitude <- function(polygon, split_lon, 
                                       west_limit = -170, 
                                       east_limit = -50,
                                       south_limit = 0,
                                       north_limit = 90) {
  
  # Make sure polygon is valid and in WGS84
  original_crs <- st_crs(polygon)
  
  if(st_crs(polygon)$epsg != 4326) {
    polygon <- st_transform(polygon, 4326)
  }
  
  polygon <- st_make_valid(polygon)
  
  # Create western box (west_limit to split_lon)
  west_box <- st_polygon(list(matrix(c(
    west_limit, south_limit,
    split_lon, south_limit,
    split_lon, north_limit,
    west_limit, north_limit,
    west_limit, south_limit
  ), ncol = 2, byrow = TRUE)))
  west_box <- st_sfc(west_box, crs = 4326)
  
  # Create eastern box (split_lon to east_limit)
  east_box <- st_polygon(list(matrix(c(
    split_lon, south_limit,
    east_limit, south_limit,
    east_limit, north_limit,
    split_lon, north_limit,
    split_lon, south_limit
  ), ncol = 2, byrow = TRUE)))
  east_box <- st_sfc(east_box, crs = 4326)
  
  # Intersect to split
  west_part <- st_intersection(polygon, west_box)
  east_part <- st_intersection(polygon, east_box)
  
  # Transform back to original CRS if needed
  if(!is.na(original_crs) && st_crs(original_crs)$epsg != 4326) {
    west_part <- st_transform(west_part, original_crs)
    east_part <- st_transform(east_part, original_crs)
  }
  
  return(list(west = west_part, east = east_part))
}

# Get geographic data
canada <- ne_countries(scale = "medium", country = "Canada", returnclass = "sf")
alaska <- ne_states(country = "United States of America", returnclass = "sf") %>%
  filter(name == "Alaska")
conus <- ne_states(country = "United States of America", returnclass = "sf") %>%
  filter(!name %in% c("Alaska", "Hawaii")) %>%
  st_union()

# Get Canada + USA as unified polygon
canada_usa <- ne_countries(scale = "medium", continent = "North America", returnclass = "sf") %>%
  filter(name %in% c("Canada", "United States of America")) %>%
  st_union()

# Remove CONUS and Hawaii to get Canada + Alaska only
us_states <- ne_states(country = "United States of America", returnclass = "sf")
conus_hawaii <- us_states %>%
  filter(!name %in% c("Alaska")) %>%
  st_union()

canada_alaska <- st_difference(canada_usa, conus_hawaii)

# Convert to sf with proper CRS
canada_alaska <- st_sf(geometry = canada_alaska, crs = 4326)

conus <- st_sf(geometry = conus, crs = 4326)

# Split both polygons at -100° longitude
canada_alaska_split <- split_polygon_at_longitude(canada_alaska, -100)
canada_alaska_west <- canada_alaska_split$west
canada_alaska_east <- canada_alaska_split$east

conus_split <- split_polygon_at_longitude(conus, -100)
conus_west <- conus_split$west
conus_east <- conus_split$east

# Transform to Albers and simplify
canada_alaska_west <- st_transform(canada_alaska_west, crs = albers_na) %>% st_simplify(dTolerance = 10000)
canada_alaska_east <- st_transform(canada_alaska_east, crs = albers_na) %>% st_simplify(dTolerance = 10000)
conus_west <- st_transform(conus_west, crs = albers_na) %>% st_simplify(dTolerance = 10000)
conus_east <- st_transform(conus_east, crs = albers_na) %>% st_simplify(dTolerance = 10000)

# Convert FPD entities to sf and project
entity_sf <- na_map_data %>%
  st_as_sf(coords = c("lon_dd_east", "lat_dd_north"), crs = 4326) %>%
  st_transform(albers_na)

# Filter ice core data to map extent and convert to sf
zhang_ice_filtered <- zhang_ice %>%
  filter(Long >= -170 & Long <= -10)

ice_sf <- zhang_ice_filtered %>%
  st_as_sf(coords = c("Long", "Lat"), crs = 4326) %>%
  st_transform(albers_na)

# Create projected coordinates for ice core labels
zhang_ice_projected <- zhang_ice_filtered %>%
  st_as_sf(coords = c("Long", "Lat"), crs = 4326) %>%
  st_transform(albers_na) %>%
  mutate(x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2])

# Create projected coordinates for legend annotation (55°N, 35°W)
legend_pos <- st_sfc(st_point(c(-35, 55)), crs = 4326) %>%
  st_transform(albers_na) %>%
  st_coordinates()
legend_x <- legend_pos[1]
legend_y <- legend_pos[2]


# Process NASEM data - filter to valid coordinates and project
nafss_filtered <- nafss %>%
  filter(!is.na(Longitude) & !is.na(Latitude)) %>%
  filter(Longitude >= -170 & Longitude <= -50) %>%
  filter(Latitude >= 14 & Latitude <= 85)

nafss_sf <- nafss_filtered %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(albers_na)


# Add proxy_type labels for unified legend
entity_sf <- entity_sf %>%
  mutate(proxy_type = paste0("Sedimentary charcoal (n=", nrow(entity_sf), ")"))

ice_sf <- ice_sf %>%
  mutate(proxy_type = paste0("Ice core (n=", nrow(ice_sf), ")"))

nafss_sf <- nafss_sf %>%
  mutate(proxy_type = paste0("Fire scar (n=", nrow(nafss_sf), ")"))

# Define the legend order and styling
proxy_levels <- c(
  paste0("Sedimentary charcoal (n=", nrow(entity_sf), ")"),
  paste0("Fire scar (n=", nrow(nafss_sf), ")"),
  paste0("Ice core (n=", nrow(ice_sf), ")")
)

proxy_colors <- c("#1f77b4", "green", "yellow")
names(proxy_colors) <- proxy_levels

proxy_shapes <- c(16, 22, 24)
names(proxy_shapes) <- proxy_levels

# Create proxy levels for first map (without NASEM)
proxy_levels_p1 <- c(
  paste0("Sedimentary charcoal (n=", nrow(entity_sf), ")"),
  paste0("Ice core (n=", nrow(ice_sf), ")")
)

# Create the first map (without NASEM)
p <- ggplot() +
  # Basemap
  geom_sf(data = north_america, fill = "gray90", color = "white", linewidth = 0.3) +
  # Lakes
  geom_sf(data = lakes, fill = "lightblue", color = "lightblue", linewidth = 0.1) +
  # FPD Entity points (circles)
  geom_sf(data = entity_sf,
          aes(fill = proxy_type, shape = proxy_type),
          color = "black", size = 2.5, alpha = 0.75, stroke = 0.3) +
  # Ice core points (triangles)
  geom_sf(data = ice_sf,
          aes(fill = proxy_type, shape = proxy_type),
          color = "black", size = 4, stroke = 1) +
  # Ice core labels
  geom_text_repel(data = zhang_ice_projected,
                  aes(x = x, y = y, label = SiteID),
                  size = 2.8, fontface = "bold",
                  color = "black",
                  bg.color = "white", bg.r = 0.1,
                  box.padding = 0.4,
                  point.padding = 0.3,
                  segment.color = "gray50",
                  segment.size = 0.3,
                  max.overlaps = 20,
                  min.segment.length = 0) +
  # Unified scales
  scale_fill_manual(
    values = c("#1f77b4", "yellow")[1:2],
    breaks = proxy_levels_p1,
    name = "Fire proxies"
  ) +
  scale_shape_manual(
    values = c(21, 24)[1:2],
    breaks = proxy_levels_p1,
    name = "Fire proxies"
  ) +
  # Use Albers Equal Area projection with limits to zoom in on continent
  coord_sf(crs = albers_na,
           xlim = c(-4500000, 3500000),
           ylim = c(200000, 7200000),
           expand = FALSE) +
  # Theme
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  guides(fill = guide_legend(override.aes = list(size = 4)),
         shape = guide_legend(override.aes = list(size = 4))) +
  labs(
    title = "Paleofire + Ice Core Data",
    subtitle = paste0("Sedimentary charcoal: ", nrow(entity_sf), " (Harrison et al. 2022) | 
                      Ice cores (Zhang et al. 2024)"),
    x = NULL,
    y = NULL
  )

# Save the map
ggsave(file.path(output_dir, "map_north_america_with_ice.png"),
       plot = p, width = 14, height = 10, dpi = 300)

ggsave(file.path(output_dir, "map_north_america_with_ice.pdf"),
       plot = p, width = 14, height = 10)


# =============================================================================
# Create second map with NASEM fire scar data
# =============================================================================

# Create proxy levels for second map (with NASEM)
proxy_levels_p2 <- c(
  paste0("Sedimentary charcoal (n=", nrow(entity_sf), ")"),
  paste0("Fire scar (n=", nrow(nafss_sf), ")"),
  paste0("Ice core (n=", nrow(ice_sf), ")")
)

p2 <- ggplot() +
  # Basemap
  geom_sf(data = north_america, fill = "gray90", color = "white", linewidth = 0.3) +
  
  # Lakes
  geom_sf(data = lakes, fill = "lightblue", color = "lightblue", linewidth = 0.1) +
  
  # Add the four bounding polygons
  geom_sf(data = canada_alaska_west, fill = NA, color = "#B8A38E", linewidth = 1) +
  geom_sf(data = canada_alaska_east, fill = NA, color = "#B2A9A9", linewidth = 1) +
  geom_sf(data = conus_west, fill = NA, color = "#A6B5A6", linewidth = 1) +
  geom_sf(data = conus_east, fill = NA, color = "#A4B0C1", linewidth = 1) +
  
  # NASEM fire scar points (small green squares)
  geom_sf(data = nafss_sf,
          aes(fill = proxy_type, shape = proxy_type),
          color = "darkgreen", size = 1.5, alpha = 0.6, stroke = 0.3) +
  # FPD Entity points (circles)
  geom_sf(data = entity_sf,
          aes(fill = proxy_type, shape = proxy_type),
          color = "black", size = 2.5, alpha = 0.75, stroke = 0.3) +
  # Ice core points (triangles)
  geom_sf(data = ice_sf,
          aes(fill = proxy_type, shape = proxy_type),
          color = "black", size = 4, stroke = 1) +

  # Ice core labels
  geom_text_repel(data = zhang_ice_projected,
                  aes(x = x, y = y, label = SiteID),
                  size = 2.8, fontface = "bold",
                  color = "black",
                  bg.color = "white", bg.r = 0.1,
                  box.padding = 0.4,
                  point.padding = 0.3,
                  segment.color = "gray50",
                  segment.size = 0.3,
                  max.overlaps = 20,
                  min.segment.length = 0) +
  # Unified scales
  scale_fill_manual(
    # values = c("#1f77b4", "green", "yellow"),
    values = c("#7B6868", "#3C7C4F", "#479EEB"),
    breaks = proxy_levels_p2,
    name = "Fire proxies"
  ) +
  scale_shape_manual(
    values = c(21, 22, 24),
    breaks = proxy_levels_p2,
    name = "Fire proxies"
  ) +
  # Use Albers Equal Area projection with limits to zoom in on continent
  coord_sf(crs = albers_na,
           xlim = c(-4500000, 3500000),
           ylim = c(-1500000, 7400000),
           expand = FALSE) +
  # Theme
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
    # legend.position = "right",
    legend.position = c(0.5, 0.95),
    legend.justification = c(0.5, 1),
    legend.background = element_rect(fill = NA, color = NA),
    legend.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  guides(fill = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 1),
         shape = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 1)) +
  labs(
    # title = "Fire Scar, Paleofire, & Ice Core Data",
    # subtitle = paste0("Sedimentary charcoal: ", nrow(entity_sf), " (Harrison et al. 2022) | ",
    #                   "Fire scars: ", nrow(nafss_sf), " (NASEM) | ",
    #                   "Ice cores (Zhang et al. 2024)"),
    x = NULL,
    y = NULL
  )

print(p2)

# Save the map with NASEM
ggsave(file.path(output_dir, "map_north_america_with_ice_nasem.png"),
       plot = p2, width = 14, height = 10, dpi = 600)

ggsave(file.path(output_dir, "map_north_america_with_ice_nasem.pdf"),
       plot = p2, width = 14, height = 10)


# # Print ice core site locations
# print(zhang_ice)



