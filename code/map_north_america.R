# =============================================================================
# Map of North American FPD Entities by Fire Proxy Type + Ice Core Data
# =============================================================================

library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggrepel)

# Set paths (relative to repo root)
  data_dir <- "data"
  output_dir <- "output"


# Load North America data (pre-filtered using fpdv1_GEOGRAPHY_GLOBAL.csv)

na_entities <- read_csv(file.path(data_dir, "fpdv1_ENTITY_NORTH_AMERICA.csv"),
                        show_col_types = FALSE)

na_data <- read_csv(file.path(data_dir, "fpdv1_DATA_NORTH_AMERICA.csv"),
                    show_col_types = FALSE)

cat(sprintf("Loaded %d North American entities and %d data records\n",
            nrow(na_entities), nrow(na_data)))

# Load Zhang ice core data
zhang_ice <- read_csv(file.path(data_dir, "zhang_ice_data.csv"),
                      show_col_types = FALSE)

cat(sprintf("Loaded %d ice core sites from Zhang data\n", nrow(zhang_ice)))

# Get the primary (most common) fire_proxy for each NA entity
entity_fire_proxy <- na_data %>%
  filter(!is.na(fire_proxy)) %>%
  count(entity_id, fire_proxy) %>%
  group_by(entity_id) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(entity_id, fire_proxy)

# Join with entity locations
na_map_data <- na_entities %>%
  left_join(entity_fire_proxy, by = "entity_id") %>%
  filter(!is.na(lat_dd_north) & !is.na(lon_dd_east) & !is.na(fire_proxy)) %>%
  filter(lat_dd_north >= 24)

cat(sprintf("North American entities with fire proxy data: %d\n", nrow(na_map_data)))

# Get North America map (include Greenland region)
north_america <- ne_countries(scale = "medium", returnclass = "sf",
                               continent = "North America")

# Get lakes for context
lakes <- ne_download(scale = "medium", type = "lakes", category = "physical",
                     returnclass = "sf")

# Convert FPD entities to sf
entity_sf <- na_map_data %>%
  st_as_sf(coords = c("lon_dd_east", "lat_dd_north"), crs = 4326)

# Filter ice core data to map extent and convert to sf
zhang_ice_filtered <- zhang_ice %>%
  filter(Long >= -170 & Long <= -10)

ice_sf <- zhang_ice_filtered %>%
  st_as_sf(coords = c("Long", "Lat"), crs = 4326)

cat(sprintf("Ice cores within map extent: %d\n", nrow(ice_sf)))

# Create summary for legend (sorted by count)
fire_proxy_summary <- na_map_data %>%
  count(fire_proxy, sort = TRUE) %>%
  mutate(label = paste0(fire_proxy, " (n=", n, ")"))

entity_sf <- entity_sf %>%
  mutate(fire_proxy_label = factor(fire_proxy,
                                   levels = fire_proxy_summary$fire_proxy,
                                   labels = fire_proxy_summary$label))

# Define colors for fire_proxy
fire_proxy_colors <- c(
  "charcoal" = "#1f77b4",
  "char" = "#ff7f0e",
  "pchar" = "#2ca02c",
  "char_pollen" = "#d62728",
  "pah" = "#9467bd",
  "unspec" = "#8c564b",
  "psug" = "#e377c2",
  "char_pres_abs" = "#7f7f7f",
  "ordinal_char" = "#bcbd22",
  "bpca" = "#17becf",
  "burnt_phyto" = "#aec7e8",
  "d13c_char" = "#ffbb78",
  "orec" = "#98df8a"
)

# Create the map
p <- ggplot() +
  # Basemap
  geom_sf(data = north_america, fill = "gray90", color = "white", linewidth = 0.3) +
  # Lakes
  geom_sf(data = lakes, fill = "lightblue", color = "lightblue", linewidth = 0.1) +
  # FPD Entity points (circles)
  geom_sf(data = entity_sf,
          aes(color = fire_proxy_label),
          shape = 16, size = 2.5, alpha = 0.75) +
  # Ice core points (triangles)
  geom_sf(data = ice_sf,
          color = "black", fill = "yellow",
          shape = 24, size = 4, stroke = 1) +
  # Ice core labels
  geom_text_repel(data = zhang_ice_filtered,
                  aes(x = Long, y = Lat, label = SiteID),
                  size = 2.8, fontface = "bold",
                  color = "black",
                  bg.color = "white", bg.r = 0.1,
                  box.padding = 0.4,
                  point.padding = 0.3,
                  segment.color = "gray50",
                  segment.size = 0.3,
                  max.overlaps = 20,
                  min.segment.length = 0) +
  # Color scale
  scale_color_manual(
    values = setNames(
      fire_proxy_colors[fire_proxy_summary$fire_proxy],
      fire_proxy_summary$label
    ),
    name = "Fire Proxy (FPD)"
  ) +
  # Limit to North America extent (extended to include Greenland)
  coord_sf(xlim = c(-170, -10), ylim = c(24, 85), expand = FALSE) +
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
  guides(color = guide_legend(override.aes = list(size = 4))) +
  # Add annotation for ice cores
  annotate("point", x = -35, y = 55, shape = 24, size = 4,
           fill = "yellow", color = "black", stroke = 1) +
  annotate("text", x = -30, y = 55, label = paste0("Ice core (n=", nrow(ice_sf), ")"),
           hjust = 0, size = 3.5) +
  labs(
    title = "Paleofire + Ice Core Data",
    subtitle = paste0("Sediment cores: ", nrow(entity_sf), " (Harrison et al. 2022) | Ice cores (Zhang et al. 2024)")
  )

# Save the map
ggsave(file.path(output_dir, "map_north_america_with_ice.png"),
       plot = p, width = 14, height = 10, dpi = 300)

ggsave(file.path(output_dir, "map_north_america_with_ice.pdf"),
       plot = p, width = 14, height = 10)

cat("\nMap saved to output/map_north_america_with_ice.png and .pdf\n")

# Print ice core site locations
cat("\nIce core sites:\n")
print(zhang_ice)
