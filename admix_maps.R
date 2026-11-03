# ADMIXTURE populations within species plotted geographically
# Potential biogeographical barriers also visualised

# load necessary libraries
library(sf)
library(ggplot2)
library(ozmaps)
library(ggspatial)  # for north arrow and scale bar
library(openxlsx)
library(ggpubr)
library(dplyr)
# needs to have columns: sample, long, lat, sp
record_data <- read.xlsx('/Users/abhishek/honours/petaurus/meta/petaurus_DPetas26-11610_meta.xlsx') 

record_data$long <- as.numeric(record_data$long)
record_data$lat <- as.numeric(record_data$lat)

record_data_nor <- filter(record_data, sp == "norfolcensis") |>
  filter(!is.na(admix_pop_k2)) 
record_data_not <- filter(record_data, sp == "notatus") |>
  filter(!is.na(admix_pop_k2)) 
record_data_bre <- filter(record_data, sp == "breviceps") |>
  filter(!is.na(admix_pop_k2)) 

# load australia map
ozmaps::abs_ste -> australia_map

xlims <- c(112, 155)
ylims <- c(-10, -43)

# simplify the australia outline
australia_map_simplified <- st_simplify(australia_map, dTolerance = 2000)

# crop both astralia map and the GPKG data using the defined bounding box (limits)
bbox <- st_bbox(c(xmin = xlims[1], xmax = xlims[2], ymin = ylims[1], ymax = ylims[2]), crs = st_crs(australia_map_simplified))
australia_map_simplified <- st_crop(australia_map_simplified, bbox)
 
#Macleay
# import the GPKG shapefile
gpkg_data_mac <- st_read("/Users/abhishek/Desktop/macleay.gpkg")

# ensure both datasets are in the same coordinate reference system
gpkg_data_mac <- st_transform(gpkg_data_mac, st_crs(australia_map_simplified))

# # crop the GPKG data by the australia map
cropped_data_mac <- st_intersection(gpkg_data_mac, australia_map_simplified)

#Clarence
# import the GPKG shapefile
gpkg_data_clar <- st_read("/Users/abhishek/Desktop/clarence_shapefile.gpkg")

# ensure both datasets are in the same coordinate reference system
gpkg_data_clar <- st_transform(gpkg_data_clar, st_crs(australia_map_simplified))

# # crop the GPKG data by the australia map
cropped_data_clar <- st_intersection(gpkg_data_clar, australia_map_simplified)


# # calculate centroids of each state
# centroids <- st_centroid(australia_map_simplified)
# 
# # filter out tas and act
# centroids <- centroids[!centroids$NAME %in% c("Tasmania", "Australian Capital Territory", "Other Territories"), ]
# 
# centroids$NAME <- c("New South\nWales", "Victoria", "Queensland", "South\nAustralia", 
#                     "Western\nAustralia", "Northern\nTerritory")
# 
# centroids$geometry[centroids$NAME == "New South\nWales"] <- 
#   centroids$geometry[centroids$NAME == "New South\nWales"] - c(1.5, -0.5)

### bound for zoom 

xlims3 <- c(150.7, 152.3)
ylims3 <- c(-29.5,-28.1)

sp_colours <- c('red','blue')
names(sp_colours) <- unique(record_data$admix_pop_k2)

# norfolcensis map
nor_admix_map <- ggplot() +
  geom_sf(data = cropped_data_mac, fill = "green", color = "transparent", alpha = 0.2) +
  geom_sf(data = cropped_data_clar, fill = "purple", color = "transparent", alpha = 0.2) +
  geom_sf(data = australia_map_simplified, fill = "transparent", color = "black", linewidth = 0.4) +
  geom_point(data = record_data_nor, mapping = aes(x = long, y = lat, fill = factor(admix_pop_k2)), size=1,
             shape=21, stroke=0.2) +
  
  coord_sf(xlim = c(140, 154), ylim = c(-39, -27)) +  #for zoomed map, delete for full map
  
  scale_fill_manual(values=sp_colours, na.value = 'transparent') +
  
  # geom_text(data = centroids, size=2.5,
  #           aes(x = st_coordinates(centroids)[,1], 
  #               y = st_coordinates(centroids)[,2], 
  #               label = toupper(NAME)), 
  #           size = 3.5, color = "black", fontface = "bold") +
  theme_void() +
  theme(
    legend.position = "none",
    legend.text = element_text(face = "italic"),
    legend.key.size = unit(0.5, "lines"),
    legend.spacing.y = unit(0.1, "lines"),
    legend.margin = margin(t = -1, unit = "pt"),        # pulls legend closer to plot
    legend.box.margin = margin(t = -15, unit = "pt")     # removes box padding above legend
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(size = 1.5),               # points in legend size 1
      title.position = "top",
      ncol = 2                                    # optional: arrange in 2 cols if many species
    ),
    fill = guide_legend(
      title = "Admixture Population"
    )
  )

nor_admix_map

ggsave('/Users/abhishek/honours/petaurus/outputs/norfolcensis/norfolcensis_admix_mapv2.png', nor_admix_map, width=15, height=15, units='cm',bg='white', dpi = 400)

# notatus map 
not_admix_map <- ggplot() +
  geom_sf(data = cropped_data_mac, fill = "green", color = "transparent", alpha = 0.2) +
  geom_sf(data = cropped_data_clar, fill = "purple", color = "transparent", alpha = 0.2) +
  geom_sf(data = australia_map_simplified, fill = "transparent", color = "black", linewidth = 0.4) +
  geom_point(data = record_data_not, mapping = aes(x = long, y = lat, fill = factor(admix_pop_k2)), size=1,
             shape=21, stroke=0.2) +
  
  coord_sf(xlim = c(140, 154), ylim = c(-39, -27)) +  #for zoomed map, delete for full map
  
  scale_fill_manual(values=sp_colours, na.value = 'transparent') +
  
  # geom_text(data = centroids, size=2.5,
  #           aes(x = st_coordinates(centroids)[,1], 
  #               y = st_coordinates(centroids)[,2], 
  #               label = toupper(NAME)), 
  #           size = 3.5, color = "black", fontface = "bold") +
  theme_void() +
  theme( 
    legend.position = "none",
    legend.text = element_text(face = "italic"),
    legend.key.size = unit(0.5, "lines"),
    legend.spacing.y = unit(0.1, "lines"),
    legend.margin = margin(t = -1, unit = "pt"),        # pulls legend closer to plot
    legend.box.margin = margin(t = -15, unit = "pt")     # removes box padding above legend
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(size = 1.5),               # points in legend size 1
      title.position = "top",
      ncol = 2                                    # optional: arrange in 2 cols if many species
    ),
    fill = guide_legend(
      title = "Admixture Population"
    )
  )

not_admix_map

ggsave('/Users/abhishek/honours/petaurus/outputs/notatus/notatus_admix_mapv2.png', not_admix_map, width=15, height=15, units='cm',bg='white', dpi = 400)

# breviceps map
bre_admix_map <- ggplot() +
  geom_sf(data = cropped_data_mac, fill = "green", color = "transparent", alpha = 0.2) +
  geom_sf(data = cropped_data_clar, fill = "purple", color = "transparent", alpha = 0.2) +
  geom_sf(data = australia_map_simplified, fill = "transparent", color = "black", linewidth = 0.4) +
  geom_point(data = record_data_bre, mapping = aes(x = long, y = lat, fill = factor(admix_pop_k2)), size=1,
             shape=21, stroke=0.2) +
  
  coord_sf(xlim = c(140, 154), ylim = c(-39, -27)) +  #for zoomed map, delete for full map
  
  scale_fill_manual(values=sp_colours, na.value = 'transparent') +
  
  # geom_text(data = centroids, size=2.5,
  #           aes(x = st_coordinates(centroids)[,1], 
  #               y = st_coordinates(centroids)[,2], 
  #               label = toupper(NAME)), 
  #           size = 3.5, color = "black", fontface = "bold") +
  theme_void() +
  theme(
    legend.position = "none",
    legend.text = element_text(face = "italic"),
    legend.key.size = unit(0.5, "lines"),
    legend.spacing.y = unit(0.1, "lines"),
    legend.margin = margin(t = -1, unit = "pt"),        # pulls legend closer to plot
    legend.box.margin = margin(t = -15, unit = "pt")     # removes box padding above legend
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(size = 1.5),               # points in legend size 1
      title.position = "top",
      ncol = 2                                    # optional: arrange in 2 cols if many species
    ),
    fill = guide_legend(
      title = "Admixture Population"
    )
  )

bre_admix_map

ggsave('/Users/abhishek/honours/petaurus/outputs/breviceps/breviceps_admix_mapv2.png', bre_admix_map, width=15, height=15, units='cm',bg='white', dpi = 400)

admix_maps_combined <- ggarrange(bre_admix_map, nor_admix_map,
                                 labels = "AUTO",
                                 font.label = list(face = "plain"),
                                 ncol = 2, nrow = 1,
                                 vjust = 6,
                                 hjust = -0.5)
admix_maps_combined

ggsave("petaurus/outputs/admix_map_panel.png", 
       admix_maps_combined, width = 25, height = 15,
       units = "cm", dpi = 300, bg = 'white')
