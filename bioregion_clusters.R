library(sf) 
library(dplyr)
library(readr)

#load lat long coordinate data 
coord_df <- read_csv("for_geodata.csv")

#remove missing values
clean_coord <- na.omit(coord_df)

#data frame converted to an sf spatial object
coord_sf <- st_as_sf(clean_coord, 
                          coords = c("Long", "Lat"), 
                          crs = 4326) # WGS84

#IBRA bioregions shapefile
ibra_regions <- st_read("/Users/abhishek/Downloads/IBRA7_subregions/ibra7_subregions.shp")

#ensure both datasets use the same Coordinate Reference System (CRS)
if(st_crs(coord_sf) != st_crs(ibra_regions)) {
  coord_sf <- st_transform(coord_sf, st_crs(ibra_regions))
}

#disable spherical geometry 
sf_use_s2(FALSE)  #not sure if correct, needed to work 

#spatial join (assigns the attributes of the bioregion to the intersecting points)
points_w_br <- st_join(coord_sf, ibra_regions, join = st_intersects)

#group data by bioregion name
dat_grouped <- points_w_br |>
  group_by(SUB_NAME_7) |>
  summarise(total_points = n()) #points per region

#view/export
print(dat_grouped)

library(writexl)
write_xlsx(dat_grouped, "bioregion_summary.xlsx")
write_xlsx(points_w_br, "bioregion_data.xlsx")
