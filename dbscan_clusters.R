install.packages("dbscan")
library(dbscan)
install.packages("geosphere")
library(geosphere)

#load data
geo_dat <- read.csv("for_geodata.csv")

#remove missing values
clean_geo <- na.omit(geo_dat)

#latitude and longitude into a matrix
coords <- as.matrix(clean_geo[, c("Long", "Lat")])

#calculate distance matrix (metres)
dist_matrix <- distm(coords, fun = distHaversine)


#smaller site clusters 
#run DBSCAN using the distance matrix
db_result <- dbscan(as.dist(dist_matrix), eps = 1000, minPts = 3)

#add clusters back to data
clean_geo$site_cluster <- db_result$cluster

#number of clusters (excluding noise)
num_clusters <- length(unique(clean_geo$site_cluster[clean_geo$site_cluster > 0]))
print(paste("Number of clusters found:", num_clusters))


#larger regional clusters
#run DBSCAN using the distance matrix 
db_result2 <- dbscan(as.dist(dist_matrix), eps = 9500, minPts = 3)

#add clusters back to data
clean_geo$reg_cluster <- db_result2$cluster

#number of clusters (excluding noise)
num_clusters <- length(unique(clean_geo$reg_cluster[clean_geo$reg_cluster > 0]))
print(paste("Number of clusters found:", num_clusters))

#export
install.packages("writexl")
library(writexl)
write_xlsx(clean_geo, "dbscan_clusters.xlsx")

#simple plot
plot(coords, col = clean_geo$reg_cluster + 1L, pch = 20)

#not sure if this works? 
#only numeric data
numeric_data <- clean_geo[sapply(clean_geo, is.numeric)]

#plot distance to kth nearest neighbour to find optimal epsilon 
kNNdistplot(numeric_data, k = 3)



