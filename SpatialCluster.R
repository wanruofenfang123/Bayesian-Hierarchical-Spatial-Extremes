### Assessing the Risk of Disruption of Wind Turbine Operations ###
### in Saudi Arabia Using Bayesian Spatial Extremes ###
### Wanfang Chen, Stefano Castruccio and Marc G. Genton ###

## Spatial clustering procedure

## install packages
library(flexclust)
library(fields)
library(sp)
library(spam)
library(MCMCpack)
library(spdep)

## Load data
cell <- readRDS("Data/cell.dat")
exc <- readRDS("Data/exc.dat")

## Spatial clustering using weighted k-means
clst <- flexclust::cclust(cbind(cell[,"lon"], cell[,"lat"], cell[,"xi"]), 
                          k=200, weights=c(0.2,0.2,0.6), method="hardcl")
kmeans <- clst@cluster
saveRDS(kmeans, "Data/kmeans.dat")
kmeans <- readRDS("Data/kmeans.dat")

# Define cluster neighbors
ID <- as.numeric(kmeans)
lon <- as.numeric(cell[,"lon"])
lat <- as.numeric(cell[,"lat"])
cell_df <- data.frame(list(lon=lon, lat=lat, id=factor(ID)),row.names =c())
cell_list <- split(cell_df, cell_df$id)
# only want lon-lats in the list, not the names
cell_list <- lapply(cell_list, function(x) { x["id"] <- NULL; x })
# create SpatialPolygons Object, convert coords to polygon
ps <- lapply(cell_list, Polygon)
# add id variable
p1 <- lapply(seq_along(ps), 
             function(i) Polygons(list(ps[[i]]),ID = names(cell_list)[i]))
# create SpatialPolygons object
my_spatial_polys <- SpatialPolygons(p1)
my_spatial_polys_df <- SpatialPolygonsDataFrame(my_spatial_polys, 
                                                data.frame(id = unique(cell_df$id),                                                           row.names = unique(cell_df$id)))
# build a neighbours list based on regions with contiguous boundaries
neighbor <- spdep::poly2nb(my_spatial_polys_df, snap=0.045, queen = FALSE)

# create datasets for cluster neighbors
for (i in 1:200) {
  ind <- which(kmeans%in%c(i,neighbor[[i]]))
  cell_look <- cell[ind,]
  exc_look <- exc[exc[,"cell"]%in%ind,]
  
  # create the adjacency matrix
  sep <- 0.045  # rounded distance of adjacent locations
  d <- as.matrix(dist(cbind(cell_look[,"lon"], cell_look[,"lat"]))) 
  adj_look <- ifelse(round(d,4)==sep,1,0)
  # exclude isolated grid cells
  no_nei <- which(apply(adj_look,1,sum)==0)
  if (length(no_nei)>0) {
    adj_look <- adj_look[-no_nei,-no_nei]
    cell_look <- cell_look[-no_nei,]
    exc_look <- exc[exc[,"cell"]%in%cell_look[,"cell"],]
  }
  
  saveRDS(cell_look,paste0("Data/Cluster_neighbors/", i, "/cell.dat"))
  saveRDS(exc_look,paste0("Data/Cluster_neighbors/", i, "/exc.dat"))
  saveRDS(adj_look,paste0("Data/Cluster_neighbors/", i, "/ADJ.dat"))
}
