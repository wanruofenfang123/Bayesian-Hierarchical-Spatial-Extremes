# Bayesian-Hierarchical-Spatial-Extremes


### R codes for ###
### Assessing the Risk of Disruption of Wind Turbine Operations ###
### in Saudi Arabia Using Bayesian Spatial Extremes ###
### Wanfang Chen, Stefano Castruccio and Marc G. Genton ###

The original WRF dataset produced by Yip (2018) is proprietary and confidential at the moment. One needs to prepare the necessary cell-wise data information (contained in "cell.dat" and "exc.dat") based on his/her own dataset in order to run the R codes to perform spatial clustering ("SpatialCluster.R") and run the Bayesian Hierarchical model with parallel computing ("RunModel.R").  

"cell.dat" is a data matrix where the columns contain the following information for each of the grid cells: 
cell: label for the grid cell  
lon: longitude
lat: latitude
elev: elevation
thold: threshold values
nObs: total number of observations 
npy: number of observations per year
nExc: number of exceedances
mu: GMLE for the GPD location parameter 
lsig: GMLE for the log GPD scale parameter
xi: GMLE for the GPD shape parameter
clon: rescaled longitude
clat: rescaled latitude
celev: rescaled elevation
The GMLEs can be obtained by using adapted functions in R (e.g., the "fpot" function in the "evd" package) with the Martins and Stedinger prior added to the log likelihood function.

"exc.dat" is a data matrix combining the threshold exceedances for each of the grid cells.
cell: label for the grid cell
ws: wind speed values exceeding the chosen threshold 

The data folder "Cluster_neighbors" contains the corresponding "cell" and "exc" information, as well as the adjacency matrix, "ADJ", for each of the cluster neighbors. "kmeans.dat" contains the grid cell labels for each of the spatial cluster. These datasets can be produced by running the code "SpatialCluster.R". 

"RunModel.R" runs the Bayesian hierarchical spatial extremes model with parallel computing, using functions from the folder "Functions".  


