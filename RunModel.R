### Assessing the Risk of Disruption of Wind Turbine Operations ###
### in Saudi Arabia Using Bayesian Spatial Extremes ###
### Wanfang Chen, Stefano Castruccio and Marc G. Genton ###

## Run the Bayesian hierarchical spatial extremes model 
## with parallel computing

library(parallel)
cl <- makeCluster(25)  # use 25 cores

f <- function (i) {
  
  source("Functions/primer.R")  #reads in data and sets up objects to store output
  primer(paste0("Data/Cluster_neighbors/", i))
  source("Functions/sourcePrograms.R")  
  
  niter <- 10000  # number of iterations
  result_look <- driver(niter, covariates = c("clon","clat","celev"))  #runs niter steps of MCMC
  
  # extract output
  kmeans <- readRDS("Data/kmeans.dat")
  ind <- which(kmeans==i)  
  # retain results only in central cluster
  keep <- which(cell[,"cell"]%in%ind) 
  n <- dim(result_look$muRecord)[1]
  
  # three PP parameters
  mu_out <- result_look$muRecord[(n/2+1):n, keep]
  lsig_out <- log(result_look$sigRecord[(n/2+1):n, keep])
  xi_out <- result_look$xiRecord[(n/2+1):n, keep]
  deviance <- mean(result_look$deviance[(n/2+1):n])
  
  # 30-year return levels
  p <- 30  # return period
  ret30 <- mu_out + exp(lsig_out)/xi_out*((-log(1-1/(p)))^(-xi_out)-1)
  # probability of 30-year return levels exceeding 25 m/s
  prob30 <- apply(ret30, 2, function(x) ifelse(sum(x>25)>0, sum(x>25)/dim(ret30)[1], 0))
  
  # store posteriors means of PP parameters and return levels, 
  # and posterior quantiles of return levels, risk and deviance 
  result <- cbind(ind, apply(mu_out,2,mean), apply(lsig_out,2,mean), 
                  apply(xi_out,2,mean), apply(ret30,2,mean),
                  apply(ret30,2,quantile,probs=0.50), 
                  apply(ret30,2,quantile,probs=0.05),
                  apply(ret30,2,quantile,probs=0.95), 
                  prob30,rep(deviance,length(ind)))
  colnames(result) <- c("cell","mu_pmean","lsig_pmean","xi_pmean","ret30_pmean",
                        "ret30_50q","ret30_5q","ret30_95q", "prob30", "deviance")
  return(result)
}

# parallel computing
result_KSA <- parLapply(cl, X=1:200, f)

