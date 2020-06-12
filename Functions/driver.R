driver <- function(numIter = 10, covariates = NULL, 
				   betaPriorMeans = NULL, paramsStart = NULL, 
				   outDir = here::here("results/"))
{	
	recordInt <- 1
	writeOutInt <- 10

	#sets up the covariates to be used
	numCov <- length("covariates")
	ones <- rep(1, numCells)
	zeros <- rep(0, numCells)
	if(is.null(covariates))
	{
	  covInfo$X <- cbind(ones)
	}
	else
	{
	  covInfo$X <- cbind(ones,
	  				rbind(as.matrix(cell[,covariates])))
	}

	#puts information about the priors into covInfo
	#this changes depending on mean structure
	if(is.null(betaPriorMeans))
	{
	  betaPriorMeans <- matrix(c(mean(cell[, "mu"]), 
	    mean(cell[, "lsig"]), mean(cell[, "xi"]),
      rep(0, 3*length(covariates))), 
	    ncol = 3, byrow = T)
	}
	covInfo$betaPriorMeans <- betaPriorMeans  
	covInfo$wishartPrior <- 1/3*diag(c(0.5, 2, 9))
	covInfo$invWishartPrior <- solve(covInfo$wishartPrior)

	#if they are not given, this sets up the initial values 
	#for the parameters, uses the individual MLE estimates
	if(is.null(paramsStart))
	{
		paramsCur <- list()
		paramsCur$mu <- cell[, "mu"]
		paramsCur$sig <- exp(cell[, "lsig"])
		paramsCur$xi <- cell[, "xi"]
		paramsCur$beta <- matrix(betaPriorMeans, ncol=3)
		paramsCur$T <- diag(c(0.5, 2, 9))
		paramsCur$tauSq <- c(4, 200, 500)
		ppMeanTemp <- covInfo$X %*% paramsCur$beta
		paramsCur$ppMean <- cbind(ppMeanTemp[1:numCells,])
		set.seed(2019)
		tempU <- matrix(c(paramsCur$mu, 
					log(paramsCur$sig), paramsCur$xi),
					ncol = 3) - paramsCur$ppMean[,1:3] +
					matrix(rnorm(3*numCells, 0, 
					rep(sqrt(paramsCur$tauSq), 1, 
					each = numCells)), ncol = 3)

		#centering U
		colMeans <- matrix(rep(apply(tempU, 2, mean), 1, 
						each = numCells), ncol = 3)
		paramsCur$U <- tempU - colMeans
	}else{
		paramsCur <- paramsStart
	}

	temp <- det.base(paramsCur$T, logarithm = T)$modulus
	attr(temp, "logarithm") <- NULL
	ppMeanTemp <- covInfo$X %*% paramsCur$beta
	paramsCur$ppMean <- cbind(ppMeanTemp[1:numCells,])
	#although they are stored in paramsCur, the jumpsize$pp
	#values and the tauSq values are not currently changing
	paramsCur$jumpsize <- list()
	paramsCur$jumpsize$pp <- matrix(rep(c(1, .5, .1), 
				numCells), ncol = numCells)	
	paramsCur$tauSq <- c(4, 200, 500)
	paramsCur$count <- list(pp = numeric(numCells))
		
	covInfo$trnsXtauSqX <- kronecker(diag(paramsCur$tauSq),
								t(covInfo$X) %*% covInfo$X)

	#sets up matrices to record output
	out <- list()
	numLines <- numIter %/% recordInt + 1
	out$muRecord <- matrix(nrow = numLines, ncol = numCells)
	out$sigRecord <- out$muRecord
	out$xiRecord <- out$muRecord
#	out$URecord <- matrix(nrow = numLines, ncol = 3*numCells)
#	out$betaRecord <- matrix(nrow = numLines, 
#			ncol = length(paramsCur$beta))
#	out$TRecord <- matrix(nrow = numLines, 
#			ncol = length(paramsCur$T))	
	out$deviance <- vector(length = numLines)

	out$muRecord[1,] <- paramsCur$mu
	out$sigRecord[1,] <- paramsCur$sig
	out$xiRecord[1,] <- paramsCur$xi
#	out$URecord[1,] <- as.vector(paramsCur$U)
#	out$betaRecord[1,] <- as.vector(paramsCur$beta)
#	out$TRecord[1,] <- as.vector(paramsCur$T)
	
	excAll <- rbind(exc)
	out$deviance[1] <- deviance(excAll, paramsCur, covInfo)

#return(paramsCur)
#return(covInfo)


	#########
	# MCMC
	#########
	print("starting MCMC")
	options(warn = 1)
	
	for(i in seq(1, numIter))
	{
#		print(i)

#print('ppParamsUpdate')
		paramsCur <- ppParamsUpdate(excList, paramsCur,
						 covInfo)

#print('UUpdate')
		paramsCur <- UUpdate(paramsCur, covInfo)

#print("betaUpdate")
		paramsCur <- betaUpdate(paramsCur, covInfo)

#print("TUpdate")
		paramsCur <- TUpdate(paramsCur, covInfo)

		
		#records output every recordInt'th step
		if(i %% recordInt == 0)
		{
			out$muRecord[i/recordInt + 1,] <- paramsCur$mu
			out$sigRecord[i/recordInt + 1,] <- paramsCur$sig
			out$xiRecord[i/recordInt + 1,] <- paramsCur$xi
#			out$URecord[i/recordInt + 1,] <- as.vector(
#											paramsCur$U)
#			out$betaRecord[i/recordInt + 1,] <- as.vector(
#											paramsCur$beta)
#			out$TRecord[i/recordInt + 1,] <- as.vector(
#											paramsCur$T)
			out$deviance[i/recordInt + 1] <- deviance(excAll, 
										paramsCur, covInfo)
		}#ends recordInt
	}
	return(out)		
}	
							