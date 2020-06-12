betaUpdate <- function(paramsCur, covInfo)
{
	#performs beta update by drawing from the appropriate
	#mv normal dist.
	
	############################
	# Priors must be assigned to the mean vector and covariace
	# matrix on beta.
	############################
	numCovariates <- dim(covInfo$X)[2] - 1	
	kappa <- 1/rep(c(.01, rep(.1, numCovariates)), 3)
	muBeta <- as.vector(covInfo$betaPriorMeans)
	theta <- matrix(c(paramsCur$mu, log(paramsCur$sig),	paramsCur$xi), ncol = 3)
	#takes advantage of kronecker and vec operators 
	thetaMinusU <- theta - rbind(paramsCur$U)
	bHere <- as.vector(t(covInfo$X) %*% thetaMinusU %*% 
				diag(paramsCur$tauSq)) + diag(kappa)%*% muBeta
	QHere <- covInfo$trnsXtauSqX + diag(kappa)

	#right from Rue and Held page 35			
	U <- chol(QHere)
	L <- t(U)
	w <- forwardsolve(L, bHere)
	mu <- backsolve(U, w)
	z <- rnorm(length(bHere))
	v <- backsolve(U, z)
	betaNew <- matrix(mu + v, ncol = 3)
	paramsCur$beta <- betaNew
	ppMeanTemp <- covInfo$X %*% paramsCur$beta
	paramsCur$ppMean <- cbind(ppMeanTemp[1:numCells,])
	
	return(paramsCur)				
}
