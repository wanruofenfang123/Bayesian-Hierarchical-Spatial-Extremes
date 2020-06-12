ppParamsUpdate <- function(excList, paramsCur, covInfo)
{
	ppMean <- paramsCur$ppMean[,1:3]

	for(i in seq(1, numCells))	
	{
		mu <- paramsCur$mu[i]
		sig <- paramsCur$sig[i]
		xi <- paramsCur$xi[i]
		thold <- covInfo$thold[i]
		nObs  <- covInfo$nObs[i]
	 	nExc  <- covInfo$nExc[i]
		npy   <- covInfo$npy[i]	
		maxVal <- covInfo$maxVal[i]

		jump <- paramsCur$jumpsize$pp[,i]

		#draws candidate from uniform distribution
		#centered at current value
		muCand <- mu + jump[1]*(runif(1) - 1/2)
 	  sigCand <- exp(log.base(sig) + 
						jump[2]*(runif(1) - 1/2))
		xiCand <- xi + jump[3]*(runif(1) - 1/2)

		#Metropolis Hastings step (has three pieces)

		#point-process pieces
		ppBot <- sum.base(-nObs/npy * (1 + xi/sig*(thold - mu))^(-1/xi)-nExc*log.base(sig))
		        -sum.base((1+1/xi) * log.base(1 + xi/sig*(excList[[i]] - mu)))
		test1 <- 1 + xiCand/sigCand*(thold - muCand) <= 0
  	test2 <- 1 + xiCand/sigCand*(maxVal - muCand) <= 0
		if(any.base(c(test1, test2)))
		ppTop <- -Inf			
		else		
		ppTop <- sum.base(-nObs/npy*(1+xiCand/sigCand*(thold-muCand))^(-1/xiCand)-nExc*log.base(sigCand))
  	        -sum.base((1+1/xiCand)*(log.base(1 + xiCand/sigCand*(excList[[i]] - muCand))))
		if(length(excList[[i]]) == 0)
		{ppTop <- ppBot <- 0}

		#stedinger prior pieces
  	stedBot <- stedLogPrior(xi)
	  stedTop <- stedLogPrior(xiCand)
	  
		#normal piece
		cellMean <- ppMean[i,] + paramsCur$U[i,]
		cellSD <- 1/sqrt(paramsCur$tauSq)
		normBot <- sum(dnorm(c(mu, log(sig), xi), cellMean, cellSD, log = T))
		normTop <- sum(dnorm(c(muCand, log(sigCand), xiCand), cellMean, cellSD, log = T))

		#accept/reject
		top <- ppTop + stedTop + normTop		
		bot <- ppBot + stedBot + normBot
		test <- top - bot
		r <- log.base(runif(1))
		if(r < test)
		{
			paramsCur$mu[i] <- muCand
			paramsCur$sig[i] <- sigCand
			paramsCur$xi[i] <- xiCand
			paramsCur$count$pp[i] <- paramsCur$count$pp[i] + 1
		}			
	}
	return(paramsCur)
}		
			
		
		
