deviance <- function(exc, paramsCur, covInfo)
{
	#calculates deviance
	#adapted from the code found in ppLogLhood.R
	
	excWS <- exc[, "ws"]
	
  mu  <- c(paramsCur$mu)
  sig <- c(paramsCur$sig)
  xi  <- c(paramsCur$xi)
  nExc <- c(covInfo$nExc)
  muXP <- rep(mu, nExc)
	sigXP <- rep(sig, nExc)
	xiXP <- rep(xi, nExc)
	
  thold <- c(covInfo$thold)
  nObs  <- c(covInfo$nObs)
  nExc  <- c(covInfo$nExc)
  npy   <- c(covInfo$npy)
  maxVal <- c(covInfo$maxVal)
    

  logLhood1 <- sum(-nObs/npy * (1 + xi/sig*(thold - mu))^(-1/xi)  - nExc*log(sig))
  logLhood2 <- -(sum((1+1/xiXP)*log(1 + xiXP/sigXP*(excWS - muXP))))
	deviance <- -2*(logLhood1 + logLhood2)
	
	return(deviance)
}