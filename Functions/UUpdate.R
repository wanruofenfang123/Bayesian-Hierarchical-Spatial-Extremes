UUpdate <- function(paramsCur, covInfo)
{
	theta <- matrix(c(paramsCur$mu, log(paramsCur$sig), paramsCur$xi), ncol = 3)
	tauSqMtx <- matrix(rep(paramsCur$tauSq, 1, each=numCells), ncol = 3)
	bHere <- as.vector(tauSqMtx * (theta - paramsCur$ppMean))
	QHere <- diag.spam(rep(paramsCur$tauSq,1,each=numCells)) 
	        + kronecker(paramsCur$T, covInfo$W)
	
	gotOutput <- F
	Upper <- chol.spam(QHere)	
	#note: because Upper is a spam object, you can use
	#forward solve with it, even though its not lower-triangular
	w <- forwardsolve(Upper, bHere)
	mu <- backsolve(Upper, w)
	z <- rnorm(length(bHere))
	v <- backsolve(Upper, z)
	UTemp <- matrix(mu + v, ncol = 3)
	
	#this forces the condition that U_mu, U_sig, U_xi all
	#sum to zero. Follows R&H p37
	A <- covInfo$A
	temp <- forwardsolve(Upper, t(A))
	V <- backsolve(Upper, temp)
	W <- A %*% V
	UpperW <- chol.base(W)
	temp2 <- forwardsolve(UpperW, t(V))
	Urh <- backsolve(UpperW, temp2)
	c <- apply(UTemp, 2, sum)
	UFinal <- UTemp - matrix(t(Urh)%*%c, ncol = 3)

	paramsCur$U <- UFinal[,1:3]

	return(paramsCur)	
}