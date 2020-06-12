TUpdate <- function(paramsCur, covInfo)
{
	U <- cbind(paramsCur$U)
	
	wishPostOne <- t(U) %*% covInfo$W %*% U
	wishPostMatrix <- wishPostOne + covInfo$invWishartPrior
	
	Tinv <- riwish(4 + numCells, wishPostMatrix)
	Upper <- chol.base(Tinv)
	Temp <- forwardsolve.base(t(Upper), diag(3))
	T <- backsolve.base(Upper, Temp)

	#this was put in because T was turning out to be not 
	#numerically symmetric
	T[2,1] <- T[1,2]
	T[3,1] <- T[1,3]
	T[3,2] <- T[2,3]
	
	paramsCur$T <- T	
	return(paramsCur)
}
	