primer <- function(dir)
{
	assign("sum.base", sum, .GlobalEnv)
	assign("chol.base", chol, .GlobalEnv)
	assign("t.base", t, .GlobalEnv)
	assign("diag.base", diag, .GlobalEnv)
	assign("forwardsolve.base", forwardsolve, .GlobalEnv)
	assign("backsolve.base", backsolve, .GlobalEnv)
	assign("any.base", any, .GlobalEnv)
	assign("matmult", base:::"%*%", .GlobalEnv)
	assign("lgamma.base", lgamma, .GlobalEnv)
	assign("log.base", log, .GlobalEnv)
	assign("c.base", c, .GlobalEnv)
	assign("subset.base", base:::"[", .GlobalEnv)
	assign("det.base", determinant, .GlobalEnv)
	
	require(spam)
	
	#this reads in exc matrix, the cell information matrix, 
	#and the adjacency matrix C and creates excList and covInfo
	
	########
	#READS IN INFO
	########
	cell <- data.matrix(readRDS(paste(dir, "/cell.dat", sep = "")))
  exc <- data.matrix(readRDS(paste(dir, "/exc.dat", sep = "")))
	C <- data.matrix(readRDS(paste(dir, "/ADJ.dat", sep = "")))
		
	#########
	#EXC LIST
	#########	
	numCells <- dim(cell)[1]
	excList <- list()
	maxVal <- numeric(numCells)

	for(i in seq(1, numCells))
	{
		excList[[i]] <- as.vector(exc[exc[,1] == cell[i, "cell"], "ws"])
		maxVal[i] <- max(c(excList[[i]], cell[i, "thold"]))
	}
	
	#########
	#MAKES W
	#########
	rowSum <- apply(C, 1, sum)
	W <- diag(rowSum) - C
		
	#########
	#COV INFO
	#########
	#this part is the same for all mean structures
	covInfo <- list()
	covInfo$thold <- cell[, "thold"]
	covInfo$nObs <- cell[, "nObs"]
	covInfo$nExc <- cell[, "nExc"]
	covInfo$npy <- cell[, "npy"]
	covInfo$maxVal <- maxVal
	covInfo$W <- as.spam(W)

	ones <- rep(1, numCells)
	zeros <- rep(0, numCells)
	A <- matrix(c(ones, zeros, zeros, zeros, ones, zeros, 
			zeros, zeros, ones), byrow = T, nrow = 3)
	covInfo$A <- A
    	
	#assigns variables to be stored in the DotRData file
	assign("covInfo", covInfo, .GlobalEnv)
	assign("cell", cell, .GlobalEnv)
	assign("excList", excList, .GlobalEnv)
	assign("exc", exc, .GlobalEnv) #deviance uses this
	
	assign("numCells", numCells, .GlobalEnv)
	assign("s", seq(1, numCells), .GlobalEnv)
	assign("nnzlmax", 200000, .GlobalEnv)
	assign("nsubmax", 200000, .GlobalEnv)
}
