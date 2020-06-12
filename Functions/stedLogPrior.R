stedLogPrior <- function(xi)
{
  if(any.base(xi < -.5 | xi > .5))
  {
    return(-1e7)
  }
  else
  {
    out1 <- lgamma.base(10) - lgamma.base(6) - lgamma.base(4)
    out2 <- 5*log.base(.5 + xi) + 3*log(.5 - xi)
    return(sum.base(out1 + out2))
  }
}