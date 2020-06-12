library(spam)
library(MCMCpack)

dir <- here::here("Rcode/Functions")

source(paste(dir, "/driver.R", sep = ""))
source(paste(dir, "/ppParamsUpdate.R", sep = ""))
source(paste(dir, "/UUpdate.R", sep = ""))
source(paste(dir, "/betaUpdate.R", sep = ""))
source(paste(dir, "/TUpdate.R", sep = ""))
source(paste(dir, "/stedLogPrior.R", sep = ""))
source(paste(dir, "/deviance.R", sep = ""))


