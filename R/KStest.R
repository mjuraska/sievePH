### Kolmogorov-Smirnov-type test of conditional independence between T and V given Z

####################################################################################
### Description: R functions to perform the test of independence between T and V ###
###              based on the Huang-Louis estimator for the joint cdf of (T,V)   ###
###              (Huang and Louis, Biometrika, 1998)                             ###
###              (the most recent update of 'ftu')                               ###
### Author:      Michal Juraska                                                  ###
### Date:        November 19, 2011                                               ###
####################################################################################

# T: survival time
# U: mark
# C: censoring time
# data: replicates of [X=min(T,C), Delta=I(T<=C), Y=U*Delta]

# input: data[size,3] each row corresponds to a data point (X,Delta,Y)
#        tu[,2] each row corresponds to the coordinates (t,u) of a point
# output: a vector of the same length as tu, cdf values at (t,u)

ftu <- function(tu,data){
  ord <- order(data[,1],-data[,2])
  size <- length(data[,1])
  prob <- numeric(size)
  surv <- 1
  for (i in 1:size){
    if (data[ord[i],2]==1){ prob[ord[i]] <- surv/(size-i+1) }
    surv <- surv-prob[ord[i]]
  }
  sapply(1:NROW(tu), function(j){
    sum(prob[data[,1]<=tu[j,1] & data[,3]<=tu[j,2]],na.rm=TRUE)
  })
}

g <- function(tu,data){
  if (is.vector(tu)) tu <- t(as.matrix(tu)) 
  tu <- na.omit(tu)
  abs(ftu(tu,data) - ftu(cbind(tu[,1],Inf),data)*ftu(cbind(Inf,tu[,2]),data[data[,2]==1,]))
}

### bootstrap to estimate the p-value for the test of independence between T and V
bootpval <- function(data, iter=1000){
  T <- max(g(data[data[,2]==1,c(1,3)],data))
  n <- NROW(data)
  resamp <- matrix(sample(1:n,n*iter,replace=TRUE),n,iter)
  obs.mark <- as.vector(data[data[,2]==1,3])
  bT <- sapply(1:iter, function(j){
    bdata <- data[resamp[,j],1:2]
    if (sum(bdata[,2])>1){
      m <- sum(bdata[,2])
      bmark <- numeric(n)
      bmark[which(bdata[,2]==1)] <- sample(obs.mark,m,replace=TRUE)
      bdata <- cbind(bdata,bmark)
      tstat <- max(g(bdata[bdata[,2]==1,c(1,3)],bdata))
    } else {
      tstat <- NULL
    }
    tstat
  })
  if (is.list(bT)){ bT <- do.call("c",lapply(bT,"[",1)) }
  mean(bT>=T)
}

indTV.pval <- bootpval(cbind(X,V,d))