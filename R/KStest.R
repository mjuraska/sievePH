### Kolmogorov-Smirnov-type test of conditional independence between T and V given Z

####################################################################################
### Description: R functions to perform the test of independence between T and V ###
###              based on the Huang-Louis estimator for the joint cdf of (T,V)   ###
###              (Huang and Louis, Biometrika, 1998)                             ###
###              (the most recent update of 'estCDF')                            ###
### Author:      Michal Juraska                                                  ###
### Date:        November 19, 2011                                               ###
####################################################################################

# T: survival time
# V: mark
# C: censoring time
# data: replicates of [X=min(T,C), Delta=I(T<=C), Y=U*Delta]

# input: data[size,3] each row corresponds to a data point (X,Delta,Y)
#        tu[,2] each row corresponds to the coordinates (t,u) of a point
# output: a vector of the same length as tu, cdf values at (t,u)


### Internal Functions

# 'estCDF' returns a vector containing the estimated CDF values at each (T, V) in 'timeMark'.
# 'timeMark' is a dataframe where each row is a case (i.e., failure) and columns 'T' and 'V' 
# specify the failure times and mark values, respectively. If calculating a marginal 
# distribution, the values of one variable are set to infinity.
# 'data' is a dataframe where each row is a subject. The following columns must be included:
#   'X' for the right-censored failure times (i.e., min of failure time and censoring time)
#   'd' for the failure indicator, 1 if failure, 0 if censored
#   'Y' for the values of the mark variable, set to 0 if no mark is observed
estCDF <- function(timeMark,data){
  ord <- order(data$X,-data$d)
  size <- length(data$X)
  prob <- numeric(size)
  surv <- 1
  for (i in 1:size){
    if (data[ord[i],"d"]==1){ prob[ord[i]] <- surv/(size-i+1) }
    surv <- surv-prob[ord[i]]
  }
  cdf <- sapply(1:NROW(timeMark), function(j){
    sum(prob[data$X <= timeMark[j,"T"] & data$Y <= timeMark[j,"V"]],na.rm=TRUE)
  })
  return(cdf)
}

# 'testStat' returns the test statistic
# 'timeMark' is a dataframe where each row is a case (i.e., failure) and columns 'T' and 'V' 
# specify the failure times and mark values, respectively
# 'data' is a dataframe where each row is a subject. The following columns must be included:
#   'X' for the right-censored failure times
#   'd' for the failure indicator, 1 if failure, 0 if censored
#   'Y' for the values of the mark variable, set to 0 if no mark is observed
testStat <- function(timeMark, data){
  timeMark <- na.omit(timeMark)
  return(max(abs(estCDF(timeMark,data) - estCDF(cbind(timeMark$T,Inf),data)*
                         estCDF(cbind(Inf,timeMark$V),data[data$d==1,]))))
}


#' Kolmogorov-Smirnov-Type Test of Conditional Independence Between Time to Failure 
#' and Mark Variable Given Treatment Group
#' 
#' \code{KStest} performs a Komogorov-Smirnov-type test on the independence of the 
#' time to failure, T, and the distribution of the mark variable, V, given the 
#' treatment group, Z. The independence of T and V given Z is a necessary assumption 
#' for parameter identifiability in the time-independent density ratio model. 
#' The function uses a bootstrap algorithm to obtain and return the p-value for the 
#' test of independence. 
#' 
#' @param data a dataframe where each row is a study subject in a given treatment group 
#' and the following columns are included: \code{X} specifying the right-censored failure 
#' times, \code{d} specifying the failure indicator (1 if failure, 0 if censored), and 
#' \code{Y} specifying the values of the mark variable (0 if no mark is observed).
#' @param iter the number of bootstrap iterations to be performed. Default is 1000.
#' 
#' @details 
#' The null hypothesis is independence between T and V given Z, and the test statistic 
#' is the supremum of the distance between the joint cdf of (T,V) given Z and the product of 
#' the conditional marginal cdfs of T and V. Estimation of the joint cdf is based on the 
#' nonparametric maximum likelihood estimator developed by Huang and Louis (1998).  
#' The estimated marginal cdf of T is one minums the Kaplan-Meier estimator of the 
#' conditional survival function of T, and the estimated marginal cdf of V is the 
#' empirical density function of the observed values of V. A bootstrap algorithm 
#' is used to estimate the critical values for the distribution under the null. 
#' 
#' @return Returns the p-value for the Kolmogorov-Smirnov test of independence between 
#' T and V given Z. 
#' 
#' @examples 
#' 
#' @export
KStest <- function(data, iter=1000){
  T <- testStat(data[data$d==1,c("X","Y")], data)
  n <- NROW(data)
  resamp <- matrix(sample(1:n, n*iter, replace=TRUE), n, iter)
  obs.mark <- as.vector(data[data$d==1,3])
  bootT <- sapply(1:iter, function(j){
    bootData <- data[resamp[,j], c("X","d")]
    if (sum(bootData$d) > 1){
      m <- sum(bootData$d)
      bootMark <- numeric(n)
      bootMark[which(bootData$d==1)] <- sample(obs.mark, m, replace=TRUE)
      bootData <- cbind(bootData, bootMark)
      tstat <- testStat(bootData[bootData$d==1,c(1,3)], bootData)
    } else {
      tstat <- NULL
    }
    tstat
  })
  if (is.list(bootT)){ bootT <- do.call("c", lapply(bootT,"[",1)) }  # if bT is a list, concatenate into a vector
  bootPval <- mean(bootT >= T)
  return(bootPval)
}