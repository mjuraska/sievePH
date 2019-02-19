# Kolmogorov-Smirnov-type test of conditional independence between T and V given Z
# R functions to perform the test of independence between T and V based on the Huang-Louis estimator for the joint cdf of (T,V) (Huang and Louis, Biometrika, 1998)

# T: survival time
# U=(U1,U2,...,Uk): multivariate mark, where k is the dimension of the mark
# C: censoring time
# data: replicates of [X=min(T,C), Delta=I(T<=C), Y=U*Delta]

# input: data[size,k+2] each row corresponds to a data point (X,Delta,Y)
#        tu[,k+1] each row corresponds to the coordinates (t,u) of a point
# output: a vector of the same length as tu, cdf values at (t,u)

ftu <- function(tu, data){
  ord <- order(data[, 1], -data[, 2])
  size <- length(data[, 1])
  prob <- numeric(size)
  surv <- 1
  for (i in 1:size){
    if (data[ord[i], 2]==1){ prob[ord[i]] <- surv / (size - i + 1) }
    surv <- surv - prob[ord[i]]
  }
  cdf <- sapply(1:NROW(tu), function(j){
    logical.mark <- sapply(2:NCOL(tu), function(col){ data[, col + 1] <= tu[j, col] })
    return(sum(prob[data[, 1] <= tu[j, 1] & Reduce("&", as.list(as.data.frame(logical.mark)))], na.rm=TRUE))
  })
  return(cdf)
}

g <- function(tu, data){
  if (is.vector(tu)) tu <- t(as.matrix(tu))
  tu <- na.omit(tu)
  return(abs(ftu(tu, data) - ftu(cbind(tu[, 1], matrix(Inf, NROW(tu), NCOL(tu) - 1)), data) * ftu(cbind(Inf, tu[, -1]), data[data[, 2]==1, ])))
}

#' Kolmogorov-Smirnov-Type Test of Conditional Independence between the Time-to-Event
#' and a Multivariate Mark Given Treatment
#'
#' A nonparametric Komogorov-Smirnov-type test of the null hypothesis that the time-to-event \eqn{T} and a possibly multivariate mark \eqn{V} are conditionally independent given treatment \eqn{Z}
#' as described in Juraska and Gilbert (2013). The conditional independence is a necessary assumption for parameter identifiability in the time-independent density ratio model. A bootstrap
#' algorithm is used to compute the p-value.
#'
#' @param data a data frame restricted to subjects in a given treatment group with the following columns (in this order): the observed right-censored time to the event of interest,
#' the event indicator (1 if event, 0 if right-censored), and the mark variable (one column for each component, if multivariate)
#' @param iter the number of bootstrap iterations (1000 by default) used for computing the p-value
#'
#' @details
#' The test statistic is the supremum of the difference between the estimated conditional joint cumulative distribution function (cdf) of \eqn{(T,V)} given \eqn{Z} and the product of
#' the estimated conditional cdfs of \eqn{T} and \eqn{V} given \eqn{Z}. The joint cdf is estimated by the nonparametric maximum likelihood estimator developed by
#' Huang and Louis (1998). The marginal cdf of \eqn{T} is estimated as one minus the Kaplan-Meier estimator for the conditional survival function of \eqn{T}, and the
#' cdf of \eqn{V} is estimated as the empirical cdf of the observed values of \eqn{V}. A bootstrap algorithm is used to compute the p-value.
#'
#' @return Returns the bootstrap p-value from the test of conditional independence between \eqn{T} and \eqn{V} given \eqn{Z}.
#'
#' @references Juraska, M. and Gilbert, P. B. (2013), Mark-specific hazard ratio model with multivariate continuous marks: an application to vaccine efficacy. \emph{Biometrics} 69(2):328–337.
#'
#' Huang, Y. and Louis, T. (1998), Nonparametric estimation of the joint distribution of survival time and mark variables. \emph{Biometrika} 85, 785–798.
#'
#' @examples
#' n <- 500
#' tx <- rep(0:1, each=n/2)
#' tm <- c(rexp(n/2, 0.2), rexp(n/2, 0.2 * exp(-0.4)))
#' cens <- runif(n, 0, 15)
#' eventTime <- pmin(tm, cens, 3)
#' eventInd <- as.numeric(tm <= pmin(cens, 3))
#' mark1 <- ifelse(eventInd==1, c(rbeta(n/2, 2, 5), rbeta(n/2, 2, 2)), NA)
#' mark2 <- ifelse(eventInd==1, c(rbeta(n/2, 1, 3), rbeta(n/2, 5, 1)), NA)
#'
#' # perform the test for a univariate mark in the placebo group
#' testIndepTimeMark(data.frame(eventTime, eventInd, mark1)[tx==0, ], iter=20)
#'
#' # perform the test for a bivariate mark in the placebo group
#' testIndepTimeMark(data.frame(eventTime, eventInd, mark1, mark2)[tx==0, ], iter=20)
#'
#' @export
testIndepTimeMark <- function(data, iter=1000){
  D <- max(g(data[data[, 2]==1, -2], data))
  n <- NROW(data)
  resamp <- matrix(sample(1:n, n*iter, replace=TRUE), n, iter)
  obs.mark <- data[data[, 2]==1, -(1:2)]
  if (is.vector(obs.mark)){ obs.mark <- as.matrix(obs.mark) }
  bD <- sapply(1:iter, function(j){
    bdata <- data[resamp[, j], 1:2]
    if (sum(bdata[, 2]) > 1){
      m <- sum(bdata[, 2])
      resamp.mark <- sample(NROW(obs.mark), m, replace=TRUE)
      bmark <- as.data.frame(matrix(0, n, NCOL(data) - 2))
      bmark[bdata[, 2]==1, ] <- obs.mark[resamp.mark, ]
      bdata <- cbind(bdata, bmark)
      tstat <- max(g(bdata[bdata[, 2]==1, -2], bdata))
    } else {
      tstat <- NULL
    }
    return(tstat)
  })
  if (is.list(bD)){ bD <- do.call("c", lapply(bD, "[", 1)) }
  return(mean(bD >= D))
}
