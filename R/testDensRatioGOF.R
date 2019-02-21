#' Goodness-of-Fit Test of the Validity of a Univariate or Multivariate Mark Density Ratio Model
#'
#' \code{testDensRatioGoF} implements the complete-case goodness-of-fit test of Qin and Zhang (1997) for evaluating the validity of the specified mark density ratio model used for modeling a component of
#' the mark-specific hazard ratio model in Juraska and Gilbert (2013). Multivariate marks are accommodated. Subjects who experienced the event of interest but their mark is missing are discarded.
#'
#' @param mark either a numeric vector specifying a univariate continuous mark or a data frame specifying a multivariate continuous mark.
#' For subjects with a right-censored time-to-event, the value(s) in \code{mark} should be set to \code{NA}.
#' @param tx a numeric vector indicating the treatment group (1 if treatment, 0 if placebo)
#' @param DRcoef a numeric vector of the coefficients \eqn{\phi} in the weight function \eqn{g(v, \phi) = \exp\{\phi^T (1, v)\}} in the density ratio model. If \code{NULL} (default), the maximum profile likelihood estimates (Qin, 1998)
#' of the coefficients are computed.
#' @param DRlambda the Lagrange multiplier in the profile score functions for \eqn{\phi} (that arises by profiling out the nuisance parameter). If \code{NULL} (default), the maximum profile likelihood estimate (Qin, 1998)
#' of the Lagrange multiplier is computed.
#' @param iter the number of bootstrap iterations (1000 by default)
#'
#' @details
#' \code{testDensRatioGoF} performs a goodness-of-fit test for the exponential form of the weight function, i.e., \eqn{g(v, \phi) = \exp\{\phi^T (1, v)\}. Other weight functions are not considered.
#'
#' @return Returns a list containing the following components:
#' \itemize{
#' \item \code{teststat}: the value of the Kolmogorov-Smirnov-type test statistic
#' \item \code{pval}: the bootstrap p-value from the test of validity of the mark density ratio model
#' \item \code{DRcoef}: the input object if different from \code{NULL} or a numeric vector of estimates of coefficients \eqn{\phi} in the weight function \eqn{g(v, \phi)} in the density ratio model
#' \item \code{DRlambda}: the input object if different from \code{NULL} or an estimate of the Lagrange multiplier in the profile score functions for \eqn{\phi}
#' }
#'
#' @references Qin, J., & Zhang, B. (1997). A goodness-of-fit test for logistic regression models based on case-control data. \emph{Biometrika}, 84(3), 609-618.
#'
#' Juraska, M. and Gilbert, P. B. (2013), Mark-specific hazard ratio model with multivariate continuous marks: an application to vaccine efficacy. \emph{Biometrics} 69(2):328-337.
#'
#' Qin, J. (1998), Inferences for case-control and semiparametric two-sample density ratio models. \emph{Biometrika} 85, 619-630.
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
#' # test goodness-of-fit for a univariate mark
#' testDensRatioGOF(mark1, tx, iter=20)
#'
#' # test goodness-of-fit for a bivariate mark
#' testDensRatioGOF(data.frame(mark1, mark2), tx, iter=20)
#'
#' @export
testDensRatioGOF <- function(mark, tx, DRcoef=NULL, DRlambda=NULL, iter=1000){
  if (is.vector(mark)) {
    mark <- mark[eventInd==1]
  } else {
    mark <- mark[eventInd==1,]
  }
  mark <- as.matrix(mark)
  tx <- tx[eventInd==1]
  
  ninf <- length(tx)   ## number of infections
  n0 <- sum(1-tx)      ## number of placebo infections
  n1 <- sum(tx)        ## number of vaccine infections
  if (any(is.null(DRcoef), is.null(DRlambda))){
    param <- densRatio(mark, tx)$coef
    DRcoef <- param[-length(param)]
    DRlambda <- param[length(param)]
  }

  g <- function(mark, DRcoef){ exp(drop(cbind(1,mark) %*% DRcoef)) }
  p <- function(mark, DRcoef, DRlambda, m){ 1/(m*(1+DRlambda*(g(mark,DRcoef)-1))) }
  F0np <- function(mark, tx){
    F0np.vector <- apply(mark, 1, function(mark.i){
      mark.vector <- sapply(1:NCOL(mark), function(col){ mark[tx==0, col] <= mark.i[col] })
      mean(Reduce("&", as.list(as.data.frame(mark.vector))))
    })
    return( F0np.vector )
  }
  F0sp <- function(mark, prob){
    F0sp.vector <- apply(mark, 1, function(mark.i){
      mark.vector <- sapply(1:NCOL(mark), function(col){ mark[, col] <= mark.i[col] })
      sum(prob*ifelse(Reduce("&", as.list(as.data.frame(mark.vector))),1,0))
    })
    return( F0sp.vector )
  }

  f0 <- p(mark,DRcoef,DRlambda,ninf)
  delta <- sqrt(ninf)*max(abs(F0np(mark,tx) - F0sp(mark,f0)))

  N0.ast <- rmultinom(iter,n0,f0)
  N1.ast <- rmultinom(iter,n1,f0*g(mark,DRcoef))
  v0.ast <- lapply(1:iter, function(iter){
    out <- cbind(sapply(1:NCOL(mark), function(col){ rep(mark[, col], N0.ast[, iter])}))
    return( out )
  })
  v1.ast <- lapply(1:iter, function(iter){
    out <- cbind(sapply(1:NCOL(mark), function(col){ rep(mark[, col], N1.ast[, iter])}))
    return( out )
  })
  z.ast <- rep(0:1,c(n0,n1))

  teststat <- sapply(1:iter, function(iter){
    v.ast.iter <- rbind(as.matrix(v0.ast[[iter]]),as.matrix(v1.ast[[iter]]))
    param <- densRatio(v.ast.iter, z.ast)$coef
    DRcoef.ast <- param[-length(param)]
    DRlambda.ast <- param[length(param)]
    sqrt(ninf)*max(abs(F0np(v.ast.iter,z.ast) - F0sp(v.ast.iter,p(v.ast.iter,DRcoef.ast,DRlambda.ast,ninf))))
  })
  pval <- mean(teststat>=delta)
  return(list(teststat=delta, pval=pval, DRcoef=DRcoef, DRlambda=DRlambda))
}
